// emi-merge: union N resolved v7 emis (same reads, possibly OVERLAPPING target sets) into one
// v7 emi for JOINT EM. Reads are unioned BY NAME (a read present in several inputs, or duplicated
// within one input, becomes ONE global read that gathers ALL its alignments); refs are unioned
// BY NAME (a target present in several inputs gets ONE global id). Each global read's
// gathered rows are then DEDUPED by target: duplicate (read,target) pairs — from duplicated read
// names or from multiple HSPs to the same target — collapse to the single best-bit_score row, so
// joint EM never double-counts a (read,target) membership. Surviving rows keep the per-read total
// order (bit desc, ref asc); alignment payload (bit/eval/identity/coords/damage/qaln/taln) and
// per-read annots pass through byte-for-byte — only read_idx and ref_idx are remapped.
//
// Streaming, like hits2emi-stream: production emis are 130-200GB, so rows are never all held
// resident. Pass1 (serial) builds the global read/ref dicts and an UPPER-BOUND per-read degree
// (sum over occurrences) that only partitions reads into blocks. Pass2 (parallel, bounded reorder
// ring) pulls each block's rows on demand via Reader::rows_for_read, dedups per (read,target),
// encodes via detail::encode_block_columns, records the true post-dedup degree, and a single
// writer drains in block order. The degree section, block row counts and num_alignments are
// therefore written AFTER pass2 (the reader is offset-driven, so section order is free).
//
// ceiling: pass2 assumes the inputs share read order (they do when both searches run on the
//   same query DB); if read orders diverge, input>0 block-cache access is non-monotonic and
//   pass2 slows (still correct). upgrade: sort each input's block reads by local id first.
#include "subcommand.hpp"
#include "dart/version.h"
#include "dart/emi_v7_io.hpp"
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include <numeric>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <map>
#include <unordered_map>

namespace dart {
namespace cli {

using namespace dart::emiv7;

namespace {
using clk = std::chrono::steady_clock;
inline double ms(clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); }
constexpr uint32_t NONE = 0xFFFFFFFFu;

// Flat name storage: name k = arena[off[k] .. off[k+1]).
struct Names { std::vector<char> arena; std::vector<uint64_t> off; };

template <class F>
inline void parallel_for(uint64_t n, int nthreads, F&& f) {
    if (n == 0) return;
    int T = std::max(1, std::min<int>(nthreads, (int)std::min<uint64_t>(n, 1u << 20)));
    std::atomic<uint64_t> next{0};
    const uint64_t chunk = std::max<uint64_t>(1, n / ((uint64_t)T * 8) + 1);
    std::vector<std::thread> pool;
    for (int t = 0; t < T; ++t) pool.emplace_back([&] {
        for (;;) { uint64_t s = next.fetch_add(chunk); if (s >= n) break; uint64_t e = std::min(n, s + chunk); for (uint64_t i = s; i < e; ++i) f(i); }
    });
    for (auto& th : pool) th.join();
}

inline uint64_t hash64(const char* p, size_t n) {           // FNV-1a 64
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; ++i) { h ^= (unsigned char)p[i]; h *= 1099511628211ULL; }
    return h;
}

// Parallel, contiguous-block decode of one name dict into a flat arena. Byte-for-byte the
// same names (and order) as NameDict::name(0..n-1), but never touches the one-block
// thread_local cache — contiguous block ranges per thread make it race-free and cache-safe.
inline Names decode_names(const detail::NameDict& d, int nthreads) {
    Names out; const uint32_t n = d.dict_size(); out.off.assign((size_t)n + 1, 0);
    if (n == 0) return out;
    const uint32_t nblk = d.dict_blocks();
    const int T = std::max(1, std::min<int>(nthreads, (int)nblk));
    std::vector<uint32_t> bnd(T + 1);
    for (int t = 0; t <= T; ++t) bnd[t] = (uint32_t)((uint64_t)nblk * t / T);
    std::vector<std::vector<char>> apart(T);
    std::vector<std::vector<uint64_t>> opart(T);
    std::vector<std::thread> pool;
    for (int t = 0; t < T; ++t) pool.emplace_back([&, t] {
        for (uint32_t b = bnd[t]; b < bnd[t + 1]; ++b) d.decode_block_bytes(b, apart[t], opart[t]);
    });
    for (auto& th : pool) th.join();
    uint64_t abase = 0, nk = 0;
    for (int t = 0; t < T; ++t) { for (uint64_t e : opart[t]) out.off[++nk] = abase + e; abase += apart[t].size(); }
    out.arena.resize(abase);
    uint64_t pos = 0;
    for (int t = 0; t < T; ++t) { std::memcpy(out.arena.data() + pos, apart[t].data(), apart[t].size()); pos += apart[t].size(); }
    return out;
}

}  // namespace

int cmd_emi_merge(int argc, char* argv[]) {
    if (argc < 5) { printf("usage: emi-merge <out.emi> <threads> <in1.emi> <in2.emi> [in3.emi ...]\n"); return 2; }
    std::string out = argv[1];
    int nthreads = atoi(argv[2]); if (nthreads < 1) nthreads = 1;
    std::vector<std::string> inpaths(argv + 3, argv + argc);
    const size_t N = inpaths.size();

    auto t0 = clk::now();
    std::vector<std::unique_ptr<Reader>> in;
    in.reserve(N);
    for (const auto& p : inpaths) {
        auto r = std::make_unique<Reader>(p);
        if (!r->has_alignment_strings())
            { fprintf(stderr, "emi-merge: input lacks alignment strings (need HAS_ALN_STRINGS): %s\n", p.c_str()); return 1; }
        in.push_back(std::move(r));
    }
    const float d_max = in[0]->d_max(), lambda = in[0]->lambda();
    for (size_t i = 1; i < N; ++i)
        if (in[i]->d_max() != d_max || in[i]->lambda() != lambda)
            { fprintf(stderr, "emi-merge: damage params differ (in0 d_max=%g lambda=%g vs in%zu d_max=%g lambda=%g)\n", d_max, lambda, i, in[i]->d_max(), in[i]->lambda()); return 1; }
    auto t1 = clk::now();

    // ---- PASS 1: union reads AND refs BY NAME (first-encounter, keep every read) ----
    // A read name seen in several inputs OR several times within one input is ONE global read that
    // gathers ALL its occurrences' rows (no read and no row is dropped). Name dicts are decoded
    // once in parallel into flat arenas (bypassing the NameDict one-block thread_local cache) and
    // re-encoded in parallel. Duplicate (read,target) pairs are removed in Pass2 (best bit_score
    // wins), so the true per-read degree is only known after Pass2 — here we build an UPPER-BOUND
    // degree (sum over occurrences) purely to partition reads into blocks.
    // merged ref names = union BY NAME across inputs (targets may OVERLAP: the two searches are
    // NOT guaranteed disjoint). ref_gid[i][local] -> global ref id, first-encounter order. Unioning
    // targets is what lets Pass2's per-read (read,target) dedup also collapse CROSS-input duplicate
    // pairs: an overlapping target gets ONE global id, so the same (read,target) from both inputs
    // lands on the same ref_idx in the read's bucket and only the best-bit_score row survives.
    std::vector<Names> fda(N);
    for (size_t i = 0; i < N; ++i) fda[i] = decode_names(in[i]->ref_dict(), nthreads);
    std::vector<std::vector<uint32_t>> ref_gid(N);
    Names refmerged; refmerged.off.push_back(0);
    { std::unordered_map<std::string_view, uint32_t> gid;
      size_t tot = 0; for (size_t i = 0; i < N; ++i) tot += (uint32_t)in[i]->num_refs();
      gid.reserve(tot * 2);
      for (size_t i = 0; i < N; ++i) {
          const uint32_t nf = (uint32_t)in[i]->num_refs();
          ref_gid[i].resize(nf);
          const char* base = fda[i].arena.data(); const uint64_t* off = fda[i].off.data();
          for (uint32_t f = 0; f < nf; ++f) {
              std::string_view sv(base + off[f], off[f + 1] - off[f]);
              auto it = gid.find(sv);
              uint32_t g;
              if (it == gid.end()) {
                  g = (uint32_t)(refmerged.off.size() - 1);
                  refmerged.arena.insert(refmerged.arena.end(), sv.begin(), sv.end());
                  refmerged.off.push_back(refmerged.arena.size());
                  gid.emplace(sv, g);
              } else g = it->second;
              ref_gid[i][f] = g;
          }
      } }
    const uint64_t num_refs = refmerged.off.size() - 1;
    std::vector<Names>().swap(fda);

    // reads: decode each input's read arena once (parallel), then serial first-encounter union.
    std::vector<Names> rda(N);
    for (size_t i = 0; i < N; ++i) rda[i] = decode_names(in[i]->read_dict(), nthreads);
    std::vector<std::vector<uint32_t>> g_of_read(N);
    Names rdmerged; rdmerged.off.push_back(0);
    std::unordered_map<std::string_view, uint32_t> read_gid;
    { size_t tot = 0; for (size_t i = 0; i < N; ++i) tot += (uint32_t)in[i]->num_reads(); read_gid.reserve(tot * 2); }
    for (size_t i = 0; i < N; ++i) {
        const uint32_t nr = (uint32_t)in[i]->num_reads();
        g_of_read[i].resize(nr);
        const char* base = rda[i].arena.data(); const uint64_t* off = rda[i].off.data();
        for (uint32_t r = 0; r < nr; ++r) {
            std::string_view sv(base + off[r], off[r + 1] - off[r]);
            auto it = read_gid.find(sv);
            uint32_t g;
            if (it == read_gid.end()) {
                g = (uint32_t)(rdmerged.off.size() - 1);
                rdmerged.arena.insert(rdmerged.arena.end(), sv.begin(), sv.end());
                rdmerged.off.push_back(rdmerged.arena.size());
                read_gid.emplace(sv, g);
            } else g = it->second;
            g_of_read[i][r] = g;
        }
    }
    const uint32_t nreads = (uint32_t)(rdmerged.off.size() - 1);
    { std::unordered_map<std::string_view, uint32_t> t; std::swap(t, read_gid); }

    // per-input global->local: primary occurrence loc[i][g], plus rare extras (duplicate read
    // names within one input). Upper-bound degree (sum of occurrences) for block partition only;
    // first-owner annot (first occurrence in argv order, read id ascending).
    std::vector<std::vector<uint32_t>> loc(N, std::vector<uint32_t>(nreads, NONE));
    std::unordered_map<uint64_t, std::vector<uint32_t>> dup_extra;   // key = i*nreads + g -> extra local ids
    std::vector<uint32_t> ub_degree(nreads, 0);
    std::vector<EmiReadAnnot> annots(nreads);
    std::vector<uint8_t> annot_set(nreads, 0);
    for (size_t i = 0; i < N; ++i) {
        const uint32_t nr = (uint32_t)in[i]->num_reads();
        for (uint32_t r = 0; r < nr; ++r) {
            uint32_t g = g_of_read[i][r];
            if (loc[i][g] == NONE) loc[i][g] = r;
            else dup_extra[(uint64_t)i * nreads + g].push_back(r);
            ub_degree[g] += in[i]->degree(r);
            if (!annot_set[g]) { annots[g] = in[i]->annot(r); annot_set[g] = 1; }
        }
    }
    std::vector<std::vector<uint32_t>>().swap(g_of_read);
    std::vector<uint8_t>().swap(annot_set);
    std::vector<Names>().swap(rda);

    // block table: identical rule to Writer (extend past target only to keep a read whole)
    struct Blk { uint32_t first_read, last_read; uint64_t start_row; uint32_t num_rows; };
    std::vector<Blk> blks;
    { uint64_t i = 0; uint32_t rd = 0; const uint32_t target = EMI_V7_BLOCK_TARGET_ROWS;
      while (rd < nreads) {
        uint64_t bstart = i; uint32_t first_read = rd, last_read = rd;
        while (rd < nreads) {
            uint64_t taken = i - bstart;
            if (taken >= target && rd != last_read) break;
            last_read = rd; i += ub_degree[rd]; rd++;
        }
        blks.push_back({first_read, last_read, bstart, (uint32_t)(i - bstart)});
      } }
    const uint32_t num_blocks = (uint32_t)blks.size();
    auto t2 = clk::now();

    // ---- sections known after pass1: read_to_block, annot, name dicts (degree AFTER pass2) ----
    std::vector<uint32_t> read_to_block(nreads, 0);
    for (uint32_t b = 0; b < num_blocks; ++b) for (uint32_t rd = blks[b].first_read; rd <= blks[b].last_read; ++rd) read_to_block[rd] = b;
    std::vector<uint32_t>().swap(ub_degree);   // partition done; true degree is computed in pass2

    std::vector<char> r2b_sec; detail::write_chunk(r2b_sec, read_to_block.data(), (size_t)nreads * sizeof(uint32_t));
    std::vector<uint32_t>().swap(read_to_block);
    std::vector<char> annot_sec; detail::write_chunk(annot_sec, annots.data(), (size_t)nreads * sizeof(EmiReadAnnot));
    std::vector<EmiReadAnnot>().swap(annots);
    std::vector<char> rdict = detail::encode_name_dict_arena(rdmerged.arena.data(), rdmerged.off.data(), nreads, EMI_V7_NAME_DICT_BLOCK, nthreads);
    { Names t; std::swap(t, rdmerged); }
    std::vector<char> fdict = detail::encode_name_dict_arena(refmerged.arena.data(), refmerged.off.data(), (uint32_t)num_refs, EMI_V7_NAME_DICT_BLOCK, nthreads);
    { Names t; std::swap(t, refmerged); }

    // Pass2 removes duplicate (read,target) pairs, so the degree section, block row counts and
    // num_alignments are only known AFTER pass2. Reader is offset-driven, so physical order is
    // free: header | block_dir | r2b | annot | block_chunks | degree | rdict | fdict | footer.
    // The reserved prefix [0, chunks_off) = header|block_dir|r2b|annot is back-patched at the end;
    // degree is appended after the chunks once actual per-read counts are known.
    const uint64_t block_dir_off = sizeof(EMIv7Header);
    const uint64_t block_dir_sz = (uint64_t)num_blocks * sizeof(BlockMeta);
    const uint64_t r2b_off = block_dir_off + block_dir_sz;
    const uint64_t annot_off = r2b_off + r2b_sec.size();
    const uint64_t chunks_off = annot_off + annot_sec.size();
    std::vector<uint32_t> actual_degree(nreads, 0);   // filled by pass2 encode_block (disjoint reads)

    // ---- PASS 2 pipeline: bounded reorder ring, workers encode, one writer drains ----
    // encode_block gathers every occurrence of each read (primary loc + dup_extra), removes
    // duplicate (read,target) pairs keeping the best bit_score per ref, then emits the per-read
    // total order (bit desc, ref asc). It records the surviving per-read counts into actual_degree
    // (blocks own disjoint read ranges, so no races) and returns the block's total surviving rows.
    struct Enc { std::vector<char> bytes; uint32_t nrows; };
    std::vector<BlockMeta> blockmeta(num_blocks);
    std::mutex mtx; std::condition_variable cv_worker, cv_writer;
    std::map<uint32_t, Enc> ready;
    uint32_t next_write = 0;
    std::atomic<uint32_t> next_claim{0};
    const uint32_t window = (uint32_t)std::min<uint64_t>(num_blocks ? num_blocks : 1, (uint64_t)nthreads * 3 + 4);

    auto encode_block = [&](uint32_t b) -> Enc {
        const Blk& blk = blks[b];
        const uint32_t nrb = blk.last_read - blk.first_read + 1;
        std::vector<std::vector<Row>> bucket(nrb);   // rows per read, gathered across all occurrences
        // gather one input at a time so each input's block cache is walked monotonically
        for (size_t i = 0; i < N; ++i) {
            const std::vector<uint32_t>& rg = ref_gid[i];   // local ref -> global (union) ref id
            for (uint32_t rd = blk.first_read; rd <= blk.last_read; ++rd) {
                auto& dst = bucket[rd - blk.first_read];
                auto gather = [&](uint32_t lr) {
                    std::vector<Row> rws = in[i]->rows_for_read(lr);
                    for (Row& row : rws) { row.read_idx = rd; row.ref_idx = rg[row.ref_idx]; dst.push_back(std::move(row)); }
                };
                uint32_t lr = loc[i][rd];
                if (lr == NONE) continue;
                gather(lr);
                auto ex = dup_extra.find((uint64_t)i * nreads + rd);   // rare duplicate read names
                if (ex != dup_extra.end()) for (uint32_t r2 : ex->second) gather(r2);
            }
        }
        std::vector<Row> rows; rows.reserve(blk.num_rows);
        uint32_t block_rows = 0;
        std::vector<uint32_t> perm, keep;
        for (uint32_t j = 0; j < nrb; ++j) {
            auto& dst = bucket[j];
            const uint32_t deg = (uint32_t)dst.size();
            perm.resize(deg); std::iota(perm.begin(), perm.end(), 0u);
            // dedup (read,target): group by ref, keep best bit_score (ties -> gather order).
            std::sort(perm.begin(), perm.end(), [&](uint32_t a, uint32_t c) {
                const Row& ra = dst[a]; const Row& rc = dst[c];
                if (ra.ref_idx != rc.ref_idx) return ra.ref_idx < rc.ref_idx;
                if (ra.bit_score != rc.bit_score) return ra.bit_score > rc.bit_score;
                return a < c;
            });
            keep.clear();
            for (uint32_t k = 0; k < deg; ++k)
                if (k == 0 || dst[perm[k]].ref_idx != dst[perm[k - 1]].ref_idx) keep.push_back(perm[k]);
            // per-read total order of survivors: bit desc, ref asc, stable by gather order
            std::sort(keep.begin(), keep.end(), [&](uint32_t a, uint32_t c) {
                const Row& ra = dst[a]; const Row& rc = dst[c];
                if (ra.bit_score != rc.bit_score) return ra.bit_score > rc.bit_score;
                if (ra.ref_idx != rc.ref_idx) return ra.ref_idx < rc.ref_idx;
                return a < c;
            });
            actual_degree[blk.first_read + j] = (uint32_t)keep.size();
            block_rows += (uint32_t)keep.size();
            for (uint32_t idx : keep) rows.push_back(std::move(dst[idx]));
        }
        Enc e; detail::encode_block_columns(rows.data(), rows.size(), e.bytes); e.nrows = block_rows;
        return e;
    };

    auto worker = [&]() {
        while (true) {
            uint32_t b = next_claim.fetch_add(1);
            if (b >= num_blocks) break;
            { std::unique_lock<std::mutex> lk(mtx); cv_worker.wait(lk, [&] { return b < next_write + window; }); }
            Enc e = encode_block(b);
            { std::unique_lock<std::mutex> lk(mtx); ready.emplace(b, std::move(e)); cv_writer.notify_one(); }
        }
    };

    std::ofstream ofs(out, std::ios::binary);
    if (!ofs) { fprintf(stderr, "cannot create %s\n", out.c_str()); return 1; }
    { std::vector<char> zeros(chunks_off, 0); ofs.write(zeros.data(), zeros.size()); }

    uint64_t running = 0;        // byte offset within block_chunks
    uint64_t running_rows = 0;   // cumulative surviving rows => block start_row + num_alignments
    std::thread writer([&] {
        for (uint32_t bid = 0; bid < num_blocks; ++bid) {
            Enc e;
            { std::unique_lock<std::mutex> lk(mtx); cv_writer.wait(lk, [&] { return ready.count(bid) != 0; }); e = std::move(ready[bid]); ready.erase(bid); }
            BlockMeta bm{}; bm.start_row = running_rows; bm.num_rows = e.nrows;
            bm.first_read_idx = blks[bid].first_read; bm.last_read_idx = blks[bid].last_read;
            bm.chunks_offset = chunks_off + running;
            blockmeta[bid] = bm;
            ofs.write(e.bytes.data(), e.bytes.size()); running += e.bytes.size(); running_rows += e.nrows;
            { std::unique_lock<std::mutex> lk(mtx); next_write = bid + 1; cv_worker.notify_all(); }
        }
    });
    std::vector<std::thread> pool; for (int t = 0; t < nthreads; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    writer.join();
    auto t3 = clk::now();

    // degree section from the ACTUAL per-read counts produced by pass2 (post dedup)
    std::vector<char> deg_sec;
    { std::vector<char> deg_raw; std::vector<uint32_t> escapes;
      for (uint32_t r = 0; r < nreads; ++r) {
        if (actual_degree[r] >= EMI_V7_DEGREE_ESCAPE) { uint16_t esc = EMI_V7_DEGREE_ESCAPE; detail::put(deg_raw, esc); escapes.push_back(actual_degree[r]); }
        else { uint16_t d = (uint16_t)actual_degree[r]; detail::put(deg_raw, d); }
      }
      detail::put<uint64_t>(deg_sec, escapes.size());
      for (uint32_t e : escapes) detail::put(deg_sec, e);
      detail::write_chunk(deg_sec, deg_raw.data(), deg_raw.size());
    }
    std::vector<uint32_t>().swap(actual_degree);

    const uint64_t total_rows = running_rows;
    const uint64_t chunks_sz = running;
    const uint64_t deg_off = chunks_off + chunks_sz;
    const uint64_t rdict_off = deg_off + deg_sec.size();
    const uint64_t fdict_off = rdict_off + rdict.size();
    ofs.write(deg_sec.data(), deg_sec.size());
    ofs.write(rdict.data(), rdict.size());
    ofs.write(fdict.data(), fdict.size());
    EMIv7Footer ftr{}; ofs.write(reinterpret_cast<const char*>(&ftr), sizeof(ftr));

    EMIv7Header h{};
    h.flags = EMI_V7_FLAG_SORTED_BY_READ | EMI_V7_FLAG_READ_ALIGNED_BLKS | EMI_V7_FLAG_DEGREE_CSR |
              EMI_V7_FLAG_SELF_CONTAINED | EMI_V7_FLAG_HAS_ALN_STRINGS | EMI_V7_FLAG_COMPRESSED_DICT;
    h.num_alignments = total_rows; h.num_reads = nreads; h.num_refs = num_refs;
    h.d_max = d_max; h.lambda = lambda; h.num_blocks = num_blocks; h.block_target = EMI_V7_BLOCK_TARGET_ROWS;
    h.sections[(int)SectionID::BLOCK_DIR]      = {block_dir_off, block_dir_sz};
    h.sections[(int)SectionID::READ_DEGREE]    = {deg_off, deg_sec.size()};
    h.sections[(int)SectionID::READ_TO_BLOCK]  = {r2b_off, r2b_sec.size()};
    h.sections[(int)SectionID::READ_ANNOT]     = {annot_off, annot_sec.size()};
    h.sections[(int)SectionID::BLOCK_CHUNKS]   = {chunks_off, chunks_sz};
    h.sections[(int)SectionID::NAME_DICT_READ] = {rdict_off, rdict.size()};
    h.sections[(int)SectionID::NAME_DICT_REF]  = {fdict_off, fdict.size()};
    ofs.seekp(0);
    ofs.write(reinterpret_cast<const char*>(&h), sizeof(h));
    ofs.write(reinterpret_cast<const char*>(blockmeta.data()), block_dir_sz);
    ofs.write(r2b_sec.data(), r2b_sec.size());
    ofs.write(annot_sec.data(), annot_sec.size());
    ofs.close();
    auto t4 = clk::now();

    std::ifstream es(out, std::ios::ate | std::ios::binary); size_t esz = es.tellg();
    printf("EMI-MERGE inputs=%zu threads=%d: open=%.0f pass1=%.0f merge+write=%.0f finalize=%.0f TOTAL=%.0fms rows=%llu reads=%u refs=%llu blocks=%u emi=%.2fMB\n",
           N, nthreads, ms(t0, t1), ms(t1, t2), ms(t2, t3), ms(t3, t4), ms(t0, t4),
           (unsigned long long)total_rows, nreads, (unsigned long long)num_refs, num_blocks, esz / 1e6);
    return 0;
}

namespace {
    struct EmiMergeRegistrar {
        EmiMergeRegistrar() {
            SubcommandRegistry::instance().register_command(
                "emi-merge",
                "Union N resolved v7 emis (same reads, disjoint targets) into one for joint EM",
                cmd_emi_merge, 53);
        }
    };
    static EmiMergeRegistrar registrar;
}

}  // namespace cli
}  // namespace dart
