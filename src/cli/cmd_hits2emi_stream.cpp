// STREAMING native merge: N (queryDB,resultDB) chunks sharing one targetDB -> one v7 emi.
// Peak resident memory INDEPENDENT of total hits: ref_of + block table + a bounded
// reorder ring. Byte-IDENTICAL output to the in-RAM merge_hits2emi: same first-encounter
// ref dict, same read-aligned 65536-row block rule, same per-read total-order comparator
// (bit desc, ref asc, orig-row-index asc), same shared detail::encode_block_columns, same
// file assembly. mmseqs DBReader<unsigned int> + Matcher are replaced by the dependency-free
// dart::mmseqs reader/parser (byte-compatible), so dart carries no mmseqs dependency.
//
// Pass1 (serial): scan chunk result DBs, build ref_of (first-encounter), per-read
//   {chunk,entry,qkey,degree}, and the block table.
// Pass2 (parallel, bounded reorder ring): each worker owns a block (a disjoint read-range),
//   expands rows, sorts each read's rows locally, ZSTD-encodes via the shared encoder, and
//   hands {block_id, blob} to one ordered writer that drains in strictly-increasing block_id.
#include "subcommand.hpp"
#include "dart/version.h"
#include "dart/mmseqs_db.hpp"
#include "dart/emi_v7_io.hpp"
#include "dart/damage_alignment_evidence.hpp"
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <map>

namespace dart {
namespace cli {

using namespace dart::emiv7;
using DB = dart::mmseqs::DB;
using AlnRec = dart::mmseqs::AlnRec;

namespace {
using clk = std::chrono::steady_clock;
inline double ms(clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); }
}  // namespace

int cmd_hits2emi_stream(int argc, char* argv[]) {
    if (argc >= 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        printf("usage: hits2emi-stream <targetDB> <out.emi> <threads> <manifest> [d_max=0.0] [lambda=0.3] [ring=0]\n\n"
               "Stream N MMseqs2 (queryDB, resultDB) chunks that share one targetDB into one EMI v7.\n"
               "  manifest  one 'queryDB resultDB' pair per line; <db>.lookup files give the names\n"
               "  d_max, lambda  damage parameters for the per-alignment damage evidence\n"
               "  ring      reorder window in blocks (0 = auto, 3*threads+4)\n");
        return 0;
    }
    if (argc < 5) { printf("usage: hits2emi-stream <targetDB> <out.emi> <threads> <manifest: 'qdb rdb' per line> [d_max=0.0] [lambda=0.3] [ring=0(auto)]\n"); return 2; }
    std::string tdb = argv[1], out = argv[2]; int nthreads = atoi(argv[3]); std::string manifest = argv[4];
    float d_max = (argc > 5) ? (float)atof(argv[5]) : 0.0f;
    float lambda = (argc > 6) ? (float)atof(argv[6]) : 0.3f;
    int ring_arg = (argc > 7) ? atoi(argv[7]) : 0;
    if (nthreads < 1) nthreads = 1;

    std::vector<std::pair<std::string, std::string>> chunks;
    { std::ifstream in(manifest); std::string q, r; while (in >> q >> r) chunks.push_back({q, r}); }
    auto t0 = clk::now();
    DB T(tdb); T.open();
    const size_t targetCard = T.getSize();
    std::vector<std::unique_ptr<DB>> Q(chunks.size()), R(chunks.size());
    for (size_t c = 0; c < chunks.size(); ++c) {
        Q[c] = std::make_unique<DB>(chunks[c].first); Q[c]->open();
        R[c] = std::make_unique<DB>(chunks[c].second); R[c]->open();
    }
    auto t1 = clk::now();

    // ---- PASS 1 (serial, no target fetch): ref dict + per-read info + block table ----
    std::vector<uint32_t> ref_of(targetCard, 0xFFFFFFFFu); std::vector<uint32_t> ref_keys;
    std::vector<uint32_t> read_chunk, read_entry, read_qkey; std::vector<uint32_t> degree;
    uint64_t total_rows = 0;
    { std::vector<AlnRec> hits;
      for (uint32_t c = 0; c < chunks.size(); ++c) {
        size_t Nc = R[c]->getSize();
        for (size_t i = 0; i < Nc; ++i) {
            hits.clear(); dart::mmseqs::parse_results(R[c]->getData(i, 0), hits);
            if (hits.empty()) continue;
            read_chunk.push_back(c); read_entry.push_back((uint32_t)i); read_qkey.push_back(R[c]->getDbKey(i));
            degree.push_back((uint32_t)hits.size());
            for (const auto& x : hits) if (ref_of[x.dbKey] == 0xFFFFFFFFu) { ref_of[x.dbKey] = (uint32_t)ref_keys.size(); ref_keys.push_back(x.dbKey); }
            total_rows += hits.size();
        }
      } }
    const uint32_t nreads = (uint32_t)degree.size();

    // block table: identical rule to Writer (extend past target only to keep a read whole)
    struct Blk { uint32_t first_read, last_read; uint64_t start_row; uint32_t num_rows; };
    std::vector<Blk> blks;
    { uint64_t i = 0; uint32_t rd = 0; const uint32_t target = EMI_V7_BLOCK_TARGET_ROWS;
      while (rd < nreads) {
        uint64_t bstart = i; uint32_t first_read = rd, last_read = rd;
        while (rd < nreads) {
            uint64_t taken = i - bstart;
            if (taken >= target && rd != last_read) break;
            last_read = rd; i += degree[rd]; rd++;
        }
        blks.push_back({first_read, last_read, bstart, (uint32_t)(i - bstart)});
      } }
    const uint32_t num_blocks = (uint32_t)blks.size();
    auto t2 = clk::now();

    // ---- sections known after pass1: degree, read_to_block, annot, name dicts ----
    std::vector<uint32_t> read_to_block(nreads, 0);
    for (uint32_t b = 0; b < num_blocks; ++b) for (uint32_t rd = blks[b].first_read; rd <= blks[b].last_read; ++rd) read_to_block[rd] = b;

    // query-name lookups (query_id -> strand for damage evidence AND read names)
    std::vector<std::vector<std::string>> qlook(chunks.size());
    for (size_t c = 0; c < chunks.size(); ++c) qlook[c] = dart::mmseqs::load_lookup_vec(chunks[c].first + ".lookup", Q[c]->getSize());
    auto tlook = dart::mmseqs::load_lookup_vec(tdb + ".lookup", targetCard);

    std::vector<std::string> read_names(nreads), ref_names(ref_keys.size());
    for (uint32_t r = 0; r < nreads; ++r) { uint32_t c = read_chunk[r], k = read_qkey[r]; read_names[r] = (k < qlook[c].size() && !qlook[c][k].empty()) ? qlook[c][k] : std::to_string(k); }
    for (size_t r = 0; r < ref_keys.size(); ++r) ref_names[r] = ref_keys[r] < tlook.size() ? tlook[ref_keys[r]] : std::to_string(ref_keys[r]);

    // degree section: u16 per read + u32 escape side array (mirror Writer)
    std::vector<char> deg_sec;
    { std::vector<char> deg_raw; std::vector<uint32_t> escapes;
      for (uint32_t r = 0; r < nreads; ++r) {
        if (degree[r] >= EMI_V7_DEGREE_ESCAPE) { uint16_t esc = EMI_V7_DEGREE_ESCAPE; detail::put(deg_raw, esc); escapes.push_back(degree[r]); }
        else { uint16_t d = (uint16_t)degree[r]; detail::put(deg_raw, d); }
      }
      detail::put<uint64_t>(deg_sec, escapes.size());
      for (uint32_t e : escapes) detail::put(deg_sec, e);
      detail::write_chunk(deg_sec, deg_raw.data(), deg_raw.size());
    }
    std::vector<uint32_t>().swap(degree);   // degree is pass1-transient — free before pass2

    std::vector<char> r2b_sec; detail::write_chunk(r2b_sec, read_to_block.data(), (size_t)nreads * sizeof(uint32_t));
    std::vector<uint32_t>().swap(read_to_block);
    std::vector<EmiReadAnnot> annots(nreads);
    std::vector<char> annot_sec; detail::write_chunk(annot_sec, annots.data(), annots.size() * sizeof(EmiReadAnnot));
    std::vector<EmiReadAnnot>().swap(annots);
    std::vector<char> rdict = detail::encode_name_dict(read_names, EMI_V7_NAME_DICT_BLOCK);
    std::vector<char> fdict = detail::encode_name_dict(ref_names, EMI_V7_NAME_DICT_BLOCK);
    std::vector<std::string>().swap(read_names); std::vector<std::string>().swap(ref_names);

    // File order: header | block_dir | degree | r2b | annot | block_chunks | rdict | fdict | footer
    const uint64_t block_dir_off = sizeof(EMIv7Header);
    const uint64_t block_dir_sz = (uint64_t)num_blocks * sizeof(BlockMeta);
    const uint64_t deg_off = block_dir_off + block_dir_sz;
    const uint64_t r2b_off = deg_off + deg_sec.size();
    const uint64_t annot_off = r2b_off + r2b_sec.size();
    const uint64_t chunks_off = annot_off + annot_sec.size();   // block payload starts here

    // ---- PASS 2 pipeline: bounded reorder ring, workers encode, one writer drains ----
    std::vector<BlockMeta> blockmeta(num_blocks);
    std::mutex mtx; std::condition_variable cv_worker, cv_writer;
    std::map<uint32_t, std::vector<char>> ready;
    uint32_t next_write = 0;
    std::atomic<uint32_t> next_claim{0};
    const uint32_t window = ring_arg > 0 ? (uint32_t)ring_arg : (uint32_t)std::min<uint64_t>(num_blocks ? num_blocks : 1, (uint64_t)nthreads * 3 + 4);

    auto encode_block = [&](uint32_t b, int tid) -> std::vector<char> {
        const Blk& blk = blks[b];
        std::vector<Row> rows; rows.reserve(blk.num_rows);
        std::vector<AlnRec> hits; std::vector<Row> tmp; std::vector<uint32_t> perm;
        for (uint32_t rd = blk.first_read; rd <= blk.last_read; ++rd) {
            uint32_t c = read_chunk[rd], entry = read_entry[rd], qkey = read_qkey[rd];
            hits.clear(); dart::mmseqs::parse_results(R[c]->getData(entry, tid), hits);
            const char* qseq = Q[c]->getData(qkey, tid);
            const std::string& qid = (qkey < qlook[c].size()) ? qlook[c][qkey] : std::string();
            const uint32_t deg = (uint32_t)hits.size();
            tmp.assign(deg, Row{});
            for (uint32_t h = 0; h < deg; ++h) {
                const auto& x = hits[h]; const char* tseq = T.getData((size_t)x.dbKey, tid);
                Row& row = tmp[h]; row.read_idx = rd; row.ref_idx = ref_of[x.dbKey];
                row.bit_score = (float)x.score; row.evalue_log10 = (x.eval > 0) ? (float)std::log10(x.eval) : -300.0f;
                row.identity_q = (uint16_t)(std::min(std::max(x.seqId, 0.0f), 1.0f) * 65535.0f);
                row.aln_len = (uint16_t)std::min<unsigned>(x.alnLength, 65535u);
                row.qstart = (uint16_t)(x.qStartPos + 1); row.qend = (uint16_t)(x.qEndPos + 1);
                row.tstart = (uint16_t)(x.dbStartPos + 1); row.tend = (uint16_t)(x.dbEndPos + 1); row.tlen = (uint16_t)x.dbLen; row.qlen = (uint16_t)x.qLen;
                int qp = x.qStartPos, tp = x.dbStartPos;
                for (char ch : x.backtrace) { if (ch == 'M') { row.qaln += qseq[qp++]; row.taln += tseq[tp++]; } else if (ch == 'I') { row.qaln += qseq[qp++]; row.taln += '-'; } else if (ch == 'D') { row.qaln += '-'; row.taln += tseq[tp++]; } }
                auto dmg = dart::compute_alignment_damage_evidence(qid, row.qaln, row.taln, row.qstart, row.qlen, d_max, lambda);
                row.dmg_k = dmg.k_hits; row.dmg_m = dmg.m_sites; row.dmg_ll_a = dmg.ll_ancient; row.dmg_ll_m = dmg.ll_modern;
            }
            // per-read local sort == global total-order restricted to this read:
            // (bit desc, ref asc, original hit-index asc). read_idx equal within a read.
            perm.resize(deg); std::iota(perm.begin(), perm.end(), 0u);
            std::sort(perm.begin(), perm.end(), [&](uint32_t a, uint32_t b2) {
                const Row& ra = tmp[a]; const Row& rb = tmp[b2];
                if (ra.bit_score != rb.bit_score) return ra.bit_score > rb.bit_score;
                if (ra.ref_idx != rb.ref_idx) return ra.ref_idx < rb.ref_idx;
                return a < b2;
            });
            for (uint32_t k = 0; k < deg; ++k) rows.push_back(std::move(tmp[perm[k]]));
        }
        std::vector<char> block_bytes; detail::encode_block_columns(rows.data(), rows.size(), block_bytes);
        return block_bytes;
    };

    auto worker = [&](int tid) {
        while (true) {
            uint32_t b = next_claim.fetch_add(1);
            if (b >= num_blocks) break;
            { std::unique_lock<std::mutex> lk(mtx); cv_worker.wait(lk, [&] { return b < next_write + window; }); }
            std::vector<char> bytes = encode_block(b, tid);
            { std::unique_lock<std::mutex> lk(mtx); ready.emplace(b, std::move(bytes)); cv_writer.notify_one(); }
        }
    };

    std::ofstream ofs(out, std::ios::binary);
    if (!ofs) { fprintf(stderr, "cannot create %s\n", out.c_str()); return 1; }
    // reserve the front region (header + block_dir + degree + r2b + annot); rewritten at the end
    { std::vector<char> zeros(chunks_off, 0); ofs.write(zeros.data(), zeros.size()); }

    uint64_t running = 0;
    std::thread writer([&] {
        for (uint32_t bid = 0; bid < num_blocks; ++bid) {
            std::vector<char> bytes;
            { std::unique_lock<std::mutex> lk(mtx); cv_writer.wait(lk, [&] { return ready.count(bid) != 0; }); bytes = std::move(ready[bid]); ready.erase(bid); }
            BlockMeta bm{}; bm.start_row = blks[bid].start_row; bm.num_rows = blks[bid].num_rows;
            bm.first_read_idx = blks[bid].first_read; bm.last_read_idx = blks[bid].last_read;
            bm.chunks_offset = chunks_off + running;
            blockmeta[bid] = bm;
            ofs.write(bytes.data(), bytes.size()); running += bytes.size();
            { std::unique_lock<std::mutex> lk(mtx); next_write = bid + 1; cv_worker.notify_all(); }
        }
    });
    std::vector<std::thread> pool; for (int t = 0; t < nthreads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();
    writer.join();
    auto t3 = clk::now();

    // tail: name dicts + footer
    const uint64_t chunks_sz = running;
    const uint64_t rdict_off = chunks_off + chunks_sz;
    const uint64_t fdict_off = rdict_off + rdict.size();
    ofs.write(rdict.data(), rdict.size());
    ofs.write(fdict.data(), fdict.size());
    EMIv7Footer ftr{}; ofs.write(reinterpret_cast<const char*>(&ftr), sizeof(ftr));

    // front: header + block_dir + degree + r2b + annot (exactly chunks_off bytes)
    EMIv7Header h{};
    h.flags = EMI_V7_FLAG_SORTED_BY_READ | EMI_V7_FLAG_READ_ALIGNED_BLKS | EMI_V7_FLAG_DEGREE_CSR |
              EMI_V7_FLAG_SELF_CONTAINED | EMI_V7_FLAG_HAS_ALN_STRINGS | EMI_V7_FLAG_COMPRESSED_DICT;
    h.num_alignments = total_rows; h.num_reads = nreads; h.num_refs = ref_keys.size();
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
    ofs.write(deg_sec.data(), deg_sec.size());
    ofs.write(r2b_sec.data(), r2b_sec.size());
    ofs.write(annot_sec.data(), annot_sec.size());
    ofs.close();
    auto t4 = clk::now();

    std::ifstream es(out, std::ios::ate | std::ios::binary); size_t esz = es.tellg();
    printf("STREAM chunks=%zu threads=%d ring=%u: open=%.0f pass1=%.0f expand+write=%.0f finalize=%.0f TOTAL=%.0fms rows=%llu reads=%u refs=%zu blocks=%u emi=%.2fMB\n",
           chunks.size(), nthreads, window, ms(t0, t1), ms(t1, t2), ms(t2, t3), ms(t3, t4), ms(t0, t4),
           (unsigned long long)total_rows, nreads, ref_keys.size(), num_blocks, esz / 1e6);
    return 0;
}

namespace {
    struct Hits2EmiStreamRegistrar {
        Hits2EmiStreamRegistrar() {
            SubcommandRegistry::instance().register_command(
                "hits2emi-stream",
                "Stream-merge N (queryDB,resultDB) chunks sharing one targetDB into one .emi",
                cmd_hits2emi_stream, 52);
        }
    };
    static Hits2EmiStreamRegistrar registrar;
}

}  // namespace cli
}  // namespace dart
