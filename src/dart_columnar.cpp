// dart_columnar.cpp — .dcx writer (write_dcx). Depends only on zstd.
// On-disk contract is dart_columnar.hpp. See that header for the byte layout;
// this file encodes a Batch into size-bounded, zstd-compressed columnar blocks
// with a trailer zone-map index, string dictionaries, and an adler32-checked
// footer. The consumer (Rust reader in convergence) mirrors read_block().

#include "dart/dart_columnar.hpp"

#include <zstd.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <future>
#include <deque>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace dcx {
namespace {

// ---- little-endian primitives (format is LE-only) ----
inline void put_u16(std::vector<uint8_t>& b, uint16_t v) {
    b.push_back(uint8_t(v)); b.push_back(uint8_t(v >> 8));
}
inline void put_u32(std::vector<uint8_t>& b, uint32_t v) {
    for (int i = 0; i < 4; ++i) b.push_back(uint8_t(v >> (8 * i)));
}
inline void put_u64(std::vector<uint8_t>& b, uint64_t v) {
    for (int i = 0; i < 8; ++i) b.push_back(uint8_t(v >> (8 * i)));
}
inline void put_bytes(std::vector<uint8_t>& b, const void* p, size_t n) {
    const auto* c = static_cast<const uint8_t*>(p);
    b.insert(b.end(), c, c + n);
}

// LEB128 unsigned varint (hoki write_varint idiom).
inline void put_varint(std::vector<uint8_t>& b, uint64_t v) {
    while (v >= 0x80) { b.push_back(uint8_t(v) | 0x80); v >>= 7; }
    b.push_back(uint8_t(v));
}
// zigzag for signed deltas (tstart resets per target; query ordinals wander).
inline uint64_t zigzag(int64_t v) { return (uint64_t(v) << 1) ^ uint64_t(v >> 63); }

// adler32 over a byte range (hoki format.hpp idiom), for the trailer index.
uint32_t adler32(const uint8_t* p, size_t n) {
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < n; ++i) { a = (a + p[i]) % 65521; b = (b + a) % 65521; }
    return (b << 16) | a;
}

// zstd frame with content checksum enabled (free XXH64 on the block payload).
std::vector<uint8_t> zstd_frame(const std::vector<uint8_t>& raw, int level) {
    ZSTD_CCtx* c = ZSTD_createCCtx();
    if (!c) throw std::runtime_error("dcx: ZSTD_createCCtx failed");
    ZSTD_CCtx_setParameter(c, ZSTD_c_compressionLevel, level);
    ZSTD_CCtx_setParameter(c, ZSTD_c_checksumFlag, 1);
    std::vector<uint8_t> out(ZSTD_compressBound(raw.size()));
    size_t n = ZSTD_compress2(c, out.data(), out.size(), raw.data(), raw.size());
    ZSTD_freeCCtx(c);
    if (ZSTD_isError(n)) throw std::runtime_error(std::string("dcx: zstd ") + ZSTD_getErrorName(n));
    out.resize(n);
    return out;
}

inline uint32_t bytes_for_bits(uint64_t nbits) { return uint32_t((nbits + 7) / 8); }

// popcount of a presence byte slice.
inline uint32_t popcount_bytes(const uint8_t* p, size_t n) {
    uint32_t s = 0;
    for (size_t i = 0; i < n; ++i) s += uint32_t(__builtin_popcount(p[i]));
    return s;
}

void push_col(std::vector<ColumnDirEntry>& dir, std::vector<uint8_t>& payload,
              ColumnId id, ColEncoding enc, uint32_t elem_count,
              const std::vector<uint8_t>& col) {
    ColumnDirEntry e{};
    e.column_id = uint16_t(id);
    e.encoding = uint16_t(enc);
    e.elem_count = elem_count;
    e.byte_offset = uint32_t(payload.size());
    e.byte_size = uint32_t(col.size());
    e.reserved = 0;
    dir.push_back(e);
    payload.insert(payload.end(), col.begin(), col.end());
}

}  // namespace


namespace {
struct EncBlock { std::vector<uint8_t> frame; BlockHeader bh; BlockIndexEntry ie; };
EncBlock encode_block(const Batch& B, const DcxWriteOptions& opt, size_t lo, size_t hi,
                      const std::vector<size_t>& pres_off, const std::vector<size_t>& obs_off) {
    const uint32_t cf = opt.col_flags | B.col_flags;
    const bool has_ga = cf & COL_GAMMA_ANCIENT;
    const bool has_em = cf & COL_EM_MARGIN;
    const bool has_neff = cf & COL_EM_REF_NEFF;
    const bool has_gnull = cf & COL_GAMMA_NULL;
    const bool has_pread = cf & COL_P_READ;
    const bool dist3_derived = cf & COL_DIST3_DERIVED;
    const bool has_strand = cf & COL_STRAND;
    const bool has_obs = cf & COL_OBS_CODON;
    const bool has_best = cf & COL_BEST_CODON;
    const bool has_crank = cf & COL_CANDIDATE_RANK;
    const bool has_fident = cf & COL_FIDENT;

        const size_t nblk = hi - lo;
        const size_t b_lo = obs_off[lo], b_hi = obs_off[hi];
        const size_t nB = b_hi - b_lo;

        // ---- encode Grain A columns ----
        std::vector<uint8_t> qid, tord, tst, spn, gq, pq, emk, strd, gaq, emq, neff, gnullq, preadq, crank, fidq;
        uint64_t pq_prev = 0, tord_prev = 0; int64_t tst_prev = 0, qid_prev = 0;
        emk.assign(bytes_for_bits(nblk), 0);
        if (has_strand) strd.assign(bytes_for_bits(nblk), 0);
        for (size_t r = lo; r < hi; ++r) {
            put_varint(qid, zigzag(int64_t(B.query_ord[r]) - qid_prev)); qid_prev = int64_t(B.query_ord[r]);
            put_varint(tord, uint64_t(B.target_ord[r]) - tord_prev); tord_prev = B.target_ord[r];
            put_varint(tst, zigzag(int64_t(B.tstart[r]) - tst_prev)); tst_prev = B.tstart[r];
            put_varint(spn, B.span[r]);
            put_u16(gq, enc_unit(B.gamma[r]));
            put_u16(pq, enc_unit(B.posterior[r]));
            if (B.em_keep[r]) emk[(r - lo) >> 3] |= uint8_t(1u << ((r - lo) & 7));  // LSB-first (Arrow validity order)
            if (has_strand && B.strand[r]) strd[(r - lo) >> 3] |= uint8_t(1u << ((r - lo) & 7));
            if (has_ga) put_u16(gaq, enc_unit(B.gamma_ancient[r]));
            if (has_em) put_u16(emq, enc_unit(B.em_margin[r]));
            if (has_neff) { float f = B.em_ref_neff[r]; put_bytes(neff, &f, 4); }
            if (has_gnull) put_u16(gnullq, enc_unit(B.gamma_null[r]));
            if (has_pread) put_u16(preadq, enc_unit(B.p_read[r]));
            if (has_crank) crank.push_back(B.candidate_rank[r]);
            if (has_fident) put_u16(fidq, enc_unit(B.fident[r]));
            (void)pq_prev;
        }
        std::vector<uint8_t> pres(B.presence.begin() + pres_off[lo], B.presence.begin() + pres_off[hi]);

        // ---- encode Grain B columns ----
        std::vector<uint8_t> refaa, obsaa, codon, ocodon, bpq, pdq, d5, d3, dcl;
        refaa.assign(B.ref_aa.begin() + b_lo, B.ref_aa.begin() + b_hi);
        obsaa.assign(B.obs_aa.begin() + b_lo, B.obs_aa.begin() + b_hi);
        if (has_best) codon.assign(B.best_codon.begin() + b_lo, B.best_codon.begin() + b_hi);
        if (has_obs) ocodon.assign(B.obs_codon.begin() + b_lo, B.obs_codon.begin() + b_hi);
        dcl.assign(bytes_for_bits(nB * 2), 0);
        for (size_t j = b_lo; j < b_hi; ++j) {
            if (has_best) put_u16(bpq, enc_unit(B.best_post[j]));
            put_u16(pdq, enc_unit(B.p_damage[j]));
            put_u16(d5, B.dist5[j]);
            if (!dist3_derived) put_u16(d3, B.dist3[j]);
            uint8_t dc = B.damage_class[j] & 0x3;
            size_t bit = (j - b_lo) * 2;
            dcl[bit >> 3] |= uint8_t(uint32_t(dc) << (bit & 7));  // 2-bit LSB-first
        }

        // ---- assemble raw payload: RawBlockHeader + dir + columns ----
        std::vector<ColumnDirEntry> dir;
        std::vector<uint8_t> cols;
        push_col(dir, cols, ColumnId::QueryId, ColEncoding::DeltaVarint, uint32_t(nblk), qid);
        push_col(dir, cols, ColumnId::TargetId, ColEncoding::DeltaVarint, uint32_t(nblk), tord);
        push_col(dir, cols, ColumnId::TStart, ColEncoding::DeltaVarint, uint32_t(nblk), tst);
        push_col(dir, cols, ColumnId::Span, ColEncoding::Varint, uint32_t(nblk), spn);
        push_col(dir, cols, ColumnId::GammaQ, ColEncoding::U16Quant, uint32_t(nblk), gq);
        push_col(dir, cols, ColumnId::PosteriorQ, ColEncoding::U16Quant, uint32_t(nblk), pq);
        push_col(dir, cols, ColumnId::EmKeep, ColEncoding::Bitset1, uint32_t(nblk), emk);
        if (has_strand) push_col(dir, cols, ColumnId::Strand, ColEncoding::Bitset1, uint32_t(nblk), strd);
        if (has_ga) push_col(dir, cols, ColumnId::GammaAncientQ, ColEncoding::U16Quant, uint32_t(nblk), gaq);
        if (has_em) push_col(dir, cols, ColumnId::EmMarginQ, ColEncoding::U16Quant, uint32_t(nblk), emq);
        if (has_neff) push_col(dir, cols, ColumnId::EmRefNeff, ColEncoding::F32, uint32_t(nblk), neff);
        if (has_gnull) push_col(dir, cols, ColumnId::GammaNull, ColEncoding::U16Quant, uint32_t(nblk), gnullq);
        if (has_pread) push_col(dir, cols, ColumnId::PReadQ, ColEncoding::U16Quant, uint32_t(nblk), preadq);
        if (has_crank) push_col(dir, cols, ColumnId::CandidateRank, ColEncoding::Raw, uint32_t(nblk), crank);
        if (has_fident) push_col(dir, cols, ColumnId::FidentQ, ColEncoding::U16Quant, uint32_t(nblk), fidq);
        push_col(dir, cols, ColumnId::Presence, ColEncoding::Bitset1, uint32_t(nblk), pres);
        push_col(dir, cols, ColumnId::RefAa, ColEncoding::Raw, uint32_t(nB), refaa);
        push_col(dir, cols, ColumnId::ObsAa, ColEncoding::Raw, uint32_t(nB), obsaa);
        if (has_best) push_col(dir, cols, ColumnId::BestCodon, ColEncoding::CodonU8, uint32_t(nB), codon);
        if (has_obs) push_col(dir, cols, ColumnId::ObsCodon, ColEncoding::CodonU8, uint32_t(nB), ocodon);
        if (has_best) push_col(dir, cols, ColumnId::BestPostQ, ColEncoding::U16Quant, uint32_t(nB), bpq);
        push_col(dir, cols, ColumnId::PDamageQ, ColEncoding::U16Quant, uint32_t(nB), pdq);
        push_col(dir, cols, ColumnId::Dist5, ColEncoding::Raw, uint32_t(nB), d5);
        if (!dist3_derived) push_col(dir, cols, ColumnId::Dist3, ColEncoding::Raw, uint32_t(nB), d3);
        push_col(dir, cols, ColumnId::DamageClass2, ColEncoding::Bits2, uint32_t(nB), dcl);

        RawBlockHeader rbh{};
        rbh.n_records_a = uint32_t(nblk);
        rbh.n_records_b = uint32_t(nB);
        rbh.n_columns = uint32_t(dir.size());
        rbh.column_dir_bytes = uint32_t(dir.size() * sizeof(ColumnDirEntry));
        std::vector<uint8_t> raw;
        put_bytes(raw, &rbh, sizeof(rbh));
        put_bytes(raw, dir.data(), rbh.column_dir_bytes);
        raw.insert(raw.end(), cols.begin(), cols.end());

        std::vector<uint8_t> frame = zstd_frame(raw, opt.zstd_level);

        // ---- block header + trailer entry ----
        uint32_t min_ts = B.tstart[lo], max_te = 0;
        for (size_t r = lo; r < hi; ++r) max_te = std::max(max_te, B.tstart[r] + B.span[r] - 1);
        BlockHeader bh{};
        std::memcpy(bh.magic, BLOCK_MAGIC, 4);
        bh.version = VERSION_MAJOR; bh.flags = uint8_t(cf); bh.reserved = 0;
        bh.compressed_sz = uint32_t(frame.size());
        bh.raw_sz = uint32_t(raw.size());
        bh.n_records_a = uint32_t(nblk);
        bh.min_target_ord = B.target_ord[lo];
        bh.max_target_ord = B.target_ord[hi - 1];

        BlockIndexEntry ie{};
        ie.compressed_sz = bh.compressed_sz; ie.raw_sz = bh.raw_sz;
        ie.n_records_a = uint32_t(nblk); ie.n_records_b = uint32_t(nB);
        ie.min_target_ord = bh.min_target_ord; ie.max_target_ord = bh.max_target_ord;
        ie.min_tstart = min_ts; ie.max_tend = max_te;

        EncBlock eb;
        eb.frame = std::move(frame);
        eb.bh = bh;
        eb.ie = ie;   // ie.offset filled during the ordered flush below
        return eb;
}

bool finish_dcx(std::ofstream& out, size_t off, const std::vector<BlockIndexEntry>& index,
                const std::vector<std::string>& target_dict, const std::vector<std::string>& query_dict,
                const DcxWriteOptions& opt, uint32_t cf, uint64_t nA, uint64_t nB) {
    // ---- trailer index + dictionaries + footer (built in a small `tail` buffer — its
    // size is O(n_blocks + dict) not O(records) — then streamed once) ----
    std::vector<uint8_t> file;   // OPT-H: now holds ONLY the trailer, not the blocks
    const uint64_t index_offset = off;
    IndexHeader ih{};
    std::memcpy(ih.magic, INDEX_MAGIC, 4);
    ih.version = VERSION_MAJOR;
    ih.n_blocks = uint32_t(index.size());
    const uint64_t entries_end = index_offset + sizeof(IndexHeader) + index.size() * sizeof(BlockIndexEntry);
    // dictionaries land right after the entries; record their offsets up front.
    auto dict_bytes = [](const std::vector<std::string>& d) {
        size_t n = 4; for (auto& s : d) n += 4 + s.size(); return n;
    };
    ih.tgt_dict_offset = entries_end;
    ih.tgt_dict_bytes = dict_bytes(target_dict);
    uint64_t after_dicts = entries_end + ih.tgt_dict_bytes;
    if (cf & COL_QID_DICT) {
        ih.qid_dict_offset = after_dicts;
        ih.qid_dict_bytes = dict_bytes(query_dict);
        after_dicts += ih.qid_dict_bytes;
    }
    ih.calib_offset = after_dicts;
    ih.calib_bytes = 2 + opt.calib_model.size() + 2 + opt.library_label.size()
                     + 4 + opt.calib_contract_blob.size();

    std::vector<uint8_t> idx_bytes;
    put_bytes(idx_bytes, &ih, sizeof(ih));
    put_bytes(idx_bytes, index.data(), index.size() * sizeof(BlockIndexEntry));
    const uint32_t idx_adler = adler32(idx_bytes.data(), idx_bytes.size());
    file.insert(file.end(), idx_bytes.begin(), idx_bytes.end());

    auto put_dict = [&](const std::vector<std::string>& d) {
        put_u32(file, uint32_t(d.size()));
        for (auto& s : d) { put_u32(file, uint32_t(s.size())); put_bytes(file, s.data(), s.size()); }
    };
    put_dict(target_dict);
    if (cf & COL_QID_DICT) put_dict(query_dict);

    // calib section: u16 model_len | model | u16 label_len | label | u32 blob_len | blob
    put_u16(file, uint16_t(opt.calib_model.size()));
    put_bytes(file, opt.calib_model.data(), opt.calib_model.size());
    put_u16(file, uint16_t(opt.library_label.size()));
    put_bytes(file, opt.library_label.data(), opt.library_label.size());
    put_u32(file, uint32_t(opt.calib_contract_blob.size()));
    put_bytes(file, opt.calib_contract_blob.data(), opt.calib_contract_blob.size());

    Footer ft{};
    ft.index_offset = index_offset;
    ft.index_adler32 = idx_adler;
    ft.reserved = 0;
    std::memcpy(ft.magic, FOOTER_MAGIC, 4);
    ft.footer_bytes = uint32_t(sizeof(Footer));
    put_bytes(file, &ft, sizeof(ft));

    // ---- patch FileHeader ----
    FileHeader fh{};
    std::memcpy(fh.magic, FILE_MAGIC, 8);
    fh.version_major = VERSION_MAJOR;
    // v2.2 (soft multi-mapping) only when the capability is actually emitted; a non-multimap
    // write stays v2.1 so it is byte-identical to the pre-2.2 writer.
    fh.version_minor = (cf & COL_SOFT_MULTIMAP) ? uint8_t(2) : VERSION_MINOR;
    fh.header_bytes = uint16_t(sizeof(FileHeader));
    fh.col_flags = cf; fh.codec = CODEC_ZSTD;
    fh.sample_flags = opt.sample_flags;
    std::memcpy(fh.calib_contract_id, opt.calib_contract_id, 16);
    fh.n_blocks = uint32_t(index.size());
    fh.n_records_a = nA; fh.n_records_b = nB;
    fh.first_block_offset = sizeof(FileHeader);

    // OPT-H: stream the trailer, then seekp(0) and overwrite the placeholder FileHeader
    // with the real one (identical bytes to the old in-place file.data() patch).
    out.write(reinterpret_cast<const char*>(file.data()), std::streamsize(file.size()));
    off += file.size();
    (void)off;
    out.seekp(0);
    out.write(reinterpret_cast<const char*>(&fh), sizeof(fh));
    return bool(out);
}

} // namespace

bool write_dcx(const char* path, const Batch& B, const DcxWriteOptions& opt) {
    const size_t nA = B.target_ord.size();
    const uint32_t cf = opt.col_flags | B.col_flags;
    const bool has_ga = cf & COL_GAMMA_ANCIENT;
    const bool has_em = cf & COL_EM_MARGIN;
    const bool has_neff = cf & COL_EM_REF_NEFF;
    const bool has_gnull = cf & COL_GAMMA_NULL;
    const bool has_pread = cf & COL_P_READ;
    const bool dist3_derived = cf & COL_DIST3_DERIVED;
    const bool has_strand = cf & COL_STRAND;
    const bool has_obs = cf & COL_OBS_CODON;
    const bool has_best = cf & COL_BEST_CODON;
    const bool has_crank = cf & COL_CANDIDATE_RANK;
    const bool has_fident = cf & COL_FIDENT;

    // Per-record presence byte offset and popcount (also the Grain B offset).
    std::vector<size_t> pres_off(nA + 1, 0), obs_off(nA + 1, 0);
    for (size_t r = 0; r < nA; ++r) {
        uint32_t pb = bytes_for_bits(B.span[r]);
        pres_off[r + 1] = pres_off[r] + pb;
        obs_off[r + 1] = obs_off[r] + popcount_bytes(&B.presence[pres_off[r]], pb);
    }
    // ref_aa is the mandatory Grain-B count reference (v2: residues, not codons).
    if (obs_off[nA] != B.ref_aa.size())
        throw std::runtime_error("dcx: Σpopcount(presence) != ref_aa size — desynced observation grain");
    if (B.obs_aa.size() != B.ref_aa.size())
        throw std::runtime_error("dcx: obs_aa size != ref_aa size — desynced observation grain");
    if (has_best && (B.best_codon.size() != B.ref_aa.size() || B.best_post.size() != B.ref_aa.size()))
        throw std::runtime_error("dcx: COL_BEST_CODON set but best_codon/best_post size != ref_aa size");
    if (has_obs && B.obs_codon.size() != B.ref_aa.size())
        throw std::runtime_error("dcx: obs_codon size != ref_aa size — desynced observation grain");
    if (has_strand && B.strand.size() != nA)
        throw std::runtime_error("dcx: strand size != n_records_a");
    if (has_pread && B.p_read.size() != nA)
        throw std::runtime_error("dcx: p_read size != n_records_a");
    if (has_crank && B.candidate_rank.size() != nA)
        throw std::runtime_error("dcx: candidate_rank size != n_records_a");
    bool calibed = false;
    for (int i = 0; i < 16; ++i) calibed |= (opt.calib_contract_id[i] != 0);
    if (!calibed)
        throw std::runtime_error("dcx: p_damage emitted without calib_contract_id — uncalibrated write");

    // OPT-H: stream blocks straight to disk instead of buffering the whole file in RAM
    // (peak was O(total compressed output)). The byte layout is IDENTICAL — same header,
    // same blocks, same trailer, in the same order — only WHEN bytes are written changes.
    // The FileHeader is written twice: a zero placeholder now, then the real one via
    // seekp(0) at the end, once n_blocks / n_records are known (as the old code patched
    // file.data() in place).
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    size_t off = 0;
    { const std::vector<uint8_t> hdr0(sizeof(FileHeader), 0);
      out.write(reinterpret_cast<const char*>(hdr0.data()), std::streamsize(hdr0.size()));
      off += hdr0.size(); }
    std::vector<BlockIndexEntry> index;

    // encode_block: pure per-block encode + zstd compress. Depends only on the block's
    // own records (delta encoders reset per block) + read-only offsets, so blocks are
    // independent and can be compressed in parallel, then flushed in order (below).

    // Greedy block cut (serial, deterministic — the only cross-block dependency): accumulate
    // records until the raw estimate exceeds the cap, then cut at the next target boundary.
    // A single oversized target is its own block. Identical cut to the pre-parallel writer.
    std::vector<std::pair<size_t, size_t>> bounds;
    { size_t lo = 0;
      while (lo < nA) {
          size_t hi = lo, est = 0;
          while (hi < nA) {
              uint32_t pb = bytes_for_bits(B.span[hi]);
              size_t nobs = obs_off[hi + 1] - obs_off[hi];
              est += 12 + pb + nobs * 8;  // ~A scalars + presence + B scalars
              ++hi;
              bool at_boundary = (hi == nA) || (B.target_ord[hi] != B.target_ord[hi - 1]);
              if (est >= opt.block_raw_cap && at_boundary) break;
          }
          bounds.emplace_back(lo, hi);
          lo = hi;
      } }

    // Parallel encode+compress, flushed strictly in block order through a bounded window of
    // K blocks. Byte-identical to the serial writer (blocks independent, zstd deterministic,
    // frames written in order); peak stays ~K blocks so OPT-H's memory cap is preserved.
#ifdef _OPENMP
    const int K = std::max(1, omp_get_max_threads());
#else
    const int K = 1;
#endif
    const size_t win_sz = std::min(size_t(K), bounds.size() ? bounds.size() : size_t(1));
    std::vector<EncBlock> win(win_sz);
    for (size_t c = 0; c < bounds.size(); c += size_t(K)) {
        const size_t ce = std::min(c + size_t(K), bounds.size());
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t bi = c; bi < ce; ++bi)
            win[bi - c] = encode_block(B, opt, bounds[bi].first, bounds[bi].second, pres_off, obs_off);
        for (size_t bi = c; bi < ce; ++bi) {
            EncBlock& eb = win[bi - c];
            eb.ie.offset = off;
            index.push_back(eb.ie);
            out.write(reinterpret_cast<const char*>(&eb.bh), sizeof(eb.bh));
            out.write(reinterpret_cast<const char*>(eb.frame.data()), std::streamsize(eb.frame.size()));
            off += sizeof(eb.bh) + eb.frame.size();
            std::vector<uint8_t>().swap(eb.frame);
        }
    }

    return finish_dcx(out, off, index, B.target_dict, B.query_dict, opt, cf, nA, B.ref_aa.size());
}

struct StreamWriter::Impl {
    DcxWriteOptions opt;
    bool auth_only, closed = false;
    std::ofstream out;
    Batch current;
    size_t estimate = 0, offset = sizeof(FileHeader), workers = 1;
    uint64_t nA = 0, nB = 0;
    std::vector<BlockIndexEntry> index;
    std::deque<std::future<EncBlock>> pending;

    Impl(const char* path, const DcxWriteOptions& options, bool auth)
        : opt(options), auth_only(auth), out(path, std::ios::binary | std::ios::trunc) {
        bool calibrated = false;
        for (uint8_t b : opt.calib_contract_id) calibrated |= b != 0;
        if (!calibrated) throw std::runtime_error("dcx: missing calib_contract_id");
        if (!out) throw std::runtime_error(std::string("dcx: cannot open ") + path);
        FileHeader placeholder{};
        out.write(reinterpret_cast<const char*>(&placeholder), sizeof(placeholder));
#ifdef _OPENMP
        workers = static_cast<size_t>(std::max(1, omp_get_max_threads()));
#endif
        current.col_flags = opt.col_flags;
    }

    void flush_one() {
        EncBlock block = pending.front().get();
        pending.pop_front();
        block.ie.offset = offset;
        index.push_back(block.ie);
        out.write(reinterpret_cast<const char*>(&block.bh), sizeof(block.bh));
        out.write(reinterpret_cast<const char*>(block.frame.data()), block.frame.size());
        offset += sizeof(block.bh) + block.frame.size();
        if (!out) throw std::runtime_error("dcx: block write failed");
    }

    void submit() {
        if (current.target_ord.empty()) return;
        if (pending.size() >= workers) flush_one();
        pending.emplace_back(std::async(std::launch::async,
            [batch = std::move(current), options = opt]() {
                const size_t n = batch.target_ord.size();
                std::vector<size_t> presence(n + 1), observations(n + 1);
                for (size_t i = 0; i < n; ++i) {
                    const size_t bytes = bytes_for_bits(batch.span[i]);
                    presence[i + 1] = presence[i] + bytes;
                    observations[i + 1] = observations[i] +
                        (bytes ? popcount_bytes(batch.presence.data() + presence[i], bytes) : 0);
                }
                if (observations.back() != batch.ref_aa.size())
                    throw std::runtime_error("dcx: presence/observation count mismatch");
                return encode_block(batch, options, 0, n, presence, observations);
            }));
        current = Batch{};
        current.col_flags = opt.col_flags;
        estimate = 0;
    }
};

StreamWriter::StreamWriter(const char* path, const DcxWriteOptions& opt, bool auth_only)
    : impl_(std::make_unique<Impl>(path, opt, auth_only)) {}
StreamWriter::~StreamWriter() = default;

void StreamWriter::append(const Record& rec, uint32_t target, uint64_t query, float posterior) {
    Impl& w = *impl_;
    if (w.closed) throw std::runtime_error("dcx: append after close");
    if (!w.current.target_ord.empty() && w.current.target_ord.back() != target &&
        w.estimate >= w.opt.block_raw_cap) w.submit();
    Batch& b = w.current;
    const uint32_t cf = w.opt.col_flags;
    b.query_ord.push_back(query); b.target_ord.push_back(target);
    b.tstart.push_back(rec.tstart); b.span.push_back(w.auth_only ? 0 : rec.span);
    b.gamma.push_back(rec.gamma); b.posterior.push_back(posterior);
    b.em_keep.push_back(rec.em_keep);
    if (cf & COL_STRAND) b.strand.push_back(rec.strand);
    if (cf & COL_GAMMA_ANCIENT) b.gamma_ancient.push_back(rec.gamma_ancient);
    if (cf & COL_EM_MARGIN) b.em_margin.push_back(rec.em_margin);
    if (cf & COL_EM_REF_NEFF) b.em_ref_neff.push_back(rec.em_ref_neff);
    if (cf & COL_GAMMA_NULL) b.gamma_null.push_back(rec.gamma_null);
    if (cf & COL_P_READ) b.p_read.push_back(rec.p_read);
    if (cf & COL_CANDIDATE_RANK) b.candidate_rank.push_back(rec.candidate_rank);
    if (cf & COL_FIDENT) b.fident.push_back(rec.fident);
    size_t observations = 0, bytes = 0;
    if (!w.auth_only) {
        bytes = bytes_for_bits(rec.span);
        observations = rec.ref_aa.size();
        if (rec.presence.size() != bytes || rec.obs_aa.size() != observations ||
            rec.p_damage.size() != observations || rec.dist5.size() != observations ||
            rec.dcls.size() != observations ||
            (!(cf & COL_DIST3_DERIVED) && rec.dist3.size() != observations) ||
            ((cf & COL_OBS_CODON) && rec.obs_codon.size() != observations) ||
            ((cf & COL_BEST_CODON) && (rec.best_codon.size() != observations || rec.best_post.size() != observations)))
            throw std::runtime_error("dcx: record column length mismatch");
        auto append = [](auto& to, const auto& from) { to.insert(to.end(), from.begin(), from.end()); };
        append(b.presence, rec.presence); append(b.ref_aa, rec.ref_aa); append(b.obs_aa, rec.obs_aa);
        append(b.damage_class, rec.dcls); append(b.p_damage, rec.p_damage); append(b.dist5, rec.dist5);
        if (!(cf & COL_DIST3_DERIVED)) append(b.dist3, rec.dist3);
        if (cf & COL_OBS_CODON) append(b.obs_codon, rec.obs_codon);
        if (cf & COL_BEST_CODON) { append(b.best_codon, rec.best_codon); append(b.best_post, rec.best_post); }
    }
    w.estimate += 12 + bytes + observations * 8;
    ++w.nA; w.nB += observations;
}

void StreamWriter::close(const std::vector<std::string>& targets, const std::vector<std::string>& queries) {
    Impl& w = *impl_;
    if (w.closed) throw std::runtime_error("dcx: already closed");
    w.submit();
    while (!w.pending.empty()) w.flush_one();
    if (!finish_dcx(w.out, w.offset, w.index, targets, queries, w.opt, w.opt.col_flags, w.nA, w.nB))
        throw std::runtime_error("dcx: trailer write failed");
    w.out.close();
    if (!w.out) throw std::runtime_error("dcx: close failed");
    w.closed = true;
}

// write_dcx_auth — the GRAIN-A AUTH SIDECAR (fold-gap fast path; fusion dual-f6f6aa14).
// A tiny companion `<lib>.dcx.auth` holding ONLY the grain-A scalars protein_auth needs
// (query_id, target_id, gamma, em_keep, p_read; the v1 alignment ordinal is the row order).
// Implemented as a grain-A-ONLY .dcx (n_records_b = 0): we reuse write_dcx verbatim with
// grain_b stripped, so the grain-A scalar VALUES are byte-identical to the full .dcx and the
// SAME reader/reduce (dcx_rs.reduce_grain_a) consumes it — zero new format, zero drift. The
// sidecar drops the entire grain-B position payload (the ~200GB bulk), so auth reads a small
// fraction. (DCX2 follow-up: also drop tstart/span/posterior for the exact grain_a.parquet
// projection; kept here for maximal reuse + guaranteed value-identity.)
bool write_dcx_auth(const char* path, const Batch& B, const DcxWriteOptions& opt) {
    Batch A = B;                       // grain-A scalars + dicts + col_flags carried verbatim
    A.span.assign(B.span.size(), 0);   // span=0 -> zero presence bytes, zero grain-B records
    A.presence.clear();
    A.ref_aa.clear();  A.obs_aa.clear();
    A.best_codon.clear();  A.obs_codon.clear();  A.best_post.clear();
    A.p_damage.clear();  A.dist5.clear();  A.dist3.clear();  A.damage_class.clear();
    return write_dcx(path, A, opt);    // Σpopcount(empty presence)=0=ref_aa.size() -> valid
}

// write_auth_stats — WRITE-TIME AUTH REDUCE output (fusion dual-66c04a26 #1). A tiny
// versioned, zstd'd per-target table {N, logLR, k@tau-grid} the annotate accumulates in its
// read-grouped pass (BEFORE the target-sort), so convergence auth reads bytes/target instead
// of re-decoding grain_a. Columnar payload: target dict + N[u64] + k[n_tau][u64] + logLR[f64].
// Header (uncompressed): magic|version|n_targets|n_tau|tau[n_tau]|lib_label|raw_sz|frame_sz.
bool write_auth_stats(const char* path, const std::string& lib_label,
                      const double* tau, uint32_t n_tau,
                      const std::vector<std::string>& targets,
                      const std::vector<uint64_t>& N,
                      const std::vector<uint64_t>& k_flat,   // n_targets*n_tau, row-major
                      const std::vector<double>& logLR, int zstd_level) {
    const uint64_t nt = targets.size();
    if (N.size() != nt || logLR.size() != nt || k_flat.size() != nt * n_tau)
        throw std::runtime_error("authstats: column length mismatch");

    std::vector<uint8_t> pay;                // columnar payload (zstd'd)
    put_u64(pay, nt);
    for (const auto& s : targets) { put_u32(pay, uint32_t(s.size())); put_bytes(pay, s.data(), s.size()); }
    for (uint64_t i = 0; i < nt; ++i) put_u64(pay, N[i]);
    for (uint32_t j = 0; j < n_tau; ++j)
        for (uint64_t i = 0; i < nt; ++i) put_u64(pay, k_flat[i * n_tau + j]);
    for (uint64_t i = 0; i < nt; ++i) { uint64_t b; std::memcpy(&b, &logLR[i], 8); put_u64(pay, b); }
    std::vector<uint8_t> frame = zstd_frame(pay, zstd_level);

    std::vector<uint8_t> file;               // header (uncompressed) + frame
    const char magic[8] = {'D','C','X','A','U','T','H','S'};
    put_bytes(file, magic, 8);
    put_u32(file, 1u);                       // format version
    put_u64(file, nt);
    put_u32(file, n_tau);
    for (uint32_t j = 0; j < n_tau; ++j) { uint64_t b; std::memcpy(&b, &tau[j], 8); put_u64(file, b); }
    put_u32(file, uint32_t(lib_label.size()));
    put_bytes(file, lib_label.data(), lib_label.size());
    put_u64(file, uint64_t(pay.size()));
    put_u64(file, uint64_t(frame.size()));
    file.insert(file.end(), frame.begin(), frame.end());

    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    out.write(reinterpret_cast<const char*>(file.data()), std::streamsize(file.size()));
    return bool(out);
}

}  // namespace dcx
