// dart_columnar.cpp — .dcx writer (write_dcx). Depends only on zstd.
// On-disk contract is dart_columnar.hpp. See that header for the byte layout;
// this file encodes a Batch into size-bounded, zstd-compressed columnar blocks
// with a trailer zone-map index, string dictionaries, and an adler32-checked
// footer. The consumer (Rust reader in convergence) mirrors read_block().

#include "dart/dart_columnar.hpp"

#include <zstd.h>

#include <cstring>
#include <fstream>
#include <stdexcept>

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
    bool calibed = false;
    for (int i = 0; i < 16; ++i) calibed |= (opt.calib_contract_id[i] != 0);
    if (!calibed)
        throw std::runtime_error("dcx: p_damage emitted without calib_contract_id — uncalibrated write");

    std::vector<uint8_t> file;
    file.resize(sizeof(FileHeader));  // FileHeader patched at the end
    std::vector<BlockIndexEntry> index;

    // Greedy block cut: accumulate records until the raw estimate exceeds the
    // cap, then cut at the next target boundary. A single oversized target is
    // its own block.
    size_t lo = 0;
    while (lo < nA) {
        size_t hi = lo;
        size_t est = 0;
        while (hi < nA) {
            uint32_t pb = bytes_for_bits(B.span[hi]);
            size_t nobs = obs_off[hi + 1] - obs_off[hi];
            est += 12 + pb + nobs * 8;  // ~A scalars + presence + B scalars
            ++hi;
            bool at_boundary = (hi == nA) || (B.target_ord[hi] != B.target_ord[hi - 1]);
            if (est >= opt.block_raw_cap && at_boundary) break;
        }
        const size_t nblk = hi - lo;
        const size_t b_lo = obs_off[lo], b_hi = obs_off[hi];
        const size_t nB = b_hi - b_lo;

        // ---- encode Grain A columns ----
        std::vector<uint8_t> qid, tord, tst, spn, gq, pq, emk, strd, gaq, emq, neff, gnullq, preadq;
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
        ie.offset = file.size();
        ie.compressed_sz = bh.compressed_sz; ie.raw_sz = bh.raw_sz;
        ie.n_records_a = uint32_t(nblk); ie.n_records_b = uint32_t(nB);
        ie.min_target_ord = bh.min_target_ord; ie.max_target_ord = bh.max_target_ord;
        ie.min_tstart = min_ts; ie.max_tend = max_te;
        index.push_back(ie);

        put_bytes(file, &bh, sizeof(bh));
        file.insert(file.end(), frame.begin(), frame.end());
        lo = hi;
    }

    // ---- trailer index + dictionaries + footer ----
    const uint64_t index_offset = file.size();
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
    ih.tgt_dict_bytes = dict_bytes(B.target_dict);
    uint64_t after_dicts = entries_end + ih.tgt_dict_bytes;
    if (cf & COL_QID_DICT) {
        ih.qid_dict_offset = after_dicts;
        ih.qid_dict_bytes = dict_bytes(B.query_dict);
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
    put_dict(B.target_dict);
    if (cf & COL_QID_DICT) put_dict(B.query_dict);

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
    fh.version_major = VERSION_MAJOR; fh.version_minor = VERSION_MINOR;
    fh.header_bytes = uint16_t(sizeof(FileHeader));
    fh.col_flags = cf; fh.codec = CODEC_ZSTD;
    fh.sample_flags = opt.sample_flags;
    std::memcpy(fh.calib_contract_id, opt.calib_contract_id, 16);
    fh.n_blocks = uint32_t(index.size());
    fh.n_records_a = nA; fh.n_records_b = B.ref_aa.size();
    fh.first_block_offset = sizeof(FileHeader);
    std::memcpy(file.data(), &fh, sizeof(fh));

    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    out.write(reinterpret_cast<const char*>(file.data()), std::streamsize(file.size()));
    return bool(out);
}

}  // namespace dcx
