// Round-trip test for the .dcx writer. Exercises the offset math that silently
// corrupts everything downstream if wrong: FileHeader (grown to 96B), footer
// adler32 over the index, the calib section offset, and the obs_codon grain.
// Compile: g++ -std=c++17 -I../include -o test_dcx_roundtrip test_dcx_roundtrip.cpp ../src/dart_columnar.cpp $(pkg-config --cflags --libs libzstd)

#include "dart/dart_columnar.hpp"

#include <zstd.h>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

using namespace dcx;

static std::vector<uint8_t> slurp(const char* p) {
    std::ifstream f(p, std::ios::binary);
    return std::vector<uint8_t>((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}
static uint32_t adler32_ref(const uint8_t* p, size_t n) {
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < n; ++i) { a = (a + p[i]) % 65521; b = (b + a) % 65521; }
    return (b << 16) | a;
}
template <class T> static T read_at(const std::vector<uint8_t>& f, size_t off) {
    T v; std::memcpy(&v, f.data() + off, sizeof(T)); return v;
}

static Batch make_batch() {
    Batch B;
    // 2 alignments, same target, spans [0,2) and [3,5).
    B.query_ord = {0, 1};
    B.target_ord = {0, 0};
    B.tstart = {0, 3};
    B.span = {2, 2};
    B.gamma = {0.5f, 0.6f};
    B.posterior = {0.9f, 0.95f};
    B.em_keep = {1, 1};
    B.strand = {0, 1};                 // COL_STRAND
    B.presence = {0x03, 0x01};         // rec0: bits 0,1 (2 obs); rec1: bit 0 (1 obs) => nB=3
    B.ref_aa = {'A', 'R', 'N'};        // MANDATORY residues (ASCII)
    B.obs_aa = {'A', 'K', 'N'};        // pos 1: R->K substitution
    B.best_codon = {10, 20, 30};       // MAP — present iff COL_BEST_CODON
    B.obs_codon = {11, 21, 31};        // COL_OBS_CODON — raw
    B.best_post = {0.9f, 0.8f, 0.7f};
    B.p_damage = {0.1f, 0.2f, 0.3f};
    B.dist5 = {1, 2, 3};
    B.dist3 = {5, 6, 7};
    B.damage_class = {1, 2, 0};
    B.target_dict = {"tgtA"};
    return B;
}

static DcxWriteOptions make_opts() {
    DcxWriteOptions o;
    o.col_flags = COL_OBS_CODON | COL_STRAND | COL_BEST_CODON;
    o.sample_flags = SAMPLE_LIB_DS;    // 0-based coords (COORD_1BASED unset)
    for (int i = 0; i < 16; ++i) o.calib_contract_id[i] = uint8_t(i + 1);
    o.calib_model = "dfit/1.4";
    o.library_label = "KapK-L03";
    o.calib_contract_blob = "lambda=0.3;d_max=30;code=abc123";
    return o;
}

int main() {
    const char* path = "/tmp/dcx_roundtrip.dcx";
    Batch B = make_batch();
    DcxWriteOptions opt = make_opts();
    assert(write_dcx(path, B, opt));

    std::vector<uint8_t> f = slurp(path);
    FileHeader fh = read_at<FileHeader>(f, 0);
    assert(std::memcmp(fh.magic, FILE_MAGIC, 8) == 0);
    assert(fh.header_bytes == sizeof(FileHeader));
    assert(fh.header_bytes == 96);
    assert(fh.version_major == VERSION_MAJOR);
    assert(fh.n_records_a == 2);
    assert(fh.n_records_b == 3);          // == Σ popcount(presence)
    assert(fh.sample_flags == SAMPLE_LIB_DS);
    assert(fh.first_block_offset == sizeof(FileHeader));
    for (int i = 0; i < 16; ++i) assert(fh.calib_contract_id[i] == uint8_t(i + 1));

    // Footer at EOF-24, then adler32 over IndexHeader+entries.
    Footer ft = read_at<Footer>(f, f.size() - sizeof(Footer));
    assert(std::memcmp(ft.magic, FOOTER_MAGIC, 4) == 0);
    assert(ft.footer_bytes == sizeof(Footer));
    IndexHeader ih = read_at<IndexHeader>(f, ft.index_offset);
    assert(std::memcmp(ih.magic, INDEX_MAGIC, 4) == 0);
    assert(ih.n_blocks == fh.n_blocks);
    size_t idx_span = sizeof(IndexHeader) + size_t(ih.n_blocks) * sizeof(BlockIndexEntry);
    assert(adler32_ref(f.data() + ft.index_offset, idx_span) == ft.index_adler32);

    // calib section round-trips.
    size_t co = ih.calib_offset;
    uint16_t mlen = read_at<uint16_t>(f, co); co += 2;
    assert(std::string(reinterpret_cast<const char*>(f.data() + co), mlen) == "dfit/1.4"); co += mlen;
    uint16_t llen = read_at<uint16_t>(f, co); co += 2;
    assert(std::string(reinterpret_cast<const char*>(f.data() + co), llen) == "KapK-L03"); co += llen;
    uint32_t blen = read_at<uint32_t>(f, co); co += 4;
    assert(std::string(reinterpret_cast<const char*>(f.data() + co), blen) == opt.calib_contract_blob);
    assert(co + blen == ih.calib_offset + ih.calib_bytes);

    // Decode block 0 and verify the raw ObsCodon column survived.
    BlockIndexEntry ie = read_at<BlockIndexEntry>(f, ft.index_offset + sizeof(IndexHeader));
    BlockHeader bh = read_at<BlockHeader>(f, ie.offset);
    assert(std::memcmp(bh.magic, BLOCK_MAGIC, 4) == 0);
    assert(ie.n_records_b == 3);
    std::vector<uint8_t> raw(bh.raw_sz);
    size_t frame_off = ie.offset + sizeof(BlockHeader);
    size_t got = ZSTD_decompress(raw.data(), raw.size(), f.data() + frame_off, bh.compressed_sz);
    assert(!ZSTD_isError(got) && got == bh.raw_sz);
    RawBlockHeader rbh = read_at<RawBlockHeader>(raw, 0);
    assert(rbh.n_records_a == 2 && rbh.n_records_b == 3);
    size_t dir_off = sizeof(RawBlockHeader);
    size_t cols_off = dir_off + rbh.column_dir_bytes;
    bool saw_obs = false, saw_strand = false, saw_best = false, saw_refaa = false, saw_obsaa = false;
    for (uint32_t c = 0; c < rbh.n_columns; ++c) {
        ColumnDirEntry e = read_at<ColumnDirEntry>(raw, dir_off + c * sizeof(ColumnDirEntry));
        const uint8_t* col = raw.data() + cols_off + e.byte_offset;
        if (e.column_id == uint16_t(ColumnId::ObsCodon)) {
            saw_obs = true;
            assert(e.elem_count == 3 && e.byte_size == 3);
            assert(col[0] == 11 && col[1] == 21 && col[2] == 31);
        }
        if (e.column_id == uint16_t(ColumnId::BestCodon)) {
            saw_best = true;
            assert(e.elem_count == 3 && col[0] == 10 && col[1] == 20 && col[2] == 30);
        }
        if (e.column_id == uint16_t(ColumnId::RefAa)) {
            saw_refaa = true;
            assert(e.elem_count == 3 && col[0] == 'A' && col[1] == 'R' && col[2] == 'N');
        }
        if (e.column_id == uint16_t(ColumnId::ObsAa)) {
            saw_obsaa = true;
            assert(e.elem_count == 3 && col[0] == 'A' && col[1] == 'K' && col[2] == 'N');
        }
        if (e.column_id == uint16_t(ColumnId::Strand)) saw_strand = true;
    }
    assert(saw_obs && saw_strand && saw_best && saw_refaa && saw_obsaa);

    // With COL_BEST_CODON unset, RefAa/ObsAa still carry residues and no codon
    // columns are emitted — the v2 codon-free write path.
    {
        Batch B2 = B; B2.best_codon.clear(); B2.best_post.clear();
        DcxWriteOptions o2 = opt; o2.col_flags = COL_STRAND;   // no OBS, no BEST
        assert(write_dcx("/tmp/dcx_nocodon.dcx", B2, o2));
        std::vector<uint8_t> g = slurp("/tmp/dcx_nocodon.dcx");
        BlockIndexEntry ie2 = read_at<BlockIndexEntry>(g, read_at<Footer>(g, g.size() - sizeof(Footer)).index_offset + sizeof(IndexHeader));
        BlockHeader bh2 = read_at<BlockHeader>(g, ie2.offset);
        std::vector<uint8_t> raw2(bh2.raw_sz);
        assert(!ZSTD_isError(ZSTD_decompress(raw2.data(), raw2.size(), g.data() + ie2.offset + sizeof(BlockHeader), bh2.compressed_sz)));
        RawBlockHeader rbh2 = read_at<RawBlockHeader>(raw2, 0);
        bool has_codon = false, has_refaa = false;
        for (uint32_t c = 0; c < rbh2.n_columns; ++c) {
            ColumnDirEntry e = read_at<ColumnDirEntry>(raw2, sizeof(RawBlockHeader) + c * sizeof(ColumnDirEntry));
            if (e.column_id == uint16_t(ColumnId::BestCodon) || e.column_id == uint16_t(ColumnId::ObsCodon)) has_codon = true;
            if (e.column_id == uint16_t(ColumnId::RefAa)) has_refaa = true;
        }
        assert(!has_codon && has_refaa);
        std::remove("/tmp/dcx_nocodon.dcx");
    }

    // Grain-A auth SIDECAR (write_dcx_auth): grain-B stripped (n_records_b==0), grain-A
    // auth columns intact. Written to a stable path for the cross-language E2E (dcx_rs
    // reads it and computes protein_auth).
    {
        Batch Ba = B;
        Ba.p_read = {0.9f, 0.3f};                 // COL_P_READ auth column
        DcxWriteOptions oa = opt;
        oa.col_flags = opt.col_flags | COL_P_READ;
        const char* apath = "/projects/caeg/scratch/kbd606/tmp/convergence_v63_rebuild/dcx_auth_test.dcx.auth";
        assert(write_dcx_auth(apath, Ba, oa));
        std::vector<uint8_t> a = slurp(apath);
        FileHeader fha = read_at<FileHeader>(a, 0);
        assert(fha.n_records_a == 2);
        assert(fha.n_records_b == 0);             // grain-B stripped -> tiny sidecar
        BlockIndexEntry iea = read_at<BlockIndexEntry>(
            a, read_at<Footer>(a, a.size() - sizeof(Footer)).index_offset + sizeof(IndexHeader));
        BlockHeader bha = read_at<BlockHeader>(a, iea.offset);
        std::vector<uint8_t> rawa(bha.raw_sz);
        assert(!ZSTD_isError(ZSTD_decompress(rawa.data(), rawa.size(),
               a.data() + iea.offset + sizeof(BlockHeader), bha.compressed_sz)));
        RawBlockHeader rbha = read_at<RawBlockHeader>(rawa, 0);
        assert(rbha.n_records_a == 2 && rbha.n_records_b == 0);
        bool saw_pread = false, saw_gamma = false, saw_keep = false, saw_qid = false, saw_tid = false;
        for (uint32_t c = 0; c < rbha.n_columns; ++c) {
            ColumnDirEntry e = read_at<ColumnDirEntry>(rawa, sizeof(RawBlockHeader) + c * sizeof(ColumnDirEntry));
            if (e.column_id == uint16_t(ColumnId::PReadQ))  { saw_pread = true; assert(e.elem_count == 2); }
            if (e.column_id == uint16_t(ColumnId::GammaQ))  { saw_gamma = true; assert(e.elem_count == 2); }
            if (e.column_id == uint16_t(ColumnId::EmKeep))    saw_keep = true;
            if (e.column_id == uint16_t(ColumnId::QueryId))   saw_qid = true;
            if (e.column_id == uint16_t(ColumnId::TargetId))  saw_tid = true;
            if (e.column_id == uint16_t(ColumnId::RefAa))    assert(e.elem_count == 0);  // grain-B stripped
        }
        assert(saw_pread && saw_gamma && saw_keep && saw_qid && saw_tid);
        std::puts("DCX_AUTH_SIDECAR_OK");
    }

    // Write-time AUTHSTATS format round-trip (write_auth_stats -> Rust read_auth_stats E2E).
    {
        std::vector<std::string> tg = {"tgtA", "tgtB"};
        std::vector<uint64_t> Nv = {3, 1};
        std::vector<uint64_t> kf = {2,2,2,1,0,  1,1,0,0,0};   // 2 targets x 5 tau, row-major
        std::vector<double>   ll = {5.189665109725213, 2.061107285416508};
        double tg_tau[5] = {0.3, 0.4, 0.5, 0.6, 0.7};
        const char* sp = "/projects/caeg/scratch/kbd606/tmp/convergence_v63_rebuild/authstats_test.authstats";
        assert(write_auth_stats(sp, "Lx", tg_tau, 5, tg, Nv, kf, ll, 19));
        std::puts("DCX_AUTHSTATS_WRITE_OK");
    }

    // The uncalibrated-write guard must fire (all-zero calib_contract_id).
    DcxWriteOptions bad = opt;
    std::memset(bad.calib_contract_id, 0, 16);
    bool threw = false;
    try { write_dcx("/tmp/dcx_bad.dcx", B, bad); } catch (const std::runtime_error&) { threw = true; }
    assert(threw);

    std::remove(path);
    std::puts("DCX_ROUNDTRIP_OK");
    return 0;
}
