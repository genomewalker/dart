#pragma once
/**
 * @file emi_v7_io.hpp
 * @brief EMI v7 writer + reader (commit 1: format round-trip + per-read random access).
 *
 * Header-only so the round-trip fixture builds standalone before CMake wiring.
 * Scope of this commit: block-aligned partitioning, degree-CSR, hot/cold column
 * split with per-block ZSTD, variable-length qaln/taln, embedded per-read annot,
 * and O(1) per-read gather. Name-dict front-coding (commit 2), native mmseqs
 * ingest (commit 4) and annotate migration (commit 5) land later.
 */

#include "dart/emi_v7.hpp"

#include <zstd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include <atomic>
#include <functional>

namespace dart {
namespace emiv7 {

// -------- in-memory row (writer input) — full annotate column set --------
struct Row {
    uint32_t read_idx = 0;      // dense read id
    uint32_t ref_idx  = 0;      // dense ref id   [HOT]
    float    bit_score = 0.0f;  //                [HOT]
    float    damage_score = 0.0f; //              [HOT]
    float    evalue_log10 = 0.0f; //              [COLD]
    uint16_t identity_q = 0;    // fident*65535   [COLD]
    uint16_t aln_len = 0;       //                [COLD]
    uint16_t qstart = 0, qend = 0;                // [COLD]
    uint16_t tstart = 0, tend = 0, tlen = 0, qlen = 0;  // [COLD]
    uint16_t dmg_k = 0, dmg_m = 0;                // [COLD] damage evidence (v5 parity)
    float    dmg_ll_a = 0.0f, dmg_ll_m = 0.0f;    // [COLD]
    std::string qaln, taln;     //                [COLD var-length]
};

// -------- tiny ZSTD helpers (dart links libzstd) --------
namespace detail {
inline std::vector<char> zc(const void* p, size_t n, int level = 3) {
    size_t bound = ZSTD_compressBound(n);
    std::vector<char> out(bound);
    size_t r = ZSTD_compress(out.data(), bound, p, n, level);
    if (ZSTD_isError(r)) throw std::runtime_error("zstd compress");
    out.resize(r);
    return out;
}
inline std::vector<char> zd(const void* p, size_t n, size_t dst) {
    std::vector<char> out(dst);
    size_t r = ZSTD_decompress(out.data(), dst, p, n);
    if (ZSTD_isError(r) || r != dst) throw std::runtime_error("zstd decompress");
    return out;
}
// append raw POD to a byte vector
template <class T> inline void put(std::vector<char>& b, const T& v) {
    const char* p = reinterpret_cast<const char*>(&v);
    b.insert(b.end(), p, p + sizeof(T));
}
template <class T> inline T get(const char* base, uint64_t& off) {
    T v; std::memcpy(&v, base + off, sizeof(T)); off += sizeof(T); return v;
}
// ---- front-coded name dictionary (ID order) ----
// Section layout: [u32 n_names][u32 n_blocks][BlockDir{u64 off,u32 raw,u32 comp}*n_blocks][zstd blocks...]
// Within a block, each name = [u16 shared_prefix_len][u16 suffix_len][suffix bytes],
// prefix taken from the PREVIOUS name in the block (first name in block: shared=0).
struct NameBlockDir { uint64_t off; uint32_t raw; uint32_t comp; };

inline std::vector<char> encode_name_dict(const std::vector<std::string>& names,
                                          uint32_t block_size) {
    std::vector<char> body;                 // concatenated zstd blocks
    std::vector<NameBlockDir> dir;
    const uint32_t n = static_cast<uint32_t>(names.size());
    for (uint32_t start = 0; start < n; start += block_size) {
        uint32_t end = std::min(start + block_size, n);
        std::vector<char> raw;
        std::string_view prev;
        for (uint32_t i = start; i < end; ++i) {
            std::string_view cur = names[i];
            uint16_t shared = 0;
            if (i != start) {
                size_t m = std::min(prev.size(), cur.size());
                while (shared < m && shared < 0xFFFF && prev[shared] == cur[shared]) shared++;
            }
            uint16_t suf = static_cast<uint16_t>(cur.size() - shared);
            put(raw, shared); put(raw, suf);
            raw.insert(raw.end(), cur.begin() + shared, cur.end());
            prev = cur;
        }
        auto c = zc(raw.data(), raw.size(), 19);   // dicts are cold: max ratio
        bool use = c.size() < raw.size();
        NameBlockDir d{ body.size(), (uint32_t)raw.size(), use ? (uint32_t)c.size() : (uint32_t)raw.size() };
        dir.push_back(d);
        if (use) body.insert(body.end(), c.begin(), c.end());
        else body.insert(body.end(), raw.begin(), raw.end());
    }
    std::vector<char> sec;
    put<uint32_t>(sec, n);
    put<uint32_t>(sec, block_size);
    put<uint32_t>(sec, (uint32_t)dir.size());
    for (auto& d : dir) { put(sec, d.off); put(sec, d.raw); put(sec, d.comp); }
    sec.insert(sec.end(), body.begin(), body.end());
    return sec;
}

// Cache identity for the thread_local caches below. A `this` pointer is NOT usable
// as a key: an object can be freed and a new one allocated at the same address, and
// the cache would then serve the previous file's data. Values are never reused.
inline uint64_t next_instance_uid() { static std::atomic<uint64_t> c{0}; return ++c; }

// Reader for a front-coded name dict section. Lazily decodes + caches one block.
class NameDict {
public:
    NameDict() = default;
    NameDict(const char* base, uint64_t off, uint64_t size) : uid_(next_instance_uid()) {
        if (size == 0) return;
        sec_ = base + off;
        uint64_t o = 0;
        n_ = get<uint32_t>(sec_, o);
        block_size_ = get<uint32_t>(sec_, o);
        uint32_t nblk = get<uint32_t>(sec_, o);
        dir_.resize(nblk);
        for (uint32_t i = 0; i < nblk; ++i) { dir_[i].off = get<uint64_t>(sec_, o); dir_[i].raw = get<uint32_t>(sec_, o); dir_[i].comp = get<uint32_t>(sec_, o); }
        body_ = sec_ + o;
    }
    uint32_t size() const { return n_; }
    // Thread-safe: decodes the owning block into a per-thread cache. Returns an
    // owned std::string (callers must NOT bind the result to a string_view).
    std::string name(uint32_t i) const {
        if (i >= n_) throw std::out_of_range("name idx");
        uint32_t b = i / block_size_;
        static thread_local std::vector<std::string> tnames;
        static thread_local uint64_t owner = 0;
        static thread_local uint32_t tblk = 0xFFFFFFFF;
        if (owner != uid_ || tblk != b) { decode_block(b, tnames); owner = uid_; tblk = b; }
        return tnames[i - b * block_size_];
    }
private:
    void decode_block(uint32_t b, std::vector<std::string>& out) const {
        const NameBlockDir& d = dir_[b];
        std::vector<char> cache;
        if (d.comp == d.raw) cache.assign(body_ + d.off, body_ + d.off + d.raw);
        else cache = zd(body_ + d.off, d.comp, d.raw);
        out.clear();
        uint64_t o = 0; std::string prev;
        while (o < cache.size()) {
            uint16_t shared = get<uint16_t>(cache.data(), o);
            uint16_t suf = get<uint16_t>(cache.data(), o);
            std::string cur = prev.substr(0, shared);
            cur.append(cache.data() + o, suf); o += suf;
            out.push_back(cur); prev = cur;
        }
    }
    const char* sec_ = nullptr; const char* body_ = nullptr;
    uint32_t n_ = 0, block_size_ = 1;
    uint64_t uid_ = 0;
    std::vector<NameBlockDir> dir_;
};

// On-disk chunk framing: [u8 codec][u64 raw_size][u64 comp_size][bytes]
inline void write_chunk(std::vector<char>& b, const void* raw, size_t raw_size) {
    if (raw_size == 0) { put<uint8_t>(b, 0); put<uint64_t>(b, 0); put<uint64_t>(b, 0); return; }
    auto c = zc(raw, raw_size);
    if (c.size() < raw_size) {
        put<uint8_t>(b, static_cast<uint8_t>(Codec::ZSTD));
        put<uint64_t>(b, raw_size); put<uint64_t>(b, c.size());
        b.insert(b.end(), c.begin(), c.end());
    } else {
        put<uint8_t>(b, static_cast<uint8_t>(Codec::NONE));
        put<uint64_t>(b, raw_size); put<uint64_t>(b, raw_size);
        const char* p = reinterpret_cast<const char*>(raw);
        b.insert(b.end(), p, p + raw_size);
    }
}
inline std::vector<char> read_chunk(const char* base, uint64_t& off) {
    uint8_t codec = get<uint8_t>(base, off);
    uint64_t raw = get<uint64_t>(base, off);
    uint64_t comp = get<uint64_t>(base, off);
    if (raw == 0) return {};
    std::vector<char> out;
    if (codec == static_cast<uint8_t>(Codec::NONE)) out.assign(base + off, base + off + raw);
    else out = zd(base + off, comp, raw);
    off += comp;
    return out;
}

// Materialize ONE read-aligned block's column chunks (hot then cold, exact chunk
// order the reader's decode_block expects) and append them to block_bytes. Shared
// by the in-RAM Writer and the streaming merge so block payload bytes are
// byte-identical between the two paths (fable byte-identity ruling). `rows` points
// at the block's first row; the caller guarantees rows[0..bn) are the block's rows
// already in final total-order.
inline void encode_block_columns(const Row* rows, uint64_t bn, std::vector<char>& block_bytes) {
    std::vector<uint32_t> c_ref(bn); std::vector<float> c_bit(bn), c_dmg(bn), c_ev(bn);
    std::vector<uint16_t> c_iq(bn), c_al(bn), c_qs(bn), c_qe(bn), c_ts(bn), c_te(bn), c_tl(bn), c_ql(bn);
    std::vector<uint16_t> c_dk(bn), c_dm(bn); std::vector<float> c_la(bn), c_lm(bn);
    std::vector<uint32_t> qoff(bn + 1, 0), toff(bn + 1, 0);
    std::string qdata, tdata;
    for (uint64_t j = 0; j < bn; ++j) {
        const Row& r = rows[j];
        c_ref[j] = r.ref_idx; c_bit[j] = r.bit_score; c_dmg[j] = r.damage_score;
        c_ev[j] = r.evalue_log10;
        c_iq[j] = r.identity_q; c_al[j] = r.aln_len; c_qs[j] = r.qstart; c_qe[j] = r.qend;
        c_ts[j] = r.tstart; c_te[j] = r.tend; c_tl[j] = r.tlen; c_ql[j] = r.qlen;
        c_dk[j] = r.dmg_k; c_dm[j] = r.dmg_m; c_la[j] = r.dmg_ll_a; c_lm[j] = r.dmg_ll_m;
        qdata += r.qaln; qoff[j + 1] = static_cast<uint32_t>(qdata.size());
        tdata += r.taln; toff[j + 1] = static_cast<uint32_t>(tdata.size());
    }
    // HOT chunks
    write_chunk(block_bytes, c_ref.data(), bn * sizeof(uint32_t));
    write_chunk(block_bytes, c_bit.data(), bn * sizeof(float));
    write_chunk(block_bytes, c_dmg.data(), bn * sizeof(float));
    // COLD chunks (order must mirror decode_block)
    write_chunk(block_bytes, c_ev.data(), bn * sizeof(float));
    write_chunk(block_bytes, c_iq.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_al.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_qs.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_qe.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_ts.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_te.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_tl.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_ql.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_dk.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_dm.data(), bn * sizeof(uint16_t));
    write_chunk(block_bytes, c_la.data(), bn * sizeof(float));
    write_chunk(block_bytes, c_lm.data(), bn * sizeof(float));
    write_chunk(block_bytes, qoff.data(), (bn + 1) * sizeof(uint32_t));
    write_chunk(block_bytes, qdata.data(), qdata.size());
    write_chunk(block_bytes, toff.data(), (bn + 1) * sizeof(uint32_t));
    write_chunk(block_bytes, tdata.data(), tdata.size());
}

// Read-only whole-file mmap with RAII teardown. Replaces the heap slurp the
// Reader used to do: a 130-200 GB .emi read into std::vector<char> is 130-200 GB
// of *anonymous* memory (resize() value-initializes, so every page is dirtied and
// resident before a single byte is read from disk). Mapped file pages are clean
// and page-cache backed, so the kernel can evict them under pressure and
// drop_cache() can release them outright.
class MappedFile {
public:
    MappedFile() = default;
    explicit MappedFile(const std::string& path) {
        fd_ = ::open(path.c_str(), O_RDONLY);
        if (fd_ < 0) throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
        struct stat st{};
        if (::fstat(fd_, &st) != 0) {
            int e = errno; ::close(fd_); fd_ = -1;
            throw std::runtime_error("cannot stat " + path + ": " + std::strerror(e));
        }
        size_ = static_cast<size_t>(st.st_size);
        if (size_ == 0) return;   // nothing to map; data() stays null
        void* p = ::mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd_, 0);
        if (p == MAP_FAILED) {
            int e = errno; ::close(fd_); fd_ = -1; size_ = 0;
            throw std::runtime_error("mmap failed for " + path + ": " + std::strerror(e));
        }
        data_ = static_cast<const char*>(p);
        // No madvise here: kernel default (MADV_NORMAL) readahead is what we want.
        // MADV_RANDOM was tried and measured — it disables readahead, so on NFS every
        // 4 KB fault becomes its own round trip, and the explicit per-block
        // MADV_WILLNEED that compensated for it bought nothing back. Same library,
        // same node, same window: with advice+prefetch 46:49.22, without either
        // 46:44.96 (0.15%), identical peak anon (50.99 GiB) and file (43.59 GiB).
    }
    ~MappedFile() {
        if (data_) ::munmap(const_cast<char*>(data_), size_);
        if (fd_ >= 0) ::close(fd_);
    }
    MappedFile(const MappedFile&) = delete;
    MappedFile& operator=(const MappedFile&) = delete;

    const char* data() const { return data_; }
    size_t size() const { return size_; }

    // Drop every resident page of the mapping. MADV_DONTNEED on a read-only
    // MAP_PRIVATE file mapping unmaps the pages (a later access re-faults them from
    // the file); POSIX_FADV_DONTNEED then evicts the underlying page cache, which
    // madvise alone does not do.
    void drop_cache() const {
        if (!data_) return;
        // data_ comes from mmap, so it is already page-aligned here.
        advise(const_cast<char*>(data_), size_, MADV_DONTNEED, "MADV_DONTNEED");
        // posix_fadvise returns the error number; it does NOT set errno.
        int rc = ::posix_fadvise(fd_, 0, 0, POSIX_FADV_DONTNEED);
        if (rc != 0) warn_once(fadv_warned_, "POSIX_FADV_DONTNEED", rc);
    }

private:
    // These are advisory hints: a failure must not abort the read, but it must not
    // pass unnoticed either. Warn once per advice kind.
    static void warn_once(std::atomic<bool>& flag, const char* what, int err) {
        bool expected = false;
        if (flag.compare_exchange_strong(expected, true))
            std::cerr << "warning: " << what << " failed: " << std::strerror(err)
                      << " (performance hint only; results unaffected)\n";
    }
    static void advise(void* p, size_t n, int adv, const char* what) {
        if (::madvise(p, n, adv) != 0) warn_once(madv_warned_, what, errno);
    }
    inline static std::atomic<bool> madv_warned_{false};
    inline static std::atomic<bool> fadv_warned_{false};

    const char* data_ = nullptr;
    size_t size_ = 0;
    int fd_ = -1;
};

}  // namespace detail

// ================= WRITER =================
class Writer {
public:
    // rows: any order; annots + read_names indexed by dense read id; ref_names by dense ref id.
    static void write(const std::string& path,
                      std::vector<Row> rows,
                      const std::vector<EmiReadAnnot>& annots,
                      const std::vector<std::string>& read_names,
                      const std::vector<std::string>& ref_names,
                      uint32_t block_target = EMI_V7_BLOCK_TARGET_ROWS,
                      float d_max = 0.0f, float lambda = 0.0f) {
        using namespace detail;
        const uint64_t num_reads = annots.size();

        // TOTAL-ORDER sort (fable DECISION A): (read_idx asc, bit_score desc,
        // ref_idx asc, original-row-index asc). The trailing original-index key is
        // globally unique, so the order is total and permutation-invariant even when
        // mmseqs emits multiple equal-bit HSPs to the SAME target (read_idx,ref_idx
        // collide). std::sort is unstable -> without this key the on-disk tie order
        // would be a non-reproducible artifact of libstdc++ internals. We sort an
        // index permutation so the comparator's final tie-break IS the original row
        // index (a<b), then materialize; no per-row index field is serialized.
        {
            const uint64_t n = rows.size();
            std::vector<uint64_t> perm(n);
            for (uint64_t k = 0; k < n; ++k) perm[k] = k;
            std::sort(perm.begin(), perm.end(), [&rows](uint64_t a, uint64_t b) {
                const Row& ra = rows[a]; const Row& rb = rows[b];
                if (ra.read_idx != rb.read_idx) return ra.read_idx < rb.read_idx;
                if (ra.bit_score != rb.bit_score) return ra.bit_score > rb.bit_score;
                if (ra.ref_idx != rb.ref_idx) return ra.ref_idx < rb.ref_idx;
                return a < b;   // original row index: guaranteed-unique final key
            });
            std::vector<Row> sorted;
            sorted.reserve(n);
            for (uint64_t idx : perm) sorted.push_back(std::move(rows[idx]));
            rows.swap(sorted);
        }

        // per-read degree (dense reads, missing reads => degree 0)
        std::vector<uint32_t> degree(num_reads, 0);
        for (const auto& r : rows) degree[r.read_idx]++;

        // read-aligned blocks: extend a block until it reaches block_target AND
        // the current read is complete (never split a read across a block).
        std::vector<BlockMeta> blocks;
        std::vector<uint32_t> read_to_block(num_reads, 0);
        std::vector<char> block_bytes;   // BLOCK_CHUNKS payload
        std::vector<uint64_t> block_chunk_off;

        uint64_t i = 0;                  // row cursor
        const uint64_t nrows = rows.size();
        while (i < nrows) {
            uint64_t bstart = i;
            uint32_t first_read = rows[i].read_idx;
            uint32_t last_read = first_read;
            // gather rows until >= target, then finish the trailing read
            while (i < nrows) {
                uint32_t cur = rows[i].read_idx;
                uint64_t taken = i - bstart;
                if (taken >= block_target && cur != last_read) break;
                last_read = cur;
                i++;
            }
            uint64_t bn = i - bstart;
            BlockMeta bm{};
            bm.start_row = bstart;
            bm.num_rows = static_cast<uint32_t>(bn);
            bm.first_read_idx = first_read;
            bm.last_read_idx = last_read;
            bm.chunks_offset = block_bytes.size();   // relative; absolutized below
            block_chunk_off.push_back(block_bytes.size());

            // materialize this block's columns (shared with streaming merge)
            encode_block_columns(&rows[bstart], bn, block_bytes);

            uint32_t bidx = static_cast<uint32_t>(blocks.size());
            for (uint32_t rd = first_read; rd <= last_read; ++rd) read_to_block[rd] = bidx;
            blocks.push_back(bm);
        }

        // ---- degree section: u16 per read + u32 escape side array ----
        std::vector<char> deg_raw;
        std::vector<uint32_t> escapes;
        for (uint64_t r = 0; r < num_reads; ++r) {
            if (degree[r] >= EMI_V7_DEGREE_ESCAPE) {
                uint16_t esc = EMI_V7_DEGREE_ESCAPE; put(deg_raw, esc);
                escapes.push_back(degree[r]);
            } else {
                uint16_t d = static_cast<uint16_t>(degree[r]); put(deg_raw, d);
            }
        }
        std::vector<char> deg_sec;
        put<uint64_t>(deg_sec, escapes.size());
        for (uint32_t e : escapes) put(deg_sec, e);
        write_chunk(deg_sec, deg_raw.data(), deg_raw.size());

        // ---- read_to_block section ----
        std::vector<char> r2b_sec;
        write_chunk(r2b_sec, read_to_block.data(), num_reads * sizeof(uint32_t));

        // ---- read_annot section ----
        std::vector<char> annot_sec;
        write_chunk(annot_sec, annots.data(), annots.size() * sizeof(EmiReadAnnot));

        // ---- name dictionaries (front-coded + ZSTD) ----
        std::vector<char> rdict = encode_name_dict(read_names, EMI_V7_NAME_DICT_BLOCK);
        std::vector<char> fdict = encode_name_dict(ref_names, EMI_V7_NAME_DICT_BLOCK);

        // ---- block dir (absolutize chunk offsets after we know layout) ----
        // File order: header|block_dir|degree|r2b|annot|block_chunks|rdict|fdict|footer
        uint64_t off = sizeof(EMIv7Header);
        uint64_t block_dir_off = off; uint64_t block_dir_sz = blocks.size() * sizeof(BlockMeta);
        off += block_dir_sz;
        uint64_t deg_off = off; off += deg_sec.size();
        uint64_t r2b_off = off; off += r2b_sec.size();
        uint64_t annot_off = off; off += annot_sec.size();
        uint64_t chunks_off = off; off += block_bytes.size();
        uint64_t rdict_off = off; off += rdict.size();
        uint64_t fdict_off = off; off += fdict.size();

        for (size_t b = 0; b < blocks.size(); ++b)
            blocks[b].chunks_offset = chunks_off + block_chunk_off[b];

        EMIv7Header h{};
        h.flags = EMI_V7_FLAG_SORTED_BY_READ | EMI_V7_FLAG_READ_ALIGNED_BLKS |
                  EMI_V7_FLAG_DEGREE_CSR | EMI_V7_FLAG_SELF_CONTAINED |
                  EMI_V7_FLAG_HAS_ALN_STRINGS | EMI_V7_FLAG_COMPRESSED_DICT;
        h.num_alignments = nrows; h.num_reads = num_reads;
        h.num_refs = ref_names.size();
        h.d_max = d_max; h.lambda = lambda;   // scoring fingerprint (parity w/ hits2emi)
        h.num_blocks = static_cast<uint32_t>(blocks.size());
        h.block_target = block_target;
        h.sections[(int)SectionID::BLOCK_DIR]      = {block_dir_off, block_dir_sz};
        h.sections[(int)SectionID::READ_DEGREE]    = {deg_off, deg_sec.size()};
        h.sections[(int)SectionID::READ_TO_BLOCK]  = {r2b_off, r2b_sec.size()};
        h.sections[(int)SectionID::READ_ANNOT]     = {annot_off, annot_sec.size()};
        h.sections[(int)SectionID::BLOCK_CHUNKS]   = {chunks_off, block_bytes.size()};
        h.sections[(int)SectionID::NAME_DICT_READ] = {rdict_off, rdict.size()};
        h.sections[(int)SectionID::NAME_DICT_REF]  = {fdict_off, fdict.size()};

        std::ofstream out(path, std::ios::binary);
        if (!out) throw std::runtime_error("cannot create " + path);
        out.write(reinterpret_cast<const char*>(&h), sizeof(h));
        out.write(reinterpret_cast<const char*>(blocks.data()), block_dir_sz);
        out.write(deg_sec.data(), deg_sec.size());
        out.write(r2b_sec.data(), r2b_sec.size());
        out.write(annot_sec.data(), annot_sec.size());
        out.write(block_bytes.data(), block_bytes.size());
        out.write(rdict.data(), rdict.size());
        out.write(fdict.data(), fdict.size());
        // Footer carries the trailing EMI_V7_MAGIC (default-initialized); the reader
        // validates it to catch truncation. header_xxh is reserved (0), not computed.
        EMIv7Footer f{};
        out.write(reinterpret_cast<const char*>(&f), sizeof(f));
    }
};

// ================= READER =================
class Reader {
public:
    explicit Reader(const std::string& path) : buf_(path) {
        const size_t sz = buf_.size();
        // Size check first: the header memcpy below indexes buf_.data(), which is
        // null for an empty file and short of 256 bytes for a stub.
        if (sz < sizeof(EMIv7Header) + sizeof(EMIv7Footer))
            throw std::runtime_error("EMI v7 file truncated (no footer): " + path);
        std::memcpy(&h_, buf_.data(), sizeof(h_));
        if (h_.magic != EMI_V7_MAGIC || h_.version != EMI_V7_VERSION)
            throw std::runtime_error("not an EMI v7 file");
        // Footer integrity: the trailing magic catches truncation/corruption.
        EMIv7Footer ftr;
        std::memcpy(&ftr, buf_.data() + sz - sizeof(EMIv7Footer), sizeof(EMIv7Footer));
        if (ftr.magic != EMI_V7_MAGIC)
            throw std::runtime_error("EMI v7 file truncated or corrupt (bad footer magic)");
        const auto& bd = h_.sections[(int)SectionID::BLOCK_DIR];
        blocks_ = reinterpret_cast<const BlockMeta*>(buf_.data() + bd.offset);
        // decode degree -> CSR
        decode_degree();
        // read_to_block
        uint64_t o = h_.sections[(int)SectionID::READ_TO_BLOCK].offset;
        auto r2b = detail::read_chunk(buf_.data(), o);
        read_to_block_.assign(reinterpret_cast<uint32_t*>(r2b.data()),
                              reinterpret_cast<uint32_t*>(r2b.data()) + h_.num_reads);
        // annot
        uint64_t ao = h_.sections[(int)SectionID::READ_ANNOT].offset;
        auto ab = detail::read_chunk(buf_.data(), ao);
        annots_.assign(reinterpret_cast<EmiReadAnnot*>(ab.data()),
                       reinterpret_cast<EmiReadAnnot*>(ab.data()) + h_.num_reads);
        // name dicts
        const auto& rd = h_.sections[(int)SectionID::NAME_DICT_READ];
        const auto& fd = h_.sections[(int)SectionID::NAME_DICT_REF];
        read_dict_ = detail::NameDict(buf_.data(), rd.offset, rd.size);
        ref_dict_  = detail::NameDict(buf_.data(), fd.offset, fd.size);
    }

    uint64_t num_alignments() const { return h_.num_alignments; }
    uint64_t num_reads() const { return h_.num_reads; }
    uint64_t num_refs() const { return h_.num_refs; }
    uint32_t num_blocks() const { return h_.num_blocks; }
    uint32_t degree(uint32_t read) const { return csr_[read + 1] - csr_[read]; }
    const EmiReadAnnot& annot(uint32_t read) const { return annots_[read]; }
    std::string read_name(uint32_t read) const { return read_dict_.name(read); }
    std::string ref_name(uint32_t ref) const { return ref_dict_.name(ref); }

    // ---- ColumnarIndexReader-compatible surface (annotate reads v7 directly) ----
    uint32_t num_row_groups() const { return h_.num_blocks; }
    uint32_t version() const { return h_.version; }
    bool is_valid() const { return h_.magic == EMI_V7_MAGIC; }
    bool has_alignment_strings() const { return (h_.flags & EMI_V7_FLAG_HAS_ALN_STRINGS) != 0; }
    bool is_sorted_by_read_then_score() const { return (h_.flags & EMI_V7_FLAG_SORTED_BY_READ) != 0; }
    bool has_row_group_read_index() const { return true; }
    float d_max() const { return h_.d_max; }
    float lambda() const { return h_.lambda; }
    uint64_t read_offset(uint32_t read) const { return csr_[read]; }
    uint32_t read_degree(uint32_t read) const { return degree(read); }
    void set_columns(std::initializer_list<ColumnID>) const {}   // v7 decodes selectively per block
    void set_all_columns() const {}
    void add_filter(const FilterPredicate&) const {}
    void clear_filters() const {}
    // Release every resident page of the mapped .emi. Callers use this between
    // passes (cmd_damage_annotate.cpp:2656, :3148) to hand the page cache back
    // before the next phase allocates.
    void flush_pages() const { buf_.drop_cache(); }
    std::vector<std::string> ref_names_pread() const {
        std::vector<std::string> v(h_.num_refs);
        for (uint32_t i = 0; i < h_.num_refs; ++i) v[i] = ref_dict_.name(i);
        return v;
    }

    using RowGroupCallback = std::function<void(
        uint32_t, uint32_t,
        const uint32_t*, const uint32_t*, const float*, const float*, const float*,
        const uint16_t*, const uint16_t*, const float*, const float*,
        const uint16_t*, const uint16_t*, const uint16_t*, const uint16_t*,
        const uint16_t*, const uint16_t*, const uint16_t*, const uint16_t*,
        const uint16_t*, const uint16_t*)>;

    // Numeric pass: skip qaln/taln decode (callbacks pull them via get_qaln, which
    // re-decodes with aln on demand). Selected pass: decode aln up front since its
    // callback calls get_qaln for every selected row.
    void parallel_scan(const RowGroupCallback& cb) const { scan_range(0, h_.num_blocks, cb, /*with_aln=*/false); }
    void parallel_scan_selected(const std::vector<uint32_t>& idx, const RowGroupCallback& cb) const {
        #pragma omp parallel for schedule(dynamic,1)
        for (size_t k = 0; k < idx.size(); ++k) invoke_block(idx[k], cb, /*with_aln=*/true);
    }

    // qaln/taln for a row within a block. The view points into the caller thread's
    // decoded-block cache (tls_decoded), so it is thread-safe inside a parallel
    // callback and valid until that thread decodes a different block. Callers use
    // it call-scoped (annotate cmd_damage_annotate.cpp:2662, :3729).
    std::string_view get_qaln(uint32_t rg, uint32_t row) const {
        const BlockDecoded& d = tls_decoded(rg, /*with_aln=*/true);
        return std::string_view(d.qdata.data() + d.qoff[row], d.qoff[row + 1] - d.qoff[row]);
    }
    std::string_view get_taln(uint32_t rg, uint32_t row) const {
        const BlockDecoded& d = tls_decoded(rg, /*with_aln=*/true);
        return std::string_view(d.tdata.data() + d.toff[row], d.toff[row + 1] - d.toff[row]);
    }

    // Random access: gather all rows for one read. Decompresses its owning block
    // and slices out that read's contiguous rows (hot + cold + var-length).
    std::vector<Row> rows_for_read(uint32_t read) const {
        std::vector<Row> out;
        if (degree(read) == 0) return out;
        uint32_t bidx = read_to_block_[read];
        const BlockDecoded& d = tls_decoded(bidx);   // per-thread cache: O(1) for reads in the same block
        // rows for `read` are contiguous within the block; slice via CSR offsets
        uint64_t lo = csr_[read] - csr_[d.first_read];   // local start within block
        uint32_t deg = csr_[read + 1] - csr_[read];
        out.reserve(deg);
        for (uint32_t k = 0; k < deg; ++k) out.push_back(d.row(static_cast<uint32_t>(lo + k)));
        return out;
    }

    // Selective HOT-only scan of a block (EM path): ref_idx + bit_score + damage_score.
    struct HotBlock { uint32_t n; std::vector<uint32_t> ref; std::vector<float> bit, dmg; };
    HotBlock hot_block(uint32_t bidx) const {
        uint64_t o = blocks_[bidx].chunks_offset;
        uint32_t n = blocks_[bidx].num_rows;
        HotBlock hb; hb.n = n;
        auto r = detail::read_chunk(buf_.data(), o);
        hb.ref.assign(reinterpret_cast<uint32_t*>(r.data()), reinterpret_cast<uint32_t*>(r.data()) + n);
        auto b = detail::read_chunk(buf_.data(), o);
        hb.bit.assign(reinterpret_cast<float*>(b.data()), reinterpret_cast<float*>(b.data()) + n);
        auto g = detail::read_chunk(buf_.data(), o);
        hb.dmg.assign(reinterpret_cast<float*>(g.data()), reinterpret_cast<float*>(g.data()) + n);
        return hb;
    }

private:
    struct BlockDecoded {
        uint32_t n = 0;
        uint32_t first_read = 0;
        std::vector<uint32_t> ref; std::vector<float> bit, dmg, ev;
        std::vector<uint16_t> iq, al, qs, qe, ts, te, tl, ql, dk, dm;
        std::vector<float> la, lm;
        std::vector<uint32_t> qoff, toff; std::string qdata, tdata;
        std::vector<uint32_t> read_id;   // per-row read id, reconstructed from CSR
        uint32_t read_of(uint32_t j, const BlockMeta&) const { return read_id[j]; }
        Row row(uint32_t j) const {
            Row r; r.read_idx = read_id[j]; r.ref_idx = ref[j];
            r.bit_score = bit[j]; r.damage_score = dmg[j]; r.evalue_log10 = ev[j];
            r.identity_q = iq[j]; r.aln_len = al[j]; r.qstart = qs[j]; r.qend = qe[j];
            r.tstart = ts[j]; r.tend = te[j]; r.tlen = tl[j]; r.qlen = ql[j];
            r.dmg_k = dk[j]; r.dmg_m = dm[j]; r.dmg_ll_a = la[j]; r.dmg_ll_m = lm[j];
            r.qaln.assign(qdata.data() + qoff[j], qdata.data() + qoff[j + 1]);
            r.taln.assign(tdata.data() + toff[j], tdata.data() + toff[j + 1]);
            return r;
        }
    };

    // Thread-safe decoded-block access. Each thread caches its most-recently
    // decoded block in a thread_local keyed by (uid_, bidx, with_aln). invoke_block
    // and get_qaln/get_taln share it, so a parallel callback decodes each block at
    // most once per thread with NO shared mutable state (no data race). If a block
    // was cached without alignment strings but aln is now required, it re-decodes.
    // uid_ (not `this`) is the key: a freed Reader's address can be reused by a new
    // Reader, and pointer keying would then serve the previous file's cached block.
    const BlockDecoded& tls_decoded(uint32_t bidx, bool with_aln = true) const {
        static thread_local BlockDecoded cache;
        static thread_local uint64_t owner = 0;
        static thread_local uint32_t cached_blk = 0xFFFFFFFFu;
        static thread_local bool cached_aln = false;
        if (owner != uid_ || cached_blk != bidx || (with_aln && !cached_aln)) {
            cache = decode_block(bidx, with_aln);
            owner = uid_; cached_blk = bidx; cached_aln = with_aln;
        }
        return cache;
    }

    // with_aln=false skips qaln/taln (var-length cold columns) for pure numeric
    // scans (annotate pass 1) — those blocks pull alignment strings on demand via
    // get_qaln, which re-decodes with aln through tls_decoded.
    BlockDecoded decode_block(uint32_t bidx, bool with_aln = true) const {
        uint64_t o = blocks_[bidx].chunks_offset;
        const BlockMeta& bm = blocks_[bidx];
        BlockDecoded d; d.n = bm.num_rows; d.first_read = bm.first_read_idx;
        auto take_u32 = [&](std::vector<uint32_t>& v, uint32_t n){ auto c=detail::read_chunk(buf_.data(),o); v.assign((uint32_t*)c.data(),(uint32_t*)c.data()+n); };
        auto take_f = [&](std::vector<float>& v, uint32_t n){ auto c=detail::read_chunk(buf_.data(),o); v.assign((float*)c.data(),(float*)c.data()+n); };
        auto take_u16 = [&](std::vector<uint16_t>& v, uint32_t n){ auto c=detail::read_chunk(buf_.data(),o); v.assign((uint16_t*)c.data(),(uint16_t*)c.data()+n); };
        auto take_str = [&](std::string& s){ auto c=detail::read_chunk(buf_.data(),o); s.assign(c.data(),c.size()); };
        take_u32(d.ref, d.n); take_f(d.bit, d.n); take_f(d.dmg, d.n);
        take_f(d.ev, d.n);
        take_u16(d.iq, d.n); take_u16(d.al, d.n); take_u16(d.qs, d.n); take_u16(d.qe, d.n);
        take_u16(d.ts, d.n); take_u16(d.te, d.n); take_u16(d.tl, d.n); take_u16(d.ql, d.n);
        take_u16(d.dk, d.n); take_u16(d.dm, d.n); take_f(d.la, d.n); take_f(d.lm, d.n);
        if (with_aln) {
            take_u32(d.qoff, d.n + 1); take_str(d.qdata);
            take_u32(d.toff, d.n + 1); take_str(d.tdata);
        }
        // reconstruct per-row read ids from CSR within this block's read span
        d.read_id.resize(d.n);
        uint64_t row = bm.start_row;
        for (uint32_t rd = bm.first_read_idx; rd <= bm.last_read_idx; ++rd) {
            uint32_t deg = csr_[rd + 1] - csr_[rd];
            for (uint32_t k = 0; k < deg; ++k) d.read_id[row - bm.start_row + k] = rd;
            row += deg;
        }
        return d;
    }

    void decode_degree() {
        uint64_t o = h_.sections[(int)SectionID::READ_DEGREE].offset;
        uint64_t nesc = detail::get<uint64_t>(buf_.data(), o);
        std::vector<uint32_t> esc(nesc);
        for (uint64_t i = 0; i < nesc; ++i) esc[i] = detail::get<uint32_t>(buf_.data(), o);
        auto degraw = detail::read_chunk(buf_.data(), o);
        const uint16_t* deg = reinterpret_cast<const uint16_t*>(degraw.data());
        csr_.assign(h_.num_reads + 1, 0);
        uint64_t ei = 0;
        for (uint64_t r = 0; r < h_.num_reads; ++r) {
            uint32_t d = deg[r];
            if (d == EMI_V7_DEGREE_ESCAPE) d = esc[ei++];
            csr_[r + 1] = csr_[r] + d;
        }
    }

    detail::MappedFile buf_;
    EMIv7Header h_{};
    const BlockMeta* blocks_ = nullptr;
    std::vector<uint32_t> read_to_block_;
    std::vector<EmiReadAnnot> annots_;
    std::vector<uint64_t> csr_;
    detail::NameDict read_dict_, ref_dict_;
    // See next_instance_uid(): keying the decoded-block cache on `this` hands a new
    // Reader the previous Reader's block when the allocator reuses the address.
    const uint64_t uid_ = detail::next_instance_uid();

    // OpenMP-safe block scan: each thread decodes its own block into a thread_local
    // cache (tls_decoded), so there is no shared mutable state across threads.
    void scan_range(uint32_t lo, uint32_t hi, const RowGroupCallback& cb, bool with_aln) const {
        #pragma omp parallel for schedule(dynamic,1)
        for (uint32_t b = lo; b < hi; ++b) invoke_block(b, cb, with_aln);
    }
    void invoke_block(uint32_t b, const RowGroupCallback& cb, bool with_aln) const {
        const BlockDecoded& d = tls_decoded(b, with_aln);   // per-thread, no shared state
        cb(b, d.n, d.read_id.data(), d.ref.data(), d.bit.data(), d.dmg.data(), d.ev.data(),
           d.dk.data(), d.dm.data(), d.la.data(), d.lm.data(), d.iq.data(), d.al.data(),
           d.qs.data(), d.qe.data(), d.ts.data(), d.te.data(), d.tl.data(), d.ql.data(),
           nullptr, nullptr);   // mismatch/gapopen not stored (annotate ignores)
    }
};

}  // namespace emiv7
}  // namespace dart
