// Columnar index implementation with parallel row group processing
//
// Key optimizations:
// - Column-oriented storage for cache efficiency
// - Row groups enable parallel decompression
// - Predicate pushdown via chunk statistics
// - Lazy loading of alignment strings

#include "dart/columnar_index.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <numeric>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef HAVE_ZSTD
#include <zstd.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dart {

namespace {

inline void decompress_chunk_to_buffer(const ColumnChunk& chunk,
                                       const char* src,
                                       std::vector<char>& out) {
    out.resize(chunk.uncompressed_size);

    switch (chunk.codec) {
        case Codec::NONE:
            std::memcpy(out.data(), src, chunk.uncompressed_size);
            return;
        case Codec::ZSTD:
#ifdef HAVE_ZSTD
        {
            size_t got = ZSTD_decompress(out.data(),
                                         static_cast<size_t>(chunk.uncompressed_size),
                                         src,
                                         static_cast<size_t>(chunk.compressed_size));
            if (ZSTD_isError(got) || got != static_cast<size_t>(chunk.uncompressed_size)) {
                throw std::runtime_error("EMI ZSTD decompression failed");
            }
            return;
        }
#else
            throw std::runtime_error("EMI chunk uses ZSTD but binary was built without ZSTD support");
#endif
        default:
            throw std::runtime_error("Unsupported EMI compression codec");
    }
}

inline void madvise_range(const char* addr, size_t len, int advice) {
    if (addr == nullptr || len == 0) return;
    long page_size = sysconf(_SC_PAGESIZE);
    if (page_size <= 0) return;
    const uintptr_t mask = static_cast<uintptr_t>(page_size - 1);
    const uintptr_t raw = reinterpret_cast<uintptr_t>(addr);
    const uintptr_t aligned_start = raw & ~mask;
    const uintptr_t aligned_end = (raw + len + mask) & ~mask;
    if (aligned_end <= aligned_start) return;
    madvise(reinterpret_cast<void*>(aligned_start),
            static_cast<size_t>(aligned_end - aligned_start),
            advice);
}

inline void madvise_willneed_range(const char* addr, size_t len) {
#ifdef MADV_WILLNEED
    madvise_range(addr, len, MADV_WILLNEED);
#else
    (void)addr;
    (void)len;
#endif
}

inline void madvise_dontneed_range(const char* addr, size_t len) {
#ifdef MADV_DONTNEED
    madvise_range(addr, len, MADV_DONTNEED);
#else
    (void)addr;
    (void)len;
#endif
}

}  // namespace

// ---------------------------------------------------------------------------
// FilterPredicate implementation
// ---------------------------------------------------------------------------

bool FilterPredicate::can_skip(const ColumnStats& stats) const {
    // Check if entire chunk can be skipped based on stats
    switch (column) {
        case ColumnID::REF_IDX:
        case ColumnID::ALN_LEN:
        case ColumnID::QSTART:
        case ColumnID::QEND:
        case ColumnID::TSTART:
        case ColumnID::TEND:
        case ColumnID::TLEN:
        case ColumnID::QLEN:
        case ColumnID::MISMATCH:
        case ColumnID::GAPOPEN:
        case ColumnID::DMG_K:
        case ColumnID::DMG_M:
        case ColumnID::IDENTITY_Q:
            // Integer columns
            switch (op) {
                case Op::GE: return stats.max_u < value_u;
                case Op::GT: return stats.max_u <= value_u;
                case Op::LE: return stats.min_u > value_u;
                case Op::LT: return stats.min_u >= value_u;
                case Op::EQ: return stats.min_u > value_u || stats.max_u < value_u;
            }
            break;
        case ColumnID::BIT_SCORE:
        case ColumnID::DAMAGE_SCORE:
        case ColumnID::EVALUE_LOG10:
        case ColumnID::DMG_LL_A:
        case ColumnID::DMG_LL_M:
            // Float columns
            switch (op) {
                case Op::GE: return stats.max_f < value_f;
                case Op::GT: return stats.max_f <= value_f;
                case Op::LE: return stats.min_f > value_f;
                case Op::LT: return stats.min_f >= value_f;
                case Op::EQ: return stats.min_f > value_f || stats.max_f < value_f;
            }
            break;
        default:
            break;
    }
    return false;
}

// ---------------------------------------------------------------------------
// ColumnarIndexWriter implementation
// ---------------------------------------------------------------------------

struct ColumnarIndexWriter::Impl {
    std::string path;
    std::ofstream out;

    // String dictionaries
    std::vector<char> read_names;
    std::vector<uint32_t> read_name_offsets;
    std::vector<char> ref_names;
    std::vector<uint32_t> ref_name_offsets;
    std::unordered_map<std::string, uint32_t> read_name_index;
    std::unordered_map<std::string, uint32_t> ref_name_index;

    // Current row group buffer
    std::vector<uint32_t> buf_read_idx;  // v2: read index per alignment
    std::vector<uint32_t> buf_ref_idx;
    std::vector<float> buf_bit_score;
    std::vector<float> buf_damage_score;
    std::vector<uint16_t> buf_dmg_k;
    std::vector<uint16_t> buf_dmg_m;
    std::vector<float> buf_dmg_ll_a;
    std::vector<float> buf_dmg_ll_m;
    std::vector<float> buf_evalue_log10;
    std::vector<uint16_t> buf_identity_q;
    std::vector<uint16_t> buf_aln_len;
    std::vector<uint16_t> buf_qstart;
    std::vector<uint16_t> buf_qend;
    std::vector<uint16_t> buf_tstart;
    std::vector<uint16_t> buf_tend;
    std::vector<uint16_t> buf_tlen;
    std::vector<uint16_t> buf_qlen;
    std::vector<uint16_t> buf_mismatch;
    std::vector<uint16_t> buf_gapopen;
    std::vector<uint32_t> buf_qaln_offsets;
    std::vector<char> buf_qaln_data;
    std::vector<uint32_t> buf_taln_offsets;
    std::vector<char> buf_taln_data;

    // Row group metadata
    std::vector<RowGroupMeta> row_groups;
    std::vector<RowGroupReadIndex> row_group_read_index;
    uint64_t emitted_rows = 0;

    // CSR offsets (per-read)
    std::vector<uint64_t> csr_offsets;
    uint32_t current_read_idx = UINT32_MAX;
    uint64_t total_alignments = 0;
    bool sorted_by_read_then_score = true;
    bool saw_alignment_strings = false;
    bool has_prev_sort_key = false;
    uint32_t prev_read_idx = 0;
    float prev_bit_score = 0.0f;

    // Temporary file for row group data
    std::string temp_path;
    std::ofstream temp_out;

    Impl(const std::string& p) : path(p), temp_path(p + ".tmp") {
        csr_offsets.push_back(0);
        reserve_buffers();
    }

    void reserve_buffers() {
        buf_read_idx.reserve(ROW_GROUP_SIZE);
        buf_ref_idx.reserve(ROW_GROUP_SIZE);
        buf_bit_score.reserve(ROW_GROUP_SIZE);
        buf_damage_score.reserve(ROW_GROUP_SIZE);
        buf_dmg_k.reserve(ROW_GROUP_SIZE);
        buf_dmg_m.reserve(ROW_GROUP_SIZE);
        buf_dmg_ll_a.reserve(ROW_GROUP_SIZE);
        buf_dmg_ll_m.reserve(ROW_GROUP_SIZE);
        buf_evalue_log10.reserve(ROW_GROUP_SIZE);
        buf_identity_q.reserve(ROW_GROUP_SIZE);
        buf_aln_len.reserve(ROW_GROUP_SIZE);
        buf_qstart.reserve(ROW_GROUP_SIZE);
        buf_qend.reserve(ROW_GROUP_SIZE);
        buf_tstart.reserve(ROW_GROUP_SIZE);
        buf_tend.reserve(ROW_GROUP_SIZE);
        buf_tlen.reserve(ROW_GROUP_SIZE);
        buf_qlen.reserve(ROW_GROUP_SIZE);
        buf_mismatch.reserve(ROW_GROUP_SIZE);
        buf_gapopen.reserve(ROW_GROUP_SIZE);
        buf_qaln_offsets.reserve(ROW_GROUP_SIZE + 1);
        buf_qaln_offsets.push_back(0);
        buf_taln_offsets.reserve(ROW_GROUP_SIZE + 1);
        buf_taln_offsets.push_back(0);
    }

    void flush_row_group() {
        if (buf_read_idx.empty()) return;

        RowGroupMeta rg{};
        rg.num_rows = static_cast<uint32_t>(buf_ref_idx.size());
        rg.first_read_idx = current_read_idx;

        // Open temp file if not open
        if (!temp_out.is_open()) {
            temp_out.open(temp_path, std::ios::binary);
        }

        // Write each column and record metadata
        auto write_column = [&](ColumnID col, const void* data, size_t elem_size, size_t count) {
            auto& chunk = rg.columns[static_cast<size_t>(col)];
            chunk.offset = temp_out.tellp();
            chunk.uncompressed_size = static_cast<uint32_t>(elem_size * count);
            chunk.compressed_size = chunk.uncompressed_size;  // No compression initially
            chunk.codec = Codec::NONE;
            temp_out.write(reinterpret_cast<const char*>(data), chunk.uncompressed_size);
        };

        // Compute stats for numeric columns
        auto compute_stats_u32 = [](const std::vector<uint32_t>& v) -> ColumnStats {
            ColumnStats s;
            s.min_u = *std::min_element(v.begin(), v.end());
            s.max_u = *std::max_element(v.begin(), v.end());
            return s;
        };

        auto compute_stats_u16 = [](const std::vector<uint16_t>& v) -> ColumnStats {
            ColumnStats s;
            s.min_u = *std::min_element(v.begin(), v.end());
            s.max_u = *std::max_element(v.begin(), v.end());
            return s;
        };

        auto compute_stats_f32 = [](const std::vector<float>& v) -> ColumnStats {
            ColumnStats s;
            s.min_f = *std::min_element(v.begin(), v.end());
            s.max_f = *std::max_element(v.begin(), v.end());
            return s;
        };

        // Write columns with stats (use ColumnID enum for correct indices)
        write_column(ColumnID::READ_IDX, buf_read_idx.data(), sizeof(uint32_t), buf_read_idx.size());
        rg.columns[static_cast<size_t>(ColumnID::READ_IDX)].stats = compute_stats_u32(buf_read_idx);

        write_column(ColumnID::REF_IDX, buf_ref_idx.data(), sizeof(uint32_t), buf_ref_idx.size());
        rg.columns[static_cast<size_t>(ColumnID::REF_IDX)].stats = compute_stats_u32(buf_ref_idx);

        write_column(ColumnID::BIT_SCORE, buf_bit_score.data(), sizeof(float), buf_bit_score.size());
        rg.columns[static_cast<size_t>(ColumnID::BIT_SCORE)].stats = compute_stats_f32(buf_bit_score);

        write_column(ColumnID::DAMAGE_SCORE, buf_damage_score.data(), sizeof(float), buf_damage_score.size());
        rg.columns[static_cast<size_t>(ColumnID::DAMAGE_SCORE)].stats = compute_stats_f32(buf_damage_score);

        write_column(ColumnID::EVALUE_LOG10, buf_evalue_log10.data(), sizeof(float), buf_evalue_log10.size());
        rg.columns[static_cast<size_t>(ColumnID::EVALUE_LOG10)].stats = compute_stats_f32(buf_evalue_log10);

        write_column(ColumnID::DMG_K, buf_dmg_k.data(), sizeof(uint16_t), buf_dmg_k.size());
        rg.columns[static_cast<size_t>(ColumnID::DMG_K)].stats = compute_stats_u16(buf_dmg_k);
        write_column(ColumnID::DMG_M, buf_dmg_m.data(), sizeof(uint16_t), buf_dmg_m.size());
        rg.columns[static_cast<size_t>(ColumnID::DMG_M)].stats = compute_stats_u16(buf_dmg_m);
        write_column(ColumnID::DMG_LL_A, buf_dmg_ll_a.data(), sizeof(float), buf_dmg_ll_a.size());
        rg.columns[static_cast<size_t>(ColumnID::DMG_LL_A)].stats = compute_stats_f32(buf_dmg_ll_a);
        write_column(ColumnID::DMG_LL_M, buf_dmg_ll_m.data(), sizeof(float), buf_dmg_ll_m.size());
        rg.columns[static_cast<size_t>(ColumnID::DMG_LL_M)].stats = compute_stats_f32(buf_dmg_ll_m);

        write_column(ColumnID::IDENTITY_Q, buf_identity_q.data(), sizeof(uint16_t), buf_identity_q.size());
        rg.columns[static_cast<size_t>(ColumnID::IDENTITY_Q)].stats = compute_stats_u16(buf_identity_q);

        write_column(ColumnID::ALN_LEN, buf_aln_len.data(), sizeof(uint16_t), buf_aln_len.size());
        rg.columns[static_cast<size_t>(ColumnID::ALN_LEN)].stats = compute_stats_u16(buf_aln_len);

        write_column(ColumnID::QSTART, buf_qstart.data(), sizeof(uint16_t), buf_qstart.size());
        rg.columns[static_cast<size_t>(ColumnID::QSTART)].stats = compute_stats_u16(buf_qstart);

        write_column(ColumnID::QEND, buf_qend.data(), sizeof(uint16_t), buf_qend.size());
        rg.columns[static_cast<size_t>(ColumnID::QEND)].stats = compute_stats_u16(buf_qend);

        write_column(ColumnID::TSTART, buf_tstart.data(), sizeof(uint16_t), buf_tstart.size());
        rg.columns[static_cast<size_t>(ColumnID::TSTART)].stats = compute_stats_u16(buf_tstart);

        write_column(ColumnID::TEND, buf_tend.data(), sizeof(uint16_t), buf_tend.size());
        rg.columns[static_cast<size_t>(ColumnID::TEND)].stats = compute_stats_u16(buf_tend);

        write_column(ColumnID::TLEN, buf_tlen.data(), sizeof(uint16_t), buf_tlen.size());
        rg.columns[static_cast<size_t>(ColumnID::TLEN)].stats = compute_stats_u16(buf_tlen);

        write_column(ColumnID::QLEN, buf_qlen.data(), sizeof(uint16_t), buf_qlen.size());
        rg.columns[static_cast<size_t>(ColumnID::QLEN)].stats = compute_stats_u16(buf_qlen);

        write_column(ColumnID::MISMATCH, buf_mismatch.data(), sizeof(uint16_t), buf_mismatch.size());
        rg.columns[static_cast<size_t>(ColumnID::MISMATCH)].stats = compute_stats_u16(buf_mismatch);

        write_column(ColumnID::GAPOPEN, buf_gapopen.data(), sizeof(uint16_t), buf_gapopen.size());
        rg.columns[static_cast<size_t>(ColumnID::GAPOPEN)].stats = compute_stats_u16(buf_gapopen);

        // Write string columns (offsets + data)
        write_column(ColumnID::QALN, buf_qaln_offsets.data(), sizeof(uint32_t), buf_qaln_offsets.size());
        temp_out.write(buf_qaln_data.data(), buf_qaln_data.size());

        write_column(ColumnID::TALN, buf_taln_offsets.data(), sizeof(uint32_t), buf_taln_offsets.size());
        temp_out.write(buf_taln_data.data(), buf_taln_data.size());

        row_groups.push_back(rg);
        row_group_read_index.push_back(RowGroupReadIndex{
            buf_read_idx.front(),
            buf_read_idx.back(),
            emitted_rows,
            emitted_rows + rg.num_rows
        });
        emitted_rows += rg.num_rows;

        // Clear buffers
        buf_read_idx.clear();
        buf_ref_idx.clear();
        buf_bit_score.clear();
        buf_damage_score.clear();
        buf_dmg_k.clear();
        buf_dmg_m.clear();
        buf_dmg_ll_a.clear();
        buf_dmg_ll_m.clear();
        buf_evalue_log10.clear();
        buf_identity_q.clear();
        buf_aln_len.clear();
        buf_qstart.clear();
        buf_qend.clear();
        buf_tstart.clear();
        buf_tend.clear();
        buf_tlen.clear();
        buf_qlen.clear();
        buf_mismatch.clear();
        buf_gapopen.clear();
        buf_qaln_offsets.clear();
        buf_qaln_offsets.push_back(0);
        buf_qaln_data.clear();
        buf_taln_offsets.clear();
        buf_taln_offsets.push_back(0);
        buf_taln_data.clear();
    }
};

ColumnarIndexWriter::ColumnarIndexWriter(const std::string& path)
    : impl_(std::make_unique<Impl>(path)) {}

ColumnarIndexWriter::~ColumnarIndexWriter() = default;

uint32_t ColumnarIndexWriter::add_ref(std::string_view name) {
    std::string key(name);
    auto it = impl_->ref_name_index.find(key);
    if (it != impl_->ref_name_index.end()) {
        return it->second;
    }
    uint32_t idx = static_cast<uint32_t>(impl_->ref_name_offsets.size());
    impl_->ref_name_offsets.push_back(static_cast<uint32_t>(impl_->ref_names.size()));
    impl_->ref_names.insert(impl_->ref_names.end(), name.begin(), name.end());
    impl_->ref_names.push_back('\0');
    impl_->ref_name_index[key] = idx;
    return idx;
}

uint32_t ColumnarIndexWriter::add_read(std::string_view name) {
    std::string key(name);
    auto it = impl_->read_name_index.find(key);
    if (it != impl_->read_name_index.end()) {
        return it->second;
    }
    uint32_t idx = static_cast<uint32_t>(impl_->read_name_offsets.size());
    impl_->read_name_offsets.push_back(static_cast<uint32_t>(impl_->read_names.size()));
    impl_->read_names.insert(impl_->read_names.end(), name.begin(), name.end());
    impl_->read_names.push_back('\0');
    impl_->read_name_index[key] = idx;

    // Update CSR
    while (impl_->csr_offsets.size() <= idx + 1) {
        impl_->csr_offsets.push_back(impl_->total_alignments);
    }
    return idx;
}

void ColumnarIndexWriter::add_alignment(const AlignmentRecord& rec) {
    // Track global sort invariant required by DART fast-path consumers.
    if (!impl_->has_prev_sort_key) {
        impl_->has_prev_sort_key = true;
    } else if (impl_->sorted_by_read_then_score) {
        if (rec.read_idx < impl_->prev_read_idx ||
            (rec.read_idx == impl_->prev_read_idx && rec.bit_score > impl_->prev_bit_score)) {
            impl_->sorted_by_read_then_score = false;
        }
    }
    impl_->prev_read_idx = rec.read_idx;
    impl_->prev_bit_score = rec.bit_score;

    // Track CSR
    if (rec.read_idx != impl_->current_read_idx) {
        impl_->current_read_idx = rec.read_idx;
        while (impl_->csr_offsets.size() <= rec.read_idx + 1) {
            impl_->csr_offsets.push_back(impl_->total_alignments);
        }
    }
    impl_->csr_offsets[rec.read_idx + 1] = impl_->total_alignments + 1;
    impl_->total_alignments++;

    // Add to current row group buffer
    impl_->buf_read_idx.push_back(rec.read_idx);
    impl_->buf_ref_idx.push_back(rec.ref_idx);
    impl_->buf_bit_score.push_back(rec.bit_score);
    impl_->buf_damage_score.push_back(rec.damage_score);
    impl_->buf_dmg_k.push_back(rec.dmg_k);
    impl_->buf_dmg_m.push_back(rec.dmg_m);
    impl_->buf_dmg_ll_a.push_back(rec.dmg_ll_a);
    impl_->buf_dmg_ll_m.push_back(rec.dmg_ll_m);
    impl_->buf_evalue_log10.push_back(rec.evalue_log10);
    impl_->buf_identity_q.push_back(static_cast<uint16_t>(std::clamp(rec.identity, 0.0f, 1.0f) * 65535.0f));
    impl_->buf_aln_len.push_back(rec.aln_len);
    impl_->buf_qstart.push_back(rec.qstart);
    impl_->buf_qend.push_back(rec.qend);
    impl_->buf_tstart.push_back(rec.tstart);
    impl_->buf_tend.push_back(rec.tend);
    impl_->buf_tlen.push_back(rec.tlen);
    impl_->buf_qlen.push_back(rec.qlen);
    impl_->buf_mismatch.push_back(rec.mismatch);
    impl_->buf_gapopen.push_back(rec.gapopen);

    // Add alignment strings
    if (!rec.qaln.empty() || !rec.taln.empty()) {
        impl_->saw_alignment_strings = true;
    }
    impl_->buf_qaln_data.insert(impl_->buf_qaln_data.end(), rec.qaln.begin(), rec.qaln.end());
    impl_->buf_qaln_offsets.push_back(static_cast<uint32_t>(impl_->buf_qaln_data.size()));

    impl_->buf_taln_data.insert(impl_->buf_taln_data.end(), rec.taln.begin(), rec.taln.end());
    impl_->buf_taln_offsets.push_back(static_cast<uint32_t>(impl_->buf_taln_data.size()));

    // Flush row group at read boundaries to avoid splitting a read's alignments
    // across groups, which would corrupt per-read normalization in streaming EM.
    // If a single read has more than 2*ROW_GROUP_SIZE alignments (pathological),
    // we allow a mid-read split rather than accumulating unboundedly.
    if (impl_->buf_ref_idx.size() >= ROW_GROUP_SIZE) {
        const bool at_read_boundary = impl_->buf_read_idx.empty() ||
                                      (rec.read_idx != impl_->buf_read_idx.back());
        if (at_read_boundary || impl_->buf_ref_idx.size() >= 2 * ROW_GROUP_SIZE) {
            impl_->flush_row_group();
        }
    }
}

void ColumnarIndexWriter::finalize(float d_max, float lambda) {
    // Flush remaining data
    impl_->flush_row_group();
    if (impl_->temp_out.is_open()) {
        impl_->temp_out.close();
    }

    // Now write the final file
    std::ofstream out(impl_->path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Cannot create output file: " + impl_->path);
    }

    // Write provisional header
    EMIHeader header;
    header.magic = EMI_MAGIC;
    header.version = EMI_VERSION;
    header.num_alignments = impl_->total_alignments;
    header.num_reads = static_cast<uint32_t>(impl_->read_name_offsets.size());
    header.num_refs = static_cast<uint32_t>(impl_->ref_name_offsets.size());
    header.num_row_groups = static_cast<uint32_t>(impl_->row_groups.size());
    header.row_group_size = ROW_GROUP_SIZE;
    header.d_max = d_max;
    header.lambda = lambda;
    if (impl_->sorted_by_read_then_score) {
        header.flags |= EMI_FLAG_SORTED_BY_READ_THEN_SCORE;
    }
    if (impl_->saw_alignment_strings) {
        header.flags |= EMI_FLAG_HAS_ALIGNMENT_STRINGS;
    }

    out.write(reinterpret_cast<const char*>(&header), sizeof(header));

    // Write row group directory
    header.row_group_dir_offset = out.tellp();
    out.write(reinterpret_cast<const char*>(impl_->row_groups.data()),
              impl_->row_groups.size() * sizeof(RowGroupMeta));

    // Write string dictionaries
    header.string_dict_offset = out.tellp();
    // Read names: offsets then data
    uint32_t read_names_count = static_cast<uint32_t>(impl_->read_name_offsets.size());
    out.write(reinterpret_cast<const char*>(&read_names_count), sizeof(read_names_count));
    out.write(reinterpret_cast<const char*>(impl_->read_name_offsets.data()),
              impl_->read_name_offsets.size() * sizeof(uint32_t));
    uint32_t read_names_size = static_cast<uint32_t>(impl_->read_names.size());
    out.write(reinterpret_cast<const char*>(&read_names_size), sizeof(read_names_size));
    out.write(impl_->read_names.data(), impl_->read_names.size());

    // Ref names: offsets then data
    uint32_t ref_names_count = static_cast<uint32_t>(impl_->ref_name_offsets.size());
    out.write(reinterpret_cast<const char*>(&ref_names_count), sizeof(ref_names_count));
    out.write(reinterpret_cast<const char*>(impl_->ref_name_offsets.data()),
              impl_->ref_name_offsets.size() * sizeof(uint32_t));
    uint32_t ref_names_size = static_cast<uint32_t>(impl_->ref_names.size());
    out.write(reinterpret_cast<const char*>(&ref_names_size), sizeof(ref_names_size));
    out.write(impl_->ref_names.data(), impl_->ref_names.size());

    // Copy row group data from temp file
    uint64_t row_group_data_offset = out.tellp();
    if (!impl_->row_groups.empty()) {
        std::ifstream temp_in(impl_->temp_path, std::ios::binary);
        out << temp_in.rdbuf();
        temp_in.close();
        std::remove(impl_->temp_path.c_str());
    }

    // Adjust row group offsets
    for (auto& rg : impl_->row_groups) {
        for (auto& col : rg.columns) {
            col.offset += row_group_data_offset;
        }
    }

    // Write CSR offsets
    if (!impl_->row_group_read_index.empty()) {
        header.flags |= EMI_FLAG_HAS_ROW_GROUP_READ_INDEX;
        header.row_group_read_index_offset = out.tellp();
        header.row_group_read_index_count = static_cast<uint32_t>(impl_->row_group_read_index.size());
        header.row_group_read_index_entry_size = sizeof(RowGroupReadIndex);
        out.write(reinterpret_cast<const char*>(impl_->row_group_read_index.data()),
                  impl_->row_group_read_index.size() * sizeof(RowGroupReadIndex));
    }

    // Write CSR offsets
    header.csr_offsets_offset = out.tellp();
    out.write(reinterpret_cast<const char*>(impl_->csr_offsets.data()),
              impl_->csr_offsets.size() * sizeof(uint64_t));

    // Re-write header with correct offsets
    out.seekp(0);
    out.write(reinterpret_cast<const char*>(&header), sizeof(header));

    // Re-write row group directory with adjusted offsets
    out.seekp(header.row_group_dir_offset);
    out.write(reinterpret_cast<const char*>(impl_->row_groups.data()),
              impl_->row_groups.size() * sizeof(RowGroupMeta));

    out.close();
}

// ---------------------------------------------------------------------------
// ColumnarIndexReader implementation
// ---------------------------------------------------------------------------

struct ColumnarIndexReader::Impl {
    void* data = nullptr;
    size_t file_size = 0;
    int fd = -1;

    const EMIHeader* header = nullptr;
    const RowGroupMeta* row_groups = nullptr;
    const uint64_t* csr_offsets = nullptr;
    const RowGroupReadIndex* row_group_read_index = nullptr;

    // String dictionaries
    const uint32_t* read_name_offsets = nullptr;
    const char* read_names = nullptr;
    uint32_t num_read_names = 0;

    const uint32_t* ref_name_offsets = nullptr;
    const char* ref_names = nullptr;
    uint32_t num_ref_names = 0;

    // File offsets for pread()-based access to ref_name data (NFS-safe, bypasses mmap).
    size_t ref_name_offsets_file_off = 0;
    size_t ref_names_file_off = 0;
    uint32_t ref_names_bytes = 0;

    // Filters
    std::vector<FilterPredicate> filters;

    // Column selection
    bool load_all_columns = true;
    bool column_needed[static_cast<size_t>(ColumnID::NUM_COLUMNS)] = {};

    void close() {
        if (data && data != MAP_FAILED) {
            munmap(data, file_size);
        }
        if (fd >= 0) {
            ::close(fd);
        }
        data = nullptr;
        fd = -1;
    }
};

ColumnarIndexReader::ColumnarIndexReader(const std::string& path)
    : impl_(std::make_unique<Impl>())
{
    impl_->fd = open(path.c_str(), O_RDONLY);
    if (impl_->fd < 0) {
        return;  // Invalid
    }

    struct stat st;
    if (fstat(impl_->fd, &st) < 0) {
        impl_->close();
        return;
    }
    impl_->file_size = static_cast<size_t>(st.st_size);

    impl_->data = mmap(nullptr, impl_->file_size, PROT_READ, MAP_PRIVATE, impl_->fd, 0);
    if (impl_->data == MAP_FAILED) {
        impl_->close();
        return;
    }

    // MADV_SEQUENTIAL: aggressive readahead keeps Pass 1 / scoring passes fast
    // and allows per-row-group MADV_DONTNEED to reclaim pages efficiently.
    // Switched to MADV_RANDOM in flush_pages() before streaming_em to prevent
    // speculative readahead from accumulating 100+ GB during 100 EM iterations.
    madvise(impl_->data, impl_->file_size, MADV_SEQUENTIAL);

    // Parse header
    impl_->header = static_cast<const EMIHeader*>(impl_->data);
    if (impl_->header->magic != EMI_MAGIC) {
        impl_->close();
        return;
    }
    if (impl_->header->version != EMI_VERSION) {
        impl_->close();
        return;
    }

    const char* base = static_cast<const char*>(impl_->data);

    // Row group directory
    impl_->row_groups = reinterpret_cast<const RowGroupMeta*>(
        base + impl_->header->row_group_dir_offset);

    // String dictionaries
    const char* dict_ptr = base + impl_->header->string_dict_offset;

    impl_->num_read_names = *reinterpret_cast<const uint32_t*>(dict_ptr);
    dict_ptr += sizeof(uint32_t);
    impl_->read_name_offsets = reinterpret_cast<const uint32_t*>(dict_ptr);
    dict_ptr += impl_->num_read_names * sizeof(uint32_t);
    uint32_t read_names_size = *reinterpret_cast<const uint32_t*>(dict_ptr);
    dict_ptr += sizeof(uint32_t);
    impl_->read_names = dict_ptr;
    dict_ptr += read_names_size;

    impl_->num_ref_names = *reinterpret_cast<const uint32_t*>(dict_ptr);
    dict_ptr += sizeof(uint32_t);
    impl_->ref_name_offsets = reinterpret_cast<const uint32_t*>(dict_ptr);
    dict_ptr += impl_->num_ref_names * sizeof(uint32_t);
    uint32_t ref_names_size = *reinterpret_cast<const uint32_t*>(dict_ptr);
    dict_ptr += sizeof(uint32_t);
    impl_->ref_names = dict_ptr;

    // Record file offsets for pread()-based access (immune to NFS page eviction).
    impl_->ref_name_offsets_file_off = static_cast<size_t>(
        reinterpret_cast<const char*>(impl_->ref_name_offsets) - base);
    impl_->ref_names_file_off = static_cast<size_t>(impl_->ref_names - base);
    impl_->ref_names_bytes = ref_names_size;

    // CSR offsets
    impl_->csr_offsets = reinterpret_cast<const uint64_t*>(
        base + impl_->header->csr_offsets_offset);

    // Optional row-group read index (v3)
    if ((impl_->header->flags & EMI_FLAG_HAS_ROW_GROUP_READ_INDEX) != 0 &&
        impl_->header->row_group_read_index_offset > 0 &&
        impl_->header->row_group_read_index_count == impl_->header->num_row_groups &&
        impl_->header->row_group_read_index_entry_size == sizeof(RowGroupReadIndex)) {
        impl_->row_group_read_index = reinterpret_cast<const RowGroupReadIndex*>(
            base + impl_->header->row_group_read_index_offset);
    }

    // Default: load all columns
    std::fill(std::begin(impl_->column_needed), std::end(impl_->column_needed), true);
}

ColumnarIndexReader::~ColumnarIndexReader() {
    if (impl_) impl_->close();
}

ColumnarIndexReader::ColumnarIndexReader(ColumnarIndexReader&&) noexcept = default;
ColumnarIndexReader& ColumnarIndexReader::operator=(ColumnarIndexReader&&) noexcept = default;

bool ColumnarIndexReader::is_valid() const {
    return impl_ && impl_->data != nullptr;
}

uint64_t ColumnarIndexReader::num_alignments() const {
    return impl_->header->num_alignments;
}

uint32_t ColumnarIndexReader::num_reads() const {
    return impl_->header->num_reads;
}

uint32_t ColumnarIndexReader::num_refs() const {
    return impl_->header->num_refs;
}

uint32_t ColumnarIndexReader::num_row_groups() const {
    return impl_->header->num_row_groups;
}

uint32_t ColumnarIndexReader::version() const {
    return impl_->header->version;
}

bool ColumnarIndexReader::has_alignment_strings() const {
    return (impl_->header->flags & EMI_FLAG_HAS_ALIGNMENT_STRINGS) != 0;
}

bool ColumnarIndexReader::is_sorted_by_read_then_score() const {
    return (impl_->header->flags & EMI_FLAG_SORTED_BY_READ_THEN_SCORE) != 0;
}

bool ColumnarIndexReader::has_row_group_read_index() const {
    return impl_->row_group_read_index != nullptr;
}

float ColumnarIndexReader::d_max() const {
    return impl_->header->d_max;
}

float ColumnarIndexReader::lambda() const {
    return impl_->header->lambda;
}

std::string_view ColumnarIndexReader::read_name(uint32_t idx) const {
    if (idx >= impl_->num_read_names) return {};
    const char* p = impl_->read_names + impl_->read_name_offsets[idx];
    return {p, strlen(p)};
}

std::string_view ColumnarIndexReader::ref_name(uint32_t idx) const {
    if (idx >= impl_->num_ref_names) return {};
    const char* p = impl_->ref_names + impl_->ref_name_offsets[idx];
    return {p, strlen(p)};
}

// Helper: pread with retry loop for partial reads (common on NFS).
static bool pread_full(int fd, void* buf, size_t count, off_t offset) {
    char* p = static_cast<char*>(buf);
    size_t remaining = count;
    while (remaining > 0) {
        ssize_t n = pread(fd, p, remaining, offset);
        if (n <= 0) return false;  // EOF or error
        p += n;
        offset += n;
        remaining -= static_cast<size_t>(n);
    }
    return true;
}

// Read all ref names via pread() — bypasses the mmap page cache entirely.
// On NFS + Linux 4.18 with MAP_PRIVATE, the kernel can evict pages at any time
// under memory pressure; re-faulting those pages from NFS returns garbage.
// pread() reads directly from the file descriptor without touching the mmap.
std::vector<std::string> ColumnarIndexReader::ref_names_pread() const {
    const uint32_t n = impl_->num_ref_names;
    std::vector<uint32_t> off_buf(n);
    if (!pread_full(impl_->fd, off_buf.data(), n * sizeof(uint32_t),
                    static_cast<off_t>(impl_->ref_name_offsets_file_off)))
        return {};
    std::vector<char> blob(impl_->ref_names_bytes);
    if (!pread_full(impl_->fd, blob.data(), impl_->ref_names_bytes,
                    static_cast<off_t>(impl_->ref_names_file_off)))
        return {};
    std::vector<std::string> result(n);
    for (uint32_t i = 0; i < n; ++i)
        result[i] = blob.data() + off_buf[i];
    return result;
}

uint64_t ColumnarIndexReader::read_offset(uint32_t read_idx) const {
    return impl_->csr_offsets[read_idx];
}

uint32_t ColumnarIndexReader::read_degree(uint32_t read_idx) const {
    return static_cast<uint32_t>(impl_->csr_offsets[read_idx + 1] - impl_->csr_offsets[read_idx]);
}

RowGroupReadIndex ColumnarIndexReader::row_group_read_index(uint32_t row_group) const {
    if (row_group >= impl_->header->num_row_groups) return {};
    if (impl_->row_group_read_index) return impl_->row_group_read_index[row_group];

    const RowGroupMeta& rg = impl_->row_groups[row_group];
    RowGroupReadIndex fallback{};
    fallback.first_read_idx = rg.first_read_idx;
    fallback.last_read_idx = rg.first_read_idx;
    fallback.start_row = static_cast<uint64_t>(row_group) * impl_->header->row_group_size;
    fallback.end_row_exclusive = fallback.start_row + rg.num_rows;
    return fallback;
}

void ColumnarIndexReader::add_filter(FilterPredicate pred) {
    impl_->filters.push_back(pred);
}

void ColumnarIndexReader::clear_filters() {
    impl_->filters.clear();
}

void ColumnarIndexReader::set_columns(std::initializer_list<ColumnID> cols) {
    impl_->load_all_columns = false;
    std::fill(std::begin(impl_->column_needed), std::end(impl_->column_needed), false);
    for (auto col : cols) {
        impl_->column_needed[static_cast<size_t>(col)] = true;
    }
}

void ColumnarIndexReader::set_all_columns() {
    impl_->load_all_columns = true;
    std::fill(std::begin(impl_->column_needed), std::end(impl_->column_needed), true);
}

void ColumnarIndexReader::flush_pages() {
    if (!impl_ || !impl_->data || impl_->file_size == 0) return;
#ifdef __linux__
    // posix_fadvise DONTNEED on the fd releases NFS-backed pages from the OS
    // page cache (more reliable than madvise alone for NFS). Call before
    // streaming_em to clear ~50-80 GB accumulated by Pass 1 + scoring passes.
    if (impl_->fd >= 0)
        posix_fadvise(impl_->fd, 0, static_cast<off_t>(impl_->file_size),
                      POSIX_FADV_DONTNEED);
    // Remove from our mmap address space too
    madvise(impl_->data, impl_->file_size, MADV_DONTNEED);
    // Switch from MADV_SEQUENTIAL to MADV_RANDOM for the streaming_em phase:
    // MADV_SEQUENTIAL's aggressive readahead is good for single-pass scans
    // (Pass 1, scoring) but accumulates 100+ GB over 100 EM iterations even
    // with per-row-group DONTNEED. MADV_RANDOM disables speculative readahead;
    // streaming_em still gets its 2-row-group explicit WILLNEED from
    // prefetch_row_group(), so throughput is not materially affected.
    madvise(impl_->data, impl_->file_size, MADV_RANDOM);
#endif
}

void ColumnarIndexReader::parallel_scan(RowGroupCallback callback) const {
    const char* base = static_cast<const char*>(impl_->data);
    const uint32_t num_rg = impl_->header->num_row_groups;
    constexpr uint32_t kPrefetchDistance = 2;

    auto prefetch_row_group = [&](uint32_t rg_idx) {
        if (rg_idx >= num_rg) return;
        const RowGroupMeta& rg = impl_->row_groups[rg_idx];
        size_t min_off = std::numeric_limits<size_t>::max();
        size_t max_end = 0;
        bool any = false;

        for (size_t c = 0; c < static_cast<size_t>(ColumnID::NUM_COLUMNS); ++c) {
            if (!impl_->load_all_columns && !impl_->column_needed[c]) continue;
            const auto& chunk = rg.columns[c];
            const size_t off = static_cast<size_t>(chunk.offset);
            const size_t sz = static_cast<size_t>(chunk.compressed_size);
            if (sz == 0) continue;
            min_off = std::min(min_off, off);
            max_end = std::max(max_end, off + sz);
            any = true;
        }
        if (any && max_end > min_off) {
            madvise_willneed_range(base + min_off, max_end - min_off);
        }
    };

    std::array<std::vector<char>, static_cast<size_t>(ColumnID::NUM_COLUMNS)> decoded_columns;

    for (uint32_t rg_idx = 0; rg_idx < num_rg; ++rg_idx) {
        const RowGroupMeta& rg = impl_->row_groups[rg_idx];
        prefetch_row_group(rg_idx + kPrefetchDistance);

        // Check if row group can be skipped via predicate pushdown
        bool skip = false;
        for (const auto& filter : impl_->filters) {
            size_t col_idx = static_cast<size_t>(filter.column);
            if (col_idx < static_cast<size_t>(ColumnID::NUM_COLUMNS)) {
                if (filter.can_skip(rg.columns[col_idx].stats)) {
                    skip = true;
                    break;
                }
            }
        }
        if (skip) continue;

        // Load column data (direct mmap pointers for uncompressed; decoded buffers for compressed).
        auto get_col = [&](ColumnID col) -> const void* {
            size_t idx = static_cast<size_t>(col);
            if (!impl_->load_all_columns && !impl_->column_needed[idx]) {
                return nullptr;
            }
            const auto& chunk = rg.columns[idx];
            const char* src = base + chunk.offset;
            if (chunk.codec == Codec::NONE) {
                return src;
            }
            decompress_chunk_to_buffer(chunk, src, decoded_columns[idx]);
            return decoded_columns[idx].data();
        };

        callback(
            rg_idx,
            rg.num_rows,
            static_cast<const uint32_t*>(get_col(ColumnID::READ_IDX)),
            static_cast<const uint32_t*>(get_col(ColumnID::REF_IDX)),
            static_cast<const float*>(get_col(ColumnID::BIT_SCORE)),
            static_cast<const float*>(get_col(ColumnID::DAMAGE_SCORE)),
            static_cast<const float*>(get_col(ColumnID::EVALUE_LOG10)),
            static_cast<const uint16_t*>(get_col(ColumnID::DMG_K)),
            static_cast<const uint16_t*>(get_col(ColumnID::DMG_M)),
            static_cast<const float*>(get_col(ColumnID::DMG_LL_A)),
            static_cast<const float*>(get_col(ColumnID::DMG_LL_M)),
            static_cast<const uint16_t*>(get_col(ColumnID::IDENTITY_Q)),
            static_cast<const uint16_t*>(get_col(ColumnID::ALN_LEN)),
            static_cast<const uint16_t*>(get_col(ColumnID::QSTART)),
            static_cast<const uint16_t*>(get_col(ColumnID::QEND)),
            static_cast<const uint16_t*>(get_col(ColumnID::TSTART)),
            static_cast<const uint16_t*>(get_col(ColumnID::TEND)),
            static_cast<const uint16_t*>(get_col(ColumnID::TLEN)),
            static_cast<const uint16_t*>(get_col(ColumnID::QLEN)),
            static_cast<const uint16_t*>(get_col(ColumnID::MISMATCH)),
            static_cast<const uint16_t*>(get_col(ColumnID::GAPOPEN))
        );

        // Release mmap pages for this row group — the callback has consumed
        // the data (copied to heap or accumulated into O(refs) sums), so these
        // pages are no longer needed until the next EM iteration re-reads them
        // from disk. Without this, all pages accumulate in RSS over multiple
        // streaming EM iterations, easily exceeding the file size.
        {
            size_t rg_min = std::numeric_limits<size_t>::max();
            size_t rg_end = 0;
            for (size_t c = 0; c < static_cast<size_t>(ColumnID::NUM_COLUMNS); ++c) {
                const auto& chunk = rg.columns[c];
                if (chunk.compressed_size == 0) continue;
                const size_t off = static_cast<size_t>(chunk.offset);
                const size_t sz  = static_cast<size_t>(chunk.compressed_size);
                rg_min = std::min(rg_min, off);
                rg_end = std::max(rg_end, off + sz);
            }
            if (rg_end > rg_min) {
#ifdef __linux__
                posix_fadvise(impl_->fd, static_cast<off_t>(rg_min),
                              static_cast<off_t>(rg_end - rg_min), POSIX_FADV_DONTNEED);
#endif
                madvise_dontneed_range(base + rg_min, rg_end - rg_min);
            }
        }
    }
}

void ColumnarIndexReader::parallel_scan_selected(
    const std::vector<uint32_t>& row_group_indices,
    RowGroupCallback callback) const {
    if (row_group_indices.empty()) return;

    const char* base = static_cast<const char*>(impl_->data);
    const uint32_t num_rg = impl_->header->num_row_groups;
    constexpr uint32_t kPrefetchDistance = 2;

    auto prefetch_selected = [&](size_t list_idx) {
        const size_t pf_idx = list_idx + kPrefetchDistance;
        if (pf_idx >= row_group_indices.size()) return;
        const uint32_t rg_idx = row_group_indices[pf_idx];
        if (rg_idx >= num_rg) return;

        const RowGroupMeta& rg = impl_->row_groups[rg_idx];
        size_t min_off = std::numeric_limits<size_t>::max();
        size_t max_end = 0;
        bool any = false;

        for (size_t c = 0; c < static_cast<size_t>(ColumnID::NUM_COLUMNS); ++c) {
            if (!impl_->load_all_columns && !impl_->column_needed[c]) continue;
            const auto& chunk = rg.columns[c];
            const size_t off = static_cast<size_t>(chunk.offset);
            const size_t sz = static_cast<size_t>(chunk.compressed_size);
            if (sz == 0) continue;
            min_off = std::min(min_off, off);
            max_end = std::max(max_end, off + sz);
            any = true;
        }
        if (any && max_end > min_off) {
            madvise_willneed_range(base + min_off, max_end - min_off);
        }
    };

    // Sequential (no OMP): each row group is read, strings copied into best_hits, then
    // DONTNEED is called before the next group loads. With OMP parallel, all threads
    // simultaneously touch different parts of the NFS-backed mmap — MADV_DONTNEED is
    // not flushed between groups, so the entire file stays resident (~114 GB for large EMI).
    // Sequential processing keeps mmap RSS at ~one-row-group at a time (~4 MB).
    for (size_t idx = 0; idx < row_group_indices.size(); ++idx) {
        const uint32_t rg_idx = row_group_indices[idx];
        if (rg_idx >= num_rg) continue;
        const RowGroupMeta& rg = impl_->row_groups[rg_idx];
        prefetch_selected(idx);

        bool skip = false;
        for (const auto& filter : impl_->filters) {
            size_t col_idx = static_cast<size_t>(filter.column);
            if (col_idx < static_cast<size_t>(ColumnID::NUM_COLUMNS)) {
                if (filter.can_skip(rg.columns[col_idx].stats)) {
                    skip = true;
                    break;
                }
            }
        }
        if (skip) continue;

        std::array<std::vector<char>, static_cast<size_t>(ColumnID::NUM_COLUMNS)> decoded_columns;

        auto get_col = [&](ColumnID col) -> const void* {
            size_t cidx = static_cast<size_t>(col);
            if (!impl_->load_all_columns && !impl_->column_needed[cidx]) {
                return nullptr;
            }
            const auto& chunk = rg.columns[cidx];
            const char* src = base + chunk.offset;
            if (chunk.codec == Codec::NONE) {
                return src;
            }
            decompress_chunk_to_buffer(chunk, src, decoded_columns[cidx]);
            return decoded_columns[cidx].data();
        };

        callback(
            rg_idx,
            rg.num_rows,
            static_cast<const uint32_t*>(get_col(ColumnID::READ_IDX)),
            static_cast<const uint32_t*>(get_col(ColumnID::REF_IDX)),
            static_cast<const float*>(get_col(ColumnID::BIT_SCORE)),
            static_cast<const float*>(get_col(ColumnID::DAMAGE_SCORE)),
            static_cast<const float*>(get_col(ColumnID::EVALUE_LOG10)),
            static_cast<const uint16_t*>(get_col(ColumnID::DMG_K)),
            static_cast<const uint16_t*>(get_col(ColumnID::DMG_M)),
            static_cast<const float*>(get_col(ColumnID::DMG_LL_A)),
            static_cast<const float*>(get_col(ColumnID::DMG_LL_M)),
            static_cast<const uint16_t*>(get_col(ColumnID::IDENTITY_Q)),
            static_cast<const uint16_t*>(get_col(ColumnID::ALN_LEN)),
            static_cast<const uint16_t*>(get_col(ColumnID::QSTART)),
            static_cast<const uint16_t*>(get_col(ColumnID::QEND)),
            static_cast<const uint16_t*>(get_col(ColumnID::TSTART)),
            static_cast<const uint16_t*>(get_col(ColumnID::TEND)),
            static_cast<const uint16_t*>(get_col(ColumnID::TLEN)),
            static_cast<const uint16_t*>(get_col(ColumnID::QLEN)),
            static_cast<const uint16_t*>(get_col(ColumnID::MISMATCH)),
            static_cast<const uint16_t*>(get_col(ColumnID::GAPOPEN))
        );

        // Release pages for this row group after processing — mirrors the sequential scan.
        // Pass 2 reads alignment strings for selected row groups; without DONTNEED the entire
        // 114 GB EMI mmap stays resident, causing ~114 GB of avoidable RSS.
        {
            size_t rg_min_off = std::numeric_limits<size_t>::max();
            size_t rg_max_end = 0;
            for (size_t c = 0; c < static_cast<size_t>(ColumnID::NUM_COLUMNS); ++c) {
                const auto& chunk = rg.columns[c];
                if (chunk.compressed_size == 0) continue;
                const size_t off = static_cast<size_t>(chunk.offset);
                const size_t sz  = static_cast<size_t>(chunk.compressed_size);
                rg_min_off = std::min(rg_min_off, off);
                rg_max_end = std::max(rg_max_end, off + sz);
            }
            if (rg_max_end > rg_min_off) {
#ifdef __linux__
                posix_fadvise(impl_->fd, static_cast<off_t>(rg_min_off),
                              static_cast<off_t>(rg_max_end - rg_min_off), POSIX_FADV_DONTNEED);
#endif
                madvise_dontneed_range(base + rg_min_off, rg_max_end - rg_min_off);
            }
        }
    }
}

std::string_view ColumnarIndexReader::get_qaln(uint32_t row_group, uint32_t row_in_group) const {
    if (!has_alignment_strings()) return {};
    if (row_group >= impl_->header->num_row_groups) return {};
    const char* base = static_cast<const char*>(impl_->data);
    const RowGroupMeta& rg = impl_->row_groups[row_group];
    if (row_in_group >= rg.num_rows) return {};

    const auto& col = rg.columns[static_cast<size_t>(ColumnID::QALN)];
    const char* chunk_bytes = nullptr;
    if (col.codec == Codec::NONE) {
        chunk_bytes = base + col.offset;
    } else {
        struct QalnCache {
            const void* impl_ptr = nullptr;
            uint64_t chunk_offset = 0;
            std::vector<char> data;
        };
        thread_local QalnCache cache;
        if (cache.impl_ptr != impl_.get() || cache.chunk_offset != col.offset ||
            cache.data.size() != static_cast<size_t>(col.uncompressed_size)) {
            decompress_chunk_to_buffer(col, base + col.offset, cache.data);
            cache.impl_ptr = impl_.get();
            cache.chunk_offset = col.offset;
        }
        chunk_bytes = cache.data.data();
    }

    const size_t offsets_bytes = (static_cast<size_t>(rg.num_rows) + 1) * sizeof(uint32_t);
    if (static_cast<size_t>(col.uncompressed_size) < offsets_bytes) return {};

    const uint32_t* offsets = reinterpret_cast<const uint32_t*>(chunk_bytes);
    const char* data = reinterpret_cast<const char*>(offsets + rg.num_rows + 1);

    uint32_t start = offsets[row_in_group];
    uint32_t end = offsets[row_in_group + 1];
    return {data + start, end - start};
}

std::string_view ColumnarIndexReader::get_taln(uint32_t row_group, uint32_t row_in_group) const {
    if (!has_alignment_strings()) return {};
    if (row_group >= impl_->header->num_row_groups) return {};
    const char* base = static_cast<const char*>(impl_->data);
    const RowGroupMeta& rg = impl_->row_groups[row_group];
    if (row_in_group >= rg.num_rows) return {};

    const auto& col = rg.columns[static_cast<size_t>(ColumnID::TALN)];
    const char* chunk_bytes = nullptr;
    if (col.codec == Codec::NONE) {
        chunk_bytes = base + col.offset;
    } else {
        struct TalnCache {
            const void* impl_ptr = nullptr;
            uint64_t chunk_offset = 0;
            std::vector<char> data;
        };
        thread_local TalnCache cache;
        if (cache.impl_ptr != impl_.get() || cache.chunk_offset != col.offset ||
            cache.data.size() != static_cast<size_t>(col.uncompressed_size)) {
            decompress_chunk_to_buffer(col, base + col.offset, cache.data);
            cache.impl_ptr = impl_.get();
            cache.chunk_offset = col.offset;
        }
        chunk_bytes = cache.data.data();
    }

    const size_t offsets_bytes = (static_cast<size_t>(rg.num_rows) + 1) * sizeof(uint32_t);
    if (static_cast<size_t>(col.uncompressed_size) < offsets_bytes) return {};

    const uint32_t* offsets = reinterpret_cast<const uint32_t*>(chunk_bytes);
    const char* data = reinterpret_cast<const char*>(offsets + rg.num_rows + 1);

    uint32_t start = offsets[row_in_group];
    uint32_t end = offsets[row_in_group + 1];
    return {data + start, end - start};
}

// ---------------------------------------------------------------------------
// EM solver on columnar index - proper per-read iteration with CSR
// ---------------------------------------------------------------------------

static double log_sum_exp(const double* vals, size_t n) {
    if (n == 0) return -std::numeric_limits<double>::infinity();
    if (n == 1) return vals[0];

    double max_val = vals[0];
    for (size_t i = 1; i < n; ++i) {
        if (vals[i] > max_val) max_val = vals[i];
    }

    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double diff = vals[i] - max_val;
        sum += (diff > -100.0) ? std::exp(diff) : 0.0;
    }
    return max_val + std::log(sum);
}

// Load all numeric columns into flat arrays for fast iteration
// v2: Now includes read_idx for correct per-read grouping
// After loading, sorts by (read_idx, bit_score DESC) and builds correct offsets
struct FlatColumnarData {
    std::vector<uint32_t> read_idx;   // v2: read index per alignment
    std::vector<uint32_t> ref_idx;
    std::vector<float> bit_score;
    std::vector<float> damage_score;
    std::vector<float> dmg_ll_a;
    std::vector<float> dmg_ll_m;
    std::vector<uint64_t> read_offsets;  // CSR offsets built from sorted read_idx

    void load(ColumnarIndexReader& reader) {
        const uint64_t n = reader.num_alignments();
        const uint32_t num_reads = reader.num_reads();
        read_idx.resize(n);
        ref_idx.resize(n);
        bit_score.resize(n);
        damage_score.resize(n);
        dmg_ll_a.resize(n);
        dmg_ll_m.resize(n);

        reader.parallel_scan([&](
            uint32_t rg_idx, uint32_t num_rows,
            const uint32_t* rd_idx, const uint32_t* rf_idx,
            const float* b_score, const float* d_score, const float*,
            const uint16_t*, const uint16_t*,
            const float* ll_a, const float* ll_m,
            const uint16_t*, const uint16_t*,
            const uint16_t*, const uint16_t*,
            const uint16_t*, const uint16_t*,
            const uint16_t*, const uint16_t*,
            const uint16_t*, const uint16_t*)
        {
            // Copy to flat arrays
            uint64_t base = static_cast<uint64_t>(rg_idx) * ROW_GROUP_SIZE;
            for (uint32_t i = 0; i < num_rows; ++i) {
                read_idx[base + i] = rd_idx[i];
                ref_idx[base + i] = rf_idx[i];
                bit_score[base + i] = b_score[i];
                damage_score[base + i] = d_score[i];
                dmg_ll_a[base + i] = ll_a ? ll_a[i] : 0.0f;
                dmg_ll_m[base + i] = ll_m ? ll_m[i] : 0.0f;
            }
        });
        // Fast path for v3+ files that are guaranteed sorted by (read_idx, bit_score DESC).
        bool already_sorted = reader.is_sorted_by_read_then_score();
        if (already_sorted && n > 1) {
            for (uint64_t i = 1; i < n; ++i) {
                if (read_idx[i] < read_idx[i - 1] ||
                    (read_idx[i] == read_idx[i - 1] && bit_score[i] > bit_score[i - 1])) {
                    already_sorted = false;
                    break;
                }
            }
        }

        if (!already_sorted) {
            // Sort by (read_idx, bit_score DESC) so best hit per read comes first.
            std::vector<uint64_t> perm(n);
            std::iota(perm.begin(), perm.end(), 0);
            std::sort(perm.begin(), perm.end(), [&](uint64_t a, uint64_t b) {
                if (read_idx[a] != read_idx[b]) return read_idx[a] < read_idx[b];
                return bit_score[a] > bit_score[b];  // Descending by score
            });

            // Apply permutation to all columns.
            auto apply_perm = [&](auto& vec) {
                using T = typename std::remove_reference<decltype(vec)>::type::value_type;
                std::vector<T> temp(vec.size());
                for (uint64_t i = 0; i < n; ++i) {
                    temp[i] = vec[perm[i]];
                }
                vec = std::move(temp);
            };
            apply_perm(read_idx);
            apply_perm(ref_idx);
            apply_perm(bit_score);
            apply_perm(damage_score);
            apply_perm(dmg_ll_a);
            apply_perm(dmg_ll_m);

            // Build CSR-style offsets from sorted read_idx.
            read_offsets.resize(num_reads + 1, 0);
            for (uint64_t i = 0; i < n; ++i) {
                if (read_idx[i] < num_reads) {
                    read_offsets[read_idx[i] + 1]++;
                }
            }
            for (uint32_t r = 1; r <= num_reads; ++r) {
                read_offsets[r] += read_offsets[r - 1];
            }
            return;
        }

        // Reuse on-disk CSR directly for already-sorted files.
        read_offsets.resize(num_reads + 1, 0);
        for (uint32_t r = 0; r <= num_reads; ++r) {
            read_offsets[r] = reader.read_offset(r);
        }
    }

    // Get read offset (replaces reader.read_offset())
    uint64_t get_read_offset(uint32_t r) const {
        return read_offsets[r];
    }
};

ColumnarEMResult em_solve_columnar(
    ColumnarIndexReader& reader,
    double lambda_b,
    uint32_t max_iters,
    double tol,
    bool use_damage,
    double alpha_prior,
    const std::vector<double>& initial_weights)
{
    const uint32_t num_refs = reader.num_refs();
    const uint64_t num_alns = reader.num_alignments();
    const uint32_t num_reads = reader.num_reads();

    // Load flat data for fast iteration
    FlatColumnarData data;
    data.load(reader);

    ColumnarEMResult result;
    if (!initial_weights.empty() && initial_weights.size() == num_refs) {
        result.weights = initial_weights;
        double sum = 0.0;
        for (double w : result.weights) sum += w;
        if (sum > 0.0) for (double& w : result.weights) w /= sum;
    } else {
        result.weights.assign(num_refs, 1.0 / static_cast<double>(num_refs));
    }
    result.gamma.resize(num_alns, 0.0);
    if (use_damage) {
        result.gamma_ancient.resize(num_alns, 0.0);
    }

    const double eps = 1e-15;
    const double min_weight = 1e-10;

    // SQUAREM buffers
    std::vector<double> w0(num_refs), w1(num_refs), w2(num_refs), r_buf(num_refs), v_buf(num_refs);

    double prev_ll = -std::numeric_limits<double>::infinity();

    for (uint32_t iter = 0; iter < max_iters; ++iter) {
        std::copy(result.weights.begin(), result.weights.end(), w0.begin());

        // --- E-step 1 ---
        double ll1 = 0.0;
        // Initialize with Dirichlet prior pseudocounts (prevents rich-get-richer)
        std::fill(w1.begin(), w1.end(), alpha_prior);

        #pragma omp parallel reduction(+:ll1)
        {
            std::vector<double> log_scores;
            std::vector<double> w_local(num_refs, 0.0);

            #pragma omp for schedule(dynamic, 1024)
            for (uint32_t r = 0; r < num_reads; ++r) {
                const uint64_t start = data.get_read_offset(r);
                const uint64_t end = data.get_read_offset(r + 1);
                const uint32_t deg = static_cast<uint32_t>(end - start);
                if (deg == 0) continue;

                // Unique mapper fast path
                if (deg == 1) {
                    const uint32_t ref = data.ref_idx[start];
                    result.gamma[start] = 1.0;
                    w_local[ref] += 1.0;
                    ll1 += std::log(std::max(w0[ref], min_weight))
                         + lambda_b * static_cast<double>(data.bit_score[start]);
                    continue;
                }

                log_scores.resize(deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    const uint32_t ref = data.ref_idx[start + j];
                    log_scores[j] = std::log(std::max(w0[ref], min_weight))
                                  + lambda_b * static_cast<double>(data.bit_score[start + j]);
                }

                double lse = log_sum_exp(log_scores.data(), deg);
                ll1 += lse;

                for (uint32_t j = 0; j < deg; ++j) {
                    double diff = log_scores[j] - lse;
                    double g = (diff > -100.0) ? std::exp(diff) : 0.0;
                    result.gamma[start + j] = g;
                    w_local[data.ref_idx[start + j]] += g;
                }
            }

            #pragma omp critical
            {
                for (uint32_t t = 0; t < num_refs; ++t) {
                    w1[t] += w_local[t];
                }
            }
        }

        // M-step 1: normalize with Dirichlet posterior mean
        // E[w_t] = (α + n_t) / (T*α + N) where N = num_reads, T = num_refs
        const double denom1 = static_cast<double>(num_reads) + num_refs * alpha_prior;
        double sum1 = 0.0;
        for (uint32_t t = 0; t < num_refs; ++t) {
            w1[t] = std::max(w1[t] / denom1, min_weight);
            sum1 += w1[t];
        }
        for (uint32_t t = 0; t < num_refs; ++t) w1[t] /= sum1;

        // --- E-step 2 ---
        // Initialize with Dirichlet prior pseudocounts (prevents rich-get-richer)
        std::fill(w2.begin(), w2.end(), alpha_prior);

        #pragma omp parallel
        {
            std::vector<double> log_scores;
            std::vector<double> w_local(num_refs, 0.0);

            #pragma omp for schedule(dynamic, 1024)
            for (uint32_t r = 0; r < num_reads; ++r) {
                const uint64_t start = data.get_read_offset(r);
                const uint64_t end = data.get_read_offset(r + 1);
                const uint32_t deg = static_cast<uint32_t>(end - start);
                if (deg == 0) continue;

                if (deg == 1) {
                    w_local[data.ref_idx[start]] += 1.0;
                    continue;
                }

                log_scores.resize(deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    const uint32_t ref = data.ref_idx[start + j];
                    log_scores[j] = std::log(std::max(w1[ref], min_weight))
                                  + lambda_b * static_cast<double>(data.bit_score[start + j]);
                }

                double lse = log_sum_exp(log_scores.data(), deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    double diff = log_scores[j] - lse;
                    double g = (diff > -100.0) ? std::exp(diff) : 0.0;
                    w_local[data.ref_idx[start + j]] += g;
                }
            }

            #pragma omp critical
            {
                for (uint32_t t = 0; t < num_refs; ++t) {
                    w2[t] += w_local[t];
                }
            }
        }

        // M-step 2: normalize with Dirichlet posterior mean
        const double denom2 = static_cast<double>(num_reads) + num_refs * alpha_prior;
        double sum2 = 0.0;
        for (uint32_t t = 0; t < num_refs; ++t) {
            w2[t] = std::max(w2[t] / denom2, min_weight);
            sum2 += w2[t];
        }
        for (uint32_t t = 0; t < num_refs; ++t) w2[t] /= sum2;

        // --- SQUAREM extrapolation ---
        double r_norm_sq = 0.0, v_norm_sq = 0.0;
        #pragma omp simd reduction(+:r_norm_sq, v_norm_sq)
        for (uint32_t t = 0; t < num_refs; ++t) {
            r_buf[t] = w1[t] - w0[t];
            v_buf[t] = (w2[t] - w1[t]) - r_buf[t];
            r_norm_sq += r_buf[t] * r_buf[t];
            v_norm_sq += v_buf[t] * v_buf[t];
        }

        double alpha = -std::sqrt(r_norm_sq / std::max(v_norm_sq, eps));
        alpha = std::clamp(alpha, -1.0, -0.01);

        #pragma omp simd
        for (uint32_t t = 0; t < num_refs; ++t) {
            w0[t] = w0[t] - 2.0 * alpha * r_buf[t] + alpha * alpha * v_buf[t];
            w0[t] = std::max(w0[t], eps);
        }

        // Re-normalize
        double sum0 = 0.0;
        for (uint32_t t = 0; t < num_refs; ++t) sum0 += w0[t];
        for (uint32_t t = 0; t < num_refs; ++t) w0[t] /= sum0;

        // --- Final E-step at accelerated point ---
        double ll = 0.0;

        #pragma omp parallel reduction(+:ll)
        {
            std::vector<double> log_scores;

            #pragma omp for schedule(dynamic, 1024)
            for (uint32_t r = 0; r < num_reads; ++r) {
                const uint64_t start = data.get_read_offset(r);
                const uint64_t end = data.get_read_offset(r + 1);
                const uint32_t deg = static_cast<uint32_t>(end - start);
                if (deg == 0) continue;

                if (deg == 1) {
                    const uint32_t ref = data.ref_idx[start];
                    result.gamma[start] = 1.0;
                    ll += std::log(std::max(w0[ref], min_weight))
                        + lambda_b * static_cast<double>(data.bit_score[start]);
                    continue;
                }

                log_scores.resize(deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    const uint32_t ref = data.ref_idx[start + j];
                    log_scores[j] = std::log(std::max(w0[ref], min_weight))
                                  + lambda_b * static_cast<double>(data.bit_score[start + j]);
                }

                double lse = log_sum_exp(log_scores.data(), deg);
                ll += lse;

                for (uint32_t j = 0; j < deg; ++j) {
                    double diff = log_scores[j] - lse;
                    result.gamma[start + j] = (diff > -100.0) ? std::exp(diff) : 0.0;
                }
            }
        }

        // Check for worsening (fall back to w2)
        if (ll < prev_ll && iter > 0) {
            std::copy(w2.begin(), w2.end(), result.weights.begin());
        } else {
            std::copy(w0.begin(), w0.end(), result.weights.begin());
        }

        result.log_likelihood = std::max(ll, prev_ll);
        result.iterations = iter + 1;

        // Convergence check
        if (iter > 0) {
            double rel_change = std::abs(ll - prev_ll) / (std::abs(prev_ll) + eps);
            if (rel_change < tol) break;
        }
        prev_ll = ll;
    }

    return result;
}

}  // namespace dart
