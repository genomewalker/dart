// Columnar index implementation with parallel row group processing
//
// Key optimizations:
// - Column-oriented storage for cache efficiency
// - Row groups enable parallel decompression
// - Predicate pushdown via chunk statistics
// - Lazy loading of alignment strings

#include "agp/columnar_index.hpp"

#include <algorithm>
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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace agp {

// ---------------------------------------------------------------------------
// FilterPredicate implementation
// ---------------------------------------------------------------------------

bool FilterPredicate::can_skip(const ColumnStats& stats) const {
    // Check if entire chunk can be skipped based on stats
    switch (column) {
        case ColumnID::REF_IDX:
        case ColumnID::ALN_LEN:
        case ColumnID::TSTART:
        case ColumnID::TEND:
        case ColumnID::TLEN:
        case ColumnID::QLEN:
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
    std::vector<float> buf_evalue_log10;
    std::vector<uint16_t> buf_identity_q;
    std::vector<uint16_t> buf_aln_len;
    std::vector<uint16_t> buf_tstart;
    std::vector<uint16_t> buf_tend;
    std::vector<uint16_t> buf_tlen;
    std::vector<uint16_t> buf_qlen;
    std::vector<uint32_t> buf_qaln_offsets;
    std::vector<char> buf_qaln_data;
    std::vector<uint32_t> buf_taln_offsets;
    std::vector<char> buf_taln_data;

    // Row group metadata
    std::vector<RowGroupMeta> row_groups;

    // CSR offsets (per-read)
    std::vector<uint64_t> csr_offsets;
    uint32_t current_read_idx = UINT32_MAX;
    uint64_t total_alignments = 0;

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
        buf_evalue_log10.reserve(ROW_GROUP_SIZE);
        buf_identity_q.reserve(ROW_GROUP_SIZE);
        buf_aln_len.reserve(ROW_GROUP_SIZE);
        buf_tstart.reserve(ROW_GROUP_SIZE);
        buf_tend.reserve(ROW_GROUP_SIZE);
        buf_tlen.reserve(ROW_GROUP_SIZE);
        buf_qlen.reserve(ROW_GROUP_SIZE);
        buf_qaln_offsets.reserve(ROW_GROUP_SIZE + 1);
        buf_qaln_offsets.push_back(0);
        buf_taln_offsets.reserve(ROW_GROUP_SIZE + 1);
        buf_taln_offsets.push_back(0);
    }

    void flush_row_group() {
        if (buf_read_idx.empty()) return;

        RowGroupMeta rg;
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

        write_column(ColumnID::IDENTITY_Q, buf_identity_q.data(), sizeof(uint16_t), buf_identity_q.size());
        rg.columns[static_cast<size_t>(ColumnID::IDENTITY_Q)].stats = compute_stats_u16(buf_identity_q);

        write_column(ColumnID::ALN_LEN, buf_aln_len.data(), sizeof(uint16_t), buf_aln_len.size());
        rg.columns[static_cast<size_t>(ColumnID::ALN_LEN)].stats = compute_stats_u16(buf_aln_len);

        write_column(ColumnID::TSTART, buf_tstart.data(), sizeof(uint16_t), buf_tstart.size());
        rg.columns[static_cast<size_t>(ColumnID::TSTART)].stats = compute_stats_u16(buf_tstart);

        write_column(ColumnID::TEND, buf_tend.data(), sizeof(uint16_t), buf_tend.size());
        rg.columns[static_cast<size_t>(ColumnID::TEND)].stats = compute_stats_u16(buf_tend);

        write_column(ColumnID::TLEN, buf_tlen.data(), sizeof(uint16_t), buf_tlen.size());
        rg.columns[static_cast<size_t>(ColumnID::TLEN)].stats = compute_stats_u16(buf_tlen);

        write_column(ColumnID::QLEN, buf_qlen.data(), sizeof(uint16_t), buf_qlen.size());
        rg.columns[static_cast<size_t>(ColumnID::QLEN)].stats = compute_stats_u16(buf_qlen);

        // Write string columns (offsets + data)
        write_column(ColumnID::QALN, buf_qaln_offsets.data(), sizeof(uint32_t), buf_qaln_offsets.size());
        temp_out.write(buf_qaln_data.data(), buf_qaln_data.size());

        write_column(ColumnID::TALN, buf_taln_offsets.data(), sizeof(uint32_t), buf_taln_offsets.size());
        temp_out.write(buf_taln_data.data(), buf_taln_data.size());

        row_groups.push_back(rg);

        // Clear buffers
        buf_read_idx.clear();
        buf_ref_idx.clear();
        buf_bit_score.clear();
        buf_damage_score.clear();
        buf_evalue_log10.clear();
        buf_identity_q.clear();
        buf_aln_len.clear();
        buf_tstart.clear();
        buf_tend.clear();
        buf_tlen.clear();
        buf_qlen.clear();
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
    impl_->buf_evalue_log10.push_back(rec.evalue_log10);
    impl_->buf_identity_q.push_back(static_cast<uint16_t>(std::clamp(rec.identity, 0.0f, 1.0f) * 65535.0f));
    impl_->buf_aln_len.push_back(rec.aln_len);
    impl_->buf_tstart.push_back(rec.tstart);
    impl_->buf_tend.push_back(rec.tend);
    impl_->buf_tlen.push_back(rec.tlen);
    impl_->buf_qlen.push_back(rec.qlen);

    // Add alignment strings
    impl_->buf_qaln_data.insert(impl_->buf_qaln_data.end(), rec.qaln.begin(), rec.qaln.end());
    impl_->buf_qaln_offsets.push_back(static_cast<uint32_t>(impl_->buf_qaln_data.size()));

    impl_->buf_taln_data.insert(impl_->buf_taln_data.end(), rec.taln.begin(), rec.taln.end());
    impl_->buf_taln_offsets.push_back(static_cast<uint32_t>(impl_->buf_taln_data.size()));

    // Flush row group if full
    if (impl_->buf_ref_idx.size() >= ROW_GROUP_SIZE) {
        impl_->flush_row_group();
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

    // String dictionaries
    const uint32_t* read_name_offsets = nullptr;
    const char* read_names = nullptr;
    uint32_t num_read_names = 0;

    const uint32_t* ref_name_offsets = nullptr;
    const char* ref_names = nullptr;
    uint32_t num_ref_names = 0;

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

    madvise(impl_->data, impl_->file_size, MADV_SEQUENTIAL);

    // Parse header
    impl_->header = static_cast<const EMIHeader*>(impl_->data);
    if (impl_->header->magic != EMI_MAGIC) {
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

    // CSR offsets
    impl_->csr_offsets = reinterpret_cast<const uint64_t*>(
        base + impl_->header->csr_offsets_offset);

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

uint64_t ColumnarIndexReader::read_offset(uint32_t read_idx) const {
    return impl_->csr_offsets[read_idx];
}

uint32_t ColumnarIndexReader::read_degree(uint32_t read_idx) const {
    return static_cast<uint32_t>(impl_->csr_offsets[read_idx + 1] - impl_->csr_offsets[read_idx]);
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

void ColumnarIndexReader::parallel_scan(RowGroupCallback callback) const {
    const char* base = static_cast<const char*>(impl_->data);
    const uint32_t num_rg = impl_->header->num_row_groups;

    #pragma omp parallel for schedule(dynamic)
    for (uint32_t rg_idx = 0; rg_idx < num_rg; ++rg_idx) {
        const RowGroupMeta& rg = impl_->row_groups[rg_idx];

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

        // Load column data (direct pointers into mmap'd region)
        auto get_col = [&](ColumnID col) -> const void* {
            size_t idx = static_cast<size_t>(col);
            return base + rg.columns[idx].offset;
        };

        callback(
            rg_idx,
            rg.num_rows,
            static_cast<const uint32_t*>(get_col(ColumnID::READ_IDX)),
            static_cast<const uint32_t*>(get_col(ColumnID::REF_IDX)),
            static_cast<const float*>(get_col(ColumnID::BIT_SCORE)),
            static_cast<const float*>(get_col(ColumnID::DAMAGE_SCORE)),
            static_cast<const float*>(get_col(ColumnID::EVALUE_LOG10)),
            static_cast<const uint16_t*>(get_col(ColumnID::IDENTITY_Q)),
            static_cast<const uint16_t*>(get_col(ColumnID::ALN_LEN)),
            static_cast<const uint16_t*>(get_col(ColumnID::TSTART)),
            static_cast<const uint16_t*>(get_col(ColumnID::TEND)),
            static_cast<const uint16_t*>(get_col(ColumnID::TLEN)),
            static_cast<const uint16_t*>(get_col(ColumnID::QLEN))
        );
    }
}

std::string_view ColumnarIndexReader::get_qaln(uint32_t row_group, uint32_t row_in_group) const {
    const char* base = static_cast<const char*>(impl_->data);
    const RowGroupMeta& rg = impl_->row_groups[row_group];

    const auto& col = rg.columns[static_cast<size_t>(ColumnID::QALN)];
    const uint32_t* offsets = reinterpret_cast<const uint32_t*>(base + col.offset);
    const char* data = reinterpret_cast<const char*>(offsets + rg.num_rows + 1);

    uint32_t start = offsets[row_in_group];
    uint32_t end = offsets[row_in_group + 1];
    return {data + start, end - start};
}

std::string_view ColumnarIndexReader::get_taln(uint32_t row_group, uint32_t row_in_group) const {
    const char* base = static_cast<const char*>(impl_->data);
    const RowGroupMeta& rg = impl_->row_groups[row_group];

    const auto& col = rg.columns[static_cast<size_t>(ColumnID::TALN)];
    const uint32_t* offsets = reinterpret_cast<const uint32_t*>(base + col.offset);
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
    std::vector<uint64_t> read_offsets;  // CSR offsets built from sorted read_idx

    void load(ColumnarIndexReader& reader) {
        const uint64_t n = reader.num_alignments();
        const uint32_t num_reads = reader.num_reads();
        read_idx.resize(n);
        ref_idx.resize(n);
        bit_score.resize(n);
        damage_score.resize(n);

        reader.parallel_scan([&](
            uint32_t rg_idx, uint32_t num_rows,
            const uint32_t* rd_idx, const uint32_t* rf_idx,
            const float* b_score, const float* d_score, const float*,
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
            }
        });

        // Sort by (read_idx, bit_score DESC) so best hit per read comes first
        std::vector<uint64_t> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::sort(perm.begin(), perm.end(), [&](uint64_t a, uint64_t b) {
            if (read_idx[a] != read_idx[b]) return read_idx[a] < read_idx[b];
            return bit_score[a] > bit_score[b];  // Descending by score
        });

        // Apply permutation to all columns
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

        // Build CSR-style offsets from sorted read_idx
        read_offsets.resize(num_reads + 1, 0);
        for (uint64_t i = 0; i < n; ++i) {
            if (read_idx[i] < num_reads) {
                read_offsets[read_idx[i] + 1]++;
            }
        }
        // Prefix sum to get offsets
        for (uint32_t r = 1; r <= num_reads; ++r) {
            read_offsets[r] += read_offsets[r - 1];
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
    double alpha_prior)
{
    const uint32_t num_refs = reader.num_refs();
    const uint64_t num_alns = reader.num_alignments();
    const uint32_t num_reads = reader.num_reads();

    // Load flat data for fast iteration
    FlatColumnarData data;
    data.load(reader);

    ColumnarEMResult result;
    result.weights.assign(num_refs, 1.0 / static_cast<double>(num_refs));
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

}  // namespace agp
