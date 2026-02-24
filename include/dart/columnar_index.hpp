#pragma once
/**
 * @file columnar_index.hpp
 * @brief Columnar binary format for EM with predicate pushdown
 *
 * Design inspired by Parquet:
 * - Row groups for parallel processing (64K rows each)
 * - Per-column compression (LZ4 for strings, none for numerics initially)
 * - Chunk statistics for predicate pushdown (min/max per column per row group)
 * - Column pruning: only load columns needed for query
 *
 * File format (.emi):
 * ┌─────────────────────────────────────────┐
 * │ Header (128 bytes)                      │
 * ├─────────────────────────────────────────┤
 * │ Row Group Directory                     │
 * │   - Per row group: offsets, sizes, stats│
 * ├─────────────────────────────────────────┤
 * │ String Dictionary                       │
 * │   - Read names (deduplicated)           │
 * │   - Ref names (deduplicated)            │
 * ├─────────────────────────────────────────┤
 * │ Row Group 0                             │
 * │   ├─ Column: ref_idx (uint32[])         │
 * │   ├─ Column: bit_score (float[])        │
 * │   ├─ Column: damage_score (float[])     │
 * │   ├─ Column: evalue_log10 (float[])     │
 * │   ├─ Column: identity_q (uint16[])      │
 * │   ├─ Column: aln_len (uint16[])         │
 * │   ├─ Column: tstart (uint16[])          │
 * │   ├─ Column: tend (uint16[])            │
 * │   ├─ Column: tlen (uint16[])            │
 * │   ├─ Column: qlen (uint16[])            │
 * │   ├─ Column: qaln_offsets (uint32[])    │
 * │   ├─ Column: qaln_data (char[])         │
 * │   ├─ Column: taln_offsets (uint32[])    │
 * │   └─ Column: taln_data (char[])         │
 * ├─────────────────────────────────────────┤
 * │ Row Group 1...N                         │
 * ├─────────────────────────────────────────┤
 * │ CSR Offsets (read -> alignment range)   │
 * └─────────────────────────────────────────┘
 */

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <functional>

namespace dart {

// Magic: "AGPEMI02" for columnar EM index format (little-endian)
// v2 adds READ_IDX column to fix CSR/row-group incompatibility
// v3 adds AGP-specific execution metadata/flags for faster consumers
// v4 adds TSV parity columns: qstart/qend/mismatch/gapopen
// v5 adds per-alignment damage evidence for reassignment:
//    dmg_k, dmg_m, dmg_ll_a, dmg_ll_m
constexpr uint64_t EMI_MAGIC = 0x3230494D45504741ULL;  // "AGPEMI02"
constexpr uint64_t EMI_MAGIC_V1 = 0x3130494D45504741ULL;  // "AGPEMI01" (legacy)
constexpr uint32_t EMI_VERSION = 5;
constexpr uint32_t ROW_GROUP_SIZE = 65536;  // 64K rows per group

// Header flags (AGP execution hints/capabilities)
constexpr uint32_t EMI_FLAG_SORTED_BY_READ_THEN_SCORE = 1u << 0;
constexpr uint32_t EMI_FLAG_HAS_ALIGNMENT_STRINGS = 1u << 1;
constexpr uint32_t EMI_FLAG_HAS_ROW_GROUP_READ_INDEX = 1u << 2;

// Column IDs for selective loading
// v2: Added READ_IDX to store read index per alignment (fixes CSR/row-group issue)
enum class ColumnID : uint8_t {
    READ_IDX = 0,      // v2: read index per alignment (enables correct grouping)
    REF_IDX = 1,
    BIT_SCORE = 2,
    DAMAGE_SCORE = 3,
    EVALUE_LOG10 = 4,
    IDENTITY_Q = 5,
    ALN_LEN = 6,
    QSTART = 7,        // v4: query start (1-based from TSV, clamped to uint16)
    QEND = 8,          // v4: query end (1-based from TSV, clamped to uint16)
    TSTART = 9,
    TEND = 10,
    TLEN = 11,
    QLEN = 12,
    MISMATCH = 13,     // v4: mismatch count (clamped to uint16)
    GAPOPEN = 14,      // v4: gap-open count (clamped to uint16)
    QALN = 15,
    TALN = 16,
    DMG_K = 17,        // v5: damage-consistent hits among opportunities
    DMG_M = 18,        // v5: damage opportunities in aligned columns
    DMG_LL_A = 19,     // v5: log P(evidence | ancient)
    DMG_LL_M = 20,     // v5: log P(evidence | modern)
    NUM_COLUMNS = 21
};

// Compression codec per column
enum class Codec : uint8_t {
    NONE = 0,
    LZ4 = 1,
    ZSTD = 2,
    DELTA = 3,     // Delta encoding for sorted/sequential values
    DICT = 4,      // Dictionary encoding
};

// Compression level (0 = fastest, higher = better compression)
constexpr int ZSTD_COMPRESSION_LEVEL = 3;  // Good balance of speed/ratio

// Per-column statistics for predicate pushdown
struct ColumnStats {
    float min_f = 0;
    float max_f = 0;
    uint32_t min_u = 0;
    uint32_t max_u = 0;
    uint32_t null_count = 0;
};

// Per-column chunk metadata
struct ColumnChunk {
    uint64_t offset = 0;        // Offset in file
    uint32_t compressed_size = 0;
    uint32_t uncompressed_size = 0;
    Codec codec = Codec::NONE;
    ColumnStats stats;
};

// Row group metadata
struct RowGroupMeta {
    uint32_t num_rows = 0;
    uint32_t first_read_idx = 0;  // For CSR mapping
    ColumnChunk columns[static_cast<size_t>(ColumnID::NUM_COLUMNS)];
};

// Optional per-row-group read span index (v3 extension)
struct RowGroupReadIndex {
    uint32_t first_read_idx = 0;
    uint32_t last_read_idx = 0;
    uint64_t start_row = 0;          // Global start row (inclusive)
    uint64_t end_row_exclusive = 0;  // Global end row (exclusive)
};

// File header (128 bytes)
struct EMIHeader {
    uint64_t magic = EMI_MAGIC;
    uint32_t version = EMI_VERSION;
    uint32_t flags = 0;
    uint64_t num_alignments = 0;
    uint32_t num_reads = 0;
    uint32_t num_refs = 0;
    uint32_t num_row_groups = 0;
    uint32_t row_group_size = ROW_GROUP_SIZE;
    uint64_t row_group_dir_offset = 0;
    uint64_t string_dict_offset = 0;
    uint64_t csr_offsets_offset = 0;
    float d_max = 0.0f;
    float lambda = 0.3f;
    // v3 extension fields (kept within former reserved space)
    uint64_t row_group_read_index_offset = 0;
    uint32_t row_group_read_index_count = 0;
    uint32_t row_group_read_index_entry_size = 0;
    uint8_t _reserved[36];
};
static_assert(sizeof(EMIHeader) == 128);

// Filter predicate for pushdown
struct FilterPredicate {
    ColumnID column;
    enum class Op { GE, GT, LE, LT, EQ } op;
    float value_f = 0;
    uint32_t value_u = 0;

    bool can_skip(const ColumnStats& stats) const;
};

// Forward declarations
class ColumnarIndexWriter;
class ColumnarIndexReader;

/**
 * @brief Alignment record for writing (expanded form)
 */
struct AlignmentRecord {
    uint32_t read_idx;
    uint32_t ref_idx;
    float bit_score;
    float damage_score;   // per-read prior p_read
    uint16_t dmg_k = 0;
    uint16_t dmg_m = 0;
    float dmg_ll_a = 0.0f;
    float dmg_ll_m = 0.0f;
    float evalue_log10;
    float identity;      // 0-1, will be quantized
    uint16_t aln_len;
    uint16_t qstart;
    uint16_t qend;
    uint16_t tstart;
    uint16_t tend;
    uint16_t tlen;
    uint16_t qlen;
    uint16_t mismatch;
    uint16_t gapopen;
    std::string_view qaln;
    std::string_view taln;
};

/**
 * @brief Writer for columnar index files
 */
class ColumnarIndexWriter {
public:
    explicit ColumnarIndexWriter(const std::string& path);
    ~ColumnarIndexWriter();

    // Add a reference name, returns index
    uint32_t add_ref(std::string_view name);

    // Add a read name, returns index
    uint32_t add_read(std::string_view name);

    // Add an alignment record
    void add_alignment(const AlignmentRecord& rec);

    // Finalize and write to disk
    void finalize(float d_max, float lambda);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/**
 * @brief Reader for columnar index files with predicate pushdown
 */
class ColumnarIndexReader {
public:
    explicit ColumnarIndexReader(const std::string& path);
    ~ColumnarIndexReader();

    // Movable but not copyable
    ColumnarIndexReader(ColumnarIndexReader&&) noexcept;
    ColumnarIndexReader& operator=(ColumnarIndexReader&&) noexcept;
    ColumnarIndexReader(const ColumnarIndexReader&) = delete;
    ColumnarIndexReader& operator=(const ColumnarIndexReader&) = delete;

    bool is_valid() const;

    // Metadata access
    uint64_t num_alignments() const;
    uint32_t num_reads() const;
    uint32_t num_refs() const;
    uint32_t num_row_groups() const;
    uint32_t version() const;
    bool has_alignment_strings() const;
    bool is_sorted_by_read_then_score() const;
    bool has_row_group_read_index() const;
    float d_max() const;
    float lambda() const;

    // Name lookup
    std::string_view read_name(uint32_t idx) const;
    std::string_view ref_name(uint32_t idx) const;

    // CSR access for per-read iteration
    uint64_t read_offset(uint32_t read_idx) const;
    uint32_t read_degree(uint32_t read_idx) const;
    RowGroupReadIndex row_group_read_index(uint32_t row_group) const;

    // Set filters for predicate pushdown (call before iterating)
    void add_filter(FilterPredicate pred);
    void clear_filters();

    // Specify which columns to load (optimization)
    void set_columns(std::initializer_list<ColumnID> cols);
    void set_all_columns();

    // Release all mmap'd file pages from the OS page cache.
    // Call before streaming_em to flush pages accumulated by Pass 1 / scoring,
    // ensuring streaming_em starts from a clean RSS baseline.
    void flush_pages();

    // Row group iteration with parallel callback
    // Callback receives (row_group_idx, num_rows, column_data_ptrs)
    // v2: Added read_idx for correct per-read grouping (fixes CSR/row-group issue)
    using RowGroupCallback = std::function<void(
        uint32_t rg_idx,
        uint32_t num_rows,
        const uint32_t* read_idx,    // v2: read index per alignment
        const uint32_t* ref_idx,
        const float* bit_score,
        const float* damage_score,
        const float* evalue_log10,
        const uint16_t* dmg_k,
        const uint16_t* dmg_m,
        const float* dmg_ll_a,
        const float* dmg_ll_m,
        const uint16_t* identity_q,
        const uint16_t* aln_len,
        const uint16_t* qstart,
        const uint16_t* qend,
        const uint16_t* tstart,
        const uint16_t* tend,
        const uint16_t* tlen,
        const uint16_t* qlen,
        const uint16_t* mismatch,
        const uint16_t* gapopen
    )>;

    // Process row groups in parallel (OpenMP)
    void parallel_scan(RowGroupCallback callback) const;
    void parallel_scan_selected(
        const std::vector<uint32_t>& row_group_indices,
        RowGroupCallback callback) const;

    // Get alignment strings for specific row (lazy load)
    std::string_view get_qaln(uint32_t row_group, uint32_t row_in_group) const;
    std::string_view get_taln(uint32_t row_group, uint32_t row_in_group) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/**
 * @brief Result of EM solve on columnar index
 */
struct ColumnarEMResult {
    std::vector<double> weights;       // Per-reference abundance
    std::vector<double> gamma;         // Per-alignment responsibility
    std::vector<double> gamma_ancient; // Ancient responsibility (if damage-aware)
    double pi = 0.1;                   // Ancient fraction
    double log_likelihood = 0.0;
    uint32_t iterations = 0;
};

/**
 * @brief Run EM on columnar index with parallel row group processing
 *
 * @param reader Columnar index reader
 * @param lambda_b Temperature for bit-score scaling (default 3.0)
 * @param max_iters Maximum iterations (default 100)
 * @param tol Convergence tolerance (default 1e-4)
 * @param use_damage Include damage-aware gamma_ancient (default true)
 * @param alpha_prior Dirichlet pseudocount for "rich get richer" prevention.
 *                    Higher values = more regularization toward uniform.
 *                    - 0.0: No prior (pure MLE, can collapse to single reference)
 *                    - 1.0: Uniform prior (recommended default)
 *                    - >1: Strong regularization for sparse data
 */
ColumnarEMResult em_solve_columnar(
    ColumnarIndexReader& reader,
    double lambda_b = 3.0,
    uint32_t max_iters = 100,
    double tol = 1e-4,
    bool use_damage = true,
    double alpha_prior = 1.0,
    const std::vector<double>& initial_weights = {});

}  // namespace dart
