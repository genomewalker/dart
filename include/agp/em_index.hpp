#pragma once
/**
 * @file em_index.hpp
 * @brief Legacy binary index format for standalone EM processing.
 *
 * Current CLI workflow uses:
 * 1. hits2emi: TSV + optional AGD -> columnar .emi
 * 2. damage-annotate: reads .emi and runs integrated EM reassignment
 *
 * File format (.emi):
 * - Header (72 bytes, see EMIHeader)
 * - String table (read names then ref names, null-separated)
 * - Alignment records (32 bytes each, see EMIAlignment)
 * - CSR offsets (8 bytes each, uint64_t)
 */

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace agp {

// Magic number: "AGPEMI01"
constexpr uint64_t EMI_MAGIC = 0x31494D455047414BULL;

struct EMIHeader {
    uint64_t magic = EMI_MAGIC;         // 8
    uint32_t version = 1;               // 4
    uint32_t flags = 0;                 // 4
    uint32_t num_reads = 0;             // 4
    uint32_t num_refs = 0;              // 4
    uint64_t num_alignments = 0;        // 8
    uint64_t string_table_offset = 0;   // 8
    uint32_t string_table_size = 0;     // 4 (max 4GB strings)
    uint32_t _pad1 = 0;                 // 4 (alignment padding)
    uint64_t alignments_offset = 0;     // 8
    uint64_t offsets_offset = 0;        // 8
    float d_max = 0.0f;                 // 4
    float lambda = 0.0f;                // 4
};                                      // Total: 72 bytes
static_assert(sizeof(EMIHeader) == 72);

// Alignment record - 32 bytes, lossless from TSV
struct EMIAlignment {
    uint32_t ref_idx;        // Reference index              4
    float bit_score;         // Alignment bit score (bits)   4
    float damage_score;      // From AGD lookup              4
    float evalue_log10;      // log10(evalue)                4
    uint16_t identity_q;     // fident * 65535               2
    uint16_t aln_len;        // Alignment length             2
    uint16_t tstart;         // Target start position        2
    uint16_t tend;           // Target end position          2
    uint16_t tlen;           // Target (reference) length    2
    uint16_t qlen;           // Query length                 2
    uint32_t _pad;           // Padding to 32 bytes          4
};                           //                         Total: 32
static_assert(sizeof(EMIAlignment) == 32);

/**
 * @brief Writer for creating .emi binary index files
 */
class EMIndexWriter {
public:
    explicit EMIndexWriter(const std::string& path);
    ~EMIndexWriter();

    // Add a read with its alignments (must be called in read order)
    void add_read(std::string_view name,
                  const std::vector<EMIAlignment>& alignments);

    // Add reference name (call before add_read for all refs)
    uint32_t add_ref(std::string_view name);

    // Finalize and write to disk
    void finalize(float d_max, float lambda);

private:
    std::string path_;
    std::vector<char> read_names_;
    std::vector<uint32_t> read_offsets_;
    std::vector<char> ref_names_;
    std::vector<uint32_t> ref_name_offsets_;
    std::vector<EMIAlignment> alignments_;
    uint32_t num_reads_ = 0;
    uint32_t num_refs_ = 0;
};

/**
 * @brief Memory-mapped reader for .emi files
 * Zero-copy access to alignment data for EM processing.
 */
class EMIndexReader {
public:
    explicit EMIndexReader(const std::string& path);
    ~EMIndexReader();

    EMIndexReader(const EMIndexReader&) = delete;
    EMIndexReader& operator=(const EMIndexReader&) = delete;

    // Header accessors
    uint32_t num_reads() const { return header_->num_reads; }
    uint32_t num_refs() const { return header_->num_refs; }
    uint64_t num_alignments() const { return header_->num_alignments; }
    float d_max() const { return header_->d_max; }
    float lambda() const { return header_->lambda; }

    // CSR access - alignments for read r are [offset(r), offset(r+1))
    uint64_t offset(uint32_t read_idx) const {
        return offsets_[read_idx];
    }

    const EMIAlignment* alignments() const { return alignments_; }
    const EMIAlignment* begin(uint32_t read_idx) const {
        return alignments_ + offsets_[read_idx];
    }
    const EMIAlignment* end(uint32_t read_idx) const {
        return alignments_ + offsets_[read_idx + 1];
    }
    uint32_t degree(uint32_t read_idx) const {
        return static_cast<uint32_t>(offsets_[read_idx + 1] - offsets_[read_idx]);
    }

    // Name lookup
    std::string_view read_name(uint32_t idx) const;
    std::string_view ref_name(uint32_t idx) const;

private:
    void* data_ = nullptr;
    size_t file_size_ = 0;
    int fd_ = -1;

    const EMIHeader* header_ = nullptr;
    const EMIAlignment* alignments_ = nullptr;
    const uint64_t* offsets_ = nullptr;
    const char* string_table_ = nullptr;

    // Pre-computed name offsets for O(1) lookup
    std::vector<uint32_t> read_name_offsets_;
    std::vector<uint32_t> ref_name_offsets_;
};

/**
 * @brief Run EM on pre-indexed data
 *
 * All data is mmap'd, no parsing needed.
 * Returns converged weights and responsibilities.
 */
struct EMResult {
    std::vector<double> weights;      // Per-reference abundance
    std::vector<double> gamma;        // Per-alignment responsibility
    std::vector<double> gamma_ancient;
    double pi = 0.1;                  // Ancient fraction
    double log_likelihood = 0.0;
    uint32_t iterations = 0;
};

EMResult em_solve_indexed(
    const EMIndexReader& index,
    double lambda_b = 3.0,
    uint32_t max_iters = 100,
    double tol = 1e-4,
    bool use_damage = true);

} // namespace agp
