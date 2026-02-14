#pragma once
/**
 * @file damage_index_reader.hpp
 * @brief Memory-mapped reader for .agd binary damage index files.
 *
 * Usage:
 *   DamageIndexReader reader("output.agd");
 *   if (auto rec = reader.find("read_name")) {
 *       // Use rec->codons_5prime, rec->nt_5prime, etc.
 *   }
 */

#include "damage_index.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace agp {

/**
 * @brief Memory-mapped reader for .agd damage index files.
 *
 * Provides O(1) lookup by read ID with zero-copy access to records.
 * Thread-safe for concurrent reads after construction.
 */
class DamageIndexReader {
public:
    /**
     * @brief Open and memory-map an .agd file.
     * @param path Path to .agd file
     * @throws std::runtime_error if file cannot be opened or is invalid
     */
    explicit DamageIndexReader(const std::string& path);

    /**
     * @brief Destructor unmaps file.
     */
    ~DamageIndexReader();

    // Non-copyable, movable
    DamageIndexReader(const DamageIndexReader&) = delete;
    DamageIndexReader& operator=(const DamageIndexReader&) = delete;
    DamageIndexReader(DamageIndexReader&& other) noexcept;
    DamageIndexReader& operator=(DamageIndexReader&& other) noexcept;

    /**
     * @brief Find record by read ID.
     * @param read_id Read ID (AGP suffix will be stripped automatically)
     * @return Pointer to record if found, nullptr otherwise
     *
     * The returned pointer is valid for the lifetime of this reader.
     */
    const AgdRecord* find(std::string_view read_id) const;

    /**
     * @brief Get sample-level d_max from header.
     */
    float d_max() const { return header_->d_max; }

    /**
     * @brief Get sample-level lambda from header.
     */
    float lambda() const { return header_->lambda; }

    /**
     * @brief Get library type (0=unknown, 1=ss, 2=ds).
     */
    uint8_t library_type() const { return header_->library_type; }

    /**
     * @brief Whether both Channel A and B agree on real damage (v3+).
     */
    bool damage_validated() const {
        return header_->version >= 3 && header_->damage_validated != 0;
    }

    /**
     * @brief Whether Channel A fired but Channel B didn't (false positive, v3+).
     */
    bool damage_artifact() const {
        return header_->version >= 3 && header_->damage_artifact != 0;
    }

    /**
     * @brief Whether sufficient data existed for Channel B (v3+).
     */
    bool channel_b_valid() const {
        return header_->version >= 3 && header_->channel_b_valid != 0;
    }

    /**
     * @brief Channel B stop-codon decay LLR (v3+, 0.0 for v2 files).
     */
    float stop_decay_llr() const {
        return header_->version >= 3 ? header_->stop_decay_llr : 0.0f;
    }

    /**
     * @brief Terminal T/(T+C) shift (v3+, 0.0 for v2 files).
     */
    float terminal_shift() const {
        return header_->version >= 3 ? header_->terminal_shift : 0.0f;
    }

    /**
     * @brief Whether this sample's damage signal is informative for scoring.
     *
     * v3+: requires damage_validated && !damage_artifact.
     * v2 fallback: uses d_max > 0 as a proxy.
     */
    bool damage_informative() const {
        if (header_->version < 3) return d_max() > 0.0f;
        return damage_validated() && !damage_artifact();
    }

    /**
     * @brief Get total number of records.
     */
    size_t record_count() const { return header_->num_records; }

    /**
     * @brief Check if file is valid and open.
     */
    bool is_valid() const { return data_ != nullptr; }

    /**
     * @brief Get record by index (for iteration).
     */
    const AgdRecord* get_record(size_t idx) const;

private:
    void* data_ = nullptr;          // mmap'd file data
    size_t file_size_ = 0;          // Total file size
    int fd_ = -1;                   // File descriptor

    // Pointers into mmap'd region
    const AgdHeader* header_ = nullptr;
    const AgdBucket* buckets_ = nullptr;
    const AgdRecord* records_ = nullptr;
    const uint32_t* chain_ = nullptr;  // Next pointers for hash collision chains

    void close();
};

/**
 * @brief Result of synonymous damage detection.
 */
struct SynonymousDamageResult {
    bool has_synonymous_damage = false;    // Any synonymous C→T/G→A detected
    int synonymous_5prime = 0;             // Count at 5' terminus
    int synonymous_3prime = 0;             // Count at 3' terminus
    int nonsynonymous_5prime = 0;          // Non-synonymous damage at 5'
    int nonsynonymous_3prime = 0;          // Non-synonymous damage at 3'

    // Detailed positions (codon index, wobble position)
    struct DamageSite {
        int codon_idx;      // 0-4 (terminal codon number)
        int nt_position;    // 0-2 within codon
        char observed_nt;   // T or A
        char expected_nt;   // C or G (before damage)
        bool is_synonymous; // True if AA unchanged
    };
    std::vector<DamageSite> sites;
};

/**
 * @brief Detect synonymous damage by comparing observed codons to expected.
 *
 * For C→T damage (5' end): checks if T at wobble positions could be from C
 * For G→A damage (3' end): checks if A at wobble positions could be from G
 *
 * @param rec Record from damage index
 * @param d_max Sample damage rate (for probability threshold)
 * @param lambda Decay rate
 * @return Detection result with counts and positions
 */
SynonymousDamageResult detect_synonymous_damage(
    const AgdRecord& rec,
    float d_max,
    float lambda);

}  // namespace agp
