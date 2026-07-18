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

#include <cmath>
#include <string>
#include <string_view>
#include <vector>

namespace dart {

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
    explicit DamageIndexReader(const std::string& path, bool prefetch = true);

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
     * @brief v5+: DamageConfidence enum value of the library's validated pi.
     * @return 0 = DETECTED; 255 for pre-v5 indices (unknown).
     */
    uint8_t pi_state() const {
        return header_->version >= 5 ? header_->pi_state : 255;
    }

    /**
     * @brief v5+: library authenticated as DETECTED damage.
     *
     * The per-read pre-gate: libtaph read_ancient_posterior only emits a per-read
     * damage claim when pi.state == DETECTED (0). Libraries that do not validate
     * are surfaced (their proteins are dropped from the damage gate), never
     * silently scored via the AA-window fallback.
     */
    bool pi_detected() const {
        return header_->version >= 5 && header_->pi_state == 0;
    }

    /**
     * @brief v5+: the library's validated damaged-fraction point (pi.point).
     * @return the point estimate, or -1.0 for pre-v5 indices.
     */
    double pi_point() const {
        return header_->version >= 5 ? header_->pi_point : -1.0;
    }

    /**
     * @brief v5+: the library's pi is USABLE for damage annotation.
     *
     * Broader than pi_detected(): a library need not validate as DETECTED(0) to be
     * annotatable. The per-read path already degrades to the prior-free
     * damage_evidence_llr ranking when read_ancient_posterior refuses, so the ONLY
     * states that must fail every protein are the abstain-equivalents with no usable
     * pi point -- NOT_DETECTED(2), UNDETERMINED(3), BELOW_FLOOR(6), or a non-positive
     * point. DETECTED(0), TRACE(1), ANCIENT_CPG(4) and LOW_ABUNDANCE(5) all carry a
     * usable positive pi.point and are annotated (advisory when not DETECTED), so a
     * library whose deamination is authenticated but whose independent-molecule count
     * lands LOW_ABUNDANCE is no longer erased.
     */
    bool pi_usable() const {
        const double p = pi_point();
        if (header_->version < 5 || !std::isfinite(p) || p <= 0.0) return false;
        switch (header_->pi_state) {
            case 0:  // DETECTED
            case 1:  // TRACE
            case 4:  // ANCIENT_CPG
            case 5:  // LOW_ABUNDANCE
                return true;
            default: // NOT_DETECTED(2), UNDETERMINED(3), BELOW_FLOOR(6)
                return false;
        }
    }

    /**
     * @brief v6+: library-level empirical-Bayes damage prior (pi_damaged = δ/d_max).
     *
     * The damaged-read FRACTION. Unlike pi_point() (the state-gated pi.point), this survives
     * the abstain states (LOW_ABUNDANCE etc.) and is the prior for per-protein p_protein_damaged.
     * A library constant. Returns -1.0 for pre-v6 indices.
     */
    double damage_prior_point() const {
        return header_->version >= 6 ? static_cast<double>(header_->damage_prior_point) : -1.0;
    }

    /**
     * @brief v6+: DamageConfidence enum value attached to the damage prior.
     * @return the enum value, or 255 for pre-v6 indices.
     */
    uint8_t damage_prior_state() const {
        return header_->version >= 6 ? header_->damage_prior_state : 255;
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

    /**
     * @brief Synchronously warm up the page cache.
     *
     * Touches every page in the mmap'd region to force it into RAM.
     * Critical for NFS performance where MADV_WILLNEED may not be effective.
     * Call this before starting random lookups.
     *
     * @return Sum of bytes read (for preventing optimizer elimination)
     */
    size_t warmup_cache() const;

private:
    void* data_ = nullptr;          // mmap'd file data
    size_t file_size_ = 0;          // Total file size
    int fd_ = -1;                   // File descriptor

    // Pointers into mmap'd or decompressed region
    const AgdHeader* header_ = nullptr;
    const AgdBucket* buckets_ = nullptr;
    mutable const AgdRecord* records_ = nullptr;
    mutable const uint32_t* chain_ = nullptr;

    // v4: ZSTD-decompressed records+chain block (heap-allocated, set by warmup_cache)
    mutable std::vector<char> decompressed_block_;

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

}  // namespace dart
