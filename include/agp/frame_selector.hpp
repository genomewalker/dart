#pragma once

#include "types.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <utility>

namespace agp {

// Forward declaration
class DamageModel;

/**
 * Sample-level damage profile computed from aggregate statistics
 */
struct SampleDamageProfile {
    // Position-specific base counts at 5' end (positions 0-14)
    // Using double to avoid float precision loss at >16M reads
    std::array<double, 15> t_freq_5prime = {};  // T count at each position
    std::array<double, 15> c_freq_5prime = {};  // C count (baseline)

    // Position-specific base counts at 3' end (positions 0-14 from end)
    std::array<double, 15> a_freq_3prime = {};  // A count at each position
    std::array<double, 15> g_freq_3prime = {};  // G count (baseline)

    // Middle-of-read baseline counts (undamaged)
    double baseline_t_freq = 0.0;
    double baseline_c_freq = 0.0;
    double baseline_a_freq = 0.0;
    double baseline_g_freq = 0.0;

    // Computed damage rates (excess over baseline)
    std::array<float, 15> damage_rate_5prime = {};  // C→T rate at each 5' position
    std::array<float, 15> damage_rate_3prime = {};  // G→A rate at each 3' position

    // Codon-position-aware damage tracking (positions 1,2,3 in codon)
    // At 5' end: T/(T+C) ratio by codon position
    std::array<float, 3> codon_pos_t_rate_5prime = {0.5f, 0.5f, 0.5f};
    // At 3' end: A/(A+G) ratio by codon position
    std::array<float, 3> codon_pos_a_rate_3prime = {0.5f, 0.5f, 0.5f};
    // Raw counts for aggregation
    std::array<size_t, 3> codon_pos_t_count_5prime = {};
    std::array<size_t, 3> codon_pos_c_count_5prime = {};
    std::array<size_t, 3> codon_pos_a_count_3prime = {};
    std::array<size_t, 3> codon_pos_g_count_3prime = {};

    // CpG context damage tracking
    float cpg_damage_rate = 0.0f;      // C→T rate in CpG context
    float non_cpg_damage_rate = 0.0f;  // C→T rate outside CpG
    size_t cpg_c_count = 0;            // C's in CpG context
    size_t cpg_t_count = 0;            // T's where C expected in CpG
    size_t non_cpg_c_count = 0;
    size_t non_cpg_t_count = 0;

    // Summary statistics
    float max_damage_5prime = 0.0f;  // Maximum C→T rate at position 0
    float max_damage_3prime = 0.0f;  // Maximum G→A rate at last position
    float sample_damage_prob = 0.0f; // Overall probability sample is ancient

    // Estimated decay constants (sample-specific)
    float lambda_5prime = 0.3f;  // Decay constant for 5' end (estimated from data)
    float lambda_3prime = 0.3f;  // Decay constant for 3' end (estimated from data)

    // Library type detection
    enum class LibraryType { UNKNOWN, DOUBLE_STRANDED, SINGLE_STRANDED };
    LibraryType library_type = LibraryType::UNKNOWN;

    size_t n_reads = 0;  // Number of reads used in computation

    bool is_valid() const { return n_reads >= 1000; }  // Need enough reads

    // Library type string for output
    const char* library_type_str() const {
        switch (library_type) {
            case LibraryType::DOUBLE_STRANDED: return "double-stranded";
            case LibraryType::SINGLE_STRANDED: return "single-stranded";
            default: return "unknown";
        }
    }
};

/**
 * Score for a single reading frame
 * Combines multiple signals to evaluate frame correctness
 */
struct FrameScore {
    int frame = 0;              // Frame offset (0, 1, 2)
    bool forward = true;        // Strand direction

    float codon_score = 0.0f;           // Codon usage bias score
    float stop_codon_penalty = 0.0f;    // Penalty for internal stops
    float aa_composition_score = 0.0f;  // Amino acid composition score
    float damage_consistency = 0.0f;    // Damage pattern consistency

    float total_score = 0.0f;           // Combined weighted score

    std::string protein;                // Translated protein sequence
    std::string corrected_protein;      // Bayesian-corrected protein (if damage-aware)
    std::string corrected_dna;          // Bayesian-corrected DNA (if damage-aware)
    size_t stops_corrected = 0;         // Number of stop codons corrected
    size_t aa_changes = 0;              // Number of AA changes from correction
};

/**
 * Frame selection for ancient DNA reads
 *
 * Combines multiple signals to select the correct reading frame:
 * 1. Codon usage bias
 * 2. Damage pattern consistency (C→T at 5', G→A at 3')
 * 3. Stop codon avoidance
 * 4. Amino acid composition
 */
class FrameSelector {
public:
    /**
     * Score a single reading frame
     *
     * @param seq DNA sequence
     * @param frame Frame offset (0, 1, or 2)
     * @param forward True for forward strand, false for reverse complement
     * @param damage Optional damage profile for damage-aware scoring
     * @return FrameScore with all scoring components
     */
    static FrameScore score_frame(
        const std::string& seq,
        int frame,
        bool forward,
        const DamageProfile* damage = nullptr);

    /**
     * Score all 6 reading frames
     *
     * @param seq DNA sequence
     * @param damage Optional damage profile
     * @return Vector of 6 FrameScores, sorted by total_score (highest first)
     */
    static std::vector<FrameScore> score_all_frames(
        const std::string& seq,
        const DamageProfile* damage = nullptr);

    /**
     * Select the best reading frame
     *
     * @param seq DNA sequence
     * @param damage Optional damage profile
     * @return FrameScore for the highest-scoring frame
     */
    static FrameScore select_best_frame(
        const std::string& seq,
        const DamageProfile* damage = nullptr);

    /**
     * Select the best reading frame from each strand
     * For short aDNA where strand cannot be determined
     *
     * @param seq DNA sequence
     * @param damage Optional damage profile
     * @return Pair of FrameScores: (best_forward, best_reverse)
     */
    static std::pair<FrameScore, FrameScore> select_best_per_strand(
        const std::string& seq,
        const DamageProfile* damage = nullptr);

    /**
     * Damage-aware frame selection with adjusted stop penalties
     *
     * Scores frames on UNCORRECTED sequences but adjusts stop codon penalties
     * based on damage probability:
     * 1. Score each frame normally (preserves natural codon patterns)
     * 2. For stop codons, calculate probability they're damage-induced
     * 3. Reduce penalty for likely damage stops (at 5' end, specific contexts)
     * 4. Select frame with best adjusted score
     *
     * This improves frame selection because:
     * - Correct frame with damage-induced stops → reduced penalties → better score
     * - Wrong frame with real stops → full penalties → worse score
     * - No sequence modification → preserves natural patterns for scoring
     *
     * @param seq DNA sequence
     * @param damage_model DamageModel with sample-level damage parameters
     * @return Pair of FrameScores: (best_forward, best_reverse)
     */
    static std::pair<FrameScore, FrameScore> select_best_per_strand_damage_aware(
        const std::string& seq,
        const DamageModel& damage_model);

    /**
     * Get reverse complement of a sequence
     */
    static std::string reverse_complement(const std::string& seq);

    /**
     * Estimate probability that a single read is ancient DNA
     * Based on T enrichment at 5' end and A enrichment at 3' end
     * (characteristic of C→T and G→A deamination damage)
     *
     * @param seq DNA sequence
     * @return Probability 0.0-1.0 that this read shows ancient damage
     */
    static float estimate_ancient_prob(const std::string& seq);

    /**
     * Compute observed damage magnitude for a read
     * Counts T's at 5' end positions 0-2 and A's at 3' end positions 0-2
     * Position-weighted: pos 0 = 1.0, pos 1 = 0.5, pos 2 = 0.25
     *
     * @param seq DNA sequence
     * @return Damage score 0.0-3.5 (sum of 5' and 3' scores)
     */
    static float compute_damage_score(const std::string& seq);

    /**
     * Compute per-read damage percentage using sample-level distribution
     * Similar to mapDamage/DamageProfiler approach:
     * - Uses sample's damage rates at each position
     * - Compares read's terminal bases to expected damage pattern
     * - Returns percentage (0-100) representing damage contribution
     *
     * @param seq DNA sequence
     * @param sample_profile Pre-computed sample-level damage statistics
     * @return Damage percentage 0-100
     */
    static float compute_damage_percentage(
        const std::string& seq,
        const SampleDamageProfile& sample_profile);

    /**
     * Enhanced damage probability using codon context
     * Exploits the fact that C→T damage in coding regions creates
     * specific codon patterns (e.g., stop codons from CAA→TAA)
     *
     * @param seq DNA sequence
     * @param frame Reading frame (0, 1, or 2)
     * @param forward True for forward strand
     * @return Probability 0.0-1.0 that this read shows ancient damage
     */
    static float estimate_ancient_prob_with_codons(
        const std::string& seq,
        int frame,
        bool forward);

    /**
     * Estimate damage probability using sample-level profile
     * Much more accurate than per-read estimation alone
     *
     * @param seq DNA sequence
     * @param sample_profile Pre-computed sample-level damage statistics
     * @return Probability 0.0-1.0 that this read shows ancient damage
     */
    static float estimate_ancient_prob_with_sample_profile(
        const std::string& seq,
        const SampleDamageProfile& sample_profile);

    /**
     * Advanced damage estimation combining all signals
     * Uses quality scores, dinucleotide context, and joint probability model
     *
     * @param seq DNA sequence
     * @param quality Quality string (Phred+33 encoded)
     * @param sample_profile Pre-computed sample-level damage statistics
     * @return Probability 0.0-1.0 that this read shows ancient damage
     */
    static float estimate_ancient_prob_advanced(
        const std::string& seq,
        const std::string& quality,
        const SampleDamageProfile& sample_profile);

    /**
     * Codon-aware damage detection using reading frame information
     *
     * Exploits key biological insight: C→T damage in coding regions
     * shows codon-position bias:
     * - Position 3 (wobble): often synonymous, more tolerated
     * - Position 1,2: usually changes amino acid
     *
     * Also detects "damaged stop codons":
     * - TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
     *
     * @param seq DNA sequence (in correct strand orientation)
     * @param frame Reading frame offset (0, 1, or 2)
     * @param sample_profile Pre-computed sample-level damage statistics
     * @return Probability 0.0-1.0 that this read shows ancient damage
     */
    static float estimate_ancient_prob_codon_aware(
        const std::string& seq,
        int frame,
        const SampleDamageProfile& sample_profile);

    /**
     * Score how well damage pattern matches a specific reading frame
     * Higher score = damage is more consistent with this frame being correct
     *
     * @param seq DNA sequence
     * @param frame Reading frame offset (0, 1, or 2)
     * @return Score indicating damage-frame consistency (0.0-1.0)
     */
    static float score_damage_frame_consistency(
        const std::string& seq,
        int frame);

    /**
     * Compute sample-level damage profile from a collection of reads
     * This should be called once on the first pass through the data
     *
     * @param sequences Vector of DNA sequences
     * @return Computed SampleDamageProfile with aggregate statistics
     */
    static SampleDamageProfile compute_sample_profile(
        const std::vector<std::string>& sequences);

    /**
     * Update sample profile incrementally with a new read
     * Thread-safe for parallel processing
     *
     * @param profile Profile to update (should be protected by mutex)
     * @param seq New sequence to add
     */
    static void update_sample_profile(
        SampleDamageProfile& profile,
        const std::string& seq);

    /**
     * Finalize sample profile after all reads have been added
     * Computes damage rates from raw frequencies
     *
     * @param profile Profile to finalize
     */
    static void finalize_sample_profile(SampleDamageProfile& profile);

    /**
     * Merge two sample profiles (for parallel aggregation)
     * Adds counts from src into dst
     *
     * @param dst Destination profile (modified in place)
     * @param src Source profile to merge from
     */
    static void merge_sample_profiles(SampleDamageProfile& dst, const SampleDamageProfile& src);

    /**
     * Update sample profile with a weighted contribution
     * Used for iterative damage refinement - weights by coding probability
     *
     * @param profile Profile to update
     * @param seq DNA sequence to add
     * @param weight Weight for this sequence (0.0-1.0, typically coding_prob)
     */
    static void update_sample_profile_weighted(
        SampleDamageProfile& profile,
        const std::string& seq,
        float weight);

    /**
     * Reset sample profile for a fresh damage estimation pass
     * Clears all counts while preserving structure
     *
     * @param profile Profile to reset
     */
    static void reset_sample_profile(SampleDamageProfile& profile);

    /**
     * Bayesian damage detection using proper likelihood ratios
     * P(seq|ancient) / P(seq|modern) with all signals
     *
     * @param seq DNA sequence
     * @param quality Quality string (empty if not available)
     * @param frame Reading frame (0, 1, or 2)
     * @param sample_profile Pre-computed sample statistics
     * @return Log-likelihood ratio (positive = ancient, negative = modern)
     */
    static float compute_damage_log_likelihood(
        const std::string& seq,
        const std::string& quality,
        int frame,
        const SampleDamageProfile& sample_profile);

    /**
     * Quality-aware damage probability
     * Low quality at damage positions strengthens the signal
     *
     * @param seq DNA sequence
     * @param quality Quality string (Phred+33)
     * @param frame Reading frame
     * @param sample_profile Pre-computed sample statistics
     * @return Probability 0.0-1.0 that this read shows ancient damage
     */
    static float estimate_ancient_prob_with_quality(
        const std::string& seq,
        const std::string& quality,
        int frame,
        const SampleDamageProfile& sample_profile);

    // Scoring helper functions (public for use by optimized frame selection)
    static float calculate_codon_usage_score(const std::vector<std::string>& codons);
    static float calculate_stop_penalty(const std::string& protein);
    static float calculate_aa_composition_score(const std::string& protein);
    static float calculate_damage_consistency(
        const std::string& original_seq,
        const std::string& frame_seq,
        int frame,
        bool forward,
        const DamageProfile& damage);

    /**
     * Hexamer-aware frame selection using GTDB frequencies
     *
     * Uses hexamer log-likelihood ratios for frame scoring with damage correction:
     * 1. Calculates LLR for each hexamer against uniform background
     * 2. Identifies damage-indicating hexamers (TGA*, TAG*, TAA* starts)
     * 3. Adjusts LLR using corrected hexamer frequencies for damage hexamers
     * 4. Returns best frame based on total hexamer evidence
     *
     * @param seq DNA sequence
     * @param terminal_length Terminal region for damage detection (default 15bp)
     * @return Pair of FrameScores: (best_forward, best_reverse)
     */
    static std::pair<FrameScore, FrameScore> select_best_per_strand_hexamer_aware(
        const std::string& seq,
        size_t terminal_length = 15);

    /**
     * Stop-only frame selection - uses ONLY stop codon count
     *
     * Analysis of frame selection accuracy on CDS fragments showed that:
     * - Stop codon count alone achieves 97.42% accuracy
     * - Combined scoring with multiple features achieves only ~46% accuracy
     * - Other features (codon usage, hexamer, etc.) add noise for CDS fragments
     *   because all 6 frames look equally "coding-like"
     *
     * This method selects frames based purely on stop codon count:
     * - Fewer stops = higher score
     * - Ties broken by frame 0 preference (standard codon position)
     *
     * @param seq DNA sequence
     * @return Pair of FrameScores: (best_forward, best_reverse)
     */
    static std::pair<FrameScore, FrameScore> select_best_per_strand_stop_only(
        const std::string& seq);

    /**
     * Stop-priority frame selection: stops as primary signal, hexamer LLR as tiebreaker
     *
     * Algorithm:
     * 1. Group frames by internal stop count
     * 2. Within each stop-count group, rank by hexamer LLR score
     * 3. Return best frame per strand
     *
     * This combines the high accuracy of stop-based selection with
     * hexamer evidence for breaking ties when multiple frames have same stop count.
     *
     * @param seq DNA sequence
     * @return Pair of FrameScores: (best_forward, best_reverse)
     */
    static std::pair<FrameScore, FrameScore> select_best_per_strand_stop_priority(
        const std::string& seq);

    /**
     * Return all 6 reading frames with full scoring information
     *
     * Used for --all-frames mode where all 6 translations are output per read.
     * Returns frames sorted by score (best first).
     *
     * @param seq DNA sequence
     * @param damage Damage profile (optional)
     * @return Vector of 6 FrameScores, sorted by total_score descending
     */
    static std::vector<FrameScore> score_all_frames_full(
        const std::string& seq,
        const DamageProfile* damage = nullptr);

    /**
     * Calculate wobble position T enrichment for strand prediction
     * In damaged DNA, C→T at wobble position (codon pos 3) is often synonymous
     * If T's are enriched at wobble position, this suggests correct strand/frame
     *
     * @param seq DNA sequence (in correct strand orientation)
     * @param frame Reading frame offset (0, 1, or 2)
     * @param terminal_region Number of bp from 5' end to analyze (default: 15)
     * @return Wobble enrichment score (0-1, higher = more T's at wobble pos)
     */
    static float calculate_wobble_enrichment(
        const std::string& seq,
        int frame,
        size_t terminal_region = 15);

    /**
     * Detect damage-induced stop codons at 5' end
     * Looks for TAA/TAG/TGA that could be from C→T damage:
     * - TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
     *
     * @param seq DNA sequence (in correct strand orientation)
     * @param frame Reading frame offset (0, 1, or 2)
     * @return Score 0-1 (0 = has damage-induced stop, 1 = no suspicious stops)
     */
    static float detect_damage_induced_stops(
        const std::string& seq,
        int frame);

    /**
     * Calculate strand confidence combining all signals
     * @param seq Original read sequence
     * @param fwd_score Best forward strand frame score
     * @param rev_score Best reverse strand frame score
     * @param fwd_frame Best forward frame
     * @param rev_frame Best reverse frame
     * @param sample_profile Sample damage statistics
     * @return Strand confidence (0.5 = uncertain, approaching 1.0 = confident)
     */
    static float calculate_strand_confidence(
        const std::string& seq,
        float fwd_score,
        float rev_score,
        int fwd_frame,
        int rev_frame,
        const SampleDamageProfile& sample_profile);

    /**
     * Compute damage codon score for iterative enrichment
     *
     * Scores reads by presence of damage-diagnostic codon patterns:
     * 1. Terminal T's in first 3 positions (C→T damage signature)
     * 2. Internal stop codons near 5' (TAA/TAG/TGA from CAA/CAG/CGA damage)
     * 3. T-rich codons at 5' end that are consistent with damaged C-containing codons
     *
     * This allows filtering for reads that show actual damage signal,
     * rather than just coding potential.
     *
     * @param seq DNA sequence
     * @return Damage codon score 0.0-1.0 (higher = more damage-indicative)
     */
    static float compute_damage_codon_score(const std::string& seq);
};

// Free functions for hexamer-based scoring and correction

/**
 * Correct stop codons using hexamer evidence
 * Fixes stop codons that are likely damage artifacts (TAA/TAG/TGA from C->T)
 * Uses GTDB hexamer frequencies to assess damage probability
 *
 * @param seq DNA sequence
 * @param frame Reading frame (0, 1, 2)
 * @param terminal_length Terminal region to consider (default 15bp)
 * @return Pair of (corrected sequence, number of corrections)
 */
std::pair<std::string, size_t> correct_stop_codons_hexamer(
    const std::string& seq,
    int frame,
    size_t terminal_length = 15);

/**
 * Calculate hexamer log-likelihood ratio score
 * Uses GTDB hexamer frequencies vs uniform background
 * Higher score = more likely to be real coding sequence
 *
 * @param seq DNA sequence
 * @param frame Reading frame (0, 1, 2)
 * @return Average log-likelihood ratio per hexamer
 */
float calculate_hexamer_llr_score(const std::string& seq, int frame);

/**
 * Damage-aware stop penalty using hexamer evidence
 * Distinguishes between likely-real stops and damage-induced stops
 *
 * @param seq DNA sequence
 * @param protein Translated protein
 * @param frame Reading frame (0, 1, 2)
 * @param terminal_length Terminal region to consider (default 15bp)
 * @return Multiplier 0-1 (1 = no real stops, 0.05 = multiple real stops)
 */
float calculate_stop_penalty_hexamer_aware(
    const std::string& seq,
    const std::string& protein,
    int frame,
    size_t terminal_length = 15);

} // namespace agp
