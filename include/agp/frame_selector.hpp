#pragma once

#include "types.hpp"
#include "hexamer_tables.hpp"  // For Domain enum
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
    // Negative control counts at 5' end (for A/(A+G) ratio)
    std::array<double, 15> a_freq_5prime = {};  // A count at each position (control)
    std::array<double, 15> g_freq_5prime = {};  // G count at each position (control)

    // Position-specific base counts at 3' end (positions 0-14 from end)
    std::array<double, 15> a_freq_3prime = {};  // A count at each position
    std::array<double, 15> g_freq_3prime = {};  // G count (baseline)
    // Negative control counts at 3' end (for T/(T+C) ratio)
    std::array<double, 15> t_freq_3prime = {};  // T count at each position (control)
    std::array<double, 15> c_freq_3prime = {};  // C count at each position (control)

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
    // Raw totals for significance testing
    std::array<double, 15> tc_total_5prime = {};  // T+C counts at 5'
    std::array<double, 15> ag_total_3prime = {};  // A+G counts at 3'

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

    // Briggs-like damage model parameters (fast closed-form estimation)
    // Model: δ(pos) = δ_s * P_overhang(pos) + δ_d * (1 - P_overhang(pos))
    // Where: P_overhang(pos) ≈ (1-λ)^pos (geometric distribution approximation)
    float delta_s_5prime = 0.0f;  // Single-stranded deamination rate at 5' end
    float delta_d_5prime = 0.0f;  // Double-stranded (background) deamination at 5'
    float delta_s_3prime = 0.0f;  // Single-stranded deamination rate at 3' end
    float delta_d_3prime = 0.0f;  // Double-stranded (background) deamination at 3'
    float r_squared_5prime = 0.0f;  // Goodness of fit for 5' model
    float r_squared_3prime = 0.0f;  // Goodness of fit for 3' model

    // Codon-position-specific damage analysis (exploits reading frame information)
    // Ancient DNA shows characteristic pattern: pos3 (wobble) > pos1 > pos2
    // because wobble position damage is often silent (synonymous)
    float codon_pos1_damage = 0.0f;  // C→T rate at codon position 1
    float codon_pos2_damage = 0.0f;  // C→T rate at codon position 2
    float codon_pos3_damage = 0.0f;  // C→T rate at codon position 3 (wobble)
    float wobble_ratio = 1.0f;       // pos3 / ((pos1 + pos2) / 2), >1 indicates ancient
    float hexamer_damage_llr = 0.0f; // Hexamer-based damage log-likelihood ratio
    float terminal_shift_5prime = 0.0f;  // terminal T/(T+C) - interior baseline
    float terminal_shift_3prime = 0.0f;  // terminal A/(A+G) - interior baseline
    float terminal_z_5prime = 0.0f;      // z-score for 5' terminal enrichment
    float terminal_z_3prime = 0.0f;      // z-score for 3' terminal enrichment
    bool terminal_inversion = false;     // true if terminals show significant depletion

    // Negative control statistics (should NOT show enrichment if damage is real)
    // 5' control: A/(A+G) at 5' end - real C→T damage shouldn't affect this
    // 3' control: T/(T+C) at 3' end - real G→A damage shouldn't affect this
    float ctrl_shift_5prime = 0.0f;   // 5' A/(A+G) terminal - interior
    float ctrl_shift_3prime = 0.0f;   // 3' T/(T+C) terminal - interior
    float ctrl_z_5prime = 0.0f;       // z-score for 5' control enrichment
    float ctrl_z_3prime = 0.0f;       // z-score for 3' control enrichment
    bool composition_bias_5prime = false;  // control comparable to damage → bias
    bool composition_bias_3prime = false;  // control comparable to damage → bias

    // Library type detection
    enum class LibraryType { UNKNOWN, DOUBLE_STRANDED, SINGLE_STRANDED };
    LibraryType library_type = LibraryType::DOUBLE_STRANDED;  // Default to double-stranded
    LibraryType forced_library_type = LibraryType::UNKNOWN;  // User override (UNKNOWN = auto-detect)

    // Inverted pattern detection (terminal T/(T+C) < interior)
    // When true, reference-free damage detection failed due to:
    // - AT-rich organisms with terminal artifacts
    // - Adapter contamination
    // - Quality trimming bias
    // Users should use metaDMG or other reference-based tools
    bool inverted_pattern_5prime = false;  // 5' terminal T/(T+C) < interior
    bool inverted_pattern_3prime = false;  // 3' terminal A/(A+G) < interior
    float terminal_gradient_5prime = 0.0f;  // pos0 - pos10-14 average (negative = inverted)
    float terminal_gradient_3prime = 0.0f;  // pos0 - pos10-14 average (negative = inverted)

    // Exponential fit parameters: p(pos) = b + A * exp(-lambda * pos)
    // Fitted baseline (asymptotic T/(T+C) or A/(A+G))
    float fit_baseline_5prime = 0.0f;   // b parameter for 5' end
    float fit_baseline_3prime = 0.0f;   // b parameter for 3' end
    // Fitted amplitude (damage signal above baseline)
    float fit_amplitude_5prime = 0.0f;  // A parameter for 5' end
    float fit_amplitude_3prime = 0.0f;  // A parameter for 3' end
    // Fit quality (RMSE of residuals)
    float fit_rmse_5prime = 0.0f;
    float fit_rmse_3prime = 0.0f;

    // Calibrated D_max values (comparable to metaDMG)
    // D = A / (1 - b), the fraction of C that became T
    float d_max_5prime = 0.0f;  // Calibrated D_max for 5' end
    float d_max_3prime = 0.0f;  // Calibrated D_max for 3' end
    float d_max_combined = 0.0f;  // Final D_max using asymmetry-aware combination
    float asymmetry = 0.0f;  // |D_5p - D_3p| / ((D_5p + D_3p) / 2)
    bool high_asymmetry = false;  // True if asymmetry > 0.5 (possible artifact)

    // Track source of d_max_combined estimate
    enum class DmaxSource { AVERAGE, MIN_ASYMMETRY, FIVE_PRIME_ONLY, THREE_PRIME_ONLY, NONE };
    DmaxSource d_max_source = DmaxSource::AVERAGE;

    const char* d_max_source_str() const {
        switch (d_max_source) {
            case DmaxSource::AVERAGE: return "average";
            case DmaxSource::MIN_ASYMMETRY: return "min_asymmetry";
            case DmaxSource::FIVE_PRIME_ONLY: return "5prime_only";
            case DmaxSource::THREE_PRIME_ONLY: return "3prime_only";
            case DmaxSource::NONE: return "none";
            default: return "unknown";
        }
    }

    size_t n_reads = 0;  // Number of reads used in computation
    size_t n_reads_gc_filtered = 0;  // Reads skipped due to low GC content
    size_t n_reads_sampled = 0;  // Total reads sampled for GC histogram

    // GC content histogram (100 bins: 0-1%, 1-2%, ..., 99-100%)
    // Used to compute adaptive GC threshold for damage detection
    std::array<size_t, 100> gc_histogram = {};
    float adaptive_gc_threshold = 0.0f;  // Computed from histogram (70th percentile)
    bool gc_threshold_computed = false;

    // Hexamer-based damage detection
    // Track hexamer counts at 5' terminal positions (first 6 bases)
    // For each hexamer, we count occurrences at terminal vs interior positions
    // C→T damage should show excess T-hexamers at terminals relative to expected
    std::array<double, 4096> hexamer_count_5prime = {};  // Hexamer counts at 5' (pos 0-5)
    std::array<double, 4096> hexamer_count_interior = {}; // Hexamer counts at interior
    size_t n_hexamers_5prime = 0;    // Total hexamers counted at 5' terminal
    size_t n_hexamers_interior = 0;  // Total hexamers counted at interior

    // Hexamer-based T/(T+C) ratios (more reliable than position 0 or 1 alone)
    // These average positions 1-6 and are less affected by first-base artifacts
    float hexamer_terminal_tc = 0.0f;   // T/(T+C) at terminal from hexamer analysis
    float hexamer_interior_tc = 0.0f;   // T/(T+C) at interior from hexamer analysis
    float hexamer_excess_tc = 0.0f;     // Terminal - interior (negative = inverted)

    // Likelihood-based model comparison (exponential decay vs constant)
    // Positive LLR = exponential fits better (real decay pattern)
    // Negative/zero LLR = constant fits better (no decay, likely composition bias)
    float decay_llr_5prime = 0.0f;      // Log-likelihood ratio for 5' damage channel
    float decay_llr_3prime = 0.0f;      // Log-likelihood ratio for 3' damage channel
    float ctrl_decay_llr_5prime = 0.0f; // Log-likelihood ratio for 5' control channel
    float ctrl_decay_llr_3prime = 0.0f; // Log-likelihood ratio for 3' control channel

    // Delta LLR: damage channel - control channel
    // Positive = damage channel has stronger decay than control = real damage
    // Near zero = both channels decay similarly = composition/trimming artifact
    float delta_llr_5prime = 0.0f;      // decay_llr - ctrl_decay_llr at 5'
    float delta_llr_3prime = 0.0f;      // decay_llr - ctrl_decay_llr at 3'

    // Channel divergence (damage channel vs control channel)
    // High divergence = real damage (only damage channel elevated)
    // Low divergence = composition bias (both channels elevated together)
    float channel_divergence_5prime = 0.0f;  // |damage_shift - control_shift| at 5'
    float channel_divergence_3prime = 0.0f;  // |damage_shift - control_shift| at 3'

    // =========================================================================
    // CHANNEL B: Convertible stop codon tracking (translation disruption signal)
    // Independent detector for real C→T damage vs compositional variation
    //
    // Tracks CAA→TAA, CAG→TAG, CGA→TGA conversions by nucleotide position
    // Real damage: should show exponential decay matching Channel A
    // Compositional variation: should be flat (no position dependence)
    // =========================================================================

    // Convertible codon pair counts at 5' end by nucleotide position (0-14)
    // Position = nucleotide position of the C/T in the codon (from read start)
    // For CAA/TAA: position of the first base (C or T)
    // Exposure = CAA + TAA, Stops = TAA
    std::array<double, 15> convertible_caa_5prime = {};  // CAA (Gln) codons
    std::array<double, 15> convertible_taa_5prime = {};  // TAA (Stop) codons
    std::array<double, 15> convertible_cag_5prime = {};  // CAG (Gln) codons
    std::array<double, 15> convertible_tag_5prime = {};  // TAG (Stop) codons
    std::array<double, 15> convertible_cga_5prime = {};  // CGA (Arg) codons
    std::array<double, 15> convertible_tga_5prime = {};  // TGA (Stop) codons

    // Total codons observed at each position (denominator for exposure)
    std::array<double, 15> total_codons_5prime = {};

    // Interior reference counts (for baseline estimation, positions 30+)
    double convertible_caa_interior = 0.0;
    double convertible_taa_interior = 0.0;
    double convertible_cag_interior = 0.0;
    double convertible_tag_interior = 0.0;
    double convertible_cga_interior = 0.0;
    double convertible_tga_interior = 0.0;
    double total_codons_interior = 0.0;

    // Computed statistics for Channel B
    float stop_conversion_rate_baseline = 0.0f;  // Interior stop/(pre+stop) ratio
    float stop_decay_llr_5prime = 0.0f;  // LLR for stop position decay (Channel B)
    float stop_amplitude_5prime = 0.0f;  // Fitted amplitude of stop excess
    bool channel_b_valid = false;  // True if sufficient data for Channel B

    // Joint evidence decision
    bool damage_validated = false;  // True if both channels agree on damage
    bool damage_artifact = false;   // True if Channel A fires but Channel B doesn't

    // Compute adaptive GC threshold from histogram
    float compute_adaptive_gc_threshold(float target_percentile = 0.70f) {
        size_t total = 0;
        for (auto c : gc_histogram) total += c;
        if (total < 1000) return 0.40f;  // Default if insufficient data

        size_t target_count = static_cast<size_t>(total * target_percentile);
        size_t cumulative = 0;
        for (int i = 0; i < 100; ++i) {
            cumulative += gc_histogram[i];
            if (cumulative >= target_count) {
                return static_cast<float>(i) / 100.0f;
            }
        }
        return 0.50f;  // Fallback
    }

    bool is_valid() const { return n_reads >= 1000; }  // Need enough reads

    // Check if damage detection is unreliable
    // Unreliable if: hexamer inverted at either end, OR composition bias at either end
    bool is_detection_unreliable() const {
        return inverted_pattern_5prime || inverted_pattern_3prime ||
               composition_bias_5prime || composition_bias_3prime;
    }

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
     * Select the best reading frame from each strand with sample-level damage info
     * Uses damage pattern asymmetry for improved strand discrimination
     *
     * @param seq DNA sequence
     * @param damage Optional per-read damage profile
     * @param sample_profile Sample-level damage statistics for strand discrimination
     * @return Pair of FrameScores: (best_forward, best_reverse)
     */
    static std::pair<FrameScore, FrameScore> select_best_per_strand(
        const std::string& seq,
        const DamageProfile* damage,
        const SampleDamageProfile* sample_profile);

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
     * Get reverse complement of a sequence (returns new string)
     */
    static std::string reverse_complement(const std::string& seq);

    /**
     * Get reverse complement into thread-local buffer (avoids allocation)
     * Returns reference to internal buffer - valid until next call from same thread
     * Use this when you only need temporary access to the RC sequence
     */
    static const std::string& reverse_complement_cached(const std::string& seq);

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

    /**
     * Bayesian frame selection using hexamer log-likelihood ratios
     *
     * Computes posterior probability over all 6 frames + noncoding hypothesis,
     * and emits 0, 1, or 2 frames based on confidence thresholds:
     * - 1 frame: posterior > tau_confident (default 0.5)
     * - 2 frames: top2 sum > tau_ambiguous AND top1 > tau_noncoding
     * - 0 frames: max coding posterior < tau_noncoding
     *
     * Returns pair where:
     * - First frame is always the best (highest posterior)
     * - Second frame is second-best if ambiguous, or invalid if confident/noncoding
     * - Check total_score >= 0 to determine if frame is valid
     *
     * @param seq DNA sequence
     * @param domain Taxonomic domain for hexamer table lookup (default: GTDB)
     * @param sample_profile Sample-level damage statistics (optional, for damage-aware scoring)
     * @return Pair of FrameScores: (best_frame, second_frame_or_invalid)
     */
    static std::pair<FrameScore, FrameScore> select_best_per_strand_bayesian(
        const std::string& seq,
        Domain domain = Domain::GTDB,
        const SampleDamageProfile* sample_profile = nullptr);
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
 * Calculate strand preference based on ancient DNA damage asymmetry
 * Forward strand: expect T-enrichment at 5' (C→T), A-enrichment at 3' (G→A complement)
 *
 * @param seq Forward strand sequence
 * @param rc_seq Reverse complement sequence
 * @param profile Sample-level damage profile
 * @return Positive for forward preference, negative for reverse
 */
float calculate_damage_strand_preference(
    const std::string& seq,
    const std::string& rc_seq,
    const SampleDamageProfile& profile);

} // namespace agp
