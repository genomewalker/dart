#pragma once

#include "types.hpp"
#include "hexamer_tables.hpp"  // For Domain enum
#include "joint_damage_model.hpp"  // Joint probabilistic damage model
#include "mixture_damage_model.hpp"  // Durbin-style mixture model
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <utility>

namespace agp {

// Forward declarations
class DamageModel;
struct UnifiedDamageContext;

/**
 * Mixture model results for organism-class heterogeneity
 * Summarizes the Durbin-style K-component model output
 */
struct MixtureResult {
    bool converged = false;        // Did EM converge?
    int K = 0;                     // Number of classes selected by BIC
    float d_population = 0.0f;     // E[δ] over all C-sites
    float d_ancient = 0.0f;        // E[δ | δ > 5%] (ancient tail)
    float d_reference = 0.0f;      // E[δ | GC > 50%] (metaDMG proxy)
    float pi_ancient = 0.0f;       // Fraction of C-sites in high-damage classes
    float bic = 0.0f;              // BIC for model selection
};

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
    bool position_0_artifact_5prime = false;  // pos0 depleted but pos1 enriched (adapter bias)
    bool position_0_artifact_3prime = false;  // pos0 depleted but pos1 enriched (adapter bias)

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
    enum class DmaxSource { AVERAGE, MIN_ASYMMETRY, FIVE_PRIME_ONLY, THREE_PRIME_ONLY, CHANNEL_B_STRUCTURAL, NONE };
    DmaxSource d_max_source = DmaxSource::AVERAGE;

    const char* d_max_source_str() const {
        switch (d_max_source) {
            case DmaxSource::AVERAGE: return "average";
            case DmaxSource::MIN_ASYMMETRY: return "min_asymmetry";
            case DmaxSource::FIVE_PRIME_ONLY: return "5prime_only";
            case DmaxSource::THREE_PRIME_ONLY: return "3prime_only";
            case DmaxSource::CHANNEL_B_STRUCTURAL: return "channel_b_structural";
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

    // Channel B: Convertible stop codon tracking (CAA→TAA, CAG→TAG, CGA→TGA)
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

    // Channel B structural d_max from multi-position stop codon conversion
    // WLS model: r_p = b0 + (1-b0) * d_max * exp(-λp)
    float d_max_from_channel_b = 0.0f;   // Structural d_max estimate from stop codons
    float channel_b_weight = 0.0f;       // Exposure weight W_B for joint likelihood
    float channel_b_slope = 0.0f;        // Raw WLS slope (positive = damage, negative = inverted)
    bool channel_b_quantifiable = false; // True if Channel B can provide d_max estimate
    bool channel_b_inverted = false;     // True if slope <= 0 (terminal stops LOWER than baseline)

    // GC-stratified damage estimation (per-bin to handle metagenome heterogeneity)
    static constexpr int N_GC_BINS = 10;  // 0-10%, 10-20%, ..., 90-100%

    struct GCBinStats {
        // Channel A: T and C counts at terminal positions (0-14)
        std::array<uint64_t, 15> t_counts = {};
        std::array<uint64_t, 15> c_counts = {};
        // Channel A: interior baseline counts
        uint64_t t_interior = 0;
        uint64_t c_interior = 0;

        // Channel B: stop codon counts at terminal positions
        std::array<uint64_t, 15> stop_counts = {};  // TAA+TAG+TGA
        std::array<uint64_t, 15> pre_counts = {};   // CAA+CAG+CGA
        // Channel B: interior baseline
        uint64_t stop_interior = 0;
        uint64_t pre_interior = 0;

        // Computed results
        float d_max = 0.0f;           // Estimated damage for this bin
        float d_max_channel_b = 0.0f; // Channel B estimate for this bin
        uint64_t n_reads = 0;         // Number of reads in this bin
        uint64_t c_sites = 0;         // Total C sites (for weighting)
        bool valid = false;           // Sufficient data for estimation

        // GC-conditional damage parameters (for per-read inference)
        float p_damaged = 0.0f;       // P(damaged) for this bin from LLR classification
        float baseline_tc = 0.5f;     // Interior T/(T+C) baseline for this bin
        float llr = 0.0f;             // Log-likelihood ratio (positive = damaged)
        bool classified_damaged = false;  // Hard classification from LLR threshold

        // Terminal observation counts (for per-read LLR update)
        uint64_t n_terminal_obs() const {
            uint64_t n = 0;
            for (int p = 0; p < 15; ++p) n += t_counts[p] + c_counts[p];
            return n;
        }
    };

    std::array<GCBinStats, N_GC_BINS> gc_bins = {};

    // Aggregated GC-stratified results
    float gc_stratified_d_max_weighted = 0.0f;  // Weighted average across bins
    float gc_stratified_d_max_peak = 0.0f;      // Max d_max across valid bins
    int gc_peak_bin = -1;                        // Which bin has peak damage
    bool gc_stratified_valid = false;            // At least one bin has valid estimate

    // GC-conditional damage summary
    float pi_damaged = 0.0f;          // Fraction of terminal obs from damaged bins
    float d_ancient = 0.0f;           // E[δ | damaged bins] - severity among damaged
    float d_population = 0.0f;        // E[δ] over all bins - average across all DNA
    int n_damaged_bins = 0;           // Number of bins classified as damaged

    // Joint evidence decision (legacy two-channel)
    bool damage_validated = false;  // True if both channels agree on damage
    bool damage_artifact = false;   // True if Channel A fires but Channel B doesn't

    // Joint probabilistic model results (Bayesian Channel A + Control + Channel B)
    float joint_delta_max = 0.0f;      // MLE estimate of damage rate
    float joint_lambda = 0.0f;         // Decay constant
    float joint_a_max = 0.0f;          // Artifact amplitude (signed)
    float joint_log_lik_m1 = 0.0f;     // Log-likelihood for M1 (damage)
    float joint_log_lik_m0 = 0.0f;     // Log-likelihood for M0 (no damage)
    float joint_delta_bic = 0.0f;      // BIC_M0 - BIC_M1 (positive = damage)
    float joint_bayes_factor = 0.0f;   // BF_10 ≈ exp(ΔBIC/2)
    float joint_p_damage = 0.0f;       // P(damage | data)
    bool joint_model_valid = false;    // Sufficient data for joint model

    // Durbin-style mixture model results (K-component over organism classes)
    MixtureResult mixture;

    // d_metamatch: Channel B-anchored estimate for metaDMG comparability
    // Formula: d_metamatch = d_global + γ × (d_channel_b - d_global)
    struct MetamatchEstimate {
        float d_metamatch = 0.0f;              // Calibrated metaDMG-comparable estimate
        float d_alignability_weighted = 0.0f;  // Raw alignability-weighted d_max
        float gamma = 0.0f;                    // Blending coefficient (0 = use d_global, 1 = use weighted)
        float mean_alignability = 0.0f;        // Mean alignability score across reads
        float alignability_damage_corr = 0.0f; // Correlation between alignability and per-read damage
    };
    MetamatchEstimate metamatch;

    // Compute adaptive GC threshold from histogram
    float compute_adaptive_gc_threshold(float target_percentile = 0.70f);

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

    // GC-conditional per-read damage helpers
    // Get GC bin index (0-9) for a sequence based on interior GC content
    static int get_gc_bin(const std::string& seq);

    // Get GC-conditional damage parameters for a read
    struct GCDamageParams {
        float p_damaged;     // Prior P(damaged) from bin classification
        float delta_s;       // Damage rate δ_s for this bin
        float baseline_tc;   // Interior T/(T+C) baseline
        float lambda;        // Decay constant (shared across bins)
        bool bin_valid;      // Whether this bin has valid estimates
    };

    GCDamageParams get_gc_params(const std::string& seq) const;

    // Get effective damage rate for a read (posterior-weighted)
    // δ_eff = P(damaged|read) * δ_s
    float get_effective_damage(const std::string& seq, float read_ancient_prob) const;
};

// Tri-state damage validation (Channel A + Channel B decision)
enum class DamageValidationState {
    VALIDATED,      // Both channels fired → full damage
    CONTRADICTED,   // Channel A fired but Channel B negative → artifact
    UNVALIDATED     // Insufficient data for Channel B → soft suppression
};

inline DamageValidationState get_damage_validation_state(const SampleDamageProfile& profile) {
    if (profile.damage_validated) {
        return DamageValidationState::VALIDATED;
    }
    if (profile.damage_artifact) {
        return DamageValidationState::CONTRADICTED;
    }
    // Channel B has sufficient data but contradicts Channel A
    if (profile.channel_b_valid && profile.stop_decay_llr_5prime < -100.0f) {
        return DamageValidationState::CONTRADICTED;
    }
    return DamageValidationState::UNVALIDATED;
}

// Returns suppression factor: 1.0 (full), 0.0 (hard suppress), 0.5 (soft suppress)
inline float get_damage_suppression_factor(DamageValidationState state) {
    switch (state) {
        case DamageValidationState::VALIDATED:    return 1.0f;
        case DamageValidationState::CONTRADICTED: return 0.0f;
        case DamageValidationState::UNVALIDATED:  return 0.5f;
    }
    return 1.0f;  // Unreachable, but satisfies compiler
}

// Convenience: get factor directly from profile
inline float get_damage_suppression_factor(const SampleDamageProfile& profile) {
    return get_damage_suppression_factor(get_damage_validation_state(profile));
}

/**
 * Score for a single reading frame
 * Combines multiple signals to evaluate frame correctness
 */
struct FrameScore {
    int frame = 0;              // Frame offset (0, 1, 2)
    bool forward = true;        // Strand direction

    float log_likelihood = 0.0f;        // Log-likelihood from Bayesian frame selection
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
     * Estimate probability that a single read is ancient DNA
     * Based on T enrichment at 5' end and A enrichment at 3' end
     * (characteristic of C→T and G→A deamination damage)
     *
     * @param seq DNA sequence
     * @return Probability 0.0-1.0 that this read shows ancient damage
     */
    static float estimate_damage_signal(const std::string& seq);

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
     * Unified damage detection using GC-conditional Bayesian inference
     *
     * Exploits key biological insight: C→T damage in coding regions
     * shows codon-position bias:
     * - Position 3 (wobble): often synonymous, more tolerated
     * - Position 1,2: usually changes amino acid
     *
     * Also detects "damaged stop codons":
     * - TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
     *
     * Uses GC-conditional parameters for proper Bayesian posterior:
     * - Prior from GC-bin specific p_damaged
     * - Damage rates from GC-bin specific delta_s
     * - Quality weighting when quality scores available
     *
     * @param seq DNA sequence (in correct strand orientation)
     * @param frame Reading frame offset (0, 1, or 2)
     * @param sample_profile Pre-computed sample-level damage statistics
     * @param quality Optional quality string (Phred+33), empty for FASTA
     * @return Probability 0.0-1.0 that this read shows ancient damage
     */
    static float estimate_damage_signal(
        const std::string& seq,
        int frame,
        const SampleDamageProfile& sample_profile,
        const std::string& quality = "");


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
     * ORF fragment from damage-aware enumeration
     */
    struct ORFFragment {
        std::string protein;          // Observed protein (stops as *, for biological output)
        std::string observed_protein; // Observed protein with X markers (for debugging)
        std::string search_protein;   // Search-optimized: terminal damage stops masked as X
        std::string corrected_nt;     // Corrected nucleotide sequence (stop codons fixed)
        size_t aa_start;              // Start position in full translation (amino acids)
        size_t aa_end;                // End position (exclusive)
        size_t nt_start;              // Start position in oriented DNA (nucleotides)
        size_t nt_end;                // End position in oriented DNA (nucleotides)
        bool is_forward;              // Strand
        int frame;                    // Reading frame (0, 1, 2)
        size_t length;                // Protein length (aa)
        size_t x_count;               // Number of X (ambiguous) positions
        size_t aa_corrections;        // Non-stop AA corrections (C→T, G→A damage fixes)
        size_t passed_stops;          // Damage-induced stops continued through (legacy, always 0)
        size_t rescued_stops;         // Damage stops X-masked in search_protein (continued through)
        size_t real_stops;            // Real (non-damage) internal stops in this frame
        float damage_evidence;        // Sum of p_damage_stop for convertible stops
        float score;                  // Ranking score (coding signal - stop penalties)
    };

    /**
     * Region of consistent reading frame (for frameshift detection)
     *
     * When a read contains a frameshift (indel), the optimal frame changes
     * mid-read. This struct represents one contiguous region with consistent frame.
     */
    struct FrameshiftRegion {
        size_t codon_start;           // Start codon index in this region
        size_t codon_end;             // End codon index (exclusive)
        size_t nt_start;              // Start nucleotide position
        size_t nt_end;                // End nucleotide position
        int frame;                    // Reading frame (0, 1, 2)
        bool forward;                 // Strand direction
        std::string protein;          // Protein sequence for this region
        float score;                  // Log-likelihood score for this region
    };

    /**
     * Result of frameshift detection via Viterbi algorithm
     */
    struct FrameshiftResult {
        std::vector<FrameshiftRegion> regions;  // Regions with consistent frame
        float viterbi_score;                    // Total Viterbi path score
        float best_single_frame_score;          // Best single-frame score (no shifts)
        int best_single_frame;                  // Best single frame index
        bool has_frameshift;                    // True if frameshift detected
        size_t frameshift_position;             // Codon index of first frameshift (if any)
    };

    /**
     * Detect frameshifts using HMM/Viterbi dynamic programming
     *
     * Uses per-codon log-likelihood scores as HMM emissions and finds the
     * optimal frame path. A frameshift is reported only if:
     * 1. The Viterbi path uses multiple frames
     * 2. Each frame region has at least min_segment_codons
     * 3. The score improvement exceeds min_score_improvement
     *
     * @param seq DNA sequence
     * @param sample_profile Sample-level damage profile
     * @param frameshift_penalty Log penalty for frame transitions (default: -5.0)
     * @param min_segment_codons Minimum codons per frame region (default: 4)
     * @param min_score_improvement Required gain over single-frame (default: 2.0)
     * @return FrameshiftResult with detected regions
     */
    static FrameshiftResult detect_frameshifts(
        const std::string& seq,
        const SampleDamageProfile& sample_profile,
        float frameshift_penalty = -5.0f,
        size_t min_segment_codons = 4,
        float min_score_improvement = 2.0f);

    /**
     * Enumerate ORF fragments with damage-aware stop handling
     *
     * Translates all 6 frames and splits at "real" stop codons only.
     * Stop codons in convertible contexts (CAA→TAA, CAG→TAG, CGA→TGA)
     * near terminals with high damage probability are marked as 'X' instead
     * of being treated as stops. The damage threshold adapts to each read's
     * estimated damage level (computed via per_read_damage or auto-computed).
     *
     * OUTPUT MODES:
     * - Default (6-frame mode): Outputs best ORF per frame (up to 6 total)
     *   This guarantees frame coverage and matches traditional 6-frame translation.
     * - Adaptive mode: Outputs ORFs within score threshold of the best
     *   Use when you want the coding model to filter candidates.
     *
     * @param seq DNA sequence
     * @param sample_profile Sample-level damage profile (for d_max, lambda)
     * @param min_aa Minimum ORF length in amino acids (default: 10)
     * @param adaptive Use adaptive score-based selection instead of 6-frame mode
     * @param per_read_damage Per-read damage prior (-1 = auto-compute from terminal bases)
     * @return Vector of ORFFragments sorted by score (highest first)
     */
    static std::vector<ORFFragment> enumerate_orf_fragments(
        const std::string& seq,
        const SampleDamageProfile& sample_profile,
        size_t min_aa = 10,
        bool adaptive = false,
        float per_read_damage = -1.0f);

    /**
     * Compute per-read damage prior using Bayesian posterior on terminal bases.
     *
     * This is a fast computation that can be done BEFORE ORF enumeration.
     * Uses sample-wide d_max/lambda and per-read terminal T/A frequencies
     * with interior composition as baseline.
     *
     * Formula: P(damaged|T) = P(T|damaged) × P(damage) / P(T)
     *
     * @param seq DNA sequence
     * @param sample_profile Sample-level damage profile (for d_max, lambda)
     * @return P(read damaged) in [0, 1], based on terminal base composition
     */
    static float compute_per_read_damage_prior(
        const std::string& seq,
        const SampleDamageProfile& sample_profile);

    /**
     * Compute damage probability for a stop codon at given position
     *
     * Checks if the stop codon could have arisen from C→T damage:
     * - TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
     *
     * Uses exponential decay model: p(damage) = d_max * exp(-lambda * min_dist)
     * where min_dist is distance to nearest terminal.
     *
     * @param codon The 3-nucleotide codon (must be TAA, TAG, or TGA)
     * @param codon_nt_start Position of first nucleotide of codon in oriented sequence
     * @param seq_len Total sequence length
     * @param sample_profile Sample-level damage profile
     * @return Probability that this stop arose from damage (0.0 if not convertible)
     */
    static float compute_stop_damage_probability(
        const char* codon,
        size_t codon_nt_start,
        size_t seq_len,
        const SampleDamageProfile& sample_profile);

    /**
     * Per-read damage evidence from Bayesian log-odds scoring
     */
    struct PerReadDamageEvidence {
        float score_logit = 0.0f;
        float p_read_damaged = 0.0f;
        float mu_aa_prior = 0.0f;
        float logbf_terminal = 0.0f;
        float logbf_stop = 0.0f;
        float logbf_prestop = 0.0f;
        uint16_t n_informative = 0;
        uint8_t n_conv_stops = 0;
        uint8_t n_prestops = 0;
    };

    /**
     * Bayesian per-read damage scoring using AA-impact prior
     *
     * Uses only the selected ORF(s) to compute damage evidence:
     * 1. AA-impact prior: expected non-synonymous damage events in terminal codons
     * 2. Terminal BF: base composition evidence from integrated Bayes factor
     * 3. Convertible stop evidence: TAA/TAG/TGA from CAA/CAG/CGA near terminals
     * 4. Pre-stop counter-evidence: surviving CAA/CAG/CGA codons
     *
     * @param seq DNA sequence
     * @param orfs Top ORFs from frame selection (best first)
     * @param profile Sample-level damage profile
     * @param evidence Optional output struct for diagnostic detail
     * @return P(read has damage) in [0, 1]
     */
    static float infer_per_read_aa_damage(
        const std::string& seq,
        const std::vector<ORFFragment>& orfs,
        const SampleDamageProfile& profile,
        PerReadDamageEvidence* evidence = nullptr);
};

// Free functions for hexamer-based scoring

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
