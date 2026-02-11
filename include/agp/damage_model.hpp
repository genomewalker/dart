#pragma once

#include "types.hpp"
#include <string>
#include <optional>
#include <array>
#include <mutex>

namespace agp {

// Forward declaration
class AdaptiveDamageCalibrator;

/**
 * Terminal region damage detection result
 * Focuses on first/last 15bp where damage signal is strongest
 */
struct TerminalDamageResult {
    float damage_score_5prime;      // Damage score from 5' terminal region (0-1)
    float damage_score_3prime;      // Damage score from 3' terminal region (0-1)
    float combined_damage_score;    // Combined damage score (0-1)
    int damaged_hexamers_5prime;    // Count of high-probability damage hexamers at 5'
    int damaged_hexamers_3prime;    // Count of high-probability damage hexamers at 3'
    int total_hexamers_5prime;      // Total hexamers analyzed at 5'
    int total_hexamers_3prime;      // Total hexamers analyzed at 3'
    bool is_likely_damaged;         // True if damage score exceeds threshold
};

/**
 * Ancient DNA damage model
 *
 * Implements position-dependent damage patterns characteristic of ancient DNA:
 * - C->T deamination at 5' ends
 * - G->A deamination at 3' ends (complement of C->T on reverse strand)
 * - Exponential decay with distance from fragment ends
 */
class DamageModel {
public:
    DamageModel() = default;

    // Copy constructor (mutex is not copyable, so we recreate cache)
    DamageModel(const DamageModel& other)
        : lambda_5prime_(other.lambda_5prime_)
        , lambda_3prime_(other.lambda_3prime_)
        , delta_max_(other.delta_max_)
        , delta_background_(other.delta_background_)
        , empirical_ct_5prime_(other.empirical_ct_5prime_)
        , empirical_ga_3prime_(other.empirical_ga_3prime_)
    {
        // Profile cache is not copied - will be rebuilt on demand
    }

    // Move constructor
    DamageModel(DamageModel&& other) noexcept
        : lambda_5prime_(other.lambda_5prime_)
        , lambda_3prime_(other.lambda_3prime_)
        , delta_max_(other.delta_max_)
        , delta_background_(other.delta_background_)
        , empirical_ct_5prime_(std::move(other.empirical_ct_5prime_))
        , empirical_ga_3prime_(std::move(other.empirical_ga_3prime_))
    {
        // Profile cache is not moved - will be rebuilt on demand
    }

    // Copy assignment
    DamageModel& operator=(const DamageModel& other) {
        if (this != &other) {
            lambda_5prime_ = other.lambda_5prime_;
            lambda_3prime_ = other.lambda_3prime_;
            delta_max_ = other.delta_max_;
            delta_background_ = other.delta_background_;
            empirical_ct_5prime_ = other.empirical_ct_5prime_;
            empirical_ga_3prime_ = other.empirical_ga_3prime_;
            // Clear cache - will be rebuilt on demand
            std::lock_guard<std::mutex> lock(cache_mutex_);
            profile_cache_ = {};
        }
        return *this;
    }

    // Move assignment
    DamageModel& operator=(DamageModel&& other) noexcept {
        if (this != &other) {
            lambda_5prime_ = other.lambda_5prime_;
            lambda_3prime_ = other.lambda_3prime_;
            delta_max_ = other.delta_max_;
            delta_background_ = other.delta_background_;
            empirical_ct_5prime_ = std::move(other.empirical_ct_5prime_);
            empirical_ga_3prime_ = std::move(other.empirical_ga_3prime_);
            // Clear cache - will be rebuilt on demand
            std::lock_guard<std::mutex> lock(cache_mutex_);
            profile_cache_ = {};
        }
        return *this;
    }

    /**
     * Create damage model from parameters
     */
    static DamageModel from_parameters(float lambda_5prime,
                                      float lambda_3prime,
                                      float delta_max,
                                      float delta_background);

    /**
     * Update damage model from sample-wide damage profile (Pass 1 results)
     * This ensures Pass 2 predictions use the actual observed damage patterns
     */
    void update_from_sample_profile(const struct SampleDamageProfile& profile);

    /**
     * Create a damage profile for a specific sequence
     */
    DamageProfile create_profile(size_t seq_len) const;

    /**
     * Create a damage profile with caching for high-throughput processing
     * Thread-safe for concurrent access
     */
    const DamageProfile& create_profile_cached(size_t seq_len) const;

    /**
     * Calculate the likelihood that a sequence exhibits ancient DNA damage
     */
    Score calculate_ancient_likelihood(const Sequence& seq) const;

    /**
     * Check if damage patterns are consistent with ancient DNA
     */
    bool is_ancient_compatible() const;

    /**
     * Correct damage in a DNA sequence
     * Reverses likely C->T at 5' end and G->A at 3' end
     * @param seq The damaged DNA sequence
     * @param confidence_threshold Only correct positions with damage probability above this
     * @return Pair of (corrected sequence, number of corrections made)
     */
    std::pair<std::string, size_t> correct_damage(
        const std::string& seq,
        float confidence_threshold = 0.3f) const;

    /**
     * Frame-aware damage correction
     * Uses codon position to adjust correction thresholds:
     * - Wobble (pos 3): more aggressive (lower threshold)
     * - Non-wobble (pos 1,2): more conservative
     *
     * @param seq DNA sequence to correct
     * @param confidence_threshold Base threshold for correction
     * @param frame Reading frame offset (0, 1, or 2), or -1 if unknown
     * @return Pair of (corrected sequence, number of corrections made)
     */
    std::pair<std::string, size_t> correct_damage(
        const std::string& seq,
        float confidence_threshold,
        int frame) const;

    /**
     * Correct damage in a protein sequence by correcting the underlying DNA
     * @param dna_seq The damaged DNA sequence
     * @param protein_seq The damaged protein sequence (for validation)
     * @return Pair of (corrected protein, number of amino acid changes)
     */
    std::pair<std::string, size_t> correct_protein_damage(
        const std::string& dna_seq,
        const std::string& protein_seq,
        float confidence_threshold = 0.3f) const;

    /**
     * Dual-strand damage correction for dsDNA libraries
     *
     * For double-stranded library preparation, reads come from both strands:
     * - Forward strand reads: C→T at 5', G→A at 3'
     * - Reverse strand reads: G→A at 5', C→T at 3'
     *
     * This method tries both correction patterns and returns the one that
     * produces better protein quality (fewer stop codons, better codon usage).
     *
     * @param seq DNA sequence to correct
     * @param confidence_threshold Base threshold for correction
     * @param frame Reading frame offset (0, 1, or 2), or -1 if unknown
     * @return Tuple of (corrected sequence, number of corrections, chosen pattern: 'F'=forward, 'R'=reverse)
     */
    std::tuple<std::string, size_t, char> correct_damage_dual(
        const std::string& seq,
        float confidence_threshold,
        int frame) const;

    /**
     * Score a DNA sequence based on protein quality
     * Used to compare different correction strategies
     * Lower score is better (fewer problems)
     */
    static float score_protein_quality(const std::string& dna, int frame);

    /**
     * Sample-profile-aware damage correction
     *
     * Uses the sample-wide damage profile from Phase 1 to make better correction
     * decisions, and scales the correction threshold by the per-read damage_signal.
     *
     * @param seq DNA sequence to correct
     * @param sample_profile Sample-wide damage profile from Phase 1
     * @param damage_signal Per-read damage probability (0-1)
     * @param frame Reading frame offset (0, 1, or 2), or -1 if unknown
     * @return Pair of (corrected sequence, number of corrections made)
     */
    std::pair<std::string, size_t> correct_damage_with_profile(
        const std::string& seq,
        const SampleDamageProfile& sample_profile,
        float damage_signal,
        int frame,
        bool allow_hexamer_correction = true) const;

    /**
     * Bayesian stop codon correction
     *
     * Uses sample-wide damage patterns to fix stop codons that are likely
     * damage artifacts. Stop codons in coding regions are rare, so when we
     * see one at a position with high damage probability, it's likely spurious.
     *
     * Damage-induced stop codons:
     * - C→T: CAA→TAA, CAG→TAG, CGA→TGA (first position)
     * - G→A: TGG→TGA (third position)
     *
     * For dsDNA libraries, we also consider reverse patterns:
     * - G→A at 5': GAA→TAA (complement), GGA→TGA
     * - C→T at 3': TGC→TGA (complement)
     *
     * @param seq DNA sequence
     * @param seq_len Sequence length
     * @param frame Reading frame (0, 1, 2)
     * @return Pair of (corrected sequence, number of stop codons fixed)
     */
    std::pair<std::string, size_t> correct_stop_codons_bayesian(
        const std::string& seq,
        int frame) const;

    /**
     * Calculate probability that a stop codon at given position is damage-induced
     * Uses Bayesian inference with position-dependent damage rates
     *
     * @param codon The observed stop codon (TAA, TAG, TGA)
     * @param codon_start Position of codon start in sequence
     * @param seq_len Total sequence length
     * @return Probability (0-1) that this stop is damage-induced
     */
    float stop_codon_damage_probability(
        const std::string& codon,
        size_t codon_start,
        size_t seq_len) const;

    /**
     * Calculate amino acid probability distribution for a codon
     * considering position-dependent damage rates.
     *
     * For each observed codon, calculates P(AA | observed_codon, position, damage)
     * by considering all possible original codons that could have been damaged
     * to produce the observed sequence.
     *
     * @param codon Observed codon (3 nucleotides)
     * @param codon_start Position of codon in sequence
     * @param seq_len Total sequence length
     * @return Array of 21 probabilities (20 amino acids + stop)
     *         Indices: A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9,
     *                  M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19, *=20
     */
    std::array<float, 21> calculate_aa_probabilities(
        const std::string& codon,
        size_t codon_start,
        size_t seq_len) const;

    /**
     * Calculate amino acid probability distribution (allocation-free overload)
     * Takes 3 chars directly instead of string to avoid substr allocation in hot loops.
     */
    std::array<float, 21> calculate_aa_probabilities(
        char c0, char c1, char c2,
        size_t codon_start,
        size_t seq_len) const;

    /**
     * Translate DNA with damage-aware probability distribution
     * Returns the most likely protein and per-position confidence
     *
     * @param seq DNA sequence
     * @param frame Reading frame (0, 1, 2)
     * @return Pair of (protein sequence, per-position confidence scores)
     */
    std::pair<std::string, std::vector<float>> translate_with_confidence(
        const std::string& seq,
        int frame) const;

    /**
     * Detect damage using GTDB hexamer frequencies in terminal regions
     *
     * Analyzes only the first and last 15bp of a sequence, where ancient DNA
     * damage signal is strongest. Uses GTDB-derived hexamer rarity to identify
     * hexamers that are likely damage artifacts (rare in real coding sequences
     * but common when their undamaged versions exist).
     *
     * Key insight: Hexamers starting with TGA, TAG, TAA (stop codons) are
     * essentially absent from real coding sequences. Their presence in the
     * terminal 15bp strongly indicates C→T damage at the first position.
     *
     * @param seq DNA sequence to analyze
     * @param terminal_length Maximum positions from each end to analyze (default 15)
     * @return TerminalDamageResult with damage scores and hexamer counts
     */
    static TerminalDamageResult detect_terminal_damage(
        const std::string& seq,
        size_t terminal_length = 15);

    /**
     * Enhanced stop codon damage probability using GTDB hexamer frequencies
     *
     * Combines position-dependent damage rates with GTDB hexamer rarity
     * to calculate probability that a stop codon is damage-induced.
     * Hexamers that are rare in GTDB (like TGA*) get higher damage probability.
     *
     * @param codon The observed stop codon (TAA, TAG, TGA)
     * @param codon_start Position of codon start in sequence
     * @param seq Full sequence (for hexamer context)
     * @return Probability (0-1) that this stop is damage-induced
     */
    float stop_codon_damage_probability_gtdb(
        const std::string& codon,
        size_t codon_start,
        const std::string& seq) const;

    /**
     * Adaptive damage correction using calibrator
     *
     * Uses the AdaptiveDamageCalibrator to make smarter correction decisions:
     * - Per-read damage LLR gating (skip low-signal reads)
     * - Correction validation (reject corrections that worsen quality)
     * - Adaptive thresholds based on batch statistics
     * - Maximum corrections per end limit
     *
     * IMPORTANT: For reverse-strand genes, the damage profile is inverted because
     * gene.sequence is the reverse complement of the original read. The damage
     * patterns are:
     * - Forward strand: C→T at 5' end, G→A at 3' end of gene.sequence
     * - Reverse strand: G→A at 5' end, C→T at 3' end of gene.sequence
     *   (because the 5' end of gene.sequence is the 3' end of the original read)
     *
     * @param seq DNA sequence to correct
     * @param sample_profile Sample-wide damage profile from Phase 1
     * @param damage_signal Per-read damage probability (0-1)
     * @param frame Reading frame offset (0, 1, or 2), or -1 if unknown
     * @param calibrator Adaptive damage calibrator with per-sample settings
     * @param is_forward True for forward-strand genes, false for reverse-strand
     * @return Pair of (corrected sequence, number of corrections made)
     */
    std::pair<std::string, size_t> correct_damage_adaptive(
        const std::string& seq,
        const SampleDamageProfile& sample_profile,
        float damage_signal,
        int frame,
        const AdaptiveDamageCalibrator& calibrator,
        bool is_forward = true) const;

    /**
     * Pre-frame-selection damage correction
     *
     * Applies conservative damage correction BEFORE frame selection to restore
     * the hexamer signal. This is simpler than full correction:
     * - Does not require frame knowledge
     * - Uses only position-based damage probability
     * - Applies stricter Bayesian posterior threshold
     * - Corrects C→T at 5' end and G→A at 3' end only
     *
     * The goal is to improve frame selection accuracy by restoring damaged
     * hexamers, not to achieve perfect correction.
     *
     * @param seq Input DNA sequence (may be damaged)
     * @param sample_profile Sample-wide damage statistics
     * @return Pair of (corrected sequence, number of corrections made)
     */
    static std::pair<std::string, size_t> correct_damage_preframe(
        const std::string& seq,
        const SampleDamageProfile& sample_profile);

    // Getters
    float lambda_5prime() const { return lambda_5prime_; }
    float lambda_3prime() const { return lambda_3prime_; }
    float delta_max() const { return delta_max_; }
    float delta_background() const { return delta_background_; }

    // Convenience getters for terminal damage rates (for positional hexamer weighting)
    // Returns estimated 5' terminal damage rate (typically delta_max for C→T)
    float get_delta_5() const { return delta_max_; }
    // Returns estimated 3' terminal damage rate (typically delta_max for G→A)
    float get_delta_3() const { return delta_max_; }

private:
    float lambda_5prime_ = 0.1f;        // Decay rate from 5' end
    float lambda_3prime_ = 0.1f;        // Decay rate from 3' end
    float delta_max_ = 0.4f;            // Maximum damage rate
    float delta_background_ = 0.01f;    // Background damage rate

    // Estimated from data
    std::optional<std::vector<float>> empirical_ct_5prime_;
    std::optional<std::vector<float>> empirical_ga_3prime_;

    // Profile cache for high-throughput processing (max 2048 different lengths)
    mutable std::array<std::optional<DamageProfile>, 2048> profile_cache_;
    mutable std::mutex cache_mutex_;  // For thread-safe caching
};

/**
 * Damage pattern analyzer
 *
 * Analyzes sequences to detect and quantify ancient DNA damage patterns
 */
class DamageAnalyzer {
public:
    /**
     * Analyze a set of sequences for damage patterns
     * Returns estimated damage parameters
     */
    static DamageModel analyze(const std::vector<std::string>& sequences);

    /**
     * Calculate C->T frequency at each position from 5' end
     */
    static std::vector<float> calculate_ct_profile(
        const std::vector<std::string>& sequences,
        size_t max_positions = 25);

    /**
     * Calculate G->A frequency at each position from 3' end
     */
    static std::vector<float> calculate_ga_profile(
        const std::vector<std::string>& sequences,
        size_t max_positions = 25);

    /**
     * Fit exponential decay model to empirical damage profile
     * Returns (lambda, delta_max, delta_background)
     */
    static std::tuple<float, float, float> fit_decay_model(
        const std::vector<float>& damage_profile);

    /**
     * Detect anomalous codon usage patterns indicative of damage
     * This can work without a damage file by detecting:
     * - Excess stop codons at fragment ends
     * - Systematic amino acid substitutions (e.g., Arg->Stop)
     */
    static DamageModel infer_from_codon_usage(
        const std::vector<std::string>& coding_sequences);

private:
    // Statistical helpers
    static float calculate_entropy(const std::vector<float>& probs);
    static float chi_square_test(const std::vector<float>& observed,
                                 const std::vector<float>& expected);
};

/**
 * Codon damage effects calculator
 *
 * Analyzes how C->T and G->A mutations affect codons:
 * - Creates stop codons (e.g., CGA->TGA, CAA->TAA)
 * - Changes amino acids (synonymous vs non-synonymous)
 * - Affects start codons
 */
class CodonDamageAnalyzer {
public:
    /**
     * Calculate the probability that damage creates a stop codon
     * given the original codon and damage rates
     */
    static Score stop_codon_probability(const std::string& codon,
                                       float ct_rate, float ga_rate);

    /**
     * Get all possible damaged versions of a codon
     * Returns vector of (damaged_codon, probability) pairs
     */
    static std::vector<std::pair<std::string, float>> get_damaged_variants(
        const std::string& codon,
        float ct_rate,
        float ga_rate);

    /**
     * Calculate amino acid substitution probabilities due to damage
     * Returns matrix[original_aa][observed_aa] = probability
     */
    static std::array<std::array<float, 21>, 21> calculate_aa_substitution_matrix(
        float ct_rate, float ga_rate);

    /**
     * Identify codons most susceptible to damage-induced stop codons
     * Returns list of (codon, stop_probability) pairs, sorted by probability
     */
    static std::vector<std::pair<std::string, float>> most_vulnerable_codons();

private:
    // Apply damage to a single nucleotide
    static std::vector<std::pair<char, float>> apply_damage_to_nt(
        char nt, float ct_rate, float ga_rate);
};

} // namespace agp
