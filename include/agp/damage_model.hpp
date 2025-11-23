#pragma once

#include "types.hpp"
#include <string>
#include <optional>
#include <array>
#include <mutex>

namespace agp {

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
     * Create damage model by estimating from sequences
     * This infers damage patterns from codon frequency anomalies
     */
    static DamageModel estimate_from_sequences(
        const std::vector<std::string>& sequences,
        const std::vector<std::vector<State>>& state_paths);

    /**
     * Create a damage profile for a specific sequence
     */
    DamageProfile create_profile(size_t seq_len) const;

    /**
     * Create a damage profile with caching for high-throughput processing
     * Thread-safe for concurrent access
     */
    const DamageProfile& create_profile_cached(size_t seq_len);

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

    // Getters
    float lambda_5prime() const { return lambda_5prime_; }
    float lambda_3prime() const { return lambda_3prime_; }
    float delta_max() const { return delta_max_; }
    float delta_background() const { return delta_background_; }

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
