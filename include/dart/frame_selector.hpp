#pragma once

#include "types.hpp"
#include "hexamer_tables.hpp"  // For Domain enum
#include <dart/sample_damage_profile.hpp>
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <utility>

namespace dart {

// Forward declarations
class DamageModel;
struct UnifiedDamageContext;

// SampleDamageProfile, DamageValidationState, get_damage_validation_state(),
// get_damage_suppression_factor(), JointDamageModel, MixtureDamageModel are
// provided by <dart/sample_damage_profile.hpp> (libdart-damage).


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

} // namespace dart
