/**
 * Frame selection for ancient DNA reads
 */

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include "agp/scoring.hpp"
#include "agp/simd_utils.hpp"
#include "agp/damage_model.hpp"
#include "agp/hexamer_tables.hpp"
#include "agp/damage_likelihood.hpp"
#include "agp/positional_hexamer.hpp"
#include "agp/dicodon_phase.hpp"  // Frame-discriminative hexamer scoring
#include "agp/forbidden_patterns.hpp"  // Forbidden pattern scoring for frame discrimination
#include "agp/gtdb_strand_hexamer.hpp"  // Sense/antisense hexamer LLR for strand discrimination
#include "agp/strand_scoring.hpp"        // Strand LLR scoring functions
#include "agp/frame_nn.hpp"  // Neural network for frame classification
#include "agp/bayesian_frame.hpp"  // Bayesian frame selection using hexamer LLR
#include <cmath>
#include <algorithm>
#include <climits>

namespace agp {

// Forward declaration for damage-adjusted stop penalty
static float calculate_damage_adjusted_stop_multiplier(
    const std::string& seq, const std::string& protein, int frame, const DamageModel& damage_model);

// ============================================================================
// Phase 1.3: Poisson model for stop codons
// Replace discrete 0/1/2+ buckets with length-normalized likelihood
// ============================================================================

// Count internal stop codons in a protein (excluding terminal stop)
static int count_internal_stops_protein(const std::string& protein) {
    int stops = 0;
    for (size_t i = 0; i + 1 < protein.length(); ++i) {
        if (protein[i] == '*') stops++;
    }
    return stops;
}

// Log-factorial using lgamma for Poisson calculation
static float log_factorial(int k) {
    if (k <= 1) return 0.0f;
    return std::lgamma(static_cast<float>(k + 1));
}

// Poisson log-probability: log P(k | lambda)
static float poisson_log_p(int k, float lambda) {
    if (lambda < 1e-10f) lambda = 1e-10f;
    return k * std::log(lambda) - lambda - log_factorial(k);
}

// Expected stop rate per amino acid under wrong-frame model
// In random sequence: P(stop) = 3/64 ≈ 0.047 (TAA, TAG, TGA)
// In coding but wrong frame: slightly lower due to codon bias
static constexpr float WRONG_FRAME_STOP_RATE = 0.035f;

// Expected stop rate per amino acid under correct-frame model
// Real coding sequences have ~0.1% internal stops (errors, selenocysteine, etc.)
static constexpr float CORRECT_FRAME_STOP_RATE = 0.001f;

// Calculate Poisson-based stop likelihood ratio
// Returns log-likelihood ratio: positive = more likely correct frame
static float calculate_poisson_stop_llr(const std::string& protein) {
    int k = count_internal_stops_protein(protein);
    int n = static_cast<int>(protein.length());
    if (n < 3) return 0.0f;

    float lambda_correct = n * CORRECT_FRAME_STOP_RATE;
    float lambda_wrong = n * WRONG_FRAME_STOP_RATE;

    // Ensure minimum lambda to avoid log(0)
    lambda_correct = std::max(lambda_correct, 0.01f);
    lambda_wrong = std::max(lambda_wrong, 0.1f);

    float log_p_correct = poisson_log_p(k, lambda_correct);
    float log_p_wrong = poisson_log_p(k, lambda_wrong);

    return log_p_correct - log_p_wrong;
}

// Convert Poisson LLR to a 0-1 score compatible with existing framework
// Uses sigmoid to bound but with much larger dynamic range than before
static float poisson_stop_to_score(float llr, size_t seq_len) {
    // Scale sigmoid gain by sequence length - longer sequences get sharper decisions
    // For 100bp (33aa): k=2, for 200bp (66aa): k=3
    float k = 1.5f + seq_len / 150.0f;
    return 1.0f / (1.0f + std::exp(-k * llr / 10.0f));
}

// ============================================================================
// Phase 1.2: Length-adaptive feature weighting
// Adjust weights based on read length to counteract signal convergence
// ============================================================================

struct LengthAdaptiveWeights {
    float stop_weight;
    float hexamer_weight;
    float unorthodox_weight;
    float dicodon_weight;
    float dipeptide_weight;
    float aa_weight;
    float gc3_weight;
    float len_bonus_weight;
    float nn_weight_multiplier;  // Multiplier for NN weight when used
};

// Get length-adaptive weights for frame scoring
// Key insight: for long reads (130-200bp), upweight frame-specific signals
static LengthAdaptiveWeights get_length_adaptive_weights(size_t seq_len) {
    LengthAdaptiveWeights w;

    if (seq_len <= 90) {
        // Short reads: standard weights work well
        w.stop_weight = 0.22f;
        w.hexamer_weight = 0.10f;
        w.unorthodox_weight = 0.18f;
        w.dicodon_weight = 0.20f;
        w.dipeptide_weight = 0.08f;
        w.aa_weight = 0.05f;
        w.gc3_weight = 0.05f;
        w.len_bonus_weight = 0.05f;
        w.nn_weight_multiplier = 1.0f;
    } else if (seq_len <= 130) {
        // Medium reads: slightly boost frame-specific signals
        w.stop_weight = 0.24f;
        w.hexamer_weight = 0.12f;
        w.unorthodox_weight = 0.20f;
        w.dicodon_weight = 0.22f;
        w.dipeptide_weight = 0.06f;
        w.aa_weight = 0.04f;
        w.gc3_weight = 0.04f;
        w.len_bonus_weight = 0.03f;
        w.nn_weight_multiplier = 1.2f;
    } else if (seq_len <= 200) {
        // Long reads (critical zone): heavily boost frame-specific signals
        // This is where accuracy drops from 83% to 57%
        w.stop_weight = 0.28f;        // +6% vs short (stop is most discriminative)
        w.hexamer_weight = 0.14f;     // +4% (hexamer is frame-specific)
        w.unorthodox_weight = 0.24f;  // +6% (unorthodox signals grow with length)
        w.dicodon_weight = 0.26f;     // +6% (dicodon phase is frame-specific)
        w.dipeptide_weight = 0.04f;   // -4% (dipeptide converges)
        w.aa_weight = 0.02f;          // -3% (aa composition converges)
        w.gc3_weight = 0.02f;         // -3% (gc3 converges)
        w.len_bonus_weight = 0.00f;   // Remove (always maxed for long reads)
        w.nn_weight_multiplier = 1.5f; // NN may help more for long reads
    } else {
        // Very long reads (200+bp): similar to long, slightly less aggressive
        w.stop_weight = 0.26f;
        w.hexamer_weight = 0.13f;
        w.unorthodox_weight = 0.22f;
        w.dicodon_weight = 0.24f;
        w.dipeptide_weight = 0.05f;
        w.aa_weight = 0.03f;
        w.gc3_weight = 0.03f;
        w.len_bonus_weight = 0.01f;
        w.nn_weight_multiplier = 1.3f;
    }

    return w;
}

// ============================================================================
// Damage-aware weight adjustment for long reads
// Key insight: In damaged ancient DNA, C→T damage creates stops in TRUE frame
// (CAA→TAA, CAG→TAG, CGA→TGA). This makes stops LESS reliable for long reads.
//
// Analysis of 130-200bp wrong predictions showed:
// - TRUE frame often has MORE stops (40.8% of cases) due to damage
// - Wrong and correct predictions have identical damage rates (~1.1%)
// - The problem is stop penalty being too dominant for long reads
//
// Solution: For long reads, aggressively reduce stop_weight since stops
// are unreliable due to damage. Boost hexamer and dicodon signals instead.
// ============================================================================
static LengthAdaptiveWeights get_damage_aware_weights(size_t seq_len, float damage_rate) {
    LengthAdaptiveWeights w = get_length_adaptive_weights(seq_len);

    // Only apply for long reads (130-200bp)
    if (seq_len < 130) {
        return w;
    }

    // For long reads: nearly eliminate stop weight
    // Ancient DNA damage creates stops in TRUE frame (CAA→TAA, CAG→TAG, CGA→TGA)
    // Analysis showed: for wrong predictions, TRUE frame has MORE stops 40.8% of time
    // Stop codons are fundamentally unreliable for damaged long reads
    // Instead, rely on hexamer and dicodon signals which are frame-specific
    float stop_reduction;
    if (seq_len >= 130) {
        // Long reads: nearly eliminate stop weight (90% reduction)
        // For 130-200bp reads, TRUE frame min stops only 23.6% of wrong predictions
        stop_reduction = 0.1f;
    } else {
        stop_reduction = 1.0f;  // No change for shorter reads
    }

    float original_stop = w.stop_weight;
    w.stop_weight *= stop_reduction;

    // Redistribute reduced stop weight to more reliable signals
    float redistributed = original_stop - w.stop_weight;
    w.dicodon_weight += redistributed * 0.45f;    // Dicodon is frame-specific
    w.hexamer_weight += redistributed * 0.35f;    // Hexamer is frame-specific
    w.unorthodox_weight += redistributed * 0.20f; // Unorthodox patterns

    return w;
}

// Thread-local buffer for reverse complement to avoid allocations
static thread_local std::string tl_rc_buffer;

// Compute reverse complement into thread-local buffer using SIMD
static const std::string& reverse_complement_fast(const std::string& seq) {
    tl_rc_buffer.resize(seq.length());

#ifdef USE_AVX512
    simd::reverse_complement_avx512(seq.data(), tl_rc_buffer.data(), seq.length());
#else
    // Scalar fallback
    for (size_t i = 0; i < seq.length(); ++i) {
        tl_rc_buffer[i] = fast_complement(seq[seq.length() - 1 - i]);
    }
#endif
    return tl_rc_buffer;
}

std::string FrameSelector::reverse_complement(const std::string& seq) {
    // Use SIMD-optimized version for large sequences
    if (seq.length() >= 16) {
        std::string rc(seq.length(), ' ');
#ifdef USE_AVX512
        simd::reverse_complement_avx512(seq.data(), rc.data(), seq.length());
#elif defined(USE_AVX2)
        // Scalar fallback for AVX2 (AVX2 RC not implemented)
        for (size_t i = 0; i < seq.length(); ++i) {
            rc[i] = fast_complement(seq[seq.length() - 1 - i]);
        }
#else
        for (size_t i = 0; i < seq.length(); ++i) {
            rc[i] = fast_complement(seq[seq.length() - 1 - i]);
        }
#endif
        return rc;
    }

    // Scalar for short sequences
    std::string rc;
    rc.reserve(seq.length());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        rc += fast_complement(*it);
    }
    return rc;
}

const std::string& FrameSelector::reverse_complement_cached(const std::string& seq) {
    // Use the existing thread-local buffer via reverse_complement_fast
    return reverse_complement_fast(seq);
}

// Pre-computed frame scores for cross-frame relative advantage calculation
struct PrecomputedFrameScores {
    float unorthodox[3];
    float dicodon[3];
    bool valid;
};

// ============================================================================
// Phase 2: Cross-frame relative stop analysis
// Key insight: When all 3 frames have similar stop counts, stops don't
// discriminate between frames. In this case, reduce stop penalty weight
// and boost secondary signals (hexamer, dicodon, unorthodox).
// ============================================================================

struct PrecomputedStopCounts {
    int counts[3];           // Internal stop counts for each frame
    float variance;          // Variance across 3 frames
    float mean;              // Mean stop count
    int min_stops;           // Minimum stops in any frame
    int max_stops;           // Maximum stops in any frame
    bool all_have_stops;     // True if min_stops >= 1
    bool valid;
};

// Pre-compute stop counts for all 3 frames on one strand
static PrecomputedStopCounts precompute_stop_counts(const std::string& seq) {
    PrecomputedStopCounts stops;
    stops.valid = true;
    stops.min_stops = 1000;
    stops.max_stops = 0;

    float sum = 0;
    for (int f = 0; f < 3; ++f) {
        std::string protein = translate_sequence(seq, f);
        stops.counts[f] = count_internal_stops_protein(protein);
        sum += stops.counts[f];
        if (stops.counts[f] < stops.min_stops) stops.min_stops = stops.counts[f];
        if (stops.counts[f] > stops.max_stops) stops.max_stops = stops.counts[f];
    }

    stops.mean = sum / 3.0f;
    stops.all_have_stops = (stops.min_stops >= 1);

    // Calculate variance
    float var_sum = 0;
    for (int f = 0; f < 3; ++f) {
        float diff = stops.counts[f] - stops.mean;
        var_sum += diff * diff;
    }
    stops.variance = var_sum / 3.0f;

    return stops;
}

// Determine if stops provide discriminative signal across frames
static bool stops_discriminate_frames(const PrecomputedStopCounts& stops, size_t seq_len) {
    // If all frames have 0 stops, no discrimination possible
    if (stops.max_stops == 0) return false;

    // If one frame has 0 stops and others have stops, strongly discriminative
    if (stops.min_stops == 0 && stops.max_stops >= 1) return true;

    // Key fix: if max - min >= 1, stops still discriminate
    // A frame with 1 stop should beat a frame with 2+ stops
    // The previous logic incorrectly disabled discrimination when all frames had stops
    if (stops.max_stops - stops.min_stops >= 1) return true;

    // Variance threshold check for edge cases
    float variance_threshold = (seq_len >= 130) ? 0.5f : 0.3f;
    if (stops.variance > variance_threshold) return true;

    // Only return false if all frames have identical stop counts
    return false;
}

// ============================================================================
// HARD STOP CONSTRAINT FOR LONG READS
// ============================================================================
// For reads ≥130bp, hexamer signals converge and become unreliable for frame
// discrimination. Stop codons are the only reliable signal:
// - Coding frames: ~0.1% internal stops
// - Wrong frames: ~3.5% internal stops (35x difference!)
//
// This function enforces a hard constraint: if one frame has significantly
// fewer stops than the current winner, override the decision.
// ============================================================================

// FrameStopInfo struct removed - no longer needed after hard stop constraint removal

// Compute unorthodox and dicodon scores for all 3 frames (called once per strand)
static PrecomputedFrameScores precompute_frame_scores(const std::string& seq) {
    PrecomputedFrameScores scores;
    scores.valid = true;

    for (int f = 0; f < 3; ++f) {
        scores.unorthodox[f] = forbidden::calculate_unorthodox_frame_score(seq.c_str(), seq.length(), f);
        if (dicodon::has_dicodon_data(active_domain())) {
            scores.dicodon[f] = dicodon::calculate_frame_score(seq, f, active_domain());
        } else {
            scores.dicodon[f] = 0.0f;
        }
    }
    return scores;
}

// Optimized frame scoring - no vector allocations
// Now uses unified positional hexamer scoring with damage awareness
// Phase 2: Added precomputed_stops for cross-frame relative stop analysis
static float score_frame_fast(
    const std::string& seq,
    int frame,
    std::string& protein_out,
    const PrecomputedFrameScores* precomputed = nullptr,
    const PrecomputedStopCounts* precomputed_stops = nullptr,
    float delta_5 = 0.0f,
    float delta_3 = 0.0f) {

    if (seq.length() < static_cast<size_t>(3 + frame)) {
        return -1000.0f;
    }

    // Translate directly
    protein_out = translate_sequence(seq, frame);

    if (protein_out.empty()) {
        return -1000.0f;
    }

    // Use damage-adjusted stop penalty to discount potential C→T damage stops
    // The hexamer-based comparison (CAA/CAG/CGA vs TAA/TAG/TGA likelihood) works
    // independently of sample damage level, so always apply it for long reads
    // where stop codon discrimination is critical.
    float stop;
    if (seq.length() >= 130 || delta_5 > 0.01f || delta_3 > 0.01f) {
        // Create minimal damage model for stop calculation
        DamageModel damage_model = DamageModel::from_parameters(0.1f, 0.1f,
            std::max(delta_5, delta_3) > 0.01f ? std::max(delta_5, delta_3) : 0.05f, 0.01f);
        stop = calculate_damage_adjusted_stop_multiplier(seq, protein_out, frame, damage_model);
    } else {
        stop = FrameSelector::calculate_stop_penalty(protein_out);
    }

    float aa_comp = FrameSelector::calculate_aa_composition_score(protein_out);
    float dipeptide = calculate_dipeptide_score(protein_out);

    // Unified positional hexamer score (replaces calculate_dicodon_score + positional_bonus)
    // Uses position-specific tables (START/INTERNAL/END) and damage-weighted scoring
    float hexamer = positional::calculate_hexamer_score_damage_aware(seq, frame, delta_5, delta_3);

    // Use pre-computed scores if available (3x faster - avoids redundant computation)
    float unorthodox_0, unorthodox_1, unorthodox_2;
    float dicodon_0, dicodon_1, dicodon_2;

    if (precomputed && precomputed->valid) {
        unorthodox_0 = precomputed->unorthodox[0];
        unorthodox_1 = precomputed->unorthodox[1];
        unorthodox_2 = precomputed->unorthodox[2];
        dicodon_0 = precomputed->dicodon[0];
        dicodon_1 = precomputed->dicodon[1];
        dicodon_2 = precomputed->dicodon[2];
    } else {
        // Fallback: compute on demand (for callers that don't pre-compute)
        unorthodox_0 = forbidden::calculate_unorthodox_frame_score(seq.c_str(), seq.length(), 0);
        unorthodox_1 = forbidden::calculate_unorthodox_frame_score(seq.c_str(), seq.length(), 1);
        unorthodox_2 = forbidden::calculate_unorthodox_frame_score(seq.c_str(), seq.length(), 2);
        if (dicodon::has_dicodon_data(active_domain())) {
            dicodon_0 = dicodon::calculate_frame_score(seq, 0, active_domain());
            dicodon_1 = dicodon::calculate_frame_score(seq, 1, active_domain());
            dicodon_2 = dicodon::calculate_frame_score(seq, 2, active_domain());
        } else {
            dicodon_0 = dicodon_1 = dicodon_2 = 0.0f;
        }
    }

    // =========================================================================
    // Phase 1.1: Unsquashed log-likelihood scoring
    // Use clipped linear transforms instead of aggressive sigmoids
    // This preserves discriminative power for long reads where differences are small
    // =========================================================================

    // Compute relative advantage for unorthodox scoring
    float this_unorthodox = (frame == 0) ? unorthodox_0 : (frame == 1) ? unorthodox_1 : unorthodox_2;
    float best_alt_unorthodox = -1000.0f;
    float second_alt_unorthodox = -1000.0f;
    for (int f = 0; f < 3; ++f) {
        if (f == frame) continue;
        float alt = (f == 0) ? unorthodox_0 : (f == 1) ? unorthodox_1 : unorthodox_2;
        if (alt > best_alt_unorthodox) {
            second_alt_unorthodox = best_alt_unorthodox;
            best_alt_unorthodox = alt;
        } else if (alt > second_alt_unorthodox) {
            second_alt_unorthodox = alt;
        }
    }

    float unorthodox_advantage = this_unorthodox - best_alt_unorthodox;

    // Phase 1.1: Use CLIPPED LINEAR instead of sigmoid for long reads
    // For short reads, keep some sigmoid for stability
    // For long reads, let the raw difference speak (clipped to avoid outliers)
    float unorthodox_relative;
    if (seq.length() >= 130) {
        // Long reads: use clipped linear transform to preserve small differences
        // Advantage range is typically [-2, +2], map to [0, 1] with clipping
        float clipped_adv = std::clamp(unorthodox_advantage, -3.0f, 3.0f);
        unorthodox_relative = 0.5f + clipped_adv / 6.0f;  // Maps [-3,3] to [0,1]
    } else {
        // Short reads: use sigmoid but with reduced compression (k=3 instead of 5)
        unorthodox_relative = 1.0f / (1.0f + std::exp(-3.0f * unorthodox_advantage));
    }

    // Compute relative advantage for dicodon phase scoring
    float dicodon_phase_score = 0.0f;
    float dicodon_raw_advantage = 0.0f;  // Keep raw for combination
    if (dicodon::has_dicodon_data(active_domain())) {
        float this_dicodon = (frame == 0) ? dicodon_0 : (frame == 1) ? dicodon_1 : dicodon_2;
        float best_alt_dicodon = -1000.0f;
        for (int f = 0; f < 3; ++f) {
            if (f == frame) continue;
            float alt = (f == 0) ? dicodon_0 : (f == 1) ? dicodon_1 : dicodon_2;
            if (alt > best_alt_dicodon) best_alt_dicodon = alt;
        }

        dicodon_raw_advantage = this_dicodon - best_alt_dicodon;

        // Phase 1.1: Less aggressive sigmoid for dicodon too
        if (seq.length() >= 130) {
            // Long reads: clipped linear
            float clipped_adv = std::clamp(dicodon_raw_advantage, -2.0f, 2.0f);
            dicodon_phase_score = 0.5f + clipped_adv / 4.0f;  // Maps [-2,2] to [0,1]
        } else {
            // Short reads: mild sigmoid (k=1.5 instead of 2)
            dicodon_phase_score = 1.0f / (1.0f + std::exp(-1.5f * dicodon_raw_advantage));
        }
    }

    // GC3
    int gc3_count = 0, total3 = 0;
    for (size_t i = frame + 2; i < seq.length(); i += 3) {
        char c = fast_upper(seq[i]);
        if (c != 'N') {
            if (c == 'G' || c == 'C') gc3_count++;
            total3++;
        }
    }
    float gc3 = total3 > 0 ? static_cast<float>(gc3_count) / total3 : 0.5f;
    // Smooth sigmoid bonus: centered at 0.5, max bonus 0.1 at gc3=0.7+
    float gc3_bonus = 0.1f / (1.0f + std::exp(-10.0f * (gc3 - 0.5f)));

    float len_bonus = std::min(0.2f, static_cast<float>(protein_out.length()) / 100.0f);

    // Damage-frame consistency score (simple T at 5' / A at 3' check)
    float damage_bonus = 0.0f;
    if (seq.length() > 0 && (delta_5 > 0.01f || delta_3 > 0.01f)) {
        float damage_consistency = FrameSelector::score_damage_frame_consistency(seq, frame);
        bool has_damage_signal = (fast_upper(seq[0]) == 'T') ||
                                  (fast_upper(seq[seq.length()-1]) == 'A');
        if (has_damage_signal) {
            damage_bonus = (damage_consistency - 0.5f) * 0.15f;  // +/- 0.075 max
        }
    }

    // Combine absolute unorthodox score with relative advantage
    // Phase 1.1: For long reads, weight the raw advantage more heavily
    float combined_unorthodox;
    if (seq.length() >= 130) {
        // Long reads: 70% relative advantage (preserves small differences)
        combined_unorthodox = 0.3f * this_unorthodox + 0.7f * unorthodox_relative;
    } else {
        combined_unorthodox = 0.5f * this_unorthodox + 0.5f * unorthodox_relative;
    }

    // =========================================================================
    // Phase 1.2: Length-adaptive feature weighting
    // Phase 1.3: Poisson stop model
    // =========================================================================

    // Get length-adaptive weights, adjusting for damage when present
    // Key insight: For damaged ancient DNA, C→T creates stops in TRUE frame,
    // making stop count less reliable for frame discrimination.
    //
    // For long reads (≥130bp): ALWAYS use damage-aware weights
    // Analysis showed this is where accuracy drops (56% vs 82%), and the
    // fundamental problem is that stops become unreliable with length due to:
    // 1. More positions for damage to hit stop-precursor codons
    // 2. Law of Large Numbers making hexamer signals converge
    // 3. Stop penalty dominating when other signals are weak
    float damage_rate = std::max(delta_5, delta_3);
    LengthAdaptiveWeights w = (seq.length() >= 130)
        ? get_damage_aware_weights(seq.length(), damage_rate)
        : get_length_adaptive_weights(seq.length());

    // Phase 1.3: Use Poisson-based stop scoring instead of discrete buckets
    // This provides much better discrimination for long reads where the old
    // 0/1/2+ bucketing caused all frames to collapse to the same penalty
    float stop_llr = calculate_poisson_stop_llr(protein_out);
    float poisson_stop = poisson_stop_to_score(stop_llr, seq.length());

    // Blend Poisson stop with traditional stop penalty
    // For short reads: traditional works well, use 50/50 blend
    // For long reads: Poisson is more discriminative, use 80/20
    float blended_stop;
    if (seq.length() >= 130) {
        blended_stop = 0.8f * poisson_stop + 0.2f * stop;
    } else {
        blended_stop = 0.5f * poisson_stop + 0.5f * stop;
    }

    // Dynamic weight adjustments when stop doesn't discriminate
    float nn_score = 0.0f;
    float nn_weight = 0.0f;
    float effective_stop_weight = w.stop_weight;
    float effective_unorthodox_weight = w.unorthodox_weight;
    float effective_dicodon_weight = w.dicodon_weight;
    float effective_hexamer_weight = w.hexamer_weight;

    // For long reads: ALWAYS compute frame NN score with damage-aware feature masking
    // Analysis showed that for 130-200bp reads:
    // 1. Hexamer signals converge (LLN)
    // 2. Stop features are ANTI-discriminative: TRUE frame has MORE stops due to damage
    // Solution: Use damage-aware NN that ignores stop features
    bool use_nn_for_long_reads = (seq.length() >= 130);
    if (use_nn_for_long_reads) {
        auto nn_features = frame_nn::extract_features(seq.c_str(), seq.length(), frame);
        // Use damage-aware scoring that masks stop-related features
        // This prevents C→T damage-induced stops from biasing toward wrong frames
        nn_score = frame_nn::score_frame_damage_aware(nn_features);
        nn_weight = 0.40f;  // Higher weight since stop-free NN is more reliable for damaged DNA
    }

    // =========================================================================
    // Phase 2: Cross-frame relative stop analysis
    // Check if Poisson stop provides discrimination WITHIN this frame
    // AND check cross-frame stop variance to see if stops discriminate ACROSS frames
    // =========================================================================
    bool stops_discriminate_single = (std::abs(stop_llr) > 2.0f) || (stop < 0.95f);

    // Cross-frame discrimination: use precomputed stop counts if available
    bool stops_discriminate_cross = true;  // Default to discriminative
    float cross_frame_stop_penalty_scale = 1.0f;

    if (precomputed_stops && precomputed_stops->valid) {
        stops_discriminate_cross = stops_discriminate_frames(*precomputed_stops, seq.length());

        // Only apply penalty scaling when stop counts are truly identical or nearly so
        // Key fix: if max - min >= 1, there's a clear difference - don't scale down
        if (seq.length() >= 130 && precomputed_stops->all_have_stops) {
            int stop_diff = precomputed_stops->max_stops - precomputed_stops->min_stops;
            if (stop_diff == 0) {
                // All frames have identical stop counts - scale down heavily
                cross_frame_stop_penalty_scale = 0.3f;
            }
            // When stop_diff >= 1, keep cross_frame_stop_penalty_scale = 1.0
            // This is the key fix: don't penalize stops when they discriminate!
        }
    }

    bool stops_discriminate = stops_discriminate_single && stops_discriminate_cross;

    if (!stops_discriminate) {
        // No stop discrimination - boost other frame-specific signals
        float stop_reduction = (seq.length() >= 130) ? 0.3f : 0.5f;
        effective_stop_weight *= stop_reduction * cross_frame_stop_penalty_scale;

        // Redistribute weight to frame-discriminative signals
        float redistributed = w.stop_weight * (1.0f - stop_reduction * cross_frame_stop_penalty_scale);
        effective_unorthodox_weight += redistributed * 0.35f;
        effective_dicodon_weight += redistributed * 0.35f;
        effective_hexamer_weight += redistributed * 0.30f;

        // Enhanced dicodon phase weighting based on signal strength
        float dicodon_signal_strength = std::abs(dicodon_raw_advantage);
        if (dicodon_signal_strength > 0.5f) {
            effective_dicodon_weight += 0.08f;
        } else if (dicodon_signal_strength > 0.2f) {
            effective_dicodon_weight += 0.04f;
        }

        // Frame NN when stops don't discriminate (only if not already computed for long reads)
        if (!use_nn_for_long_reads) {
            auto nn_features = frame_nn::extract_features(seq.c_str(), seq.length(), frame);
            nn_score = frame_nn::score_frame(nn_features);
            nn_weight = 0.25f * w.nn_weight_multiplier;
            nn_weight = std::min(nn_weight, 0.40f);
        }
        // For long reads, NN already computed with higher weight (0.35)
    }
    // Removed: no longer apply cross_frame_stop_penalty_scale when stops_discriminate
    // This was incorrectly reducing stop weight even when there was a clear difference

    float base_score = effective_stop_weight * blended_stop +
           w.aa_weight * aa_comp +
           effective_hexamer_weight * hexamer +
           effective_unorthodox_weight * combined_unorthodox +
           effective_dicodon_weight * dicodon_phase_score +
           w.dipeptide_weight * dipeptide +
           w.gc3_weight * gc3_bonus +
           w.len_bonus_weight * len_bonus +
           damage_bonus;

    // Blend with NN score when applicable
    return (1.0f - nn_weight) * base_score + nn_weight * nn_score;
}

FrameScore FrameSelector::score_frame(
    const std::string& seq,
    int frame,
    bool forward,
    const DamageProfile* damage) {

    FrameScore score;
    score.frame = frame;
    score.forward = forward;

    if (seq.length() < static_cast<size_t>(3 + frame)) {
        score.total_score = -1000.0f;
        return score;
    }

    // Get working sequence
    std::string working_seq;
    if (forward) {
        working_seq = seq;
    } else {
        working_seq = reverse_complement(seq);
    }

    // Use optimized scoring
    score.total_score = score_frame_fast(working_seq, frame, score.protein);

    if (score.total_score < -100.0f) {
        return score;
    }

    // Extract individual scores for reporting (optional)
    score.codon_score = 0.0f;  // Deprecated - codon scoring removed
    score.stop_codon_penalty = calculate_stop_penalty(score.protein);
    score.aa_composition_score = calculate_aa_composition_score(score.protein);

    if (damage) {
        score.damage_consistency = calculate_damage_consistency(
            seq, working_seq, frame, forward, *damage);
    } else {
        score.damage_consistency = 0.5f;
    }

    return score;
}

FrameScore FrameSelector::select_best_frame(
    const std::string& seq,
    const DamageProfile* damage) {

    if (seq.length() < 3) {
        FrameScore empty;
        empty.total_score = -1000.0f;
        return empty;
    }

    // Compute reverse complement ONCE using thread-local buffer (no allocation)
    const std::string& rc_seq = reverse_complement_fast(seq);

    // Extract damage parameters from profile if available
    // This enables damage-aware stop scoring which discounts damage-induced stops
    // DamageProfile uses delta_max for both ends (symmetric model)
    float delta_5 = 0.0f, delta_3 = 0.0f;
    if (damage) {
        delta_5 = damage->delta_max;
        delta_3 = damage->delta_max;
    }

    FrameScore best;
    best.total_score = -1000.0f;

    std::string protein;

    // Pre-compute unorthodox and dicodon scores once per strand (3x speedup)
    PrecomputedFrameScores fwd_precomputed = precompute_frame_scores(seq);
    PrecomputedFrameScores rev_precomputed = precompute_frame_scores(rc_seq);

    // Phase 2: Pre-compute stop counts for cross-frame analysis
    PrecomputedStopCounts fwd_stops = precompute_stop_counts(seq);
    PrecomputedStopCounts rev_stops = precompute_stop_counts(rc_seq);

    // Score all 6 frames using optimized function
    // Pass damage parameters to enable damage-aware stop scoring
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        float total = score_frame_fast(seq, frame, protein, &fwd_precomputed, &fwd_stops, delta_5, delta_3);

        if (total > best.total_score) {
            best.total_score = total;
            best.frame = frame;
            best.forward = true;
            best.codon_score = 0.0f;  // Deprecated
            best.stop_codon_penalty = calculate_stop_penalty(protein);
            best.aa_composition_score = calculate_aa_composition_score(protein);
            best.damage_consistency = 0.5f;
            best.protein = std::move(protein);
            protein.clear();
        }

        // Reverse frame
        total = score_frame_fast(rc_seq, frame, protein, &rev_precomputed, &rev_stops, delta_5, delta_3);

        if (total > best.total_score) {
            best.total_score = total;
            best.frame = frame;
            best.forward = false;
            best.codon_score = 0.0f;  // Deprecated
            best.stop_codon_penalty = calculate_stop_penalty(protein);
            best.aa_composition_score = calculate_aa_composition_score(protein);
            best.damage_consistency = 0.5f;
            best.protein = std::move(protein);
            protein.clear();
        }
    }

    return best;
}

// Helper: Calculate hexamer log-likelihood for strand comparison
// Returns sum of log(P_coding / P_uniform) for all hexamers in sequence
static float calculate_hexamer_strand_score(const std::string& seq, int frame) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    float log_prob_sum = 0.0f;
    int hexamer_count = 0;

    // Background frequency (uniform = 1/4096)
    static constexpr float LOG_BACKGROUND = -8.317766f;  // log(1/4096)

    char hexbuf[7];
    hexbuf[6] = '\0';

    for (size_t i = frame; i + 5 < seq.length(); i += 3) {
        hexbuf[0] = fast_upper(seq[i]);
        hexbuf[1] = fast_upper(seq[i + 1]);
        hexbuf[2] = fast_upper(seq[i + 2]);
        hexbuf[3] = fast_upper(seq[i + 3]);
        hexbuf[4] = fast_upper(seq[i + 4]);
        hexbuf[5] = fast_upper(seq[i + 5]);

        // Skip if any N
        if (hexbuf[0] == 'N' || hexbuf[1] == 'N' || hexbuf[2] == 'N' ||
            hexbuf[3] == 'N' || hexbuf[4] == 'N' || hexbuf[5] == 'N') continue;

        uint32_t code = encode_hexamer(hexbuf);
        if (code == UINT32_MAX) continue;

        float freq = get_hexamer_freq(code);  // Uses active domain
        if (freq < 1e-8f) freq = 1e-8f;

        log_prob_sum += std::log(freq) - LOG_BACKGROUND;
        hexamer_count++;
    }

    if (hexamer_count == 0) return 0.0f;
    return log_prob_sum;  // Return total, not average (to weight by length)
}

// Helper: Calculate sense/antisense strand LLR using strand-specific hexamer tables
// Frame-specific version: only considers hexamers at codon boundaries for the given frame
// Returns sum of log2(P_sense / P_antisense) for hexamers at codon positions
// Positive = more likely sense strand, Negative = more likely antisense strand
static float calculate_strand_sense_llr_frame(const std::string& seq, int frame) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    float total_llr = 0.0f;
    int count = 0;

    constexpr size_t BOUNDARY = 30;
    char hexbuf[7];
    hexbuf[6] = '\0';

    // Use position-aware strand LLR tables (START, INTERIOR, STOP regions)
    // Only consider hexamers at codon boundaries (every 3rd position starting from frame)
    for (size_t i = frame; i + 5 < seq.length(); i += 3) {
        hexbuf[0] = fast_upper(seq[i]);
        hexbuf[1] = fast_upper(seq[i + 1]);
        hexbuf[2] = fast_upper(seq[i + 2]);
        hexbuf[3] = fast_upper(seq[i + 3]);
        hexbuf[4] = fast_upper(seq[i + 4]);
        hexbuf[5] = fast_upper(seq[i + 5]);

        // Skip if any N
        if (hexbuf[0] == 'N' || hexbuf[1] == 'N' || hexbuf[2] == 'N' ||
            hexbuf[3] == 'N' || hexbuf[4] == 'N' || hexbuf[5] == 'N') continue;

        uint32_t code = encode_hexamer(hexbuf);
        if (code == UINT32_MAX || code >= 4096) continue;

        // Select table based on position
        float llr;
        if (i < BOUNDARY) {
            llr = STRAND_LLR_START[code];
        } else if (i + 6 > seq.length() - BOUNDARY) {
            llr = STRAND_LLR_STOP[code];
        } else {
            llr = STRAND_LLR_INTERIOR[code];
        }

        total_llr += llr;
        count++;
    }

    return (count > 0) ? total_llr / count : 0.0f;
}

// Count internal stop codons for a given frame
static int count_internal_stops_fast(const std::string& seq, int frame) {
    if (seq.length() < static_cast<size_t>(frame + 3)) return 0;
    
    int stops = 0;
    size_t last_codon_start = ((seq.length() - frame) / 3 - 1) * 3 + frame;
    
    for (size_t i = frame; i + 2 < seq.length() && i < last_codon_start; i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);
        
        // Check for stop codons: TAA, TAG, TGA
        if (c1 == 'T') {
            if ((c2 == 'A' && (c3 == 'A' || c3 == 'G')) ||
                (c2 == 'G' && c3 == 'A')) {
                stops++;
            }
        }
    }
    return stops;
}

// Compute strand decision score using hexamer log-odds and stop codons
// Returns positive if forward is preferred, negative if reverse is preferred
// Also returns confidence (magnitude of decision)
static std::pair<float, float> compute_strand_decision(const std::string& seq, const std::string& rc_seq) {
    // Compute hexamer log-odds for best frame in each orientation
    float best_fwd_hex = -1000.0f;
    float best_rev_hex = -1000.0f;
    int best_fwd_frame = 0, best_rev_frame = 0;
    
    for (int frame = 0; frame < 3; ++frame) {
        float fwd_hex = calculate_hexamer_strand_score(seq, frame);
        float rev_hex = calculate_hexamer_strand_score(rc_seq, frame);
        if (fwd_hex > best_fwd_hex) { best_fwd_hex = fwd_hex; best_fwd_frame = frame; }
        if (rev_hex > best_rev_hex) { best_rev_hex = rev_hex; best_rev_frame = frame; }
    }
    
    float hex_diff = best_fwd_hex - best_rev_hex;  // positive = fwd preferred
    
    // Compute min stop codons for each orientation
    int min_fwd_stops = 1000;
    int min_rev_stops = 1000;
    for (int frame = 0; frame < 3; ++frame) {
        int fwd_stops = count_internal_stops_fast(seq, frame);
        int rev_stops = count_internal_stops_fast(rc_seq, frame);
        if (fwd_stops < min_fwd_stops) min_fwd_stops = fwd_stops;
        if (rev_stops < min_rev_stops) min_rev_stops = rev_stops;
    }
    
    int stop_diff = min_rev_stops - min_fwd_stops;  // positive = fwd has fewer stops
    
    // Combine signals:
    // - Stop codons are very reliable when different (use as strong override)
    // - Hexamer log-odds is the primary continuous signal
    float decision;
    float confidence;
    
    if (std::abs(stop_diff) >= 1) {
        // Stop codons provide definitive signal - 96% accurate when they differ
        decision = stop_diff * 10.0f + hex_diff * 0.3f;
        confidence = std::abs(stop_diff) * 5.0f + std::abs(hex_diff) * 0.5f;
    } else {
        // No stop difference - rely on hexamer log-odds
        decision = hex_diff;
        confidence = std::abs(hex_diff);
    }
    
    return {decision, confidence};
}

std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand(
    const std::string& seq,
    const DamageProfile* damage) {

    if (seq.length() < 3) {
        FrameScore empty;
        empty.total_score = -1000.0f;
        return {empty, empty};
    }

    // Compute reverse complement ONCE using thread-local buffer (no allocation)
    const std::string& rc_seq = reverse_complement_fast(seq);

    FrameScore best_fwd, best_rev;
    best_fwd.total_score = -1000.0f;
    best_fwd.forward = true;
    best_rev.total_score = -1000.0f;
    best_rev.forward = false;

    // Pre-compute best hexamer scores for strand comparison
    // This captures the strand-specific codon usage signal
    float best_fwd_hexamer = -1000.0f;
    float best_rev_hexamer = -1000.0f;
    for (int frame = 0; frame < 3; ++frame) {
        float fwd_hex = calculate_hexamer_strand_score(seq, frame);
        float rev_hex = calculate_hexamer_strand_score(rc_seq, frame);
        if (fwd_hex > best_fwd_hexamer) best_fwd_hexamer = fwd_hex;
        if (rev_hex > best_rev_hexamer) best_rev_hexamer = rev_hex;
    }

    // Calculate strand comparison bonus
    // If one strand has much better hexamer evidence, boost its scores
    float hexamer_diff = best_fwd_hexamer - best_rev_hexamer;
    float strand_bonus_fwd = 0.0f;
    float strand_bonus_rev = 0.0f;

    // Apply sigmoid-like scaling: max bonus +/- 0.15
    // Threshold at diff > 2.0 (about 7x likelihood difference)
    if (hexamer_diff > 0.5f) {
        strand_bonus_fwd = std::min(0.15f, hexamer_diff * 0.05f);
    } else if (hexamer_diff < -0.5f) {
        strand_bonus_rev = std::min(0.15f, -hexamer_diff * 0.05f);
    }

    // =====================================================================
    // SENSE/ANTISENSE STRAND LLR (CODON-ALIGNED)
    // Uses strand-discriminative hexamer tables trained on CDS vs RC(CDS)
    // Training data was extracted at codon boundaries (frame 0) for maximum
    // discriminative power. We evaluate all 3 frames and take the best match.
    // =====================================================================
    float best_sense_fwd = -1000.0f;
    float best_sense_rev = -1000.0f;
    for (int frame = 0; frame < 3; ++frame) {
        float fwd_sense = calculate_strand_sense_llr_frame(seq, frame);
        float rev_sense = calculate_strand_sense_llr_frame(rc_seq, frame);
        if (fwd_sense > best_sense_fwd) best_sense_fwd = fwd_sense;
        if (rev_sense > best_sense_rev) best_sense_rev = rev_sense;
    }
    float sense_diff = best_sense_fwd - best_sense_rev;  // Positive = forward preferred

    // Add sense/antisense signal to strand bonus
    // Weight: 0.10 per unit LLR, max bonus 0.15
    // The sense LLR supplements the hexamer coding potential signal
    if (sense_diff > 0.1f) {
        strand_bonus_fwd += std::min(0.15f, sense_diff * 0.10f);
    } else if (sense_diff < -0.1f) {
        strand_bonus_rev += std::min(0.15f, -sense_diff * 0.10f);
    }

    std::string protein;

    // Extract damage parameters for damage-aware stop discounting
    // Critical: without these, the damage-adjusted stop penalty is never triggered!
    float delta_5 = 0.0f, delta_3 = 0.0f;
    if (damage) {
        delta_5 = damage->delta_max;
        delta_3 = damage->delta_max;
    }

    // Pre-compute unorthodox and dicodon scores once per strand (3x speedup)
    PrecomputedFrameScores fwd_precomputed = precompute_frame_scores(seq);
    PrecomputedFrameScores rev_precomputed = precompute_frame_scores(rc_seq);

    // Phase 2: Pre-compute stop counts for cross-frame analysis
    PrecomputedStopCounts fwd_stops = precompute_stop_counts(seq);
    PrecomputedStopCounts rev_stops = precompute_stop_counts(rc_seq);

    // Track all frame scores for debugging
    float fwd_scores[3] = {-1000.0f, -1000.0f, -1000.0f};
    float rev_scores[3] = {-1000.0f, -1000.0f, -1000.0f};
    std::string fwd_proteins[3];
    std::string rev_proteins[3];

    // Score all 6 frames with damage-aware stop discounting
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame - pass damage parameters!
        float total = score_frame_fast(seq, frame, protein, &fwd_precomputed, &fwd_stops, delta_5, delta_3);
        total += strand_bonus_fwd;
        fwd_scores[frame] = total;
        fwd_proteins[frame] = protein;

        if (total > best_fwd.total_score) {
            best_fwd.total_score = total;
            best_fwd.frame = frame;
            best_fwd.forward = true;
            best_fwd.codon_score = 0.0f;  // Deprecated
            best_fwd.stop_codon_penalty = calculate_stop_penalty(protein);
            best_fwd.aa_composition_score = calculate_aa_composition_score(protein);
            best_fwd.damage_consistency = 0.5f;
            best_fwd.protein = protein;
        }

        // Reverse frame - pass damage parameters!
        total = score_frame_fast(rc_seq, frame, protein, &rev_precomputed, &rev_stops, delta_5, delta_3);
        total += strand_bonus_rev;
        rev_scores[frame] = total;
        rev_proteins[frame] = protein;

        if (total > best_rev.total_score) {
            best_rev.total_score = total;
            best_rev.frame = frame;
            best_rev.forward = false;
            best_rev.codon_score = 0.0f;  // Deprecated
            best_rev.stop_codon_penalty = calculate_stop_penalty(protein);
            best_rev.aa_composition_score = calculate_aa_composition_score(protein);
            best_rev.damage_consistency = 0.5f;
            best_rev.protein = protein;
        }
    }

    // NOTE: Hard stop constraint removed - analysis showed 80.8% of frame errors
    // have TRUE frame with MORE stops (damage-induced C→T: CAA→TAA, CAG→TAG, CGA→TGA).
    // Stop-based override would pick the "cleaner" wrong frame.
    // Instead, use damage-aware stop discounting in score_frame_fast.

    return {best_fwd, best_rev};
}

// Overload with sample profile for damage-aware strand selection
std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand(
    const std::string& seq,
    const DamageProfile* damage,
    const SampleDamageProfile* sample_profile) {

    if (seq.length() < 3) {
        FrameScore empty;
        empty.total_score = -1000.0f;
        return {empty, empty};
    }

    const std::string& rc_seq = reverse_complement_fast(seq);

    FrameScore best_fwd, best_rev;
    best_fwd.total_score = -1000.0f;
    best_fwd.forward = true;
    best_rev.total_score = -1000.0f;
    best_rev.forward = false;

    // Pre-compute best hexamer scores for strand comparison
    float best_fwd_hexamer = -1000.0f;
    float best_rev_hexamer = -1000.0f;
    for (int frame = 0; frame < 3; ++frame) {
        float fwd_hex = calculate_hexamer_strand_score(seq, frame);
        float rev_hex = calculate_hexamer_strand_score(rc_seq, frame);
        if (fwd_hex > best_fwd_hexamer) best_fwd_hexamer = fwd_hex;
        if (rev_hex > best_rev_hexamer) best_rev_hexamer = rev_hex;
    }

    float hexamer_diff = best_fwd_hexamer - best_rev_hexamer;
    float strand_bonus_fwd = 0.0f;
    float strand_bonus_rev = 0.0f;

    if (hexamer_diff > 0.5f) {
        strand_bonus_fwd = std::min(0.15f, hexamer_diff * 0.05f);
    } else if (hexamer_diff < -0.5f) {
        strand_bonus_rev = std::min(0.15f, -hexamer_diff * 0.05f);
    }

    // Sense/antisense strand LLR
    float best_sense_fwd = -1000.0f;
    float best_sense_rev = -1000.0f;
    for (int frame = 0; frame < 3; ++frame) {
        float fwd_sense = calculate_strand_sense_llr_frame(seq, frame);
        float rev_sense = calculate_strand_sense_llr_frame(rc_seq, frame);
        if (fwd_sense > best_sense_fwd) best_sense_fwd = fwd_sense;
        if (rev_sense > best_sense_rev) best_sense_rev = rev_sense;
    }
    float sense_diff = best_sense_fwd - best_sense_rev;

    if (sense_diff > 0.1f) {
        strand_bonus_fwd += std::min(0.15f, sense_diff * 0.10f);
    } else if (sense_diff < -0.1f) {
        strand_bonus_rev += std::min(0.15f, -sense_diff * 0.10f);
    }

    // Damage-based strand preference when sample has significant damage
    if (sample_profile && sample_profile->max_damage_5prime > 0.05f) {
        float damage_strand_pref = calculate_damage_strand_preference(seq, rc_seq, *sample_profile);
        float damage_strand_bonus = std::clamp(damage_strand_pref * 0.15f, -0.15f, 0.15f);
        strand_bonus_fwd += damage_strand_bonus;
        strand_bonus_rev -= damage_strand_bonus;
    }

    std::string protein;

    // Extract damage parameters from sample profile for damage-aware stop discounting
    float delta_5 = 0.0f, delta_3 = 0.0f;
    if (sample_profile) {
        delta_5 = sample_profile->max_damage_5prime;
        delta_3 = sample_profile->max_damage_3prime;
    }

    // Pre-compute unorthodox and dicodon scores once per strand (3x speedup)
    PrecomputedFrameScores fwd_precomputed = precompute_frame_scores(seq);
    PrecomputedFrameScores rev_precomputed = precompute_frame_scores(rc_seq);

    // Phase 2: Pre-compute stop counts for cross-frame analysis
    PrecomputedStopCounts fwd_stops = precompute_stop_counts(seq);
    PrecomputedStopCounts rev_stops = precompute_stop_counts(rc_seq);

    // Track frame scores for debugging
    float fwd_scores[3] = {-1000.0f, -1000.0f, -1000.0f};
    float rev_scores[3] = {-1000.0f, -1000.0f, -1000.0f};
    std::string fwd_proteins[3];
    std::string rev_proteins[3];

    // Score all 6 frames with damage-aware stop discounting
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame - pass damage parameters!
        float total = score_frame_fast(seq, frame, protein, &fwd_precomputed, &fwd_stops, delta_5, delta_3);
        total += strand_bonus_fwd;
        fwd_scores[frame] = total;
        fwd_proteins[frame] = protein;

        if (total > best_fwd.total_score) {
            best_fwd.total_score = total;
            best_fwd.frame = frame;
            best_fwd.forward = true;
            best_fwd.codon_score = 0.0f;  // Deprecated
            best_fwd.stop_codon_penalty = calculate_stop_penalty(protein);
            best_fwd.aa_composition_score = calculate_aa_composition_score(protein);
            best_fwd.damage_consistency = 0.5f;
            best_fwd.protein = protein;
        }

        // Reverse frame - pass damage parameters!
        total = score_frame_fast(rc_seq, frame, protein, &rev_precomputed, &rev_stops, delta_5, delta_3);
        total += strand_bonus_rev;
        rev_scores[frame] = total;
        rev_proteins[frame] = protein;

        if (total > best_rev.total_score) {
            best_rev.total_score = total;
            best_rev.frame = frame;
            best_rev.forward = false;
            best_rev.codon_score = 0.0f;  // Deprecated
            best_rev.stop_codon_penalty = calculate_stop_penalty(protein);
            best_rev.aa_composition_score = calculate_aa_composition_score(protein);
            best_rev.damage_consistency = 0.5f;
            best_rev.protein = protein;
        }
    }

    // NOTE: Hard stop constraint removed - see comment in first select_best_per_strand overload.

    return {best_fwd, best_rev};
}

std::vector<FrameScore> FrameSelector::score_all_frames(
    const std::string& seq,
    const DamageProfile* damage) {

    std::vector<FrameScore> scores;
    scores.reserve(6);

    for (int frame = 0; frame < 3; ++frame) {
        scores.push_back(score_frame(seq, frame, true, damage));
        scores.push_back(score_frame(seq, frame, false, damage));
    }

    std::sort(scores.begin(), scores.end(),
              [](const FrameScore& a, const FrameScore& b) {
                  return a.total_score > b.total_score;
              });

    return scores;
}

float FrameSelector::calculate_wobble_enrichment(
    const std::string& seq,
    int frame,
    size_t terminal_region) {

    if (seq.length() < static_cast<size_t>(frame + 3)) {
        return 0.5f;  // Neutral if too short
    }

    // Limit analysis to terminal region
    size_t analyze_len = std::min(terminal_region, seq.length());

    // Count T's at each codon position within the terminal region
    int t_at_wobble = 0;    // T's at position 3 (wobble)
    int t_at_other = 0;     // T's at positions 1 and 2
    int total_pos3 = 0;     // Total bases at wobble position
    int total_other = 0;    // Total bases at other positions

    for (size_t i = frame; i < analyze_len; ++i) {
        char c = fast_upper(seq[i]);
        if (c == 'N') continue;

        int codon_pos = (i - frame) % 3;  // 0, 1, 2

        if (codon_pos == 2) {  // Wobble position (3rd in codon)
            total_pos3++;
            if (c == 'T') t_at_wobble++;
        } else {  // Position 1 or 2
            total_other++;
            if (c == 'T') t_at_other++;
        }
    }

    // Calculate enrichment: how much more T's at wobble vs other positions
    // Expected: if damage is evenly distributed, T freq should be similar
    // In correct frame: T's from C→T damage more likely at wobble (synonymous)
    if (total_pos3 == 0 || total_other == 0) {
        return 0.5f;
    }

    float t_freq_wobble = static_cast<float>(t_at_wobble) / total_pos3;
    float t_freq_other = static_cast<float>(t_at_other) / total_other;

    // Enrichment ratio normalized to 0-1
    // If t_freq_wobble > t_freq_other, this frame is more likely correct
    // Use logistic function to bound between 0 and 1
    float ratio = (t_freq_wobble + 0.01f) / (t_freq_other + 0.01f);
    float enrichment = 1.0f / (1.0f + std::exp(-(ratio - 1.0f) * 3.0f));

    return enrichment;
}

float FrameSelector::detect_damage_induced_stops(
    const std::string& seq,
    int frame) {

    if (seq.length() < static_cast<size_t>(frame + 3)) {
        return 1.0f;  // No stops = good
    }

    // Check first few codons for potential damage-induced stops
    // TAA could be from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
    int suspicious_stops = 0;
    int codons_checked = 0;

    // Check first 5 codons (15 bp) at 5' end
    for (size_t i = frame; i + 2 < seq.length() && codons_checked < 5; i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);
        codons_checked++;

        // Check for damage-induced stop codons
        // TAA: could be from CAA (C→T damage at position 1)
        // TAG: could be from CAG (C→T damage at position 1)
        // TGA: could be from CGA (C→T damage at position 1)
        if (c1 == 'T') {
            if ((c2 == 'A' && c3 == 'A') ||  // TAA from CAA
                (c2 == 'A' && c3 == 'G') ||  // TAG from CAG
                (c2 == 'G' && c3 == 'A')) {  // TGA from CGA
                // This stop could be damage-induced
                // Weight more heavily if at position 0 (highest damage rate)
                float pos_weight = 1.0f - static_cast<float>(i - frame) / 15.0f;
                suspicious_stops += static_cast<int>(pos_weight * 2 + 0.5f);
            }
        }
    }

    // Return score: 1.0 = no suspicious stops, 0.0 = many suspicious stops
    // Each suspicious stop reduces confidence
    float score = 1.0f / (1.0f + suspicious_stops * 0.5f);
    return score;
}

// Calculate strand preference based on ancient DNA damage asymmetry
// Forward strand: expect T-enrichment at 5' (C→T), A-enrichment at 3' (G→A complement)
// Returns positive for forward preference, negative for reverse
float calculate_damage_strand_preference(
    const std::string& seq,
    const std::string& rc_seq,
    const SampleDamageProfile& profile
) {
    if (profile.max_damage_5prime < 0.05f) {
        return 0.0f;  // Not enough damage signal
    }

    const size_t TERMINAL_BP = 10;
    size_t len = seq.length();
    if (len < TERMINAL_BP * 2) return 0.0f;

    // Count T/(T+C) at 5' and A/(A+G) at 3' for forward strand
    size_t fwd_5_t = 0, fwd_5_c = 0, fwd_3_a = 0, fwd_3_g = 0;
    for (size_t i = 0; i < TERMINAL_BP; ++i) {
        char c5 = fast_upper(seq[i]);
        char c3 = fast_upper(seq[len - 1 - i]);
        if (c5 == 'T') fwd_5_t++;
        else if (c5 == 'C') fwd_5_c++;
        if (c3 == 'A') fwd_3_a++;
        else if (c3 == 'G') fwd_3_g++;
    }

    // Same for reverse complement
    size_t rev_5_t = 0, rev_5_c = 0, rev_3_a = 0, rev_3_g = 0;
    for (size_t i = 0; i < TERMINAL_BP; ++i) {
        char c5 = fast_upper(rc_seq[i]);
        char c3 = fast_upper(rc_seq[len - 1 - i]);
        if (c5 == 'T') rev_5_t++;
        else if (c5 == 'C') rev_5_c++;
        if (c3 == 'A') rev_3_a++;
        else if (c3 == 'G') rev_3_g++;
    }

    // Calculate ratios (with Laplace smoothing)
    float fwd_5_ratio = static_cast<float>(fwd_5_t + 1) / (fwd_5_t + fwd_5_c + 2);
    float fwd_3_ratio = static_cast<float>(fwd_3_a + 1) / (fwd_3_a + fwd_3_g + 2);
    float rev_5_ratio = static_cast<float>(rev_5_t + 1) / (rev_5_t + rev_5_c + 2);
    float rev_3_ratio = static_cast<float>(rev_3_a + 1) / (rev_3_a + rev_3_g + 2);

    // Expected baseline from sample middle positions
    float baseline_tc = 0.5f;  // Neutral expectation
    float baseline_ag = 0.5f;

    // Forward strand should show: high 5' T/(T+C) and high 3' A/(A+G) in damaged DNA
    float fwd_damage_signal = (fwd_5_ratio - baseline_tc) + (fwd_3_ratio - baseline_ag);
    float rev_damage_signal = (rev_5_ratio - baseline_tc) + (rev_3_ratio - baseline_ag);

    // Scale by sample damage level for robustness
    float damage_strength = std::min(1.0f, profile.max_damage_5prime * 5.0f);

    return (fwd_damage_signal - rev_damage_signal) * damage_strength;
}

float FrameSelector::calculate_strand_confidence(
    const std::string& seq,
    float fwd_score,
    float rev_score,
    int fwd_frame,
    int rev_frame,
    const SampleDamageProfile& sample_profile) {

    // Primary signal: Strand LLR from boundary-aware hexamer tables
    // This works for ALL sequences (ancient or modern) by detecting
    // asymmetric patterns: ATG enrichment at 5', stop context at 3'
    float strand_llr_fwd, strand_llr_rev;
    bool has_strand_signal = strand::compute_both_llrs(
        seq.data(), seq.length(), strand_llr_fwd, strand_llr_rev);

    float strand_llr_diff = strand_llr_fwd - strand_llr_rev;
    float strand_conf = strand::llr_to_confidence(strand_llr_fwd, strand_llr_rev);

    // If strand LLR strongly favors one direction, use it directly
    // (LLR diff > 5 means ~30x more likely, very confident)
    if (has_strand_signal && std::abs(strand_llr_diff) > 5.0f) {
        return strand_conf;
    }

    // Secondary signal: Coding score difference
    float score_diff = std::abs(fwd_score - rev_score);
    float score_conf = 0.5f + std::min(0.5f, score_diff * 0.5f);

    // Combine strand LLR and score confidence
    float combined_conf;
    if (has_strand_signal) {
        // Weight strand LLR more heavily when it's informative
        float strand_weight = std::min(1.0f, std::abs(strand_llr_diff) / 10.0f);
        combined_conf = strand_conf * strand_weight + score_conf * (1.0f - strand_weight);
    } else {
        combined_conf = score_conf;
    }

    // Tertiary signal: Damage-based signals (when sample has damage)
    bool has_damage = sample_profile.max_damage_5prime > 0.05f ||
                      sample_profile.max_damage_3prime > 0.05f;

    if (has_damage) {
        // Compute reverse complement for reverse strand analysis
        std::string rc = reverse_complement(seq);

        // Calculate wobble enrichment for both strands
        float wobble_fwd = calculate_wobble_enrichment(seq, fwd_frame);
        float wobble_rev = calculate_wobble_enrichment(rc, rev_frame);

        // Detect damage-induced stops
        float stop_fwd = detect_damage_induced_stops(seq, fwd_frame);
        float stop_rev = detect_damage_induced_stops(rc, rev_frame);

        // Combine damage signals
        float fwd_signal = wobble_fwd * stop_fwd;
        float rev_signal = wobble_rev * stop_rev;

        float signal_ratio = (fwd_signal + 0.01f) / (rev_signal + 0.01f);
        float damage_conf;
        if (fwd_score > rev_score) {
            damage_conf = (signal_ratio > 1.0f) ? 0.5f + 0.5f * (1.0f - 1.0f/signal_ratio) : 0.5f * signal_ratio;
        } else {
            damage_conf = (signal_ratio < 1.0f) ? 0.5f + 0.5f * (1.0f - signal_ratio) : 0.5f / signal_ratio;
        }

        // Weight by damage level
        float damage_weight = std::min(0.3f, (sample_profile.max_damage_5prime + sample_profile.max_damage_3prime) * 0.3f);
        combined_conf = combined_conf * (1.0f - damage_weight) + damage_conf * damage_weight;
    }

    return std::max(0.5f, std::min(1.0f, combined_conf));
}

/**
 * Calculate damage-adjusted stop codon penalty multiplier
 * Returns a multiplier (0.0-1.0) that's more lenient for damage-induced stops
 *
 * This replaces calculate_stop_penalty() for damage-aware frame selection.
 * Instead of treating all internal stops equally, it weights them by damage probability.
 */
static float calculate_damage_adjusted_stop_multiplier(
    const std::string& seq,
    const std::string& protein,
    int frame,
    const DamageModel& damage_model) {

    if (protein.empty()) return 1.0f;

    size_t seq_len = seq.length();

    // Helper to get hexamer and calculate log probability
    auto get_hexamer_log_prob = [](const std::string& s, size_t pos) -> float {
        if (pos + 6 > s.length()) return -10.0f;
        uint32_t code = agp::encode_hexamer(s.c_str() + pos);
        if (code == UINT32_MAX) return -10.0f;
        float freq = agp::get_hexamer_freq(code);  // Uses active domain
        if (freq < 1e-8f) freq = 1e-8f;
        return std::log(freq);
    };

    // Count internal stops and their damage probabilities
    // Key insight: C→T damage creates stops (CAA→TAA, CAG→TAG, CGA→TGA) in TRUE frame.
    // The WRONG frame has different codons at the same positions, so damage creates
    // different amino acids (not necessarily stops). By discounting C→T damage stops,
    // we rescue the TRUE frame's score.
    int total_internal_stops = 0;
    float damage_weighted_stops = 0.0f;

    for (size_t i = 0; i + 1 < protein.length(); ++i) {  // i + 1 to exclude terminal stop
        if (protein[i] == '*') {
            total_internal_stops++;

            // Get codon position in DNA
            size_t codon_start = frame + (i * 3);
            if (codon_start + 3 > seq_len) {
                damage_weighted_stops += 1.0f;
                continue;
            }

            // Get stop codon
            char c1 = fast_upper(seq[codon_start]);
            char c2 = fast_upper(seq[codon_start + 1]);
            char c3 = fast_upper(seq[codon_start + 2]);

            // Check if this is a potential C→T damage stop:
            // TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
            // All stop codons start with T, so check if T could be from C→T
            bool is_taa = (c1 == 'T' && c2 == 'A' && c3 == 'A');
            bool is_tag = (c1 == 'T' && c2 == 'A' && c3 == 'G');
            bool is_tga = (c1 == 'T' && c2 == 'G' && c3 == 'A');

            float stop_weight = 1.0f;  // Default: full penalty

            // For any potential C→T damage stop, check hexamer evidence
            if (is_taa || is_tag || is_tga) {
                // Find hexamer context (codon-aligned)
                // The stop codon is at codon_start, we want the hexamer containing it
                size_t hexamer_start = codon_start;
                if (hexamer_start + 6 > seq_len) {
                    hexamer_start = (seq_len >= 6) ? seq_len - 6 : 0;
                }

                if (hexamer_start + 6 <= seq_len) {
                    // Get observed (stop-containing) hexamer probability
                    float log_p_stop = get_hexamer_log_prob(seq, hexamer_start);

                    // Create corrected hexamer: T→C to restore original codon
                    std::string corrected_hex = seq.substr(hexamer_start, 6);
                    size_t t_pos_in_hex = codon_start - hexamer_start;

                    if (t_pos_in_hex < 6) {
                        corrected_hex[t_pos_in_hex] = 'C';  // T→C: TAA→CAA, TAG→CAG, TGA→CGA

                        float log_p_corrected = get_hexamer_log_prob(corrected_hex, 0);
                        float log_likelihood_ratio = log_p_corrected - log_p_stop;

                        // If corrected hexamer is more likely, this stop is probably damage
                        // Use graduated discounting based on evidence strength:
                        // LLR > 2.0 (~7x more likely): strong evidence of damage
                        // LLR > 1.0 (~3x more likely): moderate evidence
                        // LLR > 0.5 (~1.6x more likely): weak evidence
                        if (log_likelihood_ratio > 2.0f) {
                            stop_weight = 0.3f;  // Strong evidence: heavy discount
                        } else if (log_likelihood_ratio > 1.0f) {
                            stop_weight = 0.5f;  // Moderate evidence
                        } else if (log_likelihood_ratio > 0.5f) {
                            stop_weight = 0.7f;  // Weak evidence
                        }
                        // LLR <= 0.5: stop codon is equally or more likely, keep full penalty
                    }
                }
            }

            damage_weighted_stops += stop_weight;
        }
    }

    // Return multiplier based on weighted stop count
    // Use standard penalty scale: each weighted stop has strong impact
    // 0 stops → 1.0, 1 full stop → 0.2, 2+ full stops → 0.05
    if (damage_weighted_stops < 0.01f) return 1.0f;
    if (damage_weighted_stops < 1.0f) return 1.0f - (0.8f * damage_weighted_stops);  // Linear decay
    return 0.2f / damage_weighted_stops;  // Strong penalty for multiple stops
}

std::vector<FrameScore> FrameSelector::score_all_frames_full(
    const std::string& seq,
    const DamageProfile* damage) {

    std::vector<FrameScore> scores;
    scores.reserve(6);

    if (seq.length() < 3) {
        return scores;
    }

    // Compute reverse complement ONCE
    const std::string& rc_seq = reverse_complement_fast(seq);

    std::string protein;

    // Pre-compute unorthodox and dicodon scores once per strand (3x speedup)
    PrecomputedFrameScores fwd_precomputed = precompute_frame_scores(seq);
    PrecomputedFrameScores rev_precomputed = precompute_frame_scores(rc_seq);

    // Phase 2: Pre-compute stop counts for cross-frame analysis
    PrecomputedStopCounts fwd_stops = precompute_stop_counts(seq);
    PrecomputedStopCounts rev_stops = precompute_stop_counts(rc_seq);

    // Score all 6 frames
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        FrameScore fwd;
        fwd.frame = frame;
        fwd.forward = true;
        fwd.total_score = score_frame_fast(seq, frame, protein, &fwd_precomputed, &fwd_stops);
        fwd.protein = protein;
        if (fwd.total_score > -100.0f) {
            fwd.codon_score = 0.0f;  // Deprecated - codon scoring removed
            fwd.stop_codon_penalty = calculate_stop_penalty(protein);
            fwd.aa_composition_score = calculate_aa_composition_score(protein);
            fwd.damage_consistency = score_damage_frame_consistency(seq, frame);
        }
        scores.push_back(std::move(fwd));

        // Reverse frame
        FrameScore rev;
        rev.frame = frame;
        rev.forward = false;
        rev.total_score = score_frame_fast(rc_seq, frame, protein, &rev_precomputed, &rev_stops);
        rev.protein = protein;
        if (rev.total_score > -100.0f) {
            rev.codon_score = 0.0f;  // Deprecated - codon scoring removed
            rev.stop_codon_penalty = calculate_stop_penalty(protein);
            rev.aa_composition_score = calculate_aa_composition_score(protein);
            rev.damage_consistency = score_damage_frame_consistency(rc_seq, frame);
        }
        scores.push_back(std::move(rev));
    }

    // Sort by score descending
    std::sort(scores.begin(), scores.end(),
              [](const FrameScore& a, const FrameScore& b) {
                  return a.total_score > b.total_score;
              });

    return scores;
}

// ============================================================================
// Bayesian Frame Selection
// Uses hexamer log-likelihood ratios with posterior probability model
// ============================================================================

std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand_bayesian(
    const std::string& seq,
    Domain domain,
    const SampleDamageProfile* sample_profile) {

    // Create Bayesian selector for specified domain
    BayesianFrameSelector bayesian(domain);

    // Set damage parameters from sample profile if available
    if (sample_profile && sample_profile->is_valid()) {
        bayesian.set_damage_params(
            sample_profile->max_damage_5prime,
            sample_profile->max_damage_3prime,
            sample_profile->lambda_5prime,
            sample_profile->lambda_3prime
        );
    }

    // Get Bayesian frame selection result
    BayesianFrameOutput result = bayesian.select_frames(seq);

    // Create invalid FrameScore (total_score < 0 indicates invalid)
    FrameScore invalid;
    invalid.total_score = -1000.0f;
    invalid.protein = "";

    // Convert BayesianFrameOutput to pair<FrameScore, FrameScore>
    if (result.selected_frames.empty()) {
        // No frames selected (noncoding)
        return {invalid, invalid};
    }

    // Convert first (best) frame
    const auto& best = result.selected_frames[0];
    FrameScore first;
    first.frame = best.frame;
    first.forward = !best.is_reverse;
    first.total_score = best.posterior;  // Use posterior as score (0-1)
    first.protein = best.protein;
    first.stop_codon_penalty = 0.0f;  // Not computed in Bayesian mode
    first.codon_score = best.log_likelihood;
    first.aa_composition_score = result.confidence;

    if (result.selected_frames.size() >= 2) {
        // Ambiguous: return both frames
        const auto& second = result.selected_frames[1];
        FrameScore second_frame;
        second_frame.frame = second.frame;
        second_frame.forward = !second.is_reverse;
        second_frame.total_score = second.posterior;
        second_frame.protein = second.protein;
        second_frame.codon_score = second.log_likelihood;
        return {first, second_frame};
    }

    // Confident: only first frame is valid
    return {first, invalid};
}

} // namespace agp
