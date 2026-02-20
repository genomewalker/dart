#pragma once

/**
 * UnifiedDamageContext - Unified interface for damage-aware gene prediction
 *
 * This struct provides a consistent view of damage information for all components:
 * - Frame selection (stop penalty weighting)
 * - Damage correction (threshold calculation)
 * - Confidence estimation
 *
 * Key design principles:
 * 1. Non-owning view over SampleDamageProfile (single source of truth)
 * 2. Per-read evidence (damage_signal, confidence) added at runtime
 * 3. Strand-aware: for_reverse_complement() correctly swaps 5'/3' mutation meaning
 * 4. Gates all damage-aware logic via enabled() check
 */

#include "dart/frame_selector.hpp"
#include <algorithm>
#include <cmath>
#include <array>

namespace dart {

// Per-read damage evidence computed from sequence
struct ReadDamageEvidence {
    float damage_signal = 0.5f;     // P(terminal damage | read), calibrated to sample d_max
    float damage_confidence = 1.0f; // Reliability weight (0-1), reduced for short/ambiguous reads

    static ReadDamageEvidence neutral() {
        return {0.5f, 0.0f};  // Neutral signal, no confidence
    }
};

// Unified damage context combining sample-level and per-read information
struct UnifiedDamageContext {
    const SampleDamageProfile* sample = nullptr;  // Non-owning pointer to sample profile
    ReadDamageEvidence read;                       // Per-read evidence
    bool is_reverse_complement = false;            // Whether we're analyzing RC sequence

    // Frame posterior information (optional, for correction threshold coupling)
    float frame_entropy = 0.0f;     // Entropy of frame posterior (high = uncertain)
    int selected_frame = -1;        // Selected frame (0-5, or -1 if not yet selected)

    // =========================================================================
    // Core accessors
    // =========================================================================

    /**
     * Check if damage-aware logic should be applied.
     * Returns false for artifacts, unvalidated samples, or missing data.
     */
    bool enabled() const {
        if (!sample) return false;
        if (sample->damage_artifact) return false;
        if (!sample->damage_validated && sample->d_max_combined < 0.05f) return false;
        if (read.damage_confidence < 0.01f) return false;
        return true;
    }

    /**
     * Get sample's validated d_max (0-1 scale).
     * Returns 0 if sample not available or artifact.
     */
    float d_max() const {
        if (!sample) return 0.0f;
        if (sample->damage_artifact) return 0.0f;
        return static_cast<float>(sample->d_max_combined);
    }

    /**
     * Get the effective damage gate: product of per-read signal and confidence.
     * Used to modulate how much we trust damage-based adjustments.
     */
    float damage_gate() const {
        return std::clamp(read.damage_signal, 0.0f, 1.0f) *
               std::clamp(read.damage_confidence, 0.0f, 1.0f);
    }

    // =========================================================================
    // Position-specific damage rates
    // =========================================================================

    /**
     * Get C→T damage rate at absolute position (in current sequence orientation).
     * For forward strand: uses 5' rates directly
     * For reverse complement: uses 3' rates (original 3' G→A becomes RC 5' C→T)
     */
    float ct_rate_at(size_t abs_pos, size_t seq_len) const {
        if (!sample || abs_pos >= 15) return 0.0f;

        if (!is_reverse_complement) {
            // Forward: 5' C→T damage
            float raw_rate = sample->damage_rate_5prime[abs_pos];
            // Scale to Channel B validated d_max
            float scale = (sample->max_damage_5prime > 1e-6f)
                ? (static_cast<float>(sample->d_max_5prime) / sample->max_damage_5prime)
                : 0.0f;
            return std::clamp(raw_rate * scale, 0.0f, 1.0f);
        } else {
            // Reverse complement: original 3' G→A becomes C→T on RC 5' end
            float raw_rate = sample->damage_rate_3prime[abs_pos];
            float scale = (sample->max_damage_3prime > 1e-6f)
                ? (static_cast<float>(sample->d_max_3prime) / sample->max_damage_3prime)
                : 0.0f;
            return std::clamp(raw_rate * scale, 0.0f, 1.0f);
        }
    }

    /**
     * Get G→A damage rate at absolute position (in current sequence orientation).
     * For forward strand: uses 3' rates (distance from end)
     * For reverse complement: uses 5' rates (original 5' C→T becomes RC 3' G→A)
     */
    float ga_rate_at(size_t abs_pos, size_t seq_len) const {
        if (!sample || seq_len == 0) return 0.0f;

        size_t dist_from_end = seq_len - 1 - abs_pos;
        if (dist_from_end >= 15) return 0.0f;

        if (!is_reverse_complement) {
            // Forward: 3' G→A damage
            float raw_rate = sample->damage_rate_3prime[dist_from_end];
            float scale = (sample->max_damage_3prime > 1e-6f)
                ? (static_cast<float>(sample->d_max_3prime) / sample->max_damage_3prime)
                : 0.0f;
            return std::clamp(raw_rate * scale, 0.0f, 1.0f);
        } else {
            // Reverse complement: original 5' C→T becomes G→A on RC 3' end
            float raw_rate = sample->damage_rate_5prime[dist_from_end];
            float scale = (sample->max_damage_5prime > 1e-6f)
                ? (static_cast<float>(sample->d_max_5prime) / sample->max_damage_5prime)
                : 0.0f;
            return std::clamp(raw_rate * scale, 0.0f, 1.0f);
        }
    }

    // =========================================================================
    // Context transformations
    // =========================================================================

    /**
     * Create a context for analyzing the reverse complement sequence.
     * Correctly swaps 5'/3' mutation meaning.
     */
    UnifiedDamageContext for_reverse_complement() const {
        UnifiedDamageContext ctx = *this;
        ctx.is_reverse_complement = !is_reverse_complement;
        return ctx;
    }

    /**
     * Create a context with updated per-read evidence.
     */
    UnifiedDamageContext with_read_evidence(const ReadDamageEvidence& evidence) const {
        UnifiedDamageContext ctx = *this;
        ctx.read = evidence;
        return ctx;
    }

    /**
     * Create a context with frame posterior information.
     */
    UnifiedDamageContext with_frame_info(float entropy, int frame) const {
        UnifiedDamageContext ctx = *this;
        ctx.frame_entropy = entropy;
        ctx.selected_frame = frame;
        return ctx;
    }

    // =========================================================================
    // Threshold calculations
    // =========================================================================

    /**
     * Calculate the maximum discount to apply for damage-induced stops.
     * Higher d_max + higher damage_signal → more aggressive discounting.
     * Range: 0.0 (no discount) to ~0.7 (heavy discount)
     */
    float max_stop_discount() const {
        if (!enabled()) return 0.0f;

        // Base discount scaled by sample damage level
        float base_discount = 0.7f * std::min(1.0f, d_max() * 3.0f);

        // Gate by per-read evidence
        return base_discount * damage_gate();
    }

    /**
     * Calculate correction threshold based on damage level and frame confidence.
     * Higher entropy (uncertain frame) → higher threshold (more conservative)
     * Higher d_max → lower threshold (more aggressive correction)
     * Range: 0.5 to 0.95
     */
    float correction_threshold() const {
        if (!enabled()) return 0.99f;  // Essentially never correct
        if (sample->damage_artifact) return 0.99f;

        // Base threshold: lower for high-damage samples
        float base = 0.7f - 0.2f * d_max();

        // Entropy adjustment: raise threshold when frame is uncertain
        // frame_entropy of ~1.8 (uniform over 6) → add 0.15
        // frame_entropy of ~0.5 (fairly confident) → add 0.05
        float entropy_adj = 0.08f * frame_entropy;

        return std::clamp(base + entropy_adj, 0.5f, 0.95f);
    }

    // =========================================================================
    // Factory functions
    // =========================================================================

    /**
     * Create context from sample profile (without per-read evidence).
     */
    static UnifiedDamageContext from_sample(const SampleDamageProfile& profile) {
        UnifiedDamageContext ctx;
        ctx.sample = &profile;
        ctx.read = ReadDamageEvidence::neutral();
        return ctx;
    }

    /**
     * Create context from sample profile with per-read evidence.
     */
    static UnifiedDamageContext from_sample_and_read(
        const SampleDamageProfile& profile,
        float damage_signal,
        float confidence = 1.0f
    ) {
        UnifiedDamageContext ctx;
        ctx.sample = &profile;
        ctx.read.damage_signal = damage_signal;
        ctx.read.damage_confidence = confidence;
        return ctx;
    }

    /**
     * Create a disabled/neutral context (for undamaged samples or testing).
     */
    static UnifiedDamageContext disabled() {
        return UnifiedDamageContext{};
    }
};

// =========================================================================
// Utility functions
// =========================================================================

/**
 * Compute frame posterior entropy from log-scores.
 * Higher entropy = more uncertain frame selection.
 */
inline float compute_frame_entropy(const std::array<float, 6>& log_scores) {
    // Convert to probabilities via softmax
    float max_score = *std::max_element(log_scores.begin(), log_scores.end());

    std::array<float, 6> probs;
    float sum = 0.0f;
    for (int i = 0; i < 6; ++i) {
        probs[i] = std::exp(log_scores[i] - max_score);
        sum += probs[i];
    }

    if (sum < 1e-10f) return 0.0f;

    // Compute entropy: -sum(p * log(p))
    float entropy = 0.0f;
    for (int i = 0; i < 6; ++i) {
        float p = probs[i] / sum;
        if (p > 1e-10f) {
            entropy -= p * std::log(p);
        }
    }

    return entropy;
}

/**
 * Compute frame posterior from log-scores.
 * Returns normalized probabilities for all 6 frames.
 */
inline std::array<float, 6> compute_frame_posterior(const std::array<float, 6>& log_scores) {
    float max_score = *std::max_element(log_scores.begin(), log_scores.end());

    std::array<float, 6> probs;
    float sum = 0.0f;
    for (int i = 0; i < 6; ++i) {
        probs[i] = std::exp(log_scores[i] - max_score);
        sum += probs[i];
    }

    if (sum > 1e-10f) {
        for (int i = 0; i < 6; ++i) {
            probs[i] /= sum;
        }
    }

    return probs;
}

}  // namespace dart
