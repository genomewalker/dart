#pragma once
// ThetaDCoordinator: compute a gamma-weighted MLE for d_max.
// Used by the Phase 2 unified-objective path in damage-annotate (--refine-damage).
//
// d_max MLE (lambda fixed):
//
//   d_max_hat = sum_r(gamma_r * k_r) / sum_r(gamma_r * sum_exp_decay_r)
//
// where k_r     = aa_k_hits (damage-consistent AA hits for read r)
//       sum_exp_decay_r = aa_sum_exp_decay (Σ_pos exp(-λ·dist) over susceptible positions)
//       gamma_r = gamma_ancient_r when the damage EM is active, else gamma_r (total)
//
// This is exact MLE under the per-position Bernoulli model
//   P(hit at pos | ancient) = d_max * exp(-λ·dist).
//
// Important: the numerator (sum_gamma_k) and denominator (sum_gamma_decay) are
// accumulated from different read subsets in cmd_damage_annotate.cpp:
//   - Numerator:   only reads with total_mismatches > 0 (in summaries[])
//   - Denominator: ALL reads with best-hit alignments (per_read_aa_decay[])
// This is necessary because zero-mismatch reads contribute 0 to the numerator
// but must contribute their decay to the denominator to avoid upward bias.

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace dart {

struct ThetaDSuffStats {
    double sum_gamma_k      = 0.0;  // Σ_r gamma_r * aa_k_hits_r  (numerator)
    double sum_gamma_decay  = 0.0;  // Σ_r gamma_r * aa_sum_exp_decay_r  (denominator)
    uint64_t n_reads        = 0;    // Reads contributing to numerator

    void clear() {
        sum_gamma_k     = 0.0;
        sum_gamma_decay = 0.0;
        n_reads         = 0;
    }

    // Accumulate a single read into both numerator and denominator.
    // Use this when a read contributes to both (e.g., in a single-pass loop).
    void accumulate(double gamma, uint32_t aa_k_hits, float aa_sum_exp_decay) {
        if (gamma <= 0.0 || aa_sum_exp_decay <= 0.0f) return;
        sum_gamma_k     += gamma * static_cast<double>(aa_k_hits);
        sum_gamma_decay += gamma * static_cast<double>(aa_sum_exp_decay);
        ++n_reads;
    }

    // MLE for d_max. Returns d_max_fallback when insufficient data.
    float refit(float d_max_fallback) const {
        if (sum_gamma_decay < 1e-10 || n_reads == 0) return d_max_fallback;
        float d_new = static_cast<float>(sum_gamma_k / sum_gamma_decay);
        return std::clamp(d_new, 0.0f, 1.0f);
    }
};

} // namespace dart
