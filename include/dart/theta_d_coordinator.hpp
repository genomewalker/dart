#pragma once
// ThetaDSuffStats: the gamma-weighted d_max update of damage-annotate --refine-damage.
//
// The site scorer (bayesian_damage_score.hpp) models a damage-consistent amino-acid hit at a
// susceptible position u of read r as Bernoulli with rate
//   q(u) = q_modern + AA_SCALE * d_max * exp(-lambda * dist_u),   AA_SCALE = 0.30,
// so E[k_r] = q_modern * m_r + AA_SCALE * d_max * decay_r with m_r the susceptible opportunities
// and decay_r = sum_u exp(-lambda * dist_u). Matching the gamma-weighted expected and observed hits
// over the same reads gives the moment estimator
//   d_max_hat = (sum gamma*k - q_modern * sum gamma*m) / (AA_SCALE * sum gamma*decay),
// clipped to [0, 1]. It is not the exact Bernoulli MLE (that has no closed form here) but it is
// on the scorer's scale. Every read with opportunities contributes to all three sums; reads
// without a mismatch contribute k = 0.

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace dart {

struct ThetaDSuffStats {
    static constexpr double AA_SCALE = 0.30;   // must equal the site scorer's AA_SCALE
    double sum_gamma_k      = 0.0;  // sum_r gamma_r * k_r
    double sum_gamma_m      = 0.0;  // sum_r gamma_r * m_r
    double sum_gamma_decay  = 0.0;  // sum_r gamma_r * decay_r
    uint64_t n_reads        = 0;    // reads with gamma > 0 and opportunities

    void accumulate(double gamma, uint32_t k, uint32_t m, float decay) {
        if (gamma <= 0.0 || m == 0 || decay <= 0.0f) return;
        sum_gamma_k     += gamma * static_cast<double>(k);
        sum_gamma_m     += gamma * static_cast<double>(m);
        sum_gamma_decay += gamma * static_cast<double>(decay);
        ++n_reads;
    }

    // d_max_fallback when there is too little evidence to update.
    float refit(float d_max_fallback, double q_modern) const {
        if (sum_gamma_decay < 1e-10 || n_reads == 0) return d_max_fallback;
        const double d = (sum_gamma_k - q_modern * sum_gamma_m) / (AA_SCALE * sum_gamma_decay);
        return static_cast<float>(std::clamp(d, 0.0, 1.0));
    }
};

} // namespace dart
