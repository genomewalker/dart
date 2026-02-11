// Per-read damage probability estimation for ancient DNA
//
// Uses damage-aware hexamer model:
// - Computes P(hexamer | damaged) by summing over source hexamers that could
//   produce observed hexamer via C→T (5') or G→A (3') damage
// - Combines hexamer LLR with base composition, stop codon, and AA signals
// - Prior from sample d_max and expected damage events

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include "agp/damage_hexamer.hpp"
#include "agp/types.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace agp {

namespace {

inline float safe_log(float x) {
    return std::log(std::max(x, 1e-12f));
}

// Observation for likelihood calculation
struct DamageObs {
    float d;     // damage rate at this position
    char kind;   // 'T', 'C', 'A', or 'G'
};

// Log-likelihood for scaling factor s
float compute_ll(float s, const std::vector<DamageObs>& obs, float pT, float pC, float pA, float pG) {
    s = std::clamp(s, 0.0f, 1.0f);
    float ll = 0.0f;
    for (const auto& o : obs) {
        switch (o.kind) {
            case 'T': ll += safe_log(pT + pC * s * o.d); break;
            case 'C': ll += safe_log(pC * (1.0f - s * o.d)); break;
            case 'A': ll += safe_log(pA + pG * s * o.d); break;
            case 'G': ll += safe_log(pG * (1.0f - s * o.d)); break;
        }
    }
    return ll;
}

// Gradient of log-likelihood
float compute_dll(float s, const std::vector<DamageObs>& obs, float pT, float pC, float pA, float pG) {
    s = std::clamp(s, 0.0f, 1.0f);
    float g = 0.0f;
    for (const auto& o : obs) {
        switch (o.kind) {
            case 'T': g += (pC * o.d) / std::max(pT + pC * s * o.d, 1e-12f); break;
            case 'C': g += (-o.d) / std::max(1.0f - s * o.d, 1e-12f); break;
            case 'A': g += (pG * o.d) / std::max(pA + pG * s * o.d, 1e-12f); break;
            case 'G': g += (-o.d) / std::max(1.0f - s * o.d, 1e-12f); break;
        }
    }
    return g;
}

// Hessian of log-likelihood (negative for concave)
float compute_ddll(float s, const std::vector<DamageObs>& obs, float pT, float pC, float pA, float pG) {
    s = std::clamp(s, 0.0f, 1.0f);
    float h = 0.0f;
    for (const auto& o : obs) {
        switch (o.kind) {
            case 'T': {
                float denom = std::max(pT + pC * s * o.d, 1e-12f);
                h -= (pC * o.d) * (pC * o.d) / (denom * denom);
                break;
            }
            case 'C': {
                float denom = std::max(1.0f - s * o.d, 1e-12f);
                h -= (o.d * o.d) / (denom * denom);
                break;
            }
            case 'A': {
                float denom = std::max(pA + pG * s * o.d, 1e-12f);
                h -= (pG * o.d) * (pG * o.d) / (denom * denom);
                break;
            }
            case 'G': {
                float denom = std::max(1.0f - s * o.d, 1e-12f);
                h -= (o.d * o.d) / (denom * denom);
                break;
            }
        }
    }
    return h;
}

} // anonymous namespace

// Simple heuristic damage signal (fallback when no sample profile)
float FrameSelector::estimate_damage_signal(const std::string& seq) {
    if (seq.length() < 10) return 0.5f;

    const size_t len = seq.length();
    float log_odds = 0.0f;

    // Terminal base patterns
    if (fast_upper(seq[0]) == 'T') log_odds += 0.8f;
    else if (fast_upper(seq[0]) == 'C') log_odds -= 0.3f;

    if (fast_upper(seq[len - 1]) == 'A') log_odds += 0.7f;
    else if (fast_upper(seq[len - 1]) == 'G') log_odds -= 0.25f;

    // Position 1
    if (len > 1) {
        if (fast_upper(seq[1]) == 'T') log_odds += 0.35f;
        if (fast_upper(seq[len - 2]) == 'A') log_odds += 0.3f;
    }

    // Position 2
    if (len > 2) {
        if (fast_upper(seq[2]) == 'T') log_odds += 0.15f;
        if (fast_upper(seq[len - 3]) == 'A') log_odds += 0.12f;
    }

    // Joint pattern bonus
    if (fast_upper(seq[0]) == 'T' && fast_upper(seq[len - 1]) == 'A') {
        log_odds += 0.4f;
    }

    // Convert to probability
    float prob = 1.0f / (1.0f + std::exp(-log_odds));
    return std::clamp(prob, 0.05f, 0.95f);
}

// GC-conditional damage percentage estimation
// Returns damage_pct in [0, 100] using GC-bin specific damage rate
float FrameSelector::compute_damage_percentage(
    const std::string& seq,
    const SampleDamageProfile& sample_profile) {

    constexpr int MIN_LEN = 15;
    constexpr int WINDOW = 15;

    if (seq.length() < MIN_LEN || !sample_profile.is_valid()) return 0.0f;

    // Tri-state damage validation: CONTRADICTED = hard suppress damage percentage
    if (get_damage_validation_state(sample_profile) == DamageValidationState::CONTRADICTED) {
        return 0.0f;
    }

    // Get GC-conditional parameters for this read
    auto gc_params = sample_profile.get_gc_params(seq);

    // If this GC bin has no detected damage, return 0
    if (gc_params.bin_valid && gc_params.delta_s < 0.005f) {
        return 0.0f;
    }

    const size_t len = seq.length();
    const int W = std::min(WINDOW, static_cast<int>(len));

    // Use GC-bin specific baseline
    float baseline_tc = gc_params.bin_valid ? gc_params.baseline_tc :
        static_cast<float>(sample_profile.baseline_t_freq /
            (sample_profile.baseline_t_freq + sample_profile.baseline_c_freq + 0.001));
    float baseline_ag = static_cast<float>(sample_profile.baseline_a_freq /
        (sample_profile.baseline_a_freq + sample_profile.baseline_g_freq + 0.001));

    baseline_tc = std::clamp(baseline_tc, 0.01f, 0.99f);
    baseline_ag = std::clamp(baseline_ag, 0.01f, 0.99f);

    const float pT = baseline_tc;
    const float pC = 1.0f - pT;
    const float pA = baseline_ag;
    const float pG = 1.0f - pA;

    // Collect informative observations using GC-bin specific damage rates
    std::vector<DamageObs> obs;
    obs.reserve(2 * W);

    float lambda = gc_params.lambda;  // lambda_5prime from GC params
    float lambda_3p = sample_profile.lambda_3prime > 0 ? sample_profile.lambda_3prime : lambda;
    float delta_s = gc_params.bin_valid ? gc_params.delta_s : sample_profile.d_max_combined;

    for (int i = 0; i < W; ++i) {
        // 5' end (C→T damage) - use GC-bin specific delta with decay
        float d5 = delta_s * std::exp(-lambda * i);
        d5 = std::min(d5, 0.999f);
        char b5 = fast_upper(seq[i]);
        if (b5 == 'T' || b5 == 'C') {
            obs.push_back({d5, b5});
        }

        // 3' end (G→A damage) - use separate lambda_3p for potentially asymmetric decay
        float d3 = delta_s * std::exp(-lambda_3p * i);
        d3 = std::min(d3, 0.999f);
        char b3 = fast_upper(seq[len - 1 - i]);
        if (b3 == 'A' || b3 == 'G') {
            obs.push_back({d3, b3});
        }
    }

    // No informative observations
    if (obs.empty()) return 0.0f;

    // Find s_hat via grid search + Newton refinement
    float best_s = 0.0f;
    float best_ll = compute_ll(0.0f, obs, pT, pC, pA, pG);

    // Coarse grid
    for (int k = 1; k <= 10; ++k) {
        float s = 0.1f * k;
        float ll = compute_ll(s, obs, pT, pC, pA, pG);
        if (ll > best_ll) {
            best_ll = ll;
            best_s = s;
        }
    }

    // Newton refinement
    float s = best_s;
    for (int it = 0; it < 8; ++it) {
        float g = compute_dll(s, obs, pT, pC, pA, pG);
        float h = compute_ddll(s, obs, pT, pC, pA, pG);

        if (std::fabs(g) < 1e-4f) break;
        if (h >= -1e-8f) break;

        float step = -g / h;
        float s_new = std::clamp(s + step, 0.0f, 1.0f);

        // Backtracking
        float ll_s = compute_ll(s, obs, pT, pC, pA, pG);
        float ll_new = compute_ll(s_new, obs, pT, pC, pA, pG);
        int bt = 0;
        while (ll_new < ll_s && bt++ < 10) {
            step *= 0.5f;
            s_new = std::clamp(s + step, 0.0f, 1.0f);
            ll_new = compute_ll(s_new, obs, pT, pC, pA, pG);
        }

        if (ll_new <= ll_s + 1e-8f) break;
        s = s_new;
    }

    // Convert s_hat to damage percentage using GC-bin specific delta
    // s_hat scales the bin's delta, so damage_pct = s_hat * delta_s * 100
    float damage_pct = 100.0f * std::clamp(s * delta_s, 0.0f, 1.0f);

    return damage_pct;
}

} // namespace agp
