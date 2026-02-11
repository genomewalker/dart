#pragma once

#include <array>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>

namespace agp {

/**
 * Durbin-Style Mixture Model for Metagenome Damage Quantification
 *
 * Generative story:
 * 1. Each C-site comes from latent organism class k ~ Categorical(π)
 * 2. Class k has GC distribution P(gc_bin | k) and damage rate δ_max,k
 * 3. Terminal C→T follows: P(T | C, position p, class k) = δ_max,k · exp(-λp)
 * 4. Shared artifact affects both damage and control channels
 *
 * This elegantly resolves the "reference-free vs metaDMG" gap:
 * - d_population = E[δ] over all C-sites (what we observe)
 * - d_reference = E[δ | GC > 50%] (what metaDMG approximates)
 */

static constexpr int N_GC_BINS = 10;
static constexpr int MAX_K = 4;
static constexpr int N_POSITIONS = 15;

struct MixtureDamageResult {
    int K = 0;                              // Number of classes
    std::array<float, MAX_K> pi = {};       // Mixing proportions (by C-sites)
    std::array<float, MAX_K> delta_max = {};// Damage rates per class
    std::array<float, MAX_K> mu_gc = {};    // Mean GC per class (for interpretation)

    float lambda = 0.0f;                    // Shared decay constant
    float a_max = 0.0f;                     // Shared artifact amplitude

    // Summary statistics
    float d_population = 0.0f;   // E[δ] = Σ_k π_k · δ_max,k
    float d_ancient = 0.0f;      // E[δ | δ > τ] (ancient tail, τ=5%)
    float d_reference = 0.0f;    // E[δ | GC > 50%] (metaDMG proxy)
    float pi_ancient = 0.0f;     // P(class with δ > τ)

    // Model fit
    float log_likelihood = 0.0f;
    float bic = 0.0f;
    int n_iterations = 0;
    bool converged = false;
};

struct SuperRead {
    int gc_bin = 0;  // 0-9 for 0-10%, 10-20%, ..., 90-100%

    // Channel A: T/(T+C) at 5' positions 0-14
    std::array<double, N_POSITIONS> k_tc = {};
    std::array<double, N_POSITIONS> n_tc = {};

    // Control: A/(A+G) at 5' positions 0-14
    std::array<double, N_POSITIONS> k_ag = {};
    std::array<double, N_POSITIONS> n_ag = {};

    // Channel B: stop conversion
    std::array<double, N_POSITIONS> k_stop = {};
    std::array<double, N_POSITIONS> n_stop = {};

    // Interior baselines
    double k_tc_int = 0, n_tc_int = 0;
    double k_ag_int = 0, n_ag_int = 0;
    double k_stop_int = 0, n_stop_int = 0;

    double c_sites = 0;  // Total C-sites (weight)
};

class MixtureDamageModel {
public:
    static constexpr float ANCIENT_THRESHOLD = 0.05f;  // τ for d_ancient
    static constexpr int REFERENCE_GC_MIN = 5;         // GC >= 50% for d_reference
    static constexpr int MAX_ITER = 100;
    static constexpr float CONVERGENCE_TOL = 1e-6f;
    static constexpr int N_RESTARTS = 5;
    static constexpr double PI_MIN = 1e-6;
    static constexpr float BASELINE_SHRINKAGE_ALPHA = 1000.0f;

    /**
     * Fit mixture model with automatic K selection via BIC
     */
    static MixtureDamageResult fit(const std::array<SuperRead, N_GC_BINS>& super_reads);

private:
    /**
     * Fit mixture model with fixed K
     */
    static MixtureDamageResult fit_k(
        const std::array<SuperRead, N_GC_BINS>& super_reads,
        int K,
        int restart_seed);

    /**
     * Compute binomial log-likelihood for channels
     */
    static double channel_log_likelihood(
        const SuperRead& sr,
        float delta_max,
        float lambda,
        float a_max,
        float b_tc,
        float b_ag,
        float b_stop);

    /**
     * Log-sum-exp for numerical stability
     */
    static double logsumexp(const double* values, int n) {
        double max_val = values[0];
        for (int i = 1; i < n; ++i) {
            max_val = std::max(max_val, values[i]);
        }
        if (std::isinf(max_val)) return max_val;

        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += std::exp(values[i] - max_val);
        }
        return max_val + std::log(sum);
    }

    static double binom_ll(double k, double n, double p) {
        if (n < 1.0) return 0.0;
        p = std::clamp(p, 1e-10, 1.0 - 1e-10);
        return k * std::log(p) + (n - k) * std::log(1.0 - p);
    }
};

// ============================================================================
// Implementation
// ============================================================================

inline double MixtureDamageModel::channel_log_likelihood(
    const SuperRead& sr,
    float delta_max,
    float lambda,
    float a_max,
    float b_tc,
    float b_ag,
    float b_stop)
{
    double ll = 0.0;

    for (int p = 0; p < N_POSITIONS; ++p) {
        float decay = std::exp(-lambda * p);
        float delta_p = delta_max * decay;
        float a_p = a_max * decay;

        // Channel A: π_TC(p) = b_tc + a(p) + (1 - b_tc - a(p)) · δ(p)
        float base_tc = std::clamp(b_tc + a_p, 0.01f, 0.99f);
        float pi_tc = base_tc + (1.0f - base_tc) * delta_p;
        ll += binom_ll(sr.k_tc[p], sr.n_tc[p], pi_tc);

        // Control: π_AG(p) = b_ag + a(p)
        float pi_ag = std::clamp(b_ag + a_p, 0.01f, 0.99f);
        ll += binom_ll(sr.k_ag[p], sr.n_ag[p], pi_ag);

        // Channel B: π_stop(p) = b_stop + (1 - b_stop) · δ(p)
        float pi_stop = std::clamp(b_stop + (1.0f - b_stop) * delta_p, 0.001f, 0.999f);
        ll += binom_ll(sr.k_stop[p], sr.n_stop[p], pi_stop);
    }

    return ll;
}

inline MixtureDamageResult MixtureDamageModel::fit_k(
    const std::array<SuperRead, N_GC_BINS>& super_reads,
    int K,
    int restart_seed)
{
    MixtureDamageResult result;
    result.K = K;

    // Compute total C-sites and global baselines
    double total_c_sites = 0.0;
    double global_k_tc_int = 0, global_n_tc_int = 0;
    double global_k_ag_int = 0, global_n_ag_int = 0;
    double global_k_stop_int = 0, global_n_stop_int = 0;

    for (const auto& sr : super_reads) {
        total_c_sites += sr.c_sites;
        global_k_tc_int += sr.k_tc_int;
        global_n_tc_int += sr.n_tc_int;
        global_k_ag_int += sr.k_ag_int;
        global_n_ag_int += sr.n_ag_int;
        global_k_stop_int += sr.k_stop_int;
        global_n_stop_int += sr.n_stop_int;
    }

    if (total_c_sites < 1000.0) {
        return result;  // Not enough data
    }

    float b_tc_global = global_n_tc_int > 0 ? static_cast<float>(global_k_tc_int / global_n_tc_int) : 0.5f;
    float b_ag_global = global_n_ag_int > 0 ? static_cast<float>(global_k_ag_int / global_n_ag_int) : 0.5f;
    float b_stop_global = global_n_stop_int > 0 ? static_cast<float>(global_k_stop_int / global_n_stop_int) : 0.05f;

    // Initialize parameters
    // Mixing proportions: uniform
    std::array<double, MAX_K> pi;
    for (int k = 0; k < K; ++k) {
        pi[k] = 1.0 / K;
    }

    // GC categorical: P(gc_bin | class) with Dirichlet smoothing
    std::array<std::array<double, N_GC_BINS>, MAX_K> p_gc;
    for (int k = 0; k < K; ++k) {
        for (int b = 0; b < N_GC_BINS; ++b) {
            // Initialize with slight preference for different GC ranges
            int center = (k * N_GC_BINS) / K + N_GC_BINS / (2 * K);
            double dist = std::abs(b - center);
            p_gc[k][b] = (1.0 + std::exp(-dist)) / N_GC_BINS;
        }
        // Normalize
        double sum = 0;
        for (int b = 0; b < N_GC_BINS; ++b) sum += p_gc[k][b];
        for (int b = 0; b < N_GC_BINS; ++b) p_gc[k][b] /= sum;
    }

    // Initialize delta_max from quantiles
    std::vector<std::pair<float, double>> bin_damage;  // (d_max, weight)
    for (int b = 0; b < N_GC_BINS; ++b) {
        const auto& sr = super_reads[b];
        if (sr.c_sites > 0 && sr.n_tc[0] > 0) {
            float tc_ratio = static_cast<float>(sr.k_tc[0] / sr.n_tc[0]);
            float d_est = std::max(0.0f, (tc_ratio - b_tc_global) / (1.0f - b_tc_global));
            bin_damage.push_back({d_est, sr.c_sites});
        }
    }

    std::array<float, MAX_K> delta_max;
    if (!bin_damage.empty()) {
        std::sort(bin_damage.begin(), bin_damage.end());
        for (int k = 0; k < K; ++k) {
            int idx = (k * static_cast<int>(bin_damage.size())) / K;
            delta_max[k] = std::clamp(bin_damage[idx].first + 0.01f * (restart_seed + k), 0.0f, 0.6f);
        }
    } else {
        for (int k = 0; k < K; ++k) {
            delta_max[k] = 0.1f * (k + 1);
        }
    }

    // Shared parameters
    float lambda = 0.2f;
    float a_max = 0.0f;

    // Responsibilities: w[s][k]
    std::array<std::array<double, MAX_K>, N_GC_BINS> w;

    // EM iterations
    double prev_ll = -std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        // E-step: compute responsibilities
        double obs_ll = 0.0;

        for (int s = 0; s < N_GC_BINS; ++s) {
            const auto& sr = super_reads[s];
            if (sr.c_sites < 1.0) {
                for (int k = 0; k < K; ++k) w[s][k] = 0.0;
                continue;
            }

            // Compute baselines with shrinkage
            float b_tc = static_cast<float>((sr.k_tc_int + BASELINE_SHRINKAGE_ALPHA * b_tc_global) /
                                           (sr.n_tc_int + BASELINE_SHRINKAGE_ALPHA));
            float b_ag = static_cast<float>((sr.k_ag_int + BASELINE_SHRINKAGE_ALPHA * b_ag_global) /
                                           (sr.n_ag_int + BASELINE_SHRINKAGE_ALPHA));
            float b_stop = static_cast<float>((sr.k_stop_int + BASELINE_SHRINKAGE_ALPHA * b_stop_global) /
                                             (sr.n_stop_int + BASELINE_SHRINKAGE_ALPHA));

            std::array<double, MAX_K> log_comp;
            for (int k = 0; k < K; ++k) {
                double ll_gc = std::log(std::max(p_gc[k][s], 1e-10));
                double ll_ch = channel_log_likelihood(sr, delta_max[k], lambda, a_max, b_tc, b_ag, b_stop);
                log_comp[k] = std::log(std::max(pi[k], PI_MIN)) + ll_gc + ll_ch;
            }

            double log_sum = logsumexp(log_comp.data(), K);
            obs_ll += log_sum;

            for (int k = 0; k < K; ++k) {
                w[s][k] = std::exp(log_comp[k] - log_sum);
            }
        }

        // Check convergence
        if (std::abs(obs_ll - prev_ll) < CONVERGENCE_TOL * std::abs(prev_ll)) {
            result.converged = true;
            result.n_iterations = iter + 1;
            result.log_likelihood = static_cast<float>(obs_ll);
            break;
        }
        prev_ll = obs_ll;
        result.log_likelihood = static_cast<float>(obs_ll);
        result.n_iterations = iter + 1;

        // M-step

        // Update pi_k (weighted by C-sites)
        std::array<double, MAX_K> W_k = {};
        for (int s = 0; s < N_GC_BINS; ++s) {
            for (int k = 0; k < K; ++k) {
                W_k[k] += w[s][k] * super_reads[s].c_sites;
            }
        }
        double W_total = 0;
        for (int k = 0; k < K; ++k) W_total += W_k[k];
        for (int k = 0; k < K; ++k) {
            pi[k] = std::max(W_k[k] / W_total, PI_MIN);
        }

        // Update p_gc[k][b] (categorical with Dirichlet smoothing, alpha=1)
        for (int k = 0; k < K; ++k) {
            double sum = 0;
            for (int b = 0; b < N_GC_BINS; ++b) {
                p_gc[k][b] = w[b][k] * super_reads[b].c_sites + 1.0;  // +1 Dirichlet smoothing
                sum += p_gc[k][b];
            }
            for (int b = 0; b < N_GC_BINS; ++b) {
                p_gc[k][b] /= sum;
            }
        }

        // Update delta_max_k via grid search
        for (int k = 0; k < K; ++k) {
            float best_delta = delta_max[k];
            double best_ll = -std::numeric_limits<double>::infinity();

            for (int d_idx = 0; d_idx <= 60; ++d_idx) {
                float d_test = 0.01f * d_idx;
                double ll = 0.0;

                for (int s = 0; s < N_GC_BINS; ++s) {
                    const auto& sr = super_reads[s];
                    if (w[s][k] < 1e-10 || sr.c_sites < 1.0) continue;

                    float b_tc = static_cast<float>((sr.k_tc_int + BASELINE_SHRINKAGE_ALPHA * b_tc_global) /
                                                   (sr.n_tc_int + BASELINE_SHRINKAGE_ALPHA));
                    float b_ag = static_cast<float>((sr.k_ag_int + BASELINE_SHRINKAGE_ALPHA * b_ag_global) /
                                                   (sr.n_ag_int + BASELINE_SHRINKAGE_ALPHA));
                    float b_stop = static_cast<float>((sr.k_stop_int + BASELINE_SHRINKAGE_ALPHA * b_stop_global) /
                                                     (sr.n_stop_int + BASELINE_SHRINKAGE_ALPHA));

                    ll += w[s][k] * channel_log_likelihood(sr, d_test, lambda, a_max, b_tc, b_ag, b_stop);
                }

                if (ll > best_ll) {
                    best_ll = ll;
                    best_delta = d_test;
                }
            }
            delta_max[k] = best_delta;
        }

        // Update shared a_max via golden section search
        auto eval_a = [&](float a_test) {
            double ll = 0.0;
            for (int s = 0; s < N_GC_BINS; ++s) {
                const auto& sr = super_reads[s];
                if (sr.c_sites < 1.0) continue;

                float b_tc = static_cast<float>((sr.k_tc_int + BASELINE_SHRINKAGE_ALPHA * b_tc_global) /
                                               (sr.n_tc_int + BASELINE_SHRINKAGE_ALPHA));
                float b_ag = static_cast<float>((sr.k_ag_int + BASELINE_SHRINKAGE_ALPHA * b_ag_global) /
                                               (sr.n_ag_int + BASELINE_SHRINKAGE_ALPHA));
                float b_stop = static_cast<float>((sr.k_stop_int + BASELINE_SHRINKAGE_ALPHA * b_stop_global) /
                                                 (sr.n_stop_int + BASELINE_SHRINKAGE_ALPHA));

                for (int k = 0; k < K; ++k) {
                    ll += w[s][k] * channel_log_likelihood(sr, delta_max[k], lambda, a_test, b_tc, b_ag, b_stop);
                }
            }
            return ll;
        };

        const float phi = 0.618033988749895f;
        float a = -0.2f, b = 0.2f;
        float c = b - phi * (b - a);
        float d = a + phi * (b - a);
        for (int i = 0; i < 20; ++i) {
            if (eval_a(c) > eval_a(d)) { b = d; d = c; c = b - phi * (b - a); }
            else { a = c; c = d; d = a + phi * (b - a); }
        }
        a_max = (a + b) / 2.0f;
    }

    // Store final parameters
    for (int k = 0; k < K; ++k) {
        result.pi[k] = static_cast<float>(pi[k]);
        result.delta_max[k] = delta_max[k];

        // Compute mean GC for interpretation
        double mu_gc = 0.0;
        for (int b = 0; b < N_GC_BINS; ++b) {
            mu_gc += p_gc[k][b] * (b * 10 + 5);  // bin center
        }
        result.mu_gc[k] = static_cast<float>(mu_gc);
    }
    result.lambda = lambda;
    result.a_max = a_max;

    // Compute summary statistics
    result.d_population = 0.0f;
    result.pi_ancient = 0.0f;
    float ancient_weighted_sum = 0.0f;

    for (int k = 0; k < K; ++k) {
        result.d_population += result.pi[k] * result.delta_max[k];
        if (result.delta_max[k] > ANCIENT_THRESHOLD) {
            result.pi_ancient += result.pi[k];
            ancient_weighted_sum += result.pi[k] * result.delta_max[k];
        }
    }

    result.d_ancient = result.pi_ancient > 0.01f ? ancient_weighted_sum / result.pi_ancient : 0.0f;

    // d_reference = E[δ | GC >= 50%]
    float ref_num = 0.0f, ref_den = 0.0f;
    for (int k = 0; k < K; ++k) {
        float p_gc_high = 0.0f;
        for (int b = REFERENCE_GC_MIN; b < N_GC_BINS; ++b) {
            p_gc_high += static_cast<float>(p_gc[k][b]);
        }
        ref_num += result.pi[k] * p_gc_high * result.delta_max[k];
        ref_den += result.pi[k] * p_gc_high;
    }
    result.d_reference = ref_den > 0.01f ? ref_num / ref_den : 0.0f;

    // Compute BIC
    // Parameters: (K-1) mixing + K delta_max + K*9 GC categorical + 2 shared = 11K + 1
    int n_params = 11 * K + 1;
    double n_eff = 0.0;
    for (const auto& sr : super_reads) {
        for (int p = 0; p < N_POSITIONS; ++p) {
            n_eff += sr.n_tc[p] + sr.n_ag[p] + sr.n_stop[p];
        }
    }
    result.bic = static_cast<float>(-2.0 * result.log_likelihood + n_params * std::log(n_eff));

    return result;
}

inline MixtureDamageResult MixtureDamageModel::fit(
    const std::array<SuperRead, N_GC_BINS>& super_reads)
{
    MixtureDamageResult best_result;
    best_result.bic = std::numeric_limits<float>::infinity();

    // Try K = 2, 3, 4
    for (int K = 2; K <= 4; ++K) {
        MixtureDamageResult best_for_k;
        best_for_k.bic = std::numeric_limits<float>::infinity();

        // Multiple restarts
        for (int restart = 0; restart < N_RESTARTS; ++restart) {
            auto result = fit_k(super_reads, K, restart);
            if (result.converged && result.bic < best_for_k.bic) {
                best_for_k = result;
            }
        }

        // Keep best by BIC
        if (best_for_k.converged && best_for_k.bic < best_result.bic) {
            best_result = best_for_k;
        }
    }

    return best_result;
}

} // namespace agp
