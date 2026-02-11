#pragma once

#include <array>
#include <cmath>
#include <algorithm>

namespace agp {

/**
 * Joint Probabilistic Model for Ancient DNA Damage Detection
 *
 * Unifies Channel A (nucleotide frequencies), Control channel, and Channel B
 * (stop codon conversion) under a single Bayesian framework.
 *
 * Model:
 *   δ(p) = δ_max · exp(-λp)           -- damage rate at position p
 *   a(p) = a_max · exp(-λp)           -- artifact/composition rate (shared)
 *
 * Observations (binomial):
 *   Channel A: π_TC(p) = b_TC + a(p) + (1 - b_TC - a(p)) · δ(p)
 *   Control:   π_AG(p) = b_AG + a(p)
 *   Channel B: π_stop(p) = b_stop + (1 - b_stop) · δ(p)
 *
 * The control channel shares the artifact term but NOT the damage term,
 * which identifies composition artifacts vs real damage.
 */

struct JointDamageResult {
    // Fitted parameters
    float delta_max = 0.0f;      // Maximum damage rate at terminal (δ_max)
    float lambda = 0.0f;         // Decay constant
    float a_max = 0.0f;          // Artifact amplitude (can be negative for inverted)

    // Baselines (fixed from interior data)
    float b_tc = 0.0f;           // Interior T/(T+C) baseline
    float b_ag = 0.0f;           // Interior A/(A+G) baseline
    float b_stop = 0.0f;         // Interior stop/(pre+stop) baseline

    // Model comparison
    float log_lik_m1 = 0.0f;     // Log-likelihood for M1 (δ_max > 0)
    float log_lik_m0 = 0.0f;     // Log-likelihood for M0 (δ_max = 0)
    float bic_m1 = 0.0f;         // BIC for M1
    float bic_m0 = 0.0f;         // BIC for M0
    float delta_bic = 0.0f;      // BIC_M0 - BIC_M1 (positive favors damage)
    float bayes_factor = 0.0f;   // BF_10 ≈ exp(ΔBIC/2)
    float p_damage = 0.0f;       // P(damage | data) = BF / (1 + BF)

    // Diagnostics
    float rmse = 0.0f;           // Root mean squared error of fit
    int n_positions = 15;        // Number of positions used
    uint64_t n_trials = 0;       // Total binomial trials (for BIC)
    bool valid = false;          // Sufficient data for estimation
};

struct JointDamageSuffStats {
    // Channel A: T/(T+C) at 5' end
    std::array<uint64_t, 15> k_tc = {};   // T counts (successes)
    std::array<uint64_t, 15> n_tc = {};   // T+C counts (trials)

    // Control: A/(A+G) at 5' end (same artifact, no damage)
    std::array<uint64_t, 15> k_ag = {};   // A counts
    std::array<uint64_t, 15> n_ag = {};   // A+G counts

    // Channel B: stop/(pre+stop) conversions
    std::array<uint64_t, 15> k_stop = {}; // Stop codon counts
    std::array<uint64_t, 15> n_stop = {}; // Pre-stop + stop counts

    // Interior baselines
    uint64_t k_tc_interior = 0, n_tc_interior = 0;
    uint64_t k_ag_interior = 0, n_ag_interior = 0;
    uint64_t k_stop_interior = 0, n_stop_interior = 0;

    // Compute baselines from interior
    float baseline_tc() const {
        return n_tc_interior > 0 ? static_cast<float>(k_tc_interior) / n_tc_interior : 0.5f;
    }
    float baseline_ag() const {
        return n_ag_interior > 0 ? static_cast<float>(k_ag_interior) / n_ag_interior : 0.5f;
    }
    float baseline_stop() const {
        return n_stop_interior > 0 ? static_cast<float>(k_stop_interior) / n_stop_interior : 0.05f;
    }

    uint64_t total_trials() const {
        uint64_t total = 0;
        for (int p = 0; p < 15; ++p) {
            total += n_tc[p] + n_ag[p] + n_stop[p];
        }
        return total;
    }

    bool is_valid() const {
        uint64_t min_tc = *std::min_element(n_tc.begin(), n_tc.end());
        uint64_t min_ag = *std::min_element(n_ag.begin(), n_ag.end());
        return min_tc >= 100 && min_ag >= 100 && n_tc_interior >= 1000;
    }
};

/**
 * Fit the joint damage model using grid search
 *
 * Grid searches over (λ, δ_max) and optimizes a_max analytically for each point.
 * Returns the maximum likelihood estimate and BIC-based model comparison.
 */
class JointDamageModel {
public:
    // Grid parameters
    static constexpr int N_LAMBDA = 10;
    static constexpr int N_DELTA = 100;
    static constexpr float LAMBDA_MIN = 0.05f;
    static constexpr float LAMBDA_MAX = 0.50f;
    static constexpr float DELTA_MAX_LIMIT = 0.60f;

    /**
     * Fit the joint model to sufficient statistics
     */
    static JointDamageResult fit(const JointDamageSuffStats& stats);

    /**
     * Compute log-likelihood for given parameters
     */
    static float log_likelihood(
        const JointDamageSuffStats& stats,
        float delta_max, float lambda, float a_max,
        float b_tc, float b_ag, float b_stop);

    /**
     * Optimize a_max for fixed (δ_max, λ) using Newton's method
     */
    static float optimize_a_max(
        const JointDamageSuffStats& stats,
        float delta_max, float lambda,
        float b_tc, float b_ag);

private:
    // Binomial log-likelihood: k * log(p) + (n-k) * log(1-p)
    static float binom_ll(uint64_t k, uint64_t n, float p) {
        if (n == 0) return 0.0f;
        p = std::clamp(p, 1e-10f, 1.0f - 1e-10f);
        return static_cast<float>(k) * std::log(p) +
               static_cast<float>(n - k) * std::log(1.0f - p);
    }
};

// ============================================================================
// Implementation (header-only for simplicity)
// ============================================================================

inline float JointDamageModel::log_likelihood(
    const JointDamageSuffStats& stats,
    float delta_max, float lambda, float a_max,
    float b_tc, float b_ag, float b_stop)
{
    float ll = 0.0f;

    for (int p = 0; p < 15; ++p) {
        float decay = std::exp(-lambda * p);
        float delta_p = delta_max * decay;
        float a_p = a_max * decay;

        // Channel A: π_TC(p) = b_tc + a(p) + (1 - b_tc - a(p)) · δ(p)
        float base_tc = std::clamp(b_tc + a_p, 0.01f, 0.99f);
        float pi_tc = base_tc + (1.0f - base_tc) * delta_p;
        ll += binom_ll(stats.k_tc[p], stats.n_tc[p], pi_tc);

        // Control: π_AG(p) = b_ag + a(p)
        float pi_ag = std::clamp(b_ag + a_p, 0.01f, 0.99f);
        ll += binom_ll(stats.k_ag[p], stats.n_ag[p], pi_ag);

        // Channel B: π_stop(p) = b_stop + (1 - b_stop) · δ(p)
        float pi_stop = b_stop + (1.0f - b_stop) * delta_p;
        pi_stop = std::clamp(pi_stop, 0.001f, 0.999f);
        ll += binom_ll(stats.k_stop[p], stats.n_stop[p], pi_stop);
    }

    return ll;
}

inline float JointDamageModel::optimize_a_max(
    const JointDamageSuffStats& stats,
    float delta_max, float lambda,
    float b_tc, float b_ag)
{
    // Golden section search for a_max in [-0.3, 0.3]
    // (Faster than Newton for this 1D problem)

    const float phi = 0.618033988749895f;
    float a = -0.3f, b = 0.3f;
    float c = b - phi * (b - a);
    float d = a + phi * (b - a);

    auto eval = [&](float a_max) {
        float ll = 0.0f;
        for (int p = 0; p < 15; ++p) {
            float decay = std::exp(-lambda * p);
            float delta_p = delta_max * decay;
            float a_p = a_max * decay;

            // Channel A
            float base_tc = std::clamp(b_tc + a_p, 0.01f, 0.99f);
            float pi_tc = base_tc + (1.0f - base_tc) * delta_p;
            ll += binom_ll(stats.k_tc[p], stats.n_tc[p], pi_tc);

            // Control
            float pi_ag = std::clamp(b_ag + a_p, 0.01f, 0.99f);
            ll += binom_ll(stats.k_ag[p], stats.n_ag[p], pi_ag);
        }
        return ll;
    };

    for (int iter = 0; iter < 30; ++iter) {
        if (eval(c) > eval(d)) {
            b = d;
            d = c;
            c = b - phi * (b - a);
        } else {
            a = c;
            c = d;
            d = a + phi * (b - a);
        }
    }

    return (a + b) / 2.0f;
}

inline JointDamageResult JointDamageModel::fit(const JointDamageSuffStats& stats) {
    JointDamageResult result;

    // Get baselines from interior
    result.b_tc = stats.baseline_tc();
    result.b_ag = stats.baseline_ag();
    result.b_stop = stats.baseline_stop();
    result.n_trials = stats.total_trials();
    result.valid = stats.is_valid();

    if (!result.valid) {
        return result;
    }

    // Grid search over (λ, δ_max)
    float best_ll = -1e30f;
    float best_delta = 0.0f;
    float best_lambda = 0.2f;
    float best_a = 0.0f;

    for (int i_lambda = 0; i_lambda < N_LAMBDA; ++i_lambda) {
        float lambda = LAMBDA_MIN + (LAMBDA_MAX - LAMBDA_MIN) * i_lambda / (N_LAMBDA - 1);

        for (int i_delta = 0; i_delta <= N_DELTA; ++i_delta) {
            float delta_max = DELTA_MAX_LIMIT * i_delta / N_DELTA;

            // Optimize a_max for this (λ, δ_max)
            float a_max = optimize_a_max(stats, delta_max, lambda,
                                         result.b_tc, result.b_ag);

            float ll = log_likelihood(stats, delta_max, lambda, a_max,
                                      result.b_tc, result.b_ag, result.b_stop);

            if (ll > best_ll) {
                best_ll = ll;
                best_delta = delta_max;
                best_lambda = lambda;
                best_a = a_max;
            }
        }
    }

    result.delta_max = best_delta;
    result.lambda = best_lambda;
    result.a_max = best_a;
    result.log_lik_m1 = best_ll;

    // Compute M0 likelihood (δ_max = 0)
    float a_max_m0 = optimize_a_max(stats, 0.0f, best_lambda,
                                    result.b_tc, result.b_ag);
    result.log_lik_m0 = log_likelihood(stats, 0.0f, best_lambda, a_max_m0,
                                       result.b_tc, result.b_ag, result.b_stop);

    // BIC comparison
    // M1 has 3 parameters (δ_max, λ, a_max), M0 has 2 (λ, a_max)
    // BIC = -2 * log_lik + k * log(N)
    float log_n = std::log(static_cast<float>(result.n_trials));
    result.bic_m1 = -2.0f * result.log_lik_m1 + 3.0f * log_n;
    result.bic_m0 = -2.0f * result.log_lik_m0 + 2.0f * log_n;

    // ΔBIC = BIC_M0 - BIC_M1 (positive favors M1 = damage)
    result.delta_bic = result.bic_m0 - result.bic_m1;

    // Bayes factor approximation: BF_10 ≈ exp(ΔBIC/2)
    result.bayes_factor = std::exp(result.delta_bic / 2.0f);

    // P(damage) = BF / (1 + BF) with prior P(damage) = 0.5
    result.p_damage = result.bayes_factor / (1.0f + result.bayes_factor);

    // Compute RMSE for diagnostics
    float sse = 0.0f;
    int n_obs = 0;
    for (int p = 0; p < 15; ++p) {
        if (stats.n_tc[p] > 0) {
            float decay = std::exp(-result.lambda * p);
            float base_tc = std::clamp(result.b_tc + result.a_max * decay, 0.01f, 0.99f);
            float pi_tc = base_tc + (1.0f - base_tc) * result.delta_max * decay;
            float obs_tc = static_cast<float>(stats.k_tc[p]) / stats.n_tc[p];
            sse += (pi_tc - obs_tc) * (pi_tc - obs_tc);
            ++n_obs;
        }
    }
    result.rmse = n_obs > 0 ? std::sqrt(sse / n_obs) : 0.0f;
    result.n_positions = n_obs;

    return result;
}

// ============================================================================
// TWO-COMPONENT MIXTURE MODEL FOR GC-STRATIFIED DAMAGE
// Separates damaged/undamaged subpopulations in heterogeneous metagenomes
// ============================================================================

struct DamageMixtureResult {
    float d_mean = 0.0f;         // GC-weighted average (reference-free population avg)
    float d_ancient = 0.0f;      // Damaged component mean (μ_1)
    float pi_ancient = 0.0f;     // Fraction of C-sites in damaged component
    float tau_ancient = 0.0f;    // Damaged component std dev
    int n_iterations = 0;        // EM iterations to converge
    bool converged = false;      // Did EM converge?
    bool separated = false;      // Are components well-separated? (d_ancient > 0.02)
};

struct GCBinInput {
    float d_max;      // Estimated damage for this bin
    float c_sites;    // Weight (number of C sites)
    bool valid;       // Is this bin valid?
};

class DamageMixtureModel {
public:
    // Fixed parameters
    static constexpr float MU_0 = 0.0f;       // Undamaged component mean
    static constexpr float TAU_0 = 0.01f;     // Undamaged component std dev
    static constexpr float TAU_1 = 0.10f;     // Damaged component std dev (fixed)
    static constexpr float SIGMA_FLOOR = 0.02f; // Minimum observation noise
    static constexpr int MAX_ITER = 50;
    static constexpr float CONVERGENCE_TOL = 1e-6f;

    /**
     * Fit 2-component mixture to GC-stratified damage estimates
     *
     * @param bins Array of GC bin damage estimates
     * @param n_bins Number of bins (typically 10)
     * @return Mixture model results
     */
    template<size_t N>
    static DamageMixtureResult fit(const std::array<GCBinInput, N>& bins) {
        DamageMixtureResult result;

        // Compute weighted mean (d_mean)
        float total_weight = 0.0f;
        float weighted_sum = 0.0f;
        int n_valid = 0;
        float max_d = 0.0f;

        for (size_t i = 0; i < N; ++i) {
            if (bins[i].valid && bins[i].c_sites > 0) {
                weighted_sum += bins[i].d_max * bins[i].c_sites;
                total_weight += bins[i].c_sites;
                max_d = std::max(max_d, bins[i].d_max);
                ++n_valid;
            }
        }

        if (total_weight < 1.0f || n_valid < 2) {
            return result;  // Not enough data
        }

        result.d_mean = weighted_sum / total_weight;

        // Initialize EM
        float pi_1 = 0.2f;  // Initial mixing proportion for damaged
        float mu_1 = std::max(0.05f, max_d * 0.8f);  // Initialize near max

        // Precompute observation noise σ_i for each bin
        std::array<float, N> sigma;
        for (size_t i = 0; i < N; ++i) {
            if (bins[i].valid && bins[i].c_sites > 0) {
                float d = std::clamp(bins[i].d_max, 0.001f, 0.999f);
                float sigma_binom = std::sqrt(d * (1.0f - d) / bins[i].c_sites);
                sigma[i] = std::max(SIGMA_FLOOR, sigma_binom);
            } else {
                sigma[i] = SIGMA_FLOOR;
            }
        }

        // EM iterations
        std::array<float, N> r;  // Responsibilities (posterior prob of damaged)
        float prev_ll = -1e30f;

        for (int iter = 0; iter < MAX_ITER; ++iter) {
            // E-step: compute responsibilities using log-space for stability
            float ll = 0.0f;
            for (size_t i = 0; i < N; ++i) {
                if (!bins[i].valid || bins[i].c_sites < 1.0f) {
                    r[i] = 0.0f;
                    continue;
                }

                float d = bins[i].d_max;
                float w = bins[i].c_sites;

                // Variance for each component: τ_k² + σ_i²
                float var_0 = TAU_0 * TAU_0 + sigma[i] * sigma[i];
                float var_1 = TAU_1 * TAU_1 + sigma[i] * sigma[i];

                // Log-likelihood under each component
                float ll_0 = -0.5f * std::log(var_0) - 0.5f * (d - MU_0) * (d - MU_0) / var_0;
                float ll_1 = -0.5f * std::log(var_1) - 0.5f * (d - mu_1) * (d - mu_1) / var_1;

                // Log-posterior (unnormalized)
                float log_p0 = std::log(1.0f - pi_1 + 1e-10f) + ll_0;
                float log_p1 = std::log(pi_1 + 1e-10f) + ll_1;

                // Normalize using log-sum-exp
                float log_sum = log_p0;
                if (log_p1 > log_p0) {
                    log_sum = log_p1 + std::log(1.0f + std::exp(log_p0 - log_p1));
                } else {
                    log_sum = log_p0 + std::log(1.0f + std::exp(log_p1 - log_p0));
                }

                r[i] = std::exp(log_p1 - log_sum);
                ll += w * log_sum;
            }

            // Check convergence
            if (std::abs(ll - prev_ll) < CONVERGENCE_TOL * std::abs(prev_ll)) {
                result.converged = true;
                result.n_iterations = iter + 1;
                break;
            }
            prev_ll = ll;

            // M-step: update π_1 and μ_1
            float sum_w_r = 0.0f;
            float sum_w_r_d = 0.0f;
            float sum_w = 0.0f;

            for (size_t i = 0; i < N; ++i) {
                if (bins[i].valid && bins[i].c_sites > 0) {
                    float w = bins[i].c_sites;
                    sum_w += w;
                    sum_w_r += w * r[i];
                    sum_w_r_d += w * r[i] * bins[i].d_max;
                }
            }

            if (sum_w_r > 1e-10f) {
                pi_1 = sum_w_r / sum_w;
                mu_1 = std::clamp(sum_w_r_d / sum_w_r, 0.0f, 1.0f);
            } else {
                // No mass in damaged component - declare no separation
                result.converged = true;
                result.n_iterations = iter + 1;
                result.separated = false;
                result.d_ancient = 0.0f;
                result.pi_ancient = 0.0f;
                return result;
            }

            result.n_iterations = iter + 1;
        }

        // Store results
        result.d_ancient = mu_1;
        result.pi_ancient = pi_1;
        result.tau_ancient = TAU_1;
        result.separated = (mu_1 > 0.02f && pi_1 > 0.01f);

        return result;
    }
};

} // namespace agp
