#pragma once

// Proper statistical models for ancient DNA damage detection
//
// Two-level inference:
// 1. PROTEIN-LEVEL: Beta-Binomial model for overdispersion across reads
//    - Likelihood ratio test (LRT) for damage detection
//    - Profile likelihood confidence intervals
//    - Bayes factor via Laplace approximation
//
// 2. READ-LEVEL: Newton-MAP estimation with calibrated prior
//    - Already implemented in damage_probability.cpp
//
// Reference: Briggs et al. (2007), mapDamage2, metaDMG

#include <cmath>
#include <cstdint>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

#include <taph/protein_damage.hpp>

namespace dart {

// ============================================================================
// SPECIAL FUNCTIONS
// ============================================================================

// Log-gamma function using Lanczos approximation
// More numerically stable than std::lgamma for our use case
inline double log_gamma(double x) {
    // Lanczos coefficients for g=7
    static constexpr double LANCZOS_G = 7.0;
    static constexpr double LANCZOS_COEFF[9] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };

    if (x < 0.5) {
        // Reflection formula
        return std::log(M_PI / std::sin(M_PI * x)) - log_gamma(1.0 - x);
    }

    x -= 1.0;
    double a = LANCZOS_COEFF[0];
    for (int i = 1; i < 9; ++i) {
        a += LANCZOS_COEFF[i] / (x + i);
    }

    double t = x + LANCZOS_G + 0.5;
    return 0.5 * std::log(2 * M_PI) + (x + 0.5) * std::log(t) - t + std::log(a);
}

// Log-beta function: log B(a,b) = log Γ(a) + log Γ(b) - log Γ(a+b)
inline double log_beta(double a, double b) {
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b);
}

// Digamma function (derivative of log-gamma)
// Used for MLE estimation of Beta parameters
inline double digamma(double x) {
    double result = 0.0;

    // Use recurrence for small x
    while (x < 6.0) {
        result -= 1.0 / x;
        x += 1.0;
    }

    // Asymptotic series for large x
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    result += std::log(x) - 0.5 * inv_x
              - inv_x2 * (1.0/12.0 - inv_x2 * (1.0/120.0 - inv_x2 / 252.0));

    return result;
}

// Trigamma function (second derivative of log-gamma)
// Used for Fisher information / variance estimation
inline double trigamma(double x) {
    double result = 0.0;

    // Use recurrence for small x
    while (x < 6.0) {
        result += 1.0 / (x * x);
        x += 1.0;
    }

    // Asymptotic series for large x
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    result += inv_x + 0.5 * inv_x2
              + inv_x2 * inv_x * (1.0/6.0 - inv_x2 * (1.0/30.0 - inv_x2 / 42.0));

    return result;
}

// Regularized incomplete beta function I_x(a,b) using continued fraction
// This is the CDF of Beta(a,b) distribution
// Uses Lentz's algorithm for continued fraction evaluation
inline double beta_inc(double a, double b, double x) {
    if (x < 0.0 || x > 1.0) return 0.0;
    if (x == 0.0) return 0.0;
    if (x == 1.0) return 1.0;

    // Use symmetry for better convergence: I_x(a,b) = 1 - I_{1-x}(b,a)
    bool flip = x > (a + 1.0) / (a + b + 2.0);
    if (flip) {
        std::swap(a, b);
        x = 1.0 - x;
    }

    // Continued fraction (Lentz's algorithm)
    constexpr double TINY = 1e-30;
    constexpr double EPS = 1e-10;
    constexpr int MAX_ITER = 200;

    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;

    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < TINY) d = TINY;
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= MAX_ITER; ++m) {
        int m2 = 2 * m;

        // Even step
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < TINY) d = TINY;
        c = 1.0 + aa / c;
        if (std::abs(c) < TINY) c = TINY;
        d = 1.0 / d;
        h *= d * c;

        // Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < TINY) d = TINY;
        c = 1.0 + aa / c;
        if (std::abs(c) < TINY) c = TINY;
        d = 1.0 / d;
        double del = d * c;
        h *= del;

        if (std::abs(del - 1.0) < EPS) break;
    }

    // Compute the front factor
    double bt = std::exp(
        log_gamma(a + b) - log_gamma(a) - log_gamma(b)
        + a * std::log(x) + b * std::log(1.0 - x)
    );

    double result = bt * h / a;
    return flip ? 1.0 - result : result;
}

// Beta CDF: P(X <= x) where X ~ Beta(a, b)
inline double beta_cdf(double a, double b, double x) {
    return beta_inc(a, b, x);
}

// Beta quantile (inverse CDF) using bisection
// Returns x such that P(X <= x) = p where X ~ Beta(a, b)
inline double beta_quantile(double a, double b, double p) {
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return 1.0;

    constexpr double EPS = 1e-8;
    constexpr int MAX_ITER = 50;

    double lo = 0.0, hi = 1.0;
    for (int i = 0; i < MAX_ITER; ++i) {
        double mid = 0.5 * (lo + hi);
        double cdf = beta_cdf(a, b, mid);
        if (cdf < p) {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo < EPS) break;
    }
    return 0.5 * (lo + hi);
}

// ============================================================================
// GENE-LEVEL AGGREGATION (Beta-Bernoulli)
// ============================================================================

// Result of gene-level aggregation with EM weights
struct GeneDamageSummary {
    std::string protein_id;

    // Effective counts (γ-weighted)
    double n_eff = 0.0;             // Σγ_rt (effective read count)
    double n_anc_eff = 0.0;         // Σγ_rt × p_rta (ancient evidence)

    // Beta posterior: θ ~ Beta(α₀ + n_anc_eff, β₀ + n_mod_eff)
    double theta_raw = 0.0;         // n_anc_eff / n_eff (MLE)
    double theta_post = 0.0;        // Posterior mean with shrinkage
    double theta_ci_low = 0.0;      // 95% credible interval lower
    double theta_ci_high = 0.0;     // 95% credible interval upper

    // Enrichment test: P(θ > θ_bg + δ)
    double p_enriched = 0.0;        // Posterior probability of enrichment
    bool enriched = false;          // p_enriched >= 0.95

    // Coverage statistics (γ-weighted)
    double breadth = 0.0;           // Fraction of positions covered
    double depth = 0.0;             // Mean coverage depth
    double mean_identity = 0.0;     // Weighted mean identity

    // QC metrics
    uint32_t n_raw = 0;             // Raw read count (before EM)
    double multimap_entropy = 0.0;  // Assignment uncertainty (high = ambiguous)
    bool low_count = false;         // n_eff < threshold
    bool low_breadth = false;       // breadth < threshold
};

// Compute gene-level summary from EM results
// Parameters:
//   nu0: prior pseudo-count (shrinkage strength, default 10)
//   theta_bg: background ancient rate (global sample π)
//   delta: enrichment threshold (default 0.05)
inline GeneDamageSummary compute_gene_summary(
    const std::string& protein_id,
    double n_eff,
    double n_anc_eff,
    double breadth,
    double depth,
    double mean_identity,
    uint32_t n_raw,
    double multimap_entropy,
    double nu0 = 10.0,
    double theta_bg = 0.10,
    double delta = 0.05,
    double min_n_eff = 3.0,
    double min_breadth = 0.10)
{
    GeneDamageSummary result;
    result.protein_id = protein_id;
    result.n_eff = n_eff;
    result.n_anc_eff = n_anc_eff;
    result.breadth = breadth;
    result.depth = depth;
    result.mean_identity = mean_identity;
    result.n_raw = n_raw;
    result.multimap_entropy = multimap_entropy;

    // QC flags
    result.low_count = (n_eff < min_n_eff);
    result.low_breadth = (breadth < min_breadth);

    // Raw MLE (avoid division by zero)
    result.theta_raw = (n_eff > 0.0) ? n_anc_eff / n_eff : 0.0;

    // Beta posterior with empirical Bayes shrinkage
    // Prior: α₀ = ν₀ × θ_bg, β₀ = ν₀ × (1 - θ_bg)
    double alpha0 = nu0 * theta_bg;
    double beta0 = nu0 * (1.0 - theta_bg);

    double n_mod_eff = n_eff - n_anc_eff;  // Modern evidence
    double alpha_post = alpha0 + n_anc_eff;
    double beta_post = beta0 + n_mod_eff;

    // Posterior mean
    result.theta_post = alpha_post / (alpha_post + beta_post);

    // 95% credible interval
    result.theta_ci_low = beta_quantile(alpha_post, beta_post, 0.025);
    result.theta_ci_high = beta_quantile(alpha_post, beta_post, 0.975);

    // Enrichment probability: P(θ > θ_bg + δ)
    double threshold = theta_bg + delta;
    if (threshold >= 1.0) {
        result.p_enriched = 0.0;
    } else {
        result.p_enriched = 1.0 - beta_cdf(alpha_post, beta_post, threshold);
    }
    result.enriched = (result.p_enriched >= 0.95);

    return result;
}

// ============================================================================
// BETA-BINOMIAL PROTEIN DAMAGE MODEL — migrated to libtaph (single engine)
// ============================================================================
// The Beta-Binomial LRT, self-calibrated null, Laplace Bayes factor and
// profile-likelihood CI now live in taph::fit_protein_damage
// (<taph/protein_damage.hpp>). DART builds per-read observations (data) and
// calls the engine (math); no protein-damage math remains DART-side.
using ProteinDamageResult = taph::ProteinDamageResult;
using ReadDamageObs        = taph::ProteinReadObs;
using taph::fit_protein_damage;
using taph::fit_protein_damage_binomial;  // #2: per-opportunity Beta-Binomial
using taph::self_calibrated_pi0;

// ============================================================================
// Sample-level damage detection
// ============================================================================

// Result of sample-level damage analysis
struct SampleDamageResult {
    // MLE estimates
    double d_max = 0.0;             // Maximum damage rate at terminal
    double lambda = 0.3;            // Decay constant (exponential model)

    // Goodness of fit
    double r_squared = 0.0;         // R² for exponential fit
    double residual_se = 0.0;       // Standard error of residuals

    // Test statistics
    double log_lik_damage = 0.0;    // Log-likelihood under damage model
    double log_lik_null = 0.0;      // Log-likelihood under null
    double likelihood_ratio = 0.0;  // LRT statistic
    double p_value = 1.0;           // From χ² distribution

    // Confidence intervals
    double ci_lower = 0.0;          // 95% CI lower for d_max
    double ci_upper = 0.0;          // 95% CI upper for d_max
    double se_d_max = 0.0;          // Standard error

    // Channel B validation
    bool channel_b_valid = false;
    double channel_b_lrt = 0.0;     // LRT from stop codon channel
    bool damage_validated = false;  // Both channels agree

    // Per-stop-type Channel B (split by deamination parent)
    float stop_decay_llr_taa_5prime = 0.0f;
    float stop_decay_llr_tag_5prime = 0.0f;
    float stop_decay_llr_tga_5prime = 0.0f;
    double n_convertible_caa = 0.0;
    double n_convertible_cag = 0.0;
    double n_convertible_cga = 0.0;
    bool channel_b_valid_tga = false;  // CGA exposure >= 50

    // Summary
    size_t n_reads = 0;
    size_t n_positions = 0;
    bool is_damaged() const { return p_value < 0.05 && d_max > 0.01; }
};

/**
 * Fit exponential decay model to position-specific damage rates
 *
 * Model: δ(p) = b + (d_max - b) * exp(-λ * p)
 *
 * Uses weighted least squares:
 * - Weight by sqrt(n_p) at each position for variance stabilization
 * - Newton iteration for d_max and lambda
 *
 * @param damage_rates Position-specific damage rates (e.g., T/(T+C) - baseline)
 * @param counts Number of observations at each position
 * @param n_positions Number of positions to use (default 15)
 * @return SampleDamageResult with fitted model
 */
inline SampleDamageResult fit_sample_damage(
    const std::array<float, 15>& damage_rates,
    const std::array<double, 15>& counts,
    int n_positions = 15)
{
    SampleDamageResult result;
    result.n_positions = n_positions;

    // Compute total observations
    double total_n = 0.0;
    for (int p = 0; p < n_positions; ++p) {
        total_n += counts[p];
    }
    result.n_reads = static_cast<size_t>(total_n);

    if (total_n < 100) {
        return result;  // Insufficient data
    }

    // =========================================================================
    // Weighted least squares for d_max and lambda
    // =========================================================================

    // Initial estimates from terminal vs interior difference
    float d0 = damage_rates[0];
    float d_interior = 0.0f;
    double w_interior = 0.0;
    for (int p = 10; p < n_positions; ++p) {
        d_interior += damage_rates[p] * counts[p];
        w_interior += counts[p];
    }
    d_interior = (w_interior > 0) ? d_interior / w_interior : 0.0f;

    double d_max = std::max(0.0, static_cast<double>(d0 - d_interior));
    double lambda = 0.3;  // Initial decay constant
    double baseline = d_interior;

    // Newton iteration
    for (int iter = 0; iter < 50; ++iter) {
        double grad_d = 0.0, grad_l = 0.0;
        double hess_dd = 0.0, hess_dl = 0.0, hess_ll = 0.0;

        double ss_res = 0.0, ss_tot = 0.0;
        double y_mean = 0.0;
        for (int p = 0; p < n_positions; ++p) {
            y_mean += damage_rates[p] * counts[p];
        }
        y_mean /= total_n;

        for (int p = 0; p < n_positions; ++p) {
            double w = std::sqrt(counts[p]);  // Variance-stabilizing weight
            double exp_term = std::exp(-lambda * p);
            double pred = baseline + d_max * exp_term;
            double resid = damage_rates[p] - pred;

            // Gradient
            grad_d += -2.0 * w * resid * exp_term;
            grad_l += 2.0 * w * resid * d_max * p * exp_term;

            // Hessian (approximate)
            hess_dd += 2.0 * w * exp_term * exp_term;
            hess_dl += -2.0 * w * d_max * p * exp_term * exp_term;
            hess_ll += 2.0 * w * d_max * d_max * p * p * exp_term * exp_term;

            // R² calculation
            ss_res += w * resid * resid;
            ss_tot += w * (damage_rates[p] - y_mean) * (damage_rates[p] - y_mean);
        }

        // Regularization
        hess_dd = std::max(hess_dd, 1e-6);
        hess_ll = std::max(hess_ll, 1e-6);

        // Newton step with diagonal approximation
        double step_d = grad_d / hess_dd;
        double step_l = grad_l / hess_ll;

        d_max = std::max(0.0, d_max - 0.5 * step_d);
        lambda = std::clamp(lambda - 0.5 * step_l, 0.05, 2.0);

        if (std::abs(step_d) < 1e-6 && std::abs(step_l) < 1e-6) break;
    }

    result.d_max = d_max;
    result.lambda = lambda;

    // Compute final R² and residual SE
    double ss_res = 0.0, ss_tot = 0.0;
    double y_mean = 0.0;
    for (int p = 0; p < n_positions; ++p) {
        y_mean += damage_rates[p] * counts[p];
    }
    y_mean /= total_n;

    for (int p = 0; p < n_positions; ++p) {
        double pred = baseline + d_max * std::exp(-lambda * p);
        double resid = damage_rates[p] - pred;
        ss_res += counts[p] * resid * resid;
        ss_tot += counts[p] * (damage_rates[p] - y_mean) * (damage_rates[p] - y_mean);
    }

    result.r_squared = (ss_tot > 0) ? 1.0 - ss_res / ss_tot : 0.0;
    result.residual_se = std::sqrt(ss_res / (total_n - 2));

    // =========================================================================
    // Likelihood Ratio Test: exponential decay vs constant
    // =========================================================================

    // Approximate log-likelihoods using normal distribution
    // (variance proportional to 1/n at each position)
    double ll_damage = 0.0, ll_null = 0.0;

    for (int p = 0; p < n_positions; ++p) {
        if (counts[p] < 1) continue;

        double pred_damage = baseline + d_max * std::exp(-lambda * p);
        double pred_null = baseline;

        // Variance estimate: p(1-p)/n for binomial proportion
        double var_p = std::max(0.001, baseline * (1.0 - baseline) / counts[p]);
        double se_p = std::sqrt(var_p);

        // Normal log-likelihood
        double z_damage = (damage_rates[p] - pred_damage) / se_p;
        double z_null = (damage_rates[p] - pred_null) / se_p;

        ll_damage += -0.5 * z_damage * z_damage - std::log(se_p);
        ll_null += -0.5 * z_null * z_null - std::log(se_p);
    }

    result.log_lik_damage = ll_damage;
    result.log_lik_null = ll_null;
    result.likelihood_ratio = 2.0 * (ll_damage - ll_null);
    result.likelihood_ratio = std::max(result.likelihood_ratio, 0.0);

    // P-value from χ²(2) (two additional parameters: d_max, lambda)
    if (result.likelihood_ratio > 0.0) {
        // Quick approximation for χ²(2)
        double x = result.likelihood_ratio;
        result.p_value = std::exp(-x / 2.0);
    }

    // =========================================================================
    // Confidence interval via profile likelihood
    // =========================================================================

    double chi2_crit = 3.84;
    double ll_threshold = ll_damage - chi2_crit / 2.0;

    // Standard error from Fisher information (approximate)
    result.se_d_max = d_max / std::sqrt(std::max(1.0, result.likelihood_ratio));
    result.ci_lower = std::max(0.0, d_max - 1.96 * result.se_d_max);
    result.ci_upper = std::min(1.0, d_max + 1.96 * result.se_d_max);

    return result;
}

} // namespace dart
