#pragma once
// EM multi-mapping resolution with SQUAREM acceleration
//
// Resolves ambiguous protein-to-reference alignments using an EM algorithm
// that jointly models mapping weights and ancient damage status.
//
// Core idea: reads that map to multiple references are assigned fractional
// weights based on bit score, damage evidence, and iteratively-refined
// reference abundance estimates.
//
// SQUAREM (Varadhan & Roland 2008) accelerates convergence from ~50-100
// iterations to ~10-20 by extrapolating along the EM fixed-point map.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace agp {

// EM algorithm parameters
struct EMParams {
    double lambda_b = 0.5;       // temperature for bit score: exp(lambda_b * score)
    uint32_t max_iters = 200;    // maximum EM iterations
    double tol = 1e-6;           // convergence tolerance (relative log-likelihood change)
    double min_weight = 1e-10;   // floor for reference weights (prevents log(0))
    bool use_squarem = true;     // enable SQUAREM acceleration
    bool use_damage = true;      // incorporate damage evidence in responsibilities

    // Rich-get-richer prevention
    double alpha_prior = 1.0;    // Dirichlet pseudocount (1.0 = uniform prior)
                                 // Higher values shrink weights toward uniform
                                 // 0 = no regularization (original behavior)
    bool normalize_by_length = false;  // Divide weights by sqrt(ref_length)
                                       // Prevents longer references from dominating
};

// Per-alignment record: one read mapping to one reference
struct Alignment {
    uint32_t read_idx;           // index into reads
    uint32_t ref_idx;            // index into references
    float bit_score;             // alignment bit score
    float damage_score;          // Bayesian damage posterior (0-1), from AGD index
    float fident;                // alignment identity fraction (0-1)
};

// Sparse alignment storage: CSR-like layout grouped by read
struct AlignmentData {
    std::vector<Alignment> alignments;     // all alignments, sorted by read_idx
    std::vector<uint32_t> read_offsets;    // read_offsets[r] = first alignment index for read r
                                           // read_offsets[num_reads] = alignments.size()
    std::vector<uint32_t> ref_lengths;     // ref_lengths[t] = length of reference t (for normalization)
    uint32_t num_reads = 0;
    uint32_t num_refs = 0;

    // Number of alignments for read r
    uint32_t degree(uint32_t r) const {
        return read_offsets[r + 1] - read_offsets[r];
    }

    // Alignments for read r
    const Alignment* begin(uint32_t r) const {
        return alignments.data() + read_offsets[r];
    }
    const Alignment* end(uint32_t r) const {
        return alignments.data() + read_offsets[r + 1];
    }
};

// EM state vector
struct EMState {
    std::vector<double> weights;           // w_t: reference abundance weights (sum = 1)
    std::vector<double> gamma;             // gamma_{rt}: P(ref_t | read_r), flat array [aln_idx]
    std::vector<double> gamma_ancient;     // P(ancient, ref_t | read_r) (damage extension)
    double pi = 0.10;                      // global ancient fraction estimate
    double log_likelihood = -std::numeric_limits<double>::infinity();
    uint32_t iterations = 0;
};

// Numerically stable log-sum-exp
inline double log_sum_exp(const double* vals, size_t n) {
    if (n == 0) return -std::numeric_limits<double>::infinity();
    if (n == 1) return vals[0];

    double max_val = vals[0];
    for (size_t i = 1; i < n; ++i) {
        if (vals[i] > max_val) max_val = vals[i];
    }

    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double diff = vals[i] - max_val;
        sum += (diff > -100.0) ? std::exp(diff) : 0.0;
    }
    return max_val + std::log(sum);
}

// E-step: compute responsibilities gamma_{rt} = w_t * exp(lambda * b_{rt}) / Z_r
// Returns log-likelihood sum_r log(Z_r)
inline double e_step(
    const AlignmentData& data,
    const EMParams& params,
    const double* weights,
    double* gamma)
{
    double ll = 0.0;

    // Scratch buffer for log-scores per read (reused across reads)
    // Max degree is bounded by data, but we allocate per-read to avoid global max scan
    std::vector<double> log_scores;

    for (uint32_t r = 0; r < data.num_reads; ++r) {
        const uint32_t start = data.read_offsets[r];
        const uint32_t end = data.read_offsets[r + 1];
        const uint32_t deg = end - start;
        if (deg == 0) continue;

        log_scores.resize(deg);

        // Compute unnormalized log-responsibilities
        for (uint32_t j = 0; j < deg; ++j) {
            const auto& aln = data.alignments[start + j];
            log_scores[j] = std::log(std::max(weights[aln.ref_idx], params.min_weight))
                           + params.lambda_b * static_cast<double>(aln.bit_score);
        }

        // Normalize via log-sum-exp
        double lse = log_sum_exp(log_scores.data(), deg);
        ll += lse;

        for (uint32_t j = 0; j < deg; ++j) {
            double diff = log_scores[j] - lse;
            gamma[start + j] = (diff > -100.0) ? std::exp(diff) : 0.0;
        }
    }

    return ll;
}

// M-step: update weights with Dirichlet prior and optional length normalization
// w_t ∝ (α₀ + Σ_r γ_rt) / sqrt(len_t)   [if normalize_by_length]
// w_t ∝ (α₀ + Σ_r γ_rt)                  [otherwise]
//
// Dirichlet prior prevents "rich get richer": references can't collapse to 0 weight,
// and extreme weights are shrunk toward uniform. α₀=1 is symmetric uniform prior.
//
// Length normalization prevents longer references from dominating multi-mappers
// simply because they have more opportunity for matches.
inline void m_step(
    const AlignmentData& data,
    const EMParams& params,
    const double* gamma,
    double* out_weights)
{
    const uint32_t T = data.num_refs;

    // Initialize with Dirichlet prior pseudocounts
    // α₀ = 1.0 gives uniform prior; higher values = stronger shrinkage toward uniform
    for (uint32_t t = 0; t < T; ++t) {
        out_weights[t] = params.alpha_prior;
    }

    // Accumulate responsibilities
    const size_t total = data.alignments.size();
    for (size_t i = 0; i < total; ++i) {
        out_weights[data.alignments[i].ref_idx] += gamma[i];
    }

    // Optional length normalization: divide by sqrt(ref_length)
    // This prevents longer references from attracting more multi-mappers
    if (params.normalize_by_length && !data.ref_lengths.empty()) {
        for (uint32_t t = 0; t < T; ++t) {
            double len = static_cast<double>(std::max(data.ref_lengths[t], 1u));
            out_weights[t] /= std::sqrt(len);
        }
    }

    // Apply floor and normalize to simplex
    double sum = 0.0;
    for (uint32_t t = 0; t < T; ++t) {
        out_weights[t] = std::max(out_weights[t], params.min_weight);
        sum += out_weights[t];
    }
    double inv_sum = 1.0 / sum;
    #pragma omp simd
    for (uint32_t t = 0; t < T; ++t) {
        out_weights[t] *= inv_sum;
    }
}

// Single EM fixed-point map: E-step + M-step
// Takes weights_in, writes weights_out, returns log-likelihood
inline double em_fixed_point(
    const AlignmentData& data,
    const EMParams& params,
    const double* weights_in,
    double* weights_out,
    double* gamma)
{
    double ll = e_step(data, params, weights_in, gamma);
    m_step(data, params, gamma, weights_out);
    return ll;
}

// SQUAREM acceleration step (Varadhan & Roland 2008)
// Extrapolates along the EM map direction for faster convergence.
// x: current point, x1: one EM step, x2: two EM steps
// r, v: scratch buffers (size n)
// Modifies x in-place to the accelerated point.
inline void squarem_step(
    double* x,
    const double* x1,
    const double* x2,
    double* r,
    double* v,
    size_t n)
{
    double r_norm_sq = 0.0;
    double v_norm_sq = 0.0;

    #pragma omp simd reduction(+:r_norm_sq, v_norm_sq)
    for (size_t i = 0; i < n; ++i) {
        r[i] = x1[i] - x[i];
        v[i] = (x2[i] - x1[i]) - r[i];
        r_norm_sq += r[i] * r[i];
        v_norm_sq += v[i] * v[i];
    }

    // Step length: alpha = -||r|| / ||v||
    // Clamped to [-1, -0.01] to prevent divergence
    double alpha = -std::sqrt(r_norm_sq / std::max(v_norm_sq, 1e-15));
    alpha = std::clamp(alpha, -1.0, -0.01);

    // Extrapolate: x_new = x - 2*alpha*r + alpha^2*v
    double alpha_sq = alpha * alpha;
    double neg_2alpha = -2.0 * alpha;

    #pragma omp simd
    for (size_t i = 0; i < n; ++i) {
        x[i] = x[i] + neg_2alpha * r[i] + alpha_sq * v[i];
        x[i] = std::max(x[i], 1e-15);
    }

    // Re-normalize to simplex
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) sum += x[i];
    double inv_sum = 1.0 / sum;
    #pragma omp simd
    for (size_t i = 0; i < n; ++i) {
        x[i] *= inv_sum;
    }
}

// Damage-aware E-step extension
// Computes joint responsibility: P(ref_t, ancient | read_r)
// Uses damage_score as additional evidence for ancient status
inline double e_step_damage(
    const AlignmentData& data,
    const EMParams& params,
    const double* weights,
    double pi_ancient,
    double* gamma,
    double* gamma_ancient)
{
    double ll = 0.0;
    const double eps = 1e-12;
    const double pi_a = std::clamp(pi_ancient, eps, 1.0 - eps);

    std::vector<double> log_scores;

    for (uint32_t r = 0; r < data.num_reads; ++r) {
        const uint32_t start = data.read_offsets[r];
        const uint32_t end = data.read_offsets[r + 1];
        const uint32_t deg = end - start;
        if (deg == 0) continue;

        // For damage-aware: 2*deg entries (ancient/modern for each alignment)
        log_scores.resize(2 * deg);

        for (uint32_t j = 0; j < deg; ++j) {
            const auto& aln = data.alignments[start + j];
            double log_w = std::log(std::max(weights[aln.ref_idx], params.min_weight));
            double score_contrib = params.lambda_b * static_cast<double>(aln.bit_score);

            // P(damage_score | ancient) and P(damage_score | modern)
            // Simple model: damage_score is a noisy indicator
            double ds = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
            double log_p_ds_ancient = std::log(ds);           // higher damage -> more likely ancient
            double log_p_ds_modern = std::log(1.0 - ds);      // lower damage -> more likely modern

            // Joint: P(ref_t, ancient | read_r) prop_to w_t * exp(lambda*b) * pi_a * P(ds|A)
            log_scores[2 * j]     = log_w + score_contrib + std::log(pi_a) + log_p_ds_ancient;
            log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - pi_a) + log_p_ds_modern;
        }

        double lse = log_sum_exp(log_scores.data(), 2 * deg);
        ll += lse;

        for (uint32_t j = 0; j < deg; ++j) {
            double g_ancient = 0.0, g_modern = 0.0;
            double diff_a = log_scores[2 * j] - lse;
            double diff_m = log_scores[2 * j + 1] - lse;
            g_ancient = (diff_a > -100.0) ? std::exp(diff_a) : 0.0;
            g_modern  = (diff_m > -100.0) ? std::exp(diff_m) : 0.0;

            gamma[start + j] = g_ancient + g_modern;   // marginal over ancient status
            gamma_ancient[start + j] = g_ancient;       // joint (ref_t AND ancient)
        }
    }

    return ll;
}

// Full EM solver with optional SQUAREM acceleration
inline EMState em_solve(
    const AlignmentData& data,
    const EMParams& params)
{
    const uint32_t T = data.num_refs;
    const size_t A = data.alignments.size();

    EMState state;
    state.weights.assign(T, 1.0 / static_cast<double>(T));  // uniform init
    state.gamma.resize(A, 0.0);

    if (params.use_damage) {
        state.gamma_ancient.resize(A, 0.0);
    }

    // Scratch buffers for SQUAREM
    std::vector<double> w0, w1, w2, r_buf, v_buf;
    if (params.use_squarem) {
        w0.resize(T);
        w1.resize(T);
        w2.resize(T);
        r_buf.resize(T);
        v_buf.resize(T);
    }

    double prev_ll = -std::numeric_limits<double>::infinity();

    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        double ll;

        if (params.use_squarem) {
            // Save current weights
            std::copy(state.weights.begin(), state.weights.end(), w0.begin());

            // First EM step: w0 -> w1
            if (params.use_damage) {
                e_step_damage(data, params, w0.data(), state.pi, state.gamma.data(), state.gamma_ancient.data());
            } else {
                e_step(data, params, w0.data(), state.gamma.data());
            }
            m_step(data, params, state.gamma.data(), w1.data());

            // Second EM step: w1 -> w2
            if (params.use_damage) {
                e_step_damage(data, params, w1.data(), state.pi, state.gamma.data(), state.gamma_ancient.data());
            } else {
                e_step(data, params, w1.data(), state.gamma.data());
            }
            m_step(data, params, state.gamma.data(), w2.data());

            // SQUAREM extrapolation: modifies w0 in-place
            squarem_step(w0.data(), w1.data(), w2.data(), r_buf.data(), v_buf.data(), T);

            // Evaluate at accelerated point
            if (params.use_damage) {
                ll = e_step_damage(data, params, w0.data(), state.pi, state.gamma.data(), state.gamma_ancient.data());
            } else {
                ll = e_step(data, params, w0.data(), state.gamma.data());
            }

            // If acceleration worsened likelihood, fall back to w2
            if (ll < prev_ll && iter > 0) {
                std::copy(w2.begin(), w2.end(), state.weights.begin());
                if (params.use_damage) {
                    ll = e_step_damage(data, params, state.weights.data(), state.pi, state.gamma.data(), state.gamma_ancient.data());
                } else {
                    ll = e_step(data, params, state.weights.data(), state.gamma.data());
                }
            } else {
                std::copy(w0.begin(), w0.end(), state.weights.begin());
            }
        } else {
            // Plain EM
            if (params.use_damage) {
                ll = e_step_damage(data, params, state.weights.data(), state.pi, state.gamma.data(), state.gamma_ancient.data());
            } else {
                ll = e_step(data, params, state.weights.data(), state.gamma.data());
            }
            m_step(data, params, state.gamma.data(), state.weights.data());
        }

        // Update pi (ancient fraction) from damage responsibilities
        if (params.use_damage && !state.gamma_ancient.empty()) {
            double sum_ancient = 0.0;
            double sum_total = 0.0;
            for (size_t i = 0; i < A; ++i) {
                sum_ancient += state.gamma_ancient[i];
                sum_total += state.gamma[i];
            }
            if (sum_total > 0.0) {
                state.pi = std::clamp(sum_ancient / sum_total, 0.01, 0.99);
            }
        }

        state.log_likelihood = ll;
        state.iterations = iter + 1;

        // Convergence check
        if (iter > 0) {
            double rel_change = std::abs(ll - prev_ll) / (std::abs(prev_ll) + 1e-15);
            if (rel_change < params.tol) break;
        }
        prev_ll = ll;
    }

    return state;
}

// Extract best assignment per read from converged EM state
// Returns vector of (ref_idx, responsibility) for each read
inline std::vector<std::pair<uint32_t, float>> extract_best_assignments(
    const AlignmentData& data,
    const EMState& state)
{
    std::vector<std::pair<uint32_t, float>> result(data.num_reads, {0, 0.0f});

    for (uint32_t r = 0; r < data.num_reads; ++r) {
        const uint32_t start = data.read_offsets[r];
        const uint32_t end = data.read_offsets[r + 1];

        double best_gamma = -1.0;
        uint32_t best_ref = 0;

        for (uint32_t j = start; j < end; ++j) {
            if (state.gamma[j] > best_gamma) {
                best_gamma = state.gamma[j];
                best_ref = data.alignments[j].ref_idx;
            }
        }

        result[r] = {best_ref, static_cast<float>(best_gamma)};
    }

    return result;
}

// Filter alignments: keep only those with gamma >= threshold
inline std::vector<uint32_t> filter_alignments(
    const AlignmentData& data,
    const EMState& state,
    double gamma_threshold = 0.01)
{
    std::vector<uint32_t> kept;
    kept.reserve(data.alignments.size());

    for (size_t i = 0; i < data.alignments.size(); ++i) {
        if (state.gamma[i] >= gamma_threshold) {
            kept.push_back(static_cast<uint32_t>(i));
        }
    }

    return kept;
}

// Parallel implementations (defined in em_reassign.cpp)
struct CompactAlignment;
class DamageIndexReader;
class ReferenceStatsCollector;

AlignmentData build_alignment_data(
    const CompactAlignment* alns,
    size_t n_alns,
    uint32_t num_reads,
    uint32_t num_refs,
    const DamageIndexReader* agd = nullptr);

EMState squarem_em(
    const AlignmentData& data,
    const EMParams& params,
    ReferenceStatsCollector* stats_collector = nullptr);

std::vector<std::pair<uint32_t, float>> reassign_reads(
    const AlignmentData& data,
    const EMState& state,
    double gamma_threshold = 0.01);

} // namespace agp
