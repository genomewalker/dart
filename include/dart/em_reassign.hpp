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
#include <functional>
#include <limits>
#include <vector>

namespace dart {

// EM algorithm parameters
struct EMParams {
    double lambda_b = 0.5;       // temperature for bit score: exp(lambda_b * score)
    uint32_t max_iters = 200;    // maximum EM iterations
    uint32_t min_iters = 8;      // minimum iterations before early-stop is allowed
    double tol = 1e-6;           // convergence tolerance (relative log-likelihood change)
    double min_weight = 1e-10;   // floor for reference weights (prevents log(0))
    bool use_squarem = true;     // enable SQUAREM acceleration
    uint32_t squarem_warmup_iters = 4;  // run plain EM for first N iterations, then enable SQUAREM
    int squarem_method = 3;      // SQUAREM step length: 1, 2, or 3 (paper defaults to 3)
    double squarem_step_min0 = 1.0;   // initial step.min (EM/MM contraction maps typically +1)
    double squarem_step_max0 = 1.0;   // initial step.max
    double squarem_mstep = 4.0;       // multiplicative growth/shrink for adaptive step bounds
    double squarem_objfn_inc = 1.0;   // allowed objective decrease before fallback
    bool use_damage = true;      // incorporate damage evidence in responsibilities
    bool use_alignment_damage_likelihood = false;  // if true use ll_a/ll_m with per-read prior

    // Dirichlet prior parameter (MAP): alpha=1 gives no prior effect.
    // alpha>1 shrinks toward uniform; alpha<1 promotes sparsity.
    double alpha_prior = 1.0;
    bool normalize_by_length = false;  // Apply -0.5*log(ref_length) in E-step scores
                                       // Penalizes longer references in a proper objective

    // Optional initial weights for warm-starting EM (e.g., from coverage-aware outer loop).
    // If non-empty, must have size == num_refs. EM initializes from these instead of uniform.
    std::vector<double> initial_weights;
};

// Per-alignment record: one read mapping to one reference
struct Alignment {
    uint32_t read_idx;           // index into reads
    uint32_t ref_idx;            // index into references
    float bit_score;             // alignment bit score
    float damage_score;          // per-read prior p_read (0-1), copied to alignment rows
    float damage_ll_a;           // alignment-level log P(evidence | ancient)
    float damage_ll_m;           // alignment-level log P(evidence | modern)
    float fident;                // alignment identity fraction (0-1)
    uint16_t tstart = 0;         // 0-based alignment start on reference (for coverage-EM)
    uint16_t aln_len = 0;        // alignment length on reference (for coverage-EM)
};

// Sparse alignment storage: CSR-like layout grouped by read
struct AlignmentData {
    std::vector<Alignment> alignments;     // all alignments, sorted by read_idx
    std::vector<uint32_t> read_offsets;    // read_offsets[r] = first alignment index for read r
                                           // read_offsets[num_reads] = alignments.size()
    std::vector<uint32_t> ref_lengths;     // ref_lengths[t] = length of reference t (for normalization)
    std::vector<float> ref_log_len_half;   // 0.5 * log(max(ref_length, 1))
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

inline double length_log_penalty(
    const AlignmentData& data,
    const EMParams& params,
    uint32_t ref_idx)
{
    if (!params.normalize_by_length) return 0.0;
    if (!data.ref_log_len_half.empty()) {
        return -static_cast<double>(data.ref_log_len_half[ref_idx]);
    }
    if (!data.ref_lengths.empty()) {
        const double len = static_cast<double>(std::max(data.ref_lengths[ref_idx], 1u));
        return -0.5 * std::log(len);
    }
    return 0.0;
}

// EM state vector
struct EMState {
    std::vector<double> weights;           // w_t: reference abundance weights (sum = 1)
    std::vector<double> gamma;             // gamma_{rt}: P(ref_t | read_r), flat array [aln_idx]
    std::vector<double> gamma_ancient;     // P(ancient, ref_t | read_r) (damage extension)
    double pi = 0.10;                      // global ancient fraction estimate
    double log_likelihood = -std::numeric_limits<double>::infinity();
    uint32_t iterations = 0;
};

// Per-iteration diagnostics emitted by the EM solver.
struct EMIterationDiagnostics {
    uint32_t iteration = 0;
    double log_likelihood = -std::numeric_limits<double>::infinity();
    double objective = -std::numeric_limits<double>::infinity();
    double rel_change = std::numeric_limits<double>::infinity();
    double param_change = std::numeric_limits<double>::infinity();
    bool squarem_active = false;
    bool fallback_to_em = false;
    double alpha = 1.0;
    double step_min = 0.0;
    double step_max = 0.0;
};

using EMProgressCallback = std::function<void(const EMIterationDiagnostics&)>;

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
        sum += (diff > -700.0) ? std::exp(diff) : 0.0;
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
                           + length_log_penalty(data, params, aln.ref_idx)
                           + params.lambda_b * static_cast<double>(aln.bit_score);
        }

        // Normalize via log-sum-exp
        double lse = log_sum_exp(log_scores.data(), deg);
        ll += lse;

        for (uint32_t j = 0; j < deg; ++j) {
            double diff = log_scores[j] - lse;
            gamma[start + j] = (diff > -700.0) ? std::exp(diff) : 0.0;
        }
    }

    return ll;
}

// M-step: MAP update with symmetric Dirichlet(alpha_prior)
// w_t ∝ (alpha_prior - 1 + Σ_r γ_rt), then projected to simplex with floor.
//
// alpha_prior = 1.0 reduces to standard MLE update.
inline void m_step(
    const AlignmentData& data,
    const EMParams& params,
    const double* gamma,
    double* out_weights)
{
    const uint32_t T = data.num_refs;
    const double prior_offset = params.alpha_prior - 1.0;

    // Initialize with prior offset used by MAP update.
    for (uint32_t t = 0; t < T; ++t) {
        out_weights[t] = prior_offset;
    }

    // Accumulate responsibilities
    const size_t total = data.alignments.size();
    for (size_t i = 0; i < total; ++i) {
        out_weights[data.alignments[i].ref_idx] += gamma[i];
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

// MAP-like objective used by SQUAREM safeguard/convergence:
// objective = log-likelihood + (alpha_prior - 1) * sum_t log(w_t)
// For alpha_prior <= 0 we skip prior term and fall back to log-likelihood.
inline double em_objective_from_ll(
    const EMParams& params,
    const double* weights,
    size_t n,
    double ll)
{
    if (!(params.alpha_prior > 0.0) || std::abs(params.alpha_prior - 1.0) < 1e-15) {
        return ll;
    }
    const double coeff = params.alpha_prior - 1.0;
    double prior_term = 0.0;
    for (size_t i = 0; i < n; ++i) {
        prior_term += std::log(std::max(weights[i], params.min_weight));
    }
    return ll + coeff * prior_term;
}

// SQUAREM helpers (Varadhan & Roland 2008)
// q1 = x1 - x, q2 = x2 - x1, v = q2 - q1.
inline void squarem_build_directions(
    const double* x,
    const double* x1,
    const double* x2,
    double* q1,
    double* v,
    size_t n,
    double& sr2,
    double& sv2,
    double& srv)
{
    sr2 = 0.0;
    sv2 = 0.0;
    srv = 0.0;

    #pragma omp simd reduction(+:sr2, sv2, srv)
    for (size_t i = 0; i < n; ++i) {
        q1[i] = x1[i] - x[i];
        v[i] = (x2[i] - x1[i]) - q1[i];
        sr2 += q1[i] * q1[i];
        sv2 += v[i] * v[i];
        srv += q1[i] * v[i];
    }
}

inline double squarem_step_length(int method, double sr2, double sv2, double srv) {
    constexpr double kEps = 1e-15;
    double alpha = 1.0;
    if (method == 1) {
        if (sv2 > kEps) alpha = -srv / sv2;
    } else if (method == 2) {
        if (std::abs(srv) > kEps) alpha = -sr2 / srv;
    } else {
        if (sv2 > kEps) alpha = std::sqrt(sr2 / sv2);
    }
    return std::isfinite(alpha) ? alpha : 1.0;
}

inline void squarem_extrapolate(
    double* x,
    const double* q1,
    const double* v,
    size_t n,
    double alpha,
    double min_weight)
{
    const double alpha_sq = alpha * alpha;
    const double two_alpha = 2.0 * alpha;

    #pragma omp simd
    for (size_t i = 0; i < n; ++i) {
        x[i] = x[i] + two_alpha * q1[i] + alpha_sq * v[i];
        x[i] = std::max(x[i], min_weight);
    }

    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) sum += x[i];
    const double inv_sum = 1.0 / sum;
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
            log_w += length_log_penalty(data, params, aln.ref_idx);
            double score_contrib = params.lambda_b * static_cast<double>(aln.bit_score);

            if (params.use_alignment_damage_likelihood) {
                // Principled model:
                // score_A ~ log w_t + score + log(p_read) + log P(evidence_rt | A)
                // score_M ~ log w_t + score + log(1-p_read) + log P(evidence_rt | M)
                const double p_read = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                log_scores[2 * j]     = log_w + score_contrib + std::log(p_read) + aln.damage_ll_a;
                log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - p_read) + aln.damage_ll_m;
            } else {
                // Legacy model: simple noisy-indicator interpretation of damage_score.
                double ds = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                double log_p_ds_ancient = std::log(ds);
                double log_p_ds_modern = std::log(1.0 - ds);
                log_scores[2 * j]     = log_w + score_contrib + std::log(pi_a) + log_p_ds_ancient;
                log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - pi_a) + log_p_ds_modern;
            }
        }

        double lse = log_sum_exp(log_scores.data(), 2 * deg);
        ll += lse;

        for (uint32_t j = 0; j < deg; ++j) {
            double g_ancient = 0.0, g_modern = 0.0;
            double diff_a = log_scores[2 * j] - lse;
            double diff_m = log_scores[2 * j + 1] - lse;
            g_ancient = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
            g_modern  = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;

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
    state.gamma.resize(A, 0.0);
    if (params.use_damage) {
        state.gamma_ancient.resize(A, 0.0);
    }
    if (T == 0) {
        state.log_likelihood = 0.0;
        state.iterations = 0;
        return state;
    }
    state.weights.assign(T, 1.0 / static_cast<double>(T));  // uniform init

    // Scratch buffers for SQUAREM
    std::vector<double> w0, w1, w2, r_buf, v_buf;
    if (params.use_squarem) {
        w0.resize(T);
        w1.resize(T);
        w2.resize(T);
        r_buf.resize(T);
        v_buf.resize(T);
    }

    double prev_obj = -std::numeric_limits<double>::infinity();
    double step_min = params.squarem_step_min0;
    double step_max0 = params.squarem_step_max0;
    if (step_max0 < step_min) std::swap(step_max0, step_min);
    double step_max = step_max0;
    const double mstep = std::max(params.squarem_mstep, 1.0000001);

    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        double ll;
        double objective;

        const bool squarem_active = params.use_squarem && (iter >= params.squarem_warmup_iters);
        if (squarem_active) {
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

            // SQUAREM extrapolation with textbook step rules.
            double sr2 = 0.0, sv2 = 0.0, srv = 0.0;
            squarem_build_directions(
                w0.data(), w1.data(), w2.data(), r_buf.data(), v_buf.data(), T, sr2, sv2, srv);
            double alpha = squarem_step_length(params.squarem_method, sr2, sv2, srv);
            alpha = std::clamp(alpha, step_min, step_max);
            squarem_extrapolate(w0.data(), r_buf.data(), v_buf.data(), T, alpha, params.min_weight);

            // Stabilization step from canonical SQUAREM (Varadhan & Roland 2008).
            // Always applied after extrapolation to maintain monotonicity guarantee.
            {
                if (params.use_damage) {
                    e_step_damage(data, params, w0.data(), state.pi, state.gamma.data(), state.gamma_ancient.data());
                } else {
                    e_step(data, params, w0.data(), state.gamma.data());
                }
                m_step(data, params, state.gamma.data(), w1.data());
                std::copy(w1.begin(), w1.end(), w0.begin());
            }

            // Evaluate at accelerated point
            if (params.use_damage) {
                ll = e_step_damage(data, params, w0.data(), state.pi, state.gamma.data(), state.gamma_ancient.data());
            } else {
                ll = e_step(data, params, w0.data(), state.gamma.data());
            }
            objective = em_objective_from_ll(params, w0.data(), T, ll);

            // Objective safeguard: revert to EM two-step point when acceleration is non-monotone.
            const bool non_monotone =
                (iter > 0) &&
                std::isfinite(params.squarem_objfn_inc) &&
                (objective < (prev_obj - params.squarem_objfn_inc));
            if (non_monotone) {
                std::copy(w2.begin(), w2.end(), state.weights.begin());
                if (params.use_damage) {
                    ll = e_step_damage(data, params, state.weights.data(), state.pi, state.gamma.data(), state.gamma_ancient.data());
                } else {
                    ll = e_step(data, params, state.weights.data(), state.gamma.data());
                }
                objective = em_objective_from_ll(params, state.weights.data(), T, ll);

                const bool at_upper =
                    std::abs(alpha - step_max) <= (1e-12 * (1.0 + std::abs(step_max)));
                if (at_upper) {
                    step_max = std::max(step_max0, step_max / mstep);
                }
            } else {
                std::copy(w0.begin(), w0.end(), state.weights.begin());
            }

            const bool at_upper =
                std::abs(alpha - step_max) <= (1e-12 * (1.0 + std::abs(step_max)));
            if (at_upper) {
                step_max = mstep * step_max;
            }
            const bool at_lower =
                std::abs(alpha - step_min) <= (1e-12 * (1.0 + std::abs(step_min)));
            if (step_min < 0.0 && at_lower) {
                step_min = mstep * step_min;
            }
        } else {
            // Plain EM
            if (params.use_damage) {
                ll = e_step_damage(data, params, state.weights.data(), state.pi, state.gamma.data(), state.gamma_ancient.data());
            } else {
                ll = e_step(data, params, state.weights.data(), state.gamma.data());
            }
            objective = em_objective_from_ll(params, state.weights.data(), T, ll);
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
            double rel_change =
                std::abs(objective - prev_obj) / (std::abs(prev_obj) + 1e-15);
            if (rel_change < params.tol && (iter + 1) >= params.min_iters) break;
        }
        prev_obj = objective;
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
    ReferenceStatsCollector* stats_collector = nullptr,
    const EMProgressCallback& progress_cb = EMProgressCallback());

std::vector<std::pair<uint32_t, float>> reassign_reads(
    const AlignmentData& data,
    const EMState& state,
    double gamma_threshold = 0.01);

// Forward declaration for columnar index
class ColumnarIndexReader;

// Streaming/chunked EM result - doesn't store per-alignment gamma in memory
struct StreamingEMResult {
    std::vector<double> weights;       // Per-reference abundance (T elements)
    double pi = 0.1;                   // Ancient fraction
    double log_likelihood = 0.0;
    uint32_t iterations = 0;
    uint32_t num_reads = 0;
    uint32_t num_refs = 0;
    uint64_t num_alignments = 0;
};

// Memory-efficient streaming EM that processes row groups one at a time.
// Memory: O(num_refs) + O(row_group_size) instead of O(num_alignments).
//
// For 584M alignments with 10M refs: ~80MB instead of ~18GB.
//
// The algorithm:
// 1. During iteration: E-step per row group, accumulate ref_sums, discard gamma
// 2. After convergence: weights are final, gamma computed on-demand per row group
StreamingEMResult streaming_em(
    ColumnarIndexReader& reader,
    const EMParams& params,
    const EMProgressCallback& progress_cb = EMProgressCallback());

// Callback for processing final gamma values per row group (after convergence)
// Called once per row group with the computed responsibilities
using StreamingGammaCallback = std::function<void(
    uint32_t row_group_idx,
    uint32_t num_rows,
    const uint32_t* read_idx,
    const uint32_t* ref_idx,
    const float* gamma,           // Per-alignment responsibility
    const float* gamma_ancient    // Per-alignment ancient responsibility (may be null)
)>;

// Compute final gamma values in streaming fashion after EM convergence.
// Calls callback for each row group with computed responsibilities.
void streaming_em_finalize(
    ColumnarIndexReader& reader,
    const StreamingEMResult& result,
    const EMParams& params,
    StreamingGammaCallback callback);

} // namespace dart
