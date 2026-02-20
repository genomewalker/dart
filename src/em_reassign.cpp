// EM multi-mapping resolution with SQUAREM acceleration
//
// Resolves ambiguous protein-to-reference alignments using iteratively-refined
// reference abundance estimates. Multi-mappers (reads hitting >1 reference) get
// fractional weights; unique mappers anchor the abundance profile.
//
// OpenMP parallelizes the E-step across reads; M-step uses thread-local
// accumulators merged via SIMD reduction.

#include "dart/em_reassign.hpp"
#include "dart/mmap_array.hpp"
#include "dart/reference_stats.hpp"
#include "dart/damage_index_reader.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dart {

// ---------------------------------------------------------------------------
// Build CSR-format AlignmentData from CompactAlignment array
// ---------------------------------------------------------------------------

AlignmentData build_alignment_data(
    const CompactAlignment* alns,
    size_t n_alns,
    uint32_t num_reads,
    uint32_t num_refs,
    const DamageIndexReader* agd)
{
    AlignmentData data;
    data.num_reads = num_reads;
    data.num_refs = num_refs;
    data.read_offsets.resize(num_reads + 1, 0);
    data.ref_lengths.assign(num_refs, 0);

    // Count alignments per read
    for (size_t i = 0; i < n_alns; ++i) {
        data.read_offsets[alns[i].read_idx + 1]++;
    }

    // Prefix sum
    for (uint32_t r = 1; r <= num_reads; ++r) {
        data.read_offsets[r] += data.read_offsets[r - 1];
    }

    // Fill alignments in CSR order
    data.alignments.resize(n_alns);
    std::vector<uint32_t> cursor(data.read_offsets.begin(), data.read_offsets.end() - 1);

    for (size_t i = 0; i < n_alns; ++i) {
        const auto& ca = alns[i];
        uint32_t pos = cursor[ca.read_idx]++;

        Alignment& a = data.alignments[pos];
        a.read_idx = ca.read_idx;
        a.ref_idx = ca.ref_idx;
        a.bit_score = ca.bit_score;
        a.damage_score = ca.damage_score;
        a.damage_ll_a = ca.damage_ll_a;
        a.damage_ll_m = ca.damage_ll_m;
        a.fident = static_cast<float>(ca.identity_q) / 65535.0f;
        a.tstart = ca.aln_start;
        const uint32_t aln_end32 = static_cast<uint32_t>(ca.aln_end);
        const uint32_t aln_start32 = static_cast<uint32_t>(ca.aln_start);
        a.aln_len = static_cast<uint16_t>(
            (aln_end32 > aln_start32)
                ? std::min(aln_end32 - aln_start32, uint32_t(65535))
                : 0u);
        const uint32_t tlen_hint = static_cast<uint32_t>(ca.flags);
        const uint32_t end_hint = static_cast<uint32_t>(ca.aln_end);
        const uint32_t ref_len = (tlen_hint > 0) ? tlen_hint : end_hint;
        data.ref_lengths[a.ref_idx] = std::max(data.ref_lengths[a.ref_idx], ref_len);
    }

    data.ref_log_len_half.resize(num_refs, 0.0f);
    for (uint32_t t = 0; t < num_refs; ++t) {
        const double len = static_cast<double>(std::max(data.ref_lengths[t], 1u));
        data.ref_log_len_half[t] = static_cast<float>(0.5 * std::log(len));
    }

    return data;
}

// ---------------------------------------------------------------------------
// Parallel E-step with OpenMP
// ---------------------------------------------------------------------------

static double e_step_parallel(
    const AlignmentData& data,
    const EMParams& params,
    const double* weights,
    double* gamma)
{
    double ll = 0.0;
#ifdef _OPENMP
    const int n_threads = std::max(1, omp_get_max_threads());
    std::vector<double> ll_partial(static_cast<size_t>(n_threads), 0.0);
#pragma omp parallel
    {
        std::vector<double> log_scores;
        double ll_local = 0.0;
        const int tid = omp_get_thread_num();

#pragma omp for schedule(static, 256)
        for (uint32_t r = 0; r < data.num_reads; ++r) {
            const uint32_t start = data.read_offsets[r];
            const uint32_t end = data.read_offsets[r + 1];
            const uint32_t deg = end - start;
            if (deg == 0) continue;

            // Unique mappers: skip computation
            if (deg == 1) {
                gamma[start] = 1.0;
                ll_local += std::log(std::max(weights[data.alignments[start].ref_idx], params.min_weight))
                    + length_log_penalty(data, params, data.alignments[start].ref_idx)
                    + params.lambda_b * static_cast<double>(data.alignments[start].bit_score);
                continue;
            }

            log_scores.resize(deg);

            for (uint32_t j = 0; j < deg; ++j) {
                const auto& aln = data.alignments[start + j];
                log_scores[j] = std::log(std::max(weights[aln.ref_idx], params.min_weight))
                               + length_log_penalty(data, params, aln.ref_idx)
                               + params.lambda_b * static_cast<double>(aln.bit_score);
            }

            double lse = log_sum_exp(log_scores.data(), deg);
            ll_local += lse;

            for (uint32_t j = 0; j < deg; ++j) {
                double diff = log_scores[j] - lse;
                gamma[start + j] = (diff > -700.0) ? std::exp(diff) : 0.0;
            }
        }

        ll_partial[static_cast<size_t>(tid)] = ll_local;
    }
    for (double v : ll_partial) ll += v;
#else
    std::vector<double> log_scores;
    for (uint32_t r = 0; r < data.num_reads; ++r) {
        const uint32_t start = data.read_offsets[r];
        const uint32_t end = data.read_offsets[r + 1];
        const uint32_t deg = end - start;
        if (deg == 0) continue;

        if (deg == 1) {
            gamma[start] = 1.0;
            ll += std::log(std::max(weights[data.alignments[start].ref_idx], params.min_weight))
                + length_log_penalty(data, params, data.alignments[start].ref_idx)
                + params.lambda_b * static_cast<double>(data.alignments[start].bit_score);
            continue;
        }

        log_scores.resize(deg);
        for (uint32_t j = 0; j < deg; ++j) {
            const auto& aln = data.alignments[start + j];
            log_scores[j] = std::log(std::max(weights[aln.ref_idx], params.min_weight))
                           + length_log_penalty(data, params, aln.ref_idx)
                           + params.lambda_b * static_cast<double>(aln.bit_score);
        }
        double lse = log_sum_exp(log_scores.data(), deg);
        ll += lse;
        for (uint32_t j = 0; j < deg; ++j) {
            double diff = log_scores[j] - lse;
            gamma[start + j] = (diff > -700.0) ? std::exp(diff) : 0.0;
        }
    }
#endif

    return ll;
}

// ---------------------------------------------------------------------------
// Parallel damage-aware E-step
// ---------------------------------------------------------------------------

static double e_step_damage_parallel(
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
    const double log_pi_a = std::log(pi_a);
    const double log_pi_m = std::log(1.0 - pi_a);
#ifdef _OPENMP
    const int n_threads = std::max(1, omp_get_max_threads());
    std::vector<double> ll_partial(static_cast<size_t>(n_threads), 0.0);
#pragma omp parallel
    {
        std::vector<double> log_scores;
        double ll_local = 0.0;
        const int tid = omp_get_thread_num();

#pragma omp for schedule(static, 256)
        for (uint32_t r = 0; r < data.num_reads; ++r) {
            const uint32_t start = data.read_offsets[r];
            const uint32_t end = data.read_offsets[r + 1];
            const uint32_t deg = end - start;
            if (deg == 0) continue;

            log_scores.resize(2 * deg);

            for (uint32_t j = 0; j < deg; ++j) {
                const auto& aln = data.alignments[start + j];
                double log_w = std::log(std::max(weights[aln.ref_idx], params.min_weight));
                log_w += length_log_penalty(data, params, aln.ref_idx);
                double score_contrib = params.lambda_b * static_cast<double>(aln.bit_score);

                if (params.use_alignment_damage_likelihood) {
                    const double p_read = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                    log_scores[2 * j]     = log_w + score_contrib + std::log(p_read) + aln.damage_ll_a;
                    log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - p_read) + aln.damage_ll_m;
                } else {
                    double ds = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                    double log_p_ds_ancient = std::log(ds);
                    double log_p_ds_modern = std::log(1.0 - ds);
                    log_scores[2 * j]     = log_w + score_contrib + log_pi_a + log_p_ds_ancient;
                    log_scores[2 * j + 1] = log_w + score_contrib + log_pi_m + log_p_ds_modern;
                }
            }

            double lse = log_sum_exp(log_scores.data(), 2 * deg);
            ll_local += lse;

            for (uint32_t j = 0; j < deg; ++j) {
                double diff_a = log_scores[2 * j] - lse;
                double diff_m = log_scores[2 * j + 1] - lse;
                double g_a = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
                double g_m = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;

                gamma[start + j] = g_a + g_m;
                gamma_ancient[start + j] = g_a;
            }
        }

        ll_partial[static_cast<size_t>(tid)] = ll_local;
    }
    for (double v : ll_partial) ll += v;
#else
    std::vector<double> log_scores;
    for (uint32_t r = 0; r < data.num_reads; ++r) {
        const uint32_t start = data.read_offsets[r];
        const uint32_t end = data.read_offsets[r + 1];
        const uint32_t deg = end - start;
        if (deg == 0) continue;

        log_scores.resize(2 * deg);
        for (uint32_t j = 0; j < deg; ++j) {
            const auto& aln = data.alignments[start + j];
            double log_w = std::log(std::max(weights[aln.ref_idx], params.min_weight));
            log_w += length_log_penalty(data, params, aln.ref_idx);
            double score_contrib = params.lambda_b * static_cast<double>(aln.bit_score);
            if (params.use_alignment_damage_likelihood) {
                const double p_read = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                log_scores[2 * j]     = log_w + score_contrib + std::log(p_read) + aln.damage_ll_a;
                log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - p_read) + aln.damage_ll_m;
            } else {
                double ds = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                double log_p_ds_ancient = std::log(ds);
                double log_p_ds_modern = std::log(1.0 - ds);
                log_scores[2 * j]     = log_w + score_contrib + log_pi_a + log_p_ds_ancient;
                log_scores[2 * j + 1] = log_w + score_contrib + log_pi_m + log_p_ds_modern;
            }
        }

        double lse = log_sum_exp(log_scores.data(), 2 * deg);
        ll += lse;
        for (uint32_t j = 0; j < deg; ++j) {
            double diff_a = log_scores[2 * j] - lse;
            double diff_m = log_scores[2 * j + 1] - lse;
            double g_a = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
            double g_m = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;
            gamma[start + j] = g_a + g_m;
            gamma_ancient[start + j] = g_a;
        }
    }
#endif

    return ll;
}

// ---------------------------------------------------------------------------
// Parallel M-step with thread-local accumulators
// ---------------------------------------------------------------------------

static void m_step_parallel(
    const AlignmentData& data,
    const EMParams& params,
    const double* gamma,
    double* out_weights)
{
    const uint32_t T = data.num_refs;

#ifdef _OPENMP
    const int n_threads = std::max(1, omp_get_max_threads());
#else
    const int n_threads = 1;
#endif
    std::vector<double> partial(static_cast<size_t>(n_threads) * T, 0.0);

#pragma omp parallel
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        double* local = partial.data() + static_cast<size_t>(tid) * T;

#pragma omp for schedule(static)
        for (size_t i = 0; i < data.alignments.size(); ++i) {
            local[data.alignments[i].ref_idx] += gamma[i];
        }
    }

    // MAP M-step with symmetric Dirichlet(alpha_prior).
    double sum = 0.0;
    const double prior_offset = params.alpha_prior - 1.0;
    for (uint32_t t = 0; t < T; ++t) {
        double w = prior_offset;
        for (int tid = 0; tid < n_threads; ++tid) {
            w += partial[static_cast<size_t>(tid) * T + t];
        }

        out_weights[t] = std::max(w, params.min_weight);
        sum += out_weights[t];
    }

    double inv_sum = 1.0 / sum;
#pragma omp simd
    for (uint32_t t = 0; t < T; ++t) {
        out_weights[t] *= inv_sum;
    }
}

// ---------------------------------------------------------------------------
// SQUAREM EM solver (parallel version)
// ---------------------------------------------------------------------------

EMState squarem_em(
    const AlignmentData& data,
    const EMParams& params,
    ReferenceStatsCollector* stats_collector,
    const EMProgressCallback& progress_cb)
{
    if (params.alpha_prior <= 0.0) {
        throw std::invalid_argument("EMParams::alpha_prior must be > 0 (use 1.0 for MLE)");
    }

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
    if (!params.initial_weights.empty() && params.initial_weights.size() == T) {
        state.weights = params.initial_weights;
        double sum = 0.0;
        for (double w : state.weights) sum += w;
        if (sum > 0.0) for (double& w : state.weights) w /= sum;
    } else {
        state.weights.assign(T, 1.0 / static_cast<double>(T));
    }

    // SQUAREM scratch buffers
    std::vector<double> w0(T), w1(T), w2(T), r_buf(T), v_buf(T);
    std::vector<double> prev_weights(T);
    std::vector<double> gamma_scratch(A);
    std::vector<double> gamma_ancient_scratch;
    if (params.use_damage) {
        gamma_ancient_scratch.resize(A, 0.0);
    }

    double prev_obj = -std::numeric_limits<double>::infinity();
    double step_min = params.squarem_step_min0;
    double step_max0 = params.squarem_step_max0;
    if (step_max0 < step_min) std::swap(step_max0, step_min);
    double step_max = step_max0;
    const double mstep = std::max(params.squarem_mstep, 1.0000001);
    const double param_tol = std::sqrt(std::max(params.tol, 1e-15));

    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        double ll;
        double objective;
        bool fallback_to_em = false;
        double alpha_used = 1.0;
        std::copy(state.weights.begin(), state.weights.end(), prev_weights.begin());

        const bool squarem_active = params.use_squarem && (iter >= params.squarem_warmup_iters);
        if (squarem_active) {
            // Save current weights
            std::copy(state.weights.begin(), state.weights.end(), w0.begin());

            // First EM step: w0 -> w1
            if (params.use_damage) {
                e_step_damage_parallel(data, params, w0.data(), state.pi,
                    gamma_scratch.data(), gamma_ancient_scratch.data());
            } else {
                e_step_parallel(data, params, w0.data(), gamma_scratch.data());
            }
            m_step_parallel(data, params, gamma_scratch.data(), w1.data());

            // Second EM step: w1 -> w2
            if (params.use_damage) {
                e_step_damage_parallel(data, params, w1.data(), state.pi,
                    gamma_scratch.data(), gamma_ancient_scratch.data());
            } else {
                e_step_parallel(data, params, w1.data(), gamma_scratch.data());
            }
            m_step_parallel(data, params, gamma_scratch.data(), w2.data());

            // SQUAREM extrapolation with textbook step rules.
            double sr2 = 0.0, sv2 = 0.0, srv = 0.0;
            squarem_build_directions(
                w0.data(), w1.data(), w2.data(), r_buf.data(), v_buf.data(), T, sr2, sv2, srv);
            double alpha = squarem_step_length(params.squarem_method, sr2, sv2, srv);
            alpha = std::clamp(alpha, step_min, step_max);
            alpha_used = alpha;
            squarem_extrapolate(w0.data(), r_buf.data(), v_buf.data(), T, alpha, params.min_weight);

            // Stabilization step from canonical SQUAREM (Varadhan & Roland 2008).
            // Always applied after extrapolation to maintain monotonicity guarantee.
            {
                if (params.use_damage) {
                    e_step_damage_parallel(data, params, w0.data(), state.pi,
                        gamma_scratch.data(), gamma_ancient_scratch.data());
                } else {
                    e_step_parallel(data, params, w0.data(), gamma_scratch.data());
                }
                m_step_parallel(data, params, gamma_scratch.data(), w1.data());
                std::copy(w1.begin(), w1.end(), w0.begin());
            }

            // Evaluate at accelerated point
            if (params.use_damage) {
                ll = e_step_damage_parallel(data, params, w0.data(), state.pi,
                    state.gamma.data(), state.gamma_ancient.data());
            } else {
                ll = e_step_parallel(data, params, w0.data(), state.gamma.data());
            }
            objective = em_objective_from_ll(params, w0.data(), T, ll);

            // Objective safeguard: revert to EM two-step point when acceleration is non-monotone.
            const bool non_monotone =
                (iter > 0) &&
                std::isfinite(params.squarem_objfn_inc) &&
                (objective < (prev_obj - params.squarem_objfn_inc));
            if (non_monotone) {
                fallback_to_em = true;
                std::copy(w2.begin(), w2.end(), state.weights.begin());
                if (params.use_damage) {
                    ll = e_step_damage_parallel(data, params, state.weights.data(), state.pi,
                        state.gamma.data(), state.gamma_ancient.data());
                } else {
                    ll = e_step_parallel(data, params, state.weights.data(), state.gamma.data());
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
                ll = e_step_damage_parallel(data, params, state.weights.data(), state.pi,
                    state.gamma.data(), state.gamma_ancient.data());
            } else {
                ll = e_step_parallel(data, params, state.weights.data(), state.gamma.data());
            }
            objective = em_objective_from_ll(params, state.weights.data(), T, ll);
            m_step_parallel(data, params, state.gamma.data(), state.weights.data());
        }

        // Update global ancient fraction from damage responsibilities
        if (params.use_damage && !state.gamma_ancient.empty()) {
            double sum_ancient = 0.0;
            double sum_total = 0.0;
            #pragma omp parallel for reduction(+:sum_ancient, sum_total) schedule(static)
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

        // Convergence check: relative objective change (LL + prior term when enabled)
        double rel_change = std::numeric_limits<double>::infinity();
        if (iter > 0) {
            rel_change = std::abs(objective - prev_obj) / (std::abs(prev_obj) + 1e-15);
        }
        double param_change_sq = 0.0;
#pragma omp parallel for reduction(+:param_change_sq) schedule(static)
        for (uint32_t t = 0; t < T; ++t) {
            const double d = state.weights[t] - prev_weights[t];
            param_change_sq += d * d;
        }
        const double param_change = std::sqrt(param_change_sq);

        if (progress_cb) {
            EMIterationDiagnostics diag;
            diag.iteration = iter + 1;
            diag.log_likelihood = ll;
            diag.objective = objective;
            diag.rel_change = rel_change;
            diag.param_change = param_change;
            diag.squarem_active = squarem_active;
            diag.fallback_to_em = fallback_to_em;
            diag.alpha = alpha_used;
            diag.step_min = step_min;
            diag.step_max = step_max;
            progress_cb(diag);
        }
        if (iter > 0 &&
            (iter + 1) >= params.min_iters &&
            rel_change <= params.tol &&
            param_change <= param_tol) {
            break;
        }
        prev_obj = objective;
    }

    // Populate reference stats from final responsibilities
    if (stats_collector) {
        for (size_t i = 0; i < A; ++i) {
            const auto& aln = data.alignments[i];
            float gamma_a = params.use_damage
                ? static_cast<float>(state.gamma_ancient[i])
                : 0.0f;
            stats_collector->update(
                std::to_string(aln.ref_idx),
                0,  // ref_length filled by caller
                0, 0,  // aln_start, aln_end filled by caller
                static_cast<float>(state.gamma[i]),
                gamma_a,
                aln.fident);
        }
    }

    return state;
}

// ---------------------------------------------------------------------------
// Produce reassigned output: best or weighted assignments
// ---------------------------------------------------------------------------

std::vector<std::pair<uint32_t, float>> reassign_reads(
    const AlignmentData& data,
    const EMState& state,
    double gamma_threshold)
{
    std::vector<std::pair<uint32_t, float>> result(data.num_reads, {0, 0.0f});

    #pragma omp parallel for schedule(static)
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

        float out_gamma = (best_gamma >= gamma_threshold)
            ? static_cast<float>(best_gamma) : 0.0f;
        result[r] = {best_ref, out_gamma};
    }

    return result;
}

} // namespace dart

// ---------------------------------------------------------------------------
// Streaming/Chunked EM: Memory-efficient implementation
// ---------------------------------------------------------------------------
//
// Instead of storing O(num_alignments) gamma values, we:
// 1. Process row groups one at a time
// 2. Compute E-step locally, accumulate M-step sums to ref_sums
// 3. Discard per-alignment gamma after accumulating
// 4. After convergence, compute gamma on-demand in streaming_em_finalize
//
// Memory: O(num_refs) + O(row_group_size) instead of O(num_alignments)
// For 584M alignments with 10M refs: ~80MB instead of ~18GB

#include "dart/columnar_index.hpp"
#include <atomic>

namespace {

// Cache 0.5*log(len) for uint16 lengths to avoid repeated logs in streaming hot loops.
inline double half_log_tlen(uint16_t len) {
    static const std::array<double, 65536> lut = []() {
        std::array<double, 65536> a{};
        a[0] = 0.0;
        for (size_t i = 1; i < a.size(); ++i) {
            a[i] = 0.5 * std::log(static_cast<double>(i));
        }
        return a;
    }();
    return lut[len];
}

// Per-row-group E-step: compute gamma and accumulate to ref_sums
// Returns partial log-likelihood for this row group
double streaming_e_step_rowgroup(
    uint32_t num_rows,
    const uint32_t* read_idx,
    const uint32_t* ref_idx,
    const float* bit_score,
    const float* damage_score,
    const float* dmg_ll_a,    // Alignment-level log P(evidence | ancient)
    const float* dmg_ll_m,    // Alignment-level log P(evidence | modern)
    const uint16_t* tlen,     // Reference length per alignment (for normalize_by_length)
    const double* log_weights,
    uint32_t num_refs,
    const dart::EMParams& params,
    double pi_ancient,
    std::vector<double>& ref_sums,        // Accumulated gamma per ref (output)
    std::vector<double>& ref_sums_ancient // Accumulated gamma_ancient per ref (output)
) {
    const double eps = 1e-12;
    const double pi_a = std::clamp(pi_ancient, eps, 1.0 - eps);
    double ll = 0.0;

    // Check pointers
    if (!read_idx || !ref_idx || !bit_score || !damage_score || !log_weights) {
        return 0.0;
    }

    // For alignment damage likelihood mode, we need dmg_ll_a/dmg_ll_m
    const bool use_aln_dmg = params.use_alignment_damage_likelihood && dmg_ll_a && dmg_ll_m;

    // Group alignments by read for normalization
    // Since row groups are sorted by read_idx, we can process in order
    uint32_t row = 0;
    uint32_t reads_processed = 0;
    while (row < num_rows) {
        const uint32_t current_read = read_idx[row];
        uint32_t read_start = row;

        // Find all alignments for this read in this row group
        while (row < num_rows && read_idx[row] == current_read) {
            ++row;
        }
        uint32_t deg = row - read_start;
        ++reads_processed;

        if (deg == 0) continue;

        // Compute log-scores for E-step
        thread_local std::vector<double> log_scores;
        if (params.use_damage) {
            log_scores.resize(2 * deg);
            for (uint32_t j = 0; j < deg; ++j) {
                uint32_t idx = read_start + j;
                double log_w = log_weights[ref_idx[idx]];
                if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                    log_w -= half_log_tlen(tlen[idx]);
                }
                double score_contrib = params.lambda_b * static_cast<double>(bit_score[idx]);

                if (use_aln_dmg) {
                    // Use alignment-level damage likelihoods with per-read prior
                    double p_read = std::clamp(static_cast<double>(damage_score[idx]), eps, 1.0 - eps);
                    log_scores[2 * j]     = log_w + score_contrib + std::log(p_read) + dmg_ll_a[idx];
                    log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - p_read) + dmg_ll_m[idx];
                } else {
                    // Simple mode: damage_score as probability, learn pi
                    double ds = std::clamp(static_cast<double>(damage_score[idx]), eps, 1.0 - eps);
                    log_scores[2 * j]     = log_w + score_contrib + std::log(pi_a) + std::log(ds);
                    log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - pi_a) + std::log(1.0 - ds);
                }
            }

            double lse = dart::log_sum_exp(log_scores.data(), 2 * deg);
            ll += lse;

            // Accumulate to ref_sums
            for (uint32_t j = 0; j < deg; ++j) {
                uint32_t idx = read_start + j;
                double diff_a = log_scores[2 * j] - lse;
                double diff_m = log_scores[2 * j + 1] - lse;
                double g_ancient = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
                double g_modern  = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;
                double g_total = g_ancient + g_modern;

                uint32_t ri = ref_idx[idx];
                if (ri >= num_refs) continue;
                ref_sums[ri] += g_total;
                ref_sums_ancient[ri] += g_ancient;
            }
        } else {
            log_scores.resize(deg);
            for (uint32_t j = 0; j < deg; ++j) {
                uint32_t idx = read_start + j;
                double log_w = log_weights[ref_idx[idx]];
                if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                    log_w -= half_log_tlen(tlen[idx]);
                }
                log_scores[j] = log_w + params.lambda_b * static_cast<double>(bit_score[idx]);
            }

            double lse = dart::log_sum_exp(log_scores.data(), deg);
            ll += lse;

            for (uint32_t j = 0; j < deg; ++j) {
                uint32_t idx = read_start + j;
                double diff = log_scores[j] - lse;
                double g_val = (diff > -700.0) ? std::exp(diff) : 0.0;
                ref_sums[ref_idx[idx]] += g_val;
            }
        }
    }

    return ll;
}

// M-step from accumulated ref_sums
void streaming_m_step(
    const std::vector<double>& ref_sums,
    std::vector<double>& weights,
    const dart::EMParams& params
) {
    const uint32_t T = static_cast<uint32_t>(weights.size());
    const double prior_offset = params.alpha_prior - 1.0;

    double sum = 0.0;
    for (uint32_t t = 0; t < T; ++t) {
        weights[t] = std::max(ref_sums[t] + prior_offset, params.min_weight);
        sum += weights[t];
    }

    double inv_sum = 1.0 / sum;
    for (uint32_t t = 0; t < T; ++t) {
        weights[t] *= inv_sum;
    }
}

} // anonymous namespace

namespace dart {

StreamingEMResult streaming_em(
    ColumnarIndexReader& reader,
    const EMParams& params,
    const EMProgressCallback& progress_cb)
{
    if (params.alpha_prior <= 0.0) {
        throw std::invalid_argument("EMParams::alpha_prior must be > 0 (use 1.0 for MLE)");
    }

    StreamingEMResult result;
    result.num_reads = reader.num_reads();
    result.num_refs = reader.num_refs();
    result.num_alignments = reader.num_alignments();

    const uint32_t T = result.num_refs;
    const uint32_t num_row_groups = reader.num_row_groups();

    if (T == 0 || num_row_groups == 0) {
        result.log_likelihood = 0.0;
        result.iterations = 0;
        return result;
    }

    // Initialize weights - use initial_weights if provided, otherwise uniform
    if (!params.initial_weights.empty() && params.initial_weights.size() == T) {
        result.weights = params.initial_weights;
        double sum = 0.0;
        for (double w : result.weights) sum += w;
        if (sum > 0.0) for (double& w : result.weights) w /= sum;
    } else {
        result.weights.assign(T, 1.0 / static_cast<double>(T));
    }
    result.pi = 0.1;

    // Thread-local accumulators for parallel row group processing
    // Use a reasonable upper bound that won't cause OOM
    int num_threads = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        num_threads = omp_get_num_threads();
    }
#endif
    // Cap at hardware concurrency to prevent excessive memory use
    const int hw_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (hw_threads > 0) {
        num_threads = std::min(num_threads, hw_threads);
    }

    std::vector<std::vector<double>> thread_ref_sums(num_threads);
    std::vector<std::vector<double>> thread_ref_sums_ancient(num_threads);
    for (int t = 0; t < num_threads; ++t) {
        thread_ref_sums[t].resize(T, 0.0);
        if (params.use_damage) {
            thread_ref_sums_ancient[t].resize(T, 0.0);
        }
    }

    double prev_ll = -std::numeric_limits<double>::infinity();

    // Pre-allocate buffers outside iteration loop
    std::vector<double> log_weights(T, 0.0);
    std::vector<double> ref_sums(T, 0.0);
    std::vector<double> ref_sums_ancient(T, 0.0);

    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        // Clear accumulators
        for (int t = 0; t < num_threads; ++t) {
            std::fill(thread_ref_sums[t].begin(), thread_ref_sums[t].end(), 0.0);
            if (params.use_damage) {
                std::fill(thread_ref_sums_ancient[t].begin(), thread_ref_sums_ancient[t].end(), 0.0);
            }
        }

        for (uint32_t t = 0; t < T; ++t) {
            log_weights[t] = std::log(std::max(result.weights[t], params.min_weight));
        }

        std::atomic<double> total_ll{0.0};

        // E-step: process row groups in parallel via parallel_scan
        // parallel_scan already creates OMP parallel region internally,
        // so we use thread-local storage inside the callback
        reader.parallel_scan([&](
            uint32_t /*rg_idx*/,
            uint32_t num_rows,
            const uint32_t* read_idx,
            const uint32_t* ref_idx,
            const float* bit_score,
            const float* damage_score,
            const float* /*evalue_log10*/,
            const uint16_t* /*dmg_k*/,
            const uint16_t* /*dmg_m*/,
            const float* dmg_ll_a,
            const float* dmg_ll_m,
            const uint16_t* /*identity_q*/,
            const uint16_t* /*aln_len*/,
            const uint16_t* /*qstart*/,
            const uint16_t* /*qend*/,
            const uint16_t* /*tstart*/,
            const uint16_t* /*tend*/,
            const uint16_t* tlen,
            const uint16_t* /*qlen*/,
            const uint16_t* /*mismatch*/,
            const uint16_t* /*gapopen*/
        ) {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
            if (tid >= num_threads) tid = 0;  // Safety bounds check
#else
            const int tid = 0;
#endif
            auto& local_ref_sums = thread_ref_sums[tid];
            auto& local_ref_sums_ancient = thread_ref_sums_ancient[tid];

            double local_ll = streaming_e_step_rowgroup(
                num_rows, read_idx, ref_idx, bit_score, damage_score,
                dmg_ll_a, dmg_ll_m, tlen,
                log_weights.data(), T, params, result.pi,
                local_ref_sums, local_ref_sums_ancient
            );

            total_ll.fetch_add(local_ll, std::memory_order_relaxed);
        });

        // Merge thread-local accumulators into pre-allocated buffers
        std::fill(ref_sums.begin(), ref_sums.end(), 0.0);
        if (params.use_damage) {
            std::fill(ref_sums_ancient.begin(), ref_sums_ancient.end(), 0.0);
        }
        for (int t = 0; t < num_threads; ++t) {
            for (uint32_t r = 0; r < T; ++r) {
                ref_sums[r] += thread_ref_sums[t][r];
            }
            if (params.use_damage) {
                for (uint32_t r = 0; r < T; ++r) {
                    ref_sums_ancient[r] += thread_ref_sums_ancient[t][r];
                }
            }
        }

        // M-step: update weights from accumulated sums
        streaming_m_step(ref_sums, result.weights, params);

        // Update pi from accumulated ancient sums
        if (params.use_damage) {
            double sum_ancient = 0.0, sum_total = 0.0;
            for (uint32_t r = 0; r < T; ++r) {
                sum_ancient += ref_sums_ancient[r];
                sum_total += ref_sums[r];
            }
            if (sum_total > 0.0) {
                result.pi = std::clamp(sum_ancient / sum_total, 0.01, 0.99);
            }
        }

        result.log_likelihood = total_ll.load();
        result.iterations = iter + 1;

        // Progress callback
        if (progress_cb) {
            EMIterationDiagnostics diag;
            diag.iteration = iter;
            diag.log_likelihood = result.log_likelihood;
            diag.objective = result.log_likelihood;
            diag.rel_change = std::abs(result.log_likelihood - prev_ll) / (std::abs(prev_ll) + 1e-15);
            progress_cb(diag);
        }

        // Convergence check
        if (iter > 0) {
            double rel_change = std::abs(result.log_likelihood - prev_ll) / (std::abs(prev_ll) + 1e-15);
            if (rel_change < params.tol && (iter + 1) >= params.min_iters) {
                break;
            }
        }
        prev_ll = result.log_likelihood;
    }

    return result;
}

void streaming_em_finalize(
    ColumnarIndexReader& reader,
    const StreamingEMResult& result,
    const EMParams& params,
    StreamingGammaCallback callback)
{
    const double eps = 1e-12;
    const double pi_a = std::clamp(result.pi, eps, 1.0 - eps);
    std::vector<double> log_weights(result.weights.size(), 0.0);
    for (size_t t = 0; t < result.weights.size(); ++t) {
        log_weights[t] = std::log(std::max(result.weights[t], params.min_weight));
    }

    // Process each row group and compute final gamma
    reader.parallel_scan([&](
        uint32_t rg_idx,
        uint32_t num_rows,
        const uint32_t* read_idx,
        const uint32_t* ref_idx,
        const float* bit_score,
        const float* damage_score,
        const float* /*evalue_log10*/,
        const uint16_t* /*dmg_k*/,
        const uint16_t* /*dmg_m*/,
        const float* dmg_ll_a,
        const float* dmg_ll_m,
        const uint16_t* /*identity_q*/,
        const uint16_t* /*aln_len*/,
        const uint16_t* /*qstart*/,
        const uint16_t* /*qend*/,
        const uint16_t* /*tstart*/,
        const uint16_t* /*tend*/,
        const uint16_t* tlen,
        const uint16_t* /*qlen*/,
        const uint16_t* /*mismatch*/,
        const uint16_t* /*gapopen*/
    ) {
        // Allocate per-row-group buffers (renamed to avoid conflict with std::tgamma)
        std::vector<float> g_vals(num_rows);
        std::vector<float> g_ancient_vals;
        if (params.use_damage) {
            g_ancient_vals.resize(num_rows);
        }

        // Compute gamma for each alignment
        thread_local std::vector<double> log_scores;
        uint32_t row = 0;
        while (row < num_rows) {
            const uint32_t current_read = read_idx[row];
            uint32_t read_start = row;

            while (row < num_rows && read_idx[row] == current_read) {
                ++row;
            }
            uint32_t deg = row - read_start;

            if (params.use_damage) {
                const bool use_aln_dmg = params.use_alignment_damage_likelihood && dmg_ll_a && dmg_ll_m;
                log_scores.resize(2 * deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double log_w = log_weights[ref_idx[idx]];
                    if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                        log_w -= half_log_tlen(tlen[idx]);
                    }
                    double score_contrib = params.lambda_b * static_cast<double>(bit_score[idx]);

                    if (use_aln_dmg) {
                        double p_read = std::clamp(static_cast<double>(damage_score[idx]), eps, 1.0 - eps);
                        log_scores[2 * j]     = log_w + score_contrib + std::log(p_read) + dmg_ll_a[idx];
                        log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - p_read) + dmg_ll_m[idx];
                    } else {
                        double ds = std::clamp(static_cast<double>(damage_score[idx]), eps, 1.0 - eps);
                        log_scores[2 * j]     = log_w + score_contrib + std::log(pi_a) + std::log(ds);
                        log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - pi_a) + std::log(1.0 - ds);
                    }
                }

                double lse = log_sum_exp(log_scores.data(), 2 * deg);

                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double diff_a = log_scores[2 * j] - lse;
                    double diff_m = log_scores[2 * j + 1] - lse;
                    double g_ancient = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
                    double g_modern  = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;

                    g_vals[idx] = static_cast<float>(g_ancient + g_modern);
                    g_ancient_vals[idx] = static_cast<float>(g_ancient);
                }
            } else {
                log_scores.resize(deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double log_w = log_weights[ref_idx[idx]];
                    if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                        log_w -= half_log_tlen(tlen[idx]);
                    }
                    log_scores[j] = log_w + params.lambda_b * static_cast<double>(bit_score[idx]);
                }

                double lse = log_sum_exp(log_scores.data(), deg);

                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double diff = log_scores[j] - lse;
                    g_vals[idx] = static_cast<float>((diff > -700.0) ? std::exp(diff) : 0.0);
                }
            }
        }

        // Call user callback with computed gamma
        callback(
            rg_idx, num_rows, read_idx, ref_idx,
            g_vals.data(),
            params.use_damage ? g_ancient_vals.data() : nullptr
        );
    });
}

} // namespace dart
