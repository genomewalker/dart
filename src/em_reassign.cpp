// EM multi-mapping resolution with SQUAREM acceleration
//
// Resolves ambiguous protein-to-reference alignments using iteratively-refined
// reference abundance estimates. Multi-mappers (reads hitting >1 reference) get
// fractional weights; unique mappers anchor the abundance profile.
//
// OpenMP parallelizes the E-step across reads; M-step uses thread-local
// accumulators merged via SIMD reduction.

#include "agp/em_reassign.hpp"
#include "agp/mmap_array.hpp"
#include "agp/reference_stats.hpp"
#include "agp/damage_index_reader.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace agp {

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
        a.fident = static_cast<float>(ca.identity_q) / 65535.0f;
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

    #pragma omp parallel reduction(+:ll)
    {
        std::vector<double> log_scores;

        #pragma omp for schedule(dynamic, 256)
        for (uint32_t r = 0; r < data.num_reads; ++r) {
            const uint32_t start = data.read_offsets[r];
            const uint32_t end = data.read_offsets[r + 1];
            const uint32_t deg = end - start;
            if (deg == 0) continue;

            // Unique mappers: skip computation
            if (deg == 1) {
                gamma[start] = 1.0;
                ll += std::log(std::max(weights[data.alignments[start].ref_idx], params.min_weight))
                    + params.lambda_b * static_cast<double>(data.alignments[start].bit_score);
                continue;
            }

            log_scores.resize(deg);

            for (uint32_t j = 0; j < deg; ++j) {
                const auto& aln = data.alignments[start + j];
                log_scores[j] = std::log(std::max(weights[aln.ref_idx], params.min_weight))
                               + params.lambda_b * static_cast<double>(aln.bit_score);
            }

            double lse = log_sum_exp(log_scores.data(), deg);
            ll += lse;

            for (uint32_t j = 0; j < deg; ++j) {
                double diff = log_scores[j] - lse;
                gamma[start + j] = (diff > -100.0) ? std::exp(diff) : 0.0;
            }
        }
    }

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

    #pragma omp parallel reduction(+:ll)
    {
        std::vector<double> log_scores;

        #pragma omp for schedule(dynamic, 256)
        for (uint32_t r = 0; r < data.num_reads; ++r) {
            const uint32_t start = data.read_offsets[r];
            const uint32_t end = data.read_offsets[r + 1];
            const uint32_t deg = end - start;
            if (deg == 0) continue;

            log_scores.resize(2 * deg);

            for (uint32_t j = 0; j < deg; ++j) {
                const auto& aln = data.alignments[start + j];
                double log_w = std::log(std::max(weights[aln.ref_idx], params.min_weight));
                double score_contrib = params.lambda_b * static_cast<double>(aln.bit_score);

                double ds = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                double log_p_ds_ancient = std::log(ds);
                double log_p_ds_modern = std::log(1.0 - ds);

                log_scores[2 * j]     = log_w + score_contrib + log_pi_a + log_p_ds_ancient;
                log_scores[2 * j + 1] = log_w + score_contrib + log_pi_m + log_p_ds_modern;
            }

            double lse = log_sum_exp(log_scores.data(), 2 * deg);
            ll += lse;

            for (uint32_t j = 0; j < deg; ++j) {
                double diff_a = log_scores[2 * j] - lse;
                double diff_m = log_scores[2 * j + 1] - lse;
                double g_a = (diff_a > -100.0) ? std::exp(diff_a) : 0.0;
                double g_m = (diff_m > -100.0) ? std::exp(diff_m) : 0.0;

                gamma[start + j] = g_a + g_m;
                gamma_ancient[start + j] = g_a;
            }
        }
    }

    return ll;
}

// ---------------------------------------------------------------------------
// Parallel M-step with thread-local accumulators
// ---------------------------------------------------------------------------

static void m_step_parallel(
    const AlignmentData& data,
    const double* gamma,
    double min_weight,
    double* out_weights)
{
    const uint32_t T = data.num_refs;
    const double N = static_cast<double>(data.num_reads);

    std::memset(out_weights, 0, T * sizeof(double));

    #pragma omp parallel
    {
        std::vector<double> local(T, 0.0);

        #pragma omp for schedule(static)
        for (size_t i = 0; i < data.alignments.size(); ++i) {
            local[data.alignments[i].ref_idx] += gamma[i];
        }

        #pragma omp critical
        {
            for (uint32_t t = 0; t < T; ++t) {
                out_weights[t] += local[t];
            }
        }
    }

    // Normalize
    double sum = 0.0;
    for (uint32_t t = 0; t < T; ++t) {
        out_weights[t] = std::max(out_weights[t] / N, min_weight);
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
    ReferenceStatsCollector* stats_collector)
{
    const uint32_t T = data.num_refs;
    const size_t A = data.alignments.size();

    EMState state;
    state.weights.assign(T, 1.0 / static_cast<double>(T));
    state.gamma.resize(A, 0.0);

    if (params.use_damage) {
        state.gamma_ancient.resize(A, 0.0);
    }

    // SQUAREM scratch buffers
    std::vector<double> w0(T), w1(T), w2(T), r_buf(T), v_buf(T);
    std::vector<double> gamma_scratch(A);
    std::vector<double> gamma_ancient_scratch;
    if (params.use_damage) {
        gamma_ancient_scratch.resize(A, 0.0);
    }

    double prev_ll = -std::numeric_limits<double>::infinity();

    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        double ll;

        if (params.use_squarem) {
            // Save current weights
            std::copy(state.weights.begin(), state.weights.end(), w0.begin());

            // First EM step: w0 -> w1
            if (params.use_damage) {
                e_step_damage_parallel(data, params, w0.data(), state.pi,
                    gamma_scratch.data(), gamma_ancient_scratch.data());
            } else {
                e_step_parallel(data, params, w0.data(), gamma_scratch.data());
            }
            m_step_parallel(data, gamma_scratch.data(), params.min_weight, w1.data());

            // Second EM step: w1 -> w2
            if (params.use_damage) {
                e_step_damage_parallel(data, params, w1.data(), state.pi,
                    gamma_scratch.data(), gamma_ancient_scratch.data());
            } else {
                e_step_parallel(data, params, w1.data(), gamma_scratch.data());
            }
            m_step_parallel(data, gamma_scratch.data(), params.min_weight, w2.data());

            // SQUAREM extrapolation: modifies w0 in-place
            squarem_step(w0.data(), w1.data(), w2.data(), r_buf.data(), v_buf.data(), T);

            // Evaluate at accelerated point
            if (params.use_damage) {
                ll = e_step_damage_parallel(data, params, w0.data(), state.pi,
                    state.gamma.data(), state.gamma_ancient.data());
            } else {
                ll = e_step_parallel(data, params, w0.data(), state.gamma.data());
            }

            // Safeguard: fall back to w2 if acceleration worsened likelihood
            if (ll < prev_ll && iter > 0) {
                std::copy(w2.begin(), w2.end(), state.weights.begin());
                if (params.use_damage) {
                    ll = e_step_damage_parallel(data, params, state.weights.data(), state.pi,
                        state.gamma.data(), state.gamma_ancient.data());
                } else {
                    ll = e_step_parallel(data, params, state.weights.data(), state.gamma.data());
                }
            } else {
                std::copy(w0.begin(), w0.end(), state.weights.begin());
            }
        } else {
            // Plain EM
            if (params.use_damage) {
                ll = e_step_damage_parallel(data, params, state.weights.data(), state.pi,
                    state.gamma.data(), state.gamma_ancient.data());
            } else {
                ll = e_step_parallel(data, params, state.weights.data(), state.gamma.data());
            }
            m_step_parallel(data, state.gamma.data(), params.min_weight, state.weights.data());
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

        // Convergence check: relative log-likelihood change
        if (iter > 0) {
            double rel_change = std::abs(ll - prev_ll) / (std::abs(prev_ll) + 1e-15);
            if (rel_change < params.tol) break;
        }
        prev_ll = ll;
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

        if (best_gamma >= gamma_threshold) {
            result[r] = {best_ref, static_cast<float>(best_gamma)};
        }
    }

    return result;
}

} // namespace agp
