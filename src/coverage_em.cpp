// Coverage-aware EM outer loop
//
// Alternates between inner EM (streaming_em) and per-protein coverage
// uniformity reweighting. Proteins with patchy coverage get down-weighted,
// driving the inner EM toward references with uniform read distribution.

#include "dart/coverage_em.hpp"
#include "dart/em_reassign.hpp"
#include "dart/columnar_index.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dart {

namespace {

// Cache 0.5*log(len) for uint16 lengths to avoid repeated logs in hot loops.
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

inline uint32_t accumulation_bin_count(const CoverageEMParams& p) {
    if (!p.adaptive_bins) return p.default_bins;
    return std::max({p.min_bins, p.default_bins, p.max_bins});
}

// Re-bin counts from uniform source bins to uniform destination bins by overlap.
inline void rebin_uniform_counts(
    const std::vector<float>& in,
    uint32_t out_bins,
    std::vector<float>& out)
{
    const uint32_t in_bins = static_cast<uint32_t>(in.size());
    out.assign(out_bins, 0.0f);
    if (in_bins == 0 || out_bins == 0) return;
    if (in_bins == out_bins) {
        out = in;
        return;
    }

    for (uint32_t i = 0; i < in_bins; ++i) {
        const double in_left = static_cast<double>(i) / static_cast<double>(in_bins);
        const double in_right = static_cast<double>(i + 1) / static_cast<double>(in_bins);
        for (uint32_t o = 0; o < out_bins; ++o) {
            const double out_left = static_cast<double>(o) / static_cast<double>(out_bins);
            const double out_right = static_cast<double>(o + 1) / static_cast<double>(out_bins);
            const double overlap = std::max(0.0, std::min(in_right, out_right) - std::max(in_left, out_left));
            if (overlap <= 0.0) continue;
            // Counts in each source bin are treated as uniform over bin width.
            const double frac_of_in_bin = overlap * static_cast<double>(in_bins);
            out[o] += static_cast<float>(static_cast<double>(in[i]) * frac_of_in_bin);
        }
    }
}

inline void compute_coverage_metrics_for_protein(
    ProteinCoverageStats& cs,
    const CoverageEMParams& cov_params,
    uint16_t tlen_for_completeness)
{
    if (cs.sum_gamma_sq > 0.0) {
        cs.n_eff = static_cast<float>(cs.sum_gamma * cs.sum_gamma / cs.sum_gamma_sq);
    }

    const uint32_t B_acc = static_cast<uint32_t>(cs.x_bins.size());
    if (B_acc == 0 || cs.q_bins.size() != cs.x_bins.size()) {
        cs.kl_divergence = 0.0;
        cs.coverage_deviance = 0.0f;
        cs.coverage_weight = 1.0f;
        cs.completeness_fhat = 1.0f;
        return;
    }

    uint32_t B_eval = cov_params.adaptive_bins
        ? adaptive_bin_count(cs.n_eff, cov_params)
        : cov_params.default_bins;
    if (B_eval == 0) B_eval = B_acc;
    B_eval = std::min(B_eval, B_acc);

    const std::vector<float>* x_bins_eval = &cs.x_bins;
    const std::vector<float>* q_bins_eval = &cs.q_bins;
    std::vector<float> x_rebinned;
    std::vector<float> q_rebinned;
    if (B_eval != B_acc) {
        rebin_uniform_counts(cs.x_bins, B_eval, x_rebinned);
        rebin_uniform_counts(cs.q_bins, B_eval, q_rebinned);
        x_bins_eval = &x_rebinned;
        q_bins_eval = &q_rebinned;
    }

    cs.kl_divergence = compute_kl_divergence(*x_bins_eval, *q_bins_eval, B_eval);
    cs.coverage_deviance = 2.0f * cs.n_eff * static_cast<float>(cs.kl_divergence);
    cs.coverage_weight = compute_coverage_weight(cs.n_eff, cs.kl_divergence, B_eval, cov_params);

    const float avg_aln_len = (cs.sum_gamma > 0.0)
        ? static_cast<float>(cs.sum_aln_len / cs.sum_gamma)
        : 0.0f;
    cs.completeness_fhat = estimate_completeness(
        *x_bins_eval, cs.n_eff, avg_aln_len, tlen_for_completeness, B_eval);
}

} // namespace

// ---------------------------------------------------------------------------
// coverage_accumulate_pass: streaming pass to build per-protein stats
// ---------------------------------------------------------------------------

void coverage_accumulate_pass(
    ColumnarIndexReader& reader,
    const StreamingEMResult& em_result,
    const EMParams& em_params,
    const CoverageEMParams& cov_params,
    std::unordered_map<uint32_t, ProteinCoverageStats>& stats_map)
{
    const double eps = 1e-12;
    const double pi_a = std::clamp(em_result.pi, eps, 1.0 - eps);
    const uint32_t B_acc = accumulation_bin_count(cov_params);
    if (B_acc == 0) return;
    if (em_result.weights.empty()) return;

    std::vector<double> log_weights(em_result.weights.size());
    for (size_t t = 0; t < em_result.weights.size(); ++t) {
        log_weights[t] = std::log(std::max(em_result.weights[t], em_params.min_weight));
    }

    // Determine thread count for thread-local accumulators
    int num_threads = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        num_threads = omp_get_num_threads();
    }
#endif
    num_threads = std::min(num_threads, 256);

    // Thread-local stats maps
    std::vector<std::unordered_map<uint32_t, ProteinCoverageStats>> thread_stats(num_threads);

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
        const uint16_t* aln_len,
        const uint16_t* /*qstart*/,
        const uint16_t* /*qend*/,
        const uint16_t* tstart,
        const uint16_t* /*tend*/,
        const uint16_t* tlen,
        const uint16_t* /*qlen*/,
        const uint16_t* /*mismatch*/,
        const uint16_t* /*gapopen*/
    ) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
        if (tid >= num_threads) tid = 0;
#else
        const int tid = 0;
#endif
        auto& local_stats = thread_stats[tid];

        // Recompute gamma using converged weights (mirrors streaming_em_finalize)
        thread_local std::vector<double> log_scores;
        uint32_t row = 0;
        while (row < num_rows) {
            const uint32_t current_read = read_idx[row];
            uint32_t read_start = row;
            while (row < num_rows && read_idx[row] == current_read) {
                ++row;
            }
            uint32_t deg = row - read_start;

            // Compute gamma for this read's alignments
            thread_local std::vector<double> gamma_vals;
            gamma_vals.resize(deg);

            if (em_params.use_damage) {
                const bool use_aln_dmg = em_params.use_alignment_damage_likelihood && dmg_ll_a && dmg_ll_m;
                log_scores.resize(2 * deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double log_w = log_weights[ref_idx[idx]];
                    if (em_params.normalize_by_length && tlen && tlen[idx] > 0) {
                        log_w -= half_log_tlen(tlen[idx]);
                    }
                    double score_contrib = em_params.lambda_b * static_cast<double>(bit_score[idx]);

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
                    double diff_a = log_scores[2 * j] - lse;
                    double diff_m = log_scores[2 * j + 1] - lse;
                    double g_a = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
                    double g_m = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;
                    gamma_vals[j] = g_a + g_m;
                }
            } else {
                log_scores.resize(deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double log_w = log_weights[ref_idx[idx]];
                    if (em_params.normalize_by_length && tlen && tlen[idx] > 0) {
                        log_w -= half_log_tlen(tlen[idx]);
                    }
                    log_scores[j] = log_w + em_params.lambda_b * static_cast<double>(bit_score[idx]);
                }

                double lse = log_sum_exp(log_scores.data(), deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    double diff = log_scores[j] - lse;
                    gamma_vals[j] = (diff > -700.0) ? std::exp(diff) : 0.0;
                }
            }

            // Accumulate coverage stats per protein
            for (uint32_t j = 0; j < deg; ++j) {
                double g = gamma_vals[j];
                if (g < 1e-6) continue;

                uint32_t idx = read_start + j;
                uint32_t ri = ref_idx[idx];

                auto& cs = local_stats[ri];
                if (cs.x_bins.empty()) {
                    cs.ref_idx = ri;
                    cs.x_bins.assign(B_acc, 0.0f);
                    cs.q_bins.assign(B_acc, 0.0f);
                }

                // tstart in EMI is 1-based; subtract 1 for 0-based
                uint16_t tstart_0 = (tstart && tstart[idx] > 0) ? (tstart[idx] - 1) : 0;
                uint16_t target_len = tlen ? tlen[idx] : 0;
                uint16_t alignment_len = aln_len ? aln_len[idx] : 0;

                if (target_len > 0) {
                    cs.tlen = target_len;
                }

                uint32_t bin = start_to_bin(tstart_0, target_len, B_acc);
                cs.x_bins[bin] += static_cast<float>(g);

                compute_q_bin_contribution(cs.q_bins.data(), B_acc,
                    target_len, alignment_len, static_cast<float>(g));

                cs.sum_gamma += g;
                cs.sum_gamma_sq += g * g;
                cs.sum_aln_len += g * static_cast<double>(alignment_len);
            }
        }
    });

    // Merge thread-local stats into output map
    for (int t = 0; t < num_threads; ++t) {
        for (auto& [ri, tcs] : thread_stats[t]) {
            auto& cs = stats_map[ri];
            if (cs.x_bins.empty()) {
                cs = std::move(tcs);
            } else {
                for (uint32_t b = 0; b < B_acc && b < tcs.x_bins.size(); ++b) {
                    cs.x_bins[b] += tcs.x_bins[b];
                    cs.q_bins[b] += tcs.q_bins[b];
                }
                cs.sum_gamma += tcs.sum_gamma;
                cs.sum_gamma_sq += tcs.sum_gamma_sq;
                cs.sum_aln_len += tcs.sum_aln_len;
                if (tcs.tlen > 0) cs.tlen = tcs.tlen;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// coverage_em_outer_loop: alternate inner EM + coverage reweighting
// ---------------------------------------------------------------------------

StreamingEMResult coverage_em_outer_loop(
    ColumnarIndexReader& reader,
    EMParams inner_params,
    const CoverageEMParams& cov_params,
    std::unordered_map<uint32_t, ProteinCoverageStats>& out_stats,
    bool verbose)
{
    const uint32_t T = reader.num_refs();
    StreamingEMResult result;
    bool have_updated_weights = false;

    auto recompute_coverage_stats = [&](const StreamingEMResult& em_for_stats) {
        out_stats.clear();
        coverage_accumulate_pass(reader, em_for_stats, inner_params, cov_params, out_stats);
        for (auto& [ri, cs] : out_stats) {
            (void)ri;
            compute_coverage_metrics_for_protein(cs, cov_params, cs.tlen);
        }
    };

    for (uint32_t outer_iter = 0; outer_iter < cov_params.max_outer_iters; ++outer_iter) {
        if (verbose) {
            std::cerr << "[coverage-em] outer iteration " << (outer_iter + 1)
                      << "/" << cov_params.max_outer_iters << std::endl;
        }

        // 1. Run inner EM
        result = streaming_em(reader, inner_params);

        if (verbose) {
            std::cerr << "[coverage-em]   inner EM: " << result.iterations << " iters"
                      << ", LL=" << result.log_likelihood << std::endl;
        }

        // 2-3. Recompute coverage stats from current EM solution.
        recompute_coverage_stats(result);

        // 4. Compute new weights: blend old EM weights with coverage-reweighted counts
        const std::vector<double>& prev_coverage_weights = inner_params.initial_weights;
        std::vector<double> new_weights(T, 0.0);
        double old_w_sum = 0.0;
        for (uint32_t j = 0; j < T; ++j) {
            double old_w = result.weights[j];
            double n_j = 0.0;
            double w_j = 1.0;

            auto it = out_stats.find(j);
            if (it != out_stats.end()) {
                n_j = it->second.sum_gamma;
                w_j = static_cast<double>(it->second.coverage_weight);
            }

            double raw_w = (n_j + 1e-15) * w_j;
            double log_new = (1.0 - cov_params.eta) * std::log(old_w + 1e-300)
                           + cov_params.eta * std::log(raw_w + 1e-300);
            new_weights[j] = std::exp(log_new);
            old_w_sum += new_weights[j];
        }

        // Normalize
        if (old_w_sum > 0.0) {
            double inv = 1.0 / old_w_sum;
            for (uint32_t j = 0; j < T; ++j) {
                new_weights[j] *= inv;
            }
        }

        // 5. Check outer convergence (compare to previous coverage-adjusted weights).
        bool converged = false;
        if (outer_iter > 0 && prev_coverage_weights.size() == T) {
            double max_rel_change = 0.0;
            for (uint32_t j = 0; j < T; ++j) {
                const double prev_w = prev_coverage_weights[j];
                const double rel = std::abs(new_weights[j] - prev_w) / (prev_w + 1e-15);
                if (rel > max_rel_change) max_rel_change = rel;
            }

            if (verbose) {
                std::cerr << "[coverage-em]   max_rel_change=" << max_rel_change
                          << " (tol=" << cov_params.outer_tol << ")" << std::endl;
            }

            if (max_rel_change < cov_params.outer_tol) {
                if (verbose) {
                    std::cerr << "[coverage-em]   converged after " << (outer_iter + 1)
                              << " outer iterations" << std::endl;
                }
                converged = true;
            }
        } else if (outer_iter > 0 && !prev_coverage_weights.empty() && verbose) {
            std::cerr << "[coverage-em]   warning: initial_weights size mismatch ("
                      << prev_coverage_weights.size() << " vs " << T
                      << "), skipping outer-loop convergence check\n";
        }

        // 6. Update inner_params for next iteration or final refit.
        inner_params.initial_weights = std::move(new_weights);
        have_updated_weights = true;
        if (converged) break;
    }

    // Run one final inner EM at the latest coverage-adjusted weights so returned
    // state is consistent with the final outer-loop update.
    if (have_updated_weights) {
        if (verbose) {
            std::cerr << "[coverage-em] final inner EM refit with coverage-adjusted weights\n";
        }
        result = streaming_em(reader, inner_params);
        recompute_coverage_stats(result);
    }

    return result;
}

// ---------------------------------------------------------------------------
// coverage_accumulate_from_alignments: in-memory coverage accumulation
// ---------------------------------------------------------------------------

void coverage_accumulate_from_alignments(
    const AlignmentData& data,
    const EMState& state,
    const CoverageEMParams& cov_params,
    std::unordered_map<uint32_t, ProteinCoverageStats>& stats_map)
{
    const uint32_t B_acc = accumulation_bin_count(cov_params);
    if (B_acc == 0) return;
    const size_t A = data.alignments.size();

    for (size_t i = 0; i < A; ++i) {
        const double g = state.gamma[i];
        if (g < 1e-6) continue;

        const auto& aln = data.alignments[i];
        auto& cs = stats_map[aln.ref_idx];
        if (cs.x_bins.empty()) {
            cs.ref_idx = aln.ref_idx;
            const uint32_t tlen_u32 = !data.ref_lengths.empty()
                ? data.ref_lengths[aln.ref_idx] : 0u;
            cs.tlen = static_cast<uint16_t>(std::min(tlen_u32, uint32_t(65535)));
            cs.x_bins.assign(B_acc, 0.0f);
            cs.q_bins.assign(B_acc, 0.0f);
        }

        const uint16_t target_len = cs.tlen;
        const uint16_t alignment_len = aln.aln_len;

        const uint32_t bin = start_to_bin(aln.tstart, target_len, B_acc);
        cs.x_bins[bin] += static_cast<float>(g);

        compute_q_bin_contribution(cs.q_bins.data(), B_acc,
            target_len, alignment_len, static_cast<float>(g));

        cs.sum_gamma += g;
        cs.sum_gamma_sq += g * g;
        cs.sum_aln_len += g * static_cast<double>(alignment_len);
    }
}

// ---------------------------------------------------------------------------
// coverage_em_squarem: outer loop using SQUAREM for the columnar path
// ---------------------------------------------------------------------------

EMState coverage_em_squarem(
    const AlignmentData& data,
    EMParams inner_params,
    const CoverageEMParams& cov_params,
    std::unordered_map<uint32_t, ProteinCoverageStats>& out_stats,
    const EMProgressCallback& progress_cb,
    bool verbose)
{
    const uint32_t T = data.num_refs;
    EMState em_state;
    bool have_updated_weights = false;

    auto recompute_coverage_stats = [&](const EMState& state_for_stats) {
        out_stats.clear();
        coverage_accumulate_from_alignments(data, state_for_stats, cov_params, out_stats);
        for (auto& [ri, cs] : out_stats) {
            const uint32_t tlen_u32 = (!data.ref_lengths.empty() && ri < data.ref_lengths.size())
                ? data.ref_lengths[ri]
                : uint32_t(cs.tlen);
            const uint16_t tlen_for_completeness =
                static_cast<uint16_t>(std::min(tlen_u32, uint32_t(65535)));
            compute_coverage_metrics_for_protein(cs, cov_params, tlen_for_completeness);
        }
    };

    for (uint32_t outer_iter = 0; outer_iter < cov_params.max_outer_iters; ++outer_iter) {
        if (verbose) {
            std::cerr << "[coverage-em] outer iteration " << (outer_iter + 1)
                      << "/" << cov_params.max_outer_iters << std::endl;
        }

        // 1. Run inner SQUAREM EM (warm-started via inner_params.initial_weights)
        em_state = squarem_em(data, inner_params, nullptr, progress_cb);

        if (verbose) {
            std::cerr << "[coverage-em]   inner EM: " << em_state.iterations << " iters"
                      << ", LL=" << em_state.log_likelihood << std::endl;
        }

        // 2-3. Recompute coverage stats from current EM solution.
        recompute_coverage_stats(em_state);

        // 4. Blend old EM weights with coverage-reweighted counts.
        // prev_coverage_weights holds the coverage-adjusted weights from the previous
        // outer iteration (stored in inner_params.initial_weights before SQUAREM ran).
        // We compare against these for convergence, not the unconstrained EM weights.
        const std::vector<double>& prev_coverage_weights = inner_params.initial_weights;

        std::vector<double> new_weights(T, 0.0);
        double new_w_sum = 0.0;
        for (uint32_t j = 0; j < T; ++j) {
            const double old_w = em_state.weights[j];
            double n_j = 0.0;
            double w_j = 1.0;
            auto it = out_stats.find(j);
            if (it != out_stats.end()) {
                n_j = it->second.sum_gamma;
                w_j = static_cast<double>(it->second.coverage_weight);
            }
            const double raw_w = (n_j + 1e-15) * w_j;
            const double log_new = (1.0 - cov_params.eta) * std::log(old_w + 1e-300)
                                 + cov_params.eta * std::log(raw_w + 1e-300);
            new_weights[j] = std::exp(log_new);
            new_w_sum += new_weights[j];
        }
        if (new_w_sum > 0.0) {
            const double inv = 1.0 / new_w_sum;
            for (uint32_t j = 0; j < T; ++j) new_weights[j] *= inv;
        }

        // 5. Check outer convergence: compare new coverage weights vs previous coverage weights.
        bool converged = false;
        if (outer_iter > 0 && prev_coverage_weights.size() == T) {
            double max_rel_change = 0.0;
            for (uint32_t j = 0; j < T; ++j) {
                const double prev_w = prev_coverage_weights[j];
                const double rel = std::abs(new_weights[j] - prev_w) / (prev_w + 1e-15);
                if (rel > max_rel_change) max_rel_change = rel;
            }
            if (verbose) {
                std::cerr << "[coverage-em]   max_rel_change=" << max_rel_change
                          << " (tol=" << cov_params.outer_tol << ")" << std::endl;
            }
            if (max_rel_change < cov_params.outer_tol) {
                if (verbose) {
                    std::cerr << "[coverage-em]   converged after " << (outer_iter + 1)
                              << " outer iterations" << std::endl;
                }
                converged = true;
            }
        } else if (outer_iter > 0 && !prev_coverage_weights.empty() && verbose) {
            std::cerr << "[coverage-em]   warning: initial_weights size mismatch ("
                      << prev_coverage_weights.size() << " vs " << T
                      << "), skipping outer-loop convergence check\n";
        }

        // 6. Warm-start SQUAREM for the next outer iteration
        inner_params.initial_weights = std::move(new_weights);
        have_updated_weights = true;
        if (converged) break;
    }

    // Run one final inner EM at the latest coverage-adjusted weights so returned
    // state is consistent with the final outer-loop update.
    if (have_updated_weights) {
        if (verbose) {
            std::cerr << "[coverage-em] final inner EM refit with coverage-adjusted weights\n";
        }
        em_state = squarem_em(data, inner_params, nullptr, progress_cb);
        recompute_coverage_stats(em_state);
    }

    return em_state;
}

} // namespace dart
