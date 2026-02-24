#pragma once
// Coverage-aware EM: per-protein coverage uniformity weighting
//
// Proteins with highly non-uniform read coverage (reads piling up in one region)
// get down-weighted during EM, since patchy coverage suggests spurious multi-mapping
// rather than true homology across the full protein length.
//
// Core idea: for each protein, bin alignment start positions and compare the observed
// distribution against the expected uniform fragmentation model. Proteins with large
// KL divergence between observed and expected get lower coverage_weight, reducing
// their influence in the EM M-step.
//
// The completeness estimator (fhat) uses Newton-Raphson to solve for the fraction
// of the protein that is actually present in the sample, analogous to capture-recapture.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "dart/em_reassign.hpp"

namespace dart {

// Coverage EM parameters
struct CoverageEMParams {
    uint32_t max_outer_iters = 5;
    uint32_t default_bins = 6;
    uint32_t min_bins = 3;       // adaptive: used for very low effective depth
    uint32_t max_bins = 8;       // adaptive: used for high effective depth
    float tau = 0.13f;
    float depth_gate_lo = 8.0f;
    float depth_gate_hi = 20.0f;
    float w_min = 0.2f;
    float eta = 0.5f;
    float outer_tol = 0.01f;
    bool adaptive_bins = true;
};

// Per-protein coverage statistics accumulated during the outer EM loop
struct ProteinCoverageStats {
    uint32_t ref_idx = 0;
    uint16_t tlen = 0;
    std::vector<float> x_bins;   // observed start count per bin (soft/gamma-weighted)
    std::vector<float> q_bins;   // expected start count under uniform fragmentation
    double sum_gamma = 0.0;
    double sum_gamma_sq = 0.0;
    double sum_aln_len = 0.0;
    double kl_divergence = 0.0;
    float coverage_weight = 1.0f;
    float coverage_deviance = 0.0f;
    float completeness_fhat = 1.0f;
    float n_eff = 0.0f;

    void reset(uint32_t B) {
        x_bins.assign(B, 0.0f);
        q_bins.assign(B, 0.0f);
        sum_gamma = 0.0;
        sum_gamma_sq = 0.0;
        sum_aln_len = 0.0;
        kl_divergence = 0.0;
        coverage_weight = 1.0f;
        coverage_deviance = 0.0f;
        completeness_fhat = 1.0f;
        n_eff = 0.0f;
    }
};

// Choose bin count based on effective depth.
// Low-depth proteins get fewer bins to avoid noisy KL estimates.
inline uint32_t adaptive_bin_count(float n_eff, const CoverageEMParams& p) {
    if (!std::isfinite(n_eff) || n_eff <= 0.0f) return p.default_bins;
    if (n_eff < 8.0f) return p.min_bins;
    if (n_eff < 15.0f) return p.default_bins;
    return p.max_bins;
}

// Map a 0-based alignment start position to a bin index in [0, B-1].
inline uint32_t start_to_bin(uint16_t tstart_0, uint16_t tlen, uint32_t B) {
    if (B == 0 || tlen == 0) return 0;
    uint32_t b = static_cast<uint32_t>(tstart_0) * B / static_cast<uint32_t>(tlen);
    return std::min(b, B - 1);
}

// Accumulate expected start-position contribution under uniform fragmentation.
// A read of alignment length aln_len on a protein of length tlen can start at
// positions [0 .. M-1] where M = max(1, tlen - aln_len + 1). We distribute
// gamma weight uniformly across feasible start positions, binned into B bins.
inline void compute_q_bin_contribution(float* q_bins, uint32_t B,
                                       uint16_t tlen, uint16_t aln_len,
                                       float gamma) {
    if (B == 0 || q_bins == nullptr) return;
    int M = std::max(1, static_cast<int>(tlen) - static_cast<int>(aln_len) + 1);
    float inv_M = 1.0f / static_cast<float>(M);

    for (uint32_t b = 0; b < B; ++b) {
        int bin_start = static_cast<int>(b) * static_cast<int>(tlen) / static_cast<int>(B);
        int bin_end = static_cast<int>(b + 1) * static_cast<int>(tlen) / static_cast<int>(B);
        int feasible = std::max(0, std::min(bin_end, M) - bin_start);
        q_bins[b] += gamma * static_cast<float>(feasible) * inv_M;
    }
}

// KL divergence between observed (x_bins) and expected (q_bins) distributions.
// Jeffreys smoothing: add 0.5/B to each bin before normalizing.
inline double compute_kl_divergence(const std::vector<float>& x_bins,
                                    const std::vector<float>& q_bins,
                                    uint32_t B) {
    if (B == 0 || x_bins.size() < B || q_bins.size() < B) return 0.0;

    double smooth = 0.5 / static_cast<double>(B);
    double sum_x = 0.0, sum_q = 0.0;
    for (uint32_t b = 0; b < B; ++b) {
        sum_x += static_cast<double>(x_bins[b]);
        sum_q += static_cast<double>(q_bins[b]);
    }
    double norm_x = sum_x + 0.5;
    double norm_q = sum_q + 0.5;

    double kl = 0.0;
    for (uint32_t b = 0; b < B; ++b) {
        double p_b = (static_cast<double>(x_bins[b]) + smooth) / norm_x;
        double q_b = (static_cast<double>(q_bins[b]) + smooth) / norm_q;
        kl += p_b * std::log(p_b / q_b);
    }

    if (!std::isfinite(kl) || kl < 0.0) return 0.0;
    return kl;
}

// Depth gate: linear ramp from 0 to 1 between depth_gate_lo and depth_gate_hi.
// Proteins with very few reads (n_eff < lo) get no coverage penalty.
inline float compute_depth_gate(float n_eff, const CoverageEMParams& p) {
    if (n_eff < p.depth_gate_lo) return 0.0f;
    if (n_eff >= p.depth_gate_hi) return 1.0f;
    const float span = p.depth_gate_hi - p.depth_gate_lo;
    if (span <= 0.0f) return (n_eff >= p.depth_gate_lo) ? 1.0f : 0.0f;
    return (n_eff - p.depth_gate_lo) / span;
}

// Coverage weight using null-corrected chi-squared deviance.
// D = 2·n_eff·KL ~ χ²(B-1) under H₀ (uniform coverage).
// Subtract null expectation (B-1) so proteins with truly uniform coverage
// get weight ≈ 1.0; only genuinely patchy proteins are penalized.
inline float compute_coverage_weight(float n_eff, double kl, uint32_t B,
                                     const CoverageEMParams& p) {
    if (B == 0 || !std::isfinite(kl) || kl <= 0.0 || n_eff <= 0.0f) return 1.0f;

    const double g = static_cast<double>(compute_depth_gate(n_eff, p));
    const double deviance = 2.0 * static_cast<double>(n_eff) * kl;
    const double excess = std::max(0.0, deviance - static_cast<double>(B - 1));
    const double exponent = -static_cast<double>(p.tau) * g * excess / 2.0;
    const float w = static_cast<float>(std::exp(exponent));
    return std::clamp(w, p.w_min, 1.0f);
}

// Estimate protein completeness (fraction present in sample) via Newton-Raphson.
// Solves: f * (1 - exp(-n_eff * avg_aln_len / (f * tlen * B))) = B_obs / B
// where B_obs is the number of non-empty bins.
inline float estimate_completeness(const std::vector<float>& x_bins,
                                   float n_eff, float avg_aln_len,
                                   uint16_t tlen, uint32_t B) {
    if (B == 0 || x_bins.size() < B) return 0.0f;

    // A bin counts as "observed" if it has at least 0.5 read-equivalents.
    // This threshold handles soft (gamma-weighted) assignment correctly:
    // hard assignment (gamma=1) needs 1 read; soft (gamma=0.3) needs ~2.
    uint32_t B_obs = 0;
    for (uint32_t b = 0; b < B; ++b) {
        if (x_bins[b] >= 0.5f) ++B_obs;
    }

    if (B_obs == 0 || n_eff < 1.0f || tlen == 0) return 0.0f;
    if (B_obs >= B) return 1.0f;

    float NR_LB = n_eff * avg_aln_len / (static_cast<float>(tlen) * static_cast<float>(B));
    float f = static_cast<float>(B_obs) / static_cast<float>(B);
    float B_obs_f = static_cast<float>(B_obs);
    float B_f = static_cast<float>(B);

    for (int i = 0; i < 10; ++i) {
        float f_safe = std::max(f, 0.01f);
        float e = std::exp(-NR_LB / f_safe);
        float g = f_safe * B_f * (1.0f - e) - B_obs_f;
        float g_prime = B_f * (1.0f - e) + f_safe * B_f * e * NR_LB / (f_safe * f_safe);
        if (std::abs(g_prime) < 1e-10f) break;
        f -= g / g_prime;
        f = std::clamp(f, 0.01f, 1.0f);
    }

    return f;
}

// Compact per-read entry for in-memory coverage accumulation.
// Built from best_hits + em_read_results.gamma after streaming_em_finalize.
// Sorted by ref_idx so coverage_accumulate_from_sorted can process one ref at a time.
struct CovEntry {
    uint32_t ref_idx;
    float    gamma;              // Pass-2 EM responsibility for this ref
    uint16_t tstart_0;          // 0-based alignment start on target
    uint16_t aln_len_on_target;
    uint16_t tlen;
    uint8_t  pad[2];
};  // 16 bytes

// Forward declarations
class ColumnarIndexReader;
struct StreamingEMResult;
struct EMParams;

// Accumulate per-protein coverage statistics from converged EM weights.
// Recomputes gamma from em_result.weights in a streaming pass over the EMI file.
void coverage_accumulate_pass(
    ColumnarIndexReader& reader,
    const StreamingEMResult& em_result,
    const EMParams& em_params,
    const CoverageEMParams& cov_params,
    std::unordered_map<uint32_t, ProteinCoverageStats>& stats_map);

// Coverage-aware EM outer loop (streaming path).
// Alternates between inner EM (streaming_em) and coverage reweighting until
// coverage weights stabilize or max_outer_iters is reached.
StreamingEMResult coverage_em_outer_loop(
    ColumnarIndexReader& reader,
    EMParams inner_params,
    const CoverageEMParams& cov_params,
    std::unordered_map<uint32_t, ProteinCoverageStats>& out_stats,
    bool verbose = false);

// In-memory coverage accumulation from pre-sorted CovEntry vector.
// Processes one ref at a time: no per-thread unordered_maps, no heap fragmentation.
// Also computes coverage metrics (kl_divergence, coverage_weight, completeness_fhat)
// so out_stats is ready for downstream annotation.
void coverage_accumulate_from_sorted(
    const std::vector<CovEntry>& sorted_entries,
    const CoverageEMParams& cov_params,
    std::unordered_map<uint32_t, ProteinCoverageStats>& out_stats);

// In-memory coverage accumulation from AlignmentData + converged EMState (columnar path).
// Uses gamma directly from the EM state - no disk I/O required.
void coverage_accumulate_from_alignments(
    const AlignmentData& data,
    const EMState& state,
    const CoverageEMParams& cov_params,
    std::unordered_map<uint32_t, ProteinCoverageStats>& stats_map);

// Coverage-aware EM outer loop for the in-memory columnar path.
// Uses SQUAREM internally for fast convergence; returns the final converged EMState.
EMState coverage_em_squarem(
    const AlignmentData& data,
    EMParams inner_params,
    const CoverageEMParams& cov_params,
    std::unordered_map<uint32_t, ProteinCoverageStats>& out_stats,
    const EMProgressCallback& progress_cb = EMProgressCallback(),
    bool verbose = false);

}  // namespace dart
