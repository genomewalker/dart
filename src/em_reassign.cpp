// O_DIRECT (compact_em_map_tape) needs _GNU_SOURCE before the first system header.
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
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
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

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
// Empirical null model: fit / write / read
// ---------------------------------------------------------------------------

// Method-of-moments Gumbel from (n, sum, sum_sq); n-1 sample variance.
static bool fit_gumbel_mom(uint32_t n, double sum, double sum_sq, float& mu, float& beta) {
    if (n < 2) return false;
    const double mean = sum / n;
    const double var = (sum_sq - n * mean * mean) / (n - 1.0);
    if (!(var > 0.0)) return false;
    const double b = std::sqrt(6.0) * std::sqrt(var) / M_PI;
    beta = static_cast<float>(b);
    mu = static_cast<float>(mean - 0.5772156649015329 * b);
    return true;
}

NullModel fit_null_model_moments(const std::vector<uint32_t>& n,
                                 const std::vector<double>& sum,
                                 const std::vector<double>& sum_sq,
                                 uint32_t min_n)
{
    NullModel m;
    const int nbins = static_cast<int>(n.size());
    m.mu.assign(nbins, 0.0f);
    m.beta.assign(nbins, 0.0f);
    m.n = n;
    m.sum = sum;
    m.sum_sq = sum_sq;

    // Pooled fallback fit across every bin (used for sparse bins).
    uint32_t pn = 0; double ps = 0.0, pss = 0.0;
    for (int b = 0; b < nbins; ++b) { pn += n[b]; ps += sum[b]; pss += sum_sq[b]; }
    float pooled_mu = 0.0f, pooled_beta = 0.0f;
    const bool pooled_ok = fit_gumbel_mom(pn, ps, pss, pooled_mu, pooled_beta);

    for (int b = 0; b < nbins; ++b) {
        float bmu = 0.0f, bbeta = 0.0f;
        if (n[b] >= min_n && fit_gumbel_mom(n[b], sum[b], sum_sq[b], bmu, bbeta)) {
            m.mu[b] = bmu;
            m.beta[b] = bbeta;
        } else if (pooled_ok) {
            m.mu[b] = pooled_mu;
            m.beta[b] = pooled_beta;
        }
    }
    return m;
}

NullModel fit_null_model(
    const std::vector<std::vector<double>>& per_bin_scores,
    uint32_t min_n)
{
    const size_t nbins = per_bin_scores.size();
    std::vector<uint32_t> n(nbins, 0u);
    std::vector<double> sum(nbins, 0.0), sum_sq(nbins, 0.0);
    for (size_t b = 0; b < nbins; ++b) {
        n[b] = static_cast<uint32_t>(per_bin_scores[b].size());
        for (double x : per_bin_scores[b]) { sum[b] += x; sum_sq[b] += x * x; }
    }
    return fit_null_model_moments(n, sum, sum_sq, min_n);
}

// ---------------------------------------------------------------------------
// Coherent-LR calibration model: fit f_true / write / read
// ---------------------------------------------------------------------------

// Method-of-moments single Gaussian (needs n>=2 with sd>0).
static bool fit_gaussian_mom(const std::vector<double>& xs, float& mu, float& sd) {
    const size_t n = xs.size();
    if (n < 2) return false;
    double mean = 0.0;
    for (double x : xs) mean += x;
    mean /= static_cast<double>(n);
    double var = 0.0;
    for (double x : xs) { const double d = x - mean; var += d * d; }
    var /= static_cast<double>(n - 1);
    const double s = std::sqrt(var);
    if (!(s > 0.0)) return false;
    mu = static_cast<float>(mean);
    sd = static_cast<float>(s);
    return true;
}

std::vector<TrueHitBin> fit_true_model(
    const std::vector<std::vector<double>>& per_bin_bits,
    uint32_t min_n)
{
    const int nbins = static_cast<int>(per_bin_bits.size());
    std::vector<TrueHitBin> bins(nbins);

    std::vector<double> pooled;
    for (const auto& v : per_bin_bits) pooled.insert(pooled.end(), v.begin(), v.end());
    float pmu = 0.0f, psd = 1.0f;
    const bool pooled_ok = fit_gaussian_mom(pooled, pmu, psd);

    for (int b = 0; b < nbins; ++b) {
        TrueHitBin& t = bins[b];
        t.n = static_cast<uint32_t>(per_bin_bits[b].size());
        float bmu = 0.0f, bsd = 1.0f;
        if (per_bin_bits[b].size() >= min_n && fit_gaussian_mom(per_bin_bits[b], bmu, bsd)) {
            t.w0 = 1.0f; t.mu0 = bmu; t.sd0 = bsd;
        } else if (pooled_ok) {
            t.w0 = 1.0f; t.mu0 = pmu; t.sd0 = psd;
        }
    }
    return bins;
}

// Huber-reweighted least squares of y ~ a + b*x. Robust to the heavy right tail of
// true-hit bits and the occasional spurious near-homolog. Returns residual sd.
static double robust_linfit(const std::vector<double>& xs,
                            const std::vector<double>& ys,
                            double& a, double& b) {
    const size_t n = xs.size();
    a = 0.0; b = 0.0;
    if (n < 2) { if (n == 1) a = ys[0]; return 1.0; }
    std::vector<double> w(n, 1.0);
    double resid_sd = 1.0;
    for (int iter = 0; iter < 8; ++iter) {
        double sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
        for (size_t i = 0; i < n; ++i) {
            sw += w[i]; swx += w[i] * xs[i]; swy += w[i] * ys[i];
            swxx += w[i] * xs[i] * xs[i]; swxy += w[i] * xs[i] * ys[i];
        }
        const double det = sw * swxx - swx * swx;
        if (std::abs(det) < 1e-12) { b = 0.0; a = swy / std::max(sw, 1e-12); break; }
        b = (sw * swxy - swx * swy) / det;
        a = (swy - b * swx) / sw;
        // Residual scale (MAD-ish via weighted variance) → Huber weights.
        double s2 = 0.0, ssw = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const double e = ys[i] - (a + b * xs[i]);
            s2 += e * e; ssw += 1.0;
        }
        resid_sd = std::sqrt(s2 / std::max(ssw, 1.0));
        const double c = 1.345 * std::max(resid_sd, 1e-6);
        for (size_t i = 0; i < n; ++i) {
            const double e = std::abs(ys[i] - (a + b * xs[i]));
            w[i] = (e <= c) ? 1.0 : c / e;
        }
    }
    return resid_sd;
}

bool fit_calib_continuous(
    const std::vector<std::pair<double, double>>& true_LB,
    const std::vector<std::pair<double, double>>& null_LB,
    double nu, CalibModel& m)
{
    if (true_LB.size() < 10) return false;
    std::vector<double> tL, tB;
    tL.reserve(true_LB.size()); tB.reserve(true_LB.size());
    for (const auto& p : true_LB) { tL.push_back(p.first); tB.push_back(p.second); }
    double mu_a, mu_b;
    const double s0 = robust_linfit(tL, tB, mu_a, mu_b);
    m.ft_mu_a = mu_a; m.ft_mu_b = mu_b;

    // Partial-pooled log-sd(L): anchor the level at the global residual sd and let a
    // shrunk length-slope tilt it. Slope from log(residual^2)~L, shrunk toward 0 by
    // n/(n+K) so sparse data collapses to a constant sd.
    std::vector<double> lx, ly;
    double meanL = 0.0;
    for (size_t i = 0; i < tL.size(); ++i) {
        const double e = tB[i] - (mu_a + mu_b * tL[i]);
        lx.push_back(tL[i]); ly.push_back(std::log(std::max(e * e, 1e-6)));
        meanL += tL[i];
    }
    meanL /= static_cast<double>(tL.size());
    double la, lb;
    robust_linfit(lx, ly, la, lb);
    const double K = 200.0;
    const double shrink = static_cast<double>(tL.size()) / (static_cast<double>(tL.size()) + K);
    m.ft_ls_b = 0.5 * lb * shrink;
    m.ft_ls_a = std::log(std::max(s0, 1e-3)) - m.ft_ls_b * meanL;
    m.ft_nu = nu;

    // f_0 Karlin-Altschul: fit the spurious best-hit MEAN linearly in L, convert the
    // mean/sd to a Gumbel loc/scale (mean = loc + euler*beta, beta = sqrt(6)*sd/pi).
    if (null_LB.size() >= 5) {
        std::vector<double> nL, nB;
        for (const auto& p : null_LB) { nL.push_back(p.first); nB.push_back(p.second); }
        double n_a, n_b;
        const double sn = robust_linfit(nL, nB, n_a, n_b);
        const double beta = std::max(std::sqrt(6.0) * sn / M_PI, 0.3);
        m.f0_beta = beta;
        m.f0_loc_b = n_b;
        m.f0_loc_a = n_a - 0.5772156649 * beta;
    } else {
        // No spurious calibration: place a conservative low-location null just below
        // the shortest true-hit level so only genuinely null-scale reads are rejected.
        m.f0_beta = 1.5;
        m.f0_loc_b = 0.0;
        m.f0_loc_a = std::max(0.0, static_cast<double>(m.censor) - 3.0);
    }
    m.cont = true;
    return true;
}

bool write_calib_model(const std::string& path, const CalibModel& m) {
    std::ofstream out(path);
    if (!out) return false;
    out << "#dart-calib\t" << (m.cont ? "v2" : "v1") << "\tcensor=" << m.censor << '\n';
    if (m.cont) {
        out << "#cont\tft_mu_a\tft_mu_b\tft_ls_a\tft_ls_b\tft_nu\tf0_loc_a\tf0_loc_b\tf0_beta\n";
        out << "cont\t" << m.ft_mu_a << '\t' << m.ft_mu_b << '\t' << m.ft_ls_a << '\t'
            << m.ft_ls_b << '\t' << m.ft_nu << '\t' << m.f0_loc_a << '\t'
            << m.f0_loc_b << '\t' << m.f0_beta << '\n';
    }
    // f_0 section (raw-bit Gumbel per bin)
    out << "#null\tbin\tlo\thi\tmu\tbeta\tn\n";
    const int nb = static_cast<int>(m.null.mu.size());
    for (int b = 0; b < nb && b < kNumNullBins; ++b) {
        out << "null\t" << b << '\t' << kNullBinEdges[b] << '\t' << kNullBinEdges[b + 1]
            << '\t' << m.null.mu[b] << '\t' << m.null.beta[b]
            << '\t' << (b < static_cast<int>(m.null.n.size()) ? m.null.n[b] : 0u) << '\n';
    }
    // f_true section (Gaussian mixture per bin)
    out << "#true\tbin\tw0\tmu0\tsd0\tw1\tmu1\tsd1\tn\n";
    const int nt = static_cast<int>(m.true_bins.size());
    for (int b = 0; b < nt && b < kNumNullBins; ++b) {
        const TrueHitBin& t = m.true_bins[b];
        out << "true\t" << b << '\t' << t.w0 << '\t' << t.mu0 << '\t' << t.sd0
            << '\t' << t.w1 << '\t' << t.mu1 << '\t' << t.sd1 << '\t' << t.n << '\n';
    }
    return static_cast<bool>(out);
}

bool read_calib_model(const std::string& path, CalibModel& m) {
    std::ifstream in(path);
    if (!in) return false;
    m.null = NullModel{};
    m.true_bins.clear();
    m.censor = 0.0f;
    m.cont = false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line.rfind("#dart-calib", 0) == 0) {
            const auto pos = line.find("censor=");
            if (pos != std::string::npos) m.censor = std::stof(line.substr(pos + 7));
            continue;
        }
        if (line[0] == '#') continue;  // section headers
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "cont") {
            if (ss >> m.ft_mu_a >> m.ft_mu_b >> m.ft_ls_a >> m.ft_ls_b >> m.ft_nu
                   >> m.f0_loc_a >> m.f0_loc_b >> m.f0_beta)
                m.cont = true;
        } else if (tag == "null") {
            int bin; double lo, hi, mu, beta; long n = 0;
            if (!(ss >> bin >> lo >> hi >> mu >> beta)) continue;
            ss >> n;
            m.null.mu.push_back(static_cast<float>(mu));
            m.null.beta.push_back(static_cast<float>(beta));
            m.null.n.push_back(static_cast<uint32_t>(std::max(0L, n)));
        } else if (tag == "true") {
            int bin; TrueHitBin t; long n = 0;
            if (!(ss >> bin >> t.w0 >> t.mu0 >> t.sd0 >> t.w1 >> t.mu1 >> t.sd1)) continue;
            ss >> n;
            t.n = static_cast<uint32_t>(std::max(0L, n));
            m.true_bins.push_back(t);
        }
    }
    return m.valid();
}

// ---------------------------------------------------------------------------
// Parallel E-step with OpenMP
// ---------------------------------------------------------------------------

// gamma_null (per-read, size num_reads) is optional. When non-null AND the null
// model is active, each read's mixture gains one extra "null" column whose weight
// is log_w_null_dyn (the current M-step estimate of log w_null). Real gammas are
// renormalized over deg+1 columns; gamma_null[r] holds the null responsibility.
static double e_step_parallel(
    const AlignmentData& data,
    const EMParams& params,
    const double* weights,
    double* gamma,
    double* gamma_null = nullptr,
    double log_w_null_dyn = 0.0)
{
    double ll = 0.0;
    const bool null_on = params.use_null && !params.null_mu.empty() && gamma_null != nullptr;
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

            // Unique mappers: skip computation (unless the null column competes).
            if (deg == 1 && !null_on) {
                gamma[start] = 1.0;
                ll_local += std::log(std::max(weights[data.alignments[start].ref_idx], params.min_weight))
                            + length_log_penalty(data, params, data.alignments[start].ref_idx)
                            + coverage_log_penalty(params, data.alignments[start].ref_idx)
                            + real_log_density(params, static_cast<double>(data.alignments[start].bit_score), data.alignments[start].aln_len, 1.0);
                continue;
            }

            // Per-query score normalization: scale all bit scores by 1/max so scores ∈ [0,1].
            double score_scale = 1.0;
            if (params.normalize_scores_per_query && deg > 1) {
                float max_bs = data.alignments[start].bit_score;
                for (uint32_t j = 1; j < deg; ++j)
                    max_bs = std::max(max_bs, data.alignments[start + j].bit_score);
                if (max_bs > 0.0f) score_scale = 1.0 / static_cast<double>(max_bs);
            }

            const uint32_t n_cols = null_on ? deg + 1 : deg;
            log_scores.resize(n_cols);

            for (uint32_t j = 0; j < deg; ++j) {
                const auto& aln = data.alignments[start + j];
                log_scores[j] = std::log(std::max(weights[aln.ref_idx], params.min_weight))
                               + length_log_penalty(data, params, aln.ref_idx)
                               + coverage_log_penalty(params, aln.ref_idx)
                               + real_log_density(params, static_cast<double>(aln.bit_score), aln.aln_len, score_scale);
            }
            if (null_on) {
                const double lp = read_null_log_density(params, data.alignments.data() + start, deg);
                log_scores[deg] = log_w_null_dyn + params.lambda_null * lp;
            }

            double lse = log_sum_exp(log_scores.data(), n_cols);
            ll_local += lse;

            for (uint32_t j = 0; j < deg; ++j) {
                double diff = log_scores[j] - lse;
                gamma[start + j] = (diff > -700.0) ? std::exp(diff) : 0.0;
            }
            if (null_on) {
                double diff = log_scores[deg] - lse;
                gamma_null[r] = (diff > -700.0) ? std::exp(diff) : 0.0;
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

        if (deg == 1 && !null_on) {
            gamma[start] = 1.0;
            ll += std::log(std::max(weights[data.alignments[start].ref_idx], params.min_weight))
                  + length_log_penalty(data, params, data.alignments[start].ref_idx)
                  + coverage_log_penalty(params, data.alignments[start].ref_idx)
                  + real_log_density(params, static_cast<double>(data.alignments[start].bit_score), data.alignments[start].aln_len, 1.0);
            continue;
        }

        double score_scale = 1.0;
        if (params.normalize_scores_per_query && deg > 1) {
            float max_bs = data.alignments[start].bit_score;
            for (uint32_t j = 1; j < deg; ++j)
                max_bs = std::max(max_bs, data.alignments[start + j].bit_score);
            if (max_bs > 0.0f) score_scale = 1.0 / static_cast<double>(max_bs);
        }

        const uint32_t n_cols = null_on ? deg + 1 : deg;
        log_scores.resize(n_cols);
        for (uint32_t j = 0; j < deg; ++j) {
            const auto& aln = data.alignments[start + j];
            log_scores[j] = std::log(std::max(weights[aln.ref_idx], params.min_weight))
                           + length_log_penalty(data, params, aln.ref_idx)
                           + coverage_log_penalty(params, aln.ref_idx)
                           + real_log_density(params, static_cast<double>(aln.bit_score), aln.aln_len, score_scale);
        }
        if (null_on) {
            const double lp = read_null_log_density(params, data.alignments.data() + start, deg);
            log_scores[deg] = log_w_null_dyn + params.lambda_null * lp;
        }
        double lse = log_sum_exp(log_scores.data(), n_cols);
        ll += lse;
        for (uint32_t j = 0; j < deg; ++j) {
            double diff = log_scores[j] - lse;
            gamma[start + j] = (diff > -700.0) ? std::exp(diff) : 0.0;
        }
        if (null_on) {
            double diff = log_scores[deg] - lse;
            gamma_null[r] = (diff > -700.0) ? std::exp(diff) : 0.0;
        }
    }
#endif

    return ll;
}

// ---------------------------------------------------------------------------
// Parallel damage-aware E-step
// ---------------------------------------------------------------------------

// Damage-aware E-step. gamma_null (per-read, optional) adds one extra null column
// at index 2*deg; when active it competes with the 2*deg ancient/modern columns
// and gamma_null[r] receives the marginal null responsibility.
static double e_step_damage_parallel(
    const AlignmentData& data,
    const EMParams& params,
    const double* weights,
    double pi_ancient,
    double* gamma,
    double* gamma_ancient,
    double* gamma_null = nullptr,
    double log_w_null_dyn = 0.0)
{
    double ll = 0.0;
    const double eps = 1e-12;
    const double pi_a = std::clamp(pi_ancient, eps, 1.0 - eps);
    const double log_pi_a = std::log(pi_a);
    const double log_pi_m = std::log(1.0 - pi_a);
    const bool null_on = params.use_null && !params.null_mu.empty() && gamma_null != nullptr;
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

            double score_scale = 1.0;
            if (params.normalize_scores_per_query && deg > 1) {
                float max_bs = data.alignments[start].bit_score;
                for (uint32_t j = 1; j < deg; ++j)
                    max_bs = std::max(max_bs, data.alignments[start + j].bit_score);
                if (max_bs > 0.0f) score_scale = 1.0 / static_cast<double>(max_bs);
            }

            const uint32_t n_cols = null_on ? 2 * deg + 1 : 2 * deg;
            log_scores.resize(n_cols);

            for (uint32_t j = 0; j < deg; ++j) {
                const auto& aln = data.alignments[start + j];
                double log_w = std::log(std::max(weights[aln.ref_idx], params.min_weight));
                log_w += length_log_penalty(data, params, aln.ref_idx);
                log_w += coverage_log_penalty(params, aln.ref_idx);
                double score_contrib = real_log_density(params, static_cast<double>(aln.bit_score), aln.aln_len, score_scale);

                if (params.use_alignment_damage_likelihood) {
                    const double p_read = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                    log_scores[2 * j]     = log_w + score_contrib + std::log(p_read) + aln.damage_ll_a;
                    log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - p_read) + aln.damage_ll_m;
                } else {
                    double ds = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                    log_scores[2 * j]     = log_w + score_contrib + log_pi_a + std::log(ds);
                    log_scores[2 * j + 1] = log_w + score_contrib + log_pi_m + std::log(1.0 - ds);
                }
            }
            if (null_on) {
                const double lp = read_null_log_density(params, data.alignments.data() + start, deg);
                log_scores[2 * deg] = log_w_null_dyn + params.lambda_null * lp;
            }

            double lse = log_sum_exp(log_scores.data(), n_cols);
            ll_local += lse;

            for (uint32_t j = 0; j < deg; ++j) {
                double diff_a = log_scores[2 * j] - lse;
                double diff_m = log_scores[2 * j + 1] - lse;
                double g_a = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
                double g_m = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;

                gamma[start + j] = g_a + g_m;
                gamma_ancient[start + j] = g_a;
            }
            if (null_on) {
                double diff = log_scores[2 * deg] - lse;
                gamma_null[r] = (diff > -700.0) ? std::exp(diff) : 0.0;
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

        double score_scale = 1.0;
        if (params.normalize_scores_per_query && deg > 1) {
            float max_bs = data.alignments[start].bit_score;
            for (uint32_t j = 1; j < deg; ++j)
                max_bs = std::max(max_bs, data.alignments[start + j].bit_score);
            if (max_bs > 0.0f) score_scale = 1.0 / static_cast<double>(max_bs);
        }

        const uint32_t n_cols = null_on ? 2 * deg + 1 : 2 * deg;
        log_scores.resize(n_cols);
        for (uint32_t j = 0; j < deg; ++j) {
            const auto& aln = data.alignments[start + j];
            double log_w = std::log(std::max(weights[aln.ref_idx], params.min_weight));
            log_w += length_log_penalty(data, params, aln.ref_idx);
            double score_contrib = real_log_density(params, static_cast<double>(aln.bit_score), aln.aln_len, score_scale);
            if (params.use_alignment_damage_likelihood) {
                const double p_read = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                log_scores[2 * j]     = log_w + score_contrib + std::log(p_read) + aln.damage_ll_a;
                log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - p_read) + aln.damage_ll_m;
            } else {
                double ds = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                log_scores[2 * j]     = log_w + score_contrib + log_pi_a + std::log(ds);
                log_scores[2 * j + 1] = log_w + score_contrib + log_pi_m + std::log(1.0 - ds);
            }
        }
        if (null_on) {
            const double lp = read_null_log_density(params, data.alignments.data() + start, deg);
            log_scores[2 * deg] = log_w_null_dyn + params.lambda_null * lp;
        }

        double lse = log_sum_exp(log_scores.data(), n_cols);
        ll += lse;
        for (uint32_t j = 0; j < deg; ++j) {
            double diff_a = log_scores[2 * j] - lse;
            double diff_m = log_scores[2 * j + 1] - lse;
            double g_a = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
            double g_m = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;
            gamma[start + j] = g_a + g_m;
            gamma_ancient[start + j] = g_a;
        }
        if (null_on) {
            double diff = log_scores[2 * deg] - lse;
            gamma_null[r] = (diff > -700.0) ? std::exp(diff) : 0.0;
        }
    }
#endif

    return ll;
}

// ---------------------------------------------------------------------------
// Parallel M-step with thread-local accumulators
// ---------------------------------------------------------------------------

// gamma_null (per-read, optional) contributes a null coordinate to the simplex:
//   Σ_t w_t + w_null = 1,  w_null ∝ (alpha_null - 1) + Σ_r gamma_null_r.
// When out_w_null is non-null it receives the normalized null weight.
static void m_step_parallel(
    const AlignmentData& data,
    const EMParams& params,
    const double* gamma,
    double* out_weights,
    const double* gamma_null = nullptr,
    double* out_w_null = nullptr)
{
    const uint32_t T = data.num_refs;
    const bool null_on = params.use_null && gamma_null != nullptr && out_w_null != nullptr;

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

    // Null coordinate: informative Dirichlet floor (alpha_null) prevents the null
    // from collapsing to 0 or absorbing all mass.
    double w_null = 0.0;
    if (null_on) {
        double null_mass = params.alpha_null - 1.0;
        for (uint32_t r = 0; r < data.num_reads; ++r) null_mass += gamma_null[r];
        w_null = std::max(null_mass, params.min_weight);
        sum += w_null;
    }

    double inv_sum = 1.0 / sum;
#pragma omp simd
    for (uint32_t t = 0; t < T; ++t) {
        out_weights[t] *= inv_sum;
    }
    if (null_on) *out_w_null = w_null * inv_sum;
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

    // Null component: augment the weight simplex with one extra coordinate at
    // index T holding w_null. All SQUAREM vector ops run over Tn = T+1 coords so
    // the null participates in extrapolation and renormalization. When use_null is
    // off, Tn == T and every code path is byte-identical to before.
    const bool null_on = params.use_null && !params.null_mu.empty();
    const uint32_t Tn = null_on ? T + 1 : T;

    EMState state;
    state.gamma.resize(A, 0.0);
    if (params.use_damage) {
        state.gamma_ancient.resize(A, 0.0);
    }
    if (null_on) {
        state.gamma_null.assign(data.num_reads, 0.0);
    }
    if (T == 0) {
        state.log_likelihood = 0.0;
        state.iterations = 0;
        return state;
    }
    if (!params.initial_weights.empty() && params.initial_weights.size() == T) {
        state.weights.assign(Tn, 0.0);
        double sum = 0.0;
        for (uint32_t t = 0; t < T; ++t) sum += params.initial_weights[t];
        for (uint32_t t = 0; t < T; ++t)
            state.weights[t] = (sum > 0.0) ? params.initial_weights[t] / sum : (1.0 / T);
    } else {
        state.weights.assign(Tn, 1.0 / static_cast<double>(T));
    }
    // Initialize w_null to the (renormalized) prior weight exp(log_w_null).
    if (null_on) {
        double w_null0 = std::clamp(std::exp(params.log_w_null), params.min_weight, 1.0 - 1e-6);
        for (uint32_t t = 0; t < T; ++t) state.weights[t] *= (1.0 - w_null0);
        state.weights[T] = w_null0;
    }

    // SQUAREM scratch buffers (sized Tn so the null coordinate is carried through)
    std::vector<double> w0(Tn), w1(Tn), w2(Tn), r_buf(Tn), v_buf(Tn);
    std::vector<double> prev_weights(Tn);
    std::vector<double> gamma_scratch(A);
    std::vector<double> gamma_ancient_scratch;
    if (params.use_damage) {
        gamma_ancient_scratch.resize(A, 0.0);
    }
    // Per-read null responsibility scratch (only allocated when null is active).
    std::vector<double> gamma_null_scratch;
    if (null_on) gamma_null_scratch.assign(data.num_reads, 0.0);

    // Convenience lambdas that dispatch E/M with the null coordinate wired in.
    auto do_e_step = [&](const double* w, double* gamma, double* gamma_anc, double* gnull) {
        const double log_w_null_dyn = null_on
            ? std::log(std::max(w[T], params.min_weight)) : 0.0;
        if (params.use_damage) {
            return e_step_damage_parallel(data, params, w, state.pi, gamma, gamma_anc,
                                          null_on ? gnull : nullptr, log_w_null_dyn);
        }
        return e_step_parallel(data, params, w, gamma,
                               null_on ? gnull : nullptr, log_w_null_dyn);
    };
    auto do_m_step = [&](const double* gamma, double* out_w, const double* gnull) {
        m_step_parallel(data, params, gamma, out_w,
                        null_on ? gnull : nullptr,
                        null_on ? out_w + T : nullptr);
    };

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
            do_e_step(w0.data(), gamma_scratch.data(), gamma_ancient_scratch.data(),
                      gamma_null_scratch.data());
            do_m_step(gamma_scratch.data(), w1.data(), gamma_null_scratch.data());

            // Second EM step: w1 -> w2
            do_e_step(w1.data(), gamma_scratch.data(), gamma_ancient_scratch.data(),
                      gamma_null_scratch.data());
            do_m_step(gamma_scratch.data(), w2.data(), gamma_null_scratch.data());

            // SQUAREM extrapolation with textbook step rules.
            double sr2 = 0.0, sv2 = 0.0, srv = 0.0;
            squarem_build_directions(
                w0.data(), w1.data(), w2.data(), r_buf.data(), v_buf.data(), Tn, sr2, sv2, srv);
            double alpha = squarem_step_length(params.squarem_method, sr2, sv2, srv);
            alpha = std::clamp(alpha, step_min, step_max);
            alpha_used = alpha;
            squarem_extrapolate(w0.data(), r_buf.data(), v_buf.data(), Tn, alpha, params.min_weight);

            // Stabilization step from canonical SQUAREM (Varadhan & Roland 2008).
            // Always applied after extrapolation to maintain monotonicity guarantee.
            {
                do_e_step(w0.data(), gamma_scratch.data(), gamma_ancient_scratch.data(),
                          gamma_null_scratch.data());
                do_m_step(gamma_scratch.data(), w1.data(), gamma_null_scratch.data());
                std::copy(w1.begin(), w1.end(), w0.begin());
            }

            // Evaluate at accelerated point
            ll = do_e_step(w0.data(), state.gamma.data(), state.gamma_ancient.data(),
                           state.gamma_null.data());
            objective = em_objective_from_ll(params, w0.data(), T, ll);

            // Objective safeguard: revert to EM two-step point when acceleration is non-monotone.
            const bool non_monotone =
                (iter > 0) &&
                std::isfinite(params.squarem_objfn_inc) &&
                (objective < (prev_obj - params.squarem_objfn_inc));
            if (non_monotone) {
                fallback_to_em = true;
                std::copy(w2.begin(), w2.end(), state.weights.begin());
                ll = do_e_step(state.weights.data(), state.gamma.data(),
                               state.gamma_ancient.data(), state.gamma_null.data());
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
            ll = do_e_step(state.weights.data(), state.gamma.data(),
                           state.gamma_ancient.data(), state.gamma_null.data());
            objective = em_objective_from_ll(params, state.weights.data(), T, ll);
            do_m_step(state.gamma.data(), state.weights.data(), state.gamma_null.data());
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
            fill_weight_diagnostics(diag, state.weights.data(), T, params.min_weight);
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

    // Detach the null coordinate: expose it as state.w_null and shrink weights
    // back to T so all downstream ref-indexed consumers see only real references.
    if (null_on) {
        state.w_null = state.weights[T];
        state.weights.resize(T);
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
// compact_em_map (fusion dual-b8247d43): the SINGLE EM solver. Decode-once CSR
// {ref_idx, base_a, base_m} where base = length_log_penalty + coverage_log_penalty +
// real_log_density, split into ancient/modern by em_base_pair (damage enters the
// mapping only with use_damage), all iteration-independent. Per-iter E-step adds only log(weights[ref]); fused per-ref
// s_t = Σγ; full SQUAREM; DETERMINISTIC-LANE M-step (fixed lanes, read-order per lane,
// ascending-lane reduce) -> bitwise-identical across thread counts. Frozen-damage branch
// (use_alignment_damage_likelihood): pi is post-hoc, only w_null rides the SQUAREM state.
static constexpr int COMPACT_EM_LANES = 32;
static constexpr uint32_t kEmPfDist = 16;   // gather/scatter lookahead, in alignment slots (OPT-C)

EMState compact_em_map(AlignmentData& data, const EMParams& params,
                       const TapeFinalizeCb& finalize_cb) {
    EMState state;
    const uint32_t T = data.num_refs;
    const size_t A = data.alignments.size();
    const uint32_t R = data.num_reads;
    const bool null_on = params.use_null && !params.null_mu.empty();
    const uint32_t Tn = null_on ? T + 1 : T;
    const double eps = 1e-12;

    // gamma is EM-TRANSIENT: iterations accumulate only per-ref sufficient stats
    // (lane_s) + optional per-read null mass; the full gamma is materialized ONLY in
    // the single final frozen-weight pass. When finalize_cb is supplied it is streamed
    // per read and NO O(alignments) gamma is ever held (production, ~O(num_refs) RSS);
    // otherwise the final pass fills state.gamma exactly as before (gate/back-compat).
    const bool stream = static_cast<bool>(finalize_cb);
    const bool want_anc = params.use_damage || params.posthoc_gamma_ancient;
    if (!stream) {
        state.gamma.assign(A, 0.0);
        if (want_anc) state.gamma_ancient.assign(A, 0.0);
    }
    if (null_on) state.gamma_null.assign(R, 0.0);
    if (T == 0) { state.log_likelihood = 0.0; state.iterations = 0; return state; }

    // ---- decode ONCE: base_a/base_m per alignment (iteration-independent) ----
    std::vector<uint32_t> csr_ref(A);
    std::vector<double> base_a(A), base_m(A);
    #pragma omp parallel for schedule(static)
    for (uint32_t r = 0; r < R; ++r) {
        const uint32_t s = data.read_offsets[r], e = data.read_offsets[r + 1];
        const uint32_t deg = e - s;
        if (deg == 0) continue;
        double score_scale = 1.0;
        if (params.normalize_scores_per_query && deg > 1) {
            float mx = data.alignments[s].bit_score;
            for (uint32_t j = 1; j < deg; ++j) mx = std::max(mx, data.alignments[s + j].bit_score);
            if (mx > 0.0f) score_scale = 1.0 / static_cast<double>(mx);
        }
        for (uint32_t j = s; j < e; ++j) {
            const auto& aln = data.alignments[j];
            csr_ref[j] = aln.ref_idx;
            double bc = length_log_penalty(data, params, aln.ref_idx)
                      + coverage_log_penalty(params, aln.ref_idx)
                      + real_log_density(params, static_cast<double>(aln.bit_score), aln.aln_len, score_scale);
            const double pr = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
            em_base_pair(params, bc, pr, aln.damage_ll_a, aln.damage_ll_m, base_a[j], base_m[j]);
        }
    }

    // Streaming path with the null column off never reads data.alignments again (decode
    // is done; iterations + final pass use only csr_ref/base_a/base_m). Release the
    // 32 B/aln ingestion copy NOW so steady-state RSS is ~the 20 B/aln CSR alone. The
    // non-streaming (gate/back-compat) path keeps it — callers may reuse data.
    if (stream && !null_on) std::vector<Alignment>{}.swap(data.alignments);

    // ---- init weights (+ null coordinate), same as squarem_em ----
    state.weights.assign(Tn, 1.0 / static_cast<double>(T));
    if (null_on) {
        double w0n = std::clamp(std::exp(params.log_w_null), params.min_weight, 1.0 - 1e-6);
        for (uint32_t t = 0; t < T; ++t) state.weights[t] *= (1.0 - w0n);
        state.weights[T] = w0n;
    }

    // deterministic lane read ranges (contiguous CSR ranges)
    std::array<uint32_t, COMPACT_EM_LANES + 1> lane_lo;
    for (int l = 0; l <= COMPACT_EM_LANES; ++l)
        lane_lo[l] = static_cast<uint32_t>((static_cast<uint64_t>(R) * l) / COMPACT_EM_LANES);
    std::vector<double> lane_s(static_cast<size_t>(COMPACT_EM_LANES) * T, 0.0);
    std::vector<double> gn_s(null_on ? R : 0, 0.0);

    // Fused E-step -> per-lane s_t (deterministic) + optional null mass; returns ll.
    // gamma is NOT stored here (transient): only lane_s/gn_s survive the iteration.
    auto e_step_iter = [&](const double* w) -> double {
        const double log_w_null_dyn = null_on ? std::log(std::max(w[T], params.min_weight)) : 0.0;
        std::fill(lane_s.begin(), lane_s.end(), 0.0);
        std::vector<double> ll_lane(COMPACT_EM_LANES, 0.0);
        #pragma omp parallel for schedule(dynamic, 1)
        for (int l = 0; l < COMPACT_EM_LANES; ++l) {
            double* s_t = lane_s.data() + static_cast<size_t>(l) * T;
            double ll_local = 0.0;
            std::vector<double> ls;
            for (uint32_t r = lane_lo[l]; r < lane_lo[l + 1]; ++r) {
                const uint32_t s = data.read_offsets[r], e = data.read_offsets[r + 1];
                const uint32_t deg = e - s;
                if (deg == 0) continue;
                const uint32_t ncol = null_on ? 2 * deg + 1 : 2 * deg;
                ls.resize(ncol);
                for (uint32_t j = 0; j < deg; ++j) {
                    const double lw = std::log(std::max(w[csr_ref[s + j]], params.min_weight));
                    ls[2 * j] = lw + base_a[s + j];
                    ls[2 * j + 1] = lw + base_m[s + j];
                }
                if (null_on) {
                    const double lp = read_null_log_density(params, data.alignments.data() + s, deg);
                    ls[2 * deg] = log_w_null_dyn + params.lambda_null * lp;
                }
                const double lse = log_sum_exp(ls.data(), ncol);
                ll_local += lse;
                for (uint32_t j = 0; j < deg; ++j) {
                    const double da = ls[2 * j] - lse, dm = ls[2 * j + 1] - lse;
                    const double ga = (da > -700.0) ? std::exp(da) : 0.0;
                    const double gm = (dm > -700.0) ? std::exp(dm) : 0.0;
                    s_t[csr_ref[s + j]] += ga + gm;   // read-order accumulate within lane
                }
                if (null_on) {
                    const double d = ls[2 * deg] - lse;
                    gn_s[r] = (d > -700.0) ? std::exp(d) : 0.0;
                }
            }
            ll_lane[l] = ll_local;
        }
        double ll = 0.0; for (int l = 0; l < COMPACT_EM_LANES; ++l) ll += ll_lane[l];
        return ll;
    };
    // M-step: weights[t] = (prior + Σ_lane s_lane[l][t]) ascending-lane, then normalize.
    auto m_step_fused = [&](double* out_w) {
        const double prior_offset = params.alpha_prior - 1.0;
        double tot = 0.0;
        for (uint32_t t = 0; t < T; ++t) {
            double w = prior_offset;
            for (int l = 0; l < COMPACT_EM_LANES; ++l) w += lane_s[static_cast<size_t>(l) * T + t];
            w = std::max(w, params.min_weight);
            out_w[t] = w; tot += w;
        }
        double wn = 0.0;
        if (null_on) {
            wn = prior_offset;
            for (uint32_t r = 0; r < R; ++r) wn += gn_s[r];
            wn = std::max(wn, params.min_weight); tot += wn;
        }
        if (tot > 0.0) { for (uint32_t t = 0; t < T; ++t) out_w[t] /= tot; if (null_on) out_w[T] = wn / tot; }
    };
    auto do_e = [&](const double* w) { return e_step_iter(w); };
    auto do_m = [&](double* out) { m_step_fused(out); };

    std::vector<double> w0(Tn), w1(Tn), w2(Tn), r_buf(Tn), v_buf(Tn), prev_w(Tn);
    double prev_obj = -std::numeric_limits<double>::infinity();
    double step_min = params.squarem_step_min0, step_max0 = params.squarem_step_max0;
    if (step_max0 < step_min) std::swap(step_max0, step_min);
    double step_max = step_max0;
    const double mstep = std::max(params.squarem_mstep, 1.0000001);
    const double param_tol = std::sqrt(std::max(params.tol, 1e-15));

    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        double ll, objective;
        std::copy(state.weights.begin(), state.weights.end(), prev_w.begin());
        const bool sq = params.use_squarem && (iter >= params.squarem_warmup_iters);
        if (sq) {
            std::copy(state.weights.begin(), state.weights.end(), w0.begin());
            do_e(w0.data()); do_m(w1.data());
            do_e(w1.data()); do_m(w2.data());
            double sr2 = 0, sv2 = 0, srv = 0;
            squarem_build_directions(w0.data(), w1.data(), w2.data(), r_buf.data(), v_buf.data(), Tn, sr2, sv2, srv);
            double alpha = std::clamp(squarem_step_length(params.squarem_method, sr2, sv2, srv), step_min, step_max);
            squarem_extrapolate(w0.data(), r_buf.data(), v_buf.data(), Tn, alpha, params.min_weight);
            do_e(w0.data()); do_m(w1.data()); std::copy(w1.begin(), w1.end(), w0.begin());
            ll = do_e(w0.data());
            objective = em_objective_from_ll(params, w0.data(), T, ll);
            const bool nonmono = (iter > 0) && std::isfinite(params.squarem_objfn_inc)
                                 && (objective < (prev_obj - params.squarem_objfn_inc));
            if (nonmono) {
                std::copy(w2.begin(), w2.end(), state.weights.begin());
                ll = do_e(state.weights.data());
                objective = em_objective_from_ll(params, state.weights.data(), T, ll);
                if (std::abs(alpha - step_max) <= (1e-12 * (1.0 + std::abs(step_max)))) step_max = std::max(step_max0, step_max / mstep);
            } else {
                std::copy(w0.begin(), w0.end(), state.weights.begin());
            }
            if (std::abs(alpha - step_max) <= (1e-12 * (1.0 + std::abs(step_max)))) step_max = mstep * step_max;
            if (step_min < 0.0 && std::abs(alpha - step_min) <= (1e-12 * (1.0 + std::abs(step_min)))) step_min = mstep * step_min;
        } else {
            ll = do_e(state.weights.data());
            objective = em_objective_from_ll(params, state.weights.data(), T, ll);
            do_m(state.weights.data());
        }
        state.log_likelihood = ll; state.iterations = iter + 1;
        double rel = std::numeric_limits<double>::infinity();
        if (iter > 0) rel = std::abs(objective - prev_obj) / (std::abs(prev_obj) + 1e-15);
        double pcs = 0.0; for (uint32_t t = 0; t < T; ++t) { double d = state.weights[t] - prev_w[t]; pcs += d * d; }
        if (iter > 0 && (iter + 1) >= params.min_iters && rel <= params.tol && std::sqrt(pcs) <= param_tol) break;
        prev_obj = objective;
    }

    // ---- single SERIAL final pass at frozen weights: gamma is transient (per read).
    // Recomputing at the converged weights reproduces the previous per-iteration stored
    // gamma bit-for-bit (convergence always lands in the SQUAREM branch, whose stored
    // gamma is likewise at state.weights). When streaming, hand each read to finalize_cb
    // and hold no O(A) gamma; else fill state.gamma / state.gamma_ancient as before.
    // gamma_ancient: use_damage -> ga (joint); posthoc_gamma_ancient -> gamma * p_anc with
    // p_anc = sigmoid(base_a - base_m) (bc cancels; == the CLI's post-hoc formula; the
    // damage-free split in em_base_pair keeps this difference). ----
    {
        const double log_w_null_dyn = null_on ? std::log(std::max(state.weights[T], params.min_weight)) : 0.0;
        std::vector<double> ls, gbuf, gabuf;
        std::vector<uint32_t> rbuf;
        double final_ll = 0.0, sum_a = 0.0, sum_t = 0.0;
        for (uint32_t r = 0; r < R; ++r) {
            const uint32_t s = data.read_offsets[r], e = data.read_offsets[r + 1];
            const uint32_t deg = e - s;
            if (null_on) state.gamma_null[r] = 0.0;
            if (deg == 0) { if (stream) finalize_cb(r, 0, nullptr, nullptr, nullptr); continue; }
            const uint32_t ncol = null_on ? 2 * deg + 1 : 2 * deg;
            ls.resize(ncol);
            if (stream) { gbuf.resize(deg); rbuf.resize(deg); if (want_anc) gabuf.resize(deg); }
            for (uint32_t j = 0; j < deg; ++j) {
                const uint32_t rf = csr_ref[s + j];
                if (stream) rbuf[j] = rf;
                const double lw = std::log(std::max(state.weights[rf], params.min_weight));
                ls[2 * j] = lw + base_a[s + j];
                ls[2 * j + 1] = lw + base_m[s + j];
            }
            if (null_on) {
                const double lp = read_null_log_density(params, data.alignments.data() + s, deg);
                ls[2 * deg] = log_w_null_dyn + params.lambda_null * lp;
            }
            const double lse = log_sum_exp(ls.data(), ncol);
            final_ll += lse;
            for (uint32_t j = 0; j < deg; ++j) {
                const double da = ls[2 * j] - lse, dm = ls[2 * j + 1] - lse;
                const double ga = (da > -700.0) ? std::exp(da) : 0.0;
                const double gm = (dm > -700.0) ? std::exp(dm) : 0.0;
                const double g = ga + gm;
                double ganc = 0.0;
                if (params.use_damage) {
                    ganc = ga;
                } else if (params.posthoc_gamma_ancient) {
                    const double p_anc = 1.0 / (1.0 + std::exp(base_m[s + j] - base_a[s + j]));
                    ganc = g * p_anc;
                }
                if (stream) { gbuf[j] = g; if (want_anc) gabuf[j] = ganc; }
                else { state.gamma[s + j] = g; if (want_anc) state.gamma_ancient[s + j] = ganc; }
                sum_a += ganc; sum_t += g;
            }
            if (null_on) {
                const double d = ls[2 * deg] - lse;
                state.gamma_null[r] = (d > -700.0) ? std::exp(d) : 0.0;
            }
            if (stream)
                finalize_cb(r, deg, rbuf.data(), gbuf.data(), want_anc ? gabuf.data() : nullptr);
        }
        state.log_likelihood = final_ll;
        // pi is the joint-EM ancient fraction ONLY (matches the old per-iter update);
        // in the post-hoc (damage-free mapping) path pi stays at its 0.10 prior, exactly
        // as before — the post-hoc formula consumes that prior, it does not redefine pi.
        (void)sum_a;
        if (params.use_damage && sum_t > 0.0) state.pi = std::clamp(sum_a / sum_t, 0.01, 0.99);
    }

    if (null_on) { state.w_null = state.weights[T]; state.weights.resize(T); }
    return state;
}

// ===========================================================================
// compact_em_map_tape (PI mandate 2026-07-23): decode-once, read-major "EM
// tape" on scratch NVMe, streamed per iteration with O_DIRECT. Reproduces
// compact_em_map BIT-FOR-BIT (identical decode, identical 32 contiguous-read
// lanes, identical read-order + ascending-lane reductions, identical SQUAREM),
// but never holds an O(alignments) gamma/base vector in RAM: per-iteration the
// E-step accumulates only per-ref sufficient stats (lane_s), gamma is transient
// and materialized ONLY in a single final pass. RAM = O(num_refs) + fixed I/O
// buffers => independent of alignment count.
//
// Tape = 3 SoA streams (each element indexed by CSR alignment position j):
//   <dir>/tape.ref   : uint32 ref_idx[j]
//   <dir>/tape.basea : double  base_a[j]   (== compact_em_map's base_a)
//   <dir>/tape.basem : double  base_m[j]   (== compact_em_map's base_m)
// SoA (not AoS) so every O_DIRECT request is block-aligned with zero record
// straddle: 4096 % 4 == 0 and 4096 % 8 == 0, and a chunk start aligned to a
// multiple of 1024 records makes both the 4B and 8B streams 4096-aligned.
// ===========================================================================

namespace {

constexpr size_t kODAlign = 4096;                 // O_DIRECT block alignment
constexpr uint64_t kTapeChunkRecs = 1u << 18;     // 262144 (mult. of 1024) records/refill

// Write a typed array to `path` with O_DIRECT, padding the file tail to kODAlign.
// Readers bound access by record count, so the zero padding is never consumed.
static bool tape_write_array(const std::string& path, const void* data, size_t nbytes) {
    int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_DIRECT, 0644);
    if (fd < 0) { std::perror(("tape_write open " + path).c_str()); return false; }
    const size_t BUF = size_t(1) << 22;  // 4 MiB (multiple of kODAlign)
    void* buf = nullptr;
    if (posix_memalign(&buf, kODAlign, BUF) != 0) { ::close(fd); return false; }
    bool ok = true;
    size_t off = 0;
    while (off < nbytes) {
        size_t n = std::min(BUF, nbytes - off);
        std::memcpy(buf, static_cast<const char*>(data) + off, n);
        size_t wlen = (n + kODAlign - 1) & ~(kODAlign - 1);
        if (wlen > n) std::memset(static_cast<char*>(buf) + n, 0, wlen - n);
        ssize_t w = ::pwrite(fd, buf, wlen, static_cast<off_t>(off));
        if (w < 0 || static_cast<size_t>(w) != wlen) { ok = false; std::perror("tape_write pwrite"); break; }
        off += n;  // stays kODAlign-aligned every iter except the (final) short one
    }
    free(buf);
    if (::fsync(fd) != 0) ok = false;
    ::close(fd);
    return ok;
}

// One O_DIRECT-backed SoA stream reader with a per-lane sliding record window.
// Records [lo, hi) are served by refilling a kTapeChunkRecs window whose start
// is snapped down to a 1024-record boundary (=> byte offset is kODAlign-aligned
// for both 4B and 8B element sizes). ptr(j) returns a typed pointer to record j.
template <class T>
struct TapeStream {
    int fd = -1;
    T* buf = nullptr;          // kODAlign-aligned
    uint64_t win_lo = 0;       // first record currently buffered
    uint64_t win_n = 0;        // records currently buffered
    uint64_t file_recs = 0;    // total valid records in the file

    bool open(const std::string& path, uint64_t total_recs) {
        file_recs = total_recs;
        fd = ::open(path.c_str(), O_RDONLY | O_DIRECT);
        if (fd < 0) { std::perror(("tape_read open " + path).c_str()); return false; }
        size_t bytes = size_t(kTapeChunkRecs) * sizeof(T);
        void* p = nullptr;
        if (posix_memalign(&p, kODAlign, bytes) != 0) { ::close(fd); fd = -1; return false; }
        buf = static_cast<T*>(p);
        win_lo = 0; win_n = 0;
        return true;
    }
    void close() { if (buf) { free(buf); buf = nullptr; } if (fd >= 0) { ::close(fd); fd = -1; } }

    // Ensure records [s, e) are resident; snap window start to 1024-record boundary.
    void ensure(uint64_t s, uint64_t e) {
        if (s >= win_lo && e <= win_lo + win_n) return;
        uint64_t start = s & ~uint64_t(1023);
        uint64_t want = std::min<uint64_t>(kTapeChunkRecs, file_recs - start);
        // Round read length up to kODAlign; safe because the file tail is padded.
        size_t nbytes = size_t(want) * sizeof(T);
        size_t rlen = (nbytes + kODAlign - 1) & ~(kODAlign - 1);
        off_t offb = static_cast<off_t>(start * sizeof(T));
        ssize_t got = ::pread(fd, buf, rlen, offb);
        if (got < 0) { std::perror("tape_read pread"); std::abort(); }
        win_lo = start;
        win_n = std::min<uint64_t>(want, static_cast<uint64_t>(got) / sizeof(T));
    }
    inline const T& at(uint64_t j) const { return buf[j - win_lo]; }
};

}  // anonymous namespace

// Fill base_a[]/base_m[]/csr_ref[] for one read exactly as compact_em_map's decode.
// Shared by the tape writer so the on-disk base is bit-identical to the in-RAM path.
static inline void tape_decode_read(
    const AlignmentData& data, const EMParams& params, uint32_t r,
    uint32_t* ref_out, double* basea_out, double* basem_out) {
    const uint32_t s = data.read_offsets[r], e = data.read_offsets[r + 1];
    const uint32_t deg = e - s;
    if (deg == 0) return;
    const double eps = 1e-12;
    double score_scale = 1.0;
    if (params.normalize_scores_per_query && deg > 1) {
        float mx = data.alignments[s].bit_score;
        for (uint32_t j = 1; j < deg; ++j) mx = std::max(mx, data.alignments[s + j].bit_score);
        if (mx > 0.0f) score_scale = 1.0 / static_cast<double>(mx);
    }
    for (uint32_t j = s; j < e; ++j) {
        const auto& aln = data.alignments[j];
        ref_out[j] = aln.ref_idx;
        double bc = length_log_penalty(data, params, aln.ref_idx)
                  + coverage_log_penalty(params, aln.ref_idx)
                  + real_log_density(params, static_cast<double>(aln.bit_score), aln.aln_len, score_scale);
        const double pr = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
        em_base_pair(params, bc, pr, aln.damage_ll_a, aln.damage_ll_m, basea_out[j], basem_out[j]);
    }
}

// ---- EM fixed-point checkpoint (OPT-G) ----------------------------------
// Only the SQUAREM iterations are checkpointed. The final streaming pass is NOT:
// it re-reads the whole tape and drives finalize_cb, so it always re-runs, and a
// resumed run is therefore bit-identical (final pass is a pure function of weights).
static constexpr uint64_t TAPE_CKPT_MAGIC   = 0x5441504543484B31ULL;  // "TAPECHK1"
static constexpr uint32_t TAPE_CKPT_VERSION = 2;  // 2: damage-free mapping base (em_base_pair)

struct TapeCkptHeader {
    uint64_t magic;
    uint32_t version;
    uint32_t n_weights;     // == Tn (T, or T+1 when the null column is on)
    uint64_t fingerprint;
    double   log_likelihood;
    uint32_t iterations;
    uint32_t pad;
};

static inline uint64_t tape_ckpt_fnv1a(const void* p, size_t n, uint64_t h) {
    const unsigned char* b = static_cast<const unsigned char*>(p);
    for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}

// Hashes exactly the inputs that move the EM fixed point: the tape identity
// (caller-supplied input_fingerprint) + shape (R, T, A) + the EMParams read by
// the iteration loop / stop test. use_damage / posthoc_gamma_ancient are excluded:
// they only steer the final pass, which always re-runs.
static uint64_t tape_ckpt_fingerprint(const std::string& input_fingerprint,
                                      uint32_t R, uint32_t T, uint64_t A,
                                      const EMParams& params, bool null_on) {
    uint64_t h = 14695981039346656037ULL;
    h = tape_ckpt_fnv1a(input_fingerprint.data(), input_fingerprint.size(), h);
    h = tape_ckpt_fnv1a(&R, sizeof(R), h);
    h = tape_ckpt_fnv1a(&T, sizeof(T), h);
    h = tape_ckpt_fnv1a(&A, sizeof(A), h);
    h = tape_ckpt_fnv1a(&params.min_weight, sizeof(double), h);
    h = tape_ckpt_fnv1a(&params.alpha_prior, sizeof(double), h);
    h = tape_ckpt_fnv1a(&params.max_iters, sizeof(uint32_t), h);
    h = tape_ckpt_fnv1a(&params.min_iters, sizeof(uint32_t), h);
    h = tape_ckpt_fnv1a(&params.tol, sizeof(double), h);
    const uint8_t sq = params.use_squarem ? 1 : 0;
    h = tape_ckpt_fnv1a(&sq, 1, h);
    h = tape_ckpt_fnv1a(&params.squarem_warmup_iters, sizeof(uint32_t), h);
    h = tape_ckpt_fnv1a(&params.squarem_method, sizeof(int), h);
    h = tape_ckpt_fnv1a(&params.squarem_step_min0, sizeof(double), h);
    h = tape_ckpt_fnv1a(&params.squarem_step_max0, sizeof(double), h);
    h = tape_ckpt_fnv1a(&params.squarem_mstep, sizeof(double), h);
    h = tape_ckpt_fnv1a(&params.squarem_objfn_inc, sizeof(double), h);
    const uint8_t nl = null_on ? 1 : 0;
    h = tape_ckpt_fnv1a(&nl, 1, h);
    h = tape_ckpt_fnv1a(&params.lambda_null, sizeof(double), h);
    h = tape_ckpt_fnv1a(&params.log_w_null, sizeof(double), h);
    if (null_on) {   // the null column's Gumbel table changes the converged weights
        h = tape_ckpt_fnv1a(&params.alpha_null, sizeof(double), h);
        for (float v : params.null_mu) h = tape_ckpt_fnv1a(&v, sizeof(float), h);
        for (float v : params.null_beta) h = tape_ckpt_fnv1a(&v, sizeof(float), h);
    }
    return h;
}

// Shared streaming core for the EM tape. Both the iterations and the single final
// pass run over the same 32 contiguous-read lanes in parallel (per-ref transient
// sufficient stats during iterations; per-read transient gamma in the final pass,
// streamed to finalize_cb with no O(alignments) state). Determinism (bitwise across
// thread counts) comes from the fixed lane boundaries + read-order within-lane +
// ascending-lane reduction — the lane geometry is identical to compact_em_map, but
// the final-pass scalar reductions (log_likelihood, gamma_sum, gamma_anc_sum) are
// lane-partitioned here and so differ in the last ulps from compact_em_map's serial
// running sum. Per-read gamma is bit-identical.
TapeEMResult tape_em_core(
    const std::string& p_ref, const std::string& p_basea, const std::string& p_basem,
    const uint64_t* read_offsets, uint32_t R, uint32_t T, uint64_t A,
    const EMParams& params, const double* null_lp, double* gnull_out,
    const TapeFinalizeCb& finalize_cb,
    const std::string& checkpoint_path, bool resume,
    const std::string& input_fingerprint)
{
    TapeEMResult res;
    const bool null_on = params.use_null && null_lp != nullptr;
    const uint32_t Tn = null_on ? T + 1 : T;
    if (T == 0) return res;

    // deterministic lane read ranges (identical to compact_em_map)
    std::array<uint32_t, COMPACT_EM_LANES + 1> lane_lo;
    for (int l = 0; l <= COMPACT_EM_LANES; ++l)
        lane_lo[l] = static_cast<uint32_t>((static_cast<uint64_t>(R) * l) / COMPACT_EM_LANES);
    // OPT-D: lane_s is 32*T doubles (11+ GiB at 47M refs). Value-init would zero-fill it
    // on THIS thread, first-touching every page onto one NUMA node and paying a serial
    // 11 GiB memset again on every e_step_iter call. Allocate WITHOUT value-init (new[]
    // default-inits doubles = no fill); each lane zeroes (and first-touches) its own
    // T-slice inside the parallel E-step region. Stride is exactly T — layout is
    // byte-for-byte what it was, so every l*T+t index is unchanged.
    std::unique_ptr<double[]> lane_s_buf(new double[static_cast<size_t>(COMPACT_EM_LANES) * T]);
    double* const lane_s = lane_s_buf.get();
    // gn_s: written by lane l over r in [lane_lo[l], lane_lo[l+1]) and read ascending in
    // m_step, so first-touch on the same partition. Zeroed ONCE: deg==0 reads never assign
    // gn_s[r] and must stay 0.0 across all E-steps. The CLI passes null_lp whenever a null
    // model is loaded (one double per read gathered while the tape is written).
    std::unique_ptr<double[]> gn_s_buf(null_on ? new double[R] : nullptr);
    double* const gn_s = gn_s_buf.get();
    if (null_on) {
        #pragma omp parallel for schedule(static)
        for (int l = 0; l < COMPACT_EM_LANES; ++l)
            std::fill(gn_s + lane_lo[l], gn_s + lane_lo[l + 1], 0.0);
    }
    // OPT-B: log(max(w_t, min_weight)) hoisted out of the O(A) inner gather; refilled once
    // per E-step by a static-scheduled parallel loop (stable first-touch across NUMA nodes).
    std::vector<double> logw(T);

    // init weights (+ null coordinate) — identical to compact_em_map
    std::vector<double> weights(Tn, 1.0 / static_cast<double>(T));
    if (null_on) {
        double w0n = std::clamp(std::exp(params.log_w_null), params.min_weight, 1.0 - 1e-6);
        for (uint32_t t = 0; t < T; ++t) weights[t] *= (1.0 - w0n);
        weights[T] = w0n;
    }

    // Parallel streaming E-step (iterations): transient per-ref lane_s (+ gnull), no gamma.
    auto e_step_iter = [&](const double* w) -> double {
        const double log_w_null_dyn = null_on ? std::log(std::max(w[T], params.min_weight)) : 0.0;
        // OPT-B: fill the hoisted log-weight table (schedule(static) => stable first-touch).
        #pragma omp parallel for schedule(static)
        for (int64_t t = 0; t < static_cast<int64_t>(T); ++t)
            logw[t] = std::log(std::max(w[t], params.min_weight));
        std::vector<double> ll_lane(COMPACT_EM_LANES, 0.0);
        #pragma omp parallel for schedule(dynamic, 1)
        for (int l = 0; l < COMPACT_EM_LANES; ++l) {
            TapeStream<uint32_t> ts_ref;
            TapeStream<double> ts_ba, ts_bm;
            if (!ts_ref.open(p_ref, A) || !ts_ba.open(p_basea, A) || !ts_bm.open(p_basem, A)) std::abort();
            double* s_t = lane_s + static_cast<size_t>(l) * T;
            std::fill_n(s_t, T, 0.0);   // OPT-D: parallel first-touch (was the serial std::fill)
            double ll_local = 0.0;
            std::vector<double> ls;
            for (uint32_t r = lane_lo[l]; r < lane_lo[l + 1]; ++r) {
                const uint64_t s = read_offsets[r], e = read_offsets[r + 1];
                const uint32_t deg = static_cast<uint32_t>(e - s);
                if (deg == 0) continue;
                ts_ref.ensure(s, e); ts_ba.ensure(s, e); ts_bm.ensure(s, e);
                const uint32_t ncol = null_on ? 2 * deg + 1 : 2 * deg;
                ls.resize(ncol);
                for (uint32_t j = 0; j < deg; ++j) {
                    const uint32_t pf = std::min(j + kEmPfDist, deg - 1);   // OPT-C
                    __builtin_prefetch(&logw[ts_ref.at(s + pf)], 0, 3);
                    const double lw = logw[ts_ref.at(s + j)];   // OPT-B
                    ls[2 * j] = lw + ts_ba.at(s + j);
                    ls[2 * j + 1] = lw + ts_bm.at(s + j);
                }
                if (null_on)
                    ls[2 * deg] = log_w_null_dyn + params.lambda_null * null_lp[r];
                const double lse = log_sum_exp(ls.data(), ncol);
                ll_local += lse;
                for (uint32_t j = 0; j < deg; ++j) {
                    const uint32_t pf = std::min(j + kEmPfDist, deg - 1);   // OPT-C
                    __builtin_prefetch(&s_t[ts_ref.at(s + pf)], 1, 3);
                    const double da = ls[2 * j] - lse, dm = ls[2 * j + 1] - lse;
                    const double ga = (da > -700.0) ? std::exp(da) : 0.0;
                    const double gm = (dm > -700.0) ? std::exp(dm) : 0.0;
                    s_t[ts_ref.at(s + j)] += ga + gm;   // read-order accumulate within lane
                }
                if (null_on) {
                    const double d = ls[2 * deg] - lse;
                    gn_s[r] = (d > -700.0) ? std::exp(d) : 0.0;
                }
            }
            ll_lane[l] = ll_local;
            ts_ref.close(); ts_ba.close(); ts_bm.close();
        }
        double ll = 0.0; for (int l = 0; l < COMPACT_EM_LANES; ++l) ll += ll_lane[l];
        return ll;
    };
    // M-step: ascending-lane reduction — identical to compact_em_map's m_step_fused.
    // OPT-D: parallelism is over t ONLY; the per-t reduction over lanes stays serial and
    // ascending (that ordering is the bit-identity guarantee). `tot` is a separate serial
    // ascending-t sum of the exact values stored into out_w, so it is the same left-to-right
    // sum it was when fused into the t loop. The null reduction stays serial ascending r.
    auto m_step = [&](double* out_w) {
        const double prior_offset = params.alpha_prior - 1.0;
        #pragma omp parallel for schedule(static)
        for (int64_t t = 0; t < static_cast<int64_t>(T); ++t) {
            double w = prior_offset;
            for (int l = 0; l < COMPACT_EM_LANES; ++l) w += lane_s[static_cast<size_t>(l) * T + t];
            out_w[t] = std::max(w, params.min_weight);
        }
        double tot = 0.0;
        for (uint32_t t = 0; t < T; ++t) tot += out_w[t];        // ascending t, SERIAL
        double wn = 0.0;
        if (null_on) {
            wn = prior_offset;
            for (uint32_t r = 0; r < R; ++r) wn += gn_s[r];      // ascending r, SERIAL
            wn = std::max(wn, params.min_weight); tot += wn;
        }
        if (tot > 0.0) {
            #pragma omp parallel for schedule(static)
            for (int64_t t = 0; t < static_cast<int64_t>(T); ++t) out_w[t] /= tot;
            if (null_on) out_w[T] = wn / tot;
        }
    };
    auto do_e = [&](const double* w) { return e_step_iter(w); };
    auto do_m = [&](double* out) { m_step(out); };

    // ---- OPT-G: resume from a matching EM fixed-point checkpoint ----
    const uint64_t ckpt_fp = checkpoint_path.empty()
        ? 0ULL : tape_ckpt_fingerprint(input_fingerprint, R, T, A, params, null_on);
    bool resumed = false;
    if (resume && !checkpoint_path.empty()) {
        FILE* f = std::fopen(checkpoint_path.c_str(), "rb");
        if (!f) {
            std::fprintf(stderr, "warning: EM checkpoint %s not readable (%s); running from scratch\n",
                         checkpoint_path.c_str(), std::strerror(errno));
        } else {
            TapeCkptHeader h{};
            std::vector<double> wtmp(Tn);
            if (std::fread(&h, sizeof(h), 1, f) == 1 &&
                h.magic == TAPE_CKPT_MAGIC && h.version == TAPE_CKPT_VERSION &&
                h.fingerprint == ckpt_fp && h.n_weights == Tn &&
                std::fread(wtmp.data(), sizeof(double), Tn, f) == Tn) {
                weights.swap(wtmp);
                res.iterations = h.iterations;
                res.log_likelihood = h.log_likelihood;
                resumed = true;
                std::fprintf(stderr, "EM: resumed fixed point from %s (%u iterations)\n",
                             checkpoint_path.c_str(), h.iterations);
            } else {
                std::fprintf(stderr, "warning: EM checkpoint %s stale or corrupt; running from scratch\n",
                             checkpoint_path.c_str());
            }
            std::fclose(f);
        }
    }

    // ---- SQUAREM outer loop — structurally identical to compact_em_map ----
    std::vector<double> w0(Tn), w1(Tn), w2(Tn), r_buf(Tn), v_buf(Tn), prev_w(Tn);
    double prev_obj = -std::numeric_limits<double>::infinity();
    double step_min = params.squarem_step_min0, step_max0 = params.squarem_step_max0;
    if (step_max0 < step_min) std::swap(step_max0, step_min);
    double step_max = step_max0;
    const double mstep = std::max(params.squarem_mstep, 1.0000001);
    const double param_tol = std::sqrt(std::max(params.tol, 1e-15));

    // OPT-A: true iff lane_s/gn_s currently hold E(weights) for the CURRENT `weights`.
    // Every SQUAREM iteration ends with an E-scan taken at exactly the weight vector it
    // commits to `weights` (accepted branch: scan at w0, commit w0; squeeze branch: commit
    // w2, scan at w2), so the next iteration's opening E-scan is a recomputation of an
    // unchanged pure function and can be skipped. False on entry (iter 0 / warmup pay the
    // full scans); the plain-EM branch never sets it, so it stays false there.
    bool stats_at_weights = false;

    if (!resumed)
    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        double ll, objective;
        std::copy(weights.begin(), weights.end(), prev_w.begin());
        const bool sq = params.use_squarem && (iter >= params.squarem_warmup_iters);
        if (sq) {
            std::copy(weights.begin(), weights.end(), w0.begin());
            // OPT-A: w0 is a bit-copy of `weights`; when stats_at_weights holds, lane_s/gn_s
            // already ARE E(w0), so this scan is skipped (its return was already discarded).
            if (!stats_at_weights) do_e(w0.data());
            do_m(w1.data());
            do_e(w1.data()); do_m(w2.data());
            double sr2 = 0, sv2 = 0, srv = 0;
            squarem_build_directions(w0.data(), w1.data(), w2.data(), r_buf.data(), v_buf.data(), Tn, sr2, sv2, srv);
            double alpha = std::clamp(squarem_step_length(params.squarem_method, sr2, sv2, srv), step_min, step_max);
            squarem_extrapolate(w0.data(), r_buf.data(), v_buf.data(), Tn, alpha, params.min_weight);
            do_e(w0.data()); do_m(w1.data()); std::copy(w1.begin(), w1.end(), w0.begin());
            ll = do_e(w0.data());
            objective = em_objective_from_ll(params, w0.data(), T, ll);
            const bool nonmono = (iter > 0) && std::isfinite(params.squarem_objfn_inc)
                                 && (objective < (prev_obj - params.squarem_objfn_inc));
            if (nonmono) {
                std::copy(w2.begin(), w2.end(), weights.begin());
                ll = do_e(weights.data());
                objective = em_objective_from_ll(params, weights.data(), T, ll);
                if (std::abs(alpha - step_max) <= (1e-12 * (1.0 + std::abs(step_max)))) step_max = std::max(step_max0, step_max / mstep);
            } else {
                std::copy(w0.begin(), w0.end(), weights.begin());
            }
            // OPT-A: both branches leave lane_s/gn_s equal to E(weights) — accepted: last
            // scan do_e(w0), weights<-w0; squeeze: weights<-w2, last scan do_e(weights).
            stats_at_weights = true;
            if (std::abs(alpha - step_max) <= (1e-12 * (1.0 + std::abs(step_max)))) step_max = mstep * step_max;
            if (step_min < 0.0 && std::abs(alpha - step_min) <= (1e-12 * (1.0 + std::abs(step_min)))) step_min = mstep * step_min;
        } else {
            ll = do_e(weights.data());
            objective = em_objective_from_ll(params, weights.data(), T, ll);
            do_m(weights.data());
        }
        res.log_likelihood = ll; res.iterations = iter + 1;
        double rel = std::numeric_limits<double>::infinity();
        if (iter > 0) rel = std::abs(objective - prev_obj) / (std::abs(prev_obj) + 1e-15);
        double pcs = 0.0; for (uint32_t t = 0; t < T; ++t) { double d = weights[t] - prev_w[t]; pcs += d * d; }
        if (iter > 0 && (iter + 1) >= params.min_iters && rel <= params.tol && std::sqrt(pcs) <= param_tol) break;
        prev_obj = objective;
    }

    // ---- OPT-G: persist the converged fixed point (iterations only; the final pass
    // always re-runs). Skipped when we resumed (weights already came from disk). ----
    if (!checkpoint_path.empty() && !resumed) {
        TapeCkptHeader h{};
        h.magic = TAPE_CKPT_MAGIC;
        h.version = TAPE_CKPT_VERSION;
        h.n_weights = Tn;
        h.fingerprint = ckpt_fp;
        h.log_likelihood = res.log_likelihood;
        h.iterations = res.iterations;
        h.pad = 0;
        FILE* f = std::fopen(checkpoint_path.c_str(), "wb");
        if (!f) {
            std::fprintf(stderr, "warning: cannot write EM checkpoint %s: %s\n",
                         checkpoint_path.c_str(), std::strerror(errno));
        } else {
            std::fwrite(&h, sizeof(h), 1, f);
            std::fwrite(weights.data(), sizeof(double), Tn, f);
            std::fclose(f);
        }
    }

    // ---- OPT-F: final pass at frozen weights, PARALLEL over the same 32 fixed lanes the
    // iterations use: each lane owns a contiguous read range, opens its own TapeStream trio,
    // computes per-read transient gamma, and calls finalize_cb from the worker thread.
    // finalize_cb writes only per-r slots, so its call ORDER cannot affect output. The three
    // scalars use per-lane partials folded ascending — reproducible across thread counts, but
    // NOT bit-equal to the old serial running sum (they only reach a stderr line / dead
    // fields / a pi path that is disabled whenever finalize_cb is supplied). Per-read gamma
    // is bit-identical. Recomputing at the converged weights reproduces compact_em_map's
    // gamma bit-for-bit. finalize_cb is invoked CONCURRENTLY: it must be thread-safe. ----
    {
        const double log_w_null_dyn = null_on ? std::log(std::max(weights[T], params.min_weight)) : 0.0;
        // OPT-B.4: refill the shared log-weight table at the frozen weights before the lane
        // region — one std::log per ref instead of one per alignment, same value.
        #pragma omp parallel for schedule(static)
        for (int64_t t = 0; t < static_cast<int64_t>(T); ++t)
            logw[t] = std::log(std::max(weights[t], params.min_weight));

        // gamma_ancient is delivered when the joint damage E-step ran (use_damage -> ga)
        // OR when the CLI requested a post-hoc ancient responsibility (posthoc_gamma_ancient
        // -> gamma * p_anc with p_anc = sigmoid(base_a - base_m); bc cancels). Mirrors
        // compact_em_map's final pass so the streaming tape and in-RAM paths agree bitwise.
        const bool want_anc = params.use_damage || params.posthoc_gamma_ancient;
        std::vector<double> ll_lane(COMPACT_EM_LANES, 0.0);
        std::vector<double> sa_lane(COMPACT_EM_LANES, 0.0);
        std::vector<double> st_lane(COMPACT_EM_LANES, 0.0);

        #pragma omp parallel for schedule(dynamic, 1)
        for (int l = 0; l < COMPACT_EM_LANES; ++l) {
            TapeStream<uint32_t> ts_ref;
            TapeStream<double> ts_ba, ts_bm;
            if (!ts_ref.open(p_ref, A) || !ts_ba.open(p_basea, A) || !ts_bm.open(p_basem, A)) std::abort();
            std::vector<double> ls, gbuf, gabuf;
            std::vector<uint32_t> rbuf;
            double ll_local = 0.0, sa_local = 0.0, st_local = 0.0;
            for (uint32_t r = lane_lo[l]; r < lane_lo[l + 1]; ++r) {
                const uint64_t s = read_offsets[r], e = read_offsets[r + 1];
                const uint32_t deg = static_cast<uint32_t>(e - s);
                if (gnull_out) gnull_out[r] = 0.0;
                if (deg == 0) { if (finalize_cb) finalize_cb(r, 0, nullptr, nullptr, nullptr); continue; }
                ts_ref.ensure(s, e); ts_ba.ensure(s, e); ts_bm.ensure(s, e);
                const uint32_t ncol = null_on ? 2 * deg + 1 : 2 * deg;
                ls.resize(ncol); gbuf.resize(deg); rbuf.resize(deg);
                if (want_anc) gabuf.resize(deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    rbuf[j] = ts_ref.at(s + j);
                    const double lw = logw[rbuf[j]];
                    ls[2 * j] = lw + ts_ba.at(s + j);
                    ls[2 * j + 1] = lw + ts_bm.at(s + j);
                }
                if (null_on) ls[2 * deg] = log_w_null_dyn + params.lambda_null * null_lp[r];
                const double lse = log_sum_exp(ls.data(), ncol);
                ll_local += lse;
                for (uint32_t j = 0; j < deg; ++j) {
                    const double da = ls[2 * j] - lse, dm = ls[2 * j + 1] - lse;
                    const double ga = (da > -700.0) ? std::exp(da) : 0.0;
                    const double gm = (dm > -700.0) ? std::exp(dm) : 0.0;
                    gbuf[j] = ga + gm;
                    if (params.use_damage) gabuf[j] = ga;
                    else if (params.posthoc_gamma_ancient) {
                        const double p_anc = 1.0 / (1.0 + std::exp(ts_bm.at(s + j) - ts_ba.at(s + j)));
                        gabuf[j] = (ga + gm) * p_anc;
                    }
                    sa_local += ga; st_local += ga + gm;
                }
                if (null_on && gnull_out) {
                    const double d = ls[2 * deg] - lse;
                    gnull_out[r] = (d > -700.0) ? std::exp(d) : 0.0;
                }
                if (finalize_cb)
                    finalize_cb(r, deg, rbuf.data(), gbuf.data(), want_anc ? gabuf.data() : nullptr);
            }
            ll_lane[l] = ll_local; sa_lane[l] = sa_local; st_lane[l] = st_local;
            ts_ref.close(); ts_ba.close(); ts_bm.close();
        }

        double final_ll = 0.0, sum_a = 0.0, sum_t = 0.0;
        for (int l = 0; l < COMPACT_EM_LANES; ++l) {
            final_ll += ll_lane[l]; sum_a += sa_lane[l]; sum_t += st_lane[l];
        }
        res.log_likelihood = final_ll;
        res.gamma_anc_sum = sum_a; res.gamma_sum = sum_t;
        if (params.use_damage && sum_t > 0.0) res.pi = std::clamp(sum_a / sum_t, 0.01, 0.99);
    }

    if (null_on) { res.w_null = weights[Tn - 1]; weights.resize(T); }
    res.weights = std::move(weights);
    return res;
}

// Gate/small-fixture wrapper: writes the SoA tape from an in-RAM AlignmentData, runs
// the shared streaming core, and materializes full gamma into an EMState (so
// tests/gate_compact_em can compare bit-for-bit against compact_em_map). Production
// feeds tape_em_core directly from the streaming ingestion — no AlignmentData.
EMState compact_em_map_tape(const AlignmentData& data, const EMParams& params,
                            const std::string& scratch_dir) {
    EMState state;
    const uint32_t T = data.num_refs;
    const uint64_t A = data.alignments.size();
    const uint32_t R = data.num_reads;
    const bool null_on = params.use_null && !params.null_mu.empty();
    if (T == 0) { state.log_likelihood = 0.0; state.iterations = 0; return state; }

    const std::string p_ref = scratch_dir + "/tape.ref";
    const std::string p_ba  = scratch_dir + "/tape.basea";
    const std::string p_bm  = scratch_dir + "/tape.basem";
    {
        std::vector<uint32_t> csr_ref(A);
        std::vector<double> base_a(A), base_m(A);
        #pragma omp parallel for schedule(static)
        for (uint32_t r = 0; r < R; ++r)
            tape_decode_read(data, params, r, csr_ref.data(), base_a.data(), base_m.data());
        bool ok = tape_write_array(p_ref, csr_ref.data(), A * sizeof(uint32_t))
               && tape_write_array(p_ba,  base_a.data(),  A * sizeof(double))
               && tape_write_array(p_bm,  base_m.data(),  A * sizeof(double));
        if (!ok) { std::cerr << "compact_em_map_tape: tape write failed\n"; state.iterations = 0; return state; }
    }

    std::vector<uint64_t> ro64(R + 1);
    for (uint32_t r = 0; r <= R; ++r) ro64[r] = data.read_offsets[r];

    std::vector<double> null_lp;
    if (null_on) {
        null_lp.assign(R, 0.0);
        #pragma omp parallel for schedule(static)
        for (uint32_t r = 0; r < R; ++r) {
            const uint32_t s = data.read_offsets[r], e = data.read_offsets[r + 1];
            if (e > s) null_lp[r] = read_null_log_density(params, data.alignments.data() + s, e - s);
        }
    }

    state.gamma.assign(A, 0.0);
    if (params.use_damage || params.posthoc_gamma_ancient) state.gamma_ancient.assign(A, 0.0);
    if (null_on) state.gamma_null.assign(R, 0.0);
    auto cb = [&](uint32_t r, uint32_t deg, const uint32_t*, const double* g, const double* ganc) {
        const uint64_t base = ro64[r];
        for (uint32_t j = 0; j < deg; ++j) {
            state.gamma[base + j] = g[j];
            if (ganc) state.gamma_ancient[base + j] = ganc[j];
        }
    };
    TapeEMResult res = tape_em_core(p_ref, p_ba, p_bm, ro64.data(), R, T, A, params,
                                    null_on ? null_lp.data() : nullptr,
                                    null_on ? state.gamma_null.data() : nullptr, cb);

    state.weights = std::move(res.weights);
    state.log_likelihood = res.log_likelihood;
    state.iterations = res.iterations;
    state.pi = res.pi;
    state.w_null = res.w_null;

    ::unlink(p_ref.c_str()); ::unlink(p_ba.c_str()); ::unlink(p_bm.c_str());
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
#include "dart/emi_v7_io.hpp"

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
    const uint16_t* aln_len,  // Alignment length per alignment (for coherent-LR densities)
    const double* log_weights,
    uint32_t num_refs,
    const dart::EMParams& params,
    double pi_ancient,
    bool null_on,             // Null column active (calib/null model wired)
    double log_w_null_dyn,    // Current log w_null for the null column
    std::vector<std::vector<double>>& thread_ref_sums,
    std::vector<std::vector<double>>& thread_ref_sums_ancient,
    std::vector<double>& thread_null_sums,   // Per-thread Σ gamma_null accumulator
    int num_threads
) {
    const double eps = 1e-12;
    const double pi_a = std::clamp(pi_ancient, eps, 1.0 - eps);
    double ll = 0.0;

    // Check pointers
    if (!read_idx || !ref_idx || !bit_score || !damage_score || !log_weights) {
        return 0.0;
    }
    if (num_rows == 0 || num_threads <= 0 || thread_ref_sums.empty()) {
        return 0.0;
    }

    // Group alignments by read boundaries so each OpenMP iteration processes
    // one complete read (all its alignments stay atomic).
    std::vector<uint32_t> read_starts;
    read_starts.reserve(num_rows);
    for (uint32_t row = 0; row < num_rows;) {
        read_starts.push_back(row);
        const uint32_t current_read = read_idx[row];
        while (row < num_rows && read_idx[row] == current_read) {
            ++row;
        }
    }
    if (read_starts.empty()) {
        return 0.0;
    }

    // For alignment damage likelihood mode, we need dmg_ll_a/dmg_ll_m
    const bool use_aln_dmg = params.use_alignment_damage_likelihood && dmg_ll_a && dmg_ll_m;

    const int read_count = static_cast<int>(read_starts.size());

#ifdef _OPENMP
#pragma omp parallel num_threads(num_threads) reduction(+:ll)
    {
        std::vector<double> log_scores;
        const int tid = omp_get_thread_num();
        auto& local_ref_sums = thread_ref_sums[static_cast<size_t>(tid)];
        auto& local_ref_sums_ancient = thread_ref_sums_ancient[static_cast<size_t>(tid)];
        double& local_null_sum = thread_null_sums[static_cast<size_t>(tid)];

#pragma omp for schedule(static)
        for (int read_i = 0; read_i < read_count; ++read_i) {
            const uint32_t read_start = read_starts[static_cast<size_t>(read_i)];
            const uint32_t read_end =
                (read_i + 1 < read_count)
                    ? read_starts[static_cast<size_t>(read_i + 1)]
                    : num_rows;
            const uint32_t deg = read_end - read_start;
            if (deg == 0) continue;

            if (params.use_damage) {
                const uint32_t n_cols = null_on ? 2 * deg + 1 : 2 * deg;
                log_scores.resize(n_cols);
                for (uint32_t j = 0; j < deg; ++j) {
                    const uint32_t idx = read_start + j;
                    double log_w = log_weights[ref_idx[idx]];
                    if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                        log_w -= half_log_tlen(tlen[idx]);
                    }
                    const double score_contrib = real_log_density(
                        params, static_cast<double>(bit_score[idx]),
                        aln_len ? aln_len[idx] : 0, 1.0);

                    if (use_aln_dmg) {
                        const double p_read = std::clamp(static_cast<double>(damage_score[idx]), eps, 1.0 - eps);
                        log_scores[2 * j]     = log_w + score_contrib + std::log(p_read) + dmg_ll_a[idx];
                        log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - p_read) + dmg_ll_m[idx];
                    } else {
                        const double ds = std::clamp(static_cast<double>(damage_score[idx]), eps, 1.0 - eps);
                        log_scores[2 * j]     = log_w + score_contrib + std::log(pi_a) + std::log(ds);
                        log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - pi_a) + std::log(1.0 - ds);
                    }
                }
                if (null_on) {
                    const double lp = read_null_log_density_cols(
                        params, bit_score, aln_len, read_start, deg);
                    log_scores[2 * deg] = log_w_null_dyn + params.lambda_null * lp;
                }

                const double lse = dart::log_sum_exp(log_scores.data(), n_cols);
                ll += lse;

                for (uint32_t j = 0; j < deg; ++j) {
                    const uint32_t idx = read_start + j;
                    const double diff_a = log_scores[2 * j] - lse;
                    const double diff_m = log_scores[2 * j + 1] - lse;
                    const double g_ancient = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
                    const double g_modern  = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;
                    const double g_total = g_ancient + g_modern;

                    const uint32_t ri = ref_idx[idx];
                    if (ri >= num_refs) continue;
                    local_ref_sums[ri] += g_total;
                    local_ref_sums_ancient[ri] += g_ancient;
                }
                if (null_on) {
                    const double dn = log_scores[2 * deg] - lse;
                    local_null_sum += (dn > -700.0) ? std::exp(dn) : 0.0;
                }
            } else {
                const uint32_t n_cols = null_on ? deg + 1 : deg;
                log_scores.resize(n_cols);
                for (uint32_t j = 0; j < deg; ++j) {
                    const uint32_t idx = read_start + j;
                    double log_w = log_weights[ref_idx[idx]];
                    if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                        log_w -= half_log_tlen(tlen[idx]);
                    }
                    log_scores[j] = log_w + real_log_density(
                        params, static_cast<double>(bit_score[idx]),
                        aln_len ? aln_len[idx] : 0, 1.0);
                }
                if (null_on) {
                    const double lp = read_null_log_density_cols(
                        params, bit_score, aln_len, read_start, deg);
                    log_scores[deg] = log_w_null_dyn + params.lambda_null * lp;
                }

                const double lse = dart::log_sum_exp(log_scores.data(), n_cols);
                ll += lse;

                for (uint32_t j = 0; j < deg; ++j) {
                    const uint32_t idx = read_start + j;
                    const double diff = log_scores[j] - lse;
                    const double g_val = (diff > -700.0) ? std::exp(diff) : 0.0;
                    const uint32_t ri = ref_idx[idx];
                    if (ri >= num_refs) continue;
                    local_ref_sums[ri] += g_val;
                }
                if (null_on) {
                    const double dn = log_scores[deg] - lse;
                    local_null_sum += (dn > -700.0) ? std::exp(dn) : 0.0;
                }
            }
        }
    }
#else
    auto& local_ref_sums = thread_ref_sums[0];
    auto& local_ref_sums_ancient = thread_ref_sums_ancient[0];
    double& local_null_sum = thread_null_sums[0];
    std::vector<double> log_scores;
    for (int read_i = 0; read_i < read_count; ++read_i) {
        const uint32_t read_start = read_starts[static_cast<size_t>(read_i)];
        const uint32_t read_end =
            (read_i + 1 < read_count)
                ? read_starts[static_cast<size_t>(read_i + 1)]
                : num_rows;
        const uint32_t deg = read_end - read_start;
        if (deg == 0) continue;

        if (params.use_damage) {
            const uint32_t n_cols = null_on ? 2 * deg + 1 : 2 * deg;
            log_scores.resize(n_cols);
            for (uint32_t j = 0; j < deg; ++j) {
                const uint32_t idx = read_start + j;
                double log_w = log_weights[ref_idx[idx]];
                if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                    log_w -= half_log_tlen(tlen[idx]);
                }
                const double score_contrib = real_log_density(
                    params, static_cast<double>(bit_score[idx]),
                    aln_len ? aln_len[idx] : 0, 1.0);

                if (use_aln_dmg) {
                    const double p_read = std::clamp(static_cast<double>(damage_score[idx]), eps, 1.0 - eps);
                    log_scores[2 * j]     = log_w + score_contrib + std::log(p_read) + dmg_ll_a[idx];
                    log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - p_read) + dmg_ll_m[idx];
                } else {
                    const double ds = std::clamp(static_cast<double>(damage_score[idx]), eps, 1.0 - eps);
                    log_scores[2 * j]     = log_w + score_contrib + std::log(pi_a) + std::log(ds);
                    log_scores[2 * j + 1] = log_w + score_contrib + std::log(1.0 - pi_a) + std::log(1.0 - ds);
                }
            }
            if (null_on) {
                const double lp = read_null_log_density_cols(
                    params, bit_score, aln_len, read_start, deg);
                log_scores[2 * deg] = log_w_null_dyn + params.lambda_null * lp;
            }

            const double lse = dart::log_sum_exp(log_scores.data(), n_cols);
            ll += lse;

            for (uint32_t j = 0; j < deg; ++j) {
                const uint32_t idx = read_start + j;
                const double diff_a = log_scores[2 * j] - lse;
                const double diff_m = log_scores[2 * j + 1] - lse;
                const double g_ancient = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
                const double g_modern  = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;
                const double g_total = g_ancient + g_modern;

                const uint32_t ri = ref_idx[idx];
                if (ri >= num_refs) continue;
                local_ref_sums[ri] += g_total;
                local_ref_sums_ancient[ri] += g_ancient;
            }
            if (null_on) {
                const double dn = log_scores[2 * deg] - lse;
                local_null_sum += (dn > -700.0) ? std::exp(dn) : 0.0;
            }
        } else {
            const uint32_t n_cols = null_on ? deg + 1 : deg;
            log_scores.resize(n_cols);
            for (uint32_t j = 0; j < deg; ++j) {
                const uint32_t idx = read_start + j;
                double log_w = log_weights[ref_idx[idx]];
                if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                    log_w -= half_log_tlen(tlen[idx]);
                }
                log_scores[j] = log_w + real_log_density(
                    params, static_cast<double>(bit_score[idx]),
                    aln_len ? aln_len[idx] : 0, 1.0);
            }
            if (null_on) {
                const double lp = read_null_log_density_cols(
                    params, bit_score, aln_len, read_start, deg);
                log_scores[deg] = log_w_null_dyn + params.lambda_null * lp;
            }

            const double lse = dart::log_sum_exp(log_scores.data(), n_cols);
            ll += lse;

            for (uint32_t j = 0; j < deg; ++j) {
                const uint32_t idx = read_start + j;
                const double diff = log_scores[j] - lse;
                const double g_val = (diff > -700.0) ? std::exp(diff) : 0.0;
                const uint32_t ri = ref_idx[idx];
                if (ri >= num_refs) continue;
                local_ref_sums[ri] += g_val;
            }
            if (null_on) {
                const double dn = log_scores[deg] - lse;
                local_null_sum += (dn > -700.0) ? std::exp(dn) : 0.0;
            }
        }
    }
#endif

    return ll;
}

// M-step from accumulated ref_sums
void streaming_m_step(
    const std::vector<double>& ref_sums,
    std::vector<double>& weights,
    const dart::EMParams& params,
    double null_mass,
    bool null_on,
    double& w_null_io
) {
    const uint32_t T = static_cast<uint32_t>(weights.size());
    const double prior_offset = params.alpha_prior - 1.0;

    double sum = 0.0;
    for (uint32_t t = 0; t < T; ++t) {
        weights[t] = std::max(ref_sums[t] + prior_offset, params.min_weight);
        sum += weights[t];
    }

    // Null coordinate: informative Dirichlet floor (alpha_null) mirrors the in-memory
    // m_step_parallel so streaming w_null co-evolves identically.
    double w_null = 0.0;
    if (null_on) {
        w_null = std::max(null_mass + params.alpha_null - 1.0, params.min_weight);
        sum += w_null;
    }

    double inv_sum = 1.0 / sum;
    for (uint32_t t = 0; t < T; ++t) {
        weights[t] *= inv_sum;
    }
    if (null_on) w_null_io = w_null * inv_sum;
}

} // anonymous namespace

namespace dart {

template <class ReaderT>
StreamingEMResult streaming_em(
    ReaderT& reader,
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

    // Null column: parallels the in-memory squarem_em null coordinate. When active,
    // seed w_null from exp(log_w_null) and shrink the real weights by (1 - w_null0)
    // so the simplex including the null still sums to 1.
    const bool null_on = params.use_null &&
        (!params.null_mu.empty() || params.use_calib_cont);
    if (null_on) {
        const double w_null0 = std::clamp(std::exp(params.log_w_null),
                                          params.min_weight, 1.0 - 1e-6);
        for (uint32_t t = 0; t < T; ++t) result.weights[t] *= (1.0 - w_null0);
        result.w_null = w_null0;
    }

    // Thread-local accumulators for per-row-group read-parallel processing.
    int num_threads = 1;
#ifdef _OPENMP
    num_threads = std::max(1, omp_get_max_threads());
#endif
    // Cap at hardware concurrency to prevent excessive memory use
    const int hw_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (hw_threads > 0) {
        num_threads = std::min(num_threads, hw_threads);
    }

    std::vector<std::vector<double>> thread_ref_sums(num_threads);
    std::vector<std::vector<double>> thread_ref_sums_ancient(num_threads);
    std::vector<double> thread_null_sums(num_threads, 0.0);
    for (int t = 0; t < num_threads; ++t) {
        thread_ref_sums[t].resize(T, 0.0);
        if (params.use_damage) {
            thread_ref_sums_ancient[t].resize(T, 0.0);
        }
    }

    double prev_obj = -std::numeric_limits<double>::infinity();

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
        std::fill(thread_null_sums.begin(), thread_null_sums.end(), 0.0);

        for (uint32_t t = 0; t < T; ++t) {
            log_weights[t] = std::log(std::max(result.weights[t], params.min_weight));
        }
        const double log_w_null_dyn = null_on
            ? std::log(std::max(result.w_null, params.min_weight)) : 0.0;

        double total_ll = 0.0;

        // E-step: process row groups sequentially; each row group parallelizes
        // internally over complete reads in streaming_e_step_rowgroup().
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
            const uint16_t* /*tstart*/,
            const uint16_t* /*tend*/,
            const uint16_t* tlen,
            const uint16_t* /*qlen*/,
            const uint16_t* /*mismatch*/,
            const uint16_t* /*gapopen*/
        ) {
            total_ll += streaming_e_step_rowgroup(
                num_rows, read_idx, ref_idx, bit_score, damage_score,
                dmg_ll_a, dmg_ll_m, tlen, aln_len,
                log_weights.data(), T, params, result.pi,
                null_on, log_w_null_dyn,
                thread_ref_sums, thread_ref_sums_ancient, thread_null_sums, num_threads
            );
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

        // M-step: update weights from accumulated sums (+ null coordinate)
        double total_null_mass = 0.0;
        if (null_on) {
            for (int t = 0; t < num_threads; ++t) total_null_mass += thread_null_sums[t];
        }
        streaming_m_step(ref_sums, result.weights, params,
                         total_null_mass, null_on, result.w_null);

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

        if (!std::isfinite(total_ll)) {
            break;
        }

        result.log_likelihood = total_ll;
        result.iterations = iter + 1;

        double total_obj = em_objective_from_ll(params, result.weights.data(), T, total_ll);

        // Progress callback
        if (progress_cb) {
            EMIterationDiagnostics diag;
            diag.iteration = iter + 1;
            diag.log_likelihood = total_ll;
            diag.objective = total_obj;
            diag.rel_change = std::abs(total_obj - prev_obj) / (std::abs(prev_obj) + 1e-15);
            fill_weight_diagnostics(diag, result.weights.data(), T, params.min_weight);
            progress_cb(diag);
        }

        // Convergence check
        if (iter > 0) {
            double rel_change = std::abs(total_obj - prev_obj) / (std::abs(prev_obj) + 1e-15);
            if (rel_change < params.tol && (iter + 1) >= params.min_iters) {
                break;
            }
        }
        prev_obj = total_obj;
    }

    return result;
}

template <class ReaderT>
void streaming_em_finalize(
    ReaderT& reader,
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
    // Null column against the FROZEN converged weights (matches the E-step's null col).
    const bool null_on = params.use_null &&
        (!params.null_mu.empty() || params.use_calib_cont);
    const double log_w_null_dyn = null_on
        ? std::log(std::max(result.w_null, params.min_weight)) : 0.0;

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
        const uint16_t* aln_len,
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
        // Emit gamma_ancient when the mapping itself is damage-aware (use_damage),
        // or when the mapping is damage-free but a post-hoc ancient responsibility
        // is requested against the frozen weights (posthoc_gamma_ancient). The
        // latter needs the per-alignment damage_score column to be present.
        const bool posthoc_anc = params.posthoc_gamma_ancient && damage_score;
        const bool emit_anc = params.use_damage || posthoc_anc;
        std::vector<float> g_ancient_vals;
        if (emit_anc) {
            g_ancient_vals.resize(num_rows);
        }
        std::vector<float> g_null_vals;
        if (null_on) {
            g_null_vals.resize(num_rows, 0.0f);
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
                const uint32_t n_cols = null_on ? 2 * deg + 1 : 2 * deg;
                log_scores.resize(n_cols);
                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double log_w = log_weights[ref_idx[idx]];
                    if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                        log_w -= half_log_tlen(tlen[idx]);
                    }
                    double score_contrib = real_log_density(
                        params, static_cast<double>(bit_score[idx]),
                        aln_len ? aln_len[idx] : 0, 1.0);

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
                if (null_on) {
                    double lp = read_null_log_density_cols(
                        params, bit_score, aln_len, read_start, deg);
                    log_scores[2 * deg] = log_w_null_dyn + params.lambda_null * lp;
                }

                double lse = log_sum_exp(log_scores.data(), n_cols);

                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double diff_a = log_scores[2 * j] - lse;
                    double diff_m = log_scores[2 * j + 1] - lse;
                    double g_ancient = (diff_a > -700.0) ? std::exp(diff_a) : 0.0;
                    double g_modern  = (diff_m > -700.0) ? std::exp(diff_m) : 0.0;

                    g_vals[idx] = static_cast<float>(g_ancient + g_modern);
                    g_ancient_vals[idx] = static_cast<float>(g_ancient);
                }
                if (null_on) {
                    double dn = log_scores[2 * deg] - lse;
                    float gnull = static_cast<float>((dn > -700.0) ? std::exp(dn) : 0.0);
                    for (uint32_t j = 0; j < deg; ++j) g_null_vals[read_start + j] = gnull;
                }
            } else {
                const uint32_t n_cols = null_on ? deg + 1 : deg;
                log_scores.resize(n_cols);
                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double log_w = log_weights[ref_idx[idx]];
                    if (params.normalize_by_length && tlen && tlen[idx] > 0) {
                        log_w -= half_log_tlen(tlen[idx]);
                    }
                    log_scores[j] = log_w + real_log_density(
                        params, static_cast<double>(bit_score[idx]),
                        aln_len ? aln_len[idx] : 0, 1.0);
                }
                if (null_on) {
                    double lp = read_null_log_density_cols(
                        params, bit_score, aln_len, read_start, deg);
                    log_scores[deg] = log_w_null_dyn + params.lambda_null * lp;
                }

                double lse = log_sum_exp(log_scores.data(), n_cols);

                for (uint32_t j = 0; j < deg; ++j) {
                    uint32_t idx = read_start + j;
                    double diff = log_scores[j] - lse;
                    g_vals[idx] = static_cast<float>((diff > -700.0) ? std::exp(diff) : 0.0);
                }
                if (null_on) {
                    double dn = log_scores[deg] - lse;
                    float gnull = static_cast<float>((dn > -700.0) ? std::exp(dn) : 0.0);
                    for (uint32_t j = 0; j < deg; ++j) g_null_vals[read_start + j] = gnull;
                }

                // Post-hoc ancient responsibility against the FROZEN weights: the
                // damage-free mapping gamma is left untouched; gamma_ancient is the
                // clean gamma scaled by the per-alignment P(ancient | evidence).
                // Damage never re-enters the mapping/M-step (that already converged).
                if (posthoc_anc) {
                    const bool use_aln_dmg = dmg_ll_a && dmg_ll_m;
                    for (uint32_t j = 0; j < deg; ++j) {
                        uint32_t idx = read_start + j;
                        double ds = std::clamp(static_cast<double>(damage_score[idx]), eps, 1.0 - eps);
                        double p_anc;
                        if (use_aln_dmg) {
                            double la = std::log(ds) + dmg_ll_a[idx];
                            double lm = std::log(1.0 - ds) + dmg_ll_m[idx];
                            double mx = std::max(la, lm);
                            double ea = std::exp(la - mx), em = std::exp(lm - mx);
                            p_anc = ea / (ea + em);
                        } else {
                            double num = pi_a * ds;
                            p_anc = num / (num + (1.0 - pi_a) * (1.0 - ds));
                        }
                        g_ancient_vals[idx] = static_cast<float>(g_vals[idx] * p_anc);
                    }
                }
            }
        }

        // Call user callback with computed gamma
        callback(
            rg_idx, num_rows, read_idx, ref_idx,
            g_vals.data(),
            emit_anc ? g_ancient_vals.data() : nullptr,
            null_on ? g_null_vals.data() : nullptr,
            bit_score, aln_len
        );
    });
}

// Explicit instantiations: v6 (ColumnarIndexReader) and v7 (emiv7::Reader).
template StreamingEMResult streaming_em<ColumnarIndexReader>(
    ColumnarIndexReader&, const EMParams&, const EMProgressCallback&);
template void streaming_em_finalize<ColumnarIndexReader>(
    ColumnarIndexReader&, const StreamingEMResult&, const EMParams&, StreamingGammaCallback);
template StreamingEMResult streaming_em<emiv7::Reader>(
    emiv7::Reader&, const EMParams&, const EMProgressCallback&);
template void streaming_em_finalize<emiv7::Reader>(
    emiv7::Reader&, const StreamingEMResult&, const EMParams&, StreamingGammaCallback);

} // namespace dart
