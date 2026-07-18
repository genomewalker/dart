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
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace dart {

// One length-bin's true-hit density f_true(bit): a 1- or 2-component Gaussian
// mixture over the bit scores of LABELED true-source alignments. w1 == 0 selects
// a single component. Bins share kNullBinEdges with the null Gumbel so f_true and
// f_0 describe the same statistic on the same partition.
struct TrueHitBin {
    float w0 = 1.0f, mu0 = 0.0f, sd0 = 1.0f;  // component 0 (weight, mean, sd)
    float w1 = 0.0f, mu1 = 0.0f, sd1 = 1.0f;  // component 1 (w1>0 → 2-component)
    uint32_t n = 0;                            // true-hit count used for the fit
};

// EM algorithm parameters
struct EMParams {
    // Raw-bit temperature for the LEGACY (non-calibrated) real component,
    // exp(lambda_b * bit). mmseqs bit scores are log2-odds; the principled
    // bits→nats conversion is ln2 ≈ 0.693, so that is the default. The old value
    // of 3.0 over-tempered ~4.3x and forced homolog near-ties to gamma≈1. When a
    // coherent-LR calibration model is loaded (use_calib), this knob is unused —
    // the real component becomes the normalized density log f_true (lambda ≡ 1).
    double lambda_b = M_LN2;     // 0.6931…  (was 3.0)
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
    // Damage-free mapping with post-hoc ancient responsibility. When true, the
    // mapping gamma is computed damage-free (plain E-step) but gamma_ancient is
    // still emitted as gamma * P(ancient | alignment) evaluated against the
    // FROZEN final weights — damage never feeds back into the M-step. Consumed by
    // streaming_em_finalize; inert unless set (zero behavior change when off).
    bool posthoc_gamma_ancient = false;
    double em_min_margin = 0.5;  // auth keep: require top1 - top2 responsibility >= this
    double em_max_neff = 2.0;    // auth keep: require exp(entropy(gamma)) <= this

    // Dirichlet prior parameter (MAP): alpha=1 gives no prior effect.
    // alpha>1 shrinks toward uniform; alpha<1 promotes sparsity.
    double alpha_prior = 1.0;
    bool normalize_by_length = false;  // Apply -0.5*log(ref_length) in E-step scores
                                       // Penalizes longer references in a proper objective

    // Per-query bit-score normalization (MMseqs2-style):
    // Divide each alignment's bit_score by the max bit_score for that read before applying
    // lambda_b. This makes competition between hits purely relative (all scores in [0,1]),
    // preventing high-complexity reads or long proteins from dominating.
    // Default: false (raw bit scores, existing behavior).
    bool normalize_scores_per_query = false;

    // Optional initial weights for warm-starting EM (e.g., from coverage-aware outer loop).
    // If non-empty, must have size == num_refs. EM initializes from these instead of uniform.
    std::vector<double> initial_weights;

    // Per-reference coverage log-penalty, added directly to log(w_t) in every E-step.
    // coverage_log_weights[t] = gamma_cov * log(coverage_weight_t), where coverage_weight_t
    // is in [w_min, 1.0] (computed by the outer coverage loop from KL divergence).
    // This brings coverage feedback inside the inner EM loop so every gamma computation
    // is penalized for patchy references, rather than only influencing the warm-start.
    // Empty (default) = disabled, no coverage penalty.
    std::vector<float> coverage_log_weights;

    // ---------------------------------------------------------------------
    // Null / noise component (spurious short-peptide rejection)
    //
    // Adds one extra "null" column to each read's E-step mixture. The null
    // component scores a read by an empirical per-aligned-length-bin Gumbel
    // fit from a reversed-sequence decoy search (NOT K-A / mmseqs E-values,
    // which are invalid at 5-6 residues). Real-reference weights then sum to
    // 1 - w_null, and a read dominated by the null is rejected downstream.
    //
    // All fields inert unless use_null == true → zero behavior change when off.
    // null_mu / null_beta are indexed by null length bin (see kNullBinEdges).
    // Missing/empty bins fall back to a pooled fit written by the fitter.
    std::vector<float> null_mu;    // per-bin Gumbel location (on lambda_b*bit scale)
    std::vector<float> null_beta;  // per-bin Gumbel scale
    double log_w_null = std::log(0.096);  // log soft abundance weight of the null column
    double alpha_null = 0.0;       // informative Dirichlet floor for the null coordinate (>1 = floor)
    double lambda_null = 1.0;      // tempering on the null log-density

    // ---------------------------------------------------------------------
    // Coherent likelihood-ratio calibration (f_true vs f_0)
    //
    // The legacy real component mixes an UNNORMALIZED exp(lambda_b*bit) with the
    // NORMALIZED null Gumbel, so the null underflows and the flat-DB rejection
    // guard is silently a no-op. The coherent-LR mode instead scores every EM
    // component with a NORMALIZED density of the SAME statistic (the per-candidate
    // bit score): real → log f_true(bit | Lbin), null → log f_0(bit | Lbin). f_0
    // reuses the existing per-bin Gumbel (null_mu/null_beta) on the RAW bit scale.
    // f_true is a per-length-bin 1-or-2-component Gaussian fitted ONCE from labeled
    // synthetic true-source hits and frozen into the calibration model file.
    //
    // Inert unless use_calib == true AND true_bins is non-empty → the master
    // switch preserves exact legacy behavior when off. lambda is implicitly 1 here;
    // calib_temperature is the only surviving knob, a frozen global misspecification
    // correction (1.0 = identity).
    std::vector<TrueHitBin> true_bins;  // per null-length-bin f_true params
    float  true_censor = 0.0f;     // reporting-threshold bit (provenance; see log_f_true)
    double calib_temperature = 1.0;  // frozen misspecification correction on the LR
    bool   use_calib = false;      // master switch; false → legacy raw-bit real term
    bool use_null = false;         // master switch; false → null column never added

    // v2 LENGTH-CONTINUOUS calibration (sol short-fragment fix). When use_calib_cont
    // is true these supersede the discrete true_bins / per-bin null Gumbel: both
    // f_true and f_0 carry length-continuous locations, so a short read is never
    // compared to a long-read pooled mean. Default-off; only a v2 calib model (loaded
    // via load_calib_into_params) sets use_calib_cont. See log_f_true_cont/log_f0_cont.
    bool   use_calib_cont = false;
    double ft_mu_a = 0.0, ft_mu_b = 0.0;   // f_true mean   mu(L)    = a + b*L
    double ft_ls_a = 0.0, ft_ls_b = 0.0;   // f_true log-sd sigma(L) = exp(a + b*L)
    double ft_nu   = 0.0;                   // Student-t dof (<=0 → Gaussian)
    double f0_loc_a = 0.0, f0_loc_b = 0.0;  // f_0 K-A Gumbel loc(L) = a + b*L
    double f0_beta  = 1.0;                  // f_0 K-A Gumbel scale (bit units)
};

// Aligned-length bins for the empirical null Gumbel. aln_len falls into bin b
// when kNullBinEdges[b] <= aln_len < kNullBinEdges[b+1]; the final open bin
// [25, inf) catches everything else. 6 bins → 7 edges (last = sentinel).
inline constexpr int kNumNullBins = 6;
inline constexpr uint16_t kNullBinEdges[kNumNullBins + 1] = {4, 7, 10, 13, 17, 25, 65535};

// Map an alignment length to its null bin index [0, kNumNullBins). Lengths
// below the first edge clamp to bin 0; lengths at/after the last edge → last bin.
inline int null_length_bin(uint16_t aln_len) {
    for (int b = kNumNullBins - 1; b >= 0; --b) {
        if (aln_len >= kNullBinEdges[b]) return b;
    }
    return 0;
}

// log p_null(b | L) for a Gumbel(mu, beta): the density
//   f(x) = (1/beta) exp(-(z + exp(-z))),  z = (x - mu)/beta
// so log f = -z - exp(-z) - log(beta). Guarded against beta <= 0.
inline double log_p_null_gumbel(double score, double mu, double beta) {
    if (!(beta > 0.0)) return -std::numeric_limits<double>::infinity();
    const double z = (score - mu) / beta;
    const double ez = (z < 700.0) ? std::exp(-z) : std::numeric_limits<double>::infinity();
    if (!std::isfinite(ez)) return -std::numeric_limits<double>::infinity();
    return -z - ez - std::log(beta);
}

// Forward decl: the v2 length-continuous K-A null density (defined below) is used by
// read_null_log_density, which precedes its definition.
inline double log_f0_cont(const EMParams& params, double bit, uint16_t aln_len);

// Best null log-density over a read's alignments: max_j log p_null(b_j | L_j),
// evaluated on the same tempered scale as the real scores (lambda_b * bit_score).
// Returns -inf when the model is empty (caller should treat as null disabled).
inline double read_null_log_density(
    const EMParams& params,
    const struct Alignment* alns,
    uint32_t deg);

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

inline double coverage_log_penalty(
    const EMParams& params,
    uint32_t ref_idx)
{
    if (params.coverage_log_weights.empty()) return 0.0;
    if (ref_idx >= params.coverage_log_weights.size()) return 0.0;
    return static_cast<double>(params.coverage_log_weights[ref_idx]);
}

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

// Best null log-density over a read's alignments (definition; declared above).
// Uses lambda_b * bit_score as the score (same tempered scale as real hits) and
// the per-bin Gumbel selected by null_length_bin(aln.aln_len). Falls back to the
// last (pooled) bin when a bin has no fitted parameters.
inline double read_null_log_density(
    const EMParams& params,
    const Alignment* alns,
    uint32_t deg)
{
    // v2 continuous K-A null needs no per-bin table; the discrete paths do.
    if (!params.use_calib_cont && (params.null_mu.empty() || params.null_beta.empty())) {
        return -std::numeric_limits<double>::infinity();
    }
    const int nbins = params.null_mu.empty() ? 0
        : static_cast<int>(std::min(params.null_mu.size(), params.null_beta.size()));
    double best = -std::numeric_limits<double>::infinity();
    for (uint32_t j = 0; j < deg; ++j) {
        const Alignment& a = alns[j];
        double lp;
        if (params.use_calib_cont) {
            // Length-continuous Karlin-Altschul null density on the RAW bit scale.
            lp = log_f0_cont(params, static_cast<double>(a.bit_score), a.aln_len);
        } else {
            int b = null_length_bin(a.aln_len);
            if (b >= nbins) b = nbins - 1;
            const double mu = static_cast<double>(params.null_mu[b]);
            const double beta = static_cast<double>(params.null_beta[b]);
            // Coherent-LR: f_0 lives on the RAW bit scale (same statistic as f_true);
            // legacy: on the lambda_b*bit tempered scale the null was fit against.
            const double score = params.use_calib
                ? static_cast<double>(a.bit_score)
                : params.lambda_b * static_cast<double>(a.bit_score);
            lp = log_p_null_gumbel(score, mu, beta);
        }
        if (lp > best) best = lp;
    }
    return best;
}

// Streaming twin of read_null_log_density: best null log-density over a read's
// alignments, sourced from raw columnar pointers (bit_score / aln_len) rather than
// Alignment structs. Formula is byte-identical to read_null_log_density so the
// streaming and in-memory null columns agree.
inline double read_null_log_density_cols(
    const EMParams& params,
    const float* bit_score,
    const uint16_t* aln_len,
    uint32_t start,
    uint32_t deg)
{
    if (!params.use_calib_cont && (params.null_mu.empty() || params.null_beta.empty())) {
        return -std::numeric_limits<double>::infinity();
    }
    const int nbins = params.null_mu.empty() ? 0
        : static_cast<int>(std::min(params.null_mu.size(), params.null_beta.size()));
    double best = -std::numeric_limits<double>::infinity();
    for (uint32_t j = 0; j < deg; ++j) {
        const uint32_t idx = start + j;
        const double bit = static_cast<double>(bit_score[idx]);
        const uint16_t al = aln_len ? aln_len[idx] : 0;
        double lp;
        if (params.use_calib_cont) {
            lp = log_f0_cont(params, bit, al);
        } else {
            int b = null_length_bin(al);
            if (b >= nbins) b = nbins - 1;
            const double mu = static_cast<double>(params.null_mu[b]);
            const double beta = static_cast<double>(params.null_beta[b]);
            const double score = params.use_calib ? bit : params.lambda_b * bit;
            lp = log_p_null_gumbel(score, mu, beta);
        }
        if (lp > best) best = lp;
    }
    return best;
}

// log N(x | mu, sd), guarded against sd <= 0.
inline double log_gaussian(double x, double mu, double sd) {
    if (!(sd > 0.0)) return -std::numeric_limits<double>::infinity();
    const double z = (x - mu) / sd;
    return -0.5 * z * z - std::log(sd) - 0.5 * std::log(2.0 * M_PI);
}

// log f_true(bit | Lbin): the normalized true-hit density for a candidate's bit
// score, a 1- or 2-component Gaussian selected by the alignment's length bin.
// Returns -inf when no calibration is loaded (caller must gate on use_calib).
// The reporting-threshold censoring is left un-applied: f_true is fit on observed
// (already above-threshold) true-hit bits, and the same untruncated convention is
// used for f_0, so the ratio log f_true - log f_0 is a consistent density ratio.
// ceiling: no left-truncation renormalization at true_censor; the residual
// misspecification is absorbed by calib_temperature. upgrade: fit/normalize both
// f_true and f_0 as densities left-truncated at true_censor.
inline double log_f_true(const EMParams& params, double bit, uint16_t aln_len) {
    if (params.true_bins.empty()) return -std::numeric_limits<double>::infinity();
    int b = null_length_bin(aln_len);
    const int nbins = static_cast<int>(params.true_bins.size());
    if (b >= nbins) b = nbins - 1;
    const TrueHitBin& t = params.true_bins[b];
    const double w0 = std::max(0.0, static_cast<double>(t.w0));
    const double l0 = std::log(std::max(w0, 1e-300)) + log_gaussian(bit, t.mu0, t.sd0);
    if (t.w1 <= 0.0f) return l0;
    const double l1 = std::log(static_cast<double>(t.w1)) + log_gaussian(bit, t.mu1, t.sd1);
    const double m = std::max(l0, l1);
    if (!std::isfinite(m)) return m;
    return m + std::log(std::exp(l0 - m) + std::exp(l1 - m));
}

// log Student-t(x | mu, sigma, nu). nu<=0 or sigma<=0 → -inf (caller falls back).
inline double log_student_t(double x, double mu, double sigma, double nu) {
    if (!(sigma > 0.0) || !(nu > 0.0)) return -std::numeric_limits<double>::infinity();
    const double z = (x - mu) / sigma;
    return std::lgamma(0.5 * (nu + 1.0)) - std::lgamma(0.5 * nu)
         - 0.5 * std::log(nu * M_PI) - std::log(sigma)
         - 0.5 * (nu + 1.0) * std::log1p(z * z / nu);
}

// v2 LENGTH-CONTINUOUS f_true(bit | L): robust normal/Student-t with a linear mean
// mu(L)=a+b*L and partial-pooled log-sd sigma(L)=exp(a_ls+b_ls*L). Supersedes the
// sparse discrete-bin f_true whose empty low bins fell back to the long-read pooled
// mean — the bug that made short true reads (bit≈25) lose to f_0 (mu≈20). Fit ONLY
// on real short-read true-hit calibration data; never a long-read fallback.
inline double log_f_true_cont(const EMParams& params, double bit, uint16_t aln_len) {
    const double L = static_cast<double>(aln_len);
    const double mu = params.ft_mu_a + params.ft_mu_b * L;
    const double sigma = std::exp(params.ft_ls_a + params.ft_ls_b * L);
    return (params.ft_nu > 0.0) ? log_student_t(bit, mu, sigma, params.ft_nu)
                                : log_gaussian(bit, mu, sigma);
}

// v2 LENGTH-CONTINUOUS f_0(bit | L) as the Karlin–Altschul null-maximum density.
// The K–A CDF of the best decoy score is F_0(B|L) = exp[-2^-((B-loc(L))/beta)], a
// Gumbel; its density is (1/beta)·exp(-t - e^-t) with t=(B-loc(L))/beta. loc(L) is
// the length-continuous location (∝ log2 N_eff(L)) fit empirically on the homology-
// excluded spurious best-hits after the IDENTICAL top-hit selection. beta is the
// pooled scale (bit units; ≈1/ln2 for pure log2-odds bits, but calibrated).
inline double log_f0_cont(const EMParams& params, double bit, uint16_t aln_len) {
    const double beta = (params.f0_beta > 0.0) ? params.f0_beta : 1.0;
    const double loc = params.f0_loc_a + params.f0_loc_b * static_cast<double>(aln_len);
    const double t = (bit - loc) / beta;
    const double et = (t > 700.0) ? std::numeric_limits<double>::infinity() : std::exp(-t);
    return -std::log(beta) - t - et;
}

// Real-reference log-density term for a candidate's bit score. Under the coherent
// LR (use_calib + fitted f_true) this is the NORMALIZED true-hit density scaled by
// the frozen calib_temperature (lambda ≡ 1). Otherwise it is the LEGACY raw-bit
// temperature term lambda_b*bit*score_scale. score_scale (per-query max
// normalization) applies only to the legacy branch — f_true lives on raw bits.
inline double real_log_density(const EMParams& params, double bit,
                               uint16_t aln_len, double score_scale) {
    if (params.use_calib_cont)
        return params.calib_temperature * log_f_true_cont(params, bit, aln_len);
    if (params.use_calib && !params.true_bins.empty())
        return params.calib_temperature * log_f_true(params, bit, aln_len);
    return params.lambda_b * bit * score_scale;
}

// -------------------------------------------------------------------------
// Empirical null model I/O (decoy Gumbel per aligned-length bin)
// -------------------------------------------------------------------------
struct NullModel {
    std::vector<float> mu;    // per-bin location (on lambda_b*bit scale)
    std::vector<float> beta;  // per-bin scale
    std::vector<uint32_t> n;  // per-bin decoy query count used for the fit
    bool valid() const { return !mu.empty() && mu.size() == beta.size(); }
};

// Fit the empirical null from per-query MAX (lambda_b * bit_score) values,
// grouped by aligned-length bin, via method-of-moments Gumbel:
//   beta = sqrt(6)*sd/pi,  mu = mean - 0.5772*beta.
// Inputs: per bin, a list of per-query max tempered scores. Bins with fewer
// than min_n queries fall back to the pooled (all-query) fit.
NullModel fit_null_model(
    const std::vector<std::vector<double>>& per_bin_scores,
    uint32_t min_n = 200);

// Write/read the sidecar TSV: header + one row per bin
//   bin lo hi mu beta n
bool write_null_model(const std::string& path, const NullModel& m);
bool read_null_model(const std::string& path, NullModel& m);

// -------------------------------------------------------------------------
// Coherent-LR calibration model I/O (f_true + f_0, frozen once, loaded by EM)
// -------------------------------------------------------------------------
struct CalibModel {
    NullModel               null;        // f_0: per-bin Gumbel on RAW bit scale (v1)
    std::vector<TrueHitBin> true_bins;   // f_true: per-bin Gaussian mixture (v1)
    float                   censor = 0.0f;  // reporting-threshold bit (provenance)

    // v2 length-continuous section (sol short-fragment fix). When cont is true the
    // continuous params drive f_true/f_0; the v1 discrete sections above are still
    // written for provenance but the EM uses the continuous forms.
    bool   cont = false;
    double ft_mu_a = 0.0, ft_mu_b = 0.0;
    double ft_ls_a = 0.0, ft_ls_b = 0.0;
    double ft_nu   = 0.0;
    double f0_loc_a = 0.0, f0_loc_b = 0.0;
    double f0_beta  = 1.0;

    // v1 validity requires the discrete sections; v2 requires only a positive f_0
    // scale (the continuous params are always present once cont is set).
    bool valid() const {
        if (cont) return f0_beta > 0.0;
        return null.valid() && !true_bins.empty();
    }
};

// Fit f_true per length bin from labeled true-hit bit scores via method-of-moments
// Gaussian (single component; w1 left 0). Bins with fewer than min_n hits fall back
// to the pooled fit. Bins share kNullBinEdges with the null Gumbel.
std::vector<TrueHitBin> fit_true_model(
    const std::vector<std::vector<double>>& per_bin_bits,
    uint32_t min_n = 100);

// Fit the v2 length-continuous section from raw (aln_len, bit) pairs. true_LB are
// labeled true-source hits; null_LB are per-read best spurious hits. Fills the cont
// fields on m (sets m.cont=true) via robust linear mu(L)/loc(L) and pooled scales.
// nu is the Student-t dof for f_true (<=0 → Gaussian). Returns false if too sparse.
bool fit_calib_continuous(
    const std::vector<std::pair<double, double>>& true_LB,
    const std::vector<std::pair<double, double>>& null_LB,
    double nu, CalibModel& m);

// Combined calibration file: a versioned TSV carrying both sections.
bool write_calib_model(const std::string& path, const CalibModel& m);
bool read_calib_model(const std::string& path, CalibModel& m);

// Load a calibration file into EMParams (sets f_0, f_true, and use_calib=true).
inline bool load_calib_into_params(const std::string& path, EMParams& p) {
    CalibModel m;
    if (!read_calib_model(path, m) || !m.valid()) return false;
    p.null_mu = m.null.mu;
    p.null_beta = m.null.beta;
    p.true_bins = m.true_bins;
    p.true_censor = m.censor;
    p.use_calib = true;
    p.use_null = true;
    if (m.cont) {
        p.use_calib_cont = true;
        p.ft_mu_a = m.ft_mu_a; p.ft_mu_b = m.ft_mu_b;
        p.ft_ls_a = m.ft_ls_a; p.ft_ls_b = m.ft_ls_b;
        p.ft_nu   = m.ft_nu;
        p.f0_loc_a = m.f0_loc_a; p.f0_loc_b = m.f0_loc_b;
        p.f0_beta  = m.f0_beta;
    }
    return true;
}

// EM state vector
struct EMState {
    std::vector<double> weights;           // w_t: reference abundance weights (sum = 1)
    std::vector<double> gamma;             // gamma_{rt}: P(ref_t | read_r), flat array [aln_idx]
    std::vector<double> gamma_ancient;     // P(ancient, ref_t | read_r) (damage extension)
    std::vector<double> gamma_null;        // P(null | read_r), one per read (null extension)
    double w_null = 0.0;                    // soft abundance weight of the null column
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

    // Weight-distribution telemetry
    double effective_refs = 0.0;   // exp(H(w)) — effective number of references
    double hhi = 0.0;              // Σw² — Herfindahl-Hirschman concentration
    double top1_mass = 0.0;        // weight of single largest reference
    double top5_mass = 0.0;        // combined weight of top-5 references
    uint32_t n_active_refs = 0;    // references with weight > min_weight
};

// Compute weight-distribution telemetry from a normalised weight vector.
// min_weight is the floor used in EMParams (used to define "active").
inline void fill_weight_diagnostics(
    EMIterationDiagnostics& diag,
    const double* weights,
    uint32_t T,
    double min_weight)
{
    if (T == 0) return;

    double hhi = 0.0;
    double entropy = 0.0;
    double top1 = 0.0;
    uint32_t n_active = 0;

    // Partial buffer for top-5 (insertion sort on 5 elements)
    constexpr int K = 5;
    double top_k[K] = {};

    for (uint32_t t = 0; t < T; ++t) {
        const double w = weights[t];
        hhi += w * w;
        if (w > min_weight) {
            ++n_active;
            if (w > 0.0) entropy -= w * std::log(w);
        }
        if (w > top1) top1 = w;

        // Maintain top-K via insertion into sorted small buffer
        if (w > top_k[K - 1]) {
            top_k[K - 1] = w;
            for (int i = K - 2; i >= 0 && top_k[i] < top_k[i + 1]; --i)
                std::swap(top_k[i], top_k[i + 1]);
        }
    }

    double top5 = 0.0;
    for (int i = 0; i < K; ++i) top5 += top_k[i];

    diag.effective_refs = std::exp(entropy);
    diag.hhi            = hhi;
    diag.top1_mass      = top1;
    diag.top5_mass      = top5;
    diag.n_active_refs  = n_active;
}

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

        // Per-query score normalization: divide by max bit_score so all scores are in [0,1].
        // For unique mappers the scale is always 1.0 (self is max).
        double score_scale = 1.0;
        if (params.normalize_scores_per_query && deg > 1) {
            float max_bs = data.alignments[start].bit_score;
            for (uint32_t j = 1; j < deg; ++j)
                max_bs = std::max(max_bs, data.alignments[start + j].bit_score);
            if (max_bs > 0.0f) score_scale = 1.0 / static_cast<double>(max_bs);
        }

        // Optional null column appended after the deg real columns.
        const bool null_on = params.use_null && !params.null_mu.empty();
        log_scores.resize(null_on ? deg + 1 : deg);

        // Compute unnormalized log-responsibilities
        for (uint32_t j = 0; j < deg; ++j) {
            const auto& aln = data.alignments[start + j];
            log_scores[j] = std::log(std::max(weights[aln.ref_idx], params.min_weight))
                           + length_log_penalty(data, params, aln.ref_idx)
                           + coverage_log_penalty(params, aln.ref_idx)
                           + params.lambda_b * static_cast<double>(aln.bit_score) * score_scale;
        }
        if (null_on) {
            const double lp = read_null_log_density(params, data.alignments.data() + start, deg);
            log_scores[deg] = params.log_w_null + params.lambda_null * lp;
        }

        // Normalize via log-sum-exp (over real columns + optional null column)
        const uint32_t n_cols = null_on ? deg + 1 : deg;
        double lse = log_sum_exp(log_scores.data(), n_cols);
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

        double score_scale = 1.0;
        if (params.normalize_scores_per_query && deg > 1) {
            float max_bs = data.alignments[start].bit_score;
            for (uint32_t j = 1; j < deg; ++j)
                max_bs = std::max(max_bs, data.alignments[start + j].bit_score);
            if (max_bs > 0.0f) score_scale = 1.0 / static_cast<double>(max_bs);
        }

        // For damage-aware: 2*deg entries (ancient/modern for each alignment)
        log_scores.resize(2 * deg);

        for (uint32_t j = 0; j < deg; ++j) {
            const auto& aln = data.alignments[start + j];
            double log_w = std::log(std::max(weights[aln.ref_idx], params.min_weight));
            log_w += length_log_penalty(data, params, aln.ref_idx);
            log_w += coverage_log_penalty(params, aln.ref_idx);
            double score_contrib = params.lambda_b * static_cast<double>(aln.bit_score) * score_scale;

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
    if (!(params.alpha_prior > 0.0)) {
        throw std::invalid_argument("EMParams::alpha_prior must be > 0");
    }
    if (!(params.min_weight > 0.0)) {
        throw std::invalid_argument("EMParams::min_weight must be > 0");
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

            if (!std::isfinite(ll)) {
                break;
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
            if (!std::isfinite(ll)) {
                break;
            }
            m_step(data, params, state.gamma.data(), state.weights.data());
#ifndef NDEBUG
            {
                double w_sum = 0.0;
                for (uint32_t t = 0; t < T; ++t) w_sum += state.weights[t];
                assert(std::abs(w_sum - 1.0) < 1e-6 && "M-step broke simplex invariant");
            }
#endif
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
    double w_null = 0.0;               // Converged null-column weight (0 when use_null off)
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
    const float* gamma_ancient,   // Per-alignment ancient responsibility (may be null)
    const float* gamma_null,      // Per-alignment: read's null responsibility, replicated (may be null)
    const float* bit_score,       // Per-alignment raw bit score (may be null)
    const uint16_t* aln_len       // Per-alignment alignment length (may be null)
)>;

// Compute final gamma values in streaming fashion after EM convergence.
// Calls callback for each row group with computed responsibilities.
void streaming_em_finalize(
    ColumnarIndexReader& reader,
    const StreamingEMResult& result,
    const EMParams& params,
    StreamingGammaCallback callback);

} // namespace dart
