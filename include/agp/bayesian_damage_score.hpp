#pragma once
// Bayesian damage scoring with log-odds fusion
// Calibrated ancient/modern classification
//
// Key features:
// - Empirical Bayes: all parameters derived from sample data
// - Tempered evidence: w0 < 1 downweights absence evidence for robustness
// - Stratified baseline: q_modern varies by GC content
// - Three-state output: ancient/uncertain/modern (not forced binary)
// - Channel B gating: site evidence only trusted when validated

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace agp {

// Three-state classification (not forced binary)
enum class DamageClass : uint8_t {
    Modern       = 0,  // likely undamaged
    Uncertain    = 1,  // insufficient evidence
    Ancient      = 2   // likely damaged
};

// Evidence tier for interpretability
enum class EvidenceTier : uint8_t {
    NoEvidence   = 0,  // p_read low, no damage sites
    SitesOnly    = 1,  // damage sites but low p_read
    TerminalOnly = 2,  // high p_read but no damage sites
    Convergent   = 3   // both terminal and site evidence
};

// Site evidence from alignment scan (full Bernoulli LLR model)
struct SiteEvidence {
    uint32_t m = 0;           // susceptible opportunities (C or G in target)
    uint32_t k = 0;           // observed damage hits (C→T or G→A)
    float sum_qA = 0.0f;      // sum of position-weighted ancient hazards (all positions)
    float sum_log_qA_hits = 0.0f;  // sum of log(q1/q0) at damage sites only (legacy)
    float q_eff = 0.0f;       // effective per-site rate = sum_qA / m

    // Full Bernoulli LLR: sum over ALL opportunities (hits + non-hits)
    // LLR = sum_j [ y_j*log(q1_j/q0_j) + (1-y_j)*log((1-q1_j)/(1-q0_j)) ]
    float llr_bernoulli = 0.0f;  // full position-weighted LLR
    float sum_exp_decay = 0.0f;  // sum of exp(-λ*dist) for computing effective w
    bool has_bernoulli = false;  // true if llr_bernoulli was computed
};

// Modern baseline estimate from sample (can be stratified by GC)
struct ModernBaseline {
    float q_modern = 0.005f;  // estimated baseline substitution rate
    uint64_t M0 = 0;          // pooled opportunities
    uint64_t K0 = 0;          // pooled hits
    size_t n_reads = 0;       // number of reads used
};

// Stratified baseline by GC content (4 bins: 0-25%, 25-50%, 50-75%, 75-100%)
struct StratifiedBaseline {
    std::array<ModernBaseline, 4> gc_bins{};
    ModernBaseline global{};

    // Get bin index for GC fraction
    static size_t gc_bin(float gc_frac) noexcept {
        if (gc_frac < 0.25f) return 0;
        if (gc_frac < 0.50f) return 1;
        if (gc_frac < 0.75f) return 2;
        return 3;
    }

    // Get q_modern for specific GC content (with shrinkage to global)
    float q_for_gc(float gc_frac, float shrinkage = 0.5f) const noexcept {
        size_t bin = gc_bin(gc_frac);
        const auto& b = gc_bins[bin];
        // Shrink toward global when bin has few reads
        float n_eff = static_cast<float>(b.n_reads);
        float w = n_eff / (n_eff + 100.0f);  // 100 reads = equal weight
        return w * b.q_modern + (1.0f - w) * global.q_modern;
    }
};

// Scoring parameters with tempering and gating
struct BayesianScoreParams {
    float pi = 0.10f;              // prior P(ancient)
    float pi0 = 0.10f;             // prior used when computing p_read
    float q_modern = 0.005f;       // baseline substitution rate (or use stratified)
    float w0 = 0.3f;               // tempering: weight for absence evidence [0,1]
                                   // w0=0 -> one-sided (only hits count)
                                   // w0=1 -> full Bayesian (absence penalizes)
    float terminal_threshold = 0.50f;  // p_read threshold for "strong" evidence
    float site_cap = 3.0f;         // max logBF_sites contribution (avoid single-hit overreaction)
    uint32_t min_opportunities = 3; // minimum m to trust site evidence
    bool channel_b_valid = true;   // whether sample-level damage is validated
    float eps = 1e-6f;             // numerical stability

    // Identity evidence parameters (hybrid scoring)
    // Ancient proteins have higher alignment identity (+2.8% on average)
    // because they match their reference better than spurious modern matches
    float k_identity = 20.0f;      // coefficient for identity evidence
    float identity_baseline = 0.90f;  // neutral point (identity = baseline -> no evidence)
    bool use_identity = false;     // enable identity evidence (requires fident input)

    // Damage informativeness gating
    // When false: terminal evidence is zeroed (sample has no detectable damage signal)
    // Site evidence and identity evidence still apply
    bool damage_informative = true;

    // Sample-level damage rate for fixed-average site evidence formula
    // Uses d_max * 0.5 as position-independent average damage rate
    // This matches the original proven formula (commit 0224b4a) which achieved AUC 0.82
    float d_max = 0.0f;

    // Thresholds for 3-state classification
    float ancient_threshold = 0.60f;   // posterior >= this -> Ancient
    float modern_threshold = 0.25f;    // posterior <= this -> Modern
                                       // between -> Uncertain
};

// Full decomposed output
struct BayesianScoreOutput {
    float posterior = 0.0f;        // P(ancient | all evidence)
    float logit_posterior = 0.0f;  // log-odds form
    float logBF_terminal = 0.0f;   // Bayes factor from terminal pattern
    float logBF_sites = 0.0f;      // Bayes factor from alignment sites
    float logBF_identity = 0.0f;   // Bayes factor from alignment identity
    float q_eff = 0.0f;            // effective ancient rate for this read
    uint32_t m = 0;                // opportunities
    uint32_t k = 0;                // hits
    float sum_qA = 0.0f;           // position-weighted hazard sum
    EvidenceTier tier = EvidenceTier::NoEvidence;
    DamageClass damage_class = DamageClass::Uncertain;  // 3-state output
    bool informative = true;  // false when sample damage is uninformative
};

// Precomputed exponential decay lookup table
class ExpDecayLUT {
public:
    ExpDecayLUT() = default;

    ExpDecayLUT(float lambda, size_t max_dist = 300)
        : lambda_(lambda), table_(max_dist + 1) {
        for (size_t d = 0; d <= max_dist; ++d) {
            table_[d] = std::exp(-lambda_ * static_cast<float>(d));
        }
    }

    float at(uint32_t dist) const noexcept {
        if (dist < table_.size()) return table_[dist];
        return std::exp(-lambda_ * static_cast<float>(dist));
    }

    float lambda() const noexcept { return lambda_; }

private:
    float lambda_ = 0.3f;
    std::vector<float> table_;
};

namespace detail {

// Base encoding for branch-free comparison
enum : uint8_t { B_OTHER = 0, B_A = 1, B_C = 2, B_G = 3, B_T = 4, B_GAP = 5 };

inline const std::array<uint8_t, 256>& base_lut() {
    static const std::array<uint8_t, 256> lut = [] {
        std::array<uint8_t, 256> x{};
        x.fill(B_OTHER);
        x[static_cast<uint8_t>('-')] = B_GAP;
        x[static_cast<uint8_t>('A')] = x[static_cast<uint8_t>('a')] = B_A;
        x[static_cast<uint8_t>('C')] = x[static_cast<uint8_t>('c')] = B_C;
        x[static_cast<uint8_t>('G')] = x[static_cast<uint8_t>('g')] = B_G;
        x[static_cast<uint8_t>('T')] = x[static_cast<uint8_t>('t')] = B_T;
        return x;
    }();
    return lut;
}

inline double clamp01(double x, double eps) noexcept {
    return std::min(1.0 - eps, std::max(eps, x));
}

inline double logit(double p) noexcept {
    return std::log(p) - std::log1p(-p);
}

inline double inv_logit(double z) noexcept {
    if (z >= 0.0) {
        const double e = std::exp(-z);
        return 1.0 / (1.0 + e);
    }
    const double e = std::exp(z);
    return e / (1.0 + e);
}

} // namespace detail

// Compute site evidence from alignment in single pass
// qaln/taln: aligned sequences with gaps
// qstart: 0-based position in read where alignment starts
// qlen: total read length
// d_max: sample-level damage rate
inline SiteEvidence compute_site_evidence(
    const std::string& qaln,
    const std::string& taln,
    size_t qstart,
    size_t qlen,
    float d_max,
    const ExpDecayLUT& decay_lut) noexcept
{
    SiteEvidence out{};
    if (qlen == 0 || qaln.empty() || taln.empty()) return out;

    const auto& lut = detail::base_lut();
    const size_t n = std::min(qaln.size(), taln.size());
    const uint8_t* q = reinterpret_cast<const uint8_t*>(qaln.data());
    const uint8_t* t = reinterpret_cast<const uint8_t*>(taln.data());

    uint32_t r = static_cast<uint32_t>(qstart);
    const uint32_t qlen_m1 = static_cast<uint32_t>(qlen > 0 ? qlen - 1 : 0);

    for (size_t i = 0; i < n; ++i) {
        const uint8_t qc = lut[q[i]];
        const uint8_t tc = lut[t[i]];

        // Both must be non-gap ACGT
        const bool q_valid = (qc >= detail::B_A && qc <= detail::B_T);
        const bool t_valid = (tc >= detail::B_A && tc <= detail::B_T);
        const bool both = q_valid && t_valid;

        if (both) {
            // Check if target is C or G (susceptible position)
            const bool isC = (tc == detail::B_C);
            const bool isG = (tc == detail::B_G);

            if (isC || isG) {
                out.m++;  // susceptible opportunity

                // Position-weighted hazard
                // C→T damage is 5' biased, G→A is 3' biased
                uint32_t dist = isC ? r : (qlen_m1 - r);
                float q_damage = d_max * decay_lut.at(dist);
                out.sum_qA += q_damage;

                // Check for damage-consistent substitution
                if ((isC && qc == detail::B_T) || (isG && qc == detail::B_A)) {
                    out.k++;  // damage hit
                    // Ancient rate = damage + baseline substitution
                    // This ensures BF >= 1 (damage is always more likely for ancient)
                    // q_baseline ~= 0.005 (sequencing error + natural variation)
                    constexpr float Q_BASELINE = 0.005f;
                    float q_ancient = q_damage + Q_BASELINE;
                    // log(q_ancient / q_modern) = log(1 + q_damage/q_baseline)
                    out.sum_log_qA_hits += std::log(q_ancient / Q_BASELINE);
                }
            }
        }

        // Advance read position only on non-gap query
        if (qc != detail::B_GAP) {
            r++;
        }
    }

    out.q_eff = (out.m > 0) ? (out.sum_qA / static_cast<float>(out.m)) : 0.0f;
    return out;
}

// Estimate q_modern from batch of site evidence (likely-undamaged reads)
// mk_pairs: vector of (m, k) from reads with low p_read
inline ModernBaseline estimate_modern_baseline(
    const std::vector<std::pair<uint32_t, uint32_t>>& mk_pairs,
    float alpha_prior = 1.0f,    // Beta prior alpha
    float beta_prior = 199.0f,   // Beta prior beta (mean ~0.5%)
    float q_min = 1e-4f,
    float q_max = 0.02f) noexcept
{
    ModernBaseline out;
    out.n_reads = mk_pairs.size();

    for (const auto& [m, k] : mk_pairs) {
        out.M0 += m;
        out.K0 += k;
    }

    if (out.M0 == 0) {
        out.q_modern = alpha_prior / (alpha_prior + beta_prior);
        return out;
    }

    double q = (static_cast<double>(out.K0) + alpha_prior) /
               (static_cast<double>(out.M0) + alpha_prior + beta_prior);
    out.q_modern = static_cast<float>(std::clamp(q, static_cast<double>(q_min), static_cast<double>(q_max)));
    return out;
}

// Compute final Bayesian score with decomposition
// Implements: tempered evidence, site capping, Channel B gating, 3-state output
// Optional fident: alignment identity (0-1), pass -1 to disable identity evidence
inline BayesianScoreOutput compute_bayesian_score(
    float p_read,
    const SiteEvidence& ev,
    const BayesianScoreParams& params,
    float fident = -1.0f) noexcept
{
    using namespace detail;

    const double eps = std::max(1e-12, static_cast<double>(params.eps));
    const double pi = clamp01(params.pi, eps);
    const double pi0 = clamp01(params.pi0, eps);
    const double pr = clamp01(p_read, eps);
    const double qm = clamp01(params.q_modern, eps);
    const double w0 = std::clamp(static_cast<double>(params.w0), 0.0, 1.0);

    // Terminal evidence: convert p_read posterior to Bayes factor
    // BF = P(data|A)/P(data|M) = [p_read/(1-p_read)] / [pi0/(1-pi0)]
    // When p_read is very low (no terminal signal detected), treat as no evidence (BF=1, logBF=0)
    // rather than strong evidence against damage (which would overwhelm site evidence)
    //
    // When !damage_informative: sample has no detectable damage signal (artifact or
    // unvalidated), so terminal evidence is meaningless and zeroed out.
    double logBF_terminal = 0.0;
    if (params.damage_informative && p_read > 0.01f) {
        logBF_terminal = logit(pr) - logit(pi0);
    }

    // Site evidence with tempering and gating
    double logBF_sites = 0.0;
    double q_eff_clamped = 0.0;

    // Only use site evidence if:
    // 1. Enough opportunities (m >= min_opportunities)
    // 2. Either we have hits OR (Channel B valid AND w0 > 0 for absence evidence)
    // Note: Always use positive site evidence (hits), only gate negative evidence
    const bool use_site_evidence = ev.m >= params.min_opportunities;

    if (use_site_evidence && ev.m > 0) {
        q_eff_clamped = clamp01(ev.q_eff, eps);
        const double k = static_cast<double>(ev.k);
        const double m = static_cast<double>(ev.m);

        // Tempered Bayes factor:
        // logBF = k*log(qA/qM) + w0*(m-k)*log((1-qA)/(1-qM))
        // With w0 < 1, absence evidence is downweighted for robustness

        if (ev.k > 0) {
            // Fixed-average formula for logBF_sites
            //
            // Why fixed-average (not position-weighted):
            // - Alignment AA positions don't correspond to original read positions
            // - Damage at terminal nucleotides appears at various AA positions after translation
            // - Interior alignment positions may correspond to terminal read damage
            //
            // Fix 1: AA_SCALE = 0.30 (P(observable AA event | nucleotide damage))
            // Fix 2: w = E[exp(-λ*dist)] corrected from 0.5 to data-driven value
            //
            // Data-driven w from per-read sum_exp_decay if available, else estimate:
            // With λ=0.3 and typical read lengths: w ≈ 0.10-0.15
            const double w = (ev.sum_exp_decay > 0.0f && ev.m > 0)
                ? static_cast<double>(ev.sum_exp_decay) / static_cast<double>(ev.m)
                : 0.12;  // Conservative estimate for ~50bp reads

            constexpr double AA_SCALE = 0.30;  // Nucleotide → AA damage rate conversion
            const double d_max_scaled = static_cast<double>(params.d_max) * AA_SCALE;

            if (d_max_scaled > qm) {
                // logBF_sites = k * log((d_max_scaled * w + qm) / qm)
                const double avg_q_ancient = d_max_scaled * w + qm;
                logBF_sites = k * std::log(avg_q_ancient / qm);
            } else if (ev.sum_log_qA_hits > 0.0f) {
                // Fallback: position-weighted sum (when d_max not available)
                logBF_sites = static_cast<double>(ev.sum_log_qA_hits);
            }
        }

        // Tempered negative evidence (absence of damage where expected)
        // Only apply if Channel B validates sample-level damage (otherwise no expected damage)
        if (params.channel_b_valid && w0 > 0.0 && ev.m > ev.k) {
            // Average ancient rate over all susceptible positions
            double q_ancient_avg = clamp01(ev.sum_qA / m + qm, eps);
            double log_ratio_no_damage = std::log((1.0 - q_ancient_avg) / (1.0 - qm));
            logBF_sites += w0 * (m - k) * log_ratio_no_damage;
        }

        // Cap site contribution to avoid single-hit overreaction
        logBF_sites = std::clamp(logBF_sites,
                                  -static_cast<double>(params.site_cap),
                                  static_cast<double>(params.site_cap));
    }

    // Identity evidence (hybrid scoring)
    // Ancient proteins have higher alignment identity because they're real matches
    // logBF_identity = k * (fident - baseline), where k=20, baseline=0.90
    double logBF_identity = 0.0;
    if (params.use_identity && fident >= 0.0f) {
        logBF_identity = static_cast<double>(params.k_identity) *
                         (static_cast<double>(fident) - static_cast<double>(params.identity_baseline));
        // Cap to avoid extreme values from outlier identities
        logBF_identity = std::clamp(logBF_identity, -5.0, 5.0);
    }

    // Combine in log-odds space
    const double logit_post = logit(pi) + logBF_terminal + logBF_sites + logBF_identity;
    const double post = inv_logit(logit_post);

    // Determine evidence tier
    // When uninformative, terminal signal is not usable evidence
    const bool has_sites = (ev.k > 0);
    const bool has_terminal = params.damage_informative && (p_read >= params.terminal_threshold);
    EvidenceTier tier;
    if (!params.damage_informative) {
        tier = has_sites ? EvidenceTier::SitesOnly : EvidenceTier::NoEvidence;
    } else if (has_terminal && has_sites) {
        tier = EvidenceTier::Convergent;
    } else if (has_terminal) {
        tier = EvidenceTier::TerminalOnly;
    } else if (has_sites) {
        tier = EvidenceTier::SitesOnly;
    } else {
        tier = EvidenceTier::NoEvidence;
    }

    // 3-state classification
    DamageClass damage_class;
    if (post >= params.ancient_threshold) {
        damage_class = DamageClass::Ancient;
    } else if (post <= params.modern_threshold) {
        damage_class = DamageClass::Modern;
    } else {
        damage_class = DamageClass::Uncertain;
    }

    BayesianScoreOutput out;
    out.posterior = static_cast<float>(post);
    out.logit_posterior = static_cast<float>(logit_post);
    out.logBF_terminal = static_cast<float>(logBF_terminal);
    out.logBF_sites = static_cast<float>(logBF_sites);
    out.logBF_identity = static_cast<float>(logBF_identity);
    out.q_eff = static_cast<float>(q_eff_clamped);
    out.m = ev.m;
    out.k = ev.k;
    out.sum_qA = ev.sum_qA;
    out.tier = tier;
    out.damage_class = damage_class;
    out.informative = params.damage_informative;
    return out;
}

inline const char* tier_name(EvidenceTier tier) {
    switch (tier) {
        case EvidenceTier::NoEvidence: return "none";
        case EvidenceTier::SitesOnly: return "sites";
        case EvidenceTier::TerminalOnly: return "terminal";
        case EvidenceTier::Convergent: return "convergent";
        default: return "unknown";
    }
}

inline const char* damage_class_name(DamageClass dc) {
    switch (dc) {
        case DamageClass::Modern: return "modern";
        case DamageClass::Uncertain: return "uncertain";
        case DamageClass::Ancient: return "ancient";
        default: return "unknown";
    }
}

// Estimate stratified baseline from batch of reads
// mk_gc_tuples: vector of (m, k, gc_fraction) from likely-undamaged reads
inline StratifiedBaseline estimate_stratified_baseline(
    const std::vector<std::tuple<uint32_t, uint32_t, float>>& mk_gc_tuples,
    float alpha_prior = 1.0f,
    float beta_prior = 199.0f) noexcept
{
    StratifiedBaseline out;

    // Accumulate per-bin
    for (const auto& [m, k, gc] : mk_gc_tuples) {
        size_t bin = StratifiedBaseline::gc_bin(gc);
        out.gc_bins[bin].M0 += m;
        out.gc_bins[bin].K0 += k;
        out.gc_bins[bin].n_reads++;
        out.global.M0 += m;
        out.global.K0 += k;
        out.global.n_reads++;
    }

    // Compute q_modern for each bin (Beta-Binomial posterior mean)
    auto compute_q = [alpha_prior, beta_prior](const ModernBaseline& b) -> float {
        if (b.M0 == 0) return alpha_prior / (alpha_prior + beta_prior);
        double q = (static_cast<double>(b.K0) + alpha_prior) /
                   (static_cast<double>(b.M0) + alpha_prior + beta_prior);
        return static_cast<float>(std::clamp(q, 1e-4, 0.02));
    };

    for (auto& bin : out.gc_bins) {
        bin.q_modern = compute_q(bin);
    }
    out.global.q_modern = compute_q(out.global);

    return out;
}

// Compute damage detectability score (0.0 = uninformative, 1.0 = highly informative)
// Based on AGD v3 header fields: damage_validated, damage_artifact, d_max,
// stop_decay_llr (Channel B), terminal_shift
//
// Use cases:
// - Gate Bayesian scoring: if detectability == 0, posterior is meaningless
// - Downstream filtering: weight results by sample quality
// - Reporting: flag samples where damage signal is unreliable
inline float compute_damage_detectability(
    float d_max,
    bool damage_validated,
    bool damage_artifact,
    bool channel_b_valid,
    float stop_decay_llr,
    float terminal_shift) noexcept
{
    // No damage detected at all
    if (d_max <= 0.0f) return 0.0f;

    // Channel A fired but Channel B says artifact
    if (damage_artifact) return 0.0f;

    // Channel B not valid (insufficient data) -- use d_max as weak proxy
    if (!channel_b_valid) {
        // d_max alone is unreliable, but non-zero means some signal
        // Scale: 10% d_max -> 0.3, 30% -> 0.5 (capped, never fully confident)
        return std::min(0.5f, d_max * 3.0f);
    }

    // Both channels validated: combine d_max strength with Channel B confidence
    if (damage_validated) {
        // Terminal shift strength: higher shift = clearer signal
        float shift_score = std::min(1.0f, std::abs(terminal_shift) * 5.0f);
        // Channel B LLR strength: higher LLR = more confident
        float llr_score = std::min(1.0f, std::max(0.0f, stop_decay_llr / 10000.0f));
        // d_max strength
        float dmax_score = std::min(1.0f, d_max * 5.0f);
        // Combine: geometric-ish mean, floor at 0.5 (validated = at least moderate)
        return std::max(0.5f, (shift_score + llr_score + dmax_score) / 3.0f);
    }

    // Fallback: not validated, not artifact, Channel B valid but flat
    return 0.0f;
}

// Length-binned BPA statistics for short protein z-score normalization
// Bins: 0=15-24aa, 1=25-39aa, 2=40-59aa, 3=60+aa
struct LengthBinStats {
    static constexpr int N_BINS = 4;
    struct Bin {
        double sum_bpa = 0.0;
        double sum_bpa_sq = 0.0;
        size_t count = 0;
        float mu = 0.0f;
        float sigma = 1.0f;
        bool valid = false;
    };
    std::array<Bin, N_BINS> bins{};

    static int get_bin(size_t aa_len) {
        if (aa_len < 25) return 0;
        if (aa_len < 40) return 1;
        if (aa_len < 60) return 2;
        return 3;
    }

    void add(size_t aa_len, float bpa) {
        int b = get_bin(aa_len);
        bins[b].sum_bpa += bpa;
        bins[b].sum_bpa_sq += bpa * bpa;
        bins[b].count++;
    }

    void finalize() {
        for (auto& b : bins) {
            if (b.count < 10) continue;
            b.mu = static_cast<float>(b.sum_bpa / b.count);
            double var = (b.sum_bpa_sq / b.count) - (b.mu * b.mu);
            b.sigma = static_cast<float>(std::sqrt(std::max(var, 1e-6)));
            b.valid = true;
        }
    }

    float z_score(size_t aa_len, float bpa) const {
        int b = get_bin(aa_len);
        if (!bins[b].valid) return 0.0f;
        return (bpa - bins[b].mu) / bins[b].sigma;
    }
};

// 4-class protein classification combining damage informativeness with scoring
enum class AncientClassification : uint8_t {
    AncientConfident = 0,  // damage_informative AND posterior >= threshold
    AncientLikely    = 1,  // NOT damage_informative AND source_supported
    Undetermined     = 2,  // NOT damage_informative AND no strong source signal
    ModernConfident  = 3   // damage_informative AND posterior <= modern_threshold
};

inline const char* classification_name(AncientClassification c) {
    switch (c) {
        case AncientClassification::AncientConfident: return "ancient_confident";
        case AncientClassification::AncientLikely: return "ancient_likely";
        case AncientClassification::Undetermined: return "undetermined";
        case AncientClassification::ModernConfident: return "modern_confident";
        default: return "unknown";
    }
}

// Classify a protein using damage informativeness + Bayesian posterior + quality metrics
// When damage is informative, posterior drives the classification.
// When uninformative, fall back to alignment quality heuristics.
inline AncientClassification classify_protein(
    bool damage_informative,
    float posterior,
    float fident,
    float z_bpa,
    float delta_bits,
    float ancient_threshold = 0.60f,
    float modern_threshold = 0.25f) noexcept
{
    if (damage_informative) {
        if (posterior >= ancient_threshold)
            return AncientClassification::AncientConfident;
        if (posterior <= modern_threshold)
            return AncientClassification::ModernConfident;
        return AncientClassification::Undetermined;
    }

    // Uninformative: use alignment quality as proxy
    // source_supported = high identity + good length-normalized score + clear best hit
    bool source_supported = (fident >= 0.90f) && (z_bpa >= 0.5f) && (delta_bits >= 5.0f);
    return source_supported
        ? AncientClassification::AncientLikely
        : AncientClassification::Undetermined;
}

} // namespace agp
