#pragma once
// Bayesian damage scoring with log-odds fusion
// Designed with GPT-5.3-codex for principled, calibrated ancient/modern classification
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

// Site evidence from alignment scan
struct SiteEvidence {
    uint32_t m = 0;           // susceptible opportunities (C or G in target)
    uint32_t k = 0;           // observed damage hits (C→T or G→A)
    float sum_qA = 0.0f;      // sum of position-weighted ancient hazards (all positions)
    float sum_log_qA_hits = 0.0f;  // sum of log(q) at damage sites only
    float q_eff = 0.0f;       // effective per-site rate = sum_qA / m
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
    float q_eff = 0.0f;            // effective ancient rate for this read
    uint32_t m = 0;                // opportunities
    uint32_t k = 0;                // hits
    float sum_qA = 0.0f;           // position-weighted hazard sum
    EvidenceTier tier = EvidenceTier::NoEvidence;
    DamageClass damage_class = DamageClass::Uncertain;  // 3-state output
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
inline BayesianScoreOutput compute_bayesian_score(
    float p_read,
    const SiteEvidence& ev,
    const BayesianScoreParams& params) noexcept
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
    double logBF_terminal = 0.0;
    if (p_read > 0.01f) {  // Only use terminal evidence if meaningful signal
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
            // Positive evidence from observed damage sites
            // sum_log_qA_hits = sum of log((q_damage + q_baseline) / q_baseline)
            logBF_sites = static_cast<double>(ev.sum_log_qA_hits);
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

    // Combine in log-odds space
    const double logit_post = logit(pi) + logBF_terminal + logBF_sites;
    const double post = inv_logit(logit_post);

    // Determine evidence tier
    const bool has_sites = (ev.k > 0);
    const bool has_terminal = (p_read >= params.terminal_threshold);
    EvidenceTier tier;
    if (has_terminal && has_sites) {
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
    out.q_eff = static_cast<float>(q_eff_clamped);
    out.m = ev.m;
    out.k = ev.k;
    out.sum_qA = ev.sum_qA;
    out.tier = tier;
    out.damage_class = damage_class;
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

} // namespace agp
