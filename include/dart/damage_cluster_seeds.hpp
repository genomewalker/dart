#pragma once
// Seed emission and frame scoring for damage-aware clustering
//
// Dual-seed strategy:
// - Strict seeds (k=7): no canonicalization, high precision
// - Tolerant seeds (k=5): damage-canonicalized, terminal regions only
//
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "dart/damage_cluster.hpp"

namespace dart {
namespace damage_cluster {

struct SeedHit {
    uint64_t hash = 0;
    uint16_t aa_pos = 0;
    bool terminal = false;
    bool tolerant = false;
};

namespace seed_detail {

inline char up(char c) {
    return (c >= 'a' && c <= 'z') ? static_cast<char>(c - 'a' + 'A') : c;
}

inline bool aa_code_strict(char c, uint8_t& code) {
    const uint8_t idx = aa_index(up(c));
    if (idx >= 20) return false; // reject stop(*) and unknown
    code = idx;
    return true;
}

// Conservative canonicalization to reduce false joins in mixed metagenomes.
inline char canonical_ct(char c) {
    switch (up(c)) {
        case 'W': return 'R';
        case 'Y': return 'H';
        case 'I': return 'T';
        case 'M': return 'T';
        case 'F': return 'S';
        default: return up(c);
    }
}

inline char canonical_ga(char c) {
    switch (up(c)) {
        case 'K': return 'E';
        case 'N': return 'D';
        case 'I': return 'V';
        case 'M': return 'V';
        default: return up(c);
    }
}

inline bool aa_code_ct(char c, uint8_t& code) {
    const uint8_t idx = aa_index(canonical_ct(c));
    if (idx >= 20) return false;
    code = idx;
    return true;
}

inline bool aa_code_ga(char c, uint8_t& code) {
    const uint8_t idx = aa_index(canonical_ga(c));
    if (idx >= 20) return false;
    code = idx;
    return true;
}

inline uint32_t hash_codes_32(const uint8_t* x, uint8_t n) {
    uint32_t h = 2166136261u; // FNV-1a
    for (uint8_t i = 0; i < n; ++i) {
        h ^= static_cast<uint32_t>(x[i] + 1u);
        h *= 16777619u;
    }
    return h;
}

inline uint64_t hash_codes_64(const uint8_t* x, uint8_t n) {
    uint64_t h = 1469598103934665603ULL; // FNV-1a 64
    for (uint8_t i = 0; i < n; ++i) {
        h ^= static_cast<uint64_t>(x[i] + 1u);
        h *= 1099511628211ULL;
    }
    return h;
}

inline uint8_t default_s_for_k(uint8_t k) {
    if (k <= 3) return 2;
    const uint8_t s = static_cast<uint8_t>(k / 2 + 1); // k=7->4, k=5->3
    return (s < k) ? s : static_cast<uint8_t>(k - 1);
}

template <typename CodeFn>
inline bool boundary_syncmer_hash(
    const char* kmer,
    uint8_t k,
    uint8_t s,
    CodeFn&& code_fn,
    uint64_t& out_hash
) {
    if (k == 0 || s == 0 || s > k) return false;

    // Convert full k-mer once.
    uint8_t kc[32];
    if (k > 32) return false; // guard; not expected here
    for (uint8_t i = 0; i < k; ++i) {
        if (!code_fn(kmer[i], kc[i])) return false;
    }

    // Find minimum s-mer hash and position.
    uint32_t min_h = std::numeric_limits<uint32_t>::max();
    uint8_t min_pos = 0;
    for (uint8_t i = 0; i <= static_cast<uint8_t>(k - s); ++i) {
        const uint32_t h = hash_codes_32(kc + i, s);
        if (h < min_h) {
            min_h = h;
            min_pos = i;
        }
    }

    const uint8_t right = static_cast<uint8_t>(k - s);
    const bool is_boundary = (min_pos == 0) || (min_pos == right);
    if (!is_boundary) return false;

    out_hash = hash_codes_64(kc, k);
    return true;
}

} // namespace seed_detail

// -----------------------------
// 1) Seed emission API
// -----------------------------
// Strict: no canonicalization, high precision
inline void emit_strict_seeds(
    const std::string& protein,
    std::vector<SeedHit>& out,
    uint8_t k = 7
) {
    out.clear();
    if (k < 3 || protein.size() < k) return;

    const uint8_t s = seed_detail::default_s_for_k(k);
    const size_t n = protein.size();
    out.reserve(n / 2 + 8);

    for (size_t pos = 0; pos + k <= n; ++pos) {
        uint64_t h = 0;
        if (!seed_detail::boundary_syncmer_hash(
                protein.data() + pos, k, s,
                [](char c, uint8_t& code) { return seed_detail::aa_code_strict(c, code); },
                h)) {
            continue;
        }

        out.push_back(SeedHit{
            h,
            static_cast<uint16_t>(pos),
            false,  // terminal
            false   // tolerant
        });
    }
}

inline std::vector<SeedHit> emit_strict_seeds(
    const std::string& protein,
    uint8_t k = 7
) {
    std::vector<SeedHit> out;
    emit_strict_seeds(protein, out, k);
    return out;
}

// Tolerant: canonicalized, terminal regions only
inline void emit_tolerant_seeds(
    const std::string& protein,
    std::vector<SeedHit>& out,
    uint8_t k = 5,
    uint8_t terminal_window = 12
) {
    out.clear();
    if (k < 3 || protein.size() < k) return;

    const uint8_t s = seed_detail::default_s_for_k(k);
    const size_t n = protein.size();
    const size_t tail_start = (n > terminal_window) ? (n - terminal_window) : 0;
    out.reserve(static_cast<size_t>(2) * terminal_window + 8);

    // Side tags prevent accidental 5'/3' merging for ambiguous canonicalizations.
    static constexpr uint64_t TAG_CT = 0x9E3779B97F4A7C15ULL;
    static constexpr uint64_t TAG_GA = 0xC2B2AE3D27D4EB4FULL;

    for (size_t pos = 0; pos + k <= n; ++pos) {
        const bool in_5p = (pos < terminal_window);
        const bool in_3p = ((pos + k) > tail_start);
        if (!in_5p && !in_3p) continue;

        if (in_5p) {
            uint64_t h = 0;
            if (seed_detail::boundary_syncmer_hash(
                    protein.data() + pos, k, s,
                    [](char c, uint8_t& code) { return seed_detail::aa_code_ct(c, code); },
                    h)) {
                out.push_back(SeedHit{
                    (h ^ TAG_CT),
                    static_cast<uint16_t>(pos),
                    true,
                    true
                });
            }
        }

        if (in_3p) {
            uint64_t h = 0;
            if (seed_detail::boundary_syncmer_hash(
                    protein.data() + pos, k, s,
                    [](char c, uint8_t& code) { return seed_detail::aa_code_ga(c, code); },
                    h)) {
                out.push_back(SeedHit{
                    (h ^ TAG_GA),
                    static_cast<uint16_t>(pos),
                    true,
                    true
                });
            }
        }
    }
}

inline std::vector<SeedHit> emit_tolerant_seeds(
    const std::string& protein,
    uint8_t k = 5,
    uint8_t terminal_window = 12
) {
    std::vector<SeedHit> out;
    emit_tolerant_seeds(protein, out, k, terminal_window);
    return out;
}

// -----------------------------
// 2) Frame scoring
// -----------------------------
inline float safe_z(float x, float mu, float sigma) {
    return (sigma > 1e-6f) ? ((x - mu) / sigma) : 0.0f;
}

// Frame weights container (learned from warmup, not hardcoded)
struct FrameWeights {
    float w_code = 0.0f;  // Weight for coding score
    float w_sup = 0.0f;   // Weight for cluster support
    float w_term = 0.0f;  // Weight for terminal damage LLR

    // Derivation: β_k = (μ⁺ - μ⁻) / (σ²_pooled + ε), then normalize
    // where + = high-confidence winner frames, - = runner-up frames
    static FrameWeights from_discriminative(
        float mu_plus_code, float mu_minus_code, float var_pooled_code,
        float mu_plus_sup,  float mu_minus_sup,  float var_pooled_sup,
        float mu_plus_term, float mu_minus_term, float var_pooled_term
    ) {
        constexpr float eps = 1e-4f;
        const float beta_code = (mu_plus_code - mu_minus_code) / (var_pooled_code + eps);
        const float beta_sup  = (mu_plus_sup  - mu_minus_sup)  / (var_pooled_sup  + eps);
        const float beta_term = (mu_plus_term - mu_minus_term) / (var_pooled_term + eps);

        const float sum = std::max(eps,
            std::max(0.0f, beta_code) +
            std::max(0.0f, beta_sup) +
            std::max(0.0f, beta_term));

        FrameWeights w;
        w.w_code = std::max(0.0f, beta_code) / sum;
        w.w_sup  = std::max(0.0f, beta_sup)  / sum;
        w.w_term = std::max(0.0f, beta_term) / sum;
        return w;
    }

    // Fallback defaults based on damage level (when warmup data insufficient)
    static FrameWeights defaults_for_damage(float d_max) {
        FrameWeights w;
        if (d_max < 0.05f) {
            // Low damage: coding signal dominates, terminal unreliable
            w.w_code = 0.75f; w.w_sup = 0.25f; w.w_term = 0.0f;
        } else if (d_max < 0.20f) {
            // Moderate damage: all signals contribute
            w.w_code = 0.55f; w.w_sup = 0.30f; w.w_term = 0.15f;
        } else {
            // High damage: terminal signal becomes informative
            w.w_code = 0.50f; w.w_sup = 0.28f; w.w_term = 0.22f;
        }
        return w;
    }
};

inline float score_frame(
    float code_score,
    float cluster_support,
    float terminal_llr,
    float mu_code, float sigma_code,
    float mu_sup, float sigma_sup,
    float mu_term, float sigma_term,
    const FrameWeights& weights
) {
    const float z_code = safe_z(code_score, mu_code, sigma_code);
    const float z_sup  = safe_z(cluster_support, mu_sup, sigma_sup);
    const float z_term = safe_z(terminal_llr, mu_term, sigma_term);
    return weights.w_code * z_code + weights.w_sup * z_sup + weights.w_term * z_term;
}

// Backward-compatible overload using d_max-based defaults
inline float score_frame(
    float code_score,
    float cluster_support,
    float terminal_llr,
    float mu_code, float sigma_code,
    float mu_sup, float sigma_sup,
    float mu_term, float sigma_term,
    float d_max = 0.0f
) {
    const FrameWeights w = FrameWeights::defaults_for_damage(d_max);
    return score_frame(code_score, cluster_support, terminal_llr,
                       mu_code, sigma_code, mu_sup, sigma_sup, mu_term, sigma_term, w);
}

// -----------------------------
// 3) Terminal damage LLR
// -----------------------------
namespace llr_detail {

inline bool ct_consistent(char c, char o) {
    c = seed_detail::up(c);
    o = seed_detail::up(o);
    // Conservative subset (higher precision).
    return (c == 'R' && o == 'W') ||
           (c == 'H' && o == 'Y') ||
           (c == 'T' && (o == 'I' || o == 'M')) ||
           (c == 'S' && o == 'F') ||
           (c == 'Q' && o == '*');
}

inline bool ga_consistent(char c, char o) {
    c = seed_detail::up(c);
    o = seed_detail::up(o);
    // Conservative subset (higher precision).
    return (c == 'E' && o == 'K') ||
           (c == 'D' && o == 'N') ||
           (c == 'V' && (o == 'I' || o == 'M')) ||
           (c == 'W' && o == '*');
}

inline float site_quality(float depth, float freq) {
    const float qd = std::min(1.0f, depth / 6.0f);
    const float qf = std::max(0.0f, (freq - 0.5f) / 0.5f);
    return qd * qf;
}

template <size_t MAX_W>
inline void add_site_llr(
    float& llr,
    const std::string& q,
    size_t q_idx,
    bool five_prime,
    uint8_t pos_from_end,
    const TerminalConsensus<MAX_W>& cons,
    float d_max,
    float lambda,
    float p0
) {
    const auto s = cons.site(five_prime, pos_from_end);
    if (s.depth < 1.0f) return;
    if (s.consensus == 'X') return;

    const char obs = q[q_idx];
    if (obs == s.consensus) return;

    const float q_site = site_quality(s.depth, s.freq);
    if (q_site <= 0.0f) return;

    const float pi = clampf(
        d_max * std::exp(-std::max(0.0f, lambda) * static_cast<float>(pos_from_end)),
        1e-4f, 0.95f
    );

    const bool consistent = five_prime
        ? ct_consistent(s.consensus, obs)
        : ga_consistent(s.consensus, obs);

    if (consistent) {
        llr += q_site * std::log(pi / p0);
    } else {
        llr += q_site * std::log((1.0f - pi) / (1.0f - p0));
    }
}

} // namespace llr_detail

template <size_t MAX_W>
inline float compute_terminal_damage_llr(
    const std::string& query_protein,
    const TerminalConsensus<MAX_W>& consensus,
    float d_max,
    float lambda,
    float p0 = 0.10f
) {
    if (query_protein.empty()) return 0.0f;

    const size_t n = query_protein.size();
    const uint8_t w = static_cast<uint8_t>(std::min<size_t>(consensus.window(), n));
    if (w == 0) return 0.0f;

    const float p0c = clampf(p0, 1e-4f, 0.49f);
    float llr = 0.0f;

    // 5' side
    for (uint8_t i = 0; i < w; ++i) {
        llr_detail::add_site_llr(
            llr, query_protein, static_cast<size_t>(i), true, i,
            consensus, d_max, lambda, p0c
        );
    }

    // 3' side, avoid double-counting overlap for very short proteins
    for (uint8_t i = 0; i < w; ++i) {
        const size_t q_idx = n - 1 - static_cast<size_t>(i);
        if (q_idx < static_cast<size_t>(w)) continue;
        llr_detail::add_site_llr(
            llr, query_protein, q_idx, false, i,
            consensus, d_max, lambda, p0c
        );
    }

    return llr;
}

// Exact-signature convenience overload for default TerminalConsensus<>
inline float compute_terminal_damage_llr(
    const std::string& query_protein,
    const TerminalConsensus<>& consensus,
    float d_max,
    float lambda,
    float p0 = 0.10f
) {
    return compute_terminal_damage_llr<24>(query_protein, consensus, d_max, lambda, p0);
}

// -----------------------------
// 4) Softmax frame posteriors
// -----------------------------

// Temperature estimation from score margins
// Derivation: τ = sqrt(3 * Var(Δ)) / π, where Δ = s₁ - s₂
// This gives calibrated posteriors where P(correct) matches empirical accuracy
struct SoftmaxTemperature {
    float tau = 1.0f;

    // Estimate from score margin statistics (collected during warmup)
    static SoftmaxTemperature from_margin_variance(double sum_delta, double sumsq_delta, size_t n) {
        SoftmaxTemperature t;
        if (n < 50) {
            t.tau = 1.0f;
            return t;
        }
        const double mu = sum_delta / n;
        const double var = (sumsq_delta / n) - mu * mu;
        if (var > 0) {
            // τ = sqrt(3 * Var(Δ)) / π
            t.tau = static_cast<float>(std::sqrt(3.0 * var) / 3.14159265359);
            t.tau = std::max(0.5f, std::min(3.0f, t.tau));
        }
        return t;
    }

    // Fallback based on damage level
    static SoftmaxTemperature defaults_for_damage(float d_max) {
        SoftmaxTemperature t;
        // Higher damage → more uncertainty → higher temperature
        t.tau = 1.0f + 0.5f * d_max;
        return t;
    }
};

inline void softmax_frame_posteriors(
    const float* scores,
    size_t n,
    float* posteriors,
    float tau = 1.0f  // Now defaults to 1.0, should be provided from calibration
) {
    if (n == 0) return;
    tau = std::max(1e-3f, tau);

    float mx = scores[0];
    for (size_t i = 1; i < n; ++i) mx = std::max(mx, scores[i]);

    float sum = 0.0f;
    for (size_t i = 0; i < n; ++i) {
        const float e = std::exp((scores[i] - mx) / tau);
        posteriors[i] = e;
        sum += e;
    }

    if (sum <= 0.0f) {
        const float u = 1.0f / static_cast<float>(n);
        for (size_t i = 0; i < n; ++i) posteriors[i] = u;
        return;
    }

    const float inv = 1.0f / sum;
    for (size_t i = 0; i < n; ++i) posteriors[i] *= inv;
}

inline std::array<float, 3> softmax_frame_posteriors(
    const std::array<float, 3>& scores,
    uint8_t n,
    float tau = 1.25f
) {
    std::array<float, 3> out{0.0f, 0.0f, 0.0f};
    n = static_cast<uint8_t>(std::min<uint8_t>(n, 3));
    softmax_frame_posteriors(scores.data(), n, out.data(), tau);
    return out;
}

} // namespace damage_cluster
} // namespace dart
