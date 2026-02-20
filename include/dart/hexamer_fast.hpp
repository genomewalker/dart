#pragma once

/**
 * @file hexamer_fast.hpp
 * @brief High-performance hexamer scoring with pre-computed log tables
 *
 * Optimizations:
 * 1. Pre-computed log2 lookup tables (avoid transcendentals in hot loop)
 * 2. Single-pass scoring with cached results
 * 3. SIMD-friendly memory layout
 */

#include "dart/hexamer_tables.hpp"
#include <cmath>
#include <array>
#include <cstring>

namespace dart {

// ============================================================================
// Pre-computed log2 lookup tables for each domain
// log2(freq * 4096 + 1e-10) pre-computed for all 4096 hexamers
// ============================================================================

namespace detail {

// Compute log2(freq * 4096 + 1e-10) at compile time
constexpr float compute_log2_score(float freq) {
    // Can't use std::log2 in constexpr, so we'll initialize at runtime once
    return 0.0f;  // Initialized in runtime log tables
}

// Runtime initialization of log tables
struct LogTableInitializer {
    alignas(64) float gtdb_log[4096];
    alignas(64) float fungi_log[4096];
    alignas(64) float protozoa_log[4096];
    alignas(64) float invertebrate_log[4096];
    alignas(64) float plant_log[4096];
    alignas(64) float vertebrate_mammalian_log[4096];
    alignas(64) float vertebrate_other_log[4096];
    alignas(64) float viral_log[4096];

    LogTableInitializer() {
        for (int i = 0; i < 4096; ++i) {
            gtdb_log[i] = std::log2(GTDB_HEXAMER_FREQ[i] * 4096.0f + 1e-10f);
            fungi_log[i] = std::log2(FUNGI_HEXAMER_FREQ[i] * 4096.0f + 1e-10f);
            protozoa_log[i] = std::log2(PROTOZOA_HEXAMER_FREQ[i] * 4096.0f + 1e-10f);
            invertebrate_log[i] = std::log2(INVERTEBRATE_HEXAMER_FREQ[i] * 4096.0f + 1e-10f);
            plant_log[i] = std::log2(PLANT_HEXAMER_FREQ[i] * 4096.0f + 1e-10f);
            vertebrate_mammalian_log[i] = std::log2(VERTEBRATE_MAMMALIAN_HEXAMER_FREQ[i] * 4096.0f + 1e-10f);
            vertebrate_other_log[i] = std::log2(VERTEBRATE_OTHER_HEXAMER_FREQ[i] * 4096.0f + 1e-10f);
            viral_log[i] = std::log2(VIRAL_HEXAMER_FREQ[i] * 4096.0f + 1e-10f);
        }
    }
};

inline const LogTableInitializer& get_log_tables() {
    static const LogTableInitializer tables;
    return tables;
}

} // namespace detail

// ============================================================================
// Fast hexamer encoding with lookup table
// ============================================================================

// Pre-computed base encoding table (handles both upper and lower case)
alignas(64) constexpr int8_t FAST_BASE_MAP[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  // 0-15
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  // 16-31
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  // 32-47
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  // 48-63
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,  // 64-79  (A=0, C=1, G=2)
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  // 80-95  (T=3)
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,  // 96-111 (a=0, c=1, g=2)
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  // 112-127 (t=3)
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
};

// Encode hexamer to 12-bit code (0-4095), returns UINT32_MAX on invalid
inline uint32_t encode_hexamer_fast(const char* seq) {
    int8_t b0 = FAST_BASE_MAP[static_cast<uint8_t>(seq[0])];
    int8_t b1 = FAST_BASE_MAP[static_cast<uint8_t>(seq[1])];
    int8_t b2 = FAST_BASE_MAP[static_cast<uint8_t>(seq[2])];
    int8_t b3 = FAST_BASE_MAP[static_cast<uint8_t>(seq[3])];
    int8_t b4 = FAST_BASE_MAP[static_cast<uint8_t>(seq[4])];
    int8_t b5 = FAST_BASE_MAP[static_cast<uint8_t>(seq[5])];

    // Check all valid in one comparison
    if ((b0 | b1 | b2 | b3 | b4 | b5) < 0) return UINT32_MAX;

    return (static_cast<uint32_t>(b0) << 10) |
           (static_cast<uint32_t>(b1) << 8) |
           (static_cast<uint32_t>(b2) << 6) |
           (static_cast<uint32_t>(b3) << 4) |
           (static_cast<uint32_t>(b4) << 2) |
           static_cast<uint32_t>(b5);
}

// ============================================================================
// Single-pass multi-domain scoring
// Score all 8 domains in one pass, cache results
// ============================================================================

struct HexamerScoreCache {
    alignas(64) float domain_scores[8] = {0.0f};
    int hexamer_count = 0;
    bool valid = false;

    void reset() {
        std::memset(domain_scores, 0, sizeof(domain_scores));
        hexamer_count = 0;
        valid = false;
    }
};

// Thread-local score cache to avoid repeated computation
inline HexamerScoreCache& get_thread_score_cache() {
    static thread_local HexamerScoreCache cache;
    return cache;
}

// Score all domains in single pass using pre-computed log tables
// Returns cached result - call invalidate_score_cache() when sequence changes
inline MultiDomainResult score_all_domains_fast(const char* seq, size_t len, int frame) {
    MultiDomainResult result;

    if (len < 6) return result;

    const auto& logs = detail::get_log_tables();

    alignas(32) float scores[8] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    int count = 0;

    // Single pass through all hexamers
    for (size_t i = frame; i + 5 < len; i += 3) {
        uint32_t code = encode_hexamer_fast(seq + i);
        if (code == UINT32_MAX) continue;

        // Direct table lookup - no transcendentals!
        scores[0] += logs.gtdb_log[code];
        scores[1] += logs.fungi_log[code];
        scores[2] += logs.protozoa_log[code];
        scores[3] += logs.invertebrate_log[code];
        scores[4] += logs.plant_log[code];
        scores[5] += logs.vertebrate_mammalian_log[code];
        scores[6] += logs.vertebrate_other_log[code];
        scores[7] += logs.viral_log[code];
        count++;
    }

    if (count == 0) return result;

    // Find max and best domain
    float max_score = scores[0];
    int best_idx = 0;
    for (int i = 1; i < 8; ++i) {
        if (scores[i] > max_score) {
            max_score = scores[i];
            best_idx = i;
        }
    }

    // Softmax with temperature scaling
    float sum_exp = 0.0f;
    for (int i = 0; i < 8; ++i) {
        scores[i] = std::exp((scores[i] - max_score) * 0.1f);
        sum_exp += scores[i];
    }

    if (sum_exp > 0.0f) {
        float inv_sum = 1.0f / sum_exp;
        result.gtdb_prob = scores[0] * inv_sum;
        result.fungi_prob = scores[1] * inv_sum;
        result.protozoa_prob = scores[2] * inv_sum;
        result.invertebrate_prob = scores[3] * inv_sum;
        result.plant_prob = scores[4] * inv_sum;
        result.vertebrate_mammalian_prob = scores[5] * inv_sum;
        result.vertebrate_other_prob = scores[6] * inv_sum;
        result.viral_prob = scores[7] * inv_sum;
    }

    result.best_domain = static_cast<Domain>(best_idx + 1);
    result.best_score = max_score;

    return result;
}

// Overload for std::string
inline MultiDomainResult score_all_domains_fast(const std::string& seq, int frame) {
    return score_all_domains_fast(seq.c_str(), seq.length(), frame);
}

// ============================================================================
// Fast dicodon scoring using cached domain probabilities
// ============================================================================

inline float calculate_dicodon_score_fast(const char* seq, size_t len, int frame,
                                          const MultiDomainResult& probs) {
    if (len < static_cast<size_t>(frame + 6)) return 0.5f;

    float log_prob_sum = 0.0f;
    int hexamer_count = 0;

    static constexpr float BACKGROUND_FREQ = 1.0f / 4096.0f;
    static constexpr float LOG_BACKGROUND = -12.0f;  // log(1/4096)

    for (size_t i = frame; i + 5 < len; i += 3) {
        uint32_t code = encode_hexamer_fast(seq + i);
        if (code == UINT32_MAX) continue;

        // Weighted frequency - direct table lookup
        float freq = probs.gtdb_prob * GTDB_HEXAMER_FREQ[code] +
                     probs.fungi_prob * FUNGI_HEXAMER_FREQ[code] +
                     probs.protozoa_prob * PROTOZOA_HEXAMER_FREQ[code] +
                     probs.invertebrate_prob * INVERTEBRATE_HEXAMER_FREQ[code] +
                     probs.plant_prob * PLANT_HEXAMER_FREQ[code] +
                     probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_HEXAMER_FREQ[code] +
                     probs.vertebrate_other_prob * VERTEBRATE_OTHER_HEXAMER_FREQ[code] +
                     probs.viral_prob * VIRAL_HEXAMER_FREQ[code];

        if (freq < 1e-9f) freq = BACKGROUND_FREQ * 0.01f;

        log_prob_sum += std::log(freq) - LOG_BACKGROUND;
        hexamer_count++;
    }

    if (hexamer_count == 0) return 0.5f;

    float avg_llr = log_prob_sum / hexamer_count;
    float score = 0.5f + 0.15f * avg_llr;

    return std::max(0.0f, std::min(1.0f, score));
}

inline float calculate_dicodon_score_fast(const std::string& seq, int frame,
                                          const MultiDomainResult& probs) {
    return calculate_dicodon_score_fast(seq.c_str(), seq.length(), frame, probs);
}

} // namespace dart
