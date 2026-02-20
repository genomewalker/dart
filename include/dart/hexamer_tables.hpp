#pragma once

/**
 * @file hexamer_tables.hpp
 * @brief Unified multi-domain hexamer frequency interface
 *
 * Provides domain selection and hexamer frequency lookup across all trained domains:
 * fungi, protozoa, invertebrate, plant, vertebrate_mammalian, vertebrate_other, viral
 */

#include <string>
#include <cstdint>

namespace dart {

// Encode hexamer to 12-bit code (0-4095)
// Defined here BEFORE including domain headers to avoid redefinition
inline uint32_t encode_hexamer(const char* seq) {
    static const int8_t base_map[256] = {
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
    };
    uint32_t code = 0;
    for (int i = 0; i < 6; i++) {
        int8_t b = base_map[(uint8_t)seq[i]];
        if (b < 0) return UINT32_MAX;
        code = (code << 2) | b;
    }
    return code;
}

} // namespace dart

// Include all domain-specific hexamer tables
#include "dart/gtdb_hexamer_table.hpp"
#include "dart/fungi_hexamer_table.hpp"
#include "dart/protozoa_hexamer_table.hpp"
#include "dart/invertebrate_hexamer_table.hpp"
#include "dart/plant_hexamer_table.hpp"
#include "dart/vertebrate_mammalian_hexamer_table.hpp"
#include "dart/vertebrate_other_hexamer_table.hpp"
#include "dart/viral_hexamer_table.hpp"

namespace dart {

// Domain enum for selection
enum class Domain {
    META,  // Ensemble of all domains (default for ancient metagenomes)
    GTDB,  // Bacteria + Archaea from GTDB
    FUNGI,
    PROTOZOA,
    INVERTEBRATE,
    PLANT,
    VERTEBRATE_MAMMALIAN,
    VERTEBRATE_OTHER,
    VIRAL
};

// Parse domain string to enum (defaults to META for ancient metagenomes)
inline Domain parse_domain(const std::string& domain) {
    if (domain == "meta" || domain == "metagenome" || domain == "all") return Domain::META;
    if (domain == "gtdb" || domain == "bacteria" || domain == "archaea") return Domain::GTDB;
    if (domain == "fungi") return Domain::FUNGI;
    if (domain == "protozoa") return Domain::PROTOZOA;
    if (domain == "invertebrate") return Domain::INVERTEBRATE;
    if (domain == "plant") return Domain::PLANT;
    if (domain == "mammal" || domain == "vertebrate_mammalian") return Domain::VERTEBRATE_MAMMALIAN;
    if (domain == "vertebrate" || domain == "vertebrate_other") return Domain::VERTEBRATE_OTHER;
    if (domain == "viral") return Domain::VIRAL;
    // Default to META (ensemble of all domains) for ancient metagenomes
    return Domain::META;
}

// Get domain name string
inline const char* domain_name(Domain domain) {
    switch (domain) {
        case Domain::META: return "meta";
        case Domain::GTDB: return "gtdb";
        case Domain::FUNGI: return "fungi";
        case Domain::PROTOZOA: return "protozoa";
        case Domain::INVERTEBRATE: return "invertebrate";
        case Domain::PLANT: return "plant";
        case Domain::VERTEBRATE_MAMMALIAN: return "vertebrate_mammalian";
        case Domain::VERTEBRATE_OTHER: return "vertebrate_other";
        case Domain::VIRAL: return "viral";
        default: return "meta";
    }
}

// Global active domain (thread-safe read, set at startup)
inline Domain& active_domain() {
    static Domain domain = Domain::META;
    return domain;
}

inline void set_active_domain(Domain domain) {
    active_domain() = domain;
}

inline Domain get_active_domain() {
    return active_domain();
}

// Get hexamer frequency for a specific domain
inline float get_hexamer_freq(uint32_t code, Domain domain) {
    if (code >= 4096) return 0.0f;

    switch (domain) {
        case Domain::GTDB:
            return GTDB_HEXAMER_FREQ[code];
        case Domain::FUNGI:
            return FUNGI_HEXAMER_FREQ[code];
        case Domain::PROTOZOA:
            return PROTOZOA_HEXAMER_FREQ[code];
        case Domain::INVERTEBRATE:
            return INVERTEBRATE_HEXAMER_FREQ[code];
        case Domain::PLANT:
            return PLANT_HEXAMER_FREQ[code];
        case Domain::VERTEBRATE_MAMMALIAN:
            return VERTEBRATE_MAMMALIAN_HEXAMER_FREQ[code];
        case Domain::VERTEBRATE_OTHER:
            return VERTEBRATE_OTHER_HEXAMER_FREQ[code];
        case Domain::VIRAL:
            return VIRAL_HEXAMER_FREQ[code];
        case Domain::META:
        default:
            // For META mode, use GTDB as default (bacteria-focused)
            // The ensemble weighting happens via score_all_domains()
            return GTDB_HEXAMER_FREQ[code];
    }
}

// Get hexamer frequency using active domain
inline float get_hexamer_freq(uint32_t code) {
    return get_hexamer_freq(code, active_domain());
}

// Get hexamer frequency from sequence string using active domain
inline float get_hexamer_freq(const char* seq) {
    uint32_t code = encode_hexamer(seq);
    if (code == UINT32_MAX) return 0.0f;
    return get_hexamer_freq(code, active_domain());
}

// Get hexamer frequency from sequence string for specific domain
inline float get_hexamer_freq(const char* seq, Domain domain) {
    uint32_t code = encode_hexamer(seq);
    if (code == UINT32_MAX) return 0.0f;
    return get_hexamer_freq(code, domain);
}

// Get log-ratio score (coding vs random) for active domain
inline float get_hexamer_score(const char* seq) {
    float freq = get_hexamer_freq(seq);
    if (freq <= 0.0f) return -10.0f;
    constexpr float RANDOM_FREQ = 1.0f / 4096.0f;
    return std::log2(freq / RANDOM_FREQ);
}

// Get log-ratio score for specific domain
inline float get_hexamer_score(const char* seq, Domain domain) {
    float freq = get_hexamer_freq(seq, domain);
    if (freq <= 0.0f) return -10.0f;
    constexpr float RANDOM_FREQ = 1.0f / 4096.0f;
    return std::log2(freq / RANDOM_FREQ);
}

// Multi-domain result for ensemble scoring
struct MultiDomainResult {
    float gtdb_prob = 0.0f;  // Bacteria + Archaea
    float fungi_prob = 0.0f;
    float protozoa_prob = 0.0f;
    float invertebrate_prob = 0.0f;
    float plant_prob = 0.0f;
    float vertebrate_mammalian_prob = 0.0f;
    float vertebrate_other_prob = 0.0f;
    float viral_prob = 0.0f;
    Domain best_domain = Domain::GTDB;
    float best_score = 0.0f;
};

// Thread-local ensemble mode state
inline bool& ensemble_mode_enabled() {
    static thread_local bool enabled = false;
    return enabled;
}

inline MultiDomainResult& current_domain_probs() {
    static thread_local MultiDomainResult result;
    return result;
}

inline void set_ensemble_mode(bool enabled) {
    ensemble_mode_enabled() = enabled;
}

inline bool get_ensemble_mode() {
    return ensemble_mode_enabled();
}

inline void set_domain_probs(const MultiDomainResult& result) {
    current_domain_probs() = result;
}

inline MultiDomainResult& get_domain_probs() {
    return current_domain_probs();
}

// Get hexamer frequency for specific domain from sequence string
inline float get_domain_hexamer_freq(Domain domain, const char* seq) {
    return get_hexamer_freq(seq, domain);
}

// Score sequence against all domains and return probabilities
inline MultiDomainResult score_all_domains(const std::string& seq, int frame) {
    MultiDomainResult result;

    if (seq.length() < 6) return result;

    // Score each domain (8 domains now: GTDB + 7 eukaryotic)
    float scores[8] = {0.0f};
    const char* data = seq.c_str();
    size_t len = seq.length();

    for (size_t i = frame; i + 5 < len; i += 3) {
        uint32_t code = encode_hexamer(data + i);
        if (code == UINT32_MAX) continue;

        scores[0] += std::log2(GTDB_HEXAMER_FREQ[code] * 4096.0f + 1e-10f);
        scores[1] += std::log2(FUNGI_HEXAMER_FREQ[code] * 4096.0f + 1e-10f);
        scores[2] += std::log2(PROTOZOA_HEXAMER_FREQ[code] * 4096.0f + 1e-10f);
        scores[3] += std::log2(INVERTEBRATE_HEXAMER_FREQ[code] * 4096.0f + 1e-10f);
        scores[4] += std::log2(PLANT_HEXAMER_FREQ[code] * 4096.0f + 1e-10f);
        scores[5] += std::log2(VERTEBRATE_MAMMALIAN_HEXAMER_FREQ[code] * 4096.0f + 1e-10f);
        scores[6] += std::log2(VERTEBRATE_OTHER_HEXAMER_FREQ[code] * 4096.0f + 1e-10f);
        scores[7] += std::log2(VIRAL_HEXAMER_FREQ[code] * 4096.0f + 1e-10f);
    }

    // Convert to probabilities via softmax
    float max_score = scores[0];
    int best_idx = 0;
    for (int i = 1; i < 8; ++i) {
        if (scores[i] > max_score) {
            max_score = scores[i];
            best_idx = i;
        }
    }

    float sum_exp = 0.0f;
    for (int i = 0; i < 8; ++i) {
        scores[i] = std::exp((scores[i] - max_score) * 0.1f);  // Temperature scaling
        sum_exp += scores[i];
    }

    if (sum_exp > 0.0f) {
        result.gtdb_prob = scores[0] / sum_exp;
        result.fungi_prob = scores[1] / sum_exp;
        result.protozoa_prob = scores[2] / sum_exp;
        result.invertebrate_prob = scores[3] / sum_exp;
        result.plant_prob = scores[4] / sum_exp;
        result.vertebrate_mammalian_prob = scores[5] / sum_exp;
        result.vertebrate_other_prob = scores[6] / sum_exp;
        result.viral_prob = scores[7] / sum_exp;
    }

    // Map best_idx to Domain enum (GTDB=1, FUNGI=2, etc. in enum)
    // best_idx 0=GTDB, 1=FUNGI, 2=PROTOZOA, etc.
    result.best_domain = static_cast<Domain>(best_idx + 1);  // +1 because META=0
    result.best_score = max_score;

    return result;
}

// Get ensemble-weighted hexamer frequency
inline float get_ensemble_hexamer_freq(uint32_t code) {
    if (code >= 4096) return 0.0f;

    if (!ensemble_mode_enabled()) {
        return get_hexamer_freq(code, active_domain());
    }

    const auto& probs = current_domain_probs();
    return probs.gtdb_prob * GTDB_HEXAMER_FREQ[code] +
           probs.fungi_prob * FUNGI_HEXAMER_FREQ[code] +
           probs.protozoa_prob * PROTOZOA_HEXAMER_FREQ[code] +
           probs.invertebrate_prob * INVERTEBRATE_HEXAMER_FREQ[code] +
           probs.plant_prob * PLANT_HEXAMER_FREQ[code] +
           probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_HEXAMER_FREQ[code] +
           probs.vertebrate_other_prob * VERTEBRATE_OTHER_HEXAMER_FREQ[code] +
           probs.viral_prob * VIRAL_HEXAMER_FREQ[code];
}

// Fast uppercase conversion
inline char fast_upper_hex(char c) {
    return (c >= 'a' && c <= 'z') ? (c - 32) : c;
}

// Calculate dicodon score weighted by domain probabilities
// Used in ensemble mode for metagenome scoring
inline float calculate_dicodon_score_weighted(const std::string& seq, int frame, const MultiDomainResult& probs) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    float log_prob_sum = 0.0f;
    int hexamer_count = 0;

    static constexpr float BACKGROUND_FREQ = 1.0f / 4096.0f;

    char hexbuf[7];
    hexbuf[6] = '\0';

    for (size_t i = frame; i + 5 < seq.length(); i += 3) {
        hexbuf[0] = fast_upper_hex(seq[i]);
        hexbuf[1] = fast_upper_hex(seq[i + 1]);
        hexbuf[2] = fast_upper_hex(seq[i + 2]);
        hexbuf[3] = fast_upper_hex(seq[i + 3]);
        hexbuf[4] = fast_upper_hex(seq[i + 4]);
        hexbuf[5] = fast_upper_hex(seq[i + 5]);

        if (hexbuf[0] == 'N' || hexbuf[1] == 'N' || hexbuf[2] == 'N' ||
            hexbuf[3] == 'N' || hexbuf[4] == 'N' || hexbuf[5] == 'N') continue;

        uint32_t code = encode_hexamer(hexbuf);
        if (code == UINT32_MAX) continue;

        // Weighted frequency across all domains
        float freq = probs.gtdb_prob * GTDB_HEXAMER_FREQ[code] +
                     probs.fungi_prob * FUNGI_HEXAMER_FREQ[code] +
                     probs.protozoa_prob * PROTOZOA_HEXAMER_FREQ[code] +
                     probs.invertebrate_prob * INVERTEBRATE_HEXAMER_FREQ[code] +
                     probs.plant_prob * PLANT_HEXAMER_FREQ[code] +
                     probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_HEXAMER_FREQ[code] +
                     probs.vertebrate_other_prob * VERTEBRATE_OTHER_HEXAMER_FREQ[code] +
                     probs.viral_prob * VIRAL_HEXAMER_FREQ[code];

        if (freq < 1e-9f) {
            freq = BACKGROUND_FREQ * 0.01f;
        }

        log_prob_sum += std::log(freq / BACKGROUND_FREQ);
        hexamer_count++;
    }

    if (hexamer_count == 0) return 0.5f;

    float avg_llr = log_prob_sum / hexamer_count;
    float score = 0.5f + 0.15f * avg_llr;

    return std::max(0.0f, std::min(1.0f, score));
}

} // namespace dart
