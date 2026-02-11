#pragma once

/**
 * @file positional_hexamer.hpp
 * @brief Unified multi-domain positional hexamer frequency interface
 *
 * Provides domain selection and positional hexamer frequency lookup across all trained domains:
 * fungi, protozoa, invertebrate, plant, vertebrate_mammalian, vertebrate_other, viral
 *
 * Positional hexamers track codon usage patterns at START, INTERNAL, and END positions
 * of coding sequences for improved frame detection.
 */

#include <string>
#include <cstdint>
#include <cmath>

// Include hexamer_tables.hpp for Domain enum and encode_hexamer
#include "agp/hexamer_tables.hpp"

// Include all domain-specific positional hexamer tables
#include "agp/gtdb_positional_hexamer.hpp"
#include "agp/fungi_positional_hexamer.hpp"
#include "agp/protozoa_positional_hexamer.hpp"
#include "agp/invertebrate_positional_hexamer.hpp"
#include "agp/plant_positional_hexamer.hpp"
#include "agp/vertebrate_mammalian_positional_hexamer.hpp"
#include "agp/vertebrate_other_positional_hexamer.hpp"
#include "agp/viral_positional_hexamer.hpp"

namespace agp {

// Position types for positional hexamer lookup
enum class HexamerPosition {
    START,     // First 2 codons (positions 0-5)
    INTERNAL,  // Middle codons
    END        // Last 2 codons
};

// Get START position hexamer frequency for specific domain
inline float get_start_hexamer_freq(uint32_t code, Domain domain) {
    if (code >= 4096) return 0.0f;

    switch (domain) {
        case Domain::GTDB:
            return GTDB_START_HEXAMER_FREQ[code];
        case Domain::FUNGI:
            return FUNGI_START_HEXAMER_FREQ[code];
        case Domain::PROTOZOA:
            return PROTOZOA_START_HEXAMER_FREQ[code];
        case Domain::INVERTEBRATE:
            return INVERTEBRATE_START_HEXAMER_FREQ[code];
        case Domain::PLANT:
            return PLANT_START_HEXAMER_FREQ[code];
        case Domain::VERTEBRATE_MAMMALIAN:
            return VERTEBRATE_MAMMALIAN_START_HEXAMER_FREQ[code];
        case Domain::VERTEBRATE_OTHER:
            return VERTEBRATE_OTHER_START_HEXAMER_FREQ[code];
        case Domain::VIRAL:
            return VIRAL_START_HEXAMER_FREQ[code];
        case Domain::META:
        default:
            // For META mode, use GTDB as default (bacteria-focused)
            // Ensemble weighting happens via get_ensemble_start_hexamer_freq()
            return GTDB_START_HEXAMER_FREQ[code];
    }
}

// Get INTERNAL position hexamer frequency for specific domain
inline float get_internal_hexamer_freq(uint32_t code, Domain domain) {
    if (code >= 4096) return 0.0f;

    switch (domain) {
        case Domain::GTDB:
            return GTDB_INTERNAL_HEXAMER_FREQ[code];
        case Domain::FUNGI:
            return FUNGI_INTERNAL_HEXAMER_FREQ[code];
        case Domain::PROTOZOA:
            return PROTOZOA_INTERNAL_HEXAMER_FREQ[code];
        case Domain::INVERTEBRATE:
            return INVERTEBRATE_INTERNAL_HEXAMER_FREQ[code];
        case Domain::PLANT:
            return PLANT_INTERNAL_HEXAMER_FREQ[code];
        case Domain::VERTEBRATE_MAMMALIAN:
            return VERTEBRATE_MAMMALIAN_INTERNAL_HEXAMER_FREQ[code];
        case Domain::VERTEBRATE_OTHER:
            return VERTEBRATE_OTHER_INTERNAL_HEXAMER_FREQ[code];
        case Domain::VIRAL:
            return VIRAL_INTERNAL_HEXAMER_FREQ[code];
        case Domain::META:
        default:
            // For META mode, use GTDB as default (bacteria-focused)
            // Ensemble weighting happens via get_ensemble_internal_hexamer_freq()
            return GTDB_INTERNAL_HEXAMER_FREQ[code];
    }
}

// Get END position hexamer frequency for specific domain
inline float get_end_hexamer_freq(uint32_t code, Domain domain) {
    if (code >= 4096) return 0.0f;

    switch (domain) {
        case Domain::GTDB:
            return GTDB_END_HEXAMER_FREQ[code];
        case Domain::FUNGI:
            return FUNGI_END_HEXAMER_FREQ[code];
        case Domain::PROTOZOA:
            return PROTOZOA_END_HEXAMER_FREQ[code];
        case Domain::INVERTEBRATE:
            return INVERTEBRATE_END_HEXAMER_FREQ[code];
        case Domain::PLANT:
            return PLANT_END_HEXAMER_FREQ[code];
        case Domain::VERTEBRATE_MAMMALIAN:
            return VERTEBRATE_MAMMALIAN_END_HEXAMER_FREQ[code];
        case Domain::VERTEBRATE_OTHER:
            return VERTEBRATE_OTHER_END_HEXAMER_FREQ[code];
        case Domain::VIRAL:
            return VIRAL_END_HEXAMER_FREQ[code];
        case Domain::META:
        default:
            // For META mode, use GTDB as default (bacteria-focused)
            // Ensemble weighting happens via get_ensemble_end_hexamer_freq()
            return GTDB_END_HEXAMER_FREQ[code];
    }
}

// Get positional hexamer frequency using active domain
inline float get_start_hexamer_freq(uint32_t code) {
    return get_start_hexamer_freq(code, get_active_domain());
}

inline float get_internal_hexamer_freq(uint32_t code) {
    return get_internal_hexamer_freq(code, get_active_domain());
}

inline float get_end_hexamer_freq(uint32_t code) {
    return get_end_hexamer_freq(code, get_active_domain());
}

// Get positional hexamer frequency by position type
inline float get_positional_hexamer_freq(uint32_t code, HexamerPosition pos, Domain domain) {
    switch (pos) {
        case HexamerPosition::START:
            return get_start_hexamer_freq(code, domain);
        case HexamerPosition::INTERNAL:
            return get_internal_hexamer_freq(code, domain);
        case HexamerPosition::END:
            return get_end_hexamer_freq(code, domain);
        default:
            return get_internal_hexamer_freq(code, domain);
    }
}

inline float get_positional_hexamer_freq(uint32_t code, HexamerPosition pos) {
    return get_positional_hexamer_freq(code, pos, get_active_domain());
}

// Score sequence with positional awareness
// Assumes sequence is in correct reading frame (frame 0)
inline float score_positional_hexamers(const std::string& seq, Domain domain) {
    if (seq.length() < 12) return 0.0f;  // Need at least 4 codons

    float score = 0.0f;
    const char* data = seq.c_str();
    size_t len = seq.length();
    size_t num_hexamers = 0;

    constexpr float RANDOM_FREQ = 1.0f / 4096.0f;

    // Score hexamers based on position
    for (size_t i = 0; i + 5 < len; i += 3) {
        uint32_t code = encode_hexamer(data + i);
        if (code == UINT32_MAX) continue;

        float freq;
        if (i < 6) {
            // START region (first 2 codons)
            freq = get_start_hexamer_freq(code, domain);
        } else if (i + 8 >= len) {
            // END region (last 2 codons)
            freq = get_end_hexamer_freq(code, domain);
        } else {
            // INTERNAL region
            freq = get_internal_hexamer_freq(code, domain);
        }

        if (freq > 0.0f) {
            score += std::log2(freq / RANDOM_FREQ);
            num_hexamers++;
        }
    }

    return num_hexamers > 0 ? score / num_hexamers : 0.0f;
}

inline float score_positional_hexamers(const std::string& seq) {
    return score_positional_hexamers(seq, get_active_domain());
}

// Compare positional scores across all domains (8 domains: GTDB + 7 eukaryotic)
inline MultiDomainResult score_positional_all_domains(const std::string& seq) {
    MultiDomainResult result;

    if (seq.length() < 12) return result;

    float scores[8];
    scores[0] = score_positional_hexamers(seq, Domain::GTDB);
    scores[1] = score_positional_hexamers(seq, Domain::FUNGI);
    scores[2] = score_positional_hexamers(seq, Domain::PROTOZOA);
    scores[3] = score_positional_hexamers(seq, Domain::INVERTEBRATE);
    scores[4] = score_positional_hexamers(seq, Domain::PLANT);
    scores[5] = score_positional_hexamers(seq, Domain::VERTEBRATE_MAMMALIAN);
    scores[6] = score_positional_hexamers(seq, Domain::VERTEBRATE_OTHER);
    scores[7] = score_positional_hexamers(seq, Domain::VIRAL);

    // Find best and convert to probabilities via softmax
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
        scores[i] = std::exp((scores[i] - max_score) * 0.5f);  // Temperature scaling
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

// Namespace alias for backwards compatibility with old API
namespace positional {

// Calculate hexamer score with damage awareness
// Used by frame_selector for damage-aware scoring
inline float calculate_hexamer_score_damage_aware(const std::string& seq, int frame, float delta_5, float delta_3) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    float log_prob_sum = 0.0f;
    int hexamer_count = 0;
    static constexpr float RANDOM_FREQ = 1.0f / 4096.0f;

    const char* data = seq.c_str();
    size_t len = seq.length();

    // For damage weighting:
    // delta_5: damage rate at 5' (affects start positions)
    // delta_3: damage rate at 3' (affects end positions)
    float damage_weight_5 = 1.0f - delta_5 * 0.3f;  // Reduce weight if high damage
    float damage_weight_3 = 1.0f - delta_3 * 0.3f;

    for (size_t i = frame; i + 5 < len; i += 3) {
        uint32_t code = encode_hexamer(data + i);
        if (code == UINT32_MAX) continue;

        float freq;
        float weight = 1.0f;

        if (i < static_cast<size_t>(frame + 6)) {
            // START region
            freq = get_start_hexamer_freq(code);
            weight = damage_weight_5;  // Weight by damage at 5'
        } else if (i + 8 >= len) {
            // END region
            freq = get_end_hexamer_freq(code);
            weight = damage_weight_3;  // Weight by damage at 3'
        } else {
            // INTERNAL region
            freq = get_internal_hexamer_freq(code);
        }

        if (freq > 0.0f) {
            log_prob_sum += weight * std::log2(freq / RANDOM_FREQ);
            hexamer_count++;
        }
    }

    if (hexamer_count == 0) return 0.0f;

    // Normalize and convert to 0-1 score
    float avg_log_ratio = log_prob_sum / hexamer_count;
    // Typical range is -2 to +4, map to 0-1
    float score = (avg_log_ratio + 2.0f) / 6.0f;
    return std::max(0.0f, std::min(1.0f, score));
}

} // namespace positional

// ============================================================================
// Ensemble Mode Support for META Domain
// ============================================================================

// Get ensemble-weighted START position hexamer frequency using domain probabilities
inline float get_ensemble_start_hexamer_freq(uint32_t code) {
    if (code >= 4096) return 0.0f;

    if (!ensemble_mode_enabled()) {
        return get_start_hexamer_freq(code, get_active_domain());
    }

    const auto& probs = current_domain_probs();
    return probs.gtdb_prob * GTDB_START_HEXAMER_FREQ[code] +
           probs.fungi_prob * FUNGI_START_HEXAMER_FREQ[code] +
           probs.protozoa_prob * PROTOZOA_START_HEXAMER_FREQ[code] +
           probs.invertebrate_prob * INVERTEBRATE_START_HEXAMER_FREQ[code] +
           probs.plant_prob * PLANT_START_HEXAMER_FREQ[code] +
           probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_START_HEXAMER_FREQ[code] +
           probs.vertebrate_other_prob * VERTEBRATE_OTHER_START_HEXAMER_FREQ[code] +
           probs.viral_prob * VIRAL_START_HEXAMER_FREQ[code];
}

// Get ensemble-weighted INTERNAL position hexamer frequency using domain probabilities
inline float get_ensemble_internal_hexamer_freq(uint32_t code) {
    if (code >= 4096) return 0.0f;

    if (!ensemble_mode_enabled()) {
        return get_internal_hexamer_freq(code, get_active_domain());
    }

    const auto& probs = current_domain_probs();
    return probs.gtdb_prob * GTDB_INTERNAL_HEXAMER_FREQ[code] +
           probs.fungi_prob * FUNGI_INTERNAL_HEXAMER_FREQ[code] +
           probs.protozoa_prob * PROTOZOA_INTERNAL_HEXAMER_FREQ[code] +
           probs.invertebrate_prob * INVERTEBRATE_INTERNAL_HEXAMER_FREQ[code] +
           probs.plant_prob * PLANT_INTERNAL_HEXAMER_FREQ[code] +
           probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_INTERNAL_HEXAMER_FREQ[code] +
           probs.vertebrate_other_prob * VERTEBRATE_OTHER_INTERNAL_HEXAMER_FREQ[code] +
           probs.viral_prob * VIRAL_INTERNAL_HEXAMER_FREQ[code];
}

// Get ensemble-weighted END position hexamer frequency using domain probabilities
inline float get_ensemble_end_hexamer_freq(uint32_t code) {
    if (code >= 4096) return 0.0f;

    if (!ensemble_mode_enabled()) {
        return get_end_hexamer_freq(code, get_active_domain());
    }

    const auto& probs = current_domain_probs();
    return probs.gtdb_prob * GTDB_END_HEXAMER_FREQ[code] +
           probs.fungi_prob * FUNGI_END_HEXAMER_FREQ[code] +
           probs.protozoa_prob * PROTOZOA_END_HEXAMER_FREQ[code] +
           probs.invertebrate_prob * INVERTEBRATE_END_HEXAMER_FREQ[code] +
           probs.plant_prob * PLANT_END_HEXAMER_FREQ[code] +
           probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_END_HEXAMER_FREQ[code] +
           probs.vertebrate_other_prob * VERTEBRATE_OTHER_END_HEXAMER_FREQ[code] +
           probs.viral_prob * VIRAL_END_HEXAMER_FREQ[code];
}

// Score positional hexamers using ensemble mode with domain probabilities
inline float score_positional_hexamers_weighted(const std::string& seq, const MultiDomainResult& probs) {
    if (seq.length() < 12) return 0.0f;

    float score = 0.0f;
    const char* data = seq.c_str();
    size_t len = seq.length();
    size_t num_hexamers = 0;

    constexpr float RANDOM_FREQ = 1.0f / 4096.0f;

    // Score hexamers based on position with weighted frequencies
    for (size_t i = 0; i + 5 < len; i += 3) {
        uint32_t code = encode_hexamer(data + i);
        if (code == UINT32_MAX) continue;

        float freq;
        if (i < 6) {
            // START region (first 2 codons) - weighted
            freq = probs.gtdb_prob * GTDB_START_HEXAMER_FREQ[code] +
                   probs.fungi_prob * FUNGI_START_HEXAMER_FREQ[code] +
                   probs.protozoa_prob * PROTOZOA_START_HEXAMER_FREQ[code] +
                   probs.invertebrate_prob * INVERTEBRATE_START_HEXAMER_FREQ[code] +
                   probs.plant_prob * PLANT_START_HEXAMER_FREQ[code] +
                   probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_START_HEXAMER_FREQ[code] +
                   probs.vertebrate_other_prob * VERTEBRATE_OTHER_START_HEXAMER_FREQ[code] +
                   probs.viral_prob * VIRAL_START_HEXAMER_FREQ[code];
        } else if (i + 8 >= len) {
            // END region (last 2 codons) - weighted
            freq = probs.gtdb_prob * GTDB_END_HEXAMER_FREQ[code] +
                   probs.fungi_prob * FUNGI_END_HEXAMER_FREQ[code] +
                   probs.protozoa_prob * PROTOZOA_END_HEXAMER_FREQ[code] +
                   probs.invertebrate_prob * INVERTEBRATE_END_HEXAMER_FREQ[code] +
                   probs.plant_prob * PLANT_END_HEXAMER_FREQ[code] +
                   probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_END_HEXAMER_FREQ[code] +
                   probs.vertebrate_other_prob * VERTEBRATE_OTHER_END_HEXAMER_FREQ[code] +
                   probs.viral_prob * VIRAL_END_HEXAMER_FREQ[code];
        } else {
            // INTERNAL region - weighted
            freq = probs.gtdb_prob * GTDB_INTERNAL_HEXAMER_FREQ[code] +
                   probs.fungi_prob * FUNGI_INTERNAL_HEXAMER_FREQ[code] +
                   probs.protozoa_prob * PROTOZOA_INTERNAL_HEXAMER_FREQ[code] +
                   probs.invertebrate_prob * INVERTEBRATE_INTERNAL_HEXAMER_FREQ[code] +
                   probs.plant_prob * PLANT_INTERNAL_HEXAMER_FREQ[code] +
                   probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_INTERNAL_HEXAMER_FREQ[code] +
                   probs.vertebrate_other_prob * VERTEBRATE_OTHER_INTERNAL_HEXAMER_FREQ[code] +
                   probs.viral_prob * VIRAL_INTERNAL_HEXAMER_FREQ[code];
        }

        if (freq > 0.0f) {
            score += std::log2(freq / RANDOM_FREQ);
            num_hexamers++;
        }
    }

    return num_hexamers > 0 ? score / num_hexamers : 0.0f;
}

} // namespace agp
