#pragma once

/**
 * @file damage_likelihood.hpp
 * @brief Unified multi-domain damage likelihood interface
 *
 * Provides domain selection and damage likelihood lookup across all trained domains:
 * fungi, protozoa, invertebrate, plant, vertebrate_mammalian, vertebrate_other, viral
 *
 * Damage likelihood ratios (LLR) measure how likely a hexamer pattern is to have
 * arisen from ancient DNA damage (C→T at 5' end, G→A at 3' end).
 */

#include <string>
#include <cstdint>
#include <cmath>

// Include hexamer_tables.hpp for Domain enum and encode_hexamer
#include "agp/hexamer_tables.hpp"

// Include all domain-specific damage likelihood tables
#include "agp/gtdb_damage_likelihood.hpp"
#include "agp/fungi_damage_likelihood.hpp"
#include "agp/protozoa_damage_likelihood.hpp"
#include "agp/invertebrate_damage_likelihood.hpp"
#include "agp/plant_damage_likelihood.hpp"
#include "agp/vertebrate_mammalian_damage_likelihood.hpp"
#include "agp/vertebrate_other_damage_likelihood.hpp"
#include "agp/viral_damage_likelihood.hpp"

namespace agp {

// Get 5' damage log-likelihood ratio for specific domain
// Positive values indicate damage-consistent patterns (C→T at 5' end)
inline float get_damage_llr_5prime(uint32_t code, Domain domain) {
    if (code >= 4096) return 0.0f;

    switch (domain) {
        case Domain::GTDB:
            return GTDB_DAMAGE_LLR_5PRIME[code];
        case Domain::FUNGI:
            return FUNGI_DAMAGE_LLR_5PRIME[code];
        case Domain::PROTOZOA:
            return PROTOZOA_DAMAGE_LLR_5PRIME[code];
        case Domain::INVERTEBRATE:
            return INVERTEBRATE_DAMAGE_LLR_5PRIME[code];
        case Domain::PLANT:
            return PLANT_DAMAGE_LLR_5PRIME[code];
        case Domain::VERTEBRATE_MAMMALIAN:
            return VERTEBRATE_MAMMALIAN_DAMAGE_LLR_5PRIME[code];
        case Domain::VERTEBRATE_OTHER:
            return VERTEBRATE_OTHER_DAMAGE_LLR_5PRIME[code];
        case Domain::VIRAL:
            return VIRAL_DAMAGE_LLR_5PRIME[code];
        case Domain::META:
        default:
            // For META mode, use GTDB as default (bacteria-focused)
            // Ensemble weighting happens via get_ensemble_damage_llr_5prime()
            return GTDB_DAMAGE_LLR_5PRIME[code];
    }
}

// Get 3' damage log-likelihood ratio for specific domain
// Positive values indicate damage-consistent patterns (G→A at 3' end)
inline float get_damage_llr_3prime(uint32_t code, Domain domain) {
    if (code >= 4096) return 0.0f;

    switch (domain) {
        case Domain::GTDB:
            return GTDB_DAMAGE_LLR_3PRIME[code];
        case Domain::FUNGI:
            return FUNGI_DAMAGE_LLR_3PRIME[code];
        case Domain::PROTOZOA:
            return PROTOZOA_DAMAGE_LLR_3PRIME[code];
        case Domain::INVERTEBRATE:
            return INVERTEBRATE_DAMAGE_LLR_3PRIME[code];
        case Domain::PLANT:
            return PLANT_DAMAGE_LLR_3PRIME[code];
        case Domain::VERTEBRATE_MAMMALIAN:
            return VERTEBRATE_MAMMALIAN_DAMAGE_LLR_3PRIME[code];
        case Domain::VERTEBRATE_OTHER:
            return VERTEBRATE_OTHER_DAMAGE_LLR_3PRIME[code];
        case Domain::VIRAL:
            return VIRAL_DAMAGE_LLR_3PRIME[code];
        case Domain::META:
        default:
            // For META mode, use GTDB as default (bacteria-focused)
            // Ensemble weighting happens via get_ensemble_damage_llr_3prime()
            return GTDB_DAMAGE_LLR_3PRIME[code];
    }
}

// Get damage LLR using active domain
inline float get_damage_llr_5prime(uint32_t code) {
    return get_damage_llr_5prime(code, get_active_domain());
}

inline float get_damage_llr_3prime(uint32_t code) {
    return get_damage_llr_3prime(code, get_active_domain());
}

// Get damage LLR from sequence string
inline float get_damage_llr_5prime(const char* seq, Domain domain) {
    uint32_t code = encode_hexamer(seq);
    if (code == UINT32_MAX) return 0.0f;
    return get_damage_llr_5prime(code, domain);
}

inline float get_damage_llr_3prime(const char* seq, Domain domain) {
    uint32_t code = encode_hexamer(seq);
    if (code == UINT32_MAX) return 0.0f;
    return get_damage_llr_3prime(code, domain);
}

inline float get_damage_llr_5prime(const char* seq) {
    return get_damage_llr_5prime(seq, get_active_domain());
}

inline float get_damage_llr_3prime(const char* seq) {
    return get_damage_llr_3prime(seq, get_active_domain());
}

// Calculate damage probability from sequence hexamers
// Uses hexamer context to estimate likelihood of damage
// position: 0 = 5' terminal, seq_len-1 = 3' terminal
inline float get_hexamer_damage_prob(const char* seq, size_t position, size_t seq_len, Domain domain) {
    if (seq_len < 6) return 0.0f;

    uint32_t code = encode_hexamer(seq);
    if (code == UINT32_MAX) return 0.0f;

    // Use 5' LLR for positions near 5' end, 3' LLR for positions near 3' end
    float llr;
    if (position < seq_len / 2) {
        // 5' half of the sequence
        llr = get_damage_llr_5prime(code, domain);
        // Weight by distance from 5' terminus (exponential decay)
        float distance_weight = std::exp(-static_cast<float>(position) / 10.0f);
        llr *= distance_weight;
    } else {
        // 3' half of the sequence
        llr = get_damage_llr_3prime(code, domain);
        // Weight by distance from 3' terminus
        float distance_from_3prime = static_cast<float>(seq_len - 1 - position);
        float distance_weight = std::exp(-distance_from_3prime / 10.0f);
        llr *= distance_weight;
    }

    // Convert log-likelihood ratio to probability via sigmoid
    // P(damage) = 1 / (1 + exp(-LLR))
    return 1.0f / (1.0f + std::exp(-llr));
}

inline float get_hexamer_damage_prob(const char* seq, size_t position, size_t seq_len) {
    return get_hexamer_damage_prob(seq, position, seq_len, get_active_domain());
}

// Simple hexamer damage probability lookup (1-argument version)
// Uses 5' LLR for terminal position estimate
// This is a simplified version for quick damage checks
inline float get_hexamer_damage_prob(const char* seq) {
    uint32_t code = encode_hexamer(seq);
    if (code == UINT32_MAX) return 0.0f;

    // Use 5' damage LLR as primary signal (most common damage type)
    float llr = get_damage_llr_5prime(code, get_active_domain());

    // Convert LLR to probability via sigmoid
    return 1.0f / (1.0f + std::exp(-llr));
}

// Score entire sequence for damage patterns
// Returns average damage probability across terminal regions
struct DamageScore {
    float five_prime_score = 0.0f;   // Average damage prob at 5' end
    float three_prime_score = 0.0f;  // Average damage prob at 3' end
    float combined_score = 0.0f;     // Combined damage assessment
    int five_prime_hexamers = 0;
    int three_prime_hexamers = 0;
};

inline DamageScore score_sequence_damage(const std::string& seq, Domain domain) {
    DamageScore result;

    if (seq.length() < 12) return result;

    const char* data = seq.c_str();
    size_t len = seq.length();

    // Score 5' terminal region (first ~10% or 30 bp, whichever is smaller)
    size_t five_prime_region = std::min(len / 10, size_t(30));
    for (size_t i = 0; i + 5 < five_prime_region; ++i) {
        uint32_t code = encode_hexamer(data + i);
        if (code == UINT32_MAX) continue;

        float llr = get_damage_llr_5prime(code, domain);
        result.five_prime_score += llr;
        result.five_prime_hexamers++;
    }

    // Score 3' terminal region
    size_t three_prime_start = len - std::min(len / 10, size_t(30));
    for (size_t i = three_prime_start; i + 5 < len; ++i) {
        uint32_t code = encode_hexamer(data + i);
        if (code == UINT32_MAX) continue;

        float llr = get_damage_llr_3prime(code, domain);
        result.three_prime_score += llr;
        result.three_prime_hexamers++;
    }

    // Normalize scores
    if (result.five_prime_hexamers > 0) {
        result.five_prime_score /= result.five_prime_hexamers;
    }
    if (result.three_prime_hexamers > 0) {
        result.three_prime_score /= result.three_prime_hexamers;
    }

    // Combined score: average of both termini
    int total_hexamers = result.five_prime_hexamers + result.three_prime_hexamers;
    if (total_hexamers > 0) {
        result.combined_score = (result.five_prime_score + result.three_prime_score) / 2.0f;
    }

    return result;
}

inline DamageScore score_sequence_damage(const std::string& seq) {
    return score_sequence_damage(seq, get_active_domain());
}

// ============================================================================
// Ensemble Mode Support for META Domain
// ============================================================================

// Get ensemble-weighted 5' damage LLR using domain probabilities
// Uses current_domain_probs() set by score_all_domains() in hexamer_tables.hpp
inline float get_ensemble_damage_llr_5prime(uint32_t code) {
    if (code >= 4096) return 0.0f;

    if (!ensemble_mode_enabled()) {
        return get_damage_llr_5prime(code, get_active_domain());
    }

    const auto& probs = current_domain_probs();
    return probs.gtdb_prob * GTDB_DAMAGE_LLR_5PRIME[code] +
           probs.fungi_prob * FUNGI_DAMAGE_LLR_5PRIME[code] +
           probs.protozoa_prob * PROTOZOA_DAMAGE_LLR_5PRIME[code] +
           probs.invertebrate_prob * INVERTEBRATE_DAMAGE_LLR_5PRIME[code] +
           probs.plant_prob * PLANT_DAMAGE_LLR_5PRIME[code] +
           probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_DAMAGE_LLR_5PRIME[code] +
           probs.vertebrate_other_prob * VERTEBRATE_OTHER_DAMAGE_LLR_5PRIME[code] +
           probs.viral_prob * VIRAL_DAMAGE_LLR_5PRIME[code];
}

// Get ensemble-weighted 3' damage LLR using domain probabilities
inline float get_ensemble_damage_llr_3prime(uint32_t code) {
    if (code >= 4096) return 0.0f;

    if (!ensemble_mode_enabled()) {
        return get_damage_llr_3prime(code, get_active_domain());
    }

    const auto& probs = current_domain_probs();
    return probs.gtdb_prob * GTDB_DAMAGE_LLR_3PRIME[code] +
           probs.fungi_prob * FUNGI_DAMAGE_LLR_3PRIME[code] +
           probs.protozoa_prob * PROTOZOA_DAMAGE_LLR_3PRIME[code] +
           probs.invertebrate_prob * INVERTEBRATE_DAMAGE_LLR_3PRIME[code] +
           probs.plant_prob * PLANT_DAMAGE_LLR_3PRIME[code] +
           probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_DAMAGE_LLR_3PRIME[code] +
           probs.vertebrate_other_prob * VERTEBRATE_OTHER_DAMAGE_LLR_3PRIME[code] +
           probs.viral_prob * VIRAL_DAMAGE_LLR_3PRIME[code];
}

// Score sequence damage using ensemble mode with domain probabilities
inline DamageScore score_sequence_damage_weighted(const std::string& seq, const MultiDomainResult& probs) {
    DamageScore result;

    if (seq.length() < 12) return result;

    const char* data = seq.c_str();
    size_t len = seq.length();

    // Score 5' terminal region
    size_t five_prime_region = std::min(len / 10, size_t(30));
    for (size_t i = 0; i + 5 < five_prime_region; ++i) {
        uint32_t code = encode_hexamer(data + i);
        if (code == UINT32_MAX) continue;

        // Weighted LLR across all domains
        float llr = probs.gtdb_prob * GTDB_DAMAGE_LLR_5PRIME[code] +
                    probs.fungi_prob * FUNGI_DAMAGE_LLR_5PRIME[code] +
                    probs.protozoa_prob * PROTOZOA_DAMAGE_LLR_5PRIME[code] +
                    probs.invertebrate_prob * INVERTEBRATE_DAMAGE_LLR_5PRIME[code] +
                    probs.plant_prob * PLANT_DAMAGE_LLR_5PRIME[code] +
                    probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_DAMAGE_LLR_5PRIME[code] +
                    probs.vertebrate_other_prob * VERTEBRATE_OTHER_DAMAGE_LLR_5PRIME[code] +
                    probs.viral_prob * VIRAL_DAMAGE_LLR_5PRIME[code];

        result.five_prime_score += llr;
        result.five_prime_hexamers++;
    }

    // Score 3' terminal region
    size_t three_prime_start = len - std::min(len / 10, size_t(30));
    for (size_t i = three_prime_start; i + 5 < len; ++i) {
        uint32_t code = encode_hexamer(data + i);
        if (code == UINT32_MAX) continue;

        // Weighted LLR across all domains
        float llr = probs.gtdb_prob * GTDB_DAMAGE_LLR_3PRIME[code] +
                    probs.fungi_prob * FUNGI_DAMAGE_LLR_3PRIME[code] +
                    probs.protozoa_prob * PROTOZOA_DAMAGE_LLR_3PRIME[code] +
                    probs.invertebrate_prob * INVERTEBRATE_DAMAGE_LLR_3PRIME[code] +
                    probs.plant_prob * PLANT_DAMAGE_LLR_3PRIME[code] +
                    probs.vertebrate_mammalian_prob * VERTEBRATE_MAMMALIAN_DAMAGE_LLR_3PRIME[code] +
                    probs.vertebrate_other_prob * VERTEBRATE_OTHER_DAMAGE_LLR_3PRIME[code] +
                    probs.viral_prob * VIRAL_DAMAGE_LLR_3PRIME[code];

        result.three_prime_score += llr;
        result.three_prime_hexamers++;
    }

    // Normalize scores
    if (result.five_prime_hexamers > 0) {
        result.five_prime_score /= result.five_prime_hexamers;
    }
    if (result.three_prime_hexamers > 0) {
        result.three_prime_score /= result.three_prime_hexamers;
    }

    // Combined score
    int total_hexamers = result.five_prime_hexamers + result.three_prime_hexamers;
    if (total_hexamers > 0) {
        result.combined_score = (result.five_prime_score + result.three_prime_score) / 2.0f;
    }

    return result;
}

} // namespace agp
