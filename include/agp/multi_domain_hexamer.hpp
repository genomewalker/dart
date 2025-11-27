#pragma once

// Multi-domain hexamer classifier for universal gene prediction
// Supports scoring sequences against multiple taxonomic domains
// and identifying the most likely source domain.
//
// Available domains:
// - GTDB (bacteria/archaea) - 134.7B hexamers from GTDB r220
// - FUNGI - 3.2B hexamers from RefSeq fungi CDS
// - PLANT - 4.3B hexamers from RefSeq plant CDS
// - PROTOZOA - 629M hexamers from RefSeq protozoa CDS
// - INVERTEBRATE - 8.5B hexamers from RefSeq invertebrate CDS
// - VIRAL - 165M hexamers from RefSeq viral CDS
// - VERTEBRATE_MAMMALIAN - 8.7B hexamers from RefSeq vertebrate mammalian CDS
// - VERTEBRATE_OTHER - 15.3B hexamers from RefSeq vertebrate other CDS

#include <cstdint>
#include <array>
#include <string>
#include <algorithm>
#include <cmath>

// Include all domain-specific hexamer tables
#include "agp/gtdb_hexamer_table.hpp"
#include "agp/fungi_hexamer_table.hpp"
#include "agp/plant_hexamer_table.hpp"
#include "agp/protozoa_hexamer_table.hpp"
#include "agp/invertebrate_hexamer_table.hpp"
#include "agp/viral_hexamer_table.hpp"
#include "agp/vertebrate_mammalian_hexamer_table.hpp"
#include "agp/vertebrate_other_hexamer_table.hpp"

namespace agp {

// Taxonomic domains for gene prediction
enum class Domain : uint8_t {
    GTDB = 0,              // Bacteria/Archaea (GTDB r220)
    FUNGI,                 // Fungi
    PLANT,                 // Plants
    PROTOZOA,              // Protozoa
    INVERTEBRATE,          // Invertebrates
    VIRAL,                 // Viruses
    VERTEBRATE_MAMMALIAN,  // Vertebrate mammals
    VERTEBRATE_OTHER,      // Other vertebrates (birds, fish, reptiles, amphibians)
    NUM_DOMAINS            // Number of domains (for array sizing)
};

// Domain names for output
inline const char* domain_name(Domain d) {
    switch (d) {
        case Domain::GTDB: return "bacteria_archaea";
        case Domain::FUNGI: return "fungi";
        case Domain::PLANT: return "plant";
        case Domain::PROTOZOA: return "protozoa";
        case Domain::INVERTEBRATE: return "invertebrate";
        case Domain::VIRAL: return "viral";
        case Domain::VERTEBRATE_MAMMALIAN: return "vertebrate_mammalian";
        case Domain::VERTEBRATE_OTHER: return "vertebrate_other";
        default: return "unknown";
    }
}

// Parse domain from string (case-insensitive)
inline Domain parse_domain(const std::string& name) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    if (lower == "gtdb" || lower == "bacteria" || lower == "archaea" ||
        lower == "bacteria_archaea" || lower == "prokaryote") {
        return Domain::GTDB;
    }
    if (lower == "fungi" || lower == "fungal") {
        return Domain::FUNGI;
    }
    if (lower == "plant" || lower == "plants" || lower == "viridiplantae") {
        return Domain::PLANT;
    }
    if (lower == "protozoa" || lower == "protist" || lower == "protists") {
        return Domain::PROTOZOA;
    }
    if (lower == "invertebrate" || lower == "invertebrates") {
        return Domain::INVERTEBRATE;
    }
    if (lower == "viral" || lower == "virus" || lower == "viruses") {
        return Domain::VIRAL;
    }
    if (lower == "vertebrate_mammalian" || lower == "mammal" || lower == "mammals" || lower == "mammalian") {
        return Domain::VERTEBRATE_MAMMALIAN;
    }
    if (lower == "vertebrate_other" || lower == "vertebrate" || lower == "vertebrates") {
        return Domain::VERTEBRATE_OTHER;
    }

    // Default to GTDB (bacteria/archaea) for ancient DNA
    return Domain::GTDB;
}

// Get hexamer frequency for a specific domain
inline float get_domain_hexamer_freq(Domain domain, uint32_t hex_code) {
    if (hex_code >= 4096) return 0.0f;

    switch (domain) {
        case Domain::GTDB:
            return GTDB_HEXAMER_FREQ[hex_code];
        case Domain::FUNGI:
            return FUNGI_HEXAMER_FREQ[hex_code];
        case Domain::PLANT:
            return PLANT_HEXAMER_FREQ[hex_code];
        case Domain::PROTOZOA:
            return PROTOZOA_HEXAMER_FREQ[hex_code];
        case Domain::INVERTEBRATE:
            return INVERTEBRATE_HEXAMER_FREQ[hex_code];
        case Domain::VIRAL:
            return VIRAL_HEXAMER_FREQ[hex_code];
        case Domain::VERTEBRATE_MAMMALIAN:
            return VERTEBRATE_MAMMALIAN_HEXAMER_FREQ[hex_code];
        case Domain::VERTEBRATE_OTHER:
            return VERTEBRATE_OTHER_HEXAMER_FREQ[hex_code];
        default:
            return GTDB_HEXAMER_FREQ[hex_code];
    }
}

// Get hexamer frequency for a specific domain (string version)
inline float get_domain_hexamer_freq(Domain domain, const char* hexamer) {
    uint32_t code = encode_hexamer(hexamer);
    if (code == UINT32_MAX) return 0.0f;
    return get_domain_hexamer_freq(domain, code);
}

// Result of domain classification
struct DomainScore {
    Domain domain;
    float log_likelihood;  // Log-likelihood score
    float probability;     // Normalized probability (0-1)
};

// Result of multi-domain scoring
struct MultiDomainResult {
    std::array<DomainScore, static_cast<size_t>(Domain::NUM_DOMAINS)> scores;
    Domain best_domain;
    float confidence;      // Confidence in best domain (0-1)
    int hexamer_count;     // Number of valid hexamers scored
};

// Fast uppercase conversion
inline char fast_upper_md(char c) {
    return (c >= 'a' && c <= 'z') ? c - 32 : c;
}

// Calculate log-likelihood ratio for a sequence against a specific domain
// Returns average log(P_domain / P_uniform) per hexamer
inline float calculate_domain_llr(const std::string& seq, int frame, Domain domain) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    float log_prob_sum = 0.0f;
    int hexamer_count = 0;

    // Background frequency (uniform = 1/4096)
    static constexpr float LOG_BACKGROUND = -8.317766f;  // log(1/4096)

    char hexbuf[7];
    hexbuf[6] = '\0';

    for (size_t i = frame; i + 5 < seq.length(); i += 3) {
        hexbuf[0] = fast_upper_md(seq[i]);
        hexbuf[1] = fast_upper_md(seq[i + 1]);
        hexbuf[2] = fast_upper_md(seq[i + 2]);
        hexbuf[3] = fast_upper_md(seq[i + 3]);
        hexbuf[4] = fast_upper_md(seq[i + 4]);
        hexbuf[5] = fast_upper_md(seq[i + 5]);

        // Skip if any N
        if (hexbuf[0] == 'N' || hexbuf[1] == 'N' || hexbuf[2] == 'N' ||
            hexbuf[3] == 'N' || hexbuf[4] == 'N' || hexbuf[5] == 'N') continue;

        uint32_t code = encode_hexamer(hexbuf);
        if (code == UINT32_MAX) continue;

        float freq = get_domain_hexamer_freq(domain, code);

        // Pseudocount for zero frequencies
        if (freq < 1e-9f) {
            freq = 1e-8f;  // Small pseudocount
        }

        // Log-likelihood ratio
        log_prob_sum += std::log(freq) - LOG_BACKGROUND;
        hexamer_count++;
    }

    if (hexamer_count == 0) return 0.0f;
    return log_prob_sum / hexamer_count;
}

// Score a sequence against all domains
// Returns MultiDomainResult with scores for each domain and classification
inline MultiDomainResult score_all_domains(const std::string& seq, int frame) {
    MultiDomainResult result;
    result.hexamer_count = 0;

    // Calculate LLR for each domain
    float max_llr = -1000.0f;
    Domain best = Domain::GTDB;

    for (size_t i = 0; i < static_cast<size_t>(Domain::NUM_DOMAINS); ++i) {
        Domain d = static_cast<Domain>(i);
        float llr = calculate_domain_llr(seq, frame, d);
        result.scores[i].domain = d;
        result.scores[i].log_likelihood = llr;

        if (llr > max_llr) {
            max_llr = llr;
            best = d;
        }
    }

    result.best_domain = best;

    // Convert to probabilities using softmax
    // First, shift by max for numerical stability
    float sum_exp = 0.0f;
    for (size_t i = 0; i < static_cast<size_t>(Domain::NUM_DOMAINS); ++i) {
        float shifted = result.scores[i].log_likelihood - max_llr;
        // Clamp to avoid overflow
        shifted = std::max(-20.0f, shifted);
        float exp_val = std::exp(shifted);
        result.scores[i].probability = exp_val;
        sum_exp += exp_val;
    }

    // Normalize
    if (sum_exp > 0.0f) {
        for (size_t i = 0; i < static_cast<size_t>(Domain::NUM_DOMAINS); ++i) {
            result.scores[i].probability /= sum_exp;
        }
    }

    // Confidence is the probability of the best domain
    result.confidence = result.scores[static_cast<size_t>(best)].probability;

    // Count hexamers (same for all domains)
    if (seq.length() >= static_cast<size_t>(frame + 6)) {
        for (size_t i = frame; i + 5 < seq.length(); i += 3) {
            char c0 = fast_upper_md(seq[i]);
            char c1 = fast_upper_md(seq[i+1]);
            char c2 = fast_upper_md(seq[i+2]);
            char c3 = fast_upper_md(seq[i+3]);
            char c4 = fast_upper_md(seq[i+4]);
            char c5 = fast_upper_md(seq[i+5]);
            if (c0 != 'N' && c1 != 'N' && c2 != 'N' &&
                c3 != 'N' && c4 != 'N' && c5 != 'N') {
                result.hexamer_count++;
            }
        }
    }

    return result;
}

// Calculate dicodon score using a specific domain's hexamer table
// This is the domain-aware version of calculate_dicodon_score()
inline float calculate_dicodon_score_domain(const std::string& seq, int frame, Domain domain) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    float log_prob_sum = 0.0f;
    int hexamer_count = 0;

    // Background frequency (uniform = 1/4096)
    static constexpr float BACKGROUND_FREQ = 1.0f / 4096.0f;

    char hexbuf[7];
    hexbuf[6] = '\0';

    for (size_t i = frame; i + 5 < seq.length(); i += 3) {
        hexbuf[0] = fast_upper_md(seq[i]);
        hexbuf[1] = fast_upper_md(seq[i + 1]);
        hexbuf[2] = fast_upper_md(seq[i + 2]);
        hexbuf[3] = fast_upper_md(seq[i + 3]);
        hexbuf[4] = fast_upper_md(seq[i + 4]);
        hexbuf[5] = fast_upper_md(seq[i + 5]);

        // Skip if any N
        if (hexbuf[0] == 'N' || hexbuf[1] == 'N' || hexbuf[2] == 'N' ||
            hexbuf[3] == 'N' || hexbuf[4] == 'N' || hexbuf[5] == 'N') continue;

        float freq = get_domain_hexamer_freq(domain, hexbuf);

        // Pseudocount for zero frequencies
        if (freq < 1e-9f) {
            freq = BACKGROUND_FREQ * 0.01f;
        }

        // Log-likelihood ratio: log(P_coding / P_background)
        log_prob_sum += std::log(freq / BACKGROUND_FREQ);
        hexamer_count++;
    }

    if (hexamer_count == 0) return 0.5f;

    // Normalize to 0-1 range
    float avg_llr = log_prob_sum / hexamer_count;
    float score = 0.5f + 0.15f * avg_llr;

    return std::max(0.0f, std::min(1.0f, score));
}

// Ensemble scoring: average scores across related domains
// Useful when the exact source is unknown but kingdom is known
inline float calculate_dicodon_score_ensemble(const std::string& seq, int frame,
                                               const std::vector<Domain>& domains) {
    if (domains.empty()) return 0.5f;

    float total_score = 0.0f;
    for (Domain d : domains) {
        total_score += calculate_dicodon_score_domain(seq, frame, d);
    }
    return total_score / domains.size();
}

// Probability-weighted ensemble scoring for metagenomes
// Uses domain probabilities from score_all_domains() to weight each domain's contribution
// This is robust for short reads where single-domain classification is unreliable
inline float calculate_dicodon_score_weighted(const std::string& seq, int frame,
                                               const MultiDomainResult& domain_probs) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.5f;

    float weighted_score = 0.0f;
    for (size_t i = 0; i < static_cast<size_t>(Domain::NUM_DOMAINS); ++i) {
        float prob = domain_probs.scores[i].probability;
        if (prob > 0.01f) {  // Skip negligible contributions
            float domain_score = calculate_dicodon_score_domain(seq, frame, static_cast<Domain>(i));
            weighted_score += prob * domain_score;
        }
    }
    return weighted_score;
}

// Thread-local ensemble mode flag for metagenome scoring
inline bool& get_ensemble_mode() {
    thread_local bool ensemble_mode = false;
    return ensemble_mode;
}

inline void set_ensemble_mode(bool enabled) {
    get_ensemble_mode() = enabled;
}

// Thread-local domain probabilities for weighted scoring
inline MultiDomainResult& get_domain_probs() {
    thread_local MultiDomainResult probs;
    return probs;
}

inline void set_domain_probs(const MultiDomainResult& probs) {
    get_domain_probs() = probs;
}

// Predefined domain groups for ensemble scoring
inline std::vector<Domain> get_prokaryote_domains() {
    return {Domain::GTDB};
}

inline std::vector<Domain> get_eukaryote_domains() {
    return {Domain::FUNGI, Domain::PLANT, Domain::PROTOZOA,
            Domain::INVERTEBRATE, Domain::VERTEBRATE_MAMMALIAN, Domain::VERTEBRATE_OTHER};
}

inline std::vector<Domain> get_animal_domains() {
    return {Domain::INVERTEBRATE, Domain::VERTEBRATE_MAMMALIAN, Domain::VERTEBRATE_OTHER};
}

inline std::vector<Domain> get_vertebrate_domains() {
    return {Domain::VERTEBRATE_MAMMALIAN, Domain::VERTEBRATE_OTHER};
}

inline std::vector<Domain> get_all_domains() {
    std::vector<Domain> all;
    for (size_t i = 0; i < static_cast<size_t>(Domain::NUM_DOMAINS); ++i) {
        all.push_back(static_cast<Domain>(i));
    }
    return all;
}

// Global domain setting for scoring (defaults to GTDB for ancient DNA)
// Uses thread-local storage for thread-safe per-read domain selection
inline Domain& get_active_domain() {
    thread_local Domain active = Domain::GTDB;
    return active;
}

inline void set_active_domain(Domain d) {
    get_active_domain() = d;
}

// Set global default domain (used at startup, affects all threads' initial value)
// Note: This only affects threads that haven't yet accessed get_active_domain()
inline Domain& get_global_default_domain() {
    static Domain global_default = Domain::GTDB;
    return global_default;
}

inline void set_global_default_domain(Domain d) {
    get_global_default_domain() = d;
}

// Wrapper that uses the active domain
inline float get_hexamer_freq_active(const char* hexamer) {
    return get_domain_hexamer_freq(get_active_domain(), hexamer);
}

} // namespace agp
