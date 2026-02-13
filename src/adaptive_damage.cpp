/**
 * @file adaptive_damage.cpp
 * @brief Implementation of adaptive damage calibration
 */

#include "agp/adaptive_damage.hpp"
#include "agp/hexamer_tables.hpp"
#include <cmath>
#include <algorithm>

namespace agp {

// Check if a codon is a stop codon
static inline bool is_stop_codon(const char* codon) {
    char c1 = std::toupper(static_cast<unsigned char>(codon[0]));
    char c2 = std::toupper(static_cast<unsigned char>(codon[1]));
    char c3 = std::toupper(static_cast<unsigned char>(codon[2]));

    return (c1 == 'T' && c2 == 'A' && c3 == 'A') ||  // TAA
           (c1 == 'T' && c2 == 'A' && c3 == 'G') ||  // TAG
           (c1 == 'T' && c2 == 'G' && c3 == 'A');    // TGA
}

// Core implementation for ensemble hexamer LLR calculation (const char* version)
float AdaptiveDamageCalibrator::get_ensemble_hexamer_llr(
    const char* seq, size_t len, int frame) const {

    if (len < 6) return 0.0f;

    // Compute hexamer LLRs for each domain
    // Uses same formula as score_all_domains(): log2(freq * 4096 + 1e-10)
    // This ensures consistent scaling for ensemble weight learning
    std::array<float, 8> domain_llrs = {};
    static const std::array<Domain, 8> domains = {
        Domain::GTDB, Domain::FUNGI, Domain::PROTOZOA, Domain::INVERTEBRATE,
        Domain::PLANT, Domain::VERTEBRATE_MAMMALIAN, Domain::VERTEBRATE_OTHER, Domain::VIRAL
    };

    size_t start = (frame >= 0) ? static_cast<size_t>(frame) : 0;
    size_t n_hexamers = 0;

    for (size_t i = start; i + 6 <= len; i += 3) {
        uint32_t code = encode_hexamer(seq + i);
        if (code >= 4096) continue;

        for (size_t d = 0; d < 8; ++d) {
            float freq = get_hexamer_freq(code, domains[d]);
            // Same formula as score_all_domains: log2(freq * 4096 + 1e-10)
            domain_llrs[d] += std::log2(freq * 4096.0f + 1e-10f);
        }
        n_hexamers++;
    }

    if (n_hexamers == 0) return 0.0f;

    // Normalize by hexamer count
    for (auto& llr : domain_llrs) {
        llr /= static_cast<float>(n_hexamers);
    }

    // Return ensemble-weighted LLR
    return ensemble_weights_.get_weighted_llr(domain_llrs);
}

// std::string wrapper - delegates to const char* version
float AdaptiveDamageCalibrator::get_ensemble_hexamer_llr(
    const std::string& seq, int frame) const {
    return get_ensemble_hexamer_llr(seq.data(), seq.length(), frame);
}

AdaptiveDamageCalibrator::CorrectionValidation
AdaptiveDamageCalibrator::validate_correction(
    const std::string& original_seq,
    const std::string& corrected_seq,
    size_t correction_pos,
    int frame) const {

    CorrectionValidation result;
    result.valid = true;

    if (frame < 0 || original_seq.length() < 6) {
        return result;  // Can't validate without frame info
    }

    size_t len = original_seq.length();

    // Find codon containing the correction
    // Calculate offset from frame start, handling cases where correction_pos < frame
    size_t offset_from_frame = (correction_pos >= static_cast<size_t>(frame)) 
        ? (correction_pos - static_cast<size_t>(frame)) 
        : (correction_pos + 300 - static_cast<size_t>(frame));
    size_t pos_in_codon = offset_from_frame % 3;
    
    // Ensure codon_start doesn't underflow
    if (correction_pos < pos_in_codon) {
        return result;  // Can't find valid codon start
    }
    size_t codon_start = correction_pos - pos_in_codon;
    
    if (codon_start + 3 > len) {
        return result;  // Partial codon at end
    }

    // Check if correction creates a stop codon
    bool original_is_stop = is_stop_codon(original_seq.data() + codon_start);
    bool corrected_is_stop = is_stop_codon(corrected_seq.data() + codon_start);

    if (!original_is_stop && corrected_is_stop) {
        result.creates_stop = true;
        result.valid = false;
        return result;
    }

    // Check hexamer quality before/after
    if (len >= 6) {
        float original_llr = get_ensemble_hexamer_llr(original_seq, frame);
        float corrected_llr = get_ensemble_hexamer_llr(corrected_seq, frame);

        result.quality_delta = corrected_llr - original_llr;

        // Reject if correction significantly worsens hexamer score
        if (result.quality_delta < -0.5f) {
            result.worsens_hexamer = true;
            result.valid = false;
        }
    }

    // Check upstream/downstream context for newly created stops
    // (corrections that create stops in adjacent codons)
    if (codon_start >= 3) {
        bool upstream_stop = is_stop_codon(corrected_seq.data() + codon_start - 3);
        if (upstream_stop && !is_stop_codon(original_seq.data() + codon_start - 3)) {
            result.creates_stop = true;
            result.valid = false;
        }
    }

    if (codon_start + 6 <= len) {
        bool downstream_stop = is_stop_codon(corrected_seq.data() + codon_start + 3);
        if (downstream_stop && !is_stop_codon(original_seq.data() + codon_start + 3)) {
            result.creates_stop = true;
            result.valid = false;
        }
    }

    return result;
}

} // namespace agp
