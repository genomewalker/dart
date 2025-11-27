/**
 * Frame scoring functions
 */

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include "agp/gtdb_hexamer_table.hpp"
#include "agp/hexamer_damage_lookup.hpp"
#include "agp/multi_domain_hexamer.hpp"
#include <cmath>
#include <algorithm>
#include <array>
#include <cctype>

namespace agp {

// Translate sequence to protein (no intermediate string allocations)
std::string translate_sequence(const std::string& seq, int frame) {
    std::string protein;
    size_t num_codons = (seq.length() - frame) / 3;
    protein.reserve(num_codons);

    for (size_t i = frame; i + 2 < seq.length(); i += 3) {
        char aa = translate_codon_fast(seq[i], seq[i + 1], seq[i + 2]);
        protein += aa;
    }
    return protein;
}

// Calculate codon usage score directly (no vector allocation)
float calculate_codon_usage_score_direct(const std::string& seq, int frame) {
    float log_prob_sum = 0.0f;
    int valid_codons = 0;

    for (size_t i = frame; i + 2 < seq.length(); i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        // Skip codons with N
        if (c1 == 'N' || c2 == 'N' || c3 == 'N') continue;

        float freq = get_codon_freq_fast(c1, c2, c3);
        log_prob_sum += std::log(freq + 0.01f);
        valid_codons++;
    }

    if (valid_codons == 0) return 0.0f;

    float avg_log_prob = log_prob_sum / valid_codons;
    float score = (avg_log_prob + 4.6f) / 4.6f;
    return std::max(0.0f, std::min(1.0f, score));
}

// Dicodon (hexamer) scoring using domain-specific frequencies
// Uses log-likelihood ratio against uniform background
// Now uses the active domain set via set_active_domain()
// In ensemble mode, uses weighted scoring across all domains
float calculate_dicodon_score(const std::string& seq, int frame) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    // Check for ensemble mode (metagenome weighted scoring)
    if (get_ensemble_mode()) {
        // Use probability-weighted scoring across all domains
        return calculate_dicodon_score_weighted(seq, frame, get_domain_probs());
    }

    float log_prob_sum = 0.0f;
    int hexamer_count = 0;

    // Background frequency (uniform = 1/4096)
    static constexpr float BACKGROUND_FREQ = 1.0f / 4096.0f;

    // Get the active domain for scoring
    Domain active_domain = get_active_domain();

    // Buffer for uppercase hexamer
    char hexbuf[7];
    hexbuf[6] = '\0';

    for (size_t i = frame; i + 5 < seq.length(); i += 3) {
        // Uppercase all 6 bases
        hexbuf[0] = fast_upper(seq[i]);
        hexbuf[1] = fast_upper(seq[i + 1]);
        hexbuf[2] = fast_upper(seq[i + 2]);
        hexbuf[3] = fast_upper(seq[i + 3]);
        hexbuf[4] = fast_upper(seq[i + 4]);
        hexbuf[5] = fast_upper(seq[i + 5]);

        // Skip if any N
        if (hexbuf[0] == 'N' || hexbuf[1] == 'N' || hexbuf[2] == 'N' ||
            hexbuf[3] == 'N' || hexbuf[4] == 'N' || hexbuf[5] == 'N') continue;

        // Get domain-specific frequency (uses active domain)
        float freq = get_domain_hexamer_freq(active_domain, hexbuf);

        // If hexamer has 0 frequency (rare with N's filtered),
        // use a small pseudocount
        if (freq < 1e-9f) {
            freq = BACKGROUND_FREQ * 0.01f;  // Very rare
        }

        // Log-likelihood ratio: log(P_coding / P_background)
        log_prob_sum += std::log(freq / BACKGROUND_FREQ);
        hexamer_count++;
    }

    if (hexamer_count == 0) return 0.5f;

    // Normalize: average log-likelihood ratio
    // Transform to 0-1 range using sigmoid-like scaling
    // Typical range for real coding: 1-3 bits per hexamer
    float avg_llr = log_prob_sum / hexamer_count;

    // Map typical range [-2, 4] to [0.1, 0.9]
    // avg_llr of 0 = neutral (0.5 score)
    // avg_llr of +2 = strong coding signal (0.7-0.8)
    // avg_llr of -2 = strong non-coding signal (0.2-0.3)
    float score = 0.5f + 0.15f * avg_llr;

    return std::max(0.0f, std::min(1.0f, score));
}

// Dipeptide scoring - optimized with array lookup
float calculate_dipeptide_score(const std::string& protein) {
    if (protein.length() < 2) return 0.5f;

    // Pre-computed dipeptide scores (indexed by (aa1-'A')*26 + (aa2-'A'))
    static const std::array<float, 676> DIPEPTIDE_SCORES = []() {
        std::array<float, 676> arr{};
        arr.fill(0.5f);  // Default

        auto set = [&arr](char a1, char a2, float score) {
            if (a1 >= 'A' && a1 <= 'Z' && a2 >= 'A' && a2 <= 'Z') {
                arr[(a1 - 'A') * 26 + (a2 - 'A')] = score;
            }
        };

        // Most common
        set('L', 'L', 1.00f); set('A', 'A', 0.96f); set('A', 'L', 0.96f); set('L', 'A', 0.95f);
        set('L', 'S', 0.89f); set('V', 'L', 0.88f); set('L', 'G', 0.88f); set('E', 'L', 0.87f);
        set('G', 'L', 0.87f); set('S', 'L', 0.87f); set('L', 'E', 0.87f); set('L', 'V', 0.86f);
        set('A', 'G', 0.86f); set('V', 'A', 0.85f); set('A', 'V', 0.85f); set('S', 'S', 0.85f);
        // Common
        set('G', 'G', 0.84f); set('L', 'K', 0.84f); set('E', 'A', 0.84f); set('E', 'E', 0.84f);
        set('G', 'A', 0.83f); set('L', 'R', 0.83f); set('A', 'E', 0.83f); set('R', 'L', 0.83f);
        set('T', 'L', 0.83f); set('D', 'L', 0.83f); set('L', 'D', 0.83f); set('K', 'L', 0.82f);
        set('I', 'L', 0.82f); set('L', 'T', 0.82f); set('V', 'V', 0.82f); set('L', 'I', 0.82f);
        set('A', 'S', 0.82f); set('S', 'G', 0.81f); set('G', 'V', 0.81f); set('L', 'P', 0.81f);
        set('S', 'A', 0.81f); set('I', 'A', 0.81f); set('E', 'K', 0.80f); set('G', 'S', 0.80f);
        set('V', 'E', 0.80f); set('A', 'I', 0.80f); set('A', 'R', 0.80f); set('E', 'V', 0.79f);
        set('K', 'A', 0.79f); set('V', 'G', 0.79f); set('V', 'S', 0.79f); set('A', 'K', 0.79f);
        set('K', 'E', 0.79f); set('D', 'A', 0.79f);
        // Rare (penalize)
        set('W', 'W', 0.12f); set('W', 'C', 0.12f); set('C', 'W', 0.12f); set('M', 'W', 0.12f);
        set('W', 'M', 0.12f); set('C', 'M', 0.12f); set('M', 'C', 0.12f); set('H', 'W', 0.12f);
        set('W', 'H', 0.12f); set('C', 'C', 0.13f); set('W', 'Y', 0.13f); set('Y', 'W', 0.13f);
        set('C', 'H', 0.13f); set('H', 'C', 0.13f); set('W', 'P', 0.14f); set('C', 'Y', 0.14f);

        return arr;
    }();

    float score_sum = 0.0f;
    int pair_count = 0;

    for (size_t i = 0; i + 1 < protein.length(); ++i) {
        char aa1 = protein[i];
        char aa2 = protein[i + 1];

        if (aa1 == '*' || aa2 == '*') continue;
        if (aa1 < 'A' || aa1 > 'Z' || aa2 < 'A' || aa2 > 'Z') continue;

        score_sum += DIPEPTIDE_SCORES[(aa1 - 'A') * 26 + (aa2 - 'A')];
        pair_count++;
    }

    if (pair_count == 0) return 0.5f;
    return std::max(0.0f, std::min(1.0f, score_sum / pair_count));
}

// Legacy functions that allocate vectors (for compatibility)
std::vector<std::string> extract_codons(const std::string& seq, int frame) {
    std::vector<std::string> codons;
    size_t num_codons = (seq.length() - frame) / 3;
    codons.reserve(num_codons);

    for (size_t i = frame; i + 2 < seq.length(); i += 3) {
        std::string codon;
        codon.reserve(3);
        codon += fast_upper(seq[i]);
        codon += fast_upper(seq[i + 1]);
        codon += fast_upper(seq[i + 2]);
        codons.push_back(std::move(codon));
    }
    return codons;
}

float FrameSelector::calculate_codon_usage_score(const std::vector<std::string>& codons) {
    if (codons.empty()) return 0.0f;

    float log_prob_sum = 0.0f;
    int valid_codons = 0;

    for (const auto& codon : codons) {
        if (codon.find('N') != std::string::npos) continue;
        if (codon.length() != 3) continue;

        float freq = get_codon_freq_fast(codon[0], codon[1], codon[2]);
        log_prob_sum += std::log(freq + 0.01f);
        valid_codons++;
    }

    if (valid_codons == 0) return 0.0f;

    float avg_log_prob = log_prob_sum / valid_codons;
    float score = (avg_log_prob + 4.6f) / 4.6f;
    return std::max(0.0f, std::min(1.0f, score));
}

float FrameSelector::calculate_stop_penalty(const std::string& protein) {
    if (protein.empty()) return 0.0f;

    int internal_stops = 0;
    for (size_t i = 0; i + 1 < protein.length(); ++i) {
        if (protein[i] == '*') {
            internal_stops++;
        }
    }

    if (internal_stops == 0) return 1.0f;
    if (internal_stops == 1) return 0.2f;
    return 0.05f;
}

float FrameSelector::calculate_aa_composition_score(const std::string& protein) {
    if (protein.empty()) return 0.0f;

    // Use array instead of map for AA counts
    std::array<int, 26> aa_counts{};
    int total = 0;

    for (char aa : protein) {
        if (aa != '*' && aa >= 'A' && aa <= 'Z') {
            aa_counts[aa - 'A']++;
            total++;
        }
    }

    if (total == 0) return 0.0f;

    float divergence = 0.0f;
    const char amino_acids[] = "ACDEFGHIKLMNPQRSTVWY";

    for (const char* p = amino_acids; *p != '\0'; ++p) {
        char aa = *p;
        float expected = get_aa_freq_fast(aa);
        float observed = static_cast<float>(aa_counts[aa - 'A']) / total;
        float diff = observed - expected;
        divergence += (diff * diff) / (expected + 0.01f);
    }

    float score = 1.0f - std::min(1.0f, divergence / 3.0f);

    float total_f = static_cast<float>(total);
    if (aa_counts['C' - 'A'] / total_f > 0.15f) score -= 0.1f;
    if (aa_counts['W' - 'A'] / total_f > 0.15f) score -= 0.1f;
    if (aa_counts['M' - 'A'] / total_f > 0.15f) score -= 0.1f;

    int common_present = 0;
    const char common_aas[] = {'L', 'A', 'E', 'G', 'V', 'S', 'K', 'R'};
    for (char aa : common_aas) {
        if (aa_counts[aa - 'A'] / total_f > 0.03f) common_present++;
    }
    if (common_present >= 4) score += 0.1f;

    return std::max(0.0f, std::min(1.0f, score));
}

float FrameSelector::calculate_damage_consistency(
    const std::string& original_seq,
    const std::string& frame_seq,
    int frame,
    bool forward,
    const DamageProfile& damage) {

    if (frame_seq.length() < 10) return 0.5f;

    int t_count_5prime = 0, c_count_5prime = 0;
    int a_count_3prime = 0, g_count_3prime = 0;

    size_t check_len = std::min(size_t(5), frame_seq.length());

    for (size_t i = 0; i < check_len; ++i) {
        char base = fast_upper(frame_seq[i]);
        if (base == 'T') t_count_5prime++;
        else if (base == 'C') c_count_5prime++;

        base = fast_upper(frame_seq[frame_seq.length() - 1 - i]);
        if (base == 'A') a_count_3prime++;
        else if (base == 'G') g_count_3prime++;
    }

    float t_ratio_5 = (t_count_5prime + c_count_5prime > 0) ?
                      static_cast<float>(t_count_5prime) / (t_count_5prime + c_count_5prime) : 0.5f;
    float a_ratio_3 = (a_count_3prime + g_count_3prime > 0) ?
                      static_cast<float>(a_count_3prime) / (a_count_3prime + g_count_3prime) : 0.5f;

    float score = 0.5f;
    if (forward) {
        if (t_ratio_5 > 0.6f) score += 0.15f;
        if (a_ratio_3 > 0.6f) score += 0.15f;
        if (t_ratio_5 < 0.3f && a_ratio_3 < 0.3f) score -= 0.1f;
    } else {
        if (t_ratio_5 < 0.4f) score += 0.1f;
        if (a_ratio_3 < 0.4f) score += 0.1f;
    }

    return std::max(0.0f, std::min(1.0f, score));
}

// Hexamer-based stop codon correction
// Corrects stop codons that are likely damage artifacts based on GTDB hexamer evidence
// Returns: (corrected sequence, number of corrections made)
std::pair<std::string, size_t> correct_stop_codons_hexamer(
    const std::string& seq,
    int frame,
    size_t terminal_length) {

    std::string corrected = seq;
    size_t corrections = 0;

    if (seq.length() < static_cast<size_t>(frame + 6)) {
        return {corrected, 0};
    }

    // Scan for stop codons in-frame
    for (size_t i = frame; i + 2 < seq.length(); i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        // Check if this is a stop codon (TAA, TAG, TGA)
        bool is_stop = false;
        char original_base = '\0';
        size_t correction_pos = 0;

        if (c1 == 'T' && c2 == 'A' && c3 == 'A') {
            // TAA - could come from CAA (C->T damage at position 0)
            is_stop = true;
            original_base = 'C';
            correction_pos = i;
        } else if (c1 == 'T' && c2 == 'A' && c3 == 'G') {
            // TAG - could come from CAG (C->T damage at position 0)
            is_stop = true;
            original_base = 'C';
            correction_pos = i;
        } else if (c1 == 'T' && c2 == 'G' && c3 == 'A') {
            // TGA - could come from CGA (C->T damage at position 0)
            is_stop = true;
            original_base = 'C';
            correction_pos = i;
        }

        if (!is_stop) continue;

        // Check position - only correct in terminal regions where damage is likely
        bool in_5prime_terminal = (i < terminal_length);
        bool in_3prime_terminal = (i >= seq.length() - terminal_length - 3);

        if (!in_5prime_terminal && !in_3prime_terminal) continue;

        // Get hexamer context if available (codon + next 3 bases)
        float damage_prob = 0.0f;
        if (i + 5 < seq.length()) {
            char hexbuf[7];
            hexbuf[0] = fast_upper(seq[i]);
            hexbuf[1] = fast_upper(seq[i + 1]);
            hexbuf[2] = fast_upper(seq[i + 2]);
            hexbuf[3] = fast_upper(seq[i + 3]);
            hexbuf[4] = fast_upper(seq[i + 4]);
            hexbuf[5] = fast_upper(seq[i + 5]);
            hexbuf[6] = '\0';

            // Get damage probability from hexamer lookup
            damage_prob = get_hexamer_damage_prob(hexbuf);

            // Also check GTDB frequency - stop-starting hexamers should be very rare
            float gtdb_freq = get_hexamer_freq(hexbuf);

            // If GTDB frequency is very low (< 1e-6), this hexamer is rare in coding
            // which supports it being a damage artifact
            if (gtdb_freq < 1e-6f && damage_prob < 0.5f) {
                // Boost damage probability for very rare hexamers
                damage_prob = std::max(damage_prob, 0.7f);
            }
        } else {
            // Short sequence, use position-based estimate
            // 5' terminal damage rate peaks at ~30-40% for ancient DNA
            if (in_5prime_terminal) {
                float pos_factor = 1.0f - static_cast<float>(i) / terminal_length;
                damage_prob = 0.3f * pos_factor;  // Conservative estimate
            }
        }

        // Correct if damage probability is high enough
        if (damage_prob >= 0.5f) {
            corrected[correction_pos] = original_base;
            corrections++;
        }
    }

    return {corrected, corrections};
}

// Calculate hexamer log-likelihood score for frame selection
// Higher score = more likely to be real coding sequence
// Uses GTDB hexamer frequencies vs uniform background
float calculate_hexamer_llr_score(const std::string& seq, int frame) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    float log_prob_sum = 0.0f;
    int hexamer_count = 0;

    // Background frequency (uniform = 1/4096)
    static constexpr float BACKGROUND_FREQ = 1.0f / 4096.0f;
    static constexpr float LOG_BACKGROUND = -8.317766f;  // log(1/4096)

    char hexbuf[7];
    hexbuf[6] = '\0';

    for (size_t i = frame; i + 5 < seq.length(); i += 3) {
        hexbuf[0] = fast_upper(seq[i]);
        hexbuf[1] = fast_upper(seq[i + 1]);
        hexbuf[2] = fast_upper(seq[i + 2]);
        hexbuf[3] = fast_upper(seq[i + 3]);
        hexbuf[4] = fast_upper(seq[i + 4]);
        hexbuf[5] = fast_upper(seq[i + 5]);

        // Skip if any N
        if (hexbuf[0] == 'N' || hexbuf[1] == 'N' || hexbuf[2] == 'N' ||
            hexbuf[3] == 'N' || hexbuf[4] == 'N' || hexbuf[5] == 'N') continue;

        float freq = get_hexamer_freq(hexbuf);

        // Apply damage correction: if hexamer has high damage probability,
        // use the frequency of the corrected (C instead of T) version
        float damage_prob = get_hexamer_damage_prob(hexbuf);
        if (damage_prob > 0.5f && hexbuf[0] == 'T') {
            // Try corrected version (T->C at position 0)
            char corrected_hex[7];
            corrected_hex[0] = 'C';
            for (int j = 1; j < 6; j++) corrected_hex[j] = hexbuf[j];
            corrected_hex[6] = '\0';

            float corrected_freq = get_hexamer_freq(corrected_hex);
            // Use weighted average based on damage probability
            freq = (1.0f - damage_prob) * freq + damage_prob * corrected_freq;
        }

        // Pseudocount for zero frequencies
        if (freq < 1e-9f) {
            freq = BACKGROUND_FREQ * 0.01f;
        }

        // Log-likelihood ratio
        log_prob_sum += std::log(freq) - LOG_BACKGROUND;
        hexamer_count++;
    }

    if (hexamer_count == 0) return 0.0f;

    // Return average log-likelihood ratio
    return log_prob_sum / hexamer_count;
}

// Damage-aware stop penalty using hexamer evidence
// Returns a multiplier (0-1) where higher = fewer likely-real stops
float calculate_stop_penalty_hexamer_aware(
    const std::string& seq,
    const std::string& protein,
    int frame,
    size_t terminal_length) {

    if (protein.empty()) return 0.0f;

    int real_internal_stops = 0;
    int damage_induced_stops = 0;

    // Scan for internal stop codons (not the last one)
    size_t codon_idx = 0;
    for (size_t i = frame; i + 2 < seq.length() && codon_idx + 1 < protein.length(); i += 3, codon_idx++) {
        if (protein[codon_idx] != '*') continue;

        // This is an internal stop - check if it's likely damage
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        // Only damage-correctable stops start with T (from C->T)
        bool could_be_damage = (c1 == 'T');

        // Check position
        bool in_5prime = (i < terminal_length);
        bool in_3prime = (i >= seq.length() - terminal_length - 3);
        bool in_terminal = in_5prime || in_3prime;

        // Get hexamer damage probability if available
        float damage_prob = 0.0f;
        if (could_be_damage && i + 5 < seq.length()) {
            char hexbuf[7];
            hexbuf[0] = c1; hexbuf[1] = c2; hexbuf[2] = c3;
            hexbuf[3] = fast_upper(seq[i + 3]);
            hexbuf[4] = fast_upper(seq[i + 4]);
            hexbuf[5] = fast_upper(seq[i + 5]);
            hexbuf[6] = '\0';

            damage_prob = get_hexamer_damage_prob(hexbuf);

            // Boost for terminal regions
            if (in_terminal && damage_prob < 0.3f) {
                float pos_factor = in_5prime ? (1.0f - static_cast<float>(i) / terminal_length)
                                            : (static_cast<float>(i - (seq.length() - terminal_length)) / terminal_length);
                damage_prob = std::max(damage_prob, 0.2f * pos_factor);
            }
        }

        // Classify the stop
        if (damage_prob >= 0.5f) {
            damage_induced_stops++;
        } else {
            real_internal_stops++;
        }
    }

    // Calculate penalty
    // Real stops are heavily penalized, damage stops are lightly penalized
    if (real_internal_stops == 0 && damage_induced_stops == 0) return 1.0f;
    if (real_internal_stops == 0) {
        // Only damage-induced stops - light penalty
        if (damage_induced_stops == 1) return 0.8f;
        return 0.6f;
    }
    if (real_internal_stops == 1) return 0.2f;
    return 0.05f;
}

} // namespace agp
