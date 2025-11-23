/**
 * Frame scoring functions
 */

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include <cmath>
#include <algorithm>
#include <array>

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

// Dicodon (hexamer) scoring - optimized
float calculate_dicodon_score(const std::string& seq, int frame) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    float score_sum = 0.0f;
    int dicodon_count = 0;

    for (size_t i = frame; i + 5 < seq.length(); i += 3) {
        // Uppercase all 6 bases at once
        char c0 = fast_upper(seq[i]);
        char c1 = fast_upper(seq[i + 1]);
        char c2 = fast_upper(seq[i + 2]);
        char c3 = fast_upper(seq[i + 3]);
        char c4 = fast_upper(seq[i + 4]);
        char c5 = fast_upper(seq[i + 5]);

        // Skip if any N
        if (c0 == 'N' || c1 == 'N' || c2 == 'N' ||
            c3 == 'N' || c4 == 'N' || c5 == 'N') continue;

        // ATG start patterns (high)
        if (c0 == 'A' && c1 == 'T' && c2 == 'G') {
            score_sum += 0.9f;
        }
        // CTG (common Leu codon)
        else if (c0 == 'C' && c1 == 'T' && c2 == 'G') {
            score_sum += 0.8f;
        }
        // GAA patterns (common)
        else if (c0 == 'G' && c1 == 'A' && c2 == 'A') {
            score_sum += 0.85f;
        }
        // Stop-stop (penalize)
        else if ((c0 == 'T' && c1 == 'A' && c2 == 'A') ||
                 (c0 == 'T' && c1 == 'G' && c2 == 'A') ||
                 (c0 == 'T' && c1 == 'A' && c2 == 'G')) {
            if ((c3 == 'T') || (c4 == 'A' && c5 == 'A')) {
                score_sum += 0.05f;
            } else {
                score_sum += 0.3f;
            }
        }
        else {
            score_sum += 0.3f;
        }
        dicodon_count++;
    }

    if (dicodon_count == 0) return 0.5f;
    return std::max(0.0f, std::min(1.0f, score_sum / dicodon_count));
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

} // namespace agp
