// Damage probability calculations (ancient likelihood, AA distributions)

#include "agp/damage_model.hpp"
#include "agp/codon_tables.hpp"
#include <cmath>
#include <algorithm>
#include <cctype>

namespace agp {

Score DamageModel::calculate_ancient_likelihood(const Sequence& seq) const {
    if (seq.empty()) return 0.0f;

    size_t seq_len = seq.length();

    // Get expected frequencies under the ancient model (H1)
    const auto& profile = create_profile_cached(seq_len);

    // We compare two hypotheses:
    // H1 (Ancient): Observed damage patterns (C->T at 5', G->A at 3') follow decay model
    // H0 (Modern): Observed patterns follow flat background rate

    float log_likelihood_ratio = 0.0f;

    // Check up to 15bp at each end, but avoid overlap
    size_t check_len = std::min(size_t(15), seq_len / 2);
    if (check_len == 0 && seq_len > 0) check_len = 1;

    // Clamp background probability to avoid log(0) or division by zero
    float p_bg = std::clamp(delta_background_, 0.001f, 0.999f);

    // 5' end: Check for C->T damage signal
    for (size_t i = 0; i < check_len; ++i) {
        char base = seq[i];

        // Only look at Pyrimidines (C/T)
        if (base == 'C' || base == 'c') {
            // Observed C (No damage)
            float p_ancient = std::clamp(profile.ct_prob_5prime[i], 0.001f, 0.999f);
            log_likelihood_ratio += std::log((1.0f - p_ancient) / (1.0f - p_bg));
        } else if (base == 'T' || base == 't') {
            // Observed T (Potential damage)
            float p_ancient = std::clamp(profile.ct_prob_5prime[i], 0.001f, 0.999f);
            log_likelihood_ratio += std::log((1.0f + p_ancient) / (1.0f + p_bg));
        }
    }

    // 3' end: Check for G->A damage signal
    for (size_t i = 0; i < check_len; ++i) {
        size_t pos = seq_len - 1 - i;
        char base = seq[pos];

        // Only look at Purines (G/A)
        if (base == 'G' || base == 'g') {
            // Observed G (No damage)
            float p_ancient = std::clamp(profile.ga_prob_3prime[pos], 0.001f, 0.999f);
            log_likelihood_ratio += std::log((1.0f - p_ancient) / (1.0f - p_bg));
        } else if (base == 'A' || base == 'a') {
            // Observed A (Potential damage)
            float p_ancient = std::clamp(profile.ga_prob_3prime[pos], 0.001f, 0.999f);
            log_likelihood_ratio += std::log((1.0f + p_ancient) / (1.0f + p_bg));
        }
    }

    // Convert LLR to probability (sigmoid)
    return 1.0f / (1.0f + std::exp(-log_likelihood_ratio));
}

// Calculate amino acid probability distribution for a codon considering damage
// String overload - delegates to 3-char version
std::array<float, 21> DamageModel::calculate_aa_probabilities(
    const std::string& codon,
    size_t codon_start,
    size_t seq_len) const {
    if (codon.length() != 3) return std::array<float, 21>{};
    return calculate_aa_probabilities(codon[0], codon[1], codon[2], codon_start, seq_len);
}

// Allocation-free overload - takes 3 chars directly
std::array<float, 21> DamageModel::calculate_aa_probabilities(
    char c0, char c1, char c2,
    size_t codon_start,
    size_t seq_len) const {

    // Initialize all probabilities to zero
    // Indices: A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9,
    //          M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19, *=20
    std::array<float, 21> aa_probs{};
    aa_probs.fill(0.0f);

    // Get damage profile for this sequence length (cached to avoid O(N²))
    const auto& profile = create_profile_cached(seq_len);

    // Positions of each nucleotide
    size_t pos0 = codon_start;
    size_t pos1 = codon_start + 1;
    size_t pos2 = codon_start + 2;

    // Get position-dependent damage rates
    // C→T at 5' end (forward strand damage)
    float ct_rate_0 = (pos0 < profile.ct_prob_5prime.size()) ? profile.ct_prob_5prime[pos0] : delta_background_;
    float ct_rate_1 = (pos1 < profile.ct_prob_5prime.size()) ? profile.ct_prob_5prime[pos1] : delta_background_;
    float ct_rate_2 = (pos2 < profile.ct_prob_5prime.size()) ? profile.ct_prob_5prime[pos2] : delta_background_;

    // G→A at 3' end (forward strand damage)
    float ga_rate_0 = (pos0 < profile.ga_prob_3prime.size()) ? profile.ga_prob_3prime[pos0] : delta_background_;
    float ga_rate_1 = (pos1 < profile.ga_prob_3prime.size()) ? profile.ga_prob_3prime[pos1] : delta_background_;
    float ga_rate_2 = (pos2 < profile.ga_prob_3prime.size()) ? profile.ga_prob_3prime[pos2] : delta_background_;

    // Helper to convert amino acid char to index
    auto aa_to_idx = [](char aa) -> int {
        switch (aa) {
            case 'A': return 0;  case 'C': return 1;  case 'D': return 2;
            case 'E': return 3;  case 'F': return 4;  case 'G': return 5;
            case 'H': return 6;  case 'I': return 7;  case 'K': return 8;
            case 'L': return 9;  case 'M': return 10; case 'N': return 11;
            case 'P': return 12; case 'Q': return 13; case 'R': return 14;
            case 'S': return 15; case 'T': return 16; case 'V': return 17;
            case 'W': return 18; case 'Y': return 19; case '*': return 20;
            default: return -1;
        }
    };

    // For each position, determine possible original nucleotides
    auto get_possible_originals = [](char observed, float ct_rate, float ga_rate)
        -> std::vector<std::pair<char, float>> {
        std::vector<std::pair<char, float>> originals;
        char obs = std::toupper(observed);

        if (obs == 'T') {
            float denom = 1.0f + ct_rate;
            originals.push_back({'T', 1.0f / denom});
            originals.push_back({'C', ct_rate / denom});
        }
        else if (obs == 'A') {
            float denom = 1.0f + ga_rate;
            originals.push_back({'A', 1.0f / denom});
            originals.push_back({'G', ga_rate / denom});
        }
        else if (obs == 'C') {
            originals.push_back({'C', 1.0f});
        }
        else if (obs == 'G') {
            originals.push_back({'G', 1.0f});
        }
        else {
            // N or unknown - uniform distribution
            originals.push_back({'A', 0.25f});
            originals.push_back({'C', 0.25f});
            originals.push_back({'G', 0.25f});
            originals.push_back({'T', 0.25f});
        }
        return originals;
    };

    // Get possible original nucleotides for each position
    auto orig0 = get_possible_originals(c0, ct_rate_0, ga_rate_0);
    auto orig1 = get_possible_originals(c1, ct_rate_1, ga_rate_1);
    auto orig2 = get_possible_originals(c2, ct_rate_2, ga_rate_2);

    // Enumerate all possible original codons and their probabilities
    float total_prob = 0.0f;
    for (const auto& [nt0, p0] : orig0) {
        for (const auto& [nt1, p1] : orig1) {
            for (const auto& [nt2, p2] : orig2) {
                float codon_prob = p0 * p1 * p2;
                char aa = CodonTable::translate_codon(nt0, nt1, nt2);
                int idx = aa_to_idx(aa);
                if (idx >= 0 && idx < 21) {
                    aa_probs[idx] += codon_prob;
                    total_prob += codon_prob;
                }
            }
        }
    }

    // Normalize probabilities
    if (total_prob > 1e-9f) {
        for (auto& p : aa_probs) {
            p /= total_prob;
        }
    } else {
        float uniform_prob = 1.0f / 21.0f;
        for (auto& p : aa_probs) {
            p = uniform_prob;
        }
    }

    return aa_probs;
}

// Translate DNA with damage-aware probability distribution
std::pair<std::string, std::vector<float>> DamageModel::translate_with_confidence(
    const std::string& seq,
    int frame) const {

    std::string protein;
    std::vector<float> confidences;

    if (seq.empty()) return {protein, confidences};

    size_t seq_len = seq.length();
    size_t start = (frame >= 0) ? static_cast<size_t>(frame) : 0;

    // Index to amino acid character
    auto idx_to_aa = [](int idx) -> char {
        const char* aas = "ACDEFGHIKLMNPQRSTVWY*";
        return (idx >= 0 && idx < 21) ? aas[idx] : 'X';
    };

    for (size_t i = start; i + 2 < seq_len; i += 3) {
        // Calculate AA probability distribution (allocation-free - no substr)
        auto aa_probs = calculate_aa_probabilities(seq[i], seq[i+1], seq[i+2], i, seq_len);

        // Find the most likely amino acid
        int best_idx = 0;
        float best_prob = aa_probs[0];
        for (int j = 1; j < 21; ++j) {
            if (aa_probs[j] > best_prob) {
                best_prob = aa_probs[j];
                best_idx = j;
            }
        }

        protein += idx_to_aa(best_idx);
        confidences.push_back(best_prob);
    }

    return {protein, confidences};
}

}  // namespace agp
