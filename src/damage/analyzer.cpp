// DamageAnalyzer and CodonDamageAnalyzer implementations

#include "dart/damage_model.hpp"
#include "dart/codon_tables.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace dart {

// DamageAnalyzer implementation
DamageModel DamageAnalyzer::analyze(const std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return DamageModel::from_parameters(0.1f, 0.1f, 0.4f, 0.01f);
    }

    // Calculate C->T profile from 5' end
    auto ct_profile = calculate_ct_profile(sequences, 25);

    // Calculate G->A profile from 3' end
    auto ga_profile = calculate_ga_profile(sequences, 25);

    // Fit exponential decay model
    auto [lambda_5, delta_max_5, delta_bg_5] = fit_decay_model(ct_profile);
    auto [lambda_3, delta_max_3, delta_bg_3] = fit_decay_model(ga_profile);

    // Average parameters
    float lambda = (lambda_5 + lambda_3) / 2.0f;
    float delta_max = (delta_max_5 + delta_max_3) / 2.0f;
    float delta_bg = (delta_bg_5 + delta_bg_3) / 2.0f;

    return DamageModel::from_parameters(lambda, lambda, delta_max, delta_bg);
}

std::vector<float> DamageAnalyzer::calculate_ct_profile(
    const std::vector<std::string>& sequences,
    size_t max_positions) {

    std::vector<size_t> t_count(max_positions, 0);
    std::vector<size_t> total_count(max_positions, 0);

    for (const auto& seq : sequences) {
        for (size_t i = 0; i < std::min(max_positions, seq.length()); ++i) {
            char c = seq[i];
            if (c == 'T' || c == 't' || c == 'C' || c == 'c') {
                total_count[i]++;
                if (c == 'T' || c == 't') {
                    t_count[i]++;
                }
            }
        }
    }

    std::vector<float> profile(max_positions);
    for (size_t i = 0; i < max_positions; ++i) {
        if (total_count[i] > 0) {
            profile[i] = static_cast<float>(t_count[i]) / total_count[i];
        } else {
            profile[i] = 0.25f;  // Background
        }
    }

    return profile;
}

std::vector<float> DamageAnalyzer::calculate_ga_profile(
    const std::vector<std::string>& sequences,
    size_t max_positions) {

    std::vector<size_t> a_count(max_positions, 0);
    std::vector<size_t> total_count(max_positions, 0);

    for (const auto& seq : sequences) {
        size_t len = seq.length();
        for (size_t i = 0; i < std::min(max_positions, len); ++i) {
            size_t pos = len - i - 1;
            char c = seq[pos];
            if (c == 'A' || c == 'a' || c == 'G' || c == 'g') {
                total_count[i]++;
                if (c == 'A' || c == 'a') {
                    a_count[i]++;
                }
            }
        }
    }

    std::vector<float> profile(max_positions);
    for (size_t i = 0; i < max_positions; ++i) {
        if (total_count[i] > 0) {
            profile[i] = static_cast<float>(a_count[i]) / total_count[i];
        } else {
            profile[i] = 0.25f;
        }
    }

    return profile;
}

std::tuple<float, float, float> DamageAnalyzer::fit_decay_model(
    const std::vector<float>& damage_profile) {

    if (damage_profile.empty()) {
        return {0.1f, 0.4f, 0.01f};
    }

    // Simple fitting using first and last values
    float delta_max = damage_profile[0] - 0.25f;  // Subtract background base frequency
    delta_max = std::max(0.0f, delta_max);

    float delta_bg = damage_profile.back();

    // Estimate lambda from half-life of exponential component
    // δ(p) = δ_max * e^(-λp) + δ_bg, half-life when e^(-λp) = 0.5
    float half_damage = delta_max / 2.0f + delta_bg;
    size_t half_pos = 0;
    for (size_t i = 0; i < damage_profile.size(); ++i) {
        if (damage_profile[i] <= half_damage) {
            half_pos = i;
            break;
        }
    }

    float lambda = 0.693f / std::max(1.0f, static_cast<float>(half_pos));

    return {lambda, delta_max, delta_bg};
}

DamageModel DamageAnalyzer::infer_from_codon_usage(
    const std::vector<std::string>& coding_sequences) {

    // Count stop codons at different positions
    std::vector<size_t> stop_count(25, 0);
    std::vector<size_t> total_codons(25, 0);

    for (const auto& seq : coding_sequences) {
        for (size_t i = 0; i + 2 < seq.length(); i += 3) {
            size_t codon_pos = i / 3;
            if (codon_pos >= 25) break;

            total_codons[codon_pos]++;

            if (CodonTable::is_stop_codon(seq[i], seq[i+1], seq[i+2])) {
                stop_count[codon_pos]++;
            }
        }
    }

    // Stop codons should be rare in coding regions
    std::vector<float> stop_freq(25);
    for (size_t i = 0; i < 25; ++i) {
        if (total_codons[i] > 0) {
            stop_freq[i] = static_cast<float>(stop_count[i]) / total_codons[i];
        }
    }

    auto [lambda, delta_max, delta_bg] = fit_decay_model(stop_freq);

    return DamageModel::from_parameters(lambda, lambda, delta_max, delta_bg);
}

// CodonDamageAnalyzer implementation
Score CodonDamageAnalyzer::stop_codon_probability(const std::string& codon,
                                                 float ct_rate, float ga_rate) {
    if (codon.length() != 3) return 0.0f;

    float prob_creates_stop = 0.0f;

    for (int pos = 0; pos < 3; ++pos) {
        char original = codon[pos];
        std::string damaged = codon;

        // C->T damage
        if (original == 'C') {
            damaged[pos] = 'T';
            if (CodonTable::is_stop_codon(damaged[0], damaged[1], damaged[2])) {
                prob_creates_stop += ct_rate;
            }
        }

        // G->A damage
        damaged = codon;
        if (original == 'G') {
            damaged[pos] = 'A';
            if (CodonTable::is_stop_codon(damaged[0], damaged[1], damaged[2])) {
                prob_creates_stop += ga_rate;
            }
        }
    }

    return std::min(1.0f, prob_creates_stop);
}

std::vector<std::pair<std::string, float>> CodonDamageAnalyzer::get_damaged_variants(
    const std::string& codon,
    float ct_rate,
    float ga_rate) {

    std::vector<std::pair<std::string, float>> variants;

    if (codon.length() != 3) return variants;

    // Original codon (no damage)
    float no_damage_prob = 1.0f;
    for (char c : codon) {
        if (c == 'C') no_damage_prob *= (1.0f - ct_rate);
        if (c == 'G') no_damage_prob *= (1.0f - ga_rate);
    }
    variants.push_back({codon, no_damage_prob});

    // Single position damage
    for (int pos = 0; pos < 3; ++pos) {
        char original = codon[pos];

        // C->T
        if (original == 'C') {
            std::string damaged = codon;
            damaged[pos] = 'T';
            float prob = ct_rate;
            for (int i = 0; i < 3; ++i) {
                if (i != pos) {
                    if (codon[i] == 'C') prob *= (1.0f - ct_rate);
                    if (codon[i] == 'G') prob *= (1.0f - ga_rate);
                }
            }
            variants.push_back({damaged, prob});
        }

        // G->A
        if (original == 'G') {
            std::string damaged = codon;
            damaged[pos] = 'A';
            float prob = ga_rate;
            for (int i = 0; i < 3; ++i) {
                if (i != pos) {
                    if (codon[i] == 'C') prob *= (1.0f - ct_rate);
                    if (codon[i] == 'G') prob *= (1.0f - ga_rate);
                }
            }
            variants.push_back({damaged, prob});
        }
    }

    return variants;
}

std::vector<std::pair<std::string, float>> CodonDamageAnalyzer::most_vulnerable_codons() {
    std::vector<std::pair<std::string, float>> vulnerable;

    const std::vector<std::string> sensitive_codons = {
        "CGA",  // ->TGA (stop)
        "CAA",  // ->TAA (stop)
        "CAG",  // ->TAG (stop)
        "TGG",  // ->TGA (stop) via G->A
    };

    float typical_damage = 0.3f;

    for (const auto& codon : sensitive_codons) {
        float prob = stop_codon_probability(codon, typical_damage, typical_damage);
        vulnerable.push_back({codon, prob});
    }

    std::sort(vulnerable.begin(), vulnerable.end(),
             [](const auto& a, const auto& b) { return a.second > b.second; });

    return vulnerable;
}

}  // namespace dart
