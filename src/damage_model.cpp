#include "agp/damage_model.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace agp {

DamageModel DamageModel::from_parameters(float lambda_5prime,
                                        float lambda_3prime,
                                        float delta_max,
                                        float delta_background) {
    DamageModel model;
    model.lambda_5prime_ = lambda_5prime;
    model.lambda_3prime_ = lambda_3prime;
    model.delta_max_ = delta_max;
    model.delta_background_ = delta_background;
    return model;
}

DamageModel DamageModel::estimate_from_sequences(
    const std::vector<std::string>& sequences,
    const std::vector<std::vector<State>>& state_paths) {

    return DamageAnalyzer::analyze(sequences);
}

DamageProfile DamageModel::create_profile(size_t seq_len) const {
    DamageProfile profile;
    profile.lambda_5prime = lambda_5prime_;
    profile.lambda_3prime = lambda_3prime_;
    profile.delta_max = delta_max_;
    profile.delta_background = delta_background_;

    // Precompute damage probabilities for each position
    profile.ct_prob_5prime.resize(seq_len);
    profile.ct_prob_3prime.resize(seq_len);
    profile.ga_prob_5prime.resize(seq_len);
    profile.ga_prob_3prime.resize(seq_len);

    for (size_t pos = 0; pos < seq_len; ++pos) {
        // C->T from 5' end
        float dist_5prime = static_cast<float>(pos);
        profile.ct_prob_5prime[pos] = delta_max_ * std::exp(-lambda_5prime_ * dist_5prime)
                                    + delta_background_;

        // G->A from 3' end (reverse complement of C->T)
        float dist_3prime = static_cast<float>(seq_len - pos - 1);
        profile.ga_prob_3prime[pos] = delta_max_ * std::exp(-lambda_3prime_ * dist_3prime)
                                    + delta_background_;
    }

    // Use empirical data if available
    if (empirical_ct_5prime_) {
        for (size_t i = 0; i < std::min(empirical_ct_5prime_->size(), seq_len); ++i) {
            profile.ct_prob_5prime[i] = (*empirical_ct_5prime_)[i];
        }
    }

    if (empirical_ga_3prime_) {
        for (size_t i = 0; i < std::min(empirical_ga_3prime_->size(), seq_len); ++i) {
            profile.ga_prob_3prime[seq_len - i - 1] = (*empirical_ga_3prime_)[i];
        }
    }

    return profile;
}

const DamageProfile& DamageModel::create_profile_cached(size_t seq_len) {
    // Clamp to cache size
    if (seq_len >= profile_cache_.size()) {
        seq_len = profile_cache_.size() - 1;
    }

    // Check cache first (double-checked locking for performance)
    if (profile_cache_[seq_len].has_value()) {
        return *profile_cache_[seq_len];
    }

    // Lock and create profile if not cached
    std::lock_guard<std::mutex> lock(cache_mutex_);

    // Double-check after acquiring lock
    if (!profile_cache_[seq_len].has_value()) {
        profile_cache_[seq_len] = create_profile(seq_len);
    }

    return *profile_cache_[seq_len];
}

Score DamageModel::calculate_ancient_likelihood(const Sequence& seq) const {
    if (seq.empty()) return 0.0f;

    float log_likelihood = 0.0f;
    size_t seq_len = seq.length();

    // Check for C->T enrichment at 5' end
    size_t ct_count_5prime = 0;
    size_t c_count_5prime = 0;

    // Check first 10 positions
    size_t check_len = std::min(size_t(10), seq_len);

    for (size_t i = 0; i < check_len; ++i) {
        if (seq[i] == 'C' || seq[i] == 'c') {
            c_count_5prime++;
        } else if (seq[i] == 'T' || seq[i] == 't') {
            ct_count_5prime++;
        }
    }

    // Expected damage rate at 5' end
    float expected_ct = delta_max_ * 0.5f;  // Average over first few positions

    if (c_count_5prime + ct_count_5prime > 0) {
        float observed_ct = static_cast<float>(ct_count_5prime) /
                           (c_count_5prime + ct_count_5prime);

        // Log likelihood ratio
        log_likelihood += (observed_ct - delta_background_) / (expected_ct + 0.1f);
    }

    return 1.0f / (1.0f + std::exp(-log_likelihood));
}

bool DamageModel::is_ancient_compatible() const {
    return delta_max_ > 0.1f && delta_max_ < 0.6f &&
           lambda_5prime_ > 0.01f && lambda_5prime_ < 1.0f;
}

std::pair<std::string, size_t> DamageModel::correct_damage(
    const std::string& seq,
    float confidence_threshold,
    int frame) const {

    if (seq.empty()) return {seq, 0};

    std::string corrected = seq;
    size_t corrections = 0;
    size_t seq_len = seq.length();

    // Create damage profile for this sequence length
    auto profile = create_profile(seq_len);

    // Frame-aware thresholds:
    // - Wobble position (codon pos 3): lower threshold (more likely to correct)
    //   because synonymous mutations are tolerated, so observed T is likely damage
    // - Non-wobble (pos 1, 2): higher threshold (be more careful)
    //   because non-synonymous mutations are selected against
    auto get_threshold = [&](size_t pos, float base_threshold) -> float {
        if (frame < 0) return base_threshold;  // No frame info
        int codon_pos = (pos - frame + 300) % 3;
        if (codon_pos == 2) {
            // Wobble position: more aggressive correction
            return base_threshold * 0.7f;
        } else if (codon_pos == 0) {
            // First position: most conservative
            return base_threshold * 1.2f;
        }
        return base_threshold;
    };

    // CpG context detection
    auto is_cpg_context = [&](size_t pos, bool check_5prime) -> bool {
        if (check_5prime) {
            // Check if position is followed by G (CpG or TpG from damage)
            return (pos + 1 < seq_len && (seq[pos + 1] == 'G' || seq[pos + 1] == 'g'));
        } else {
            // Check if position is preceded by C (CpG context for G→A)
            return (pos > 0 && (seq[pos - 1] == 'C' || seq[pos - 1] == 'c'));
        }
    };

    // Only correct positions very close to the ends where damage is most likely
    // C->T at 5' end (first few positions)
    size_t max_5prime_pos = std::min(size_t(7), seq_len);  // Extended to 7
    for (size_t i = 0; i < max_5prime_pos; ++i) {
        float ct_prob = profile.ct_prob_5prime[i];
        float threshold = get_threshold(i, confidence_threshold);

        // CpG sites have higher damage rate, so lower threshold
        if (is_cpg_context(i, true)) {
            threshold *= 0.8f;
        }

        // Only correct if damage probability exceeds threshold
        if (ct_prob > threshold && ct_prob > 0.15f) {
            if (seq[i] == 'T' || seq[i] == 't') {
                // Determine codon boundaries based on frame
                size_t codon_start = (frame >= 0) ?
                    ((i / 3) * 3 + frame) % 3 + ((i - frame + 300) / 3) * 3 - 300 + frame :
                    (i / 3) * 3;
                codon_start = std::min(codon_start, i);

                if (codon_start + 2 < seq_len) {
                    std::string test_codon;
                    for (size_t j = 0; j < 3 && codon_start + j < seq_len; j++) {
                        test_codon += corrected[codon_start + j];
                    }
                    if (test_codon.length() == 3) {
                        size_t pos_in_codon = i - codon_start;
                        char original = test_codon[pos_in_codon];
                        test_codon[pos_in_codon] = (original == 'T') ? 'C' : 'c';

                        // Check if correction creates a valid (non-stop) codon
                        if (!CodonTable::is_stop_codon(test_codon[0], test_codon[1], test_codon[2])) {
                            corrected[i] = (original == 'T') ? 'C' : 'c';
                            corrections++;
                        }
                    }
                } else if (i + 2 < seq_len) {
                    // Fallback for edge cases
                    char c1 = (seq[i] == 'T') ? 'C' : 'c';
                    if (!CodonTable::is_stop_codon(c1, seq[i + 1], seq[i + 2])) {
                        corrected[i] = c1;
                        corrections++;
                    }
                }
            }
        }
    }

    // G->A at 3' end (last few positions)
    size_t min_3prime_pos = seq_len > 7 ? seq_len - 7 : 0;  // Extended to 7
    for (size_t i = min_3prime_pos; i < seq_len; ++i) {
        float ga_prob = profile.ga_prob_3prime[i];
        float threshold = get_threshold(i, confidence_threshold);

        // CpG context (G preceded by C)
        if (is_cpg_context(i, false)) {
            threshold *= 0.8f;
        }

        if (ga_prob > threshold && ga_prob > 0.15f) {
            if (seq[i] == 'A' || seq[i] == 'a') {
                // Frame-aware codon boundary
                size_t codon_start;
                if (frame >= 0) {
                    codon_start = frame + ((i - frame) / 3) * 3;
                    if (codon_start > i) codon_start -= 3;
                } else {
                    codon_start = (i / 3) * 3;
                }

                if (codon_start + 2 < seq_len) {
                    std::string test_codon = corrected.substr(codon_start, 3);
                    size_t pos_in_codon = i - codon_start;
                    test_codon[pos_in_codon] = (seq[i] == 'A') ? 'G' : 'g';

                    if (!CodonTable::is_stop_codon(test_codon[0], test_codon[1], test_codon[2])) {
                        corrected[i] = (seq[i] == 'A') ? 'G' : 'g';
                        corrections++;
                    }
                }
            }
        }
    }

    return {corrected, corrections};
}

// Backward-compatible overload without frame
std::pair<std::string, size_t> DamageModel::correct_damage(
    const std::string& seq,
    float confidence_threshold) const {
    return correct_damage(seq, confidence_threshold, -1);  // No frame info
}

std::pair<std::string, size_t> DamageModel::correct_protein_damage(
    const std::string& dna_seq,
    const std::string& protein_seq,
    float confidence_threshold) const {

    if (dna_seq.empty()) return {protein_seq, 0};

    // First correct the DNA sequence
    auto [corrected_dna, dna_corrections] = correct_damage(dna_seq, confidence_threshold);

    // Translate the corrected DNA to protein
    std::string corrected_protein;
    size_t aa_changes = 0;

    for (size_t i = 0; i + 2 < corrected_dna.length(); i += 3) {
        char c1 = corrected_dna[i];
        char c2 = corrected_dna[i + 1];
        char c3 = corrected_dna[i + 2];

        // Translate codon (using CodonTable)
        char aa = CodonTable::translate_codon(c1, c2, c3);

        // Compare with original protein if available
        size_t aa_pos = i / 3;
        if (aa_pos < protein_seq.length() && protein_seq[aa_pos] != aa) {
            aa_changes++;
        }

        corrected_protein += aa;
    }

    return {corrected_protein, aa_changes};
}

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

    // Estimate lambda from half-life
    float half_damage = (delta_max + delta_bg) / 2.0f;
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

            // Check if it's a stop codon
            if (CodonTable::is_stop_codon(seq[i], seq[i+1], seq[i+2])) {
                stop_count[codon_pos]++;
            }
        }
    }

    // Stop codons should be rare in coding regions
    // If enriched at 5' end, suggests C->T damage
    std::vector<float> stop_freq(25);
    for (size_t i = 0; i < 25; ++i) {
        if (total_codons[i] > 0) {
            stop_freq[i] = static_cast<float>(stop_count[i]) / total_codons[i];
        }
    }

    // Estimate damage from stop codon frequency
    // Codons like CGA->TGA, CAA->TAA, CAG->TAG are damage-sensitive
    float estimated_damage = stop_freq[0] * 3.0f;  // Rough estimate

    auto [lambda, delta_max, delta_bg] = fit_decay_model(stop_freq);

    return DamageModel::from_parameters(lambda, lambda, delta_max, delta_bg);
}

// CodonDamageAnalyzer implementation
Score CodonDamageAnalyzer::stop_codon_probability(const std::string& codon,
                                                 float ct_rate, float ga_rate) {
    if (codon.length() != 3) return 0.0f;

    // Check all possible damage scenarios
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

    // Hard-coded list of codons that become stop codons with single C->T or G->A
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

    // Sort by probability
    std::sort(vulnerable.begin(), vulnerable.end(),
             [](const auto& a, const auto& b) { return a.second > b.second; });

    return vulnerable;
}

} // namespace agp
