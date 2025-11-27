#include "agp/damage_model.hpp"
#include "agp/frame_selector.hpp"
#include "agp/hexamer_damage_lookup.hpp"
#include "agp/gtdb_hexamer_table.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cctype>
#include <iostream>
#include <atomic>

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

void DamageModel::update_from_sample_profile(const SampleDamageProfile& profile) {
    // Update decay rates from Pass 1 observations
    lambda_5prime_ = profile.lambda_5prime;
    lambda_3prime_ = profile.lambda_3prime;

    // Update maximum damage rates from terminal positions
    delta_max_ = std::max(profile.max_damage_5prime, profile.max_damage_3prime);

    // Update background rate from internal positions (use minimum observed)
    // Internal positions should have minimal damage
    float avg_internal_5 = 0.0f;
    float avg_internal_3 = 0.0f;
    for (size_t i = 10; i < 15; ++i) {
        if (profile.t_freq_5prime[i] + profile.c_freq_5prime[i] > 0) {
            avg_internal_5 += profile.t_freq_5prime[i] / (profile.t_freq_5prime[i] + profile.c_freq_5prime[i]);
        }
        if (profile.a_freq_3prime[i] + profile.g_freq_3prime[i] > 0) {
            avg_internal_3 += profile.a_freq_3prime[i] / (profile.a_freq_3prime[i] + profile.g_freq_3prime[i]);
        }
    }
    delta_background_ = std::min(avg_internal_5 / 5.0f, avg_internal_3 / 5.0f);

    // Clear profile cache so new profiles reflect updated parameters
    std::lock_guard<std::mutex> lock(cache_mutex_);
    profile_cache_ = {};
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

const DamageProfile& DamageModel::create_profile_cached(size_t seq_len) const {
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

// Sample-profile-aware damage correction
std::pair<std::string, size_t> DamageModel::correct_damage_with_profile(
    const std::string& seq,
    const SampleDamageProfile& sample_profile,
    float ancient_prob,
    int frame) const {

    if (seq.empty() || !sample_profile.is_valid()) {
        // Fallback to standard correction if no profile
        return correct_damage(seq, 0.3f, frame);
    }

    std::string corrected = seq;
    size_t corrections = 0;
    size_t seq_len = seq.length();

    // Adaptive threshold based on ancient_prob
    // Higher ancient_prob → more confidence in damage → lower threshold
    float base_threshold = 0.5f - (ancient_prob * 0.25f);  // Range: 0.25-0.5
    base_threshold = std::clamp(base_threshold, 0.20f, 0.50f);

    // Frame-aware threshold adjustment
    auto get_threshold = [&](size_t pos, float base) -> float {
        if (frame < 0) return base;
        int codon_pos = (pos - frame + 300) % 3;
        if (codon_pos == 2) {
            // Wobble position: more aggressive (synonymous changes)
            return base * 0.7f;
        } else if (codon_pos == 0) {
            // First codon position: more conservative (non-synonymous)
            return base * 1.2f;
        }
        return base;
    };

    // Get sample baseline rates for comparison
    double baseline_tc_total = sample_profile.baseline_t_freq + sample_profile.baseline_c_freq;
    double baseline_ag_total = sample_profile.baseline_a_freq + sample_profile.baseline_g_freq;

    float baseline_t_rate = (baseline_tc_total > 0) ?
        static_cast<float>(sample_profile.baseline_t_freq / baseline_tc_total) : 0.5f;
    float baseline_a_rate = (baseline_ag_total > 0) ?
        static_cast<float>(sample_profile.baseline_a_freq / baseline_ag_total) : 0.5f;


    // =========================================================================
    // 5' END: C→T damage correction using hexamer likelihood model
    // Strategy: Calculate P(corrected_hexamer) / P(observed_hexamer)
    // Only correct if corrected codon is significantly more likely
    // =========================================================================

    // Helper lambda to get hexamer context (6bp = 2 codons)
    auto get_hexamer = [&](const std::string& s, size_t pos) -> std::string {
        if (pos + 6 > s.length()) return "";
        return s.substr(pos, 6);
    };

    // Helper to calculate log probability from GTDB frequencies
    auto log_prob_hexamer = [](const std::string& hex) -> float {
        if (hex.length() != 6) return -10.0f;  // Very low probability
        uint32_t code = agp::encode_hexamer(hex.c_str());
        if (code == UINT32_MAX) return -10.0f;
        float freq = agp::GTDB_HEXAMER_FREQ[code];
        if (freq < 1e-8f) freq = 1e-8f;  // Floor to avoid log(0)
        return std::log(freq);
    };

    size_t max_5prime_pos = std::min(size_t(5), seq_len);

    for (size_t i = 0; i < max_5prime_pos; ++i) {
        char base = std::toupper(corrected[i]);
        if (base != 'T') continue;  // Only correct T→C

        // Require sample-level damage signal
        if (sample_profile.damage_rate_5prime[i] < 0.08f) {
            continue;  // Skip if sample shows < 8% damage
        }

        // Check if in CpG context (TG from CG)
        bool is_cpg = (i + 1 < seq_len && std::toupper(seq[i + 1]) == 'G');

        // Get hexamer context (align to codon boundary if possible)
        size_t hexamer_start = (i >= 3) ? ((i / 3) * 3) : 0;
        if (hexamer_start + 6 > seq_len) {
            if (seq_len < 6) continue;
            hexamer_start = seq_len - 6;
        }

        std::string observed_hex = get_hexamer(corrected, hexamer_start);
        if (observed_hex.empty()) continue;

        // Create corrected hexamer (T→C at position i)
        std::string corrected_hex = observed_hex;
        size_t pos_in_hex = i - hexamer_start;
        if (pos_in_hex >= 6) continue;
        corrected_hex[pos_in_hex] = 'C';

        // Calculate log likelihood ratio
        float log_p_observed = log_prob_hexamer(observed_hex);
        float log_p_corrected = log_prob_hexamer(corrected_hex);
        float log_likelihood_ratio = log_p_corrected - log_p_observed;

        // Damage prior: position-dependent + ancient_prob
        float position_weight = 1.0f - (i / 10.0f);  // Decay: 1.0→0.5
        float damage_prior = ancient_prob * position_weight;
        if (is_cpg) {
            damage_prior *= 1.5f;  // CpG boost
        }

        // Posterior = likelihood ratio + log prior
        float log_posterior = log_likelihood_ratio + std::log(damage_prior + 0.01f);

        // Correct if posterior favors correction
        // Require strong evidence: LLR > 1.0 (3x more likely) or very high ancient_prob
        bool high_confidence = (log_likelihood_ratio > 1.0f) || (ancient_prob > 0.9f && log_likelihood_ratio > 0.5f);

        if (high_confidence && log_posterior > -2.0f) {
            corrected[i] = (corrected[i] == 'T') ? 'C' : 'c';
            corrections++;
        }
    }

    // =========================================================================
    // 3' END: G→A damage correction using hexamer likelihood model
    // Strategy: Same as 5' end - likelihood-based correction
    // =========================================================================
    size_t max_3prime_pos = std::min(size_t(5), seq_len);

    for (size_t i = 0; i < max_3prime_pos; ++i) {
        size_t pos = seq_len - 1 - i;
        char base = std::toupper(corrected[pos]);
        if (base != 'A') continue;  // Only correct A→G

        // Require sample-level damage signal
        if (sample_profile.damage_rate_3prime[i] < 0.08f) {
            continue;
        }

        // Check if in CpG context (CA from CG on complement)
        bool is_cpg = (pos > 0 && std::toupper(seq[pos - 1]) == 'C');

        // Get hexamer context
        size_t hexamer_start;
        if (pos >= 5) {
            hexamer_start = ((pos - 2) / 3) * 3;  // Align to codon
            if (hexamer_start + 6 > seq_len) hexamer_start = seq_len - 6;
        } else {
            if (seq_len < 6) continue;
            hexamer_start = 0;
        }

        std::string observed_hex = get_hexamer(corrected, hexamer_start);
        if (observed_hex.empty()) continue;

        // Create corrected hexamer (A→G at position pos)
        std::string corrected_hex = observed_hex;
        size_t pos_in_hex = pos - hexamer_start;
        if (pos_in_hex >= 6) continue;
        corrected_hex[pos_in_hex] = 'G';

        // Calculate log likelihood ratio
        float log_p_observed = log_prob_hexamer(observed_hex);
        float log_p_corrected = log_prob_hexamer(corrected_hex);
        float log_likelihood_ratio = log_p_corrected - log_p_observed;

        // Damage prior
        float position_weight = 1.0f - (i / 10.0f);
        float damage_prior = ancient_prob * position_weight;
        if (is_cpg) {
            damage_prior *= 1.5f;
        }

        float log_posterior = log_likelihood_ratio + std::log(damage_prior + 0.01f);

        // Require strong evidence
        bool high_confidence = (log_likelihood_ratio > 1.0f) || (ancient_prob > 0.9f && log_likelihood_ratio > 0.5f);

        if (high_confidence && log_posterior > -2.0f) {
            corrected[pos] = (corrected[pos] == 'A') ? 'G' : 'g';
            corrections++;
        }
    }

    return {corrected, corrections};
}

// Score protein quality - lower is better
float DamageModel::score_protein_quality(const std::string& dna, int frame) {
    if (dna.empty()) return 1000.0f;

    float score = 0.0f;
    size_t stop_count = 0;
    size_t total_codons = 0;

    // Start position based on frame
    size_t start = (frame >= 0) ? static_cast<size_t>(frame) : 0;

    for (size_t i = start; i + 2 < dna.length(); i += 3) {
        char c1 = dna[i], c2 = dna[i + 1], c3 = dna[i + 2];
        total_codons++;

        // Heavy penalty for stop codons (especially not at the end)
        if (CodonTable::is_stop_codon(c1, c2, c3)) {
            stop_count++;
            // Heavier penalty for early stop codons
            float position_weight = 1.0f - (float(i) / dna.length());
            score += 100.0f * position_weight;
        }

        // Small bonus for common codons (indicates better protein quality)
        char aa = CodonTable::translate_codon(c1, c2, c3);
        // Common amino acids in proteins: L, A, G, V, E, S
        if (aa == 'L' || aa == 'A' || aa == 'G' || aa == 'V' || aa == 'E' || aa == 'S') {
            score -= 0.5f;
        }
        // Rare amino acids: W, C, M
        if (aa == 'W' || aa == 'C' || aa == 'M') {
            score += 0.2f;
        }
    }

    // Normalize by length
    if (total_codons > 0) {
        score = score / total_codons;
    }

    return score;
}

// Reverse strand correction: A→G at 5' end, T→C at 3' end
static std::pair<std::string, size_t> correct_damage_reverse_impl(
    const std::string& seq,
    const DamageProfile& profile,
    float confidence_threshold,
    int frame) {

    if (seq.empty()) return {seq, 0};

    std::string corrected = seq;
    size_t corrections = 0;
    size_t seq_len = seq.length();

    // Frame-aware thresholds
    auto get_threshold = [&](size_t pos, float base_threshold) -> float {
        if (frame < 0) return base_threshold;
        int codon_pos = (pos - frame + 300) % 3;
        if (codon_pos == 2) {
            return base_threshold * 0.7f;  // Wobble: more aggressive
        } else if (codon_pos == 0) {
            return base_threshold * 1.2f;  // First position: conservative
        }
        return base_threshold;
    };

    // For reverse strand reads: G→A damage appears at 5' end
    // So we correct A→G at 5' end
    size_t max_5prime_pos = std::min(size_t(7), seq_len);
    for (size_t i = 0; i < max_5prime_pos; ++i) {
        // Use the G→A profile at the 5' end (reversed pattern)
        float ga_prob = profile.ga_prob_3prime[seq_len - i - 1];  // Use 3' profile reversed
        float threshold = get_threshold(i, confidence_threshold);

        if (ga_prob > threshold && ga_prob > 0.15f) {
            if (seq[i] == 'A' || seq[i] == 'a') {
                // Check codon boundaries
                size_t codon_start;
                if (frame >= 0) {
                    codon_start = frame + ((i - frame + 300) / 3) * 3;
                    if (codon_start > i) codon_start = (codon_start >= 3) ? codon_start - 3 : 0;
                } else {
                    codon_start = (i / 3) * 3;
                }

                if (codon_start + 2 < seq_len) {
                    std::string test_codon = corrected.substr(codon_start, 3);
                    size_t pos_in_codon = i - codon_start;
                    if (pos_in_codon < 3) {
                        test_codon[pos_in_codon] = (seq[i] == 'A') ? 'G' : 'g';

                        if (!CodonTable::is_stop_codon(test_codon[0], test_codon[1], test_codon[2])) {
                            corrected[i] = (seq[i] == 'A') ? 'G' : 'g';
                            corrections++;
                        }
                    }
                }
            }
        }
    }

    // For reverse strand reads: C→T damage appears at 3' end
    // So we correct T→C at 3' end
    size_t min_3prime_pos = seq_len > 7 ? seq_len - 7 : 0;
    for (size_t i = min_3prime_pos; i < seq_len; ++i) {
        // Use the C→T profile at the 3' end (reversed pattern)
        float ct_prob = profile.ct_prob_5prime[seq_len - i - 1];  // Use 5' profile reversed
        float threshold = get_threshold(i, confidence_threshold);

        if (ct_prob > threshold && ct_prob > 0.15f) {
            if (seq[i] == 'T' || seq[i] == 't') {
                size_t codon_start;
                if (frame >= 0) {
                    codon_start = frame + ((i - frame) / 3) * 3;
                    if (codon_start > i) codon_start = (codon_start >= 3) ? codon_start - 3 : 0;
                } else {
                    codon_start = (i / 3) * 3;
                }

                if (codon_start + 2 < seq_len) {
                    std::string test_codon = corrected.substr(codon_start, 3);
                    size_t pos_in_codon = i - codon_start;
                    if (pos_in_codon < 3) {
                        test_codon[pos_in_codon] = (seq[i] == 'T') ? 'C' : 'c';

                        if (!CodonTable::is_stop_codon(test_codon[0], test_codon[1], test_codon[2])) {
                            corrected[i] = (seq[i] == 'T') ? 'C' : 'c';
                            corrections++;
                        }
                    }
                }
            }
        }
    }

    return {corrected, corrections};
}

// Dual-strand correction: try both patterns and pick the better one
std::tuple<std::string, size_t, char> DamageModel::correct_damage_dual(
    const std::string& seq,
    float confidence_threshold,
    int frame) const {

    if (seq.empty()) return {seq, 0, 'N'};

    // Get damage profile
    auto profile = create_profile(seq.length());

    // Try forward strand correction (standard: T→C at 5', A→G at 3')
    auto [forward_seq, forward_corrections] = correct_damage(seq, confidence_threshold, frame);

    // Try reverse strand correction (A→G at 5', T→C at 3')
    auto [reverse_seq, reverse_corrections] = correct_damage_reverse_impl(
        seq, profile, confidence_threshold, frame);

    // If neither made any corrections, return original
    if (forward_corrections == 0 && reverse_corrections == 0) {
        return {seq, 0, 'N'};
    }

    // If only one made corrections, use that one
    if (forward_corrections > 0 && reverse_corrections == 0) {
        return {forward_seq, forward_corrections, 'F'};
    }
    if (reverse_corrections > 0 && forward_corrections == 0) {
        return {reverse_seq, reverse_corrections, 'R'};
    }

    // Both made corrections - score each result
    float original_score = score_protein_quality(seq, frame);
    float forward_score = score_protein_quality(forward_seq, frame);
    float reverse_score = score_protein_quality(reverse_seq, frame);

    // Choose the correction that improves protein quality the most
    float forward_improvement = original_score - forward_score;
    float reverse_improvement = original_score - reverse_score;

    if (forward_improvement >= reverse_improvement && forward_score <= original_score) {
        return {forward_seq, forward_corrections, 'F'};
    } else if (reverse_score <= original_score) {
        return {reverse_seq, reverse_corrections, 'R'};
    } else if (forward_score <= original_score) {
        return {forward_seq, forward_corrections, 'F'};
    }

    // If both make things worse, return original
    return {seq, 0, 'N'};
}

// Calculate probability that a stop codon at given position is damage-induced
float DamageModel::stop_codon_damage_probability(
    const std::string& codon,
    size_t codon_start,
    size_t seq_len) const {

    if (codon.length() != 3) return 0.0f;

    // Get damage profile for this sequence length (cached to avoid O(N²))
    const auto& profile = create_profile_cached(seq_len);

    // Position of each nucleotide in the codon
    size_t pos1 = codon_start;
    size_t pos2 = codon_start + 1;
    size_t pos3 = codon_start + 2;

    char c1 = std::toupper(codon[0]);
    char c2 = std::toupper(codon[1]);
    char c3 = std::toupper(codon[2]);

    // Get damage probabilities at each position
    // For dsDNA, damage can come from either direction
    float ct_5prime_1 = (pos1 < profile.ct_prob_5prime.size()) ? profile.ct_prob_5prime[pos1] : delta_background_;
    float ct_5prime_3 = (pos3 < profile.ct_prob_5prime.size()) ? profile.ct_prob_5prime[pos3] : delta_background_;
    float ga_3prime_1 = (pos1 < profile.ga_prob_3prime.size()) ? profile.ga_prob_3prime[pos1] : delta_background_;
    float ga_3prime_3 = (pos3 < profile.ga_prob_3prime.size()) ? profile.ga_prob_3prime[pos3] : delta_background_;

    // Prior: stop codons are rare in coding regions (~0.5% expected vs 4.7% random)
    const float prior_real_stop = 0.005f;

    // Calculate P(damage | stop) for each stop codon type
    float p_damage = 0.0f;

    if (c1 == 'T' && c2 == 'A' && c3 == 'A') {
        // TAA: could come from CAA (C→T at pos 1) or GAA (G→A complement at 5')
        // C→T at position 1 (from 5' end damage)
        float p_caa_to_taa = ct_5prime_1 * 0.25f;  // P(C at pos1) * P(C→T)
        // G→A at position 1 (from 3' end on reverse strand - appears as complement)
        float p_gaa_to_taa = ga_3prime_1 * 0.25f;
        p_damage = std::max(p_caa_to_taa, p_gaa_to_taa);
    }
    else if (c1 == 'T' && c2 == 'A' && c3 == 'G') {
        // TAG: could come from CAG (C→T at pos 1)
        float p_cag_to_tag = ct_5prime_1 * 0.25f;
        p_damage = p_cag_to_tag;
    }
    else if (c1 == 'T' && c2 == 'G' && c3 == 'A') {
        // TGA: could come from CGA (C→T at pos 1) or TGG (G→A at pos 3)
        float p_cga_to_tga = ct_5prime_1 * 0.25f;
        float p_tgg_to_tga = ga_3prime_3 * 0.25f;
        p_damage = std::max(p_cga_to_tga, p_tgg_to_tga);
    }

    // Bayesian: P(damage | stop) = P(stop | damage) * P(damage) / P(stop)
    // P(stop) ≈ prior_real_stop + P(damage creates stop)
    // Simplification: if damage prob is high and position is near ends, likely damage
    float p_stop = prior_real_stop + p_damage * (1 - prior_real_stop);
    float p_damage_given_stop = (p_stop > 0) ? (p_damage / p_stop) : 0.0f;

    return std::min(1.0f, p_damage_given_stop);
}

// Bayesian stop codon correction
std::pair<std::string, size_t> DamageModel::correct_stop_codons_bayesian(
    const std::string& seq,
    int frame) const {

    if (seq.empty()) return {seq, 0};

    std::string corrected = seq;
    size_t corrections = 0;
    size_t seq_len = seq.length();

    // Threshold for correction (conservative)
    const float correction_threshold = 0.5f;

    // Start from frame offset
    size_t start = (frame >= 0) ? static_cast<size_t>(frame) : 0;

    for (size_t i = start; i + 2 < seq_len; i += 3) {
        char c1 = std::toupper(corrected[i]);
        char c2 = std::toupper(corrected[i + 1]);
        char c3 = std::toupper(corrected[i + 2]);

        // Check if this is a stop codon
        if (!CodonTable::is_stop_codon(c1, c2, c3)) continue;

        std::string codon = {c1, c2, c3};

        // Calculate probability this stop is damage-induced
        float p_damage = stop_codon_damage_probability(codon, i, seq_len);

        if (p_damage > correction_threshold) {
            // Correct the stop codon back to most likely original
            if (c1 == 'T' && c2 == 'A' && c3 == 'A') {
                // TAA → CAA (Gln) - fix C→T at position 1
                corrected[i] = (corrected[i] == 'T') ? 'C' : 'c';
                corrections++;
            }
            else if (c1 == 'T' && c2 == 'A' && c3 == 'G') {
                // TAG → CAG (Gln) - fix C→T at position 1
                corrected[i] = (corrected[i] == 'T') ? 'C' : 'c';
                corrections++;
            }
            else if (c1 == 'T' && c2 == 'G' && c3 == 'A') {
                // TGA → CGA (Arg) or TGG (Trp)
                // CGA→TGA is more common than TGG→TGA, so prefer CGA
                // But check position: if near 5' end, likely C→T; if near 3' end, likely G→A
                float dist_5prime = static_cast<float>(i);
                float dist_3prime = static_cast<float>(seq_len - i - 3);

                if (dist_5prime < dist_3prime) {
                    // Closer to 5' end: C→T damage, fix TGA → CGA
                    corrected[i] = (corrected[i] == 'T') ? 'C' : 'c';
                } else {
                    // Closer to 3' end: G→A damage, fix TGA → TGG
                    corrected[i + 2] = (corrected[i + 2] == 'A') ? 'G' : 'g';
                }
                corrections++;
            }
        }
    }

    return {corrected, corrections};
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

// Calculate amino acid probability distribution for a codon considering damage
std::array<float, 21> DamageModel::calculate_aa_probabilities(
    const std::string& codon,
    size_t codon_start,
    size_t seq_len) const {

    // Initialize all probabilities to zero
    // Indices: A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9,
    //          M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19, *=20
    std::array<float, 21> aa_probs{};

    if (codon.length() != 3) return aa_probs;

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
    // If we see T, it could have been C (damaged) or T (original)
    // If we see A, it could have been G (damaged) or A (original)
    // If we see C or G, unlikely to be damaged (damage goes C→T, G→A)

    auto get_possible_originals = [](char observed, float ct_rate, float ga_rate)
        -> std::vector<std::pair<char, float>> {
        std::vector<std::pair<char, float>> originals;
        char obs = std::toupper(observed);

        if (obs == 'T') {
            // Bayesian: P(original=X | observed=T) ∝ P(observed=T | original=X) * P(X)
            // P(observed=T | original=T) = 1 (T stays T)
            // P(observed=T | original=C) = ct_rate (C becomes T)
            // With uniform prior P(T) = P(C) = 0.25:
            // P(original=T | observed=T) = 1 / (1 + ct_rate)
            // P(original=C | observed=T) = ct_rate / (1 + ct_rate)
            float denom = 1.0f + ct_rate;
            originals.push_back({'T', 1.0f / denom});
            originals.push_back({'C', ct_rate / denom});
        }
        else if (obs == 'A') {
            // P(observed=A | original=A) = 1 (A stays A)
            // P(observed=A | original=G) = ga_rate (G becomes A)
            // P(original=A | observed=A) = 1 / (1 + ga_rate)
            // P(original=G | observed=A) = ga_rate / (1 + ga_rate)
            float denom = 1.0f + ga_rate;
            originals.push_back({'A', 1.0f / denom});
            originals.push_back({'G', ga_rate / denom});
        }
        else if (obs == 'C') {
            // C is unlikely to be damaged from something else
            // (damage goes C→T, not X→C)
            originals.push_back({'C', 1.0f});
        }
        else if (obs == 'G') {
            // G is unlikely to be damaged from something else
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
    auto orig0 = get_possible_originals(codon[0], ct_rate_0, ga_rate_0);
    auto orig1 = get_possible_originals(codon[1], ct_rate_1, ga_rate_1);
    auto orig2 = get_possible_originals(codon[2], ct_rate_2, ga_rate_2);

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
    if (total_prob > 0) {
        for (auto& p : aa_probs) {
            p /= total_prob;
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
        std::string codon = seq.substr(i, 3);

        // Calculate AA probability distribution
        auto aa_probs = calculate_aa_probabilities(codon, i, seq_len);

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

// Terminal region damage detection using GTDB hexamer frequencies
TerminalDamageResult DamageModel::detect_terminal_damage(
    const std::string& seq,
    size_t terminal_length) {

    TerminalDamageResult result = {};
    result.damage_score_5prime = 0.0f;
    result.damage_score_3prime = 0.0f;
    result.combined_damage_score = 0.0f;
    result.damaged_hexamers_5prime = 0;
    result.damaged_hexamers_3prime = 0;
    result.total_hexamers_5prime = 0;
    result.total_hexamers_3prime = 0;
    result.is_likely_damaged = false;

    if (seq.length() < 6) return result;

    // Uppercase buffer for hexamer
    char hexbuf[7];
    hexbuf[6] = '\0';

    // Analyze 5' terminal region (first terminal_length bases)
    // Look for hexamers that are rare in GTDB - suggests damage
    float damage_sum_5prime = 0.0f;
    size_t max_5prime = std::min(terminal_length, seq.length() - 5);

    for (size_t i = 0; i < max_5prime; ++i) {
        // Copy and uppercase hexamer
        bool valid = true;
        for (int j = 0; j < 6; ++j) {
            char c = std::toupper(seq[i + j]);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                valid = false;
                break;
            }
            hexbuf[j] = c;
        }
        if (!valid) continue;

        result.total_hexamers_5prime++;

        // Get GTDB damage probability for this hexamer
        float damage_prob = get_hexamer_damage_prob(hexbuf);
        if (damage_prob > 0.5f) {
            result.damaged_hexamers_5prime++;
            damage_sum_5prime += damage_prob;
        }

        // Also check if hexamer is extremely rare in GTDB
        // Rare hexamers in terminal regions suggest damage
        float gtdb_freq = get_hexamer_freq(hexbuf);
        if (gtdb_freq < 1e-6f) {
            // Very rare hexamer - check if undamaged version exists
            // Check T at position 0 (could be C→T)
            if (hexbuf[0] == 'T') {
                hexbuf[0] = 'C';
                float undamaged_freq = get_hexamer_freq(hexbuf);
                hexbuf[0] = 'T';  // restore
                if (undamaged_freq > 1e-4f) {
                    // Undamaged version is common - this is likely damage
                    damage_sum_5prime += 0.8f;
                    result.damaged_hexamers_5prime++;
                }
            }
        }
    }

    // Analyze 3' terminal region (last terminal_length bases)
    // Look for G→A damage patterns
    float damage_sum_3prime = 0.0f;
    size_t seq_len = seq.length();
    size_t start_3prime = (seq_len > terminal_length + 5) ? (seq_len - terminal_length - 5) : 0;

    for (size_t i = start_3prime; i + 5 < seq_len; ++i) {
        // Copy and uppercase hexamer
        bool valid = true;
        for (int j = 0; j < 6; ++j) {
            char c = std::toupper(seq[i + j]);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                valid = false;
                break;
            }
            hexbuf[j] = c;
        }
        if (!valid) continue;

        result.total_hexamers_3prime++;

        // Get GTDB damage probability
        float damage_prob = get_hexamer_damage_prob(hexbuf);
        if (damage_prob > 0.5f) {
            result.damaged_hexamers_3prime++;
            damage_sum_3prime += damage_prob;
        }

        // Check for rare hexamers with A that could be G→A damage
        float gtdb_freq = get_hexamer_freq(hexbuf);
        if (gtdb_freq < 1e-6f) {
            // Check A positions (could be G→A)
            for (int pos = 0; pos < 6; ++pos) {
                if (hexbuf[pos] == 'A') {
                    hexbuf[pos] = 'G';
                    float undamaged_freq = get_hexamer_freq(hexbuf);
                    hexbuf[pos] = 'A';  // restore
                    if (undamaged_freq > 1e-4f) {
                        damage_sum_3prime += 0.6f;
                        result.damaged_hexamers_3prime++;
                        break;  // Only count once per hexamer
                    }
                }
            }
        }
    }

    // Calculate scores
    if (result.total_hexamers_5prime > 0) {
        result.damage_score_5prime = damage_sum_5prime / result.total_hexamers_5prime;
    }
    if (result.total_hexamers_3prime > 0) {
        result.damage_score_3prime = damage_sum_3prime / result.total_hexamers_3prime;
    }

    // Combined score: weight 5' more heavily (C→T is more common than G→A)
    result.combined_damage_score = 0.6f * result.damage_score_5prime +
                                   0.4f * result.damage_score_3prime;

    // Threshold for "likely damaged"
    result.is_likely_damaged = (result.damaged_hexamers_5prime >= 1) ||
                               (result.combined_damage_score > 0.3f);

    return result;
}

// Enhanced stop codon damage probability using GTDB hexamer frequencies
float DamageModel::stop_codon_damage_probability_gtdb(
    const std::string& codon,
    size_t codon_start,
    const std::string& seq) const {

    if (codon.length() != 3) return 0.0f;
    if (codon_start + 5 >= seq.length()) {
        // Can't form hexamer - fall back to standard method
        return stop_codon_damage_probability(codon, codon_start, seq.length());
    }

    // Extract hexamer starting at codon position
    char hexbuf[7];
    for (int i = 0; i < 6; ++i) {
        hexbuf[i] = std::toupper(seq[codon_start + i]);
    }
    hexbuf[6] = '\0';

    // Validate hexamer
    for (int i = 0; i < 6; ++i) {
        char c = hexbuf[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            return stop_codon_damage_probability(codon, codon_start, seq.length());
        }
    }

    // Get GTDB damage probability for this hexamer
    float hexamer_damage_prob = get_hexamer_damage_prob(hexbuf);

    // Get position-dependent damage probability
    float pos_damage_prob = stop_codon_damage_probability(codon, codon_start, seq.length());

    // Combine: if hexamer indicates damage AND position is terminal, very high confidence
    // The hexamer lookup already captures that TGA*, TAG*, TAA* are rare in real coding
    if (hexamer_damage_prob > 0.9f) {
        // Strong hexamer signal - use it directly, boost with position
        return std::min(1.0f, hexamer_damage_prob + 0.1f * pos_damage_prob);
    }

    // Check hexamer rarity as additional signal
    float gtdb_freq = get_hexamer_freq(hexbuf);
    if (gtdb_freq < 1e-7f) {
        // Extremely rare hexamer starting with stop codon
        // Almost certainly damage-induced
        return 0.95f;
    }
    if (gtdb_freq < 1e-5f) {
        // Very rare - likely damage
        return 0.8f;
    }

    // Fall back to position-based estimate with hexamer boost
    return std::max(pos_damage_prob, hexamer_damage_prob);
}

} // namespace agp
