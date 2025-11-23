/**
 * Ancient DNA damage detection
 *
 * Position-weighted analysis of C->T at 5' end and G->A at 3' end
 * Sample-level aggregation for improved accuracy
 */

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include <cmath>
#include <algorithm>
#include <array>

namespace agp {

// Position-dependent damage probability (exponential decay from ends)
static float damage_prob_at_position(size_t pos, size_t seq_len, bool is_5prime) {
    const float base_rate = 0.127f;  // From gargammel analysis
    const float lambda = 0.693f;     // Decay constant

    size_t distance;
    if (is_5prime) {
        distance = pos;
    } else {
        distance = seq_len - 1 - pos;
    }

    return base_rate * std::exp(-lambda * static_cast<float>(distance));
}

float FrameSelector::estimate_ancient_prob(const std::string& seq) {
    if (seq.length() < 15) return 0.5f;

    const size_t len = seq.length();
    const size_t analyze_len = std::min(size_t(15), len / 3);
    const float lambda = 0.3f;

    // Compute baseline T/C and A/G ratios from the MIDDLE of the sequence
    size_t mid_start = len / 3;
    size_t mid_end = 2 * len / 3;
    int mid_t = 0, mid_c = 0, mid_a = 0, mid_g = 0;

    for (size_t i = mid_start; i < mid_end; ++i) {
        char b = fast_upper(seq[i]);
        switch (b) {
            case 'T': mid_t++; break;
            case 'C': mid_c++; break;
            case 'A': mid_a++; break;
            case 'G': mid_g++; break;
        }
    }

    float baseline_tc = (mid_t + mid_c > 2) ? static_cast<float>(mid_t) / (mid_t + mid_c) : 0.5f;
    float baseline_ag = (mid_a + mid_g > 2) ? static_cast<float>(mid_a) / (mid_a + mid_g) : 0.5f;

    // Analyze 5' end with position weighting
    float weighted_t_5 = 0.0f, weighted_c_5 = 0.0f;
    float total_weight_5 = 0.0f;

    for (size_t i = 0; i < analyze_len; ++i) {
        float weight = std::exp(-lambda * static_cast<float>(i));
        char b = fast_upper(seq[i]);
        if (b == 'T') weighted_t_5 += weight;
        else if (b == 'C') weighted_c_5 += weight;

        if (b == 'T' || b == 'C') total_weight_5 += weight;
    }

    // Analyze 3' end with position weighting
    float weighted_a_3 = 0.0f, weighted_g_3 = 0.0f;
    float total_weight_3 = 0.0f;

    for (size_t i = 0; i < analyze_len; ++i) {
        float weight = std::exp(-lambda * static_cast<float>(i));
        size_t pos = len - 1 - i;
        char b = fast_upper(seq[pos]);
        if (b == 'A') weighted_a_3 += weight;
        else if (b == 'G') weighted_g_3 += weight;

        if (b == 'A' || b == 'G') total_weight_3 += weight;
    }

    // Calculate weighted ratios at ends
    float t_ratio_5 = (total_weight_5 > 0.5f) ? weighted_t_5 / total_weight_5 : baseline_tc;
    float a_ratio_3 = (total_weight_3 > 0.5f) ? weighted_a_3 / total_weight_3 : baseline_ag;

    // Excess over baseline (damage signal)
    float excess_t_5 = std::max(0.0f, t_ratio_5 - baseline_tc);
    float excess_a_3 = std::max(0.0f, a_ratio_3 - baseline_ag);

    const float damage_threshold = 0.05f;

    float damage_score_5 = (excess_t_5 > damage_threshold) ?
                           std::min(1.0f, excess_t_5 / 0.25f) : 0.0f;
    float damage_score_3 = (excess_a_3 > damage_threshold) ?
                           std::min(1.0f, excess_a_3 / 0.25f) : 0.0f;

    // Terminal position check
    float terminal_bonus = 0.0f;

    int terminal_t = 0, terminal_c = 0;
    for (size_t i = 0; i < std::min(size_t(3), len); ++i) {
        char b = fast_upper(seq[i]);
        if (b == 'T') terminal_t++;
        else if (b == 'C') terminal_c++;
    }

    int terminal_a = 0, terminal_g = 0;
    for (size_t i = 0; i < std::min(size_t(3), len); ++i) {
        char b = fast_upper(seq[len - 1 - i]);
        if (b == 'A') terminal_a++;
        else if (b == 'G') terminal_g++;
    }

    if (terminal_t >= 2 && terminal_a >= 2) {
        terminal_bonus = 0.2f;
    } else if (terminal_t >= 2 || terminal_a >= 2) {
        terminal_bonus = 0.1f;
    }

    float evidence_against = 0.0f;
    if (len > 0 && fast_upper(seq[0]) == 'C') evidence_against += 0.1f;
    if (len > 0 && fast_upper(seq[len-1]) == 'G') evidence_against += 0.1f;

    float combined_damage;
    if (damage_score_5 > 0.1f && damage_score_3 > 0.1f) {
        combined_damage = 0.5f * (damage_score_5 + damage_score_3) + 0.1f;
    } else {
        combined_damage = std::max(damage_score_5, damage_score_3);
    }

    float prob = 0.4f + 0.35f * combined_damage + terminal_bonus - evidence_against;

    return std::clamp(prob, 0.1f, 0.95f);
}

float FrameSelector::compute_damage_score(const std::string& seq) {
    if (seq.length() < 6) return 0.0f;

    const size_t len = seq.length();
    float score_5prime = 0.0f;
    float score_3prime = 0.0f;

    // Position weights: pos 0 = 1.0, pos 1 = 0.5, pos 2 = 0.25
    const float weights[3] = {1.0f, 0.5f, 0.25f};

    // 5' end: count T's (from C→T damage)
    for (int i = 0; i < 3 && i < static_cast<int>(len); ++i) {
        if (fast_upper(seq[i]) == 'T') {
            score_5prime += weights[i];
        }
    }

    // 3' end: count A's (from G→A damage)
    for (int i = 0; i < 3 && i < static_cast<int>(len); ++i) {
        if (fast_upper(seq[len - 1 - i]) == 'A') {
            score_3prime += weights[i];
        }
    }

    return score_5prime + score_3prime;
}

float FrameSelector::compute_damage_percentage(
    const std::string& seq,
    const SampleDamageProfile& sample_profile) {

    if (seq.length() < 6 || !sample_profile.is_valid()) return 0.0f;

    const size_t len = seq.length();

    // Compute a damage score that compares this read to the sample distribution
    // Uses log-likelihood ratio: P(bases | damaged) / P(bases | undamaged)

    float log_lr = 0.0f;
    int informative_positions = 0;

    // 5' end: C→T damage
    // For each position, if we see T where C could be damaged:
    //   P(T | damaged) = baseline_T + (1 - baseline_T) * damage_rate
    //   P(T | undamaged) = baseline_T
    //   LR = P(T|dmg) / P(T|undmg) = 1 + damage_rate * (1 - baseline_T) / baseline_T
    // If we see C:
    //   P(C | damaged) = baseline_C * (1 - damage_rate)
    //   P(C | undamaged) = baseline_C
    //   LR = (1 - damage_rate)

    float baseline_t = static_cast<float>(sample_profile.baseline_t_freq);
    float baseline_c = static_cast<float>(sample_profile.baseline_c_freq);
    float baseline_a = static_cast<float>(sample_profile.baseline_a_freq);
    float baseline_g = static_cast<float>(sample_profile.baseline_g_freq);

    // Ensure baselines are reasonable
    if (baseline_t < 0.1f) baseline_t = 0.25f;
    if (baseline_c < 0.1f) baseline_c = 0.25f;
    if (baseline_a < 0.1f) baseline_a = 0.25f;
    if (baseline_g < 0.1f) baseline_g = 0.25f;

    // 5' end positions (C→T damage)
    for (int i = 0; i < 5 && i < static_cast<int>(len); ++i) {
        float dmg_rate = sample_profile.damage_rate_5prime[i];
        if (dmg_rate < 0.001f) dmg_rate = 0.001f;  // Minimum for numerical stability

        char base = fast_upper(seq[i]);
        if (base == 'T') {
            // T observed - evidence for damage
            // P(T|damaged) / P(T|undamaged)
            float p_t_damaged = baseline_t + baseline_c * dmg_rate;
            float p_t_undamaged = baseline_t;
            log_lr += std::log(p_t_damaged / p_t_undamaged);
            informative_positions++;
        } else if (base == 'C') {
            // C observed - evidence against damage at this position
            float p_c_damaged = baseline_c * (1.0f - dmg_rate);
            float p_c_undamaged = baseline_c;
            log_lr += std::log(p_c_damaged / p_c_undamaged);
            informative_positions++;
        }
    }

    // 3' end positions (G→A damage)
    for (int i = 0; i < 5 && i < static_cast<int>(len); ++i) {
        size_t pos = len - 1 - i;
        float dmg_rate = sample_profile.damage_rate_3prime[i];
        if (dmg_rate < 0.001f) dmg_rate = 0.001f;

        char base = fast_upper(seq[pos]);
        if (base == 'A') {
            // A observed - evidence for damage
            float p_a_damaged = baseline_a + baseline_g * dmg_rate;
            float p_a_undamaged = baseline_a;
            log_lr += std::log(p_a_damaged / p_a_undamaged);
            informative_positions++;
        } else if (base == 'G') {
            // G observed - evidence against damage
            float p_g_damaged = baseline_g * (1.0f - dmg_rate);
            float p_g_undamaged = baseline_g;
            log_lr += std::log(p_g_damaged / p_g_undamaged);
            informative_positions++;
        }
    }

    if (informative_positions == 0) return 50.0f;  // No informative positions

    // Convert log-likelihood ratio to a percentage
    // log_lr > 0 means evidence for damage
    // log_lr < 0 means evidence against damage

    // Scale: use sample's max damage rate to calibrate
    // A "fully damaged" read at max sample damage would have log_lr around:
    float max_expected_lr = 0.0f;
    for (int i = 0; i < 5; ++i) {
        float dmg_5 = sample_profile.damage_rate_5prime[i];
        float dmg_3 = sample_profile.damage_rate_3prime[i];
        if (dmg_5 > 0.001f) {
            max_expected_lr += std::log(1.0f + dmg_5 * baseline_c / baseline_t);
        }
        if (dmg_3 > 0.001f) {
            max_expected_lr += std::log(1.0f + dmg_3 * baseline_g / baseline_a);
        }
    }

    // Normalize: 0% = max evidence against, 50% = neutral, 100% = max evidence for
    float normalized;
    if (max_expected_lr > 0.01f) {
        normalized = 50.0f + 50.0f * (log_lr / max_expected_lr);
    } else {
        // Low-damage sample: use simpler scaling
        normalized = 50.0f + log_lr * 20.0f;  // Each log unit = 20%
    }

    return std::clamp(normalized, 0.0f, 100.0f);
}

float FrameSelector::estimate_ancient_prob_with_codons(
    const std::string& seq,
    int frame,
    bool forward) {

    if (seq.length() < 15) return 0.5f;

    const std::string& working_seq = seq;
    size_t len = working_seq.length();

    float score = 0.0f;

    // FEATURE 1: Single position features at read ends
    char base_0 = fast_upper(working_seq[0]);
    char base_end = fast_upper(working_seq[len - 1]);

    if (base_0 == 'T') score += 0.35f;
    else if (base_0 == 'C') score -= 0.12f;

    if (base_end == 'A') score += 0.20f;
    else if (base_end == 'G') score -= 0.10f;

    // FEATURE 2: First codon pattern
    if (len >= 3) {
        std::string first_codon;
        first_codon.reserve(3);
        for (int i = 0; i < 3; ++i) {
            first_codon += fast_upper(working_seq[i]);
        }

        if (first_codon == "TGC") score += 0.25f;
        else if (first_codon == "TCG") score += 0.24f;
        else if (first_codon == "TGA") score += 0.20f;
        else if (first_codon == "TGT") score += 0.18f;
        else if (first_codon == "TCT") score += 0.17f;
        else if (first_codon == "TTG") score += 0.17f;
        else if (first_codon == "TTC") score += 0.15f;
        else if (first_codon == "TCA") score += 0.15f;
        else if (first_codon == "TAG") score += 0.15f;
        else if (first_codon == "CAA") score -= 0.30f;
        else if (first_codon == "CCA") score -= 0.30f;
        else if (first_codon == "CAT") score -= 0.22f;
        else if (first_codon == "CCT") score -= 0.15f;
        else if (first_codon[0] == 'T') score += 0.08f;
        else if (first_codon[0] == 'C') score -= 0.08f;
    }

    // FEATURE 3: Count T's in first 5 positions
    int t_count_5 = 0, c_count_5 = 0;
    for (size_t i = 0; i < std::min(size_t(5), len); ++i) {
        char b = fast_upper(working_seq[i]);
        if (b == 'T') t_count_5++;
        else if (b == 'C') c_count_5++;
    }
    score += 0.04f * t_count_5;
    score -= 0.025f * c_count_5;

    // FEATURE 4: Count A's in last 5 positions
    int a_count_3 = 0, g_count_3 = 0;
    for (size_t i = 0; i < std::min(size_t(5), len); ++i) {
        char b = fast_upper(working_seq[len - 1 - i]);
        if (b == 'A') a_count_3++;
        else if (b == 'G') g_count_3++;
    }
    score += 0.03f * a_count_3;
    score -= 0.02f * g_count_3;

    // FEATURE 5: Last codon pattern
    if (len >= 3) {
        std::string last_codon;
        last_codon.reserve(3);
        for (size_t i = 0; i < 3; ++i) {
            last_codon += fast_upper(working_seq[len - 3 + i]);
        }

        if (last_codon[2] == 'A') score += 0.08f;
        else if (last_codon[2] == 'G') score -= 0.05f;
    }

    // Convert score to probability using sigmoid
    float prob = 1.0f / (1.0f + std::exp(-2.5f * score));
    prob = 0.15f + 0.77f * prob;

    return std::clamp(prob, 0.1f, 0.95f);
}

// ============================================================================
// SAMPLE-LEVEL DAMAGE AGGREGATION
// ============================================================================

void FrameSelector::update_sample_profile(
    SampleDamageProfile& profile,
    const std::string& seq) {

    if (seq.length() < 30) return;

    size_t len = seq.length();

    // Count bases at 5' end positions
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.t_freq_5prime[i] += 1.0;
        else if (base == 'C') profile.c_freq_5prime[i] += 1.0;
    }

    // Count bases at 3' end positions
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        if (base == 'A') profile.a_freq_3prime[i] += 1.0;
        else if (base == 'G') profile.g_freq_3prime[i] += 1.0;
    }

    // Count bases in middle (undamaged baseline)
    size_t mid_start = len / 3;
    size_t mid_end = 2 * len / 3;
    for (size_t i = mid_start; i < mid_end; ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.baseline_t_freq += 1.0;
        else if (base == 'C') profile.baseline_c_freq += 1.0;
        else if (base == 'A') profile.baseline_a_freq += 1.0;
        else if (base == 'G') profile.baseline_g_freq += 1.0;
    }

    // Codon-position-aware counting at 5' end (first 15 bases)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        int codon_pos = i % 3;  // 0, 1, 2
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.codon_pos_t_count_5prime[codon_pos]++;
        else if (base == 'C') profile.codon_pos_c_count_5prime[codon_pos]++;
    }

    // Codon-position-aware counting at 3' end (last 15 bases)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        int codon_pos = (len - 1 - i) % 3;
        char base = fast_upper(seq[pos]);
        if (base == 'A') profile.codon_pos_a_count_3prime[codon_pos]++;
        else if (base == 'G') profile.codon_pos_g_count_3prime[codon_pos]++;
    }

    // CpG context damage tracking (5' end, first 5 bases)
    // We use the same range for both contexts to avoid asymmetric counting
    for (size_t i = 0; i < std::min(size_t(5), len - 1); ++i) {
        char base = fast_upper(seq[i]);
        char next = fast_upper(seq[i + 1]);

        // Check for CpG context: C followed by G, or T followed by G (damaged CpG)
        if (next == 'G') {
            if (base == 'C') {
                profile.cpg_c_count++;
            } else if (base == 'T') {
                profile.cpg_t_count++;  // Likely C→T in CpG
            }
        } else {
            // Non-CpG context (same position range as CpG)
            if (base == 'C') {
                profile.non_cpg_c_count++;
            } else if (base == 'T') {
                profile.non_cpg_t_count++;
            }
        }
    }

    profile.n_reads++;
}

void FrameSelector::finalize_sample_profile(SampleDamageProfile& profile) {
    if (profile.n_reads == 0) return;

    // Normalize baseline frequencies (using double to preserve precision)
    double mid_total = profile.baseline_t_freq + profile.baseline_c_freq +
                       profile.baseline_a_freq + profile.baseline_g_freq;
    if (mid_total > 0) {
        profile.baseline_t_freq /= mid_total;
        profile.baseline_c_freq /= mid_total;
        profile.baseline_a_freq /= mid_total;
        profile.baseline_g_freq /= mid_total;
    }

    // Compute damage rates as C→T transition rate at 5' and G→A at 3'
    // For ancient DNA: damage converts C→T, so we measure T's where C was expected
    // Damage rate = mismatches / (matches + mismatches) where mismatch = T at C site

    // Normalize position-specific frequencies and compute damage rates
    for (int i = 0; i < 15; ++i) {
        // 5' end: C→T damage
        // The damage rate is the fraction of C+T sites that show T
        // For aDNA with 30% damage: if original had 50% C, now 35% C + 15% T from damage
        // T/(T+C) = (original_T + damaged_C) / (original_T + original_C)
        // We report T/(T+C) - 0.5 (assuming balanced genome) as the damage signal
        double tc_total = profile.t_freq_5prime[i] + profile.c_freq_5prime[i];
        if (tc_total > 0) {
            double t_freq = profile.t_freq_5prime[i] / tc_total;
            // Report C→T rate: fraction of T's at C/T sites above baseline 50%
            // For highly damaged reads, this will be high (e.g., 0.65 - 0.5 = 0.15 = 15%)
            profile.damage_rate_5prime[i] = std::max(0.0f, static_cast<float>(t_freq - 0.5));

            profile.t_freq_5prime[i] = t_freq;
            profile.c_freq_5prime[i] = 1.0 - t_freq;
        }

        // 3' end: G→A damage
        double ag_total = profile.a_freq_3prime[i] + profile.g_freq_3prime[i];
        if (ag_total > 0) {
            double a_freq = profile.a_freq_3prime[i] / ag_total;
            // Report G→A rate: fraction of A's at A/G sites above baseline 50%
            profile.damage_rate_3prime[i] = std::max(0.0f, static_cast<float>(a_freq - 0.5));

            profile.a_freq_3prime[i] = a_freq;
            profile.g_freq_3prime[i] = 1.0 - a_freq;
        }
    }

    // Compute codon-position-aware damage rates
    double baseline_tc = profile.baseline_t_freq /
                        (profile.baseline_t_freq + profile.baseline_c_freq + 0.001);
    double baseline_ag = profile.baseline_a_freq /
                        (profile.baseline_a_freq + profile.baseline_g_freq + 0.001);

    for (int p = 0; p < 3; p++) {
        // 5' end codon position rates
        size_t tc_total = profile.codon_pos_t_count_5prime[p] + profile.codon_pos_c_count_5prime[p];
        if (tc_total > 0) {
            profile.codon_pos_t_rate_5prime[p] = static_cast<float>(profile.codon_pos_t_count_5prime[p]) / tc_total;
        }

        // 3' end codon position rates
        size_t ag_total = profile.codon_pos_a_count_3prime[p] + profile.codon_pos_g_count_3prime[p];
        if (ag_total > 0) {
            profile.codon_pos_a_rate_3prime[p] = static_cast<float>(profile.codon_pos_a_count_3prime[p]) / ag_total;
        }
    }

    // Compute CpG damage rates
    size_t cpg_total = profile.cpg_c_count + profile.cpg_t_count;
    if (cpg_total > 10) {  // Need enough data
        profile.cpg_damage_rate = static_cast<float>(profile.cpg_t_count) / cpg_total;
    }

    size_t non_cpg_total = profile.non_cpg_c_count + profile.non_cpg_t_count;
    if (non_cpg_total > 10) {
        profile.non_cpg_damage_rate = static_cast<float>(profile.non_cpg_t_count) / non_cpg_total;
    }

    // Summary statistics
    profile.max_damage_5prime = profile.damage_rate_5prime[0];
    profile.max_damage_3prime = profile.damage_rate_3prime[0];

    // Estimate decay constants (lambda) from the damage profile
    // Using half-life method: lambda = ln(2) / half_position
    // Find position where damage drops to 50% of max
    auto estimate_lambda = [](const std::array<float, 15>& rates, float max_rate) -> float {
        if (max_rate < 0.02f) return 0.3f;  // Not enough signal, use default

        float half_rate = max_rate / 2.0f;
        int half_pos = 1;  // Default if decay is very steep

        for (int i = 1; i < 15; ++i) {
            if (rates[i] <= half_rate) {
                // Interpolate for more precision
                if (i > 0 && rates[i - 1] > half_rate) {
                    float frac = (rates[i - 1] - half_rate) / (rates[i - 1] - rates[i] + 0.001f);
                    half_pos = i - 1;
                    float interp_pos = static_cast<float>(half_pos) + frac;
                    return 0.693f / std::max(0.5f, interp_pos);  // ln(2) / half_pos
                }
                half_pos = i;
                break;
            }
        }

        return 0.693f / std::max(1.0f, static_cast<float>(half_pos));
    };

    profile.lambda_5prime = estimate_lambda(profile.damage_rate_5prime, profile.max_damage_5prime);
    profile.lambda_3prime = estimate_lambda(profile.damage_rate_3prime, profile.max_damage_3prime);

    // Clamp to reasonable range [0.1, 1.0]
    profile.lambda_5prime = std::clamp(profile.lambda_5prime, 0.1f, 1.0f);
    profile.lambda_3prime = std::clamp(profile.lambda_3prime, 0.1f, 1.0f);

    // Library type detection
    // Double-stranded: C→T at 5' AND G→A at 3' (complementary pattern, ~symmetric)
    // Single-stranded: Damage at only ONE end (asymmetric pattern)
    //   - 5' C→T only (no 3' G→A) - one ss orientation
    //   - 3' G→A only (no 5' C→T) - opposite ss orientation
    const float damage_threshold = 0.02f;  // 2% minimum for detection
    const float ratio_threshold = 0.3f;    // asymmetric if ratio < 0.3

    bool has_5prime = profile.max_damage_5prime > damage_threshold;
    bool has_3prime = profile.max_damage_3prime > damage_threshold;

    if (has_5prime && has_3prime) {
        // Both ends show damage - check symmetry
        float ratio = std::min(profile.max_damage_5prime, profile.max_damage_3prime) /
                      (std::max(profile.max_damage_5prime, profile.max_damage_3prime) + 0.001f);
        if (ratio > ratio_threshold) {
            // Symmetric damage - double-stranded
            profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
        } else {
            // Asymmetric - likely single-stranded
            profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
        }
    } else if (has_5prime && !has_3prime) {
        // Only 5' damage - single-stranded
        profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
    } else if (!has_5prime && has_3prime) {
        // Only 3' damage - single-stranded (opposite orientation)
        profile.library_type = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
    }
    // else: not enough damage at either end, stays UNKNOWN

    // Improved sample classification using multiple signals
    float damage_signal = (profile.max_damage_5prime + profile.max_damage_3prime) / 2.0f;

    // CpG and wobble boosts are SUPPORTIVE signals only - require baseline damage first
    // Without visible damage signal, these could trigger on natural sequence variation
    float cpg_boost = 0.0f;
    float wobble_boost = 0.0f;

    // Only consider boosts if there's already some damage signal (>1%)
    if (damage_signal > 0.01f) {
        // CpG boost: ancient DNA shows elevated CpG damage
        if (profile.cpg_damage_rate > profile.non_cpg_damage_rate + 0.05f) {
            cpg_boost = 0.03f;  // Reduced boost - supportive evidence only
        }

        // Wobble position boost: ancient coding DNA shows pos3 enrichment
        float wobble_enrichment = profile.codon_pos_t_rate_5prime[2] -
                                  (profile.codon_pos_t_rate_5prime[0] + profile.codon_pos_t_rate_5prime[1]) / 2.0f;
        if (wobble_enrichment > 0.05f) {
            wobble_boost = 0.02f;  // Reduced boost - supportive evidence only
        }
    }

    float total_signal = damage_signal + cpg_boost + wobble_boost;

    if (total_signal > 0.12f) {
        profile.sample_damage_prob = 0.95f;
    } else if (total_signal > 0.06f) {
        profile.sample_damage_prob = 0.80f;
    } else if (total_signal > 0.03f) {
        profile.sample_damage_prob = 0.50f;
    } else {
        profile.sample_damage_prob = 0.20f;
    }
}

void FrameSelector::merge_sample_profiles(SampleDamageProfile& dst, const SampleDamageProfile& src) {
    // Merge position-specific counts (before normalization)
    for (int i = 0; i < 15; ++i) {
        dst.t_freq_5prime[i] += src.t_freq_5prime[i];
        dst.c_freq_5prime[i] += src.c_freq_5prime[i];
        dst.a_freq_3prime[i] += src.a_freq_3prime[i];
        dst.g_freq_3prime[i] += src.g_freq_3prime[i];
    }

    // Merge baseline counts
    dst.baseline_t_freq += src.baseline_t_freq;
    dst.baseline_c_freq += src.baseline_c_freq;
    dst.baseline_a_freq += src.baseline_a_freq;
    dst.baseline_g_freq += src.baseline_g_freq;

    // Merge codon position counts
    for (int p = 0; p < 3; ++p) {
        dst.codon_pos_t_count_5prime[p] += src.codon_pos_t_count_5prime[p];
        dst.codon_pos_c_count_5prime[p] += src.codon_pos_c_count_5prime[p];
        dst.codon_pos_a_count_3prime[p] += src.codon_pos_a_count_3prime[p];
        dst.codon_pos_g_count_3prime[p] += src.codon_pos_g_count_3prime[p];
    }

    // Merge CpG counts
    dst.cpg_c_count += src.cpg_c_count;
    dst.cpg_t_count += src.cpg_t_count;
    dst.non_cpg_c_count += src.non_cpg_c_count;
    dst.non_cpg_t_count += src.non_cpg_t_count;

    // Merge read count
    dst.n_reads += src.n_reads;
}

SampleDamageProfile FrameSelector::compute_sample_profile(
    const std::vector<std::string>& sequences) {

    SampleDamageProfile profile;

    for (const auto& seq : sequences) {
        update_sample_profile(profile, seq);
    }

    finalize_sample_profile(profile);
    return profile;
}

float FrameSelector::estimate_ancient_prob_with_sample_profile(
    const std::string& seq,
    const SampleDamageProfile& sample_profile) {

    if (seq.length() < 15) return 0.5f;
    if (!sample_profile.is_valid()) {
        return estimate_ancient_prob_with_codons(seq, 0, true);
    }

    size_t len = seq.length();
    float log_odds = 0.0f;

    // Sample-calibrated position weights
    float damage_weight_5prime = 0.0f;
    for (size_t i = 0; i < std::min(size_t(5), len); ++i) {
        char base = fast_upper(seq[i]);
        float pos_damage_rate = sample_profile.damage_rate_5prime[i];

        if (pos_damage_rate > 0.02f) {
            if (base == 'T') {
                float weight = pos_damage_rate / (sample_profile.max_damage_5prime + 0.01f);
                damage_weight_5prime += weight * 0.3f;
            }
            else if (base == 'C') {
                float weight = pos_damage_rate / (sample_profile.max_damage_5prime + 0.01f);
                damage_weight_5prime -= weight * 0.1f;
            }
        }
    }

    float damage_weight_3prime = 0.0f;
    for (size_t i = 0; i < std::min(size_t(5), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        float pos_damage_rate = sample_profile.damage_rate_3prime[i];

        if (pos_damage_rate > 0.02f) {
            if (base == 'A') {
                float weight = pos_damage_rate / (sample_profile.max_damage_3prime + 0.01f);
                damage_weight_3prime += weight * 0.3f;
            }
            else if (base == 'G') {
                float weight = pos_damage_rate / (sample_profile.max_damage_3prime + 0.01f);
                damage_weight_3prime -= weight * 0.1f;
            }
        }
    }

    log_odds += damage_weight_5prime + damage_weight_3prime;

    // Simple position signals
    if (len > 0 && fast_upper(seq[0]) == 'T') log_odds += 0.7f;
    if (len > 1 && fast_upper(seq[1]) == 'T') log_odds += 0.25f;
    if (len > 0 && fast_upper(seq[len - 1]) == 'A') log_odds += 0.6f;
    if (len > 1 && fast_upper(seq[len - 2]) == 'A') log_odds += 0.2f;

    // Negative signals
    if (len > 0 && fast_upper(seq[len - 1]) == 'G') log_odds -= 0.2f;
    if (len > 0 && fast_upper(seq[0]) == 'C') log_odds -= 0.2f;

    float prob = 1.0f / (1.0f + std::exp(-log_odds));

    return std::clamp(prob, 0.05f, 0.95f);
}

float FrameSelector::estimate_ancient_prob_advanced(
    const std::string& seq,
    const std::string& quality,
    const SampleDamageProfile& sample_profile) {

    if (seq.length() < 15) return 0.5f;

    size_t len = seq.length();
    bool has_quality = !quality.empty() && quality.length() == seq.length();

    float log_odds = 0.0f;

    // PRIMARY SIGNAL: Position-based damage detection
    if (len > 0 && fast_upper(seq[0]) == 'T') log_odds += 0.7f;
    if (len > 1 && fast_upper(seq[1]) == 'T') log_odds += 0.25f;
    if (len > 0 && fast_upper(seq[len - 1]) == 'A') log_odds += 0.6f;
    if (len > 1 && fast_upper(seq[len - 2]) == 'A') log_odds += 0.2f;

    // Negative signals
    if (len > 0 && fast_upper(seq[0]) == 'C') log_odds -= 0.2f;
    if (len > 0 && fast_upper(seq[len - 1]) == 'G') log_odds -= 0.15f;

    // FEATURE 1: Quality-weighted base observations
    auto phred_to_prob = [](char q) -> float {
        int phred = static_cast<int>(q) - 33;
        phred = std::clamp(phred, 0, 40);
        return 1.0f - std::pow(10.0f, -phred / 10.0f);
    };

    for (size_t i = 0; i < std::min(size_t(3), len); ++i) {
        char base = fast_upper(seq[i]);
        float quality_weight = has_quality ? phred_to_prob(quality[i]) : 0.9f;

        if (base == 'T') {
            float quality_bonus = (quality_weight < 0.95f) ? 0.1f : 0.0f;
            log_odds += quality_bonus;
        }
        else if (base == 'C' && quality_weight > 0.99f) {
            float pos_weight = 1.0f - (float)i / 3.0f;
            log_odds -= pos_weight * 0.05f;
        }
    }

    for (size_t i = 0; i < std::min(size_t(3), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        float quality_weight = has_quality ? phred_to_prob(quality[pos]) : 0.9f;

        if (base == 'A') {
            float quality_bonus = (quality_weight < 0.95f) ? 0.08f : 0.0f;
            log_odds += quality_bonus;
        }
        else if (base == 'G' && quality_weight > 0.99f) {
            float pos_weight = 1.0f - (float)i / 3.0f;
            log_odds -= pos_weight * 0.04f;
        }
    }

    // FEATURE 2: Dinucleotide context
    for (size_t i = 0; i < std::min(size_t(4), len - 1); ++i) {
        char b1 = fast_upper(seq[i]);
        char b2 = fast_upper(seq[i + 1]);

        if (b1 == 'T' && b2 == 'G') log_odds += 0.3f;
        else if (b1 == 'T' && b2 == 'A') log_odds += 0.15f;
    }

    for (size_t i = len - 1; i > 0 && i > len - 5; --i) {
        char b1 = fast_upper(seq[i - 1]);
        char b2 = fast_upper(seq[i]);

        if (b1 == 'C' && b2 == 'A') log_odds += 0.25f;
    }

    // FEATURE 3: Joint probability
    bool t_at_5prime = (len > 0 && fast_upper(seq[0]) == 'T');
    bool a_at_3prime = (len > 0 && fast_upper(seq[len - 1]) == 'A');

    if (t_at_5prime && a_at_3prime) log_odds += 0.3f;

    int t_count_5prime = 0;
    for (size_t i = 0; i < std::min(size_t(3), len); ++i) {
        if (fast_upper(seq[i]) == 'T') t_count_5prime++;
    }
    if (t_count_5prime >= 2) log_odds += 0.25f;

    int a_count_3prime = 0;
    for (size_t i = 0; i < std::min(size_t(3), len); ++i) {
        if (fast_upper(seq[len - 1 - i]) == 'A') a_count_3prime++;
    }
    if (a_count_3prime >= 2) log_odds += 0.2f;

    // FEATURE 4: Sample-calibrated weights
    if (sample_profile.is_valid() && sample_profile.max_damage_5prime > 0.02f) {
        float sample_damage_scale = sample_profile.max_damage_5prime / 0.1f;
        sample_damage_scale = std::clamp(sample_damage_scale, 0.5f, 2.0f);
        log_odds *= sample_damage_scale;
    }

    // FEATURE 5: Negative evidence
    bool no_damage_signature = true;
    for (size_t i = 0; i < std::min(size_t(3), len); ++i) {
        if (fast_upper(seq[i]) == 'T') no_damage_signature = false;
        if (fast_upper(seq[len - 1 - i]) == 'A') no_damage_signature = false;
    }
    if (no_damage_signature) log_odds -= 0.3f;

    float prob = 1.0f / (1.0f + std::exp(-log_odds));

    return std::clamp(prob, 0.05f, 0.95f);
}

// ============================================================================
// CODON-AWARE DAMAGE DETECTION
// ============================================================================

// Check if a codon could be a damaged stop codon
// Returns: 0 = not a damaged stop, 1 = possible, 2 = likely
static int is_damaged_stop_codon(char c1, char c2, char c3) {
    // TAA could be from CAA (Gln)
    if (c1 == 'T' && c2 == 'A' && c3 == 'A') return 2;
    // TAG could be from CAG (Gln)
    if (c1 == 'T' && c2 == 'A' && c3 == 'G') return 2;
    // TGA could be from CGA (Arg)
    if (c1 == 'T' && c2 == 'G' && c3 == 'A') return 2;
    // TGG could be from CGG (Arg) - not stop but damage pattern
    if (c1 == 'T' && c2 == 'G' && c3 == 'G') return 1;
    return 0;
}

// Check if codon shows C→T damage at position 1
static bool has_ct_damage_pos1(char c1, char c2, char c3) {
    // T at position 1 where C would make sense
    if (c1 != 'T') return false;

    // TNN patterns that could be CNN
    // Common amino acids from CNN: Pro (CCN), Leu (CTN), His/Gln (CAN), Arg (CGN)
    // If we see TNN that's unusual, likely damage
    return true;
}

float FrameSelector::score_damage_frame_consistency(
    const std::string& seq,
    int frame) {

    if (seq.length() < 15) return 0.5f;

    size_t len = seq.length();
    float score = 0.0f;

    // Count T's at each codon position in first 15 bases
    int t_at_pos1 = 0, t_at_pos2 = 0, t_at_pos3 = 0;
    int c_at_pos1 = 0, c_at_pos2 = 0, c_at_pos3 = 0;
    int total_codons = 0;

    for (size_t i = frame; i + 2 < std::min(size_t(18), len); i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        if (c1 == 'T') t_at_pos1++;
        else if (c1 == 'C') c_at_pos1++;

        if (c2 == 'T') t_at_pos2++;
        else if (c2 == 'C') c_at_pos2++;

        if (c3 == 'T') t_at_pos3++;
        else if (c3 == 'C') c_at_pos3++;

        total_codons++;
    }

    if (total_codons == 0) return 0.5f;

    // In coding regions with damage:
    // - T at position 3 (wobble) is often tolerated (synonymous)
    // - T at position 1 or 2 usually changes amino acid

    // Calculate T/(T+C) ratio at each position
    float tc_ratio_pos1 = (t_at_pos1 + c_at_pos1 > 0) ?
        (float)t_at_pos1 / (t_at_pos1 + c_at_pos1) : 0.5f;
    float tc_ratio_pos3 = (t_at_pos3 + c_at_pos3 > 0) ?
        (float)t_at_pos3 / (t_at_pos3 + c_at_pos3) : 0.5f;

    // If T is enriched at position 3 relative to position 1,
    // this frame is more consistent with real coding + damage
    // (because wobble position tolerates C→T better)
    float wobble_bias = tc_ratio_pos3 - tc_ratio_pos1;

    // Check for damaged stop codons at 5' end
    int damaged_stops = 0;
    for (size_t i = frame; i + 2 < std::min(size_t(12), len); i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        damaged_stops += is_damaged_stop_codon(c1, c2, c3);
    }

    // Similarly check 3' end for G→A damage patterns
    // A at position 3 where G would be expected
    int a_at_pos3_3prime = 0, g_at_pos3_3prime = 0;

    size_t end_start = (len > 18) ? len - 18 : 0;
    // Align to frame
    while ((end_start - frame) % 3 != 0 && end_start < len) end_start++;

    for (size_t i = end_start; i + 2 < len; i += 3) {
        char c3 = fast_upper(seq[i + 2]);
        if (c3 == 'A') a_at_pos3_3prime++;
        else if (c3 == 'G') g_at_pos3_3prime++;
    }

    float ag_ratio_pos3_3prime = (a_at_pos3_3prime + g_at_pos3_3prime > 0) ?
        (float)a_at_pos3_3prime / (a_at_pos3_3prime + g_at_pos3_3prime) : 0.5f;

    // Combine signals
    score = 0.5f;
    score += wobble_bias * 0.3f;  // Wobble position bias
    score += damaged_stops * 0.15f;  // Damaged stop codons
    score += (ag_ratio_pos3_3prime - 0.5f) * 0.2f;  // 3' A enrichment at wobble

    return std::clamp(score, 0.0f, 1.0f);
}

float FrameSelector::estimate_ancient_prob_codon_aware(
    const std::string& seq,
    int frame,
    const SampleDamageProfile& sample_profile) {

    if (seq.length() < 15) return 0.5f;

    size_t len = seq.length();
    float log_odds = 0.0f;

    // =========================================================================
    // FEATURE 1: Terminal base patterns (standard damage signal)
    // =========================================================================
    if (fast_upper(seq[0]) == 'T') log_odds += 0.6f;
    else if (fast_upper(seq[0]) == 'C') log_odds -= 0.15f;

    if (fast_upper(seq[len - 1]) == 'A') log_odds += 0.5f;
    else if (fast_upper(seq[len - 1]) == 'G') log_odds -= 0.12f;

    // =========================================================================
    // FEATURE 1b: Terminal DINUCLEOTIDE patterns (stronger damage signal)
    // CT→TT, CG→TG at 5' end; GA→AA, CA→TA at 3' end
    // =========================================================================
    if (len >= 2) {
        char b0 = fast_upper(seq[0]);
        char b1 = fast_upper(seq[1]);

        // 5' dinucleotide patterns indicative of C→T damage
        if (b0 == 'T' && b1 == 'T') log_odds += 0.2f;      // TT from CT
        if (b0 == 'T' && b1 == 'G') log_odds += 0.35f;     // TG from CG (CpG damage!)
        if (b0 == 'T' && b1 == 'A') log_odds += 0.15f;     // TA from CA

        char bm1 = fast_upper(seq[len - 1]);
        char bm2 = fast_upper(seq[len - 2]);

        // 3' dinucleotide patterns indicative of G→A damage
        if (bm1 == 'A' && bm2 == 'A') log_odds += 0.2f;    // AA from GA
        if (bm1 == 'A' && bm2 == 'C') log_odds += 0.3f;    // CA from CG on complement
        if (bm1 == 'A' && bm2 == 'T') log_odds += 0.15f;   // TA from TG on complement
    }

    // =========================================================================
    // FEATURE 1c: Joint 5' AND 3' damage - BONUS only (no penalty)
    // Only ~6% of ancient reads have damage at BOTH ends, so don't penalize
    // =========================================================================
    bool has_5prime_signal = (fast_upper(seq[0]) == 'T');
    bool has_3prime_signal = (fast_upper(seq[len-1]) == 'A');

    if (has_5prime_signal && has_3prime_signal) {
        log_odds += 0.4f;  // Damage at BOTH ends - very strong signal
    }
    // No penalty for single-end or no damage - too common in real ancient DNA

    // =========================================================================
    // FEATURE 1d: T/C ratio in first 5 bases (gradient damage)
    // Damage decays with distance from ends - check positions 0-4
    // =========================================================================
    if (len >= 5) {
        int t_5prime = 0, c_5prime = 0;
        for (size_t i = 0; i < 5; ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'T') t_5prime++;
            else if (b == 'C') c_5prime++;
        }
        // Elevated T/(T+C) ratio suggests damage even if position 0 isn't T
        if (t_5prime + c_5prime >= 2) {
            float tc_ratio = (float)t_5prime / (t_5prime + c_5prime);
            if (tc_ratio > 0.65f) log_odds += 0.25f;  // >65% T suggests damage
            else if (tc_ratio > 0.55f) log_odds += 0.1f;  // Modest elevation
        }

        // Same for 3' end (A/G ratio)
        int a_3prime = 0, g_3prime = 0;
        for (size_t i = len - 5; i < len; ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'A') a_3prime++;
            else if (b == 'G') g_3prime++;
        }
        if (a_3prime + g_3prime >= 2) {
            float ag_ratio = (float)a_3prime / (a_3prime + g_3prime);
            if (ag_ratio > 0.65f) log_odds += 0.2f;
            else if (ag_ratio > 0.55f) log_odds += 0.08f;
        }
    }

    // =========================================================================
    // FEATURE 2: Codon-position-aware damage at 5' end
    // C→T damage should show codon position pattern in real coding regions
    // =========================================================================

    // Count T and C at each codon position in first ~15 bases
    std::array<int, 3> t_count = {0, 0, 0};
    std::array<int, 3> c_count = {0, 0, 0};

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        int codon_pos = (i - frame + 300) % 3;  // +300 to handle negative
        char base = fast_upper(seq[i]);

        if (base == 'T') t_count[codon_pos]++;
        else if (base == 'C') c_count[codon_pos]++;
    }

    // Calculate T/(T+C) ratio at each codon position
    std::array<float, 3> tc_ratio;
    for (int p = 0; p < 3; p++) {
        tc_ratio[p] = (t_count[p] + c_count[p] > 0) ?
            (float)t_count[p] / (t_count[p] + c_count[p]) : 0.5f;
    }

    // Key insight: In real coding regions with C→T damage,
    // position 3 (wobble) should show MORE T enrichment because
    // synonymous changes are tolerated/preserved
    float wobble_enrichment = tc_ratio[2] - (tc_ratio[0] + tc_ratio[1]) / 2.0f;

    if (wobble_enrichment > 0.1f) {
        log_odds += 0.4f * wobble_enrichment;  // Strong signal of coding + damage
    }

    // =========================================================================
    // FEATURE 3: Damaged stop codon detection
    // TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
    // =========================================================================

    int damaged_stop_score = 0;
    for (size_t i = frame; i + 2 < std::min(size_t(15), len); i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        int ds = is_damaged_stop_codon(c1, c2, c3);
        if (ds > 0 && i < 9) {  // Stronger signal if near 5' end
            damaged_stop_score += ds * (3 - i / 3);  // Weight by position
        }
    }

    log_odds += damaged_stop_score * 0.2f;

    // =========================================================================
    // FEATURE 4: T at codon position 1 patterns
    // TNN where CNN would be common (Pro, Leu, His, Gln, Arg)
    // =========================================================================

    int t_pos1_damage_patterns = 0;
    for (size_t i = frame; i + 2 < std::min(size_t(12), len); i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        if (c1 == 'T') {
            // TCA/TCG/TCC/TCT could be from SCA/SCG... (Ser from damage)
            // But also check for patterns that suggest C→T
            if (c2 == 'C' || c2 == 'G' || c2 == 'A') {
                t_pos1_damage_patterns++;
            }
        }
    }

    if (t_pos1_damage_patterns >= 2) {
        log_odds += 0.15f;
    }

    // =========================================================================
    // FEATURE 5: 3' end G→A damage at codon position 3 (wobble)
    // =========================================================================

    std::array<int, 3> a_count_3prime = {0, 0, 0};
    std::array<int, 3> g_count_3prime = {0, 0, 0};

    size_t end_region = (len > 15) ? len - 15 : 0;
    for (size_t i = end_region; i < len; ++i) {
        int codon_pos = (i - frame + 300) % 3;
        char base = fast_upper(seq[i]);

        if (base == 'A') a_count_3prime[codon_pos]++;
        else if (base == 'G') g_count_3prime[codon_pos]++;
    }

    // A enrichment at position 3 (wobble) at 3' end
    float ag_ratio_wobble = (a_count_3prime[2] + g_count_3prime[2] > 0) ?
        (float)a_count_3prime[2] / (a_count_3prime[2] + g_count_3prime[2]) : 0.5f;

    if (ag_ratio_wobble > 0.6f) {
        log_odds += 0.25f;
    }

    // =========================================================================
    // FEATURE 6: Sample-level calibration (CRITICAL for reducing false positives)
    // =========================================================================

    if (sample_profile.is_valid()) {
        // Scale by observed sample damage level
        float sample_scale = sample_profile.max_damage_5prime / 0.1f;
        sample_scale = std::clamp(sample_scale, 0.1f, 2.0f);  // Allow lower minimum

        // Scale the codon-specific signals by sample damage level
        log_odds *= sample_scale;

        // Sample-level prior: CRUCIAL for false positive reduction
        if (sample_profile.sample_damage_prob > 0.7f) {
            log_odds += 0.3f;  // High-damage sample: boost ancient signal
        } else if (sample_profile.sample_damage_prob < 0.3f) {
            // MODERN sample: apply strong penalty
            // Per-read signals without sample-level support are likely noise
            log_odds -= 0.8f;
        } else {
            // Mixed/uncertain sample: slight penalty
            log_odds -= 0.2f;
        }
    }

    // =========================================================================
    // FEATURE 7: Joint 5' + 3' pattern
    // =========================================================================

    bool strong_5prime = (fast_upper(seq[0]) == 'T' && t_count[2] > t_count[0]);
    bool strong_3prime = (fast_upper(seq[len-1]) == 'A' && ag_ratio_wobble > 0.5f);

    if (strong_5prime && strong_3prime) {
        log_odds += 0.35f;  // Both ends show frame-consistent damage
    }

    // Convert to probability
    float prob = 1.0f / (1.0f + std::exp(-log_odds));

    return std::clamp(prob, 0.05f, 0.95f);
}

// ============================================================================
// BAYESIAN LOG-LIKELIHOOD MODEL
// ============================================================================

float FrameSelector::compute_damage_log_likelihood(
    const std::string& seq,
    const std::string& quality,
    int frame,
    const SampleDamageProfile& sample_profile) {

    if (seq.length() < 15) return 0.0f;

    size_t len = seq.length();
    bool has_quality = !quality.empty() && quality.length() == seq.length();

    float log_lr = 0.0f;  // Log likelihood ratio: log(P(seq|ancient) / P(seq|modern))

    // Get damage rates from sample profile
    float p_damage_5prime = sample_profile.is_valid() ?
        sample_profile.max_damage_5prime : 0.15f;  // Default ancient rate
    float p_damage_3prime = sample_profile.is_valid() ?
        sample_profile.max_damage_3prime : 0.12f;

    // Background rates (modern DNA)
    float p_modern = 0.01f;  // Very low damage in modern DNA

    // Position decay constants (use sample-specific if available)
    const float lambda_5 = sample_profile.is_valid() ? sample_profile.lambda_5prime : 0.3f;
    const float lambda_3 = sample_profile.is_valid() ? sample_profile.lambda_3prime : 0.25f;

    // Phred to error probability
    auto phred_to_error = [](char q) -> float {
        int phred = static_cast<int>(q) - 33;
        phred = std::clamp(phred, 0, 40);
        return std::pow(10.0f, -phred / 10.0f);
    };

    // =========================================================================
    // 5' END: C→T damage likelihood
    // =========================================================================
    for (size_t i = 0; i < std::min(size_t(10), len); ++i) {
        char base = fast_upper(seq[i]);
        float pos_weight = std::exp(-lambda_5 * static_cast<float>(i));
        float p_ancient = p_damage_5prime * pos_weight + p_modern;
        p_ancient = std::clamp(p_ancient, 0.01f, 0.5f);

        // Quality adjustment: low quality = less informative, high quality C = strong evidence against damage
        float quality_weight = 1.0f;
        if (has_quality) {
            float error_prob = phred_to_error(quality[i]);
            // Low quality bases are less informative (could be miscall)
            if (error_prob > 0.1f) {
                quality_weight = 0.5f;  // Reduce weight for unreliable bases
            }
            // High quality C is strong evidence against damage
            if (base == 'C' && error_prob < 0.01f) {
                quality_weight = 1.5f;  // Strengthen negative evidence
            }
        }

        if (base == 'T') {
            // T observed: P(T|ancient,C→T) vs P(T|modern)
            // In CpG context, damage is more likely
            bool in_cpg = (i + 1 < len && fast_upper(seq[i + 1]) == 'G');
            float cpg_boost = in_cpg ? 1.5f : 1.0f;

            // Codon position: wobble (pos 3) tolerates damage better
            int codon_pos = (i - frame + 300) % 3;
            float codon_boost = (codon_pos == 2) ? 1.2f : 1.0f;  // Wobble

            float lr = (p_ancient * cpg_boost * codon_boost) /
                       (0.25f + p_modern);  // 0.25 = background T frequency
            log_lr += std::log(lr) * pos_weight * quality_weight;
        }
        else if (base == 'C') {
            // C observed: evidence AGAINST damage at this position
            // High quality C is especially strong evidence
            float lr = (1.0f - p_ancient) / (1.0f - p_modern);
            log_lr += std::log(lr) * pos_weight * quality_weight;
        }
    }

    // =========================================================================
    // 3' END: G→A damage likelihood
    // =========================================================================
    for (size_t i = 0; i < std::min(size_t(10), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        float pos_weight = std::exp(-lambda_3 * static_cast<float>(i));
        float p_ancient = p_damage_3prime * pos_weight + p_modern;
        p_ancient = std::clamp(p_ancient, 0.01f, 0.5f);

        // Quality adjustment: low quality = less informative, high quality G = strong evidence against damage
        float quality_weight = 1.0f;
        if (has_quality) {
            float error_prob = phred_to_error(quality[pos]);
            if (error_prob > 0.1f) {
                quality_weight = 0.5f;  // Reduce weight for unreliable bases
            }
            if (base == 'G' && error_prob < 0.01f) {
                quality_weight = 1.5f;  // Strengthen negative evidence
            }
        }

        if (base == 'A') {
            // A observed at 3' end
            int codon_pos = (pos - frame + 300) % 3;
            float codon_boost = (codon_pos == 2) ? 1.2f : 1.0f;

            float lr = (p_ancient * codon_boost) / (0.25f + p_modern);
            log_lr += std::log(lr) * pos_weight * quality_weight;
        }
        else if (base == 'G') {
            // G observed: evidence AGAINST damage at this position
            float lr = (1.0f - p_ancient) / (1.0f - p_modern);
            log_lr += std::log(lr) * pos_weight * quality_weight;
        }
    }

    // =========================================================================
    // JOINT 5' + 3' PATTERN (both ends damaged is strong evidence)
    // =========================================================================
    if (fast_upper(seq[0]) == 'T' && fast_upper(seq[len - 1]) == 'A') {
        log_lr += 0.5f;  // Bonus for joint pattern
    }

    // =========================================================================
    // SAMPLE-LEVEL PRIOR
    // =========================================================================
    if (sample_profile.is_valid()) {
        // Add prior based on sample classification
        float prior_odds = sample_profile.sample_damage_prob /
                          (1.0f - sample_profile.sample_damage_prob + 0.01f);
        log_lr += 0.3f * std::log(prior_odds);
    }

    return log_lr;
}

// ============================================================================
// QUALITY-AWARE DAMAGE DETECTION
// ============================================================================

float FrameSelector::estimate_ancient_prob_with_quality(
    const std::string& seq,
    const std::string& quality,
    int frame,
    const SampleDamageProfile& sample_profile) {

    if (seq.length() < 15) return 0.5f;

    // Use Bayesian log-likelihood
    float log_lr = compute_damage_log_likelihood(seq, quality, frame, sample_profile);

    // Convert log-likelihood ratio to probability
    // P(ancient|seq) = 1 / (1 + exp(-log_lr))
    float prob = 1.0f / (1.0f + std::exp(-log_lr));

    return std::clamp(prob, 0.05f, 0.95f);
}

} // namespace agp
