/**
 * Ancient DNA damage detection
 *
 * Position-weighted analysis of C->T at 5' end and G->A at 3' end
 * Sample-level aggregation for improved accuracy
 */

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include "agp/hexamer_damage_lookup.hpp"
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

    // Compute baseline T/C and A/G ratios from the middle of reads (already normalized)
    // This is more accurate than assuming 50% as genomes have varying GC content
    double baseline_tc = profile.baseline_t_freq /
                        (profile.baseline_t_freq + profile.baseline_c_freq + 0.001);
    double baseline_ag = profile.baseline_a_freq /
                        (profile.baseline_a_freq + profile.baseline_g_freq + 0.001);

    // Normalize position-specific frequencies and compute damage rates
    for (int i = 0; i < 15; ++i) {
        // 5' end: C→T damage
        // The damage rate is the fraction of C+T sites that show T
        // For aDNA with 30% damage: if original had 50% C, now 35% C + 15% T from damage
        // T/(T+C) = (original_T + damaged_C) / (original_T + original_C)
        // We report T/(T+C) - baseline_tc as the damage signal (using actual genome baseline)
        double tc_total = profile.t_freq_5prime[i] + profile.c_freq_5prime[i];
        if (tc_total > 0) {
            double t_freq = profile.t_freq_5prime[i] / tc_total;
            // Report C→T rate: fraction of T's at C/T sites above the sample baseline
            // For highly damaged reads, this will be high (e.g., 0.65 - 0.5 = 0.15 = 15%)
            profile.damage_rate_5prime[i] = std::max(0.0f, static_cast<float>(t_freq - baseline_tc));

            profile.t_freq_5prime[i] = t_freq;
            profile.c_freq_5prime[i] = 1.0 - t_freq;
        }

        // 3' end: G→A damage
        double ag_total = profile.a_freq_3prime[i] + profile.g_freq_3prime[i];
        if (ag_total > 0) {
            double a_freq = profile.a_freq_3prime[i] / ag_total;
            // Report G→A rate: fraction of A's at A/G sites above the sample baseline
            profile.damage_rate_3prime[i] = std::max(0.0f, static_cast<float>(a_freq - baseline_ag));

            profile.a_freq_3prime[i] = a_freq;
            profile.g_freq_3prime[i] = 1.0 - a_freq;
        }
    }

    // Compute codon-position-aware damage rates (using baseline_tc and baseline_ag from above)
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

void FrameSelector::update_sample_profile_weighted(
    SampleDamageProfile& profile,
    const std::string& seq,
    float weight) {

    if (seq.length() < 30 || weight < 0.001f) return;

    size_t len = seq.length();

    // Count bases at 5' end positions (weighted)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.t_freq_5prime[i] += weight;
        else if (base == 'C') profile.c_freq_5prime[i] += weight;
    }

    // Count bases at 3' end positions (weighted)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        if (base == 'A') profile.a_freq_3prime[i] += weight;
        else if (base == 'G') profile.g_freq_3prime[i] += weight;
    }

    // Count bases in middle (undamaged baseline) - weighted
    size_t mid_start = len / 3;
    size_t mid_end = 2 * len / 3;
    for (size_t i = mid_start; i < mid_end; ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.baseline_t_freq += weight;
        else if (base == 'C') profile.baseline_c_freq += weight;
        else if (base == 'A') profile.baseline_a_freq += weight;
        else if (base == 'G') profile.baseline_g_freq += weight;
    }

    // Codon-position-aware counting at 5' end (weighted as integer approximation)
    // For simplicity, round weight to get integer count contribution
    size_t weight_count = std::max(size_t(1), static_cast<size_t>(weight * 10 + 0.5));
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        int codon_pos = i % 3;
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.codon_pos_t_count_5prime[codon_pos] += weight_count;
        else if (base == 'C') profile.codon_pos_c_count_5prime[codon_pos] += weight_count;
    }

    // Codon-position-aware counting at 3' end (weighted)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        int codon_pos = (len - 1 - i) % 3;
        char base = fast_upper(seq[pos]);
        if (base == 'A') profile.codon_pos_a_count_3prime[codon_pos] += weight_count;
        else if (base == 'G') profile.codon_pos_g_count_3prime[codon_pos] += weight_count;
    }

    // CpG context damage tracking (weighted)
    for (size_t i = 0; i < std::min(size_t(5), len - 1); ++i) {
        char base = fast_upper(seq[i]);
        char next = fast_upper(seq[i + 1]);

        if (next == 'G') {
            if (base == 'C') {
                profile.cpg_c_count += weight_count;
            } else if (base == 'T') {
                profile.cpg_t_count += weight_count;
            }
        } else {
            if (base == 'C') {
                profile.non_cpg_c_count += weight_count;
            } else if (base == 'T') {
                profile.non_cpg_t_count += weight_count;
            }
        }
    }

    profile.n_reads++;
}

void FrameSelector::reset_sample_profile(SampleDamageProfile& profile) {
    // Reset all position-specific counts
    for (int i = 0; i < 15; ++i) {
        profile.t_freq_5prime[i] = 0.0;
        profile.c_freq_5prime[i] = 0.0;
        profile.a_freq_3prime[i] = 0.0;
        profile.g_freq_3prime[i] = 0.0;
        profile.damage_rate_5prime[i] = 0.0f;
        profile.damage_rate_3prime[i] = 0.0f;
    }

    // Reset baseline counts
    profile.baseline_t_freq = 0.0;
    profile.baseline_c_freq = 0.0;
    profile.baseline_a_freq = 0.0;
    profile.baseline_g_freq = 0.0;

    // Reset codon position counts
    for (int p = 0; p < 3; ++p) {
        profile.codon_pos_t_count_5prime[p] = 0;
        profile.codon_pos_c_count_5prime[p] = 0;
        profile.codon_pos_a_count_3prime[p] = 0;
        profile.codon_pos_g_count_3prime[p] = 0;
        profile.codon_pos_t_rate_5prime[p] = 0.5f;
        profile.codon_pos_a_rate_3prime[p] = 0.5f;
    }

    // Reset CpG counts
    profile.cpg_c_count = 0;
    profile.cpg_t_count = 0;
    profile.non_cpg_c_count = 0;
    profile.non_cpg_t_count = 0;
    profile.cpg_damage_rate = 0.0f;
    profile.non_cpg_damage_rate = 0.0f;

    // Reset summary statistics
    profile.max_damage_5prime = 0.0f;
    profile.max_damage_3prime = 0.0f;
    profile.sample_damage_prob = 0.0f;
    profile.lambda_5prime = 0.3f;
    profile.lambda_3prime = 0.3f;
    profile.library_type = SampleDamageProfile::LibraryType::UNKNOWN;
    profile.n_reads = 0;
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

    // ==========================================================================
    // BAYESIAN LIKELIHOOD-BASED DAMAGE DETECTION
    //
    // Uses sample profile's position-specific damage rates to calculate:
    //   log_LR = sum over positions of log(P(base|ancient) / P(base|modern))
    //
    // Key insight: Ancient DNA damage shows EXPONENTIAL DECAY from termini
    // Position 0 has highest damage, decaying with distance according to lambda
    // Sample profile gives us actual damage rates at each position
    // ==========================================================================

    const size_t len = seq.length();
    float log_lr = 0.0f;  // Log-likelihood ratio

    // Get baseline frequencies from sample profile (middle of reads)
    // These represent the UNDAMAGED base composition
    float baseline_t = 0.25f, baseline_c = 0.25f;
    float baseline_a = 0.25f, baseline_g = 0.25f;

    if (sample_profile.is_valid()) {
        baseline_t = static_cast<float>(sample_profile.baseline_t_freq);
        baseline_c = static_cast<float>(sample_profile.baseline_c_freq);
        baseline_a = static_cast<float>(sample_profile.baseline_a_freq);
        baseline_g = static_cast<float>(sample_profile.baseline_g_freq);

        // Ensure reasonable bounds
        baseline_t = std::max(0.15f, std::min(0.35f, baseline_t));
        baseline_c = std::max(0.15f, std::min(0.35f, baseline_c));
        baseline_a = std::max(0.15f, std::min(0.35f, baseline_a));
        baseline_g = std::max(0.15f, std::min(0.35f, baseline_g));
    }

    // ==========================================================================
    // 5' END LIKELIHOOD: C→T DAMAGE
    // Analyze first 10 positions with position-specific damage rates
    // ==========================================================================

    const size_t analyze_5prime = std::min(size_t(10), len);

    for (size_t i = 0; i < analyze_5prime; ++i) {
        char base = fast_upper(seq[i]);

        // Get damage rate at this position
        float dmg_rate = 0.0f;
        if (sample_profile.is_valid() && i < 15) {
            dmg_rate = sample_profile.damage_rate_5prime[i];
        } else {
            // Default exponential decay: ~15% at pos 0, decaying with lambda=0.3
            dmg_rate = 0.15f * std::exp(-0.3f * static_cast<float>(i));
        }

        // Codon position modulation: wobble position (3rd) more tolerant
        int codon_pos = (static_cast<int>(i) - frame + 300) % 3;
        float codon_mod = 1.0f;
        if (codon_pos == 2) {
            // Wobble position: damage more likely preserved (synonymous)
            codon_mod = 1.2f;
        } else if (codon_pos == 0) {
            // First codon position: damage less likely preserved
            codon_mod = 0.9f;
        }

        float effective_dmg = dmg_rate * codon_mod;

        if (base == 'T') {
            // T observed: could be original T or C→T damage
            // P(T|ancient) = baseline_T + baseline_C * damage_rate
            // P(T|modern) = baseline_T
            float p_t_ancient = baseline_t + baseline_c * effective_dmg;
            float p_t_modern = baseline_t;

            if (p_t_modern > 0.01f) {
                log_lr += std::log(p_t_ancient / p_t_modern);
            }
        } else if (base == 'C') {
            // C observed: C that survived (wasn't damaged)
            // P(C|ancient) = baseline_C * (1 - damage_rate)
            // P(C|modern) = baseline_C
            float p_c_ancient = baseline_c * (1.0f - effective_dmg);
            float p_c_modern = baseline_c;

            if (p_c_modern > 0.01f) {
                log_lr += std::log(p_c_ancient / p_c_modern);
            }
        }
        // A and G at 5' end are uninformative for C→T damage
    }

    // ==========================================================================
    // 3' END LIKELIHOOD: G→A DAMAGE
    // Analyze last 10 positions with position-specific damage rates
    // ==========================================================================

    const size_t analyze_3prime = std::min(size_t(10), len);

    for (size_t i = 0; i < analyze_3prime; ++i) {
        size_t pos = len - 1 - i;  // Position from end of sequence
        char base = fast_upper(seq[pos]);

        // Get damage rate at this position (distance from 3' end)
        float dmg_rate = 0.0f;
        if (sample_profile.is_valid() && i < 15) {
            dmg_rate = sample_profile.damage_rate_3prime[i];
        } else {
            // Default exponential decay
            dmg_rate = 0.12f * std::exp(-0.25f * static_cast<float>(i));
        }

        // Codon position modulation
        int codon_pos = (static_cast<int>(pos) - frame + 300) % 3;
        float codon_mod = 1.0f;
        if (codon_pos == 2) {
            codon_mod = 1.2f;  // Wobble more tolerant
        } else if (codon_pos == 0) {
            codon_mod = 0.9f;
        }

        float effective_dmg = dmg_rate * codon_mod;

        if (base == 'A') {
            // A observed: could be original A or G→A damage
            // P(A|ancient) = baseline_A + baseline_G * damage_rate
            // P(A|modern) = baseline_A
            float p_a_ancient = baseline_a + baseline_g * effective_dmg;
            float p_a_modern = baseline_a;

            if (p_a_modern > 0.01f) {
                log_lr += std::log(p_a_ancient / p_a_modern);
            }
        } else if (base == 'G') {
            // G observed: G that survived
            // P(G|ancient) = baseline_G * (1 - damage_rate)
            // P(G|modern) = baseline_G
            float p_g_ancient = baseline_g * (1.0f - effective_dmg);
            float p_g_modern = baseline_g;

            if (p_g_modern > 0.01f) {
                log_lr += std::log(p_g_ancient / p_g_modern);
            }
        }
    }

    // ==========================================================================
    // DAMAGED STOP CODON DETECTION
    // Stop codons near 5' end that could be from C→T damage are strong signals:
    //   TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
    // ==========================================================================

    for (size_t i = static_cast<size_t>(frame); i + 2 < std::min(size_t(15), len); i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        int ds = is_damaged_stop_codon(c1, c2, c3);
        if (ds > 0 && i < 12) {
            // Weight by position: closer to 5' = more likely damage
            float pos_weight = std::exp(-0.2f * static_cast<float>(i));
            log_lr += ds * 0.5f * pos_weight;
        }
    }

    // ==========================================================================
    // WOBBLE POSITION T ENRICHMENT (5' end)
    // In coding regions with C→T damage, wobble position should show more T
    // because synonymous C→T changes are preserved by selection
    // ==========================================================================

    std::array<int, 3> t_count = {0, 0, 0};
    std::array<int, 3> c_count = {0, 0, 0};

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        int codon_pos = (static_cast<int>(i) - frame + 300) % 3;
        char base = fast_upper(seq[i]);

        if (base == 'T') t_count[codon_pos]++;
        else if (base == 'C') c_count[codon_pos]++;
    }

    // Calculate T/(T+C) at wobble vs other positions
    float tc_wobble = (t_count[2] + c_count[2] > 0) ?
        static_cast<float>(t_count[2]) / (t_count[2] + c_count[2]) : 0.5f;
    float tc_other = 0.5f;
    int other_total = t_count[0] + c_count[0] + t_count[1] + c_count[1];
    if (other_total > 0) {
        tc_other = static_cast<float>(t_count[0] + t_count[1]) / other_total;
    }

    float wobble_enrichment = tc_wobble - tc_other;
    if (wobble_enrichment > 0.1f) {
        // Wobble enrichment is strong evidence for coding + damage
        log_lr += 0.6f * wobble_enrichment;
    }

    // ==========================================================================
    // SAMPLE-LEVEL PRIOR
    // Use the sample's overall damage probability as a prior
    // ==========================================================================

    float prior_log_odds = 0.0f;
    if (sample_profile.is_valid()) {
        // Convert sample_damage_prob to log-odds
        float sample_prior = std::clamp(sample_profile.sample_damage_prob, 0.1f, 0.9f);
        prior_log_odds = std::log(sample_prior / (1.0f - sample_prior));

        // Weight the prior based on sample size (more reads = more confidence)
        float prior_weight = std::min(1.0f, static_cast<float>(sample_profile.n_reads) / 10000.0f);
        prior_log_odds *= prior_weight * 0.5f;  // Dampen to not overwhelm read evidence
    }

    // Combine likelihood and prior
    float total_log_odds = log_lr + prior_log_odds;

    // Convert log-odds to probability
    float prob = 1.0f / (1.0f + std::exp(-total_log_odds));

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

    // Use fixed baseline of 0.5 for damage likelihood calculation
    // This creates better contrast for per-read classification
    // The sample-specific rates are used for the damage at each position
    const float baseline_T = 0.50f;
    const float baseline_A = 0.50f;

    // Default damage rates (used if sample profile not available)
    const float default_damage_5prime = 0.15f;
    const float default_damage_3prime = 0.12f;
    const float default_lambda_5 = 0.3f;
    const float default_lambda_3 = 0.25f;

    // Phred to error probability
    auto phred_to_error = [](char q) -> float {
        int phred = static_cast<int>(q) - 33;
        phred = std::clamp(phred, 0, 40);
        return std::pow(10.0f, -phred / 10.0f);
    };

    // =========================================================================
    // 5' END: FEATURE-BASED DAMAGE DETECTION
    // Use fixed bonuses for specific damage signatures (more discriminative than
    // Bayesian likelihood ratios which have tiny contrast in AT-rich samples)
    // =========================================================================

    // Terminal base pattern (position 0)
    if (fast_upper(seq[0]) == 'T') {
        log_lr += 0.7f;  // T at position 0 - primary damage signal
    } else if (fast_upper(seq[0]) == 'C') {
        log_lr -= 0.2f;  // C at position 0 - evidence against damage
    }

    // Position 1 pattern
    if (len > 1 && fast_upper(seq[1]) == 'T') {
        log_lr += 0.3f;  // T at position 1 - secondary damage signal
    }

    // Dinucleotide patterns at 5' end (strong damage indicators)
    if (len >= 2) {
        char b0 = fast_upper(seq[0]);
        char b1 = fast_upper(seq[1]);
        if (b0 == 'T' && b1 == 'G') log_lr += 0.5f;  // TG from CG (CpG damage!)
        if (b0 == 'T' && b1 == 'T') log_lr += 0.2f;  // TT from CT
        if (b0 == 'T' && b1 == 'A') log_lr += 0.15f; // TA from CA
    }

    // Quality-weighted position scoring for first 5 positions
    float decay_weights[] = {1.0f, 0.6f, 0.4f, 0.25f, 0.15f};
    for (size_t i = 2; i < std::min(size_t(5), len); ++i) {
        char base = fast_upper(seq[i]);
        float weight = decay_weights[i];

        float quality_mod = 1.0f;
        if (has_quality) {
            float error_prob = phred_to_error(quality[i]);
            if (error_prob > 0.1f) quality_mod = 0.5f;
        }

        if (base == 'T') {
            log_lr += 0.15f * weight * quality_mod;
        } else if (base == 'C') {
            log_lr -= 0.08f * weight * quality_mod;
        }
    }

    // =========================================================================
    // 3' END: FEATURE-BASED DAMAGE DETECTION
    // =========================================================================

    // Terminal base pattern (last position)
    if (fast_upper(seq[len - 1]) == 'A') {
        log_lr += 0.6f;  // A at last position - primary 3' damage signal
    } else if (fast_upper(seq[len - 1]) == 'G') {
        log_lr -= 0.15f;  // G at last position - evidence against damage
    }

    // Position len-2 pattern
    if (len > 1 && fast_upper(seq[len - 2]) == 'A') {
        log_lr += 0.25f;
    }

    // Dinucleotide patterns at 3' end
    if (len >= 2) {
        char bm1 = fast_upper(seq[len - 1]);
        char bm2 = fast_upper(seq[len - 2]);
        if (bm1 == 'A' && bm2 == 'C') log_lr += 0.4f;  // CA from CG complement
        if (bm1 == 'A' && bm2 == 'A') log_lr += 0.2f;  // AA from GA
        if (bm1 == 'A' && bm2 == 'T') log_lr += 0.15f; // TA from TG complement
    }

    // Quality-weighted position scoring for last 5 positions
    for (size_t i = 2; i < std::min(size_t(5), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        float weight = decay_weights[i];

        float quality_mod = 1.0f;
        if (has_quality) {
            float error_prob = phred_to_error(quality[pos]);
            if (error_prob > 0.1f) quality_mod = 0.5f;
        }

        if (base == 'A') {
            log_lr += 0.12f * weight * quality_mod;
        } else if (base == 'G') {
            log_lr -= 0.06f * weight * quality_mod;
        }
    }

    // =========================================================================
    // HEXAMER DAMAGE LOOKUP (5' end)
    // Use GTDB-derived hexamer damage probabilities
    // =========================================================================
    if (len >= 6) {
        // Check first few hexamers at 5' end
        for (size_t i = 0; i < std::min(size_t(3), len - 5); ++i) {
            float hex_prob = get_hexamer_damage_prob(seq.c_str() + i);
            if (hex_prob > 0.9f) {
                // High-confidence damage hexamer (e.g., TGA from CGA)
                float weight = 1.0f / (1 << i);  // 1.0, 0.5, 0.25
                log_lr += std::log(hex_prob / 0.1f) * weight * 0.3f;
            }
        }
    }

    // =========================================================================
    // DAMAGED STOP CODON DETECTION
    // TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
    // These are strong signals when found near 5' end
    // =========================================================================
    for (size_t i = frame; i + 2 < std::min(size_t(15), len); i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        // Check for damaged stop codons
        bool is_damaged_stop = false;
        if (c1 == 'T' && c2 == 'A' && c3 == 'A') is_damaged_stop = true;  // TAA from CAA
        if (c1 == 'T' && c2 == 'A' && c3 == 'G') is_damaged_stop = true;  // TAG from CAG
        if (c1 == 'T' && c2 == 'G' && c3 == 'A') is_damaged_stop = true;  // TGA from CGA

        if (is_damaged_stop && i < 12) {
            // Weight by position (closer to 5' = more likely damage)
            float pos_weight = std::exp(-0.3f * static_cast<float>(i));
            log_lr += 0.8f * pos_weight;  // Strong signal
        }
    }

    // =========================================================================
    // WOBBLE POSITION T ENRICHMENT (5' end)
    // In coding regions with C→T damage, wobble position (pos 3) should
    // show more T because synonymous changes are better tolerated
    // =========================================================================
    std::array<int, 3> t_count = {0, 0, 0};
    std::array<int, 3> c_count = {0, 0, 0};

    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        int codon_pos = (static_cast<int>(i) - frame + 300) % 3;
        char base = fast_upper(seq[i]);
        if (base == 'T') t_count[codon_pos]++;
        else if (base == 'C') c_count[codon_pos]++;
    }

    // Calculate T/(T+C) ratio at wobble position vs others
    float tc_wobble = (t_count[2] + c_count[2] > 0) ?
        static_cast<float>(t_count[2]) / (t_count[2] + c_count[2]) : 0.5f;
    float tc_other = 0.5f;
    int other_total = t_count[0] + c_count[0] + t_count[1] + c_count[1];
    if (other_total > 0) {
        tc_other = static_cast<float>(t_count[0] + t_count[1]) / other_total;
    }

    float wobble_enrichment = tc_wobble - tc_other;
    if (wobble_enrichment > 0.1f) {
        log_lr += 0.5f * wobble_enrichment;  // Coding + damage signal
    }

    // =========================================================================
    // JOINT 5' + 3' PATTERN (both ends damaged is strong evidence)
    // =========================================================================
    if (fast_upper(seq[0]) == 'T' && fast_upper(seq[len - 1]) == 'A') {
        log_lr += 0.5f;  // Bonus for joint pattern
    }

    // =========================================================================
    // WITHIN-READ GRADIENT: terminal vs interior contrast
    // Damaged reads should have higher T/(T+C) at terminals than interior
    // This helps discriminate in AT-rich samples
    // =========================================================================
    if (len >= 20) {
        // 5' gradient: compare positions 0-4 to positions 10-14
        int t_term5 = 0, c_term5 = 0;
        int t_int5 = 0, c_int5 = 0;
        for (int i = 0; i < 5; ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'T') t_term5++;
            else if (b == 'C') c_term5++;
        }
        for (int i = 10; i < 15 && i < static_cast<int>(len); ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'T') t_int5++;
            else if (b == 'C') c_int5++;
        }

        // Compare terminal to interior T/(T+C)
        if (t_term5 + c_term5 >= 2 && t_int5 + c_int5 >= 2) {
            float tc_term = static_cast<float>(t_term5) / (t_term5 + c_term5);
            float tc_int = static_cast<float>(t_int5) / (t_int5 + c_int5);
            float gradient_5 = tc_term - tc_int;  // Positive = more T at terminal
            // Lower threshold and higher weight for better sensitivity
            if (gradient_5 > 0.10f) {
                log_lr += 1.0f * gradient_5;  // Strong gradient signal
            } else if (gradient_5 < -0.10f) {
                log_lr += 0.5f * gradient_5;  // Less T at terminal = evidence against damage
            }
        }

        // 3' gradient: compare positions len-5 to len-1 vs len-15 to len-11
        int a_term3 = 0, g_term3 = 0;
        int a_int3 = 0, g_int3 = 0;
        for (size_t i = len - 5; i < len; ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'A') a_term3++;
            else if (b == 'G') g_term3++;
        }
        for (size_t i = len - 15; i < len - 10; ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'A') a_int3++;
            else if (b == 'G') g_int3++;
        }

        if (a_term3 + g_term3 >= 2 && a_int3 + g_int3 >= 2) {
            float ag_term = static_cast<float>(a_term3) / (a_term3 + g_term3);
            float ag_int = static_cast<float>(a_int3) / (a_int3 + g_int3);
            float gradient_3 = ag_term - ag_int;
            // Lower threshold and higher weight for better sensitivity
            if (gradient_3 > 0.10f) {
                log_lr += 0.8f * gradient_3;
            } else if (gradient_3 < -0.10f) {
                log_lr += 0.4f * gradient_3;
            }
        }
    }

    // =========================================================================
    // SAMPLE-LEVEL CALIBRATION
    // Don't add a sample-wide prior as it reduces per-read discrimination
    // Instead, use sample info to scale feature weights based on expected damage level
    // =========================================================================
    if (sample_profile.is_valid()) {
        // In high-damage samples, we expect more damage patterns, so they're less surprising
        // Scale down the bonus slightly in high-damage samples to maintain discrimination
        float damage_level = (sample_profile.max_damage_5prime + sample_profile.max_damage_3prime) / 2.0f;
        if (damage_level > 0.15f) {
            // High damage sample: the signal is expected, so dampen slightly
            log_lr *= 0.9f;
        } else if (damage_level < 0.05f) {
            // Low damage sample: damage signal is more surprising when present
            log_lr *= 1.1f;
        }
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

// ============================================================================
// DAMAGE CODON SCORE FOR ITERATIVE ENRICHMENT
// ============================================================================

float FrameSelector::compute_damage_codon_score(const std::string& seq) {
    if (seq.length() < 15) return 0.0f;

    const size_t len = seq.length();
    float score = 0.0f;
    float max_possible = 0.0f;

    // =========================================================================
    // FEATURE 1: Terminal T enrichment at 5' end (positions 0-4)
    // C→T damage creates T's at terminal positions
    // Weight by position (closer to 5' = more indicative of damage)
    // =========================================================================
    const float pos_weights[] = {1.0f, 0.7f, 0.5f, 0.35f, 0.25f};

    for (size_t i = 0; i < std::min(size_t(5), len); ++i) {
        char base = fast_upper(seq[i]);
        float weight = pos_weights[i];
        max_possible += weight;

        if (base == 'T') {
            score += weight;  // T at 5' = damage signal
        } else if (base == 'C') {
            score -= weight * 0.3f;  // C at 5' = evidence against damage
        }
    }

    // =========================================================================
    // FEATURE 2: Terminal A enrichment at 3' end (positions len-1 to len-5)
    // G→A damage creates A's at terminal positions
    // =========================================================================
    for (size_t i = 0; i < std::min(size_t(5), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        float weight = pos_weights[i] * 0.8f;  // 3' signal slightly weaker
        max_possible += weight;

        if (base == 'A') {
            score += weight;  // A at 3' = damage signal
        } else if (base == 'G') {
            score -= weight * 0.3f;  // G at 3' = evidence against damage
        }
    }

    // =========================================================================
    // FEATURE 3: Damage-susceptible dinucleotide patterns at 5' end
    // CpG context shows elevated damage - TG from CG is strong signal
    // =========================================================================
    for (size_t i = 0; i < std::min(size_t(4), len - 1); ++i) {
        char b1 = fast_upper(seq[i]);
        char b2 = fast_upper(seq[i + 1]);
        float pos_mod = 1.0f - (i * 0.2f);  // Position decay

        if (b1 == 'T' && b2 == 'G') {
            // TG from CpG damage - very strong signal
            score += 0.8f * pos_mod;
            max_possible += 0.8f * pos_mod;
        } else if (b1 == 'T' && b2 == 'A') {
            // TA from CA - moderate signal
            score += 0.3f * pos_mod;
            max_possible += 0.3f * pos_mod;
        } else if (b1 == 'C' && b2 == 'G') {
            // Intact CpG at 5' - evidence against damage
            score -= 0.4f * pos_mod;
            max_possible += 0.4f * pos_mod;
        }
    }

    // =========================================================================
    // FEATURE 4: Damage-induced stop codons near 5' end
    // TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
    // These are VERY strong damage indicators in coding sequences
    // =========================================================================
    // Check all 3 frames for damage stops in first 15bp
    for (int frame = 0; frame < 3; ++frame) {
        for (size_t i = frame; i + 2 < std::min(size_t(15), len); i += 3) {
            char c1 = fast_upper(seq[i]);
            char c2 = fast_upper(seq[i + 1]);
            char c3 = fast_upper(seq[i + 2]);

            float pos_mod = std::exp(-0.2f * static_cast<float>(i));

            // Check for damaged stop codons
            bool is_damaged_stop = false;
            float stop_weight = 0.0f;

            if (c1 == 'T' && c2 == 'A' && c3 == 'A') {
                // TAA from CAA (Gln) - very common damage pattern
                is_damaged_stop = true;
                stop_weight = 1.5f;
            } else if (c1 == 'T' && c2 == 'A' && c3 == 'G') {
                // TAG from CAG (Gln)
                is_damaged_stop = true;
                stop_weight = 1.5f;
            } else if (c1 == 'T' && c2 == 'G' && c3 == 'A') {
                // TGA from CGA (Arg)
                is_damaged_stop = true;
                stop_weight = 1.2f;
            }

            if (is_damaged_stop) {
                score += stop_weight * pos_mod;
                max_possible += stop_weight * pos_mod;
            }

            // Also check for TGG from CGG (Arg) - not stop but damage pattern
            if (c1 == 'T' && c2 == 'G' && c3 == 'G') {
                score += 0.5f * pos_mod;
                max_possible += 0.5f * pos_mod;
            }
        }
    }

    // =========================================================================
    // FEATURE 5: T-rich codons at 5' end consistent with C→T damage
    // TNT, TTN patterns from CNC, CTN codons
    // =========================================================================
    for (int frame = 0; frame < 3; ++frame) {
        for (size_t i = frame; i + 2 < std::min(size_t(12), len); i += 3) {
            char c1 = fast_upper(seq[i]);
            char c2 = fast_upper(seq[i + 1]);
            char c3 = fast_upper(seq[i + 2]);

            float pos_mod = std::exp(-0.25f * static_cast<float>(i));
            int t_count = (c1 == 'T') + (c2 == 'T') + (c3 == 'T');
            int c_count = (c1 == 'C') + (c2 == 'C') + (c3 == 'C');

            // Multiple T's in a codon at 5' = potential damage
            if (t_count >= 2 && c_count == 0) {
                score += 0.3f * pos_mod;
                max_possible += 0.3f * pos_mod;
            }
            // Multiple C's preserved at 5' = unlikely to be damaged
            else if (c_count >= 2 && t_count == 0) {
                score -= 0.2f * pos_mod;
                max_possible += 0.2f * pos_mod;
            }
        }
    }

    // =========================================================================
    // FEATURE 6: Within-read T gradient (terminal vs interior)
    // Damaged reads should have higher T/(T+C) at terminals than interior
    // This is the key discriminator for AT-rich samples
    // =========================================================================
    if (len >= 20) {
        // 5' gradient: compare positions 0-4 to positions 10-14
        int t_term = 0, c_term = 0;
        int t_int = 0, c_int = 0;

        for (int i = 0; i < 5; ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'T') t_term++;
            else if (b == 'C') c_term++;
        }
        for (int i = 10; i < 15 && i < static_cast<int>(len); ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'T') t_int++;
            else if (b == 'C') c_int++;
        }

        // Calculate T/(T+C) gradient
        if (t_term + c_term >= 2 && t_int + c_int >= 2) {
            float tc_term = static_cast<float>(t_term) / (t_term + c_term);
            float tc_int = static_cast<float>(t_int) / (t_int + c_int);
            float gradient = tc_term - tc_int;

            // Positive gradient = more T at terminal = damage signal
            if (gradient > 0.15f) {
                score += gradient * 2.0f;  // Strong weight for gradient
                max_possible += 0.5f;  // Cap contribution
            } else if (gradient < -0.15f) {
                // Negative gradient = MORE C at terminal than interior
                // This is anti-damage signal
                score += gradient * 1.0f;  // Negative contribution
                max_possible += 0.3f;
            }
        }

        // 3' gradient: compare positions len-5 to len-1 vs len-15 to len-11
        int a_term3 = 0, g_term3 = 0;
        int a_int3 = 0, g_int3 = 0;

        for (size_t i = len - 5; i < len; ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'A') a_term3++;
            else if (b == 'G') g_term3++;
        }
        for (size_t i = len - 15; i < len - 10 && i < len; ++i) {
            char b = fast_upper(seq[i]);
            if (b == 'A') a_int3++;
            else if (b == 'G') g_int3++;
        }

        if (a_term3 + g_term3 >= 2 && a_int3 + g_int3 >= 2) {
            float ag_term = static_cast<float>(a_term3) / (a_term3 + g_term3);
            float ag_int = static_cast<float>(a_int3) / (a_int3 + g_int3);
            float gradient = ag_term - ag_int;

            if (gradient > 0.15f) {
                score += gradient * 1.5f;
                max_possible += 0.4f;
            } else if (gradient < -0.15f) {
                score += gradient * 0.8f;
                max_possible += 0.25f;
            }
        }
    }

    // =========================================================================
    // FEATURE 7: Joint 5' + 3' terminal pattern
    // T at 5' AND A at 3' together is stronger evidence than either alone
    // =========================================================================
    if (fast_upper(seq[0]) == 'T' && fast_upper(seq[len - 1]) == 'A') {
        score += 0.6f;  // Bonus for joint pattern
        max_possible += 0.6f;
    }

    // Normalize score to 0.0-1.0 range
    if (max_possible > 0.0f) {
        float normalized = (score + max_possible) / (2.0f * max_possible);
        return std::clamp(normalized, 0.0f, 1.0f);
    }

    return 0.5f;  // Default neutral score
}

} // namespace agp
