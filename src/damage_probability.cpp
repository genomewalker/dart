// Per-read damage probability estimation for ancient DNA

#include "agp/frame_selector.hpp"
#include <algorithm>
#include <cmath>

namespace agp {

static inline char fast_upper(char c) {
    return (c >= 'a' && c <= 'z') ? c - 32 : c;
}

float FrameSelector::estimate_damage_signal(const std::string& seq) {
    if (seq.length() < 10) return 0.5f;  // Too short to estimate

    size_t len = seq.length();
    float log_odds = 0.0f;

    // PRIMARY SIGNAL: Terminal base patterns
    // T at 5' end (C→T damage) - strongest signal
    if (fast_upper(seq[0]) == 'T') {
        log_odds += 0.8f;  // Strong positive evidence
    } else if (fast_upper(seq[0]) == 'C') {
        log_odds -= 0.3f;  // Moderate negative evidence (survived without damage)
    }

    // A at 3' end (G→A damage) - also strong
    if (fast_upper(seq[len - 1]) == 'A') {
        log_odds += 0.7f;
    } else if (fast_upper(seq[len - 1]) == 'G') {
        log_odds -= 0.25f;
    }

    // SECONDARY SIGNAL: Position 1 (second base)
    if (len > 1) {
        if (fast_upper(seq[1]) == 'T') log_odds += 0.35f;
        if (fast_upper(seq[len - 2]) == 'A') log_odds += 0.3f;
    }

    // TERTIARY SIGNAL: Position 2
    if (len > 2) {
        if (fast_upper(seq[2]) == 'T') log_odds += 0.15f;
        if (fast_upper(seq[len - 3]) == 'A') log_odds += 0.12f;
    }

    // JOINT SIGNAL: Both ends damaged is stronger evidence
    if (fast_upper(seq[0]) == 'T' && fast_upper(seq[len - 1]) == 'A') {
        log_odds += 0.4f;  // Bonus for joint pattern
    }

    // NEGATIVE SIGNAL: Multiple C's at 5' or G's at 3'
    int c_count_5 = 0, g_count_3 = 0;
    for (size_t i = 0; i < std::min(size_t(3), len); ++i) {
        if (fast_upper(seq[i]) == 'C') c_count_5++;
        if (fast_upper(seq[len - 1 - i]) == 'G') g_count_3++;
    }
    if (c_count_5 >= 2) log_odds -= 0.3f;
    if (g_count_3 >= 2) log_odds -= 0.25f;

    // Convert log-odds to probability
    float prob = 1.0f / (1.0f + std::exp(-log_odds));

    return std::clamp(prob, 0.05f, 0.95f);
}

float FrameSelector::compute_damage_percentage(
    const std::string& seq,
    const SampleDamageProfile& sample_profile) {

    if (seq.length() < 15 || !sample_profile.is_valid()) return 0.0f;

    // =========================================================================
    // CHANNEL B SAMPLE-LEVEL CHECK
    // If sample has no real damage, return 0 for all reads
    // =========================================================================
    if (sample_profile.damage_artifact) {
        return 0.0f;  // Sample is compositional artifact
    }
    if (!sample_profile.damage_validated && sample_profile.channel_b_valid &&
        sample_profile.stop_decay_llr_5prime < -100.0f) {
        return 0.0f;  // Channel B strongly contradicts damage
    }

    const size_t len = seq.length();
    const size_t analyze_len = std::min<size_t>(15, len);  // extend window to match sample profile

    // Calculate log-likelihood ratio using sample's damage rates
    float log_lr = 0.0f;

    // Get baseline frequencies from sample
    float baseline_tc = static_cast<float>(sample_profile.baseline_t_freq /
        (sample_profile.baseline_t_freq + sample_profile.baseline_c_freq + 0.001));
    float baseline_ag = static_cast<float>(sample_profile.baseline_a_freq /
        (sample_profile.baseline_a_freq + sample_profile.baseline_g_freq + 0.001));

    // Light GC-normalization to avoid extreme baselines suppressing signal
    baseline_tc = std::clamp(baseline_tc, 0.05f, 0.95f);
    baseline_ag = std::clamp(baseline_ag, 0.05f, 0.95f);

    const float baseline_c = 1.0f - baseline_tc;
    const float baseline_g = 1.0f - baseline_ag;

    // 5' end likelihood contribution
    for (size_t i = 0; i < analyze_len; ++i) {
        char base = fast_upper(seq[i]);
        float dmg_rate = (i < 15) ? sample_profile.damage_rate_5prime[i] : 0.0f;

        if (base == 'T') {
            // T observed: could be original or C→T
            // P(T|ancient) = baseline_T + baseline_C * damage_rate
            // P(T|modern) = baseline_T
            float p_ancient = baseline_tc + (1.0f - baseline_tc) * dmg_rate;
            float p_modern = baseline_tc;
            if (p_modern > 0.01f) {
                log_lr += std::log(p_ancient / p_modern);
            }
        } else if (base == 'C') {
            // C observed: survived without damage
            float p_ancient = (1.0f - baseline_tc) * (1.0f - dmg_rate);
            float p_modern = (1.0f - baseline_tc);
            if (p_modern > 0.01f) {
                log_lr += std::log(p_ancient / p_modern);
            }
        }
    }

    // 3' end likelihood contribution
    for (size_t i = 0; i < analyze_len; ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        float dmg_rate = (i < 15) ? sample_profile.damage_rate_3prime[i] : 0.0f;

        if (base == 'A') {
            float p_ancient = baseline_ag + (1.0f - baseline_ag) * dmg_rate;
            float p_modern = baseline_ag;
            if (p_modern > 0.01f) {
                log_lr += std::log(p_ancient / p_modern);
            }
        } else if (base == 'G') {
            float p_ancient = (1.0f - baseline_ag) * (1.0f - dmg_rate);
            float p_modern = (1.0f - baseline_ag);
            if (p_modern > 0.01f) {
                log_lr += std::log(p_ancient / p_modern);
            }
        }
    }

    // Calibrate against the maximum expected log-LR for this sample's damage rates.
    // This maps 0 = no evidence, 100 = per-read signal matching sample D_max.
    float max_expected_lr = 0.0f;
    for (size_t i = 0; i < analyze_len; ++i) {
        float dmg_5 = (i < 15) ? sample_profile.damage_rate_5prime[i] : 0.0f;
        float dmg_3 = (i < 15) ? sample_profile.damage_rate_3prime[i] : 0.0f;
        if (dmg_5 > 0.0f && baseline_tc > 0.001f && baseline_c > 0.001f) {
            max_expected_lr += std::log(1.0f + dmg_5 * baseline_c / baseline_tc);
        }
        if (dmg_3 > 0.0f && baseline_ag > 0.001f && baseline_g > 0.001f) {
            max_expected_lr += std::log(1.0f + dmg_3 * baseline_g / baseline_ag);
        }
    }

    float percentage = 0.0f;
    if (max_expected_lr > 1e-3f) {
        percentage = 100.0f * (log_lr / max_expected_lr);
    }

    // Negative evidence should not produce negative percentages
    percentage = std::clamp(percentage, 0.0f, 100.0f);

    return std::clamp(percentage, 0.0f, 100.0f);
}

float FrameSelector::estimate_damage_signal_with_codons(
    const std::string& seq,
    int frame,
    bool forward) {

    if (seq.length() < 15) return 0.5f;

    size_t len = seq.length();
    float log_odds = 0.0f;

    // Basic terminal signals (as in estimate_damage_signal)
    if (fast_upper(seq[0]) == 'T') log_odds += 0.7f;
    if (fast_upper(seq[len - 1]) == 'A') log_odds += 0.6f;
    if (len > 1 && fast_upper(seq[1]) == 'T') log_odds += 0.3f;
    if (len > 1 && fast_upper(seq[len - 2]) == 'A') log_odds += 0.25f;

    // Negative signals
    if (fast_upper(seq[0]) == 'C') log_odds -= 0.2f;
    if (fast_upper(seq[len - 1]) == 'G') log_odds -= 0.15f;

    // CODON-SPECIFIC SIGNALS

    // Damaged stop codon detection at 5' end
    // TAA from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
    for (size_t i = frame; i + 2 < std::min(size_t(15), len); i += 3) {
        char c1 = fast_upper(seq[i]);
        char c2 = fast_upper(seq[i + 1]);
        char c3 = fast_upper(seq[i + 2]);

        if (c1 == 'T' && c2 == 'A' && c3 == 'A') {
            float pos_weight = std::exp(-0.2f * static_cast<float>(i));
            log_odds += 0.6f * pos_weight;  // TAA from CAA
        } else if (c1 == 'T' && c2 == 'A' && c3 == 'G') {
            float pos_weight = std::exp(-0.2f * static_cast<float>(i));
            log_odds += 0.6f * pos_weight;  // TAG from CAG
        } else if (c1 == 'T' && c2 == 'G' && c3 == 'A') {
            float pos_weight = std::exp(-0.2f * static_cast<float>(i));
            log_odds += 0.5f * pos_weight;  // TGA from CGA
        }
    }

    // T at wobble position (codon pos 3) - often synonymous
    int t_at_wobble = 0;
    int t_at_other = 0;
    for (size_t i = 0; i < std::min(size_t(12), len); ++i) {
        int codon_pos = (static_cast<int>(i) - frame + 300) % 3;
        if (fast_upper(seq[i]) == 'T') {
            if (codon_pos == 2) t_at_wobble++;
            else t_at_other++;
        }
    }

    // Wobble enrichment is a positive signal for damage
    if (t_at_wobble > t_at_other && t_at_wobble >= 2) {
        log_odds += 0.3f;
    }

    // Joint terminal pattern bonus
    if (fast_upper(seq[0]) == 'T' && fast_upper(seq[len - 1]) == 'A') {
        log_odds += 0.3f;
    }

    float prob = 1.0f / (1.0f + std::exp(-log_odds));
    return std::clamp(prob, 0.05f, 0.95f);
}

float FrameSelector::estimate_damage_signal_with_sample_profile(
    const std::string& seq,
    const SampleDamageProfile& sample_profile) {

    if (seq.length() < 15) return 0.5f;
    if (!sample_profile.is_valid()) {
        return estimate_damage_signal_with_codons(seq, 0, true);
    }

    // =========================================================================
    // CHANNEL B SAMPLE-LEVEL PRIOR
    // If the sample was identified as having no real damage (artifact or no signal),
    // individual reads shouldn't get high damage probabilities even if they happen
    // to have T at terminal positions - those T's are compositional, not damage.
    // =========================================================================
    float sample_prior = 0.0f;  // Log-odds prior from sample-level evidence

    if (sample_profile.damage_artifact) {
        // Sample identified as compositional artifact by Channel B
        // Strong negative prior - reads are very unlikely to be damaged
        sample_prior = -2.0f;
    } else if (!sample_profile.damage_validated && sample_profile.channel_b_valid &&
               sample_profile.stop_decay_llr_5prime < -100.0f) {
        // Channel B strongly contradicts damage (no stop excess)
        // Moderate negative prior
        sample_prior = -1.5f;
    } else if (sample_profile.damage_validated) {
        // Both channels agree on real damage
        // Positive prior scaled by damage level
        float damage_level = std::clamp(sample_profile.d_max_combined, 0.0f, 0.5f);
        sample_prior = damage_level * 2.0f;  // 0 to 1.0 log-odds boost
    }
    // Else: uncertain, use neutral prior (0.0)

    size_t len = seq.length();
    float log_odds = sample_prior;  // Start with sample-level prior

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

} // namespace agp
