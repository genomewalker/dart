// Damage scoring and frame consistency for ancient DNA

#include "agp/frame_selector.hpp"
#include "agp/damage_likelihood.hpp"
#include <algorithm>
#include <cmath>
#include <array>

namespace agp {

static inline char fast_upper(char c) {
    return (c >= 'a' && c <= 'z') ? c - 32 : c;
}

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
    // CHANNEL B SAMPLE-LEVEL PRIOR
    // If the sample was identified as having no real damage (artifact or no signal),
    // individual reads shouldn't get high damage probabilities.
    // ==========================================================================
    if (sample_profile.damage_artifact) {
        // Sample identified as compositional artifact by Channel B
        // Return low probability - reads are very unlikely to be damaged
        return 0.15f;
    }
    if (!sample_profile.damage_validated && sample_profile.channel_b_valid &&
        sample_profile.stop_decay_llr_5prime < -100.0f) {
        // Channel B strongly contradicts damage (no stop excess)
        return 0.20f;
    }

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

    // ==========================================================================
    // CHANNEL B SAMPLE-LEVEL PRIOR
    // If the sample was identified as having no real damage, return low probability
    // ==========================================================================
    if (sample_profile.damage_artifact) {
        return 0.15f;  // Sample is compositional artifact
    }
    if (!sample_profile.damage_validated && sample_profile.channel_b_valid &&
        sample_profile.stop_decay_llr_5prime < -100.0f) {
        return 0.20f;  // Channel B strongly contradicts damage
    }

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
            int t_cnt = (c1 == 'T') + (c2 == 'T') + (c3 == 'T');
            int c_cnt = (c1 == 'C') + (c2 == 'C') + (c3 == 'C');

            // Multiple T's in a codon at 5' = potential damage
            if (t_cnt >= 2 && c_cnt == 0) {
                score += 0.3f * pos_mod;
                max_possible += 0.3f * pos_mod;
            }
            // Multiple C's preserved at 5' = unlikely to be damaged
            else if (c_cnt >= 2 && t_cnt == 0) {
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
