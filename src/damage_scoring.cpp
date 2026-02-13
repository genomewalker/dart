// Damage scoring and frame consistency for ancient DNA

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include "agp/damage_likelihood.hpp"
#include "agp/hexamer_tables.hpp"
#include <algorithm>
#include <cmath>
#include <array>

namespace agp {

// ============================================================================
// ADAPTIVE CALIBRATION FOR TERMINAL DAMAGE PROBABILITY
//
// Key insight: Raw scores cluster around ~0.65 regardless of actual damage rate.
// We need to recalibrate so that:
//   - Mean output ≈ sample's validated d_max
//   - Scores spread appropriately around that mean
//
// Method: Sigmoid scaling based on sample d_max
//   calibrated = d_max * sigmoid(k * (raw - threshold))
// where threshold is chosen so mean output ≈ d_max
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

float FrameSelector::estimate_damage_signal(
    const std::string& seq,
    int frame,
    const SampleDamageProfile& sample_profile,
    const std::string& quality) {

    if (seq.length() < 15) return 0.5f;

    const bool has_quality = !quality.empty() && quality.length() == seq.length();

    // ==========================================================================
    // TRI-STATE DAMAGE VALIDATION
    // Uses shared get_damage_suppression_factor() for consistent semantics:
    //   - VALIDATED: factor = 1.0 (full confidence in per-read scores)
    //   - CONTRADICTED: factor = 0.0 (hard suppress for artifacts)
    //   - UNVALIDATED: factor = 0.5 (soft suppression)
    // ==========================================================================
    const float channel_b_factor = get_damage_suppression_factor(sample_profile);

    // Early exit for contradicted samples (artifact identified)
    if (channel_b_factor == 0.0f) {
        return 0.0f;  // Hard suppress for artifacts
    }

    // ==========================================================================
    // GC-CONDITIONAL BAYESIAN LIKELIHOOD-BASED DAMAGE DETECTION
    //
    // Uses GC-bin specific damage rates and priors for per-read inference:
    //   1. Get GC-conditional parameters for this read
    //   2. Compute per-read LLR using bin-specific delta_s
    //   3. Combine with GC-bin prior: posterior = sigmoid(logit(p0) + LLR)
    //
    // Key insight: Ancient DNA damage shows EXPONENTIAL DECAY from termini
    // Different GC bins may have different damage levels due to organism heterogeneity
    // ==========================================================================

    const size_t len = seq.length();
    float log_lr = 0.0f;  // Log-likelihood ratio

    // Get GC-conditional parameters for this read
    auto gc_params = sample_profile.get_gc_params(seq);

    // Use GC-bin prior if available (p_damaged), else neutral prior
    float p0 = gc_params.bin_valid ? gc_params.p_damaged : 0.5f;

    // Get damage rate: use GC-bin specific delta_s if available
    float delta_s = gc_params.bin_valid ? gc_params.delta_s : sample_profile.d_max_combined;
    float lambda = sample_profile.lambda_5prime;

    // Use GC-bin specific baseline if available
    float baseline_tc = gc_params.bin_valid ? gc_params.baseline_tc : 0.5f;
    if (!gc_params.bin_valid && sample_profile.is_valid()) {
        baseline_tc = static_cast<float>(sample_profile.baseline_t_freq /
            (sample_profile.baseline_t_freq + sample_profile.baseline_c_freq + 0.001));
    }
    baseline_tc = std::clamp(baseline_tc, 0.15f, 0.85f);

    // Get baseline frequencies for 3' end (global)
    float baseline_ag = 0.5f;
    if (sample_profile.is_valid()) {
        baseline_ag = static_cast<float>(sample_profile.baseline_a_freq /
            (sample_profile.baseline_a_freq + sample_profile.baseline_g_freq + 0.001));
    }
    baseline_ag = std::clamp(baseline_ag, 0.15f, 0.85f);

    // For undamaged bins, use sample-level damage rate as delta_s fallback
    // Don't shortcut - still compute per-read LLR to detect individual damaged reads
    if (gc_params.bin_valid && gc_params.delta_s < 0.005f) {
        // Use global damage rate instead of bin-specific (which is 0)
        delta_s = sample_profile.d_max_combined > 0.01f ?
                  static_cast<float>(sample_profile.d_max_combined) : 0.05f;
    }

    // ==========================================================================
    // 5' END LIKELIHOOD: C→T DAMAGE
    // Analyze first 10 positions with GC-bin specific damage rates
    // ==========================================================================

    const size_t analyze_5prime = std::min(size_t(10), len);

    for (size_t i = 0; i < analyze_5prime; ++i) {
        char base = fast_upper(seq[i]);

        // Quality-based confidence weight (1.0 for FASTA or high quality)
        float quality_weight = 1.0f;
        if (has_quality) {
            float error_prob = phred_to_error_prob(quality[i]);
            // Low quality at damage-prone positions: reduce weight
            // High error rate means less confident in the base call
            quality_weight = 1.0f - std::min(error_prob * 2.0f, 0.5f);
        }

        // Compute damage rate at this position using GC-bin delta_s
        float decay = std::exp(-lambda * static_cast<float>(i));
        float dmg_rate = delta_s * decay;

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

        float effective_dmg = std::min(dmg_rate * codon_mod, 0.95f);

        if (base == 'T') {
            // T observed: could be original T or C→T damage
            // P(T|ancient) = baseline_T + baseline_C * damage_rate
            // P(T|modern) = baseline_T
            float p_t_ancient = baseline_tc + (1.0f - baseline_tc) * effective_dmg;
            float p_t_modern = baseline_tc;

            if (p_t_modern > 0.01f) {
                log_lr += std::log(p_t_ancient / p_t_modern) * quality_weight;
            }
        } else if (base == 'C') {
            // C observed: C that survived (wasn't damaged)
            // P(C|ancient) = baseline_C * (1 - damage_rate)
            // P(C|modern) = baseline_C
            float p_c_ancient = (1.0f - baseline_tc) * (1.0f - effective_dmg);
            float p_c_modern = 1.0f - baseline_tc;

            if (p_c_modern > 0.01f) {
                log_lr += std::log(p_c_ancient / p_c_modern) * quality_weight;
            }
        }
        // A and G at 5' end are uninformative for C→T damage
    }

    // ==========================================================================
    // 3' END LIKELIHOOD: G→A DAMAGE
    // Analyze last 10 positions with GC-bin specific damage rates
    // ==========================================================================

    const size_t analyze_3prime = std::min(size_t(10), len);
    float lambda_3p = sample_profile.lambda_3prime > 0 ? sample_profile.lambda_3prime : lambda;

    for (size_t i = 0; i < analyze_3prime; ++i) {
        size_t pos = len - 1 - i;  // Position from end of sequence
        char base = fast_upper(seq[pos]);

        // Quality-based confidence weight
        float quality_weight = 1.0f;
        if (has_quality) {
            float error_prob = phred_to_error_prob(quality[pos]);
            quality_weight = 1.0f - std::min(error_prob * 2.0f, 0.5f);
        }

        // Compute damage rate at this position using GC-bin delta_s
        float decay = std::exp(-lambda_3p * static_cast<float>(i));
        float dmg_rate = delta_s * decay;

        // Codon position modulation
        int codon_pos = (static_cast<int>(pos) - frame + 300) % 3;
        float codon_mod = 1.0f;
        if (codon_pos == 2) {
            codon_mod = 1.2f;  // Wobble more tolerant
        } else if (codon_pos == 0) {
            codon_mod = 0.9f;
        }

        float effective_dmg = std::min(dmg_rate * codon_mod, 0.95f);

        if (base == 'A') {
            // A observed: could be original A or G→A damage
            // P(A|ancient) = baseline_A + baseline_G * damage_rate
            // P(A|modern) = baseline_A
            float p_a_ancient = baseline_ag + (1.0f - baseline_ag) * effective_dmg;
            float p_a_modern = baseline_ag;

            if (p_a_modern > 0.01f) {
                log_lr += std::log(p_a_ancient / p_a_modern) * quality_weight;
            }
        } else if (base == 'G') {
            // G observed: G that survived
            // P(G|ancient) = baseline_G * (1 - damage_rate)
            // P(G|modern) = baseline_G
            float p_g_ancient = (1.0f - baseline_ag) * (1.0f - effective_dmg);
            float p_g_modern = 1.0f - baseline_ag;

            if (p_g_modern > 0.01f) {
                log_lr += std::log(p_g_ancient / p_g_modern) * quality_weight;
            }
        }
    }

    // ==========================================================================
    // MULTI-SIGNAL DAMAGE DETECTION
    // Three signals: ORF extension, context hexamer, hexamer improvement
    // ==========================================================================

    int terminal_stop_pos = -1;
    int terminal_stop_codon_idx = -1;
    for (int codon_idx = 0; codon_idx < 4; ++codon_idx) {
        size_t pos = static_cast<size_t>(frame) + codon_idx * 3;
        if (pos + 3 > len) break;

        int ds = is_damaged_stop_codon(fast_upper(seq[pos]),
                                        fast_upper(seq[pos + 1]),
                                        fast_upper(seq[pos + 2]));
        if (ds == 2) {
            terminal_stop_pos = static_cast<int>(pos);
            terminal_stop_codon_idx = codon_idx;
            break;
        }
    }

    if (terminal_stop_pos >= 0) {
        float damage_evidence = 0.0f;
        float no_damage_evidence = 0.0f;
        float pos_weight = (terminal_stop_codon_idx == 0) ? 1.0f :
                          (terminal_stop_codon_idx == 1) ? 0.8f :
                          (terminal_stop_codon_idx == 2) ? 0.6f : 0.4f;

        // Signal 1: ORF extension
        int next_stop_pos = -1;
        for (size_t i = terminal_stop_pos + 3; i + 3 <= len; i += 3) {
            char c1 = fast_upper(seq[i]);
            char c2 = fast_upper(seq[i + 1]);
            char c3 = fast_upper(seq[i + 2]);
            if ((c1 == 'T' && c2 == 'A' && c3 == 'A') ||
                (c1 == 'T' && c2 == 'A' && c3 == 'G') ||
                (c1 == 'T' && c2 == 'G' && c3 == 'A')) {
                next_stop_pos = static_cast<int>(i);
                break;
            }
        }
        int orf_extension = (next_stop_pos >= 0) ?
            (next_stop_pos - terminal_stop_pos) : (static_cast<int>(len) - terminal_stop_pos);
        if (orf_extension > 60) damage_evidence += 1.5f;
        else if (orf_extension > 30) damage_evidence += 0.8f;
        else if (orf_extension < 15) no_damage_evidence += 0.5f;

        // Signal 2: Context hexamer
        if (terminal_stop_pos >= 3 && len >= 6) {
            size_t ctx_start = (terminal_stop_pos >= 6) ? terminal_stop_pos - 6 : 0;
            if (ctx_start + 6 <= len) {
                char hexbuf[7];
                for (int i = 0; i < 6; i++) hexbuf[i] = fast_upper(seq[ctx_start + i]);
                hexbuf[6] = '\0';
                float ctx_freq = get_hexamer_freq(hexbuf);
                if (ctx_freq > 0.001f) damage_evidence += 0.6f;
                else if (ctx_freq < 0.0001f) no_damage_evidence += 0.3f;
            }
        }

        // Signal 3: Hexamer improvement
        if (len >= 6) {
            std::string corrected = seq;
            corrected[terminal_stop_pos] = 'C';
            float improvement = calculate_hexamer_llr_score(corrected, frame) -
                               calculate_hexamer_llr_score(seq, frame);
            if (improvement > 1.0f) damage_evidence += 0.5f * (improvement - 1.0f);
            else if (improvement < 0.3f) no_damage_evidence += 0.4f;
        }

        log_lr += std::clamp((damage_evidence - no_damage_evidence) * pos_weight, -1.5f, 3.0f);
    }

    // Non-stop T-first codons
    if (terminal_stop_pos < 0) {
        for (int codon_idx = 0; codon_idx < 3; ++codon_idx) {
            size_t pos = static_cast<size_t>(frame) + codon_idx * 3;
            if (pos + 6 > len) break;
            if (fast_upper(seq[pos]) != 'T') continue;

            std::string corrected = seq;
            corrected[pos] = 'C';
            float improvement = calculate_hexamer_llr_score(corrected, frame) -
                               calculate_hexamer_llr_score(seq, frame);
            float pos_weight = std::exp(-0.3f * static_cast<float>(codon_idx));
            if (improvement > 0.3f) log_lr += std::min(0.8f, improvement * 0.5f) * pos_weight;
            else if (improvement < -0.2f) log_lr += improvement * 0.3f * pos_weight;
        }
    }

    // Surviving pre-stops
    for (int codon_idx = 0; codon_idx < 3; ++codon_idx) {
        size_t pos = static_cast<size_t>(frame) + codon_idx * 3;
        if (pos + 3 > len) break;
        char c1 = fast_upper(seq[pos]);
        char c2 = fast_upper(seq[pos + 1]);
        char c3 = fast_upper(seq[pos + 2]);
        if ((c1 == 'C' && c2 == 'A' && c3 == 'A') ||
            (c1 == 'C' && c2 == 'A' && c3 == 'G') ||
            (c1 == 'C' && c2 == 'G' && c3 == 'A')) {
            log_lr -= 0.2f * std::exp(-0.3f * static_cast<float>(codon_idx));
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
    // PRIOR FOR PER-READ CLASSIFICATION
    // Use a balanced prior (0.5) so per-read LLR determines classification.
    // The GC-bin prior (p0) is too strong for individual read classification -
    // it treats ALL reads in a bin the same, ignoring per-read features.
    // ==========================================================================

    // Use neutral prior - let LLR do the discrimination
    // The GC-conditional sample-level statistics inform delta_s but not the prior
    float prior_clamped = 0.5f;
    float prior_log_odds = 0.0f;  // log(0.5 / 0.5) = 0

    // Combine likelihood and prior: posterior = sigmoid(logit(prior) + LLR)
    float total_log_odds = log_lr + prior_log_odds;

    // Convert log-odds to probability - this is the GC-conditional posterior
    float posterior = 1.0f / (1.0f + std::exp(-total_log_odds));

    // ==========================================================================
    // CHANNEL B MODULATION
    // Apply Channel B factor to ensure sample-level validation affects per-read scores
    // - Validated damage samples: full posterior
    // - Artifact samples: heavily suppressed (already handled in early exit)
    // - Uncertain samples: partial confidence
    // ==========================================================================
    posterior *= channel_b_factor;

    return std::clamp(posterior, 0.0f, 1.0f);
}

} // namespace agp
