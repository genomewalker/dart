// Sample-level damage profile management

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include "agp/hexamer_tables.hpp"
#include <algorithm>
#include <cmath>
#include <array>

namespace agp {

// Compute binomial log-likelihood for a single observation
// k = successes (e.g., T count), n = total trials (e.g., T+C count), p = probability
static inline double binomial_ll(double k, double n, double p) {
    if (n < 1 || p <= 0 || p >= 1) return 0.0;
    // Log-likelihood: k*log(p) + (n-k)*log(1-p) + constant (ignored for LLR)
    return k * std::log(p) + (n - k) * std::log(1.0 - p);
}

// Compute log-likelihood ratio: exponential decay model vs constant model
// Returns LLR > 0 if exponential fits better (real decay pattern)
// Returns LLR ~0 if constant fits as well (no pattern)
// Returns LLR < 0 if data is inverted (terminal lower than interior)
// Uses positions 1-9 for fitting (excludes position 0 which has artifacts)
// Computes its own best-fit amplitude directly from data (without clamping)
static float compute_decay_llr(
    const std::array<double, 15>& freq,      // T/(T+C) or A/(A+G) at each position (normalized)
    const std::array<double, 15>& total,     // T+C or A+G counts at each position
    float baseline,                          // Middle-of-read baseline
    float /*amplitude*/,                     // Not used - we compute directly
    float lambda) {                          // Decay constant for exponential model

    const double MIN_COVERAGE = 100.0;

    // First, compute best-fit amplitude from data (without clamping to positive)
    // This lets us detect both positive decay (damage) and negative decay (inverted)
    double sum_signal = 0.0, sum_weight = 0.0;
    for (int i = 1; i < 10; ++i) {
        if (total[i] < MIN_COVERAGE) continue;
        double weight = std::exp(-lambda * i);  // Weight by expected decay contribution
        double excess = freq[i] - baseline;
        sum_signal += total[i] * excess / weight;  // Infer amplitude from each position
        sum_weight += total[i];
    }
    double raw_amplitude = (sum_weight > 0) ? sum_signal / sum_weight : 0.0;

    // Now compute likelihoods
    double ll_exp = 0.0;   // Log-likelihood under exponential model
    double ll_const = 0.0; // Log-likelihood under constant model

    // Sum log-likelihoods over positions 1-9 (exclude pos 0)
    for (int i = 1; i < 10; ++i) {
        if (total[i] < MIN_COVERAGE) continue;

        double n = total[i];
        double k = freq[i] * n;  // freq is already normalized to 0-1

        // Exponential model: p = baseline + amplitude * exp(-lambda * i)
        double p_exp = baseline + raw_amplitude * std::exp(-lambda * i);
        p_exp = std::clamp(p_exp, 0.001, 0.999);
        ll_exp += binomial_ll(k, n, p_exp);

        // Constant model: p = baseline (no decay)
        double p_const = std::clamp(static_cast<double>(baseline), 0.001, 0.999);
        ll_const += binomial_ll(k, n, p_const);
    }

    // Log-likelihood ratio: positive means exponential fits better
    // For INVERTED patterns (raw_amplitude < 0), the exponential would fit INVERTED decay
    // We want positive LLR only for POSITIVE decay, so we penalize inverted patterns
    float llr = static_cast<float>(ll_exp - ll_const);

    // If amplitude is negative (inverted pattern), negate the LLR
    // This makes decay_llr negative for inverted patterns, positive for true damage
    if (raw_amplitude < 0) {
        return -std::abs(llr);  // Negative LLR for inverted patterns
    }
    return llr;
}

// Fit exponential decay model: p(pos) = b + A * exp(-lambda * pos)
// Uses weighted least squares with coverage-based weights
// Returns: {baseline (b), amplitude (A), lambda, rmse}
// If external_baseline >= 0, use it instead of estimating from positions 10-14
// This allows using middle-of-read baseline which is more reliable
static std::array<float, 4> fit_exponential_decay(
    const std::array<double, 15>& freq,      // T/(T+C) or A/(A+G) at each position
    const std::array<double, 15>& coverage,  // T+C or A+G counts at each position
    float lambda_init = 0.2f,
    float external_baseline = -1.0f) {       // If >= 0, use this as baseline

    // Minimum coverage to include a position
    const double MIN_COVERAGE = 100.0;

    // Count valid positions
    int n_valid = 0;
    for (int i = 0; i < 15; ++i) {
        if (coverage[i] >= MIN_COVERAGE) n_valid++;
    }
    if (n_valid < 5) {
        // Not enough data - return zeros
        return {0.0f, 0.0f, 0.0f, 1.0f};
    }

    // Use external baseline (middle-of-read) if provided, otherwise estimate from positions 10-14
    // CRITICAL: Middle-of-read baseline is more reliable because positions 10-14 may still
    // have read-end composition artifacts that inflate the apparent "damage" signal
    float b;
    if (external_baseline >= 0.0f) {
        b = external_baseline;
    } else {
        // Fallback: estimate baseline from positions 10-14
        double baseline_sum = 0.0, baseline_weight = 0.0;
        for (int i = 10; i < 15; ++i) {
            if (coverage[i] >= MIN_COVERAGE) {
                baseline_sum += freq[i] * coverage[i];
                baseline_weight += coverage[i];
            }
        }
        b = (baseline_weight > 0) ?
            static_cast<float>(baseline_sum / baseline_weight) : 0.5f;
    }

    // CRITICAL: Estimate amplitude from position 1, NOT position 0
    // Position 0 is uniquely prone to first-cycle/ligation/trimming artifacts
    // that can arbitrarily bias the estimate. Position 1 is more reliable.
    float A = (coverage[1] >= MIN_COVERAGE) ?
              std::max(0.0f, static_cast<float>(freq[1]) - b) : 0.0f;

    // Refine lambda using linear regression on log(p - b)
    // CRITICAL: Exclude position 0 from fit - it can have arbitrary bias
    // Use positions 1-9 where decay is significant
    float lambda = lambda_init;
    if (A > 0.01f) {
        double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_w = 0.0;
        for (int i = 1; i < 10; ++i) {  // Start from 1, not 0
            if (coverage[i] < MIN_COVERAGE) continue;
            double excess = freq[i] - b;
            if (excess > 0.005) {  // Only use points above baseline
                double w = coverage[i];
                double y = std::log(excess);
                sum_x += w * i;
                sum_y += w * y;
                sum_xy += w * i * y;
                sum_xx += w * i * i;
                sum_w += w;
            }
        }
        if (sum_w > 0 && (sum_w * sum_xx - sum_x * sum_x) > 1e-6) {
            // Weighted linear regression: y = log(A) - lambda * x
            double slope = (sum_w * sum_xy - sum_x * sum_y) /
                          (sum_w * sum_xx - sum_x * sum_x);
            lambda = std::max(0.05f, std::min(0.5f, static_cast<float>(-slope)));
            // Re-estimate A from intercept
            double intercept = (sum_y - slope * sum_x) / sum_w;
            float A_new = static_cast<float>(std::exp(intercept));
            if (A_new > 0.0f && A_new < 1.0f - b) {
                A = A_new;
            }
        }
    }

    // Constrain A to valid range [0, 1-b]
    A = std::max(0.0f, std::min(A, 1.0f - b - 0.001f));

    // Compute RMSE of fit (excluding position 0 for consistency)
    double sse = 0.0, weight_sum = 0.0;
    for (int i = 1; i < 15; ++i) {  // Start from 1, not 0
        if (coverage[i] >= MIN_COVERAGE) {
            double pred = b + A * std::exp(-lambda * i);
            double resid = freq[i] - pred;
            sse += resid * resid * coverage[i];
            weight_sum += coverage[i];
        }
    }
    float rmse = (weight_sum > 0) ?
                 static_cast<float>(std::sqrt(sse / weight_sum)) : 1.0f;

    return {b, A, lambda, rmse};
}

void FrameSelector::update_sample_profile(
    SampleDamageProfile& profile,
    const std::string& seq) {

    if (seq.length() < 30) return;  // Too short for reliable statistics

    size_t len = seq.length();

    // Count bases at 5' end positions (first 15 bases)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = fast_upper(seq[i]);
        // Damage signal: T/(T+C) for C→T
        if (base == 'T') {
            profile.t_freq_5prime[i]++;
            profile.tc_total_5prime[i]++;
        } else if (base == 'C') {
            profile.c_freq_5prime[i]++;
            profile.tc_total_5prime[i]++;
        }
        // Negative control: A/(A+G) - should NOT be elevated by C→T damage
        if (base == 'A') {
            profile.a_freq_5prime[i]++;
        } else if (base == 'G') {
            profile.g_freq_5prime[i]++;
        }
    }

    // Count bases at 3' end positions (last 15 bases)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        // Damage signal: A/(A+G) for G→A
        if (base == 'A') {
            profile.a_freq_3prime[i]++;
            profile.ag_total_3prime[i]++;
        } else if (base == 'G') {
            profile.g_freq_3prime[i]++;
            profile.ag_total_3prime[i]++;
        }
        // Negative control: T/(T+C) - should NOT be elevated by G→A damage
        if (base == 'T') {
            profile.t_freq_3prime[i]++;
        } else if (base == 'C') {
            profile.c_freq_3prime[i]++;
        }
    }

    // Count bases in middle third (undamaged baseline)
    size_t mid_start = len / 3;
    size_t mid_end = 2 * len / 3;
    for (size_t i = mid_start; i < mid_end; ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'T') profile.baseline_t_freq++;
        else if (base == 'C') profile.baseline_c_freq++;
        else if (base == 'A') profile.baseline_a_freq++;
        else if (base == 'G') profile.baseline_g_freq++;
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

    // =========================================================================
    // CHANNEL B: Convertible stop codon tracking
    // Track CAA→TAA, CAG→TAG, CGA→TGA pairs by nucleotide position
    // This is BEFORE frame selection, so no circularity issue
    // For forward frames (f=0,1,2), the C/T position is at f + 3*k
    // =========================================================================
    if (len >= 18) {
        // Track convertible codons for all 3 forward frames at 5' end
        for (int frame = 0; frame < 3; ++frame) {
            // Scan codons from 5' end up to position 14
            for (size_t k = 0; ; ++k) {
                size_t codon_start = frame + 3 * k;
                if (codon_start + 3 > len || codon_start > 14) break;

                // Position of the first base (C or T in CAx/TAx codons)
                size_t p = codon_start;
                if (p >= 15) break;

                // Extract codon
                char b0 = fast_upper(seq[codon_start]);
                char b1 = fast_upper(seq[codon_start + 1]);
                char b2 = fast_upper(seq[codon_start + 2]);

                // Skip if any base is ambiguous
                if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                    (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                    (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                    continue;
                }

                // Count total codons at this position
                profile.total_codons_5prime[p]++;

                // Check for convertible pairs (CAA/TAA, CAG/TAG, CGA/TGA)
                // CAA (Gln) → TAA (Stop) via C→T at position 0
                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'C') profile.convertible_caa_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_taa_5prime[p]++;
                }
                // CAG (Gln) → TAG (Stop) via C→T at position 0
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'C') profile.convertible_cag_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tag_5prime[p]++;
                }
                // CGA (Arg) → TGA (Stop) via C→T at position 0
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'C') profile.convertible_cga_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tga_5prime[p]++;
                }
            }
        }

        // Track interior convertible codons (positions 30+ from start)
        // This gives us the baseline stop conversion rate
        for (int frame = 0; frame < 3; ++frame) {
            for (size_t k = 0; ; ++k) {
                size_t codon_start = frame + 3 * k;
                // Interior region: 30 to len-30 (away from both ends)
                if (codon_start < 30 || codon_start + 3 > len - 30) {
                    if (codon_start + 3 > len - 30) break;
                    continue;
                }

                char b0 = fast_upper(seq[codon_start]);
                char b1 = fast_upper(seq[codon_start + 1]);
                char b2 = fast_upper(seq[codon_start + 2]);

                if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                    (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                    (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                    continue;
                }

                profile.total_codons_interior++;

                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'C') profile.convertible_caa_interior++;
                    else if (b0 == 'T') profile.convertible_taa_interior++;
                }
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'C') profile.convertible_cag_interior++;
                    else if (b0 == 'T') profile.convertible_tag_interior++;
                }
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'C') profile.convertible_cga_interior++;
                    else if (b0 == 'T') profile.convertible_tga_interior++;
                }
            }
        }
    }

    // =========================================================================
    // Hexamer-based damage detection
    // Collect hexamers from ALL reads - works regardless of sample composition
    // Using position 0 (standard); inversion correction handles unusual cases
    // =========================================================================
    if (len >= 18) {  // Need at least 18 bases for hexamers
        // 5' terminal hexamer starting at position 0
        char hex_5prime[7];
        bool valid_5prime = true;
        for (int i = 0; i < 6; ++i) {
            char b = fast_upper(seq[i]);
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') {
                valid_5prime = false;
                break;
            }
            hex_5prime[i] = b;
        }
        hex_5prime[6] = '\0';

        if (valid_5prime) {
            uint32_t code_5prime = encode_hexamer(hex_5prime);
            if (code_5prime < 4096) {
                profile.hexamer_count_5prime[code_5prime] += 1.0;
                profile.n_hexamers_5prime++;
            }
        }

        // Interior hexamer (starting at len/2 - 3, approximately middle)
        size_t interior_start = len / 2 - 3;
        char hex_interior[7];
        bool valid_interior = true;
        for (int i = 0; i < 6; ++i) {
            char b = fast_upper(seq[interior_start + i]);
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') {
                valid_interior = false;
                break;
            }
            hex_interior[i] = b;
        }
        hex_interior[6] = '\0';

        if (valid_interior) {
            uint32_t code_interior = encode_hexamer(hex_interior);
            if (code_interior < 4096) {
                profile.hexamer_count_interior[code_interior] += 1.0;
                profile.n_hexamers_interior++;
            }
        }
    }

    profile.n_reads++;
}

void FrameSelector::finalize_sample_profile(SampleDamageProfile& profile) {
    if (profile.n_reads == 0) return;

    // Capture raw counts before normalization for statistical tests
    double base_tc_total = profile.baseline_t_freq + profile.baseline_c_freq;
    double base_ag_total = profile.baseline_a_freq + profile.baseline_g_freq;

    auto compute_terminal_stats = [](double term_t, double term_c, double base_t, double base_c) {
        double term_total = term_t + term_c;
        double base_total = base_t + base_c;
        if (term_total < 1.0 || base_total < 1.0) {
            return std::pair<float, float>{0.0f, 0.0f};
        }
        double p_term = term_t / term_total;
        double p_base = base_t / base_total;
        double p_pool = (term_t + base_t) / (term_total + base_total);
        double se = std::sqrt(std::max(1e-9, p_pool * (1.0 - p_pool) * (1.0 / term_total + 1.0 / base_total)));
        double z = se > 0.0 ? (p_term - p_base) / se : 0.0;
        return std::pair<float, float>{static_cast<float>(p_term - p_base), static_cast<float>(z)};
    };

    // Terminal enrichment/depletion statistics (position 0)
    auto stats_5 = compute_terminal_stats(profile.t_freq_5prime[0], profile.c_freq_5prime[0],
                                          profile.baseline_t_freq, profile.baseline_c_freq);
    profile.terminal_shift_5prime = stats_5.first;
    profile.terminal_z_5prime = stats_5.second;

    auto stats_3 = compute_terminal_stats(profile.a_freq_3prime[0], profile.g_freq_3prime[0],
                                          profile.baseline_a_freq, profile.baseline_g_freq);
    profile.terminal_shift_3prime = stats_3.first;
    profile.terminal_z_3prime = stats_3.second;

    // Negative control statistics (should NOT show enrichment if damage is real)
    // 5' control: A/(A+G) at 5' end vs interior
    auto ctrl_5 = compute_terminal_stats(profile.a_freq_5prime[0], profile.g_freq_5prime[0],
                                         profile.baseline_a_freq, profile.baseline_g_freq);
    profile.ctrl_shift_5prime = ctrl_5.first;
    profile.ctrl_z_5prime = ctrl_5.second;

    // 3' control: T/(T+C) at 3' end vs interior
    auto ctrl_3 = compute_terminal_stats(profile.t_freq_3prime[0], profile.c_freq_3prime[0],
                                         profile.baseline_t_freq, profile.baseline_c_freq);
    profile.ctrl_shift_3prime = ctrl_3.first;
    profile.ctrl_z_3prime = ctrl_3.second;

    // Channel divergence: difference between damage shift and control shift
    // High divergence = real damage (only damage channel elevated)
    // Low divergence = composition bias (both channels elevated together)
    profile.channel_divergence_5prime = std::abs(profile.terminal_shift_5prime - profile.ctrl_shift_5prime);
    profile.channel_divergence_3prime = std::abs(profile.terminal_shift_3prime - profile.ctrl_shift_3prime);

    // Normalize baseline frequencies (using double to preserve precision)
    double mid_total = profile.baseline_t_freq + profile.baseline_c_freq +
                       profile.baseline_a_freq + profile.baseline_g_freq;
    if (mid_total > 0) {
        profile.baseline_t_freq /= mid_total;
        profile.baseline_c_freq /= mid_total;
        profile.baseline_a_freq /= mid_total;
        profile.baseline_g_freq /= mid_total;
    }

    // Compute baseline T/C and A/G ratios from the middle of reads (for inverted pattern detection)
    double baseline_tc = profile.baseline_t_freq /
                        (profile.baseline_t_freq + profile.baseline_c_freq + 0.001);
    double baseline_ag = profile.baseline_a_freq /
                        (profile.baseline_a_freq + profile.baseline_g_freq + 0.001);

    // Step 1: Normalize position-specific frequencies (but don't compute damage rates yet)
    for (int i = 0; i < 15; ++i) {
        double tc_total = profile.t_freq_5prime[i] + profile.c_freq_5prime[i];
        if (tc_total > 0) {
            profile.t_freq_5prime[i] = profile.t_freq_5prime[i] / tc_total;
            profile.c_freq_5prime[i] = 1.0 - profile.t_freq_5prime[i];
        }
        double ag_total = profile.a_freq_3prime[i] + profile.g_freq_3prime[i];
        if (ag_total > 0) {
            profile.a_freq_3prime[i] = profile.a_freq_3prime[i] / ag_total;
            profile.g_freq_3prime[i] = 1.0 - profile.a_freq_3prime[i];
        }
    }

    // Step 2: Compute exponential fit using MIDDLE-OF-READ baseline
    // p(pos) = b + A * exp(-lambda * pos)
    // CRITICAL: Use middle-of-read baseline instead of positions 10-14
    // Positions 10-14 are still at read ends and may have composition artifacts
    auto fit_5p = fit_exponential_decay(profile.t_freq_5prime, profile.tc_total_5prime, 0.2f,
                                        static_cast<float>(baseline_tc));
    profile.fit_baseline_5prime = fit_5p[0];
    profile.fit_amplitude_5prime = fit_5p[1];

    // DEBUG: Verify baseline values match
    fprintf(stderr, "  [DEBUG-BASELINE] baseline_tc=%.4f, fit_baseline_5prime=%.4f, t_freq[0]=%.4f\n",
            baseline_tc, profile.fit_baseline_5prime, profile.t_freq_5prime[0]);
    float fit_lambda_5p = fit_5p[2];
    profile.fit_rmse_5prime = fit_5p[3];

    auto fit_3p = fit_exponential_decay(profile.a_freq_3prime, profile.ag_total_3prime, 0.2f,
                                        static_cast<float>(baseline_ag));
    profile.fit_baseline_3prime = fit_3p[0];
    profile.fit_amplitude_3prime = fit_3p[1];
    float fit_lambda_3p = fit_3p[2];
    profile.fit_rmse_3prime = fit_3p[3];

    // Compute decay log-likelihood ratio (exponential vs constant model)
    // Positive LLR = exponential fits better (real decay pattern)
    // This helps distinguish real damage (exponential decay) from composition bias (uniform)
    profile.decay_llr_5prime = compute_decay_llr(
        profile.t_freq_5prime, profile.tc_total_5prime,
        profile.fit_baseline_5prime, profile.fit_amplitude_5prime, fit_lambda_5p);
    profile.decay_llr_3prime = compute_decay_llr(
        profile.a_freq_3prime, profile.ag_total_3prime,
        profile.fit_baseline_3prime, profile.fit_amplitude_3prime, fit_lambda_3p);

    // Compute control channel decay LLR
    // Control channels: A/(A+G) at 5', T/(T+C) at 3'
    // If control channel also shows decay, it's likely composition/trimming artifact, not damage
    {
        // Normalize control channel frequencies
        std::array<double, 15> ctrl_freq_5p = {};  // A/(A+G) at 5'
        std::array<double, 15> ctrl_total_5p = {}; // A+G at 5'
        std::array<double, 15> ctrl_freq_3p = {};  // T/(T+C) at 3'
        std::array<double, 15> ctrl_total_3p = {}; // T+C at 3'

        for (int i = 0; i < 15; ++i) {
            double ag_5p = profile.a_freq_5prime[i] + profile.g_freq_5prime[i];
            ctrl_total_5p[i] = ag_5p;
            ctrl_freq_5p[i] = (ag_5p > 0) ? profile.a_freq_5prime[i] / ag_5p : 0.5;

            double tc_3p = profile.t_freq_3prime[i] + profile.c_freq_3prime[i];
            ctrl_total_3p[i] = tc_3p;
            ctrl_freq_3p[i] = (tc_3p > 0) ? profile.t_freq_3prime[i] / tc_3p : 0.5;
        }

        // Compute control channel baselines from middle of read
        double ctrl_baseline_5p = baseline_ag;  // A/(A+G) baseline
        double ctrl_baseline_3p = baseline_tc;  // T/(T+C) baseline

        // Compute control channel decay LLR (using same lambda as damage channel)
        profile.ctrl_decay_llr_5prime = compute_decay_llr(
            ctrl_freq_5p, ctrl_total_5p,
            static_cast<float>(ctrl_baseline_5p), 0.0f, fit_lambda_5p);
        profile.ctrl_decay_llr_3prime = compute_decay_llr(
            ctrl_freq_3p, ctrl_total_3p,
            static_cast<float>(ctrl_baseline_3p), 0.0f, fit_lambda_3p);

        // Delta LLR: damage channel - control channel
        // Positive = damage channel has stronger decay = real damage
        // Near zero = both channels decay similarly = composition artifact
        profile.delta_llr_5prime = profile.decay_llr_5prime - profile.ctrl_decay_llr_5prime;
        profile.delta_llr_3prime = profile.decay_llr_3prime - profile.ctrl_decay_llr_3prime;
    }

    // =========================================================================
    // CHANNEL B: Convertible stop codon decay analysis
    // Compute stop conversion rate decay and compare to Channel A
    // This is the "smoking gun" test: real C→T damage MUST create stops
    // in damage-susceptible contexts (CAA→TAA, CAG→TAG, CGA→TGA)
    // =========================================================================
    {
        // Compute interior baseline: stop / (pre-image + stop) ratio
        double total_pre_interior = profile.convertible_caa_interior +
                                   profile.convertible_cag_interior +
                                   profile.convertible_cga_interior;
        double total_stop_interior = profile.convertible_taa_interior +
                                    profile.convertible_tag_interior +
                                    profile.convertible_tga_interior;
        double total_convertible_interior = total_pre_interior + total_stop_interior;

        if (total_convertible_interior > 100) {
            profile.stop_conversion_rate_baseline = static_cast<float>(
                total_stop_interior / total_convertible_interior);
            profile.channel_b_valid = true;
        }

        fprintf(stderr, "  [CHANNEL-B] Interior: pre=%.0f, stop=%.0f, baseline=%.4f\n",
                total_pre_interior, total_stop_interior, profile.stop_conversion_rate_baseline);

        // Compute stop conversion rate at each terminal position
        std::array<double, 15> stop_rate = {};
        std::array<double, 15> stop_exposure = {};  // pre + stop

        for (int p = 0; p < 15; ++p) {
            double pre = profile.convertible_caa_5prime[p] +
                        profile.convertible_cag_5prime[p] +
                        profile.convertible_cga_5prime[p];
            double stop = profile.convertible_taa_5prime[p] +
                         profile.convertible_tag_5prime[p] +
                         profile.convertible_tga_5prime[p];
            stop_exposure[p] = pre + stop;
            if (stop_exposure[p] > 10) {
                stop_rate[p] = stop / stop_exposure[p];
            } else {
                stop_rate[p] = profile.stop_conversion_rate_baseline;  // Use baseline if no data
            }
        }

        // Debug: print stop rates at first few positions
        fprintf(stderr, "  [CHANNEL-B] Terminal stop rates: p0=%.4f, p1=%.4f, p2=%.4f, p3=%.4f\n",
                stop_rate[0], stop_rate[1], stop_rate[2], stop_rate[3]);
        fprintf(stderr, "  [CHANNEL-B] Terminal exposure: p0=%.0f, p1=%.0f, p2=%.0f, p3=%.0f\n",
                stop_exposure[0], stop_exposure[1], stop_exposure[2], stop_exposure[3]);

        // Compute stop decay LLR using same lambda from Channel A
        // If both channels have same decay shape, it confirms real damage
        // BUG FIX: Use fit_lambda_5p directly since profile.lambda_5prime isn't updated yet
        if (profile.channel_b_valid) {
            float baseline_b = profile.stop_conversion_rate_baseline;
            float lambda_b = fit_lambda_5p;  // Use fit lambda from Channel A (NOT profile.lambda_5prime which isn't set yet)

            // Estimate amplitude from positions 1-4 (exclude pos 0 like Channel A does)
            // BUG FIX: Exclude position 0 which has first-base artifacts
            double sum_excess = 0.0, sum_weight = 0.0;
            for (int i = 1; i < 5; ++i) {  // Start from 1, not 0
                if (stop_exposure[i] > 50) {
                    double weight = std::exp(-lambda_b * i);
                    double excess = stop_rate[i] - baseline_b;
                    sum_excess += stop_exposure[i] * excess / weight;
                    sum_weight += stop_exposure[i];
                }
            }
            float amplitude_b = (sum_weight > 0) ? static_cast<float>(sum_excess / sum_weight) : 0.0f;
            profile.stop_amplitude_5prime = std::max(0.0f, amplitude_b);

            fprintf(stderr, "  [CHANNEL-B] Fitted amplitude=%.4f (excess stops at terminal, lambda=%.3f)\n",
                    profile.stop_amplitude_5prime, lambda_b);

            // Compute log-likelihood ratio for exponential decay vs constant
            // BUG FIX: Use actual stop counts, not stop_rate * n
            // BUG FIX: Exclude position 0 for consistency with Channel A
            double ll_exp = 0.0, ll_const = 0.0;
            for (int p = 1; p < 10; ++p) {  // Start from 1, not 0
                if (stop_exposure[p] < 50) continue;

                double n = stop_exposure[p];
                // BUG FIX: Use actual stop count directly
                double pre = profile.convertible_caa_5prime[p] +
                            profile.convertible_cag_5prime[p] +
                            profile.convertible_cga_5prime[p];
                double stop = profile.convertible_taa_5prime[p] +
                             profile.convertible_tag_5prime[p] +
                             profile.convertible_tga_5prime[p];
                double k = stop;  // Actual stop count

                // Exponential model: p = baseline + amplitude * exp(-lambda * p)
                double p_exp = baseline_b + amplitude_b * std::exp(-lambda_b * p);
                p_exp = std::clamp(p_exp, 0.001, 0.999);
                ll_exp += binomial_ll(k, n, p_exp);

                // Constant model: p = baseline
                double p_const = std::clamp(static_cast<double>(baseline_b), 0.001, 0.999);
                ll_const += binomial_ll(k, n, p_const);
            }

            profile.stop_decay_llr_5prime = static_cast<float>(ll_exp - ll_const);

            // If amplitude is negative (inverted), negate LLR
            if (amplitude_b < 0) {
                profile.stop_decay_llr_5prime = -std::abs(profile.stop_decay_llr_5prime);
            }

            fprintf(stderr, "  [CHANNEL-B] Stop decay LLR=%.2f (exponential vs constant)\n",
                    profile.stop_decay_llr_5prime);

            // =========================================================================
            // JOINT EVIDENCE DECISION
            // Compare Channel A (nucleotide T/(T+C)) vs Channel B (stop conversion)
            // Real damage: both channels show positive decay
            // Compositional artifact: Channel A positive, Channel B flat/negative
            //
            // BUG FIX: Use delta_llr (control-adjusted) instead of raw decay_llr
            // This accounts for compositional variation in both channels
            // =========================================================================
            float channel_a_llr = profile.delta_llr_5prime;  // Control-adjusted Channel A
            float channel_b_llr = profile.stop_decay_llr_5prime;

            fprintf(stderr, "  [JOINT] Channel A (delta_llr) = %.2f, Channel B (stop_llr) = %.2f\n",
                    channel_a_llr, channel_b_llr);

            // Thresholds for deciding damage is real
            // Use lower threshold for delta_llr since it's already control-adjusted
            const float CHANNEL_A_THRESHOLD = 50.0f;    // Delta LLR threshold (lower than raw)
            const float CHANNEL_B_THRESHOLD = 10.0f;    // Moderate stop decay (weaker signal)

            // BUG FIX: Per spec, "Channel A fires but Channel B flat" should be ARTIFACT
            // Previously only flagged when Channel B was actively negative
            if (channel_a_llr > CHANNEL_A_THRESHOLD && channel_b_llr > CHANNEL_B_THRESHOLD) {
                // Both channels agree: real damage
                profile.damage_validated = true;
                profile.damage_artifact = false;
                fprintf(stderr, "  [JOINT] DECISION: REAL DAMAGE (both channels fire)\n");
            } else if (channel_a_llr > CHANNEL_A_THRESHOLD && channel_b_llr <= CHANNEL_B_THRESHOLD) {
                // Channel A fires but Channel B is flat or negative: artifact
                // BUG FIX: Treat "flat" as artifact, not just "negative"
                profile.damage_validated = false;
                profile.damage_artifact = true;
                fprintf(stderr, "  [JOINT] DECISION: COMPOSITIONAL ARTIFACT (Channel A only, B=%.2f)\n",
                        channel_b_llr);
            } else {
                // No strong signal in Channel A
                profile.damage_validated = false;
                profile.damage_artifact = false;
                fprintf(stderr, "  [JOINT] DECISION: NO/LOW DAMAGE SIGNAL (A=%.2f, B=%.2f)\n",
                        channel_a_llr, channel_b_llr);
            }
        }
    }

    // Update lambda from fit if reasonable
    if (fit_lambda_5p > 0.05f && fit_lambda_5p < 0.5f) {
        profile.lambda_5prime = fit_lambda_5p;
    }
    if (fit_lambda_3p > 0.05f && fit_lambda_3p < 0.5f) {
        profile.lambda_3prime = fit_lambda_3p;
    }

    // Step 3: Compute per-position damage rates using FIT baseline (not middle-of-read)
    // This produces rates that are comparable to metaDMG and consistent with d_max
    float fit_baseline_c_frac_5p = 1.0f - profile.fit_baseline_5prime;
    float fit_baseline_g_frac_3p = 1.0f - profile.fit_baseline_3prime;

    for (int i = 0; i < 15; ++i) {
        // 5' end: C→T damage using fit baseline
        if (fit_baseline_c_frac_5p > 0.1f) {
            float raw_signal = static_cast<float>(profile.t_freq_5prime[i]) - profile.fit_baseline_5prime;
            profile.damage_rate_5prime[i] = std::max(0.0f, raw_signal / fit_baseline_c_frac_5p);
        } else {
            profile.damage_rate_5prime[i] = 0.0f;
        }

        // 3' end: G→A damage using fit baseline
        if (fit_baseline_g_frac_3p > 0.1f) {
            float raw_signal = static_cast<float>(profile.a_freq_3prime[i]) - profile.fit_baseline_3prime;
            profile.damage_rate_3prime[i] = std::max(0.0f, raw_signal / fit_baseline_g_frac_3p);
        } else {
            profile.damage_rate_3prime[i] = 0.0f;
        }
    }

    // DEBUG: Print damage_rate[0] calculation
    fprintf(stderr, "  [DEBUG-DAMAGE] 5': t_freq[0]=%.4f, baseline=%.4f, c_frac=%.4f, rate[0]=%.4f\n",
            profile.t_freq_5prime[0], profile.fit_baseline_5prime, fit_baseline_c_frac_5p, profile.damage_rate_5prime[0]);
    fprintf(stderr, "  [DEBUG-DAMAGE] 3': a_freq[0]=%.4f, baseline=%.4f, g_frac=%.4f, rate[0]=%.4f\n",
            profile.a_freq_3prime[0], profile.fit_baseline_3prime, fit_baseline_g_frac_3p, profile.damage_rate_3prime[0]);

    // =========================================================================
    // Inverted pattern detection
    // Detect when terminal T/(T+C) < baseline from middle of reads (opposite of damage pattern)
    // This indicates reference-free detection failure due to:
    // - AT-rich organisms with terminal artifacts
    // - Adapter contamination
    // - Quality trimming bias
    //
    // NOTE: We compare terminal position 0 against TRUE baseline (middle of reads),
    // NOT against positions 10-14 which are still in the terminal region and may
    // share composition bias with position 0.
    // =========================================================================
    {
        // Compute terminal gradients for reporting (but DON'T set inverted patterns here)
        // Position 0 has artifacts - we use hexamer-based detection (positions 1-6) instead.
        // The hexamer detection at the end of this function sets inverted_pattern_5prime.

        // 5' end: terminal T/(T+C) - baseline (for reporting only)
        double terminal_tc_5 = profile.t_freq_5prime[0];  // Already normalized
        profile.terminal_gradient_5prime = static_cast<float>(terminal_tc_5 - baseline_tc);

        // 3' end: terminal A/(A+G) - baseline (for reporting only)
        double terminal_ag_3 = profile.a_freq_3prime[0];  // Already normalized
        profile.terminal_gradient_3prime = static_cast<float>(terminal_ag_3 - baseline_ag);

        // NOTE: inverted_pattern flags are set later by hexamer-based detection,
        // which uses positions 1-6 and is more reliable than position 0.
    }

    // Compute codon-position-aware damage rates
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

    // Compute codon-position-specific damage for JSON output (using fit baseline)
    {
        float pos1_t_rate = profile.codon_pos_t_rate_5prime[0];
        float pos2_t_rate = profile.codon_pos_t_rate_5prime[1];
        float pos3_t_rate = profile.codon_pos_t_rate_5prime[2];  // Wobble position

        float baseline_tc_f = profile.fit_baseline_5prime;
        float baseline_c_f = 1.0f - profile.fit_baseline_5prime;

        if (baseline_c_f > 0.1f) {
            profile.codon_pos1_damage = std::max(0.0f, (pos1_t_rate - baseline_tc_f) / baseline_c_f);
            profile.codon_pos2_damage = std::max(0.0f, (pos2_t_rate - baseline_tc_f) / baseline_c_f);
            profile.codon_pos3_damage = std::max(0.0f, (pos3_t_rate - baseline_tc_f) / baseline_c_f);
        } else {
            profile.codon_pos1_damage = 0.0f;
            profile.codon_pos2_damage = 0.0f;
            profile.codon_pos3_damage = 0.0f;
        }

        // Calculate wobble ratio: pos3 / ((pos1 + pos2) / 2)
        float avg_pos12_damage = (profile.codon_pos1_damage + profile.codon_pos2_damage) / 2.0f;
        if (avg_pos12_damage > 0.005f) {
            profile.wobble_ratio = profile.codon_pos3_damage / avg_pos12_damage;
        } else if (profile.codon_pos3_damage > 0.005f) {
            profile.wobble_ratio = 2.0f;  // Cap at 2x
        } else {
            profile.wobble_ratio = 1.0f;
        }

        profile.wobble_ratio = std::clamp(profile.wobble_ratio, 0.5f, 3.0f);
    }

    // Compute CpG damage rates
    size_t cpg_total = profile.cpg_c_count + profile.cpg_t_count;
    if (cpg_total > 10) {
        profile.cpg_damage_rate = static_cast<float>(profile.cpg_t_count) / cpg_total;
    }

    size_t non_cpg_total = profile.non_cpg_c_count + profile.non_cpg_t_count;
    if (non_cpg_total > 10) {
        profile.non_cpg_damage_rate = static_cast<float>(profile.non_cpg_t_count) / non_cpg_total;
    }

    // Summary statistics
    profile.max_damage_5prime = profile.damage_rate_5prime[0];
    profile.max_damage_3prime = profile.damage_rate_3prime[0];

    // Briggs-like damage model parameter estimation (closed-form)
    auto estimate_briggs_params = [](const std::array<float, 15>& rates, float max_rate,
                                     float& delta_s, float& delta_d, float& lambda, float& r_squared) {
        delta_s = 0.0f;
        delta_d = 0.0f;
        lambda = 0.3f;
        r_squared = 0.0f;

        if (max_rate < 0.02f) return;

        delta_s = rates[0];

        float sum_background = 0.0f;
        for (int i = 10; i < 15; ++i) {
            sum_background += rates[i];
        }
        delta_d = sum_background / 5.0f;

        if (delta_s <= delta_d + 0.01f) {
            lambda = 0.3f;
            return;
        }

        // Weighted log-linear regression for lambda estimation
        float sum_xy = 0.0f;
        float sum_x2 = 0.0f;
        float sum_y = 0.0f;
        float sum_x = 0.0f;
        float sum_w = 0.0f;
        float sum_y2 = 0.0f;

        for (int pos = 0; pos < 15; ++pos) {
            float normalized = (rates[pos] - delta_d) / (delta_s - delta_d + 0.001f);
            normalized = std::clamp(normalized, 0.01f, 1.0f);

            float y = std::log(normalized);
            float x = static_cast<float>(pos);
            float w = normalized * normalized;

            sum_xy += w * x * y;
            sum_x2 += w * x * x;
            sum_y += w * y;
            sum_x += w * x;
            sum_w += w;
            sum_y2 += w * y * y;
        }

        float denom = sum_w * sum_x2 - sum_x * sum_x;
        if (std::abs(denom) < 0.001f) {
            lambda = 0.3f;
            return;
        }

        float slope = (sum_w * sum_xy - sum_x * sum_y) / denom;

        if (slope >= 0.0f) {
            lambda = 0.3f;
            return;
        }

        lambda = 1.0f - std::exp(slope);
        lambda = std::clamp(lambda, 0.05f, 0.95f);

        // Calculate R²
        float ss_tot = sum_y2 - (sum_y * sum_y) / sum_w;
        float intercept = (sum_y - slope * sum_x) / sum_w;
        float ss_res = 0.0f;
        for (int pos = 0; pos < 15; ++pos) {
            float normalized = (rates[pos] - delta_d) / (delta_s - delta_d + 0.001f);
            normalized = std::clamp(normalized, 0.01f, 1.0f);
            float y = std::log(normalized);
            float y_pred = slope * static_cast<float>(pos) + intercept;
            float w = normalized * normalized;
            ss_res += w * (y - y_pred) * (y - y_pred);
        }
        r_squared = (ss_tot > 0.001f) ? std::max(0.0f, 1.0f - ss_res / ss_tot) : 0.0f;
    };

    estimate_briggs_params(profile.damage_rate_5prime, profile.max_damage_5prime,
                           profile.delta_s_5prime, profile.delta_d_5prime,
                           profile.lambda_5prime, profile.r_squared_5prime);

    estimate_briggs_params(profile.damage_rate_3prime, profile.max_damage_3prime,
                           profile.delta_s_3prime, profile.delta_d_3prime,
                           profile.lambda_3prime, profile.r_squared_3prime);

    profile.lambda_5prime = std::clamp(profile.lambda_5prime, 0.1f, 1.0f);
    profile.lambda_3prime = std::clamp(profile.lambda_3prime, 0.1f, 1.0f);

    // Library type handling: default to double-stranded unless user forces single-stranded
    if (profile.forced_library_type != SampleDamageProfile::LibraryType::UNKNOWN) {
        profile.library_type = profile.forced_library_type;
    } else {
        profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
    }

    // Flag inversions: terminals depleted relative to interior with statistical support
    const bool inversion_5 = (profile.terminal_shift_5prime < -0.01f && profile.terminal_z_5prime < -2.0f);
    const bool inversion_3 = (profile.terminal_shift_3prime < -0.01f && profile.terminal_z_3prime < -2.0f);
    profile.terminal_inversion = inversion_5 || inversion_3;

    // =========================================================================
    // Hexamer-based damage detection (reference-independent)
    // Compare hexamer frequencies at 5' terminal vs interior positions
    // =========================================================================
    fprintf(stderr, "  [DEBUG] Hexamer collection (pos 1-6): 5' count=%zu, interior count=%zu\n",
            profile.n_hexamers_5prime, profile.n_hexamers_interior);

    if (profile.n_hexamers_5prime >= 1000 && profile.n_hexamers_interior >= 1000) {
        double llr_sum = 0.0;
        double weight_sum = 0.0;

        // For each pair of hexamers differing only at position 0 (C vs T):
        for (uint32_t base_code = 0; base_code < 1024; ++base_code) {
            uint32_t c_hex = 0x400 | base_code;  // C at position 0
            uint32_t t_hex = 0xC00 | base_code;  // T at position 0

            float expected_c = get_hexamer_freq(c_hex);
            float expected_t = get_hexamer_freq(t_hex);

            if (expected_c < 1e-6f || expected_t < 1e-6f) continue;

            double term_c = profile.hexamer_count_5prime[c_hex];
            double term_t = profile.hexamer_count_5prime[t_hex];
            double int_c = profile.hexamer_count_interior[c_hex];
            double int_t = profile.hexamer_count_interior[t_hex];

            double total_term = term_c + term_t;
            double total_int = int_c + int_t;
            if (total_term < 5 || total_int < 5) continue;

            double ratio_term = term_t / total_term;
            double ratio_int = int_t / total_int;
            double excess = ratio_term - ratio_int;

            double weight = (expected_c + expected_t) * std::sqrt(total_term * total_int);
            llr_sum += excess * weight;
            weight_sum += weight;
        }

        if (weight_sum > 0) {
            float raw_llr = static_cast<float>(llr_sum / weight_sum);

            // Store raw hexamer LLR without inversion correction
            // The raw value is the actual terminal-vs-interior hexamer composition difference
            // Positive = T-starting hexamers enriched at terminal (damage pattern)
            // Negative = C-starting hexamers enriched at terminal (AT-rich composition bias)
            profile.hexamer_damage_llr = raw_llr;

            if (raw_llr < -0.01f && profile.terminal_inversion) {
                fprintf(stderr, "  [DEBUG] Inverted hexamer pattern: raw_llr=%.4f (NOT using as damage proxy)\n",
                        raw_llr);
            }
        }

        // Compute hexamer-based T/(T+C) ratios (more reliable than position 0 or 1 alone)
        double total_term_t = 0, total_term_c = 0;
        double total_int_t = 0, total_int_c = 0;
        for (uint32_t base_code = 0; base_code < 1024; ++base_code) {
            uint32_t c_hex = 0x400 | base_code;
            uint32_t t_hex = 0xC00 | base_code;
            total_term_c += profile.hexamer_count_5prime[c_hex];
            total_term_t += profile.hexamer_count_5prime[t_hex];
            total_int_c += profile.hexamer_count_interior[c_hex];
            total_int_t += profile.hexamer_count_interior[t_hex];
        }
        double term_t_ratio = total_term_t / (total_term_t + total_term_c + 1e-10);
        double int_t_ratio = total_int_t / (total_int_t + total_int_c + 1e-10);

        // Store hexamer-based T/(T+C) ratios for use in amplitude calculation
        profile.hexamer_terminal_tc = static_cast<float>(term_t_ratio);
        profile.hexamer_interior_tc = static_cast<float>(int_t_ratio);
        profile.hexamer_excess_tc = static_cast<float>(term_t_ratio - int_t_ratio);

        fprintf(stderr, "  [DEBUG] C-starting hexamers: terminal=%.0f, interior=%.0f\n",
                total_term_c, total_int_c);
        fprintf(stderr, "  [DEBUG] T-starting hexamers: terminal=%.0f, interior=%.0f\n",
                total_term_t, total_int_t);
        fprintf(stderr, "  [DEBUG] T/(T+C) from hexamers: terminal=%.1f%%, interior=%.1f%%, excess=%.2f%%\n",
                term_t_ratio * 100, int_t_ratio * 100, (term_t_ratio - int_t_ratio) * 100);
        fprintf(stderr, "  [DEBUG] Hexamer damage LLR (weighted): %.4f\n", profile.hexamer_damage_llr);

        // CRITICAL: If hexamer analysis shows terminal T/(T+C) < interior,
        // this overrides position-0 based detection.
        // Position 0 often has sequencing artifacts (adapter ligation, first-base bias)
        // that don't represent real damage. Hexamer analysis (positions 1-6) is more reliable.
        if (profile.hexamer_damage_llr < -0.02f) {
            // Terminal T/(T+C) is LOWER than interior - no damage signal at 5' end
            profile.inverted_pattern_5prime = true;
            fprintf(stderr, "  [DEBUG] Hexamer inversion detected (LLR=%.4f < -0.02): setting inverted_pattern_5prime=true\n",
                    profile.hexamer_damage_llr);
        }
    }

    // =========================================================================
    // Composition bias detection using negative controls
    // =========================================================================
    // If the negative control (A/(A+G) at 5', T/(T+C) at 3') shows comparable
    // enrichment to the "damage" signal, it's likely composition bias, not damage.
    // Decision rule from GPT-5.2:
    //   Flag as bias if: |ctrl_shift| >= max(0.005, 0.5 * |damage_shift|)

    float damage_shift_5 = profile.terminal_shift_5prime;
    float ctrl_shift_5 = std::abs(profile.ctrl_shift_5prime);
    float threshold_5 = std::max(0.005f, 0.5f * std::abs(damage_shift_5));
    if (ctrl_shift_5 >= threshold_5 && damage_shift_5 > 0.01f) {
        profile.composition_bias_5prime = true;
        fprintf(stderr, "  [DEBUG] 5' control shift (%.3f) >= threshold (%.3f): composition_bias_5prime=true\n",
                ctrl_shift_5, threshold_5);
    }

    float damage_shift_3 = profile.terminal_shift_3prime;
    float ctrl_shift_3 = std::abs(profile.ctrl_shift_3prime);
    float threshold_3 = std::max(0.005f, 0.5f * std::abs(damage_shift_3));
    if (ctrl_shift_3 >= threshold_3 && damage_shift_3 > 0.01f) {
        profile.composition_bias_3prime = true;
        fprintf(stderr, "  [DEBUG] 3' control shift (%.3f) >= threshold (%.3f): composition_bias_3prime=true\n",
                ctrl_shift_3, threshold_3);
    }

    // Composition bias flags unreliable detection, but does NOT set inverted pattern.
    // Inverted pattern is reserved for hexamer-based inversion (terminal < interior).
    // Composition bias just means the signal is confounded - we can still estimate d_max
    // from the position-based analysis, but it's less reliable.
    if (profile.composition_bias_5prime) {
        fprintf(stderr, "  [DEBUG] Composition bias at 5': detection unreliable but estimate available\n");
    }
    if (profile.composition_bias_3prime) {
        fprintf(stderr, "  [DEBUG] Composition bias at 3': detection unreliable but estimate available\n");
    }

    // Sample classification
    float damage_signal = (profile.max_damage_5prime + profile.max_damage_3prime) / 2.0f;

    // Hexamer-based damage estimation
    // For non-inverted samples: use positive hexamer LLR as damage boost
    // For inverted samples: use z-score asymmetry to detect real damage vs composition bias
    //   - Real damage: 3' G→A signal should be relatively stronger (z_ratio < 1.0)
    //   - Composition bias: 5' T enrichment dominates (z_ratio > 1.5)
    float hexamer_boost = 0.0f;
    if (profile.hexamer_damage_llr > 0.02f && !profile.terminal_inversion) {
        // Normal sample with clear hexamer signal
        hexamer_boost = profile.hexamer_damage_llr * 8.0f;
        fprintf(stderr, "  [DEBUG] Using hexamer-based damage boost: %.1f%% (from LLR %.4f)\n",
                hexamer_boost * 100, profile.hexamer_damage_llr);
    } else if (profile.terminal_inversion) {
        // Inverted sample: check z-score asymmetry to distinguish real damage from composition bias
        float z5_abs = std::abs(profile.terminal_z_5prime);
        float z3_abs = std::abs(profile.terminal_z_3prime);
        float z_ratio = (z3_abs > 0) ? z5_abs / z3_abs : 10.0f;

        fprintf(stderr, "  [DEBUG] Inverted sample z-scores: 5'=%.1f, 3'=%.1f, ratio=%.2f\n",
                profile.terminal_z_5prime, profile.terminal_z_3prime, z_ratio);

        // Use absolute hexamer LLR magnitude for inverted samples
        float abs_llr = std::abs(profile.hexamer_damage_llr);

        if (z_ratio < 1.2f && abs_llr > 0.02f) {
            // 3' signal is relatively strong → likely real G→A damage
            // Use hexamer estimate but with more conservative scaling
            hexamer_boost = abs_llr * 8.0f;
            fprintf(stderr, "  [DEBUG] Inverted + asymmetric (z_ratio=%.2f < 1.2): likely REAL damage, hexamer boost=%.1f%%\n",
                    z_ratio, hexamer_boost * 100);
        } else if (z_ratio > 1.5f) {
            // 5' signal dominates → likely composition bias, NOT damage
            fprintf(stderr, "  [DEBUG] Inverted + 5'-dominant (z_ratio=%.2f > 1.5): likely COMPOSITION BIAS, no hexamer estimate\n",
                    z_ratio);
        } else {
            // Ambiguous case: use conservative estimate
            hexamer_boost = abs_llr * 4.0f;  // Half the normal scaling
            fprintf(stderr, "  [DEBUG] Inverted + ambiguous (z_ratio=%.2f): conservative hexamer boost=%.1f%%\n",
                    z_ratio, hexamer_boost * 100);
        }

        // Apply hexamer boost for inverted samples
        if (hexamer_boost > 0.01f && damage_signal < 0.01f) {
            damage_signal = hexamer_boost;
            profile.max_damage_5prime = hexamer_boost;
            profile.max_damage_3prime = hexamer_boost;
        }
    }

    float cpg_boost = 0.0f;
    float wobble_boost = 0.0f;

    if (damage_signal > 0.01f) {
        if (profile.cpg_damage_rate > profile.non_cpg_damage_rate + 0.05f) {
            cpg_boost = 0.03f;
        }

        float wobble_enrichment = profile.codon_pos_t_rate_5prime[2] -
                                  (profile.codon_pos_t_rate_5prime[0] + profile.codon_pos_t_rate_5prime[1]) / 2.0f;
        if (wobble_enrichment > 0.05f) {
            wobble_boost = 0.02f;
        }
    }

    float total_signal = damage_signal + cpg_boost + wobble_boost + hexamer_boost;

    if (total_signal > 0.12f) {
        profile.sample_damage_prob = 0.95f;
    } else if (total_signal > 0.06f) {
        profile.sample_damage_prob = 0.80f;
    } else if (total_signal > 0.03f) {
        profile.sample_damage_prob = 0.50f;
    } else {
        profile.sample_damage_prob = 0.20f;
    }

    // =========================================================================
    // D_max estimation - JOINT EVIDENCE from Channel A + Channel B
    //
    // NEW APPROACH: Use Channel B (stop conversion) as independent validator
    // - If Channel A fires AND Channel B fires: real damage → report d_max
    // - If Channel A fires BUT Channel B is flat: compositional artifact → d_max = 0
    // - If neither fires: no damage → d_max = 0
    //
    // This solves the fundamental limitation of reference-free detection:
    // T/(T+C) elevation can be from composition OR damage, but stop codons
    // appearing in damage-susceptible contexts can ONLY be from real C→T damage.
    // =========================================================================
    {
        // First compute raw d_max values from Channel A (nucleotide frequencies)
        float raw_d_max_5prime = profile.damage_rate_5prime[0];
        float raw_d_max_3prime = profile.damage_rate_3prime[0];

        // Clamp to valid range [0, 1]
        raw_d_max_5prime = std::clamp(raw_d_max_5prime, 0.0f, 1.0f);
        raw_d_max_3prime = std::clamp(raw_d_max_3prime, 0.0f, 1.0f);

        fprintf(stderr, "  [D_MAX] Raw from Channel A: 5'=%.2f%%, 3'=%.2f%%\n",
                raw_d_max_5prime * 100, raw_d_max_3prime * 100);

        // Compute asymmetry: |D_5p - D_3p| / ((D_5p + D_3p) / 2)
        float d_sum = raw_d_max_5prime + raw_d_max_3prime;
        if (d_sum > 0.01f) {
            profile.asymmetry = std::abs(raw_d_max_5prime - raw_d_max_3prime) / (d_sum / 2.0f);
        } else {
            profile.asymmetry = 0.0f;
        }
        profile.high_asymmetry = (profile.asymmetry > 0.5f);

        // =========================================================================
        // JOINT EVIDENCE DECISION FOR D_MAX
        // =========================================================================
        if (profile.damage_artifact) {
            // Channel A fires but Channel B identified it as compositional artifact
            // This is the KEY FIX: set d_max to 0 for false positives
            profile.d_max_5prime = 0.0f;
            profile.d_max_3prime = 0.0f;
            profile.d_max_combined = 0.0f;
            profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
            fprintf(stderr, "  [D_MAX] ARTIFACT DETECTED by Channel B: setting d_max = 0\n");
        } else if (profile.damage_validated) {
            // Both channels agree: use Channel A's raw estimates
            // BUT if raw values are 0 due to inverted pattern, use fit_amplitude instead
            // This handles the case where position 0 is below baseline but the
            // exponential fit and Channel B both confirm real damage
            float effective_d_max_5prime = raw_d_max_5prime;
            float effective_d_max_3prime = raw_d_max_3prime;

            // If raw estimate is 0 but fit amplitude shows decay, use amplitude
            if (raw_d_max_5prime < 0.001f && profile.fit_amplitude_5prime > 0.02f) {
                effective_d_max_5prime = profile.fit_amplitude_5prime;
                fprintf(stderr, "  [D_MAX] 5' raw=0 but fit_amplitude=%.1f%%, using amplitude\n",
                        profile.fit_amplitude_5prime * 100);
            }
            if (raw_d_max_3prime < 0.001f && profile.fit_amplitude_3prime > 0.02f) {
                effective_d_max_3prime = profile.fit_amplitude_3prime;
                fprintf(stderr, "  [D_MAX] 3' raw=0 but fit_amplitude=%.1f%%, using amplitude\n",
                        profile.fit_amplitude_3prime * 100);
            }

            profile.d_max_5prime = effective_d_max_5prime;
            profile.d_max_3prime = effective_d_max_3prime;

            // Combine using asymmetry-aware logic
            if (profile.high_asymmetry) {
                profile.d_max_combined = std::min(profile.d_max_5prime, profile.d_max_3prime);
                profile.d_max_source = SampleDamageProfile::DmaxSource::MIN_ASYMMETRY;
            } else {
                profile.d_max_combined = (profile.d_max_5prime + profile.d_max_3prime) / 2.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
            fprintf(stderr, "  [D_MAX] VALIDATED by Channel B: d_max = %.2f%%\n",
                    profile.d_max_combined * 100);
        } else if (profile.channel_b_valid && profile.stop_amplitude_5prime > 0.005f) {
            // Channel B has data and shows some stop excess
            // Use Channel B's stop amplitude as a LOWER BOUND for d_max
            // r(p) = r0 + (1-r0) * delta(p), so delta(p) = (r(p) - r0) / (1 - r0)
            // At p=0: delta_max = stop_amplitude / (1 - baseline)
            float channel_b_dmax = profile.stop_amplitude_5prime /
                                  (1.0f - profile.stop_conversion_rate_baseline + 0.01f);
            channel_b_dmax = std::clamp(channel_b_dmax, 0.0f, 1.0f);

            // Use the MINIMUM of Channel A and Channel B estimates
            // This is conservative: Channel A can overestimate due to composition
            profile.d_max_5prime = std::min(raw_d_max_5prime, channel_b_dmax * 2.0f);
            profile.d_max_3prime = std::min(raw_d_max_3prime, channel_b_dmax * 2.0f);

            if (profile.high_asymmetry) {
                profile.d_max_combined = std::min(profile.d_max_5prime, profile.d_max_3prime);
                profile.d_max_source = SampleDamageProfile::DmaxSource::MIN_ASYMMETRY;
            } else {
                profile.d_max_combined = (profile.d_max_5prime + profile.d_max_3prime) / 2.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
            fprintf(stderr, "  [D_MAX] Channel B estimate: %.2f%%, combined: %.2f%%\n",
                    channel_b_dmax * 100, profile.d_max_combined * 100);
        } else {
            // No strong evidence either way - use raw values but with caution
            profile.d_max_5prime = raw_d_max_5prime;
            profile.d_max_3prime = raw_d_max_3prime;

            if (profile.high_asymmetry) {
                profile.d_max_combined = std::min(profile.d_max_5prime, profile.d_max_3prime);
                profile.d_max_source = SampleDamageProfile::DmaxSource::MIN_ASYMMETRY;
            } else if (d_sum > 0.01f) {
                profile.d_max_combined = (profile.d_max_5prime + profile.d_max_3prime) / 2.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else {
                profile.d_max_combined = 0.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
            }
            fprintf(stderr, "  [D_MAX] No Channel B validation: d_max = %.2f%% (use with caution)\n",
                    profile.d_max_combined * 100);
        }

        // Additional reliability check: if either end has inverted pattern, use the other
        // BUT skip this check if Channel B validated damage - Channel B is the ground truth
        if (!profile.damage_validated) {
            if (profile.inverted_pattern_5prime && !profile.inverted_pattern_3prime) {
                profile.d_max_combined = profile.d_max_3prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::THREE_PRIME_ONLY;
            } else if (profile.inverted_pattern_3prime && !profile.inverted_pattern_5prime) {
                profile.d_max_combined = profile.d_max_5prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::FIVE_PRIME_ONLY;
            } else if (profile.inverted_pattern_5prime && profile.inverted_pattern_3prime) {
                // Both ends inverted - can't reliably estimate damage
                profile.d_max_combined = 0.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
            }
        }

        fprintf(stderr, "  [D_MAX] Final: d_max_5p=%.2f%%, d_max_3p=%.2f%%, combined=%.2f%% (source=%s)\n",
                profile.d_max_5prime * 100, profile.d_max_3prime * 100,
                profile.d_max_combined * 100, profile.d_max_source_str());
    }
}

void FrameSelector::merge_sample_profiles(SampleDamageProfile& dst, const SampleDamageProfile& src) {
    // Merge position-specific counts (before normalization)
    for (int i = 0; i < 15; ++i) {
        // Damage signal counts
        dst.t_freq_5prime[i] += src.t_freq_5prime[i];
        dst.c_freq_5prime[i] += src.c_freq_5prime[i];
        dst.a_freq_3prime[i] += src.a_freq_3prime[i];
        dst.g_freq_3prime[i] += src.g_freq_3prime[i];
        dst.tc_total_5prime[i] += src.tc_total_5prime[i];
        dst.ag_total_3prime[i] += src.ag_total_3prime[i];
        // Negative control counts
        dst.a_freq_5prime[i] += src.a_freq_5prime[i];
        dst.g_freq_5prime[i] += src.g_freq_5prime[i];
        dst.t_freq_3prime[i] += src.t_freq_3prime[i];
        dst.c_freq_3prime[i] += src.c_freq_3prime[i];
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

    // Merge hexamer counts
    for (uint32_t i = 0; i < 4096; ++i) {
        dst.hexamer_count_5prime[i] += src.hexamer_count_5prime[i];
        dst.hexamer_count_interior[i] += src.hexamer_count_interior[i];
    }
    dst.n_hexamers_5prime += src.n_hexamers_5prime;
    dst.n_hexamers_interior += src.n_hexamers_interior;

    // Merge Channel B convertible codon counts
    for (int i = 0; i < 15; ++i) {
        dst.convertible_caa_5prime[i] += src.convertible_caa_5prime[i];
        dst.convertible_taa_5prime[i] += src.convertible_taa_5prime[i];
        dst.convertible_cag_5prime[i] += src.convertible_cag_5prime[i];
        dst.convertible_tag_5prime[i] += src.convertible_tag_5prime[i];
        dst.convertible_cga_5prime[i] += src.convertible_cga_5prime[i];
        dst.convertible_tga_5prime[i] += src.convertible_tga_5prime[i];
        dst.total_codons_5prime[i] += src.total_codons_5prime[i];
    }
    dst.convertible_caa_interior += src.convertible_caa_interior;
    dst.convertible_taa_interior += src.convertible_taa_interior;
    dst.convertible_cag_interior += src.convertible_cag_interior;
    dst.convertible_tag_interior += src.convertible_tag_interior;
    dst.convertible_cga_interior += src.convertible_cga_interior;
    dst.convertible_tga_interior += src.convertible_tga_interior;
    dst.total_codons_interior += src.total_codons_interior;

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
        if (base == 'T') {
            profile.t_freq_5prime[i] += weight;
            profile.tc_total_5prime[i] += weight;
        } else if (base == 'C') {
            profile.c_freq_5prime[i] += weight;
            profile.tc_total_5prime[i] += weight;
        }
    }

    // Count bases at 3' end positions (weighted)
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        size_t pos = len - 1 - i;
        char base = fast_upper(seq[pos]);
        if (base == 'A') {
            profile.a_freq_3prime[i] += weight;
            profile.ag_total_3prime[i] += weight;
        } else if (base == 'G') {
            profile.g_freq_3prime[i] += weight;
            profile.ag_total_3prime[i] += weight;
        }
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
        profile.tc_total_5prime[i] = 0.0;
        profile.ag_total_3prime[i] = 0.0;
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
    profile.terminal_shift_5prime = 0.0f;
    profile.terminal_shift_3prime = 0.0f;
    profile.terminal_z_5prime = 0.0f;
    profile.terminal_z_3prime = 0.0f;
    profile.terminal_inversion = false;
    // Default to double-stranded unless user forces single-stranded
    profile.library_type = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;

    // Reset hexamer counts
    profile.hexamer_count_5prime.fill(0.0);
    profile.hexamer_count_interior.fill(0.0);
    profile.n_hexamers_5prime = 0;
    profile.n_hexamers_interior = 0;
    profile.hexamer_damage_llr = 0.0f;

    // Reset Channel B convertible codon counts
    profile.convertible_caa_5prime.fill(0.0);
    profile.convertible_taa_5prime.fill(0.0);
    profile.convertible_cag_5prime.fill(0.0);
    profile.convertible_tag_5prime.fill(0.0);
    profile.convertible_cga_5prime.fill(0.0);
    profile.convertible_tga_5prime.fill(0.0);
    profile.total_codons_5prime.fill(0.0);
    profile.convertible_caa_interior = 0.0;
    profile.convertible_taa_interior = 0.0;
    profile.convertible_cag_interior = 0.0;
    profile.convertible_tag_interior = 0.0;
    profile.convertible_cga_interior = 0.0;
    profile.convertible_tga_interior = 0.0;
    profile.total_codons_interior = 0.0;
    profile.stop_conversion_rate_baseline = 0.0f;
    profile.stop_decay_llr_5prime = 0.0f;
    profile.stop_amplitude_5prime = 0.0f;
    profile.channel_b_valid = false;
    profile.damage_validated = false;
    profile.damage_artifact = false;

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

} // namespace agp
