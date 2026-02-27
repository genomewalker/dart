// Sample-level damage profile management

#include "dart/frame_selector.hpp"
#include "dart/codon_tables.hpp"
#include "dart/hexamer_tables.hpp"
#include <algorithm>
#include <cmath>
#include <array>

namespace dart {

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

    // Use external baseline (middle-of-read) if provided, otherwise estimate from positions 10-14
    // CRITICAL: Middle-of-read baseline is more reliable because positions 10-14 may still
    // have read-end composition artifacts that inflate the apparent "damage" signal.
    //
    // IMPORTANT: Even when coverage is low (e.g. small benchmarks), we still need a sane
    // baseline for damage-rate computation. Returning b=0 here inflates damage to the raw
    // T/(T+C) ratio and can cascade into 0-gene failure modes.
    float b;
    if (external_baseline >= 0.0f) {
        b = std::clamp(external_baseline, 0.001f, 0.999f);
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
        b = std::clamp(b, 0.001f, 0.999f);
    }

    // Count valid positions
    int n_valid = 0;
    for (int i = 0; i < 15; ++i) {
        if (coverage[i] >= MIN_COVERAGE) n_valid++;
    }

    // CRITICAL: Estimate amplitude from position 1, NOT position 0
    // Position 0 is uniquely prone to first-cycle/ligation/trimming artifacts
    // that can arbitrarily bias the estimate. Position 1 is more reliable.
    float A = (coverage[1] >= MIN_COVERAGE) ?
              std::max(0.0f, static_cast<float>(freq[1]) - b) : 0.0f;

    // If we don't have enough positions to fit a decay curve, return baseline + a
    // conservative amplitude estimate and keep the initial lambda.
    if (n_valid < 5) {
        return {b, A, lambda_init, 1.0f};
    }

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
        // Guard: len >= 63 required to have valid interior region [30, len-30)
        // Without this, len - 30 underflows for short reads causing OOB access
        constexpr size_t INTERIOR_MIN_LEN = 63;  // 30 + 3 + 30
        if (len >= INTERIOR_MIN_LEN) {
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
    }

    // =========================================================================
    // CHANNEL C: Oxidative stop codon tracking (G→T transversions)
    // Track GAG→TAG, GAA→TAA, GGA→TGA pairs by nucleotide position
    // Unlike deamination, oxidative damage is UNIFORM across read length
    // =========================================================================
    if (len >= 18) {
        // Track oxidative convertible codons at 5' end
        for (int frame = 0; frame < 3; ++frame) {
            for (size_t k = 0; ; ++k) {
                size_t codon_start = frame + 3 * k;
                if (codon_start + 3 > len || codon_start > 14) break;

                size_t p = codon_start;
                if (p >= 15) break;

                char b0 = fast_upper(seq[codon_start]);
                char b1 = fast_upper(seq[codon_start + 1]);
                char b2 = fast_upper(seq[codon_start + 2]);

                if ((b0 != 'A' && b0 != 'C' && b0 != 'G' && b0 != 'T') ||
                    (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') ||
                    (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T')) {
                    continue;
                }

                // GAG (Glu) → TAG (Stop) via G→T at position 0
                if (b1 == 'A' && b2 == 'G') {
                    if (b0 == 'G') profile.convertible_gag_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tag_ox_5prime[p]++;
                }
                // GAA (Glu) → TAA (Stop) via G→T at position 0
                if (b1 == 'A' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gaa_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_taa_ox_5prime[p]++;
                }
                // GGA (Gly) → TGA (Stop) via G→T at position 0
                if (b1 == 'G' && b2 == 'A') {
                    if (b0 == 'G') profile.convertible_gga_5prime[p]++;
                    else if (b0 == 'T') profile.convertible_tga_ox_5prime[p]++;
                }
            }
        }

        // Track interior oxidative convertible codons (positions 30+)
        constexpr size_t INTERIOR_MIN_LEN_OX = 63;
        if (len >= INTERIOR_MIN_LEN_OX) {
            for (int frame = 0; frame < 3; ++frame) {
                for (size_t k = 0; ; ++k) {
                    size_t codon_start = frame + 3 * k;
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

                    // Track GAG, GAA, GGA (oxidation pre-images) and TAG, TAA, TGA (stops)
                    if (b1 == 'A' && b2 == 'G') {
                        if (b0 == 'G') profile.convertible_gag_interior++;
                        else if (b0 == 'T') profile.convertible_tag_ox_interior++;
                    }
                    if (b1 == 'A' && b2 == 'A') {
                        if (b0 == 'G') profile.convertible_gaa_interior++;
                        else if (b0 == 'T') profile.convertible_taa_ox_interior++;
                    }
                    if (b1 == 'G' && b2 == 'A') {
                        if (b0 == 'G') profile.convertible_gga_interior++;
                        else if (b0 == 'T') profile.convertible_tga_ox_interior++;
                    }
                }
            }
        }
    }

    // =========================================================================
    // CHANNEL D: G→T / C→A transversion tracking (oxidative damage)
    // Track G and T counts at each position for G→T rate calculation
    // Also track asymmetry: G→T vs T→G (real oxidation shows G→T > T→G)
    // =========================================================================
    // Count G and T at 5' end positions for G→T tracking
    for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
        char base = fast_upper(seq[i]);
        if (base == 'G') {
            profile.g_count_5prime[i]++;
        }
        // Note: T counts for oxidation tracking come from T where we'd expect G
        // This requires reference alignment, so for reference-free we use asymmetry instead
    }

    // =========================================================================
    // CHANNEL E: Depurination detection (purine enrichment at termini)
    // Depurination creates strand breaks preferentially at purine sites
    // Detection: terminal purine (A+G) rate vs interior purine rate
    // =========================================================================
    // Purine tracking is already done via a_freq/g_freq arrays
    // We'll compute purine enrichment in finalize_sample_profile()

    // =========================================================================
    // GC-STRATIFIED DAMAGE ACCUMULATION
    // Bin reads by interior GC content (positions 5+ to avoid terminal damage)
    // This handles metagenome heterogeneity where different organisms have
    // different GC content and damage levels
    // =========================================================================
    if (len >= 30) {
        // Compute interior GC from positions 5 to end-5 (avoid both terminal regions)
        size_t interior_start = 5;
        size_t interior_end = len > 10 ? len - 5 : len;
        uint64_t gc_count = 0;
        uint64_t at_count = 0;
        for (size_t i = interior_start; i < interior_end; ++i) {
            char base = fast_upper(seq[i]);
            if (base == 'G' || base == 'C') gc_count++;
            else if (base == 'A' || base == 'T') at_count++;
        }

        if (gc_count + at_count > 0) {
            float gc_frac = static_cast<float>(gc_count) / (gc_count + at_count);
            int bin_idx = std::min(static_cast<int>(gc_frac * SampleDamageProfile::N_GC_BINS),
                                   SampleDamageProfile::N_GC_BINS - 1);
            auto& bin = profile.gc_bins[bin_idx];

            // Accumulate Channel A counts (T and C at terminal positions)
            for (size_t i = 0; i < std::min(size_t(15), len); ++i) {
                char base = fast_upper(seq[i]);
                if (base == 'T') bin.t_counts[i]++;
                else if (base == 'C') bin.c_counts[i]++;
            }

            // Accumulate Channel A interior baseline (middle third)
            size_t mid_start = len / 3;
            size_t mid_end = 2 * len / 3;
            for (size_t i = mid_start; i < mid_end; ++i) {
                char base = fast_upper(seq[i]);
                if (base == 'T') bin.t_interior++;
                else if (base == 'C') bin.c_interior++;
            }

            // Accumulate Channel B counts (stop codons at terminal positions)
            for (int frame = 0; frame < 3; ++frame) {
                for (size_t k = 0; ; ++k) {
                    size_t codon_start = frame + 3 * k;
                    if (codon_start + 3 > len || codon_start > 14) break;
                    size_t p = codon_start;

                    char b0 = fast_upper(seq[codon_start]);
                    char b1 = fast_upper(seq[codon_start + 1]);
                    char b2 = fast_upper(seq[codon_start + 2]);

                    // CAA/TAA, CAG/TAG, CGA/TGA
                    if (b1 == 'A' && b2 == 'A') {
                        if (b0 == 'C') bin.pre_counts[p]++;
                        else if (b0 == 'T') bin.stop_counts[p]++;
                    }
                    if (b1 == 'A' && b2 == 'G') {
                        if (b0 == 'C') bin.pre_counts[p]++;
                        else if (b0 == 'T') bin.stop_counts[p]++;
                    }
                    if (b1 == 'G' && b2 == 'A') {
                        if (b0 == 'C') bin.pre_counts[p]++;
                        else if (b0 == 'T') bin.stop_counts[p]++;
                    }
                }
            }

            // Accumulate Channel B interior baseline
            // Guard: len >= 63 required to prevent unsigned underflow in len - 30
            if (len >= 63) {
                for (int frame = 0; frame < 3; ++frame) {
                    for (size_t k = 0; ; ++k) {
                        size_t codon_start = frame + 3 * k;
                        if (codon_start < 30 || codon_start + 3 > len - 30) {
                            if (codon_start + 3 > len - 30) break;
                            continue;
                        }

                        char b0 = fast_upper(seq[codon_start]);
                        char b1 = fast_upper(seq[codon_start + 1]);
                        char b2 = fast_upper(seq[codon_start + 2]);

                        if (b1 == 'A' && b2 == 'A') {
                            if (b0 == 'C') bin.pre_interior++;
                            else if (b0 == 'T') bin.stop_interior++;
                        }
                        if (b1 == 'A' && b2 == 'G') {
                            if (b0 == 'C') bin.pre_interior++;
                            else if (b0 == 'T') bin.stop_interior++;
                        }
                        if (b1 == 'G' && b2 == 'A') {
                            if (b0 == 'C') bin.pre_interior++;
                            else if (b0 == 'T') bin.stop_interior++;
                        }
                    }
                }
            }

            bin.n_reads++;
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

    // Terminal enrichment/depletion statistics
    // Check for position-0 artifact: if pos0 shows depletion but pos1 shows enrichment,
    // use pos1 instead (common adapter ligation bias pattern)
    auto stats_5_pos0 = compute_terminal_stats(profile.t_freq_5prime[0], profile.c_freq_5prime[0],
                                               profile.baseline_t_freq, profile.baseline_c_freq);
    auto stats_5_pos1 = compute_terminal_stats(profile.t_freq_5prime[1], profile.c_freq_5prime[1],
                                               profile.baseline_t_freq, profile.baseline_c_freq);

    // Detect position-0 artifact: two criteria
    // 1. Classic: pos0 depleted but pos1 enriched (relative to baseline)
    // 2. Jump pattern: Large pos0-to-pos1 difference (>3%) with pos1 enriched (>2%)
    //    This catches cases where baseline elevation masks the pos0 depletion
    bool pos0_artifact_5 = false;
    if (stats_5_pos0.first < -0.005f && stats_5_pos1.first > 0.01f) {
        // Classic pattern: pos0 below baseline, pos1 above
        pos0_artifact_5 = true;
    } else if (stats_5_pos1.first > 0.02f && (stats_5_pos1.first - stats_5_pos0.first) > 0.03f) {
        // Large pos0-to-pos1 jump pattern: pos1 clearly elevated but pos0 isn't
        // This indicates adapter artifact masking damage at pos0
        pos0_artifact_5 = true;
    }

    if (pos0_artifact_5) {
        profile.terminal_shift_5prime = stats_5_pos1.first;
        profile.terminal_z_5prime = stats_5_pos1.second;
        profile.position_0_artifact_5prime = true;
    } else {
        profile.terminal_shift_5prime = stats_5_pos0.first;
        profile.terminal_z_5prime = stats_5_pos0.second;
        profile.position_0_artifact_5prime = false;
    }

    auto stats_3_pos0 = compute_terminal_stats(profile.a_freq_3prime[0], profile.g_freq_3prime[0],
                                               profile.baseline_a_freq, profile.baseline_g_freq);
    auto stats_3_pos1 = compute_terminal_stats(profile.a_freq_3prime[1], profile.g_freq_3prime[1],
                                               profile.baseline_a_freq, profile.baseline_g_freq);

    // Same check for 3' end
    bool pos0_artifact_3 = false;
    if (stats_3_pos0.first < -0.005f && stats_3_pos1.first > 0.01f) {
        pos0_artifact_3 = true;
    } else if (stats_3_pos1.first > 0.02f && (stats_3_pos1.first - stats_3_pos0.first) > 0.03f) {
        pos0_artifact_3 = true;
    }

    if (pos0_artifact_3) {
        profile.terminal_shift_3prime = stats_3_pos1.first;
        profile.terminal_z_3prime = stats_3_pos1.second;
        profile.position_0_artifact_3prime = true;
    } else {
        profile.terminal_shift_3prime = stats_3_pos0.first;
        profile.terminal_z_3prime = stats_3_pos0.second;
        profile.position_0_artifact_3prime = false;
    }

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

    // =========================================================================
    // JOINT PROBABILISTIC MODEL
    // Fit unified model BEFORE normalization (need raw counts)
    // =========================================================================
    {
        JointDamageSuffStats jstats;

        // Channel A: T/(T+C) at 5' end (raw counts still in t_freq_5prime, c_freq_5prime)
        for (int p = 0; p < 15; ++p) {
            jstats.k_tc[p] = static_cast<uint64_t>(profile.t_freq_5prime[p]);
            jstats.n_tc[p] = static_cast<uint64_t>(profile.t_freq_5prime[p] + profile.c_freq_5prime[p]);
        }

        // Control: A/(A+G) at 5' end
        for (int p = 0; p < 15; ++p) {
            jstats.k_ag[p] = static_cast<uint64_t>(profile.a_freq_5prime[p]);
            jstats.n_ag[p] = static_cast<uint64_t>(profile.a_freq_5prime[p] + profile.g_freq_5prime[p]);
        }

        // Channel B: stop codon conversions (TAA+TAG+TGA) / (CAA+CAG+CGA + TAA+TAG+TGA)
        for (int p = 0; p < 15; ++p) {
            double stops = profile.convertible_taa_5prime[p] +
                          profile.convertible_tag_5prime[p] +
                          profile.convertible_tga_5prime[p];
            double pre = profile.convertible_caa_5prime[p] +
                        profile.convertible_cag_5prime[p] +
                        profile.convertible_cga_5prime[p];
            jstats.k_stop[p] = static_cast<uint64_t>(stops);
            jstats.n_stop[p] = static_cast<uint64_t>(pre + stops);
        }

        // Interior baselines (from middle of reads)
        jstats.k_tc_interior = static_cast<uint64_t>(base_tc_total * baseline_tc);
        jstats.n_tc_interior = static_cast<uint64_t>(base_tc_total);
        jstats.k_ag_interior = static_cast<uint64_t>(base_ag_total * baseline_ag);
        jstats.n_ag_interior = static_cast<uint64_t>(base_ag_total);

        // Channel B interior
        double stop_int = profile.convertible_taa_interior +
                         profile.convertible_tag_interior +
                         profile.convertible_tga_interior;
        double pre_int = profile.convertible_caa_interior +
                        profile.convertible_cag_interior +
                        profile.convertible_cga_interior;
        jstats.k_stop_interior = static_cast<uint64_t>(stop_int);
        jstats.n_stop_interior = static_cast<uint64_t>(pre_int + stop_int);

        // Fit the joint model
        JointDamageResult jresult = JointDamageModel::fit(jstats);

        // Store results in profile
        profile.joint_delta_max = jresult.delta_max;
        profile.joint_lambda = jresult.lambda;
        profile.joint_a_max = jresult.a_max;
        profile.joint_log_lik_m1 = jresult.log_lik_m1;
        profile.joint_log_lik_m0 = jresult.log_lik_m0;
        profile.joint_delta_bic = jresult.delta_bic;
        profile.joint_bayes_factor = jresult.bayes_factor;
        profile.joint_p_damage = jresult.p_damage;
        profile.joint_model_valid = jresult.valid;

        // Use joint model for damage decision
        if (jresult.valid) {
            // Check for strong evidence of damage
            // P(damage) > 0.95 OR ΔBIC > 10 (handles overflow when BF=inf)
            bool strong_damage_evidence = jresult.p_damage > 0.95f ||
                                          jresult.delta_bic > 10.0f ||
                                          std::isinf(jresult.bayes_factor);

            if (strong_damage_evidence) {
                profile.damage_validated = true;
                profile.damage_artifact = false;
            } else if (jresult.p_damage < 0.05f && jresult.a_max > 0.02f) {
                // Strong evidence for artifact (composition bias detected)
                profile.damage_validated = false;
                profile.damage_artifact = true;
            } else if (jresult.delta_bic < -10.0f) {
                // Strong evidence against damage (M0 fits much better)
                profile.damage_validated = false;
                profile.damage_artifact = (jresult.a_max > 0.02f);
            } else {
                // Uncertain - could go either way
                profile.damage_validated = false;
                profile.damage_artifact = false;
            }
        }
    }

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

        // Compute LOCAL baseline from positions 5-14 (same reads, past damage zone)
        // This avoids bias from interior baseline which comes from longer reads only
        double local_pre = 0.0, local_stop = 0.0;
        for (int p = 5; p < 15; ++p) {
            local_pre += profile.convertible_caa_5prime[p] +
                        profile.convertible_cag_5prime[p] +
                        profile.convertible_cga_5prime[p];
            local_stop += profile.convertible_taa_5prime[p] +
                         profile.convertible_tag_5prime[p] +
                         profile.convertible_tga_5prime[p];
        }
        float local_baseline = (local_pre + local_stop > 100)
            ? static_cast<float>(local_stop / (local_pre + local_stop))
            : profile.stop_conversion_rate_baseline;

        // Compute stop decay LLR using same lambda from Channel A
        // If both channels have same decay shape, it confirms real damage
        if (profile.channel_b_valid) {
            // Use local baseline instead of interior baseline
            float baseline_b = local_baseline;
            float lambda_b = fit_lambda_5p;  // Use fit lambda from Channel A (NOT profile.lambda_5prime which isn't set yet)

            // Estimate amplitude from positions 0-4
            // Include position 0 because for steep damage decay (AT-rich samples),
            // the signal is concentrated at position 0 and excluding it causes false negatives
            double sum_excess = 0.0, sum_weight = 0.0;
            for (int i = 0; i < 5; ++i) {
                if (stop_exposure[i] > 50) {
                    double weight = std::exp(-lambda_b * i);
                    double excess = stop_rate[i] - baseline_b;
                    sum_excess += stop_exposure[i] * excess / weight;
                    sum_weight += stop_exposure[i];
                }
            }
            float amplitude_b = (sum_weight > 0) ? static_cast<float>(sum_excess / sum_weight) : 0.0f;
            profile.stop_amplitude_5prime = std::max(0.0f, amplitude_b);

            // Compute log-likelihood ratio for exponential decay vs constant
            // Include position 0 for same reason as amplitude fitting above
            double ll_exp = 0.0, ll_const = 0.0;
            for (int p = 0; p < 10; ++p) {
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

            // =========================================================================
            // CHANNEL B STRUCTURAL QUANTIFICATION
            // Compute d_max directly from stop codon conversion rate at position 0.
            // Formula: d_max_B = stops_excess / convertible_original
            //
            // This directly measures the C→T damage rate at terminal positions.
            // Convertible codons (CAA, CAG, CGA) only become stops when their
            // first-position C is damaged, giving us a direct damage measurement.
            // Since d_max is a RATE (not a count), no multiplication factor needed.
            // =========================================================================
            {
                // =====================================================================
                // JOINT WLS FIT FOR b0 AND d_max (corrects baseline contamination)
                // Model: r_p = b0 + (1-b0) * d_max * x_p, where x_p = exp(-λp)
                // Fit via weighted least squares: y_p = a + c*x_p
                // Then: b0 = a, d_max = c / (1 - b0)
                // Use per-sample fitted lambda from Channel A (NOT fixed 0.3)
                // =====================================================================
                const double lambda = std::clamp(static_cast<double>(fit_lambda_5p), 0.1, 0.5);
                const int N_POSITIONS = 15;  // Use all available positions

                // Weighted least squares accumulators
                double S_w = 0, S_x = 0, S_xx = 0, S_y = 0, S_xy = 0;
                double total_exposure = 0;

                for (int p = 0; p < N_POSITIONS; ++p) {
                    double x_p = std::exp(-lambda * p);

                    double stops_p = profile.convertible_taa_5prime[p] +
                                    profile.convertible_tag_5prime[p] +
                                    profile.convertible_tga_5prime[p];
                    double pre_p = profile.convertible_caa_5prime[p] +
                                  profile.convertible_cag_5prime[p] +
                                  profile.convertible_cga_5prime[p];
                    double exposure_p = stops_p + pre_p;

                    if (exposure_p < 100) continue;  // Skip low-coverage positions

                    double y_p = stops_p / exposure_p;  // Stop fraction
                    double w_p = exposure_p;            // Weight by exposure

                    S_w  += w_p;
                    S_x  += w_p * x_p;
                    S_xx += w_p * x_p * x_p;
                    S_y  += w_p * y_p;
                    S_xy += w_p * x_p * y_p;
                    total_exposure += exposure_p;
                }

                // Solve weighted least squares: y = a + c*x
                double denom = S_w * S_xx - S_x * S_x;

                if (total_exposure > 1000 && std::abs(denom) > 1e-10) {
                    double c = (S_w * S_xy - S_x * S_y) / denom;  // Slope
                    double a = (S_y - c * S_x) / S_w;             // Intercept

                    // Store raw slope for diagnostics
                    profile.channel_b_slope = static_cast<float>(c);

                    // Slope sign is natural boundary: c > 0 = damage, c <= 0 = inverted
                    if (c > 0) {
                        double b0 = std::clamp(a, 0.01, 0.99);
                        double d_max_b = std::clamp(c / (1.0 - b0), 0.0, 1.0);

                        profile.d_max_from_channel_b = static_cast<float>(d_max_b);
                        profile.channel_b_weight = static_cast<float>(S_w);
                        profile.channel_b_quantifiable = true;
                        profile.channel_b_inverted = false;
                    } else {
                        // Inverted pattern: terminal stops LOWER than baseline
                        profile.d_max_from_channel_b = 0.0f;
                        profile.channel_b_weight = 0.0f;
                        profile.channel_b_quantifiable = false;
                        profile.channel_b_inverted = true;
                    }
                } else {
                    profile.d_max_from_channel_b = 0.0f;
                    profile.channel_b_weight = 0.0f;
                    profile.channel_b_quantifiable = false;
                    profile.channel_b_slope = 0.0f;
                }
            }

            // =========================================================================
            // LEGACY LLR DIAGNOSTICS (kept for debug output)
            // Decision now made by joint probabilistic model earlier
            // =========================================================================
            (void)profile.delta_llr_5prime;  // Used in debug output
            (void)profile.stop_decay_llr_5prime;
        }
    }

    // =========================================================================
    // CHANNEL C: Oxidative stop codon analysis (G→T transversions)
    // Unlike deamination (terminal decay), oxidation is UNIFORM across reads
    // Real oxidation: terminal rate ≈ interior rate (uniformity ratio ≈ 1)
    // =========================================================================
    {
        // Compute interior baseline for oxidative stops
        double ox_pre_interior = profile.convertible_gag_interior +
                                 profile.convertible_gaa_interior +
                                 profile.convertible_gga_interior;
        double ox_stop_interior = profile.convertible_tag_ox_interior +
                                  profile.convertible_taa_ox_interior +
                                  profile.convertible_tga_ox_interior;
        double ox_total_interior = ox_pre_interior + ox_stop_interior;

        if (ox_total_interior > 100) {
            profile.ox_stop_conversion_rate_baseline = static_cast<float>(
                ox_stop_interior / ox_total_interior);
            profile.channel_c_valid = true;
        }

        // Compute terminal (pos 0-4) vs interior (pos 5-14) stop rates
        double ox_pre_terminal = 0.0, ox_stop_terminal = 0.0;
        double ox_pre_mid = 0.0, ox_stop_mid = 0.0;

        for (int p = 0; p < 5; ++p) {
            ox_pre_terminal += profile.convertible_gag_5prime[p] +
                              profile.convertible_gaa_5prime[p] +
                              profile.convertible_gga_5prime[p];
            ox_stop_terminal += profile.convertible_tag_ox_5prime[p] +
                               profile.convertible_taa_ox_5prime[p] +
                               profile.convertible_tga_ox_5prime[p];
        }
        for (int p = 5; p < 15; ++p) {
            ox_pre_mid += profile.convertible_gag_5prime[p] +
                         profile.convertible_gaa_5prime[p] +
                         profile.convertible_gga_5prime[p];
            ox_stop_mid += profile.convertible_tag_ox_5prime[p] +
                          profile.convertible_taa_ox_5prime[p] +
                          profile.convertible_tga_ox_5prime[p];
        }

        double ox_total_terminal = ox_pre_terminal + ox_stop_terminal;
        double ox_total_mid = ox_pre_mid + ox_stop_mid;

        if (ox_total_terminal > 50 && ox_total_mid > 50) {
            profile.ox_stop_rate_terminal = static_cast<float>(ox_stop_terminal / ox_total_terminal);
            profile.ox_stop_rate_interior = static_cast<float>(ox_stop_mid / ox_total_mid);

            // Uniformity ratio: terminal/interior (≈1 for uniform oxidation)
            if (profile.ox_stop_rate_interior > 0.001f) {
                profile.ox_uniformity_ratio = profile.ox_stop_rate_terminal / profile.ox_stop_rate_interior;
            }
        }

        // Detect oxidative damage: elevated stop rate AND uniform distribution
        // Oxidation is characterized by:
        // 1. Elevated stop conversion rate (above baseline)
        // 2. Uniform distribution (uniformity ratio 0.8-1.2)
        // 3. NOT correlated with deamination damage (which also elevates terminal regions)
        float ox_stop_excess = profile.ox_stop_rate_terminal - profile.ox_stop_conversion_rate_baseline;

        // Use stricter threshold if deamination is present (correlated damage)
        // Deamination enriches ALL terminal damage, including G→T via adjacent effects
        float threshold = 0.02f;  // 2% excess required (was 0.5%)
        if (profile.d_max_combined > 5.0f || profile.damage_validated) {
            threshold = 0.05f;  // 5% excess if deamination present
        }

        bool elevated = ox_stop_excess > threshold;
        bool uniform = profile.ox_uniformity_ratio > 0.85f && profile.ox_uniformity_ratio < 1.15f;

        if (profile.channel_c_valid && elevated && uniform) {
            profile.ox_damage_detected = true;
            profile.ox_d_max = std::max(0.0f, ox_stop_excess * 100.0f);  // Convert to percentage
        }

        // Check for artifact pattern: terminal much higher than interior (like deamination)
        // This suggests contamination with C→T damage signal, not true G→T oxidation
        if (profile.ox_uniformity_ratio > 1.5f) {
            profile.ox_is_artifact = true;
            profile.ox_damage_detected = false;  // Override - not true oxidation
        }

        // Also flag as artifact if terminal is LOWER than interior (inverted pattern)
        // This suggests library prep bias or compositional artifact
        if (profile.ox_uniformity_ratio < 0.7f) {
            profile.ox_is_artifact = true;
            profile.ox_damage_detected = false;
        }
    }

    // =========================================================================
    // CHANNEL D: G→T / C→A transversion asymmetry analysis
    // Real oxidation: G→T rate > T→G rate (asymmetric)
    // Sequencing error: G→T ≈ T→G (symmetric)
    // =========================================================================
    {
        // G→T asymmetry is already captured in the control channel (A/(A+G))
        // For reference-free analysis, we use the existing a_freq/g_freq arrays
        // and check if G→T shows asymmetry relative to the complement

        // Compute G→T rate from g_count_5prime (tracked during accumulation)
        // Note: Without a reference, we can't directly measure G→T vs T→G
        // Instead, we infer from the oxidative stop codon pattern

        // For now, set asymmetry based on Channel C results
        if (profile.channel_c_valid) {
            // Use stop codon data to infer asymmetry
            // If oxidative stops are elevated, G→T is the driver
            profile.ox_gt_asymmetry = profile.ox_uniformity_ratio;
        }
    }

    // =========================================================================
    // CHANNEL E: Depurination detection (purine enrichment at termini)
    // Depurination creates strand breaks at purine sites
    // Detection: terminal purine (A+G) rate vs interior purine rate
    // =========================================================================
    {
        // Compute purine rate at terminal positions (5' end, pos 0-4)
        double purine_terminal = 0.0, total_terminal = 0.0;
        for (int p = 0; p < 5; ++p) {
            // Use the raw counts before normalization
            // a_freq_5prime and g_freq_5prime are now normalized (0-1 ratios)
            // We need to use tc_total_5prime for scaling
            double ag = profile.a_freq_5prime[p] + profile.g_freq_5prime[p];
            double tc = profile.t_freq_5prime[p] + profile.c_freq_5prime[p];
            // After normalization, a+g ≈ 1 and t+c ≈ 1, so total is meaningless
            // Use the raw totals instead
            purine_terminal += profile.tc_total_5prime[p] * (1.0 - profile.t_freq_5prime[p] - profile.c_freq_5prime[p]);
            total_terminal += profile.tc_total_5prime[p];
        }

        // Compute purine rate in middle of reads (from baseline)
        double purine_baseline = profile.baseline_a_freq + profile.baseline_g_freq;
        // baseline values are already normalized to sum to 1

        if (total_terminal > 100 && purine_baseline > 0.01) {
            // Use baseline A+G ratio as reference
            profile.purine_rate_interior = static_cast<float>(purine_baseline);

            // Terminal purine rate requires raw A+G counts
            // Since we don't have separate AG totals at 5', estimate from existing data
            // The a_freq_5prime and g_freq_5prime were captured before normalization
            // but are now ratios. We can reconstruct approximate purine enrichment
            // by comparing A/(A+G) and G totals

            // Simpler approach: use the existing terminal shift metrics
            // If both T/(T+C) AND A/(A+G) are elevated at termini, it's composition bias
            // If only T/(T+C) is elevated, it's damage
            // If A/(A+G) is elevated more than T/(T+C), check for depurination

            profile.purine_enrichment_5prime = profile.ctrl_shift_5prime;  // A/(A+G) shift at 5'
            profile.purine_enrichment_3prime = profile.ctrl_shift_3prime;  // T/(T+C) shift at 3'

            // Depurination detected if:
            // 1. Purine enrichment at 5' is significantly positive (>2%)
            // 2. AND it's not explained by composition bias (divergence from damage channel)
            if (profile.purine_enrichment_5prime > 0.02f &&
                profile.channel_divergence_5prime > 0.01f) {
                profile.depurination_detected = true;
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
    // Use both shift threshold AND z-score override for very strong statistical signals
    const bool inversion_5 = (profile.terminal_shift_5prime < -0.01f && profile.terminal_z_5prime < -2.0f)
                           || profile.terminal_z_5prime < -10.0f;  // Very strong depletion signal
    const bool inversion_3 = (profile.terminal_shift_3prime < -0.01f && profile.terminal_z_3prime < -2.0f)
                           || profile.terminal_z_3prime < -10.0f;  // Very strong depletion signal
    profile.terminal_inversion = inversion_5 || inversion_3;

    // Set individual inverted pattern flags (used for asymmetric pattern handling)
    // Note: inverted_pattern_5prime may also be set later by hexamer-based detection
    if (inversion_3) {
        profile.inverted_pattern_3prime = true;
    }
    if (inversion_5) {
        profile.inverted_pattern_5prime = true;
    }

    // =========================================================================
    // Hexamer-based damage detection (reference-independent)
    // Compare hexamer frequencies at 5' terminal vs interior positions
    // =========================================================================
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

            (void)raw_llr;  // Used for hexamer_damage_llr
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

        // If hexamer analysis shows terminal T/(T+C) < interior AND we didn't detect
        // a position-0 artifact, mark as inverted.
        // When position-0 artifact is detected, the hexamer signal may be corrupted too
        // since hexamers starting at position 0 include the artifact.
        if (profile.hexamer_damage_llr < -0.02f && !profile.position_0_artifact_5prime) {
            // Terminal T/(T+C) is LOWER than interior - no damage signal at 5' end
            profile.inverted_pattern_5prime = true;
        }
    }

    // =========================================================================
    // Composition bias detection using negative controls
    // =========================================================================
    // If the negative control (A/(A+G) at 5', T/(T+C) at 3') shows comparable
    // enrichment to the "damage" signal, it's likely composition bias, not damage.
    // Decision rule:
    //   Flag as bias if: |ctrl_shift| >= max(0.005, 0.5 * |damage_shift|)

    float damage_shift_5 = profile.terminal_shift_5prime;
    float ctrl_shift_5 = std::abs(profile.ctrl_shift_5prime);
    float threshold_5 = std::max(0.005f, 0.5f * std::abs(damage_shift_5));
    if (ctrl_shift_5 >= threshold_5 && damage_shift_5 > 0.01f) {
        profile.composition_bias_5prime = true;
    }

    float damage_shift_3 = profile.terminal_shift_3prime;
    float ctrl_shift_3 = std::abs(profile.ctrl_shift_3prime);
    float threshold_3 = std::max(0.005f, 0.5f * std::abs(damage_shift_3));
    if (ctrl_shift_3 >= threshold_3 && damage_shift_3 > 0.01f) {
        profile.composition_bias_3prime = true;
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
    } else if (profile.terminal_inversion) {
        // Inverted sample: check z-score asymmetry to distinguish real damage from composition bias
        float z5_abs = std::abs(profile.terminal_z_5prime);
        float z3_abs = std::abs(profile.terminal_z_3prime);
        float z_ratio = (z3_abs > 0) ? z5_abs / z3_abs : 10.0f;

        // Use absolute hexamer LLR magnitude for inverted samples
        float abs_llr = std::abs(profile.hexamer_damage_llr);

        if (z_ratio < 1.2f && abs_llr > 0.02f) {
            // 3' signal is relatively strong → likely real G→A damage
            // Use hexamer estimate but with more conservative scaling
            hexamer_boost = abs_llr * 8.0f;
        } else if (z_ratio > 1.5f) {
            // 5' signal dominates → likely composition bias, NOT damage
            // No boost applied
        } else {
            // Ambiguous case: use conservative estimate
            hexamer_boost = abs_llr * 4.0f;  // Half the normal scaling
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


        // Compute asymmetry: |D_5p - D_3p| / ((D_5p + D_3p) / 2)
        float d_sum = raw_d_max_5prime + raw_d_max_3prime;
        if (d_sum > 0.01f) {
            profile.asymmetry = std::abs(raw_d_max_5prime - raw_d_max_3prime) / (d_sum / 2.0f);
        } else {
            profile.asymmetry = 0.0f;
        }
        profile.high_asymmetry = (profile.asymmetry > 0.5f);

        // =========================================================================
        // GC-STRATIFIED DAMAGE CALCULATION
        // Calculate d_max for each GC bin, then aggregate
        // =========================================================================
        {
            const uint64_t MIN_C_SITES = 10000;  // Minimum C sites for valid estimate
            float weighted_sum = 0.0f;
            float weight_sum = 0.0f;
            float peak_dmax = 0.0f;
            int peak_bin = -1;

            for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
                auto& b = profile.gc_bins[bin];

                // Sum up C sites for this bin
                b.c_sites = b.c_interior;
                for (int p = 0; p < 15; ++p) {
                    b.c_sites += b.c_counts[p];
                }

                if (b.n_reads < 1000 || b.c_sites < MIN_C_SITES) {
                    continue;  // Skip bins with insufficient data
                }

                // Channel A: calculate d_max = (t_terminal - t_baseline) / c_baseline
                double t_baseline = static_cast<double>(b.t_interior) /
                                   (b.t_interior + b.c_interior + 1);
                double c_baseline = 1.0 - t_baseline;
                double t_terminal = static_cast<double>(b.t_counts[0]) /
                                   (b.t_counts[0] + b.c_counts[0] + 1);

                if (c_baseline > 0.1) {
                    b.d_max = std::clamp(static_cast<float>((t_terminal - t_baseline) / c_baseline),
                                         0.0f, 1.0f);
                }

                // Channel B: calculate d_max from stop codon conversion
                double stop_baseline = static_cast<double>(b.stop_interior) /
                                      (b.stop_interior + b.pre_interior + 1);
                double stop_terminal = static_cast<double>(b.stop_counts[0]) /
                                      (b.stop_counts[0] + b.pre_counts[0] + 1);

                if (stop_baseline < 0.99) {
                    b.d_max_channel_b = std::clamp(
                        static_cast<float>((stop_terminal - stop_baseline) / (1.0 - stop_baseline)),
                        0.0f, 1.0f);
                }

                b.valid = true;

                // Use max of Channel A and B (they measure same thing)
                float bin_dmax = std::max(b.d_max, b.d_max_channel_b);

                // Weight by C sites (more C = more signal)
                float weight = static_cast<float>(b.c_sites);
                weighted_sum += bin_dmax * weight;
                weight_sum += weight;

                // Track peak
                if (bin_dmax > peak_dmax) {
                    peak_dmax = bin_dmax;
                    peak_bin = bin;
                }

                int gc_low = bin * 10;
                int gc_high = gc_low + 10;
            }

            if (weight_sum > 0) {
                profile.gc_stratified_d_max_weighted = weighted_sum / weight_sum;
                profile.gc_stratified_d_max_peak = peak_dmax;
                profile.gc_peak_bin = peak_bin;
                profile.gc_stratified_valid = true;


                // =========================================================================
                // GC-CONDITIONAL DAMAGE CLASSIFICATION
                // Per-bin LLR classification and aggregate metrics (π_damaged, d_ancient)
                // This is the principled approach recommended by OpenCode discussion
                // =========================================================================
                {
                    constexpr float LLR_THRESHOLD = 10.0f;  // Log-likelihood ratio threshold for classification
                    constexpr float MIN_DMAX_THRESHOLD = 0.01f;  // Minimum d_max to consider as damaged

                    uint64_t total_obs = 0;      // Total terminal observations
                    uint64_t damaged_obs = 0;    // Observations from damaged bins
                    float damaged_weighted_d = 0.0f;
                    uint64_t damaged_weight = 0;
                    int n_damaged = 0;

                    for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
                        auto& b = profile.gc_bins[bin];
                        if (!b.valid) continue;

                        // Compute per-bin baseline
                        double baseline_tc = static_cast<double>(b.t_interior) /
                                            std::max(1.0, static_cast<double>(b.t_interior + b.c_interior));
                        b.baseline_tc = static_cast<float>(std::clamp(baseline_tc, 0.01, 0.99));

                        // Compute terminal observation count
                        uint64_t n_obs = b.n_terminal_obs();
                        total_obs += n_obs;

                        // Compute per-bin LLR: damaged vs undamaged model
                        // LL(damaged) - LL(undamaged) at terminal positions
                        float ll_damaged = 0.0f;
                        float ll_undamaged = 0.0f;
                        float lambda = profile.lambda_5prime;

                        for (int p = 0; p < 15; ++p) {
                            float decay = std::exp(-lambda * p);
                            float delta_p = b.d_max * decay;

                            // Undamaged model: P(T) = baseline_tc
                            float pi_undamaged = b.baseline_tc;
                            // Damaged model: P(T) = baseline + (1-baseline) * delta_p
                            float pi_damaged = b.baseline_tc + (1.0f - b.baseline_tc) * delta_p;

                            pi_undamaged = std::clamp(pi_undamaged, 0.001f, 0.999f);
                            pi_damaged = std::clamp(pi_damaged, 0.001f, 0.999f);

                            double k = static_cast<double>(b.t_counts[p]);
                            double n = static_cast<double>(b.t_counts[p] + b.c_counts[p]);
                            if (n > 0) {
                                ll_damaged += static_cast<float>(k * std::log(pi_damaged) + (n - k) * std::log(1.0f - pi_damaged));
                                ll_undamaged += static_cast<float>(k * std::log(pi_undamaged) + (n - k) * std::log(1.0f - pi_undamaged));
                            }
                        }

                        b.llr = ll_damaged - ll_undamaged;

                        // Classify bin as damaged if LLR > threshold AND d_max > minimum
                        b.classified_damaged = (b.llr > LLR_THRESHOLD) && (b.d_max > MIN_DMAX_THRESHOLD);

                        // Convert LLR to soft probability (logistic)
                        // p_damaged = sigmoid(LLR - threshold) scaled to [0, 1]
                        float llr_centered = b.llr - LLR_THRESHOLD;
                        b.p_damaged = 1.0f / (1.0f + std::exp(-0.5f * llr_centered));

                        // Require minimum d_max for non-trivial p_damaged
                        if (b.d_max < MIN_DMAX_THRESHOLD) {
                            b.p_damaged = 0.0f;
                        }

                        // Accumulate for aggregates
                        if (b.classified_damaged) {
                            damaged_obs += n_obs;
                            damaged_weighted_d += b.d_max * static_cast<float>(n_obs);
                            damaged_weight += n_obs;
                            ++n_damaged;
                        }
                    }

                    // Compute aggregate metrics
                    if (total_obs > 0) {
                        profile.pi_damaged = static_cast<float>(damaged_obs) / static_cast<float>(total_obs);
                    }
                    if (damaged_weight > 0) {
                        profile.d_ancient = damaged_weighted_d / static_cast<float>(damaged_weight);
                    }

                    // d_population = E[δ] across all bins (weighted by observations)
                    float pop_weighted_d = 0.0f;
                    for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
                        const auto& b = profile.gc_bins[bin];
                        if (!b.valid) continue;
                        pop_weighted_d += b.d_max * static_cast<float>(b.n_terminal_obs());
                    }
                    if (total_obs > 0) {
                        profile.d_population = pop_weighted_d / static_cast<float>(total_obs);
                    }

                    profile.n_damaged_bins = n_damaged;

                }

                // =========================================================================
                // DURBIN-STYLE K-COMPONENT MIXTURE MODEL
                // Models organism-class heterogeneity for principled damage quantification
                // d_population = what we observe, d_reference = metaDMG proxy
                // =========================================================================
                double total_c_sites = static_cast<double>(weight_sum);
                std::array<SuperRead, N_GC_BINS> super_reads;
                for (int bin = 0; bin < N_GC_BINS; ++bin) {
                    const auto& b = profile.gc_bins[bin];
                    super_reads[bin].gc_bin = bin;
                    super_reads[bin].c_sites = static_cast<double>(b.c_sites);

                    // Channel A: T/(T+C) - use the GC bin's terminal counts
                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_tc[p] = static_cast<double>(b.t_counts[p]);
                        super_reads[bin].n_tc[p] = static_cast<double>(b.t_counts[p] + b.c_counts[p]);
                    }

                    // Control: A/(A+G) at 5' - approximate from global profile scaled by bin weight
                    // (We don't track A/G per GC bin, so use global proportionally)
                    double bin_fraction = b.c_sites > 0 ? static_cast<double>(b.c_sites) / total_c_sites : 0.0;
                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_ag[p] = profile.a_freq_5prime[p] * bin_fraction * base_ag_total;
                        super_reads[bin].n_ag[p] = (profile.a_freq_5prime[p] + profile.g_freq_5prime[p]) * bin_fraction * base_ag_total;
                    }

                    // Channel B: stop conversion
                    for (int p = 0; p < N_POSITIONS; ++p) {
                        super_reads[bin].k_stop[p] = static_cast<double>(b.stop_counts[p]);
                        super_reads[bin].n_stop[p] = static_cast<double>(b.stop_counts[p] + b.pre_counts[p]);
                    }

                    // Interior baselines
                    super_reads[bin].k_tc_int = static_cast<double>(b.t_interior);
                    super_reads[bin].n_tc_int = static_cast<double>(b.t_interior + b.c_interior);
                    super_reads[bin].k_stop_int = static_cast<double>(b.stop_interior);
                    super_reads[bin].n_stop_int = static_cast<double>(b.stop_interior + b.pre_interior);

                    // Control interior - approximate
                    super_reads[bin].k_ag_int = profile.baseline_a_freq * bin_fraction * base_ag_total;
                    super_reads[bin].n_ag_int = (profile.baseline_a_freq + profile.baseline_g_freq) * bin_fraction * base_ag_total;
                }

                auto mixture_result = MixtureDamageModel::fit(super_reads);
                profile.mixture_K = mixture_result.K;
                profile.mixture_d_population = mixture_result.d_population;
                profile.mixture_d_ancient = mixture_result.d_ancient;
                profile.mixture_d_reference = mixture_result.d_reference;
                profile.mixture_pi_ancient = mixture_result.pi_ancient;
                profile.mixture_bic = mixture_result.bic;
                profile.mixture_converged = mixture_result.converged;

                if (mixture_result.converged) {
                    // Print per-class details
                    for (int k = 0; k < mixture_result.K; ++k) {
                    }
                } else {
                }
            }
        }

        // =========================================================================
        // JOINT MODEL FOR CLASSIFICATION + MIXTURE MODEL FOR QUANTIFICATION
        // Joint model P(damage) for hypothesis testing
        // Mixture model d_reference for metaDMG-comparable magnitude
        // =========================================================================
        profile.d_max_5prime = raw_d_max_5prime;  // Keep raw estimates for reference
        profile.d_max_3prime = raw_d_max_3prime;

        if (profile.damage_artifact) {
            // Artifact detected (P(damage) < 0.05 with positive artifact term)
            profile.d_max_5prime = 0.0f;
            profile.d_max_3prime = 0.0f;
            profile.d_max_combined = 0.0f;
            profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
        } else if (profile.damage_validated) {
            // Damage validated (P(damage) > 0.95)
            // Check for inversion patterns or position-0 artifacts: Channel A may be unreliable
            bool channel_a_unreliable = (profile.inverted_pattern_5prime && profile.inverted_pattern_3prime)
                                      || profile.position_0_artifact_5prime
                                      || profile.position_0_artifact_3prime;

            if (channel_a_unreliable && profile.channel_b_quantifiable && profile.d_max_from_channel_b > 0.01f) {
                // Channel A corrupted by inversion or position-0 artifact - use Channel B
                profile.d_max_combined = profile.d_max_from_channel_b;
                profile.d_max_source = SampleDamageProfile::DmaxSource::CHANNEL_B_STRUCTURAL;
            } else if (profile.inverted_pattern_3prime && !profile.inverted_pattern_5prime) {
                // 3' inverted but 5' valid - use 5' raw estimate
                profile.d_max_combined = raw_d_max_5prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::FIVE_PRIME_ONLY;
            } else if (profile.inverted_pattern_5prime && !profile.inverted_pattern_3prime) {
                // 5' inverted but 3' valid - use 3' raw estimate
                profile.d_max_combined = raw_d_max_3prime;
                profile.d_max_source = SampleDamageProfile::DmaxSource::THREE_PRIME_ONLY;
            } else if (profile.mixture_converged && profile.mixture_d_reference > 0.01f) {
                // No inversion - use mixture model d_reference (metaDMG proxy: E[δ | GC > 50%])
                profile.d_max_combined = profile.mixture_d_reference;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else if (profile.gc_stratified_valid) {
                // Fallback to GC-weighted average
                profile.d_max_combined = profile.gc_stratified_d_max_weighted;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else {
                // Fallback to joint model
                profile.d_max_combined = profile.joint_delta_max;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
        } else if (profile.joint_model_valid && profile.joint_p_damage > 0.5f) {
            // Uncertain but leaning toward damage
            if (profile.mixture_converged && profile.mixture_d_reference > 0.01f) {
                profile.d_max_combined = profile.mixture_d_reference;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else if (profile.gc_stratified_valid) {
                profile.d_max_combined = profile.gc_stratified_d_max_weighted;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            } else {
                profile.d_max_combined = profile.joint_delta_max;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
        } else if (profile.joint_model_valid) {
            // Joint model valid and does NOT support damage.
            // EXCEPTION: For single-stranded libraries, Channel B tracks 5' C→T stops which
            // are absent in ss libraries (ss damage is G→A at 3'). Fall back to Channel A.
            const bool is_ss = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                               (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
            if (is_ss) {
                // Use Channel A for ss libraries since Channel B is not applicable
                profile.d_max_5prime = raw_d_max_5prime;
                profile.d_max_3prime = raw_d_max_3prime;
                profile.d_max_combined = std::max(raw_d_max_5prime, raw_d_max_3prime);
                profile.d_max_source = SampleDamageProfile::DmaxSource::MAX_SS_ASYMMETRY;
            } else {
                // For ds libraries: check if BOTH channels show positive signal
                // This is the key insight: if terminal_shift > 0 AND stop_decay_llr > 0,
                // we have concordant evidence from two independent channels = real damage.
                // Only zero out if channels disagree or both are negative.
                bool channel_a_positive = profile.terminal_shift_5prime > 0.005f || profile.terminal_shift_3prime > 0.005f;
                bool channel_b_positive = profile.stop_decay_llr_5prime > 0.0f;

                if (channel_a_positive && channel_b_positive) {
                    // Concordant positive signal from both channels = real damage
                    // Use Channel B d_max if quantifiable, otherwise use raw Channel A
                    if (profile.channel_b_quantifiable && profile.d_max_from_channel_b > 0.001f) {
                        profile.d_max_5prime = raw_d_max_5prime;
                        profile.d_max_3prime = raw_d_max_3prime;
                        profile.d_max_combined = profile.d_max_from_channel_b;
                        profile.d_max_source = SampleDamageProfile::DmaxSource::CHANNEL_B_STRUCTURAL;
                    } else {
                        // Use raw Channel A values (symmetric for ds)
                        profile.d_max_5prime = raw_d_max_5prime;
                        profile.d_max_3prime = raw_d_max_3prime;
                        profile.d_max_combined = (raw_d_max_5prime + raw_d_max_3prime) / 2.0f;
                        profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
                    }
                } else {
                    // Channels disagree or both negative - compositional artifact
                    profile.d_max_5prime = 0.0f;
                    profile.d_max_3prime = 0.0f;
                    profile.d_max_combined = 0.0f;
                    profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
                }
            }
        } else if (profile.channel_b_quantifiable) {
            // Channel B WLS fit succeeded - use its d_max directly
            // This path is for samples with some stop signal but not validated
            profile.d_max_5prime = raw_d_max_5prime;
            profile.d_max_3prime = raw_d_max_3prime;
            profile.d_max_combined = profile.d_max_from_channel_b;
            profile.d_max_source = SampleDamageProfile::DmaxSource::CHANNEL_B_STRUCTURAL;
        } else {
            // No Channel B quantification - use Channel A raw values
            profile.d_max_5prime = raw_d_max_5prime;
            profile.d_max_3prime = raw_d_max_3prime;

            if (profile.high_asymmetry) {
                // For single-stranded libraries, asymmetry is EXPECTED (damage at one terminus only)
                // Use max() to capture the damaged terminus, not min() which would give 0%
                const bool is_ss = (profile.forced_library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) ||
                                   (profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED);
                if (is_ss) {
                    profile.d_max_combined = std::max(profile.d_max_5prime, profile.d_max_3prime);
                    profile.d_max_source = SampleDamageProfile::DmaxSource::MAX_SS_ASYMMETRY;
                } else {
                    // For ds libraries, high asymmetry suggests artifact - use conservative min
                    profile.d_max_combined = std::min(profile.d_max_5prime, profile.d_max_3prime);
                    profile.d_max_source = SampleDamageProfile::DmaxSource::MIN_ASYMMETRY;
                }
            } else {
                profile.d_max_combined = (profile.d_max_5prime + profile.d_max_3prime) / 2.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::AVERAGE;
            }
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
                profile.d_max_5prime = 0.0f;
                profile.d_max_3prime = 0.0f;
                profile.d_max_combined = 0.0f;
                profile.d_max_source = SampleDamageProfile::DmaxSource::NONE;
            }
        }

    }

    // =========================================================================
    // D_METAMATCH CALCULATION: Channel B-anchored damage estimate
    //
    // metaDMG uses aligned reads (selection bias toward well-preserved sequences).
    // DART's d_global uses ALL reads. The gap arises because alignable reads
    // tend to show cleaner damage signal.
    //
    // APPROACH: Use Channel B (stop codon conversion) as the primary anchor.
    // Channel B directly measures C→T damage via CAA→TAA, CAG→TAG, CGA→TGA
    // conversions, which is biologically equivalent to what metaDMG measures.
    //
    // Formula:
    //   d_metamatch = d_global + γ × (d_channel_b - d_global)
    //
    // Where γ is based on Channel B confidence (higher LLR = more trust in
    // Channel B estimate). For highly damaged, validated samples, this pulls
    // d_global toward the Channel B structural estimate.
    // =========================================================================
    {
        // Step 1: Compute GC-weighted d_max for diagnostics (kept for reporting)
        double weighted_sum = 0.0;
        double weight_sum = 0.0;

        for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
            const auto& b = profile.gc_bins[bin];
            if (!b.valid || b.n_reads < 100) continue;

            float gc_center = (bin + 0.5f) / static_cast<float>(SampleDamageProfile::N_GC_BINS);
            float gc_deviation = gc_center - 0.50f;
            float alignability_weight = std::exp(-gc_deviation * gc_deviation / (2.0f * 0.15f * 0.15f));
            double read_weight = std::sqrt(static_cast<double>(b.n_reads));
            double total_weight = alignability_weight * read_weight;

            weighted_sum += b.d_max * total_weight;
            weight_sum += total_weight;
        }

        if (weight_sum > 0) {
            profile.d_alignability_weighted = static_cast<float>(weighted_sum / weight_sum);
        } else {
            profile.d_alignability_weighted = profile.d_max_combined;
        }

        // Step 2: Compute confidence coefficient (γ) for Channel B blending
        // High LLR = confident damage → trust Channel B estimate
        // Low/negative LLR = uncertain → stay close to d_global
        float channel_b_llr = profile.stop_decay_llr_5prime;

        // Sigmoid scaling: γ ranges from 0 to 1
        // LLR = 0 → γ = 0.5, LLR = 50000 → γ ≈ 0.92, LLR = 100000 → γ ≈ 0.99
        float gamma_raw = 1.0f / (1.0f + std::exp(-channel_b_llr / 50000.0f));

        // Step 3: Determine blending strategy based on validation state
        float d_global = profile.d_max_combined;
        float d_channel_b = profile.d_max_from_channel_b;

        if (!profile.damage_validated || profile.damage_artifact) {
            // No damage or artifact: d_metamatch = d_global (typically 0)
            profile.metamatch_gamma = 0.0f;
            profile.d_metamatch = d_global;
        } else if (profile.channel_b_quantifiable && d_channel_b > 0.01f) {
            // Channel B is quantifiable: blend toward it
            // Apply asymmetric blending:
            // - If Channel B > d_global: use γ to pull UP (metaDMG usually higher)
            // - If Channel B < d_global: use weaker pull to avoid under-estimation
            if (d_channel_b > d_global) {
                profile.metamatch_gamma = gamma_raw;
            } else {
                // More conservative when Channel B suggests LOWER damage
                profile.metamatch_gamma = 0.3f * gamma_raw;
            }
            profile.d_metamatch = d_global + profile.metamatch_gamma * (d_channel_b - d_global);
        } else {
            // Channel B not quantifiable: use GC-weighted estimate with softer blending
            profile.metamatch_gamma = 0.5f * gamma_raw;
            profile.d_metamatch = d_global + profile.metamatch_gamma * (profile.d_alignability_weighted - d_global);
        }

        // Ensure d_metamatch is non-negative and bounded
        profile.d_metamatch = std::clamp(profile.d_metamatch, 0.0f, 1.0f);

        // Compute mean alignability for diagnostics
        double alignability_total = 0.0;
        double n_total = 0.0;
        for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
            const auto& b = profile.gc_bins[bin];
            if (b.n_reads == 0) continue;

            float gc_center = (bin + 0.5f) / static_cast<float>(SampleDamageProfile::N_GC_BINS);
            float gc_deviation = gc_center - 0.50f;
            float alignability = std::exp(-gc_deviation * gc_deviation / (2.0f * 0.15f * 0.15f));

            alignability_total += alignability * b.n_reads;
            n_total += b.n_reads;
        }
        profile.mean_alignability = (n_total > 0) ? static_cast<float>(alignability_total / n_total) : 0.5f;
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

    // Merge Channel C oxidative codon counts (G→T transversions)
    for (int i = 0; i < 15; ++i) {
        dst.convertible_gag_5prime[i] += src.convertible_gag_5prime[i];
        dst.convertible_tag_ox_5prime[i] += src.convertible_tag_ox_5prime[i];
        dst.convertible_gaa_5prime[i] += src.convertible_gaa_5prime[i];
        dst.convertible_taa_ox_5prime[i] += src.convertible_taa_ox_5prime[i];
        dst.convertible_gga_5prime[i] += src.convertible_gga_5prime[i];
        dst.convertible_tga_ox_5prime[i] += src.convertible_tga_ox_5prime[i];
        // Channel D: G count for asymmetry tracking
        dst.g_count_5prime[i] += src.g_count_5prime[i];
    }
    dst.convertible_gag_interior += src.convertible_gag_interior;
    dst.convertible_tag_ox_interior += src.convertible_tag_ox_interior;
    dst.convertible_gaa_interior += src.convertible_gaa_interior;
    dst.convertible_taa_ox_interior += src.convertible_taa_ox_interior;
    dst.convertible_gga_interior += src.convertible_gga_interior;
    dst.convertible_tga_ox_interior += src.convertible_tga_ox_interior;

    // Merge GC-stratified bin counts
    for (int bin = 0; bin < SampleDamageProfile::N_GC_BINS; ++bin) {
        auto& db = dst.gc_bins[bin];
        const auto& sb = src.gc_bins[bin];
        for (int p = 0; p < 15; ++p) {
            db.t_counts[p] += sb.t_counts[p];
            db.c_counts[p] += sb.c_counts[p];
            db.stop_counts[p] += sb.stop_counts[p];
            db.pre_counts[p] += sb.pre_counts[p];
        }
        db.t_interior += sb.t_interior;
        db.c_interior += sb.c_interior;
        db.stop_interior += sb.stop_interior;
        db.pre_interior += sb.pre_interior;
        db.n_reads += sb.n_reads;
    }

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

    // Joint probabilistic model
    profile.joint_delta_max = 0.0f;
    profile.joint_lambda = 0.0f;
    profile.joint_a_max = 0.0f;
    profile.joint_log_lik_m1 = 0.0f;
    profile.joint_log_lik_m0 = 0.0f;
    profile.joint_delta_bic = 0.0f;
    profile.joint_bayes_factor = 0.0f;
    profile.joint_p_damage = 0.0f;
    profile.joint_model_valid = false;

    // Durbin-style mixture model
    profile.mixture_K = 0;
    profile.mixture_d_population = 0.0f;
    profile.mixture_d_ancient = 0.0f;
    profile.mixture_d_reference = 0.0f;
    profile.mixture_pi_ancient = 0.0f;
    profile.mixture_bic = 0.0f;
    profile.mixture_converged = false;

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

} // namespace dart
