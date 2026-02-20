#pragma once

/**
 * @file adaptive_damage.hpp
 * @brief Adaptive damage detection and correction calibration
 *
 * Implements advanced damage detection that separates true ancient DNA damage
 * from composition bias, with online calibration and validation.
 */

#include "frame_selector.hpp"
#include "hexamer_tables.hpp"
#include <array>
#include <vector>
#include <atomic>
#include <cmath>
#include <algorithm>
#include <optional>

namespace dart {

/**
 * Terminal k-mer statistics for damage vs composition separation
 */
struct TerminalKmerStats {
    // 5' terminal (positions 0-5): T/(T+C) ratios
    std::array<double, 6> terminal_5prime_t = {};
    std::array<double, 6> terminal_5prime_c = {};

    // 3' terminal (positions 0-5 from end): A/(A+G) ratios
    std::array<double, 6> terminal_3prime_a = {};
    std::array<double, 6> terminal_3prime_g = {};

    // Interior reference (positions 15-25): baseline composition
    double interior_t = 0.0;
    double interior_c = 0.0;
    double interior_a = 0.0;
    double interior_g = 0.0;

    size_t n_reads = 0;

    // Compute decay slope (real damage shows decay, AT-bias is flat)
    float compute_decay_slope_5prime() const {
        if (n_reads < 100) return 0.0f;

        // Linear regression: T/(T+C) vs position
        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
        int n = 0;
        for (int i = 0; i < 6; ++i) {
            double total = terminal_5prime_t[i] + terminal_5prime_c[i];
            if (total < 10) continue;
            double y = terminal_5prime_t[i] / total;
            sum_x += i;
            sum_y += y;
            sum_xy += i * y;
            sum_x2 += i * i;
            n++;
        }
        if (n < 4) return 0.0f;

        double denom = n * sum_x2 - sum_x * sum_x;
        if (std::abs(denom) < 1e-10) return 0.0f;

        return static_cast<float>((n * sum_xy - sum_x * sum_y) / denom);
    }

    float compute_decay_slope_3prime() const {
        if (n_reads < 100) return 0.0f;

        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
        int n = 0;
        for (int i = 0; i < 6; ++i) {
            double total = terminal_3prime_a[i] + terminal_3prime_g[i];
            if (total < 10) continue;
            double y = terminal_3prime_a[i] / total;
            sum_x += i;
            sum_y += y;
            sum_xy += i * y;
            sum_x2 += i * i;
            n++;
        }
        if (n < 4) return 0.0f;

        double denom = n * sum_x2 - sum_x * sum_x;
        if (std::abs(denom) < 1e-10) return 0.0f;

        return static_cast<float>((n * sum_xy - sum_x * sum_y) / denom);
    }

    // Check if terminal enrichment is flat (AT-bias, not damage)
    bool is_flat_enrichment() const {
        float slope_5 = compute_decay_slope_5prime();
        float slope_3 = compute_decay_slope_3prime();

        // Real damage: slope should be negative (decay from terminus)
        // AT-bias: slope near zero (flat)
        return slope_5 > -0.005f && slope_3 > -0.005f;
    }

    // Compute mixture probability: P(damage | terminal pattern)
    float damage_posterior() const {
        if (n_reads < 100) return 0.5f;

        // Compute terminal excess over interior
        double interior_tc = interior_t / (interior_t + interior_c + 1);
        double terminal_tc = terminal_5prime_t[0] / (terminal_5prime_t[0] + terminal_5prime_c[0] + 1);

        double interior_ag = interior_a / (interior_a + interior_g + 1);
        double terminal_ag = terminal_3prime_a[0] / (terminal_3prime_a[0] + terminal_3prime_g[0] + 1);

        double excess_5 = terminal_tc - interior_tc;
        double excess_3 = terminal_ag - interior_ag;

        // Decay slope (negative = real damage)
        float slope_5 = compute_decay_slope_5prime();
        float slope_3 = compute_decay_slope_3prime();

        // Combine signals into posterior
        // Strong negative slope + positive excess → high damage probability
        float signal = 0.0f;

        if (excess_5 > 0.05 && slope_5 < -0.01) signal += 0.3f;
        if (excess_3 > 0.05 && slope_3 < -0.01) signal += 0.3f;
        if (excess_5 > 0.10) signal += 0.15f;
        if (excess_3 > 0.10) signal += 0.15f;
        if (slope_5 < -0.02) signal += 0.05f;
        if (slope_3 < -0.02) signal += 0.05f;

        return std::clamp(signal, 0.0f, 1.0f);
    }
};

/**
 * Per-batch statistics for online calibration
 */
struct BatchStats {
    size_t n_reads = 0;
    size_t n_corrected = 0;
    size_t total_corrections = 0;
    size_t stop_restorations = 0;

    // Coding quality before/after correction
    double total_coding_score_before = 0.0;
    double total_coding_score_after = 0.0;
    double total_hexamer_llr_before = 0.0;
    double total_hexamer_llr_after = 0.0;

    // Damage statistics
    double sum_damage_llr = 0.0;  // Sum of per-read damage LLRs
    float max_damage_rate = 0.0f;

    // Per-domain statistics for ensemble calibration
    std::array<double, 8> domain_hexamer_sums = {};  // Sum of hexamer LLRs per domain
    std::array<size_t, 8> domain_counts = {};        // Reads per domain

    void reset() {
        n_reads = n_corrected = total_corrections = stop_restorations = 0;
        total_coding_score_before = total_coding_score_after = 0.0;
        total_hexamer_llr_before = total_hexamer_llr_after = 0.0;
        sum_damage_llr = 0.0;
        max_damage_rate = 0.0f;
        domain_hexamer_sums.fill(0.0);
        domain_counts.fill(0);
    }

    // Check if corrections are improving quality
    bool corrections_improving() const {
        if (n_corrected < 10) return true;  // Not enough data
        return total_coding_score_after >= total_coding_score_before * 0.95 &&
               total_hexamer_llr_after >= total_hexamer_llr_before * 0.90;
    }

    // Correction rate
    float correction_rate() const {
        return n_reads > 0 ? static_cast<float>(n_corrected) / n_reads : 0.0f;
    }
};

/**
 * Sliding window for on-the-fly decay monitoring
 */
struct DecayMonitor {
    static constexpr size_t WINDOW_SIZE = 10000;

    // Circular buffer of per-read decay estimates
    std::vector<float> decay_estimates_5prime;
    std::vector<float> decay_estimates_3prime;
    size_t write_pos = 0;
    size_t n_samples = 0;

    DecayMonitor() : decay_estimates_5prime(WINDOW_SIZE, 0.0f),
                     decay_estimates_3prime(WINDOW_SIZE, 0.0f) {}

    void add_sample(float decay_5, float decay_3) {
        decay_estimates_5prime[write_pos] = decay_5;
        decay_estimates_3prime[write_pos] = decay_3;
        write_pos = (write_pos + 1) % WINDOW_SIZE;
        if (n_samples < WINDOW_SIZE) n_samples++;
    }

    // Compute mean decay slope over window
    std::pair<float, float> mean_decay() const {
        if (n_samples < 100) return {0.0f, 0.0f};

        double sum_5 = 0, sum_3 = 0;
        size_t count = std::min(n_samples, WINDOW_SIZE);
        for (size_t i = 0; i < count; ++i) {
            sum_5 += decay_estimates_5prime[i];
            sum_3 += decay_estimates_3prime[i];
        }
        return {static_cast<float>(sum_5 / count),
                static_cast<float>(sum_3 / count)};
    }

    // Check if decay signal is collapsing
    bool is_decay_collapsing() const {
        auto [mean_5, mean_3] = mean_decay();
        // Real damage should show negative slopes
        // Collapsing = slopes approaching zero or positive
        return mean_5 > -0.005f && mean_3 > -0.005f;
    }
};

/**
 * Adaptive thresholds that adjust based on batch statistics
 */
struct AdaptiveThresholds {
    // Higher base threshold (0.50) to be more conservative
    // Most T's at terminal positions are NOT damage even with high sample damage
    float base_correction_threshold = 0.50f;
    float current_correction_threshold = 0.50f;

    float base_damage_llr_threshold = 0.5f;
    float current_damage_llr_threshold = 0.5f;

    int max_corrections_per_end = 3;
    int current_max_corrections = 3;

    bool corrections_enabled = true;
    bool stop_restoration_enabled = true;

    // Tighten thresholds (more conservative)
    void tighten() {
        current_correction_threshold = std::min(0.6f, current_correction_threshold + 0.05f);
        current_damage_llr_threshold = std::min(1.0f, current_damage_llr_threshold + 0.1f);
        current_max_corrections = std::max(1, current_max_corrections - 1);
    }

    // Loosen thresholds (more aggressive)
    void loosen() {
        current_correction_threshold = std::max(0.2f, current_correction_threshold - 0.02f);
        current_damage_llr_threshold = std::max(0.3f, current_damage_llr_threshold - 0.05f);
        current_max_corrections = std::min(5, current_max_corrections + 1);
    }

    // Reset to base values
    void reset() {
        current_correction_threshold = base_correction_threshold;
        current_damage_llr_threshold = base_damage_llr_threshold;
        current_max_corrections = max_corrections_per_end;
    }

    // Disable corrections entirely
    void disable_corrections() {
        corrections_enabled = false;
        stop_restoration_enabled = false;
    }
};

/**
 * Ensemble domain weights for META mode
 */
struct EnsembleWeights {
    std::array<float, 8> weights = {0.125f, 0.125f, 0.125f, 0.125f,
                                    0.125f, 0.125f, 0.125f, 0.125f};
    bool weights_calibrated = false;
    float uncertainty = 1.0f;  // High uncertainty = flat weights

    // Update weights based on observed hexamer LLRs
    void update_from_batch(const BatchStats& stats) {
        double total = 0.0;
        for (size_t i = 0; i < 8; ++i) {
            if (stats.domain_counts[i] > 0) {
                double mean_llr = stats.domain_hexamer_sums[i] / stats.domain_counts[i];
                // Convert LLR to weight (higher LLR = better fit)
                weights[i] = static_cast<float>(std::exp(mean_llr * 0.1));
                total += weights[i];
            }
        }

        if (total > 0) {
            for (auto& w : weights) w /= static_cast<float>(total);
            weights_calibrated = true;

            // Compute uncertainty: entropy of weight distribution
            float entropy = 0.0f;
            for (auto w : weights) {
                if (w > 0.001f) entropy -= w * std::log(w);
            }
            // Max entropy for 8 domains = ln(8) ≈ 2.08
            uncertainty = entropy / 2.08f;
        }
    }

    // Check if weights are too flat (uncertain)
    bool is_uncertain() const {
        return !weights_calibrated || uncertainty > 0.8f;
    }

    // Get ensemble-weighted hexamer LLR
    float get_weighted_llr(const std::array<float, 8>& domain_llrs) const {
        float weighted = 0.0f;
        for (size_t i = 0; i < 8; ++i) {
            weighted += weights[i] * domain_llrs[i];
        }
        return weighted;
    }
};

/**
 * Per-read damage likelihood ratio
 */
struct ReadDamageLLR {
    float llr_5prime = 0.0f;  // Log-likelihood ratio for 5' damage
    float llr_3prime = 0.0f;  // Log-likelihood ratio for 3' damage
    float combined_llr = 0.0f;

    bool is_low_signal() const {
        return combined_llr < 0.1f;  // Below threshold
    }

    static ReadDamageLLR compute(const std::string& seq,
                                  const SampleDamageProfile& profile) {
        ReadDamageLLR result;
        if (seq.length() < 20) return result;

        size_t len = seq.length();

        // 5' end: compare observed T/(T+C) to expected under damage vs no-damage
        double t_count = 0, c_count = 0;
        for (size_t i = 0; i < std::min<size_t>(5, len); ++i) {
            char b = std::toupper(static_cast<unsigned char>(seq[i]));
            if (b == 'T') t_count++;
            else if (b == 'C') c_count++;
        }

        if (t_count + c_count > 0) {
            double obs_ratio = t_count / (t_count + c_count);
            double base_ratio = profile.baseline_t_freq /
                               (profile.baseline_t_freq + profile.baseline_c_freq + 0.001);
            double damage_ratio = base_ratio + profile.max_damage_5prime * (1.0 - base_ratio);

            // LLR = log(P(obs|damage) / P(obs|no-damage))
            double p_damage = std::pow(damage_ratio, t_count) * std::pow(1-damage_ratio, c_count);
            double p_no_damage = std::pow(base_ratio, t_count) * std::pow(1-base_ratio, c_count);

            if (p_no_damage > 1e-100) {
                result.llr_5prime = static_cast<float>(std::log(p_damage / p_no_damage + 1e-100));
            }
        }

        // 3' end: similar for A/(A+G)
        double a_count = 0, g_count = 0;
        for (size_t i = len - std::min<size_t>(5, len); i < len; ++i) {
            char b = std::toupper(static_cast<unsigned char>(seq[i]));
            if (b == 'A') a_count++;
            else if (b == 'G') g_count++;
        }

        if (a_count + g_count > 0) {
            double obs_ratio = a_count / (a_count + g_count);
            double base_ratio = profile.baseline_a_freq /
                               (profile.baseline_a_freq + profile.baseline_g_freq + 0.001);
            double damage_ratio = base_ratio + profile.max_damage_3prime * (1.0 - base_ratio);

            double p_damage = std::pow(damage_ratio, a_count) * std::pow(1-damage_ratio, g_count);
            double p_no_damage = std::pow(base_ratio, a_count) * std::pow(1-base_ratio, g_count);

            if (p_no_damage > 1e-100) {
                result.llr_3prime = static_cast<float>(std::log(p_damage / p_no_damage + 1e-100));
            }
        }

        result.combined_llr = (result.llr_5prime + result.llr_3prime) / 2.0f;
        return result;
    }
};

/**
 * Main adaptive damage calibration class
 */
class AdaptiveDamageCalibrator {
public:
    AdaptiveDamageCalibrator() = default;

    // Initialize from sample profile
    void initialize(const SampleDamageProfile& profile) {
        sample_profile_ = profile;
        thresholds_.reset();

        // Set base thresholds based on damage level
        float d_max = std::max(profile.max_damage_5prime, profile.max_damage_3prime);
        if (d_max > 0.30f) {
            // High damage: more aggressive correction
            thresholds_.base_correction_threshold = 0.25f;
            thresholds_.max_corrections_per_end = 4;
        } else if (d_max > 0.10f) {
            // Moderate damage: balanced
            thresholds_.base_correction_threshold = 0.35f;
            thresholds_.max_corrections_per_end = 3;
        } else {
            // Low damage: conservative
            thresholds_.base_correction_threshold = 0.45f;
            thresholds_.max_corrections_per_end = 2;
        }
        thresholds_.reset();
    }

    // Update with terminal k-mer stats from a read
    void update_terminal_stats(const std::string& seq) {
        if (seq.length() < 30) return;

        size_t len = seq.length();

        // 5' terminal (positions 0-5)
        for (size_t i = 0; i < 6 && i < len; ++i) {
            char b = std::toupper(static_cast<unsigned char>(seq[i]));
            if (b == 'T') terminal_stats_.terminal_5prime_t[i]++;
            else if (b == 'C') terminal_stats_.terminal_5prime_c[i]++;
        }

        // 3' terminal (positions 0-5 from end)
        for (size_t i = 0; i < 6 && i < len; ++i) {
            char b = std::toupper(static_cast<unsigned char>(seq[len - 1 - i]));
            if (b == 'A') terminal_stats_.terminal_3prime_a[i]++;
            else if (b == 'G') terminal_stats_.terminal_3prime_g[i]++;
        }

        // Interior (positions 15-25)
        for (size_t i = 15; i < 25 && i < len; ++i) {
            char b = std::toupper(static_cast<unsigned char>(seq[i]));
            if (b == 'T') terminal_stats_.interior_t++;
            else if (b == 'C') terminal_stats_.interior_c++;
            else if (b == 'A') terminal_stats_.interior_a++;
            else if (b == 'G') terminal_stats_.interior_g++;
        }

        terminal_stats_.n_reads++;
    }

    // Check if this read should have corrections applied
    bool should_correct_read(const std::string& seq, float damage_signal, int frame = -1) const {
        if (!thresholds_.corrections_enabled) return false;

        // Relaxed length guard: allow down to 20bp for high-prob ancient reads
        size_t min_len = (damage_signal > 0.7f) ? 20 : 24;
        if (seq.length() < min_len) return false;

        // Sample-aware damage threshold:
        // - If Channel B validated real damage: lower threshold (0.15) since we're confident
        // - If Channel B identified artifact: block all corrections
        // - Otherwise: moderate threshold (0.25)
        float damage_threshold = 0.25f;
        if (sample_profile_.damage_validated) {
            damage_threshold = 0.15f;  // Confident sample has real damage
        } else if (sample_profile_.damage_artifact) {
            return false;  // Sample damage is artifact, don't correct
        }

        if (damage_signal < damage_threshold) return false;

        // Compute per-read damage LLR (terminal-based)
        auto llr = ReadDamageLLR::compute(seq, sample_profile_);
        if (llr.is_low_signal()) return false;

        // Check ensemble hexamer LLR for coding quality validation
        // High coding potential means corrections are more likely to be valid
        if (seq.length() >= 30) {
            float ensemble_llr = get_ensemble_hexamer_llr(seq, frame);
            // Require minimum coding signal (ensemble LLR > -2 means better than random)
            if (ensemble_llr < -2.0f) {
                // Poor coding signal - be more conservative with corrections
                return llr.combined_llr > thresholds_.current_damage_llr_threshold * 1.5f;
            }
        }

        // Check if terminal pattern is flat (composition bias)
        if (terminal_stats_.n_reads > 1000 && terminal_stats_.is_flat_enrichment()) {
             #ifdef VERBOSE_DEBUG
             if (llr.combined_llr <= thresholds_.current_damage_llr_threshold) {
                // Log only if rejecting due to flat enrichment + low LLR
                // std::cerr << "REJECTED: flat enrichment and low LLR" << std::endl;
             }
             #endif
            // Only correct if this specific read shows strong damage signal
            return llr.combined_llr > thresholds_.current_damage_llr_threshold;
        }

        return true;
    }

    // Get maximum corrections allowed for a read
    int get_max_corrections(const std::string& seq) const {
        // Shorter reads get fewer corrections
        if (seq.length() < 40) return std::min(2, thresholds_.current_max_corrections);
        if (seq.length() < 60) return std::min(3, thresholds_.current_max_corrections);
        return thresholds_.current_max_corrections;
    }

    // Get correction threshold for a position
    float get_correction_threshold(size_t pos, size_t seq_len, int frame) const {
        float base = thresholds_.current_correction_threshold;

        // Frame-aware adjustment
        if (frame >= 0) {
            // Safe calculation avoiding underflow when pos < frame
            size_t frame_u = static_cast<size_t>(frame);
            size_t offset = (pos >= frame_u) ? (pos - frame_u) : (pos + 300 - frame_u);
            int codon_pos = static_cast<int>(offset % 3);
            if (codon_pos == 2) base *= 0.8f;  // Wobble: more aggressive
            else if (codon_pos == 0) base *= 1.15f;  // First pos: conservative
        }

        // Position-based adjustment (more confident at terminus)
        if (pos < 3 || pos >= seq_len - 3) {
            base *= 0.9f;  // Very terminal: slightly more aggressive
        }

        return base;
    }

    // Update batch statistics
    void update_batch_stats(size_t n_reads, size_t n_corrected,
                           size_t total_corrections, size_t stop_restorations,
                           double coding_before, double coding_after,
                           double hexamer_before, double hexamer_after) {
        batch_stats_.n_reads += n_reads;
        batch_stats_.n_corrected += n_corrected;
        batch_stats_.total_corrections += total_corrections;
        batch_stats_.stop_restorations += stop_restorations;
        batch_stats_.total_coding_score_before += coding_before;
        batch_stats_.total_coding_score_after += coding_after;
        batch_stats_.total_hexamer_llr_before += hexamer_before;
        batch_stats_.total_hexamer_llr_after += hexamer_after;
    }

    // Update per-domain hexamer stats for ensemble weight calibration
    void update_domain_stats(const std::array<double, 8>& domain_sums,
                             const std::array<size_t, 8>& domain_counts) {
        for (size_t i = 0; i < 8; ++i) {
            batch_stats_.domain_hexamer_sums[i] += domain_sums[i];
            batch_stats_.domain_counts[i] += domain_counts[i];
        }
    }

    // Called after each batch to adjust thresholds
    void calibrate_after_batch() {
        if (batch_stats_.n_reads < 1000) return;

        // Check if corrections are improving quality
        if (!batch_stats_.corrections_improving()) {
            thresholds_.tighten();
        }

        // Check if correction rate is too high (suspicious)
        if (batch_stats_.correction_rate() > 0.5f) {
            thresholds_.tighten();
        }

        // Check decay monitor
        if (decay_monitor_.n_samples > 1000 && decay_monitor_.is_decay_collapsing()) {
            thresholds_.tighten();
            // If decay completely collapsed, disable corrections
            auto [mean_5, mean_3] = decay_monitor_.mean_decay();
            if (mean_5 > 0.0f && mean_3 > 0.0f) {
                thresholds_.disable_corrections();
            }
        }

        // Update ensemble weights
        ensemble_weights_.update_from_batch(batch_stats_);

        // Reset batch stats for next batch
        batch_stats_.reset();
    }

    // Get damage posterior from terminal k-mer mixture model
    float get_damage_posterior() const {
        return terminal_stats_.damage_posterior();
    }

    // Check if we should use conservative mode (flat/uncertain ensemble)
    bool use_conservative_mode() const {
        return ensemble_weights_.is_uncertain() ||
               (terminal_stats_.n_reads > 1000 && terminal_stats_.is_flat_enrichment());
    }

    // Get ensemble-weighted hexamer LLR (string version)
    float get_ensemble_hexamer_llr(const std::string& seq, int frame) const;

    // Get ensemble-weighted hexamer LLR (const char* version - avoids allocation)
    float get_ensemble_hexamer_llr(const char* seq, size_t len, int frame) const;

    // Validate a proposed correction
    struct CorrectionValidation {
        bool valid = true;
        float quality_delta = 0.0f;
        bool creates_stop = false;
        bool worsens_hexamer = false;
    };

    CorrectionValidation validate_correction(
        const std::string& original_seq,
        const std::string& corrected_seq,
        size_t correction_pos,
        int frame) const;

    // Accessors
    const AdaptiveThresholds& thresholds() const { return thresholds_; }
    const TerminalKmerStats& terminal_stats() const { return terminal_stats_; }
    const BatchStats& batch_stats() const { return batch_stats_; }
    const EnsembleWeights& ensemble_weights() const { return ensemble_weights_; }
    const DecayMonitor& decay_monitor() const { return decay_monitor_; }

    // Add decay sample from a read
    void add_decay_sample(float decay_5, float decay_3) {
        decay_monitor_.add_sample(decay_5, decay_3);
    }

    // Merge batch terminal stats into calibrator's cumulative stats
    void merge_terminal_stats(const TerminalKmerStats& batch_stats) {
        for (size_t i = 0; i < 6; ++i) {
            terminal_stats_.terminal_5prime_t[i] += batch_stats.terminal_5prime_t[i];
            terminal_stats_.terminal_5prime_c[i] += batch_stats.terminal_5prime_c[i];
            terminal_stats_.terminal_3prime_a[i] += batch_stats.terminal_3prime_a[i];
            terminal_stats_.terminal_3prime_g[i] += batch_stats.terminal_3prime_g[i];
        }
        terminal_stats_.interior_t += batch_stats.interior_t;
        terminal_stats_.interior_c += batch_stats.interior_c;
        terminal_stats_.interior_a += batch_stats.interior_a;
        terminal_stats_.interior_g += batch_stats.interior_g;
        terminal_stats_.n_reads += batch_stats.n_reads;
    }

private:
    SampleDamageProfile sample_profile_;
    TerminalKmerStats terminal_stats_;
    BatchStats batch_stats_;
    AdaptiveThresholds thresholds_;
    EnsembleWeights ensemble_weights_;
    DecayMonitor decay_monitor_;
};

} // namespace dart
