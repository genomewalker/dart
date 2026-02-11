// FDR estimation using target-decoy approach
// Provides sample-specific FDR calibration for gene predictions

#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace agp {

struct FDRResult {
    float threshold;        // Score threshold used
    float fdr;              // Estimated FDR at this threshold
    size_t n_targets;       // Number of targets passing threshold
    size_t n_decoys;        // Number of decoys passing threshold
    float sensitivity;      // Fraction of targets passing
};

struct GeneWithScore {
    float score;
    size_t length;
    float damage_level;
    bool is_decoy;
};

class FDREstimator {
public:
    // Estimate FDR at a given score threshold using target-decoy competition
    // FDR = n_decoys / (n_targets + n_decoys) at threshold
    static float estimate_fdr(
        float threshold,
        const std::vector<float>& target_scores,
        const std::vector<float>& decoy_scores
    ) {
        size_t n_targets = 0, n_decoys = 0;

        for (float s : target_scores) {
            if (s >= threshold) ++n_targets;
        }
        for (float s : decoy_scores) {
            if (s >= threshold) ++n_decoys;
        }

        if (n_targets + n_decoys == 0) return 0.0f;

        // Conservative FDR estimate
        return static_cast<float>(n_decoys) / static_cast<float>(n_targets + n_decoys);
    }

    // Find threshold that achieves target FDR (e.g., 0.01 for 1% FDR)
    // Auto-determines score range from data if not specified
    static float find_threshold_for_fdr(
        float target_fdr,
        const std::vector<float>& target_scores,
        const std::vector<float>& decoy_scores,
        float min_threshold = 0.0f,
        float max_threshold = -1.0f,  // -1 = auto from data
        float step = 0.01f
    ) {
        // Auto-determine max threshold from actual score range
        if (max_threshold < 0.0f) {
            max_threshold = 0.0f;
            for (float s : target_scores) {
                if (s > max_threshold) max_threshold = s;
            }
            for (float s : decoy_scores) {
                if (s > max_threshold) max_threshold = s;
            }
            max_threshold = std::ceil(max_threshold * 100.0f) / 100.0f;  // Round up
        }

        float best_threshold = max_threshold;

        for (float t = max_threshold; t >= min_threshold; t -= step) {
            float fdr = estimate_fdr(t, target_scores, decoy_scores);
            if (fdr <= target_fdr) {
                best_threshold = t;
            } else {
                break;  // FDR exceeded, stop lowering threshold
            }
        }

        return best_threshold;
    }

    // Calculate q-values (minimum FDR at which each score would pass)
    static std::vector<float> calculate_q_values(
        const std::vector<float>& target_scores,
        const std::vector<float>& decoy_scores
    ) {
        // Combine and sort all scores
        std::vector<std::pair<float, bool>> all_scores;
        all_scores.reserve(target_scores.size() + decoy_scores.size());

        for (float s : target_scores) {
            all_scores.push_back({s, false});  // false = target
        }
        for (float s : decoy_scores) {
            all_scores.push_back({s, true});   // true = decoy
        }

        // Sort descending by score
        std::sort(all_scores.begin(), all_scores.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        // Calculate FDR at each score threshold
        std::vector<float> fdrs;
        fdrs.reserve(all_scores.size());

        size_t cum_targets = 0, cum_decoys = 0;
        for (const auto& [score, is_decoy] : all_scores) {
            if (is_decoy) {
                ++cum_decoys;
            } else {
                ++cum_targets;
            }
            float fdr = (cum_targets + cum_decoys > 0) ?
                static_cast<float>(cum_decoys) / (cum_targets + cum_decoys) : 0.0f;
            fdrs.push_back(fdr);
        }

        // Convert FDRs to q-values (monotonically increasing from bottom)
        std::vector<float> q_values(fdrs.size());
        float min_fdr = 1.0f;
        for (int i = fdrs.size() - 1; i >= 0; --i) {
            min_fdr = std::min(min_fdr, fdrs[i]);
            q_values[i] = min_fdr;
        }

        // Map q-values back to target scores only
        std::vector<float> target_q_values;
        target_q_values.reserve(target_scores.size());

        // Create index mapping
        std::vector<size_t> target_indices;
        for (size_t i = 0; i < all_scores.size(); ++i) {
            if (!all_scores[i].second) {
                target_q_values.push_back(q_values[i]);
            }
        }

        return target_q_values;
    }

    // Calculate Posterior Error Probability (PEP) using mixture model
    // PEP = P(incorrect | score) = P(score | incorrect) * P(incorrect) / P(score)
    static std::vector<float> calculate_pep(
        const std::vector<float>& target_scores,
        const std::vector<float>& decoy_scores
    ) {
        if (target_scores.empty() || decoy_scores.empty()) {
            return std::vector<float>(target_scores.size(), 0.5f);
        }

        // Estimate null distribution (incorrect) from decoys
        float null_mean = 0.0f, null_var = 0.0f;
        for (float s : decoy_scores) null_mean += s;
        null_mean /= decoy_scores.size();

        for (float s : decoy_scores) {
            float diff = s - null_mean;
            null_var += diff * diff;
        }
        null_var /= decoy_scores.size();
        float null_sd = std::sqrt(null_var + 1e-6f);

        // Estimate prior P(incorrect) from proportion of decoys in top hits
        size_t top_n = std::min(decoy_scores.size(), target_scores.size() / 10);
        if (top_n < 10) top_n = std::min(decoy_scores.size(), size_t(100));

        std::vector<float> all_scores = target_scores;
        all_scores.insert(all_scores.end(), decoy_scores.begin(), decoy_scores.end());
        std::sort(all_scores.begin(), all_scores.end(), std::greater<float>());

        float score_threshold = (top_n > 0 && top_n <= all_scores.size()) ?
            all_scores[top_n - 1] : 0.0f;

        size_t decoys_in_top = 0;
        for (float s : decoy_scores) {
            if (s >= score_threshold) ++decoys_in_top;
        }
        float pi0 = static_cast<float>(decoys_in_top) / top_n;
        pi0 = std::max(0.01f, std::min(0.99f, pi0));  // Bound pi0

        // Calculate PEP for each target
        std::vector<float> pep_values;
        pep_values.reserve(target_scores.size());

        for (float score : target_scores) {
            // P(score | incorrect) ~ N(null_mean, null_sd)
            float z = (score - null_mean) / null_sd;
            float p_score_null = std::exp(-0.5f * z * z) / (null_sd * 2.507f);

            // P(score | correct) ~ empirical from high-scoring targets
            // Approximation: use exponential tail
            float p_score_correct = std::exp(2.0f * (score - 1.0f));
            p_score_correct = std::max(1e-10f, p_score_correct);

            // Bayes' rule
            float p_score = pi0 * p_score_null + (1.0f - pi0) * p_score_correct;
            float pep = (p_score > 0) ? (pi0 * p_score_null) / p_score : 0.5f;
            pep = std::max(0.0f, std::min(1.0f, pep));

            pep_values.push_back(pep);
        }

        return pep_values;
    }

    // Length-stratified FDR estimation
    static std::vector<FDRResult> estimate_fdr_by_length(
        const std::vector<GeneWithScore>& genes,
        const std::vector<size_t>& length_bins = {50, 100, 150, 200, 300, 500}
    ) {
        std::vector<FDRResult> results;

        for (size_t i = 0; i <= length_bins.size(); ++i) {
            size_t min_len = (i == 0) ? 0 : length_bins[i - 1];
            size_t max_len = (i < length_bins.size()) ? length_bins[i] : SIZE_MAX;

            std::vector<float> bin_targets, bin_decoys;
            for (const auto& g : genes) {
                if (g.length >= min_len && g.length < max_len) {
                    if (g.is_decoy) {
                        bin_decoys.push_back(g.score);
                    } else {
                        bin_targets.push_back(g.score);
                    }
                }
            }

            if (!bin_targets.empty() && !bin_decoys.empty()) {
                float threshold = find_threshold_for_fdr(0.01f, bin_targets, bin_decoys);
                float fdr = estimate_fdr(threshold, bin_targets, bin_decoys);

                size_t n_passing = 0;
                for (float s : bin_targets) {
                    if (s >= threshold) ++n_passing;
                }

                results.push_back({
                    threshold,
                    fdr,
                    n_passing,
                    0,  // Decoys passing at this threshold
                    static_cast<float>(n_passing) / bin_targets.size()
                });
            }
        }

        return results;
    }

    // Generate FDR calibration report
    struct FDRReport {
        size_t n_targets;
        size_t n_decoys;
        float empirical_null_mean;
        float empirical_null_sd;
        std::vector<std::pair<float, float>> fdr_curve;  // (threshold, fdr)
        float threshold_01;   // Threshold for 1% FDR
        float threshold_05;   // Threshold for 5% FDR
        float threshold_10;   // Threshold for 10% FDR
        size_t genes_01;      // Genes passing 1% FDR
        size_t genes_05;      // Genes passing 5% FDR
        size_t genes_10;      // Genes passing 10% FDR
    };

    static FDRReport generate_report(
        const std::vector<float>& target_scores,
        const std::vector<float>& decoy_scores
    ) {
        FDRReport report;
        report.n_targets = target_scores.size();
        report.n_decoys = decoy_scores.size();

        // Null distribution statistics
        float sum = 0.0f;
        for (float s : decoy_scores) sum += s;
        report.empirical_null_mean = sum / decoy_scores.size();

        float var = 0.0f;
        for (float s : decoy_scores) {
            float diff = s - report.empirical_null_mean;
            var += diff * diff;
        }
        report.empirical_null_sd = std::sqrt(var / decoy_scores.size());

        // Determine actual max score for FDR curve
        float max_score = 0.0f;
        for (float s : target_scores) {
            if (s > max_score) max_score = s;
        }
        for (float s : decoy_scores) {
            if (s > max_score) max_score = s;
        }
        max_score = std::ceil(max_score * 20.0f) / 20.0f;  // Round up to nearest 0.05

        // FDR curve at various thresholds (covering full score range)
        for (float t = 0.0f; t <= max_score; t += 0.05f) {
            float fdr = estimate_fdr(t, target_scores, decoy_scores);
            report.fdr_curve.push_back({t, fdr});
        }

        // Key thresholds (auto-determine max from data)
        report.threshold_01 = find_threshold_for_fdr(0.01f, target_scores, decoy_scores);
        report.threshold_05 = find_threshold_for_fdr(0.05f, target_scores, decoy_scores);
        report.threshold_10 = find_threshold_for_fdr(0.10f, target_scores, decoy_scores);

        // Count genes at each threshold
        report.genes_01 = report.genes_05 = report.genes_10 = 0;
        for (float s : target_scores) {
            if (s >= report.threshold_01) ++report.genes_01;
            if (s >= report.threshold_05) ++report.genes_05;
            if (s >= report.threshold_10) ++report.genes_10;
        }

        return report;
    }
};

}  // namespace agp
