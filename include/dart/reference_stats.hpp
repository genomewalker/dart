#pragma once
// Per-reference coverage and damage statistics for EM reassignment
//
// Tracks breadth, depth, evenness, and EM-derived damage decomposition
// for each reference protein. Thread-safe accumulation via per-thread
// collectors merged in finalize().

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iosfwd>
#include <mutex>
#include <numeric>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dart {

// Per-reference summary statistics (output of finalize)
struct ReferenceStats {
    std::string ref_id;
    uint32_t ref_length = 0;

    // Coverage metrics
    float breadth = 0.0f;          // nonzero_positions / ref_length
    float depth_mean = 0.0f;       // mean depth across covered positions
    float depth_std = 0.0f;        // std dev of depth
    float depth_evenness = 0.0f;   // 1 - (std/mean), measures uniformity

    // Positional diversity (better for short proteins than entropy)
    // High = reads start at different positions (authentic)
    // Low = all reads start at same position (spurious motif match)
    float start_diversity = 0.0f;  // n_unique_starts / min(n_reads, ref_len)
    float start_span = 0.0f;       // (max_start - min_start + 1) / ref_length
    float positional_score = 0.0f; // sqrt(diversity * span) [0-1]

    // Terminal coverage (detects middle-only spurious matches)
    // Real genes have reads at both ends, not just conserved middle
    float terminal_ratio = 0.0f;   // (5'+3' coverage) / (2*middle coverage)

    // Counts
    uint32_t total_covered_bases = 0;
    uint64_t total_depth = 0;
    uint32_t n_reads = 0;
    uint32_t n_unique_starts = 0;  // count of distinct read start positions

    // Quality metrics
    float avg_alignment_length = 0.0f;
    float avg_read_length = 0.0f;
    float avg_identity = 0.0f;

    // Damage decomposition (from EM)
    double n_effective = 0.0;      // sum gamma_{rt} (fractional assignment)
    double n_ancient = 0.0;        // sum gamma_{rt,ancient}
    double n_modern = 0.0;         // sum gamma_{rt,modern}
    float damage_enrichment = 0.0f; // log2(n_ancient/n_modern) - log2(pi/(1-pi))
};

// Filtering thresholds for reference selection
struct FilterThresholds {
    uint32_t min_reads = 3;
    float min_breadth = 0.10f;
    float min_depth = 0.5f;
    float min_identity = 0.5f;
    float min_enrichment = 0.0f;   // positive = ancient-enriched

    // Positional diversity filtering (for short proteins)
    // Filters spurious matches where all reads hit same position
    float min_positional_score = 0.0f;  // sqrt(diversity * span), 0-1
    float min_terminal_ratio = 0.0f;    // terminal vs middle coverage

    bool passes(const ReferenceStats& s) const noexcept {
        return s.n_reads >= min_reads
            && s.breadth >= min_breadth
            && s.depth_mean >= min_depth
            && s.avg_identity >= min_identity
            && s.damage_enrichment >= min_enrichment
            && s.positional_score >= min_positional_score
            && s.terminal_ratio >= min_terminal_ratio;
    }
};

// Streaming coverage accumulator for a single reference
// Tracks per-position weighted depth and quality sums
class CoverageAccumulator {
public:
    CoverageAccumulator() = default;

    explicit CoverageAccumulator(const std::string& ref_id, uint32_t ref_length)
        : ref_id_(ref_id)
        , ref_length_(ref_length)
        , coverage_(ref_length, 0.0f) {}

    void add_alignment(uint32_t start, uint32_t end, float weight,
                       float identity, uint32_t aln_len, uint32_t read_len) {
        uint32_t clamped_end = std::min(end, ref_length_);
        for (uint32_t p = start; p < clamped_end; ++p) {
            coverage_[p] += weight;
        }
        sum_identity_ += identity;
        sum_aln_len_ += aln_len;
        sum_read_len_ += read_len;
        ++n_reads_;

        // Track unique start positions for positional diversity
        unique_starts_.insert(start);
        if (n_reads_ == 1) {
            min_start_ = max_start_ = start;
        } else {
            min_start_ = std::min(min_start_, start);
            max_start_ = std::max(max_start_, start);
        }
    }

    void add_em_weights(double gamma, double gamma_ancient) {
        n_effective_ += gamma;
        n_ancient_ += gamma_ancient;
        n_modern_ += (gamma - gamma_ancient);
    }

    ReferenceStats finalize(float pi = 0.10f) const {
        ReferenceStats s;
        s.ref_id = ref_id_;
        s.ref_length = ref_length_;
        s.n_reads = n_reads_;

        if (ref_length_ == 0 || n_reads_ == 0) return s;

        // Coverage metrics
        uint32_t nonzero = 0;
        double sum = 0.0;
        double sum_sq = 0.0;
        for (uint32_t p = 0; p < ref_length_; ++p) {
            float d = coverage_[p];
            if (d > 0.0f) {
                ++nonzero;
                sum += d;
                sum_sq += static_cast<double>(d) * d;
            }
        }

        s.total_covered_bases = nonzero;
        s.total_depth = static_cast<uint64_t>(sum + 0.5);
        s.breadth = static_cast<float>(nonzero) / static_cast<float>(ref_length_);

        if (nonzero > 0) {
            s.depth_mean = static_cast<float>(sum / nonzero);
            double var = (sum_sq / nonzero) - (sum / nonzero) * (sum / nonzero);
            s.depth_std = static_cast<float>(std::sqrt(std::max(0.0, var)));
            s.depth_evenness = (s.depth_mean > 0.0f)
                ? 1.0f - s.depth_std / s.depth_mean
                : 0.0f;
        }

        // Positional diversity (better for short proteins than entropy)
        s.n_unique_starts = static_cast<uint32_t>(unique_starts_.size());
        if (n_reads_ > 0 && ref_length_ > 0) {
            // diversity: fraction of possible start positions used
            float max_possible = static_cast<float>(std::min(n_reads_, ref_length_));
            s.start_diversity = static_cast<float>(s.n_unique_starts) / max_possible;

            // span: how much of the reference do start positions cover
            uint32_t span_len = (max_start_ >= min_start_)
                ? (max_start_ - min_start_ + 1) : 1;
            s.start_span = static_cast<float>(span_len) / static_cast<float>(ref_length_);

            // Combined score: geometric mean (both must be high)
            s.positional_score = std::sqrt(s.start_diversity * s.start_span);
        }

        // Terminal coverage ratio
        // Real genes have reads covering both ends, not just conserved middle
        if (ref_length_ >= 6) {
            // 5' terminal: first 2 positions (or 10% of length)
            // 3' terminal: last 2 positions (or 10% of length)
            uint32_t term_size = std::max(2u, ref_length_ / 10);
            term_size = std::min(term_size, ref_length_ / 3);

            double term_5prime = 0.0, term_3prime = 0.0, middle = 0.0;
            for (uint32_t p = 0; p < ref_length_; ++p) {
                double d = coverage_[p];
                if (p < term_size) {
                    term_5prime += d;
                } else if (p >= ref_length_ - term_size) {
                    term_3prime += d;
                } else {
                    middle += d;
                }
            }

            // Normalize by region size
            double term_avg = (term_5prime + term_3prime) / (2.0 * term_size);
            uint32_t middle_size = ref_length_ - 2 * term_size;
            double middle_avg = (middle_size > 0) ? middle / middle_size : 0.0;

            // Ratio: high = good terminal coverage, low = middle-only (spurious)
            // Add small epsilon to avoid division by zero
            s.terminal_ratio = static_cast<float>(
                (term_avg + 0.01) / (middle_avg + 0.01));
        } else {
            s.terminal_ratio = 1.0f;  // Too short to measure
        }

        // Quality metrics
        float nr = static_cast<float>(n_reads_);
        s.avg_identity = sum_identity_ / nr;
        s.avg_alignment_length = static_cast<float>(sum_aln_len_) / nr;
        s.avg_read_length = static_cast<float>(sum_read_len_) / nr;

        // EM damage decomposition
        s.n_effective = n_effective_;
        s.n_ancient = n_ancient_;
        s.n_modern = n_modern_;

        // Damage enrichment: log2(ancient/modern) - log2(pi/(1-pi))
        // Positive = more ancient than expected by prior
        if (n_ancient_ > 0.0 && n_modern_ > 0.0 && pi > 0.0f && pi < 1.0f) {
            double obs_ratio = n_ancient_ / n_modern_;
            double prior_ratio = static_cast<double>(pi) / (1.0 - static_cast<double>(pi));
            s.damage_enrichment = static_cast<float>(
                std::log2(obs_ratio) - std::log2(prior_ratio));
        }

        return s;
    }

private:
    std::string ref_id_;
    uint32_t ref_length_ = 0;
    std::vector<float> coverage_;

    uint32_t n_reads_ = 0;
    double sum_identity_ = 0.0;
    uint64_t sum_aln_len_ = 0;
    uint64_t sum_read_len_ = 0;

    double n_effective_ = 0.0;
    double n_ancient_ = 0.0;
    double n_modern_ = 0.0;

    // Positional diversity tracking
    std::unordered_set<uint32_t> unique_starts_;
    uint32_t min_start_ = 0;
    uint32_t max_start_ = 0;
};

// Thread-safe collector for parallel read processing
// Each thread calls update(); finalize() merges and computes stats.
class ReferenceStatsCollector {
public:
    void reserve(size_t n_refs) {
        accumulators_.reserve(n_refs);
    }

    // Register a reference before any updates (not thread-safe, call from single thread)
    void add_reference(const std::string& ref_id, uint32_t ref_length) {
        accumulators_.emplace(ref_id, CoverageAccumulator(ref_id, ref_length));
    }

    // Thread-safe update: add alignment hit to a reference
    void update(const std::string& ref_id, uint32_t ref_len,
                uint32_t aln_start, uint32_t aln_end,
                float gamma, float gamma_ancient,
                float identity, uint32_t aln_len = 0, uint32_t read_len = 0) {
        std::lock_guard<std::mutex> lock(mu_);
        auto it = accumulators_.find(ref_id);
        if (it == accumulators_.end()) {
            it = accumulators_.emplace(ref_id,
                CoverageAccumulator(ref_id, ref_len)).first;
        }
        it->second.add_alignment(aln_start, aln_end, gamma,
                                 identity, aln_len, read_len);
        it->second.add_em_weights(gamma, gamma_ancient);
    }

    // Finalize all references and return stats (call after all updates)
    std::vector<ReferenceStats> finalize(float pi = 0.10f) const {
        std::vector<ReferenceStats> out;
        out.reserve(accumulators_.size());
        for (const auto& [id, acc] : accumulators_) {
            out.push_back(acc.finalize(pi));
        }
        return out;
    }

    // Finalize with filtering
    std::vector<ReferenceStats> finalize_filtered(
            float pi, const FilterThresholds& thresh) const {
        std::vector<ReferenceStats> out;
        out.reserve(accumulators_.size());
        for (const auto& [id, acc] : accumulators_) {
            auto s = acc.finalize(pi);
            if (thresh.passes(s)) {
                out.push_back(std::move(s));
            }
        }
        return out;
    }

    size_t size() const {
        std::lock_guard<std::mutex> lock(mu_);
        return accumulators_.size();
    }

private:
    mutable std::mutex mu_;
    std::unordered_map<std::string, CoverageAccumulator> accumulators_;
};

// Output formatting (implemented in reference_stats.cpp)
void write_reference_stats_tsv(std::ostream& out,
                               const std::vector<ReferenceStats>& stats);
void write_reference_stats_json(std::ostream& out,
                                const std::vector<ReferenceStats>& stats,
                                float pi);

} // namespace dart
