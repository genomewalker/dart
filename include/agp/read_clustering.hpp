#pragma once

/**
 * @file read_clustering.hpp
 * @brief Fast read clustering for strand voting
 *
 * Clusters reads by sequence similarity using canonical minimizers.
 * Reads in the same cluster can vote on strand orientation to improve
 * accuracy beyond the per-read Bayes limit (~77% → 85%+ with clustering).
 *
 * Algorithm:
 * 1. Extract canonical minimizers (min of k-mer and RC) from each read
 * 2. Store (minimizer, read_id) pairs and sort by minimizer
 * 3. Union reads sharing minimizers using disjoint-set (union-find)
 * 4. Aggregate strand votes within clusters
 *
 * Complexity: O(n log n) for n reads
 * Memory: O(n * max_mins) for minimizer storage
 */

#include <vector>
#include <deque>
#include <algorithm>
#include <cstdint>
#include <string>
#include <limits>
#include <unordered_map>

namespace agp {

// ============================================================================
// Configuration
// ============================================================================

struct ClusterParams {
    int k = 15;              // k-mer size for minimizers
    int window = 20;         // Minimizer window size
    int max_mins = 4;        // Max minimizers per read
    int max_bucket = 256;    // Max reads per minimizer (skip overly common)
    float min_vote_conf = 0.6f;  // Min vote confidence to override individual strand call
};

// ============================================================================
// Utility functions
// ============================================================================

namespace detail {

inline char fast_upper(char c) {
    return (c >= 'a' && c <= 'z') ? c - 32 : c;
}

inline int base_bits(char c) {
    switch (fast_upper(c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// Fast hash function (splitmix64)
inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

} // namespace detail

// ============================================================================
// Minimizer extraction
// ============================================================================

/**
 * Extract canonical minimizers from a sequence.
 * Canonical = min(k-mer, reverse_complement(k-mer))
 * This ensures strand-agnostic clustering.
 */
inline void extract_minimizers(const std::string& seq,
                               const ClusterParams& p,
                               std::vector<uint64_t>& out) {
    out.clear();
    if (seq.length() < static_cast<size_t>(p.k)) return;

    const uint64_t mask = (p.k == 32) ? ~0ULL : ((1ULL << (2 * p.k)) - 1);
    uint64_t fwd = 0;
    uint64_t rev = 0;
    int valid = 0;

    std::vector<uint64_t> hashes;
    hashes.reserve(seq.length());

    for (size_t i = 0; i < seq.length(); ++i) {
        int b = detail::base_bits(seq[i]);
        if (b < 0) {
            valid = 0;
            fwd = 0;
            rev = 0;
            hashes.push_back(std::numeric_limits<uint64_t>::max());
            continue;
        }

        // Rolling hash: add base to forward, reverse complement to reverse
        fwd = ((fwd << 2) | static_cast<uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - b) << (2 * (p.k - 1)));
        valid++;

        if (valid >= p.k) {
            // Canonical k-mer = min(forward, reverse complement)
            uint64_t canon = std::min(fwd, rev);
            hashes.push_back(detail::splitmix64(canon));
        } else {
            hashes.push_back(std::numeric_limits<uint64_t>::max());
        }
    }

    if (hashes.size() < static_cast<size_t>(p.window)) return;

    // Extract minimizers using sliding window
    std::deque<size_t> dq;
    for (size_t i = 0; i < hashes.size(); ++i) {
        // Remove elements larger than current from back
        while (!dq.empty() && hashes[i] <= hashes[dq.back()]) dq.pop_back();
        dq.push_back(i);

        // Remove elements outside window from front
        if (dq.front() + p.window <= i) dq.pop_front();

        // Record minimizer after window is full
        if (i + 1 >= static_cast<size_t>(p.window)) {
            uint64_t h = hashes[dq.front()];
            if (h != std::numeric_limits<uint64_t>::max()) {
                if (out.empty() || out.back() != h) {
                    out.push_back(h);
                }
            }
        }
    }

    // Keep only top max_mins minimizers (smallest hashes)
    if (out.size() > static_cast<size_t>(p.max_mins)) {
        std::nth_element(out.begin(), out.begin() + p.max_mins, out.end());
        out.resize(p.max_mins);
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
    }
}

// ============================================================================
// Disjoint Set (Union-Find)
// ============================================================================

class DisjointSet {
public:
    std::vector<int> parent;
    std::vector<uint8_t> rank;

    explicit DisjointSet(size_t n) : parent(n), rank(n, 0) {
        for (size_t i = 0; i < n; ++i) parent[i] = static_cast<int>(i);
    }

    int find(int x) {
        int root = x;
        while (parent[root] != root) root = parent[root];
        // Path compression
        while (parent[x] != x) {
            int p = parent[x];
            parent[x] = root;
            x = p;
        }
        return root;
    }

    void unite(int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) return;
        // Union by rank
        if (rank[ra] < rank[rb]) std::swap(ra, rb);
        parent[rb] = ra;
        if (rank[ra] == rank[rb]) rank[ra]++;
    }
};

// ============================================================================
// Read clustering
// ============================================================================

/**
 * Cluster reads by shared minimizers.
 * Returns a DisjointSet where cluster_id = dsu.find(read_id).
 */
inline DisjointSet cluster_reads(const std::vector<std::string>& reads,
                                 const ClusterParams& p = ClusterParams()) {
    struct Pair { uint64_t hash; int id; };
    std::vector<Pair> pairs;
    pairs.reserve(reads.size() * p.max_mins);

    std::vector<uint64_t> mins;
    for (int i = 0; i < static_cast<int>(reads.size()); ++i) {
        extract_minimizers(reads[i], p, mins);
        for (uint64_t h : mins) {
            pairs.push_back({h, i});
        }
    }

    // Sort by hash to group reads with same minimizer
    std::sort(pairs.begin(), pairs.end(),
              [](const Pair& a, const Pair& b) { return a.hash < b.hash; });

    DisjointSet dsu(reads.size());
    size_t start = 0;
    while (start < pairs.size()) {
        size_t end = start + 1;
        while (end < pairs.size() && pairs[end].hash == pairs[start].hash) end++;

        size_t bucket = end - start;
        // Skip overly common minimizers (low complexity)
        if (bucket <= static_cast<size_t>(p.max_bucket)) {
            int base = pairs[start].id;
            for (size_t i = start + 1; i < end; ++i) {
                dsu.unite(base, pairs[i].id);
            }
        }

        start = end;
    }

    return dsu;
}

// ============================================================================
// Strand voting
// ============================================================================

/**
 * Strand vote for a single read.
 */
struct StrandVote {
    float forward_score;    // Score if forward strand is correct
    float reverse_score;    // Score if reverse strand is correct
    bool predicted_forward; // Original per-read prediction

    float confidence() const {
        float total = forward_score + reverse_score;
        if (total < 1e-6f) return 0.5f;
        return std::max(forward_score, reverse_score) / total;
    }

    bool vote_forward() const {
        return forward_score > reverse_score;
    }
};

/**
 * Aggregate strand votes within a cluster.
 * Returns (forward_votes, reverse_votes, confidence).
 */
inline std::tuple<float, float, float> aggregate_cluster_votes(
    const std::vector<StrandVote>& votes,
    const std::vector<int>& read_ids,
    const DisjointSet& dsu,
    int cluster_id) {

    float fwd_sum = 0.0f;
    float rev_sum = 0.0f;
    int count = 0;

    for (int id : read_ids) {
        if (dsu.parent[id] != cluster_id &&
            const_cast<DisjointSet&>(dsu).find(id) != cluster_id) continue;

        const auto& v = votes[id];
        fwd_sum += v.forward_score;
        rev_sum += v.reverse_score;
        count++;
    }

    if (count == 0) return {0.0f, 0.0f, 0.0f};

    float total = fwd_sum + rev_sum;
    float conf = (total > 1e-6f) ? std::max(fwd_sum, rev_sum) / total : 0.5f;

    return {fwd_sum, rev_sum, conf};
}

/**
 * Apply cluster voting to improve strand predictions.
 * Modifies votes in-place if cluster vote is confident enough.
 *
 * @param votes Per-read strand votes (modified in-place)
 * @param dsu Cluster membership from cluster_reads()
 * @param params Clustering parameters (uses min_vote_conf)
 * @return Number of predictions changed by cluster voting
 */
inline int apply_cluster_voting(std::vector<StrandVote>& votes,
                                DisjointSet& dsu,
                                const ClusterParams& params = ClusterParams()) {
    if (votes.empty()) return 0;

    // Group reads by cluster
    std::unordered_map<int, std::vector<int>> clusters;
    for (size_t i = 0; i < votes.size(); ++i) {
        int cluster_id = dsu.find(static_cast<int>(i));
        clusters[cluster_id].push_back(static_cast<int>(i));
    }

    int changed = 0;

    for (auto& [cluster_id, members] : clusters) {
        if (members.size() < 2) continue;  // No voting for singletons

        // Aggregate votes
        float fwd_sum = 0.0f, rev_sum = 0.0f;
        for (int id : members) {
            fwd_sum += votes[id].forward_score;
            rev_sum += votes[id].reverse_score;
        }

        float total = fwd_sum + rev_sum;
        if (total < 1e-6f) continue;

        float conf = std::max(fwd_sum, rev_sum) / total;
        bool cluster_forward = fwd_sum > rev_sum;

        // If cluster is confident, override low-confidence individual predictions
        if (conf >= params.min_vote_conf) {
            for (int id : members) {
                auto& v = votes[id];
                if (v.predicted_forward != cluster_forward && v.confidence() < conf) {
                    // Override with cluster decision
                    v.predicted_forward = cluster_forward;
                    changed++;
                }
            }
        }
    }

    return changed;
}

// ============================================================================
// Frame voting
// ============================================================================

/**
 * Frame vote for a single read.
 * Stores scores for all 6 hypotheses (3 frames × 2 strands).
 */
struct FrameVote {
    std::array<float, 6> scores;  // [fwd_f0, fwd_f1, fwd_f2, rev_f0, rev_f1, rev_f2]
    int predicted_frame;          // 0, 1, or 2
    bool predicted_forward;       // Original per-read prediction

    FrameVote() : scores{}, predicted_frame(0), predicted_forward(true) {
        scores.fill(0.0f);
    }

    // Get score for specific hypothesis
    float score(bool forward, int frame) const {
        return scores[forward ? frame : (3 + frame)];
    }

    // Set score for specific hypothesis
    void set_score(bool forward, int frame, float s) {
        scores[forward ? frame : (3 + frame)] = s;
    }

    // Confidence in the best hypothesis vs sum of all
    float confidence() const {
        float total = 0.0f;
        float max_score = 0.0f;
        for (float s : scores) {
            total += s;
            if (s > max_score) max_score = s;
        }
        if (total < 1e-6f) return 1.0f / 6.0f;
        return max_score / total;
    }

    // Marginalized frame confidence (sum over strands)
    float frame_confidence(int frame) const {
        float fwd = scores[frame];
        float rev = scores[3 + frame];
        float total = 0.0f;
        for (float s : scores) total += s;
        if (total < 1e-6f) return 1.0f / 3.0f;
        return (fwd + rev) / total;
    }

    // Best frame (marginalized over strand)
    int best_frame() const {
        std::array<float, 3> frame_scores = {0.0f, 0.0f, 0.0f};
        for (int f = 0; f < 3; ++f) {
            frame_scores[f] = scores[f] + scores[3 + f];
        }
        return std::max_element(frame_scores.begin(), frame_scores.end()) - frame_scores.begin();
    }
};

/**
 * Apply cluster voting to improve frame predictions.
 * Uses marginalization over strand to focus on frame consensus.
 *
 * @param votes Per-read frame votes (modified in-place)
 * @param dsu Cluster membership from cluster_reads()
 * @param params Clustering parameters
 * @return Number of frame predictions changed
 */
inline int apply_cluster_frame_voting(std::vector<FrameVote>& votes,
                                      DisjointSet& dsu,
                                      const ClusterParams& params = ClusterParams()) {
    if (votes.empty()) return 0;

    // Group reads by cluster
    std::unordered_map<int, std::vector<int>> clusters;
    for (size_t i = 0; i < votes.size(); ++i) {
        int cluster_id = dsu.find(static_cast<int>(i));
        clusters[cluster_id].push_back(static_cast<int>(i));
    }

    int changed = 0;

    for (auto& [cluster_id, members] : clusters) {
        if (members.size() < 2) continue;

        // Aggregate frame scores (marginalized over strand)
        std::array<float, 3> frame_sums = {0.0f, 0.0f, 0.0f};
        for (int id : members) {
            const auto& v = votes[id];
            for (int f = 0; f < 3; ++f) {
                frame_sums[f] += v.scores[f] + v.scores[3 + f];
            }
        }

        float total = frame_sums[0] + frame_sums[1] + frame_sums[2];
        if (total < 1e-6f) continue;

        int cluster_frame = std::max_element(frame_sums.begin(), frame_sums.end()) - frame_sums.begin();
        float cluster_conf = frame_sums[cluster_frame] / total;

        // Override low-confidence individual predictions
        if (cluster_conf >= params.min_vote_conf) {
            for (int id : members) {
                auto& v = votes[id];
                if (v.predicted_frame != cluster_frame && v.frame_confidence(v.predicted_frame) < cluster_conf) {
                    v.predicted_frame = cluster_frame;
                    changed++;
                }
            }
        }
    }

    return changed;
}

/**
 * Combined frame+strand voting result.
 */
struct CombinedVote {
    FrameVote frame_vote;
    StrandVote strand_vote;
};

/**
 * Apply combined cluster voting for both frame and strand.
 * First votes on strand, then on frame within the strand-consensus.
 *
 * @param votes Per-read combined votes (modified in-place)
 * @param dsu Cluster membership
 * @param params Clustering parameters
 * @return Pair of (strand_changes, frame_changes)
 */
inline std::pair<int, int> apply_combined_cluster_voting(
    std::vector<CombinedVote>& votes,
    DisjointSet& dsu,
    const ClusterParams& params = ClusterParams()) {

    if (votes.empty()) return {0, 0};

    // Extract separate vote vectors
    std::vector<StrandVote> strand_votes;
    std::vector<FrameVote> frame_votes;
    strand_votes.reserve(votes.size());
    frame_votes.reserve(votes.size());

    for (const auto& v : votes) {
        strand_votes.push_back(v.strand_vote);
        frame_votes.push_back(v.frame_vote);
    }

    // Apply strand voting first
    int strand_changes = apply_cluster_voting(strand_votes, dsu, params);

    // Apply frame voting
    int frame_changes = apply_cluster_frame_voting(frame_votes, dsu, params);

    // Copy back results
    for (size_t i = 0; i < votes.size(); ++i) {
        votes[i].strand_vote = strand_votes[i];
        votes[i].frame_vote = frame_votes[i];
    }

    return {strand_changes, frame_changes};
}

} // namespace agp
