#pragma once
// VTML-aware damage scoring for protein clustering
//
// Uses VTML20 substitution matrix to:
// 1. Generate damage-tolerant k-mer variants for clustering
// 2. Weight consensus building (high VTML = equivalent, low = damage)
// 3. Score damage probability at each mismatch position
//
// Key insight: VTML encodes biochemical similarity. Damage produces
// conservative substitutions (high VTML), so we can distinguish:
// - High VTML mismatch: natural variation OR damage (ambiguous)
// - Low VTML mismatch: more likely damage (unusual substitution)

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <cmath>

namespace agp {

// VTML20 substitution matrix (20x20 for standard amino acids)
// Indexed by AA code: A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9,
//                     M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19
// Score 127 = invalid pair
class VTML20 {
public:
    static constexpr int8_t INVALID = 127;

    // Convert AA letter to index (0-19), returns -1 for invalid
    static constexpr int aa_to_idx(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'D': case 'd': return 2;
            case 'E': case 'e': return 3;
            case 'F': case 'f': return 4;
            case 'G': case 'g': return 5;
            case 'H': case 'h': return 6;
            case 'I': case 'i': return 7;
            case 'K': case 'k': return 8;
            case 'L': case 'l': return 9;
            case 'M': case 'm': return 10;
            case 'N': case 'n': return 11;
            case 'P': case 'p': return 12;
            case 'Q': case 'q': return 13;
            case 'R': case 'r': return 14;
            case 'S': case 's': return 15;
            case 'T': case 't': return 16;
            case 'V': case 'v': return 17;
            case 'W': case 'w': return 18;
            case 'Y': case 'y': return 19;
            default: return -1;
        }
    }

    // Get VTML20 score for substitution (from → to)
    static int8_t score(char from, char to) {
        int i = aa_to_idx(from);
        int j = aa_to_idx(to);
        if (i < 0 || j < 0) return INVALID;
        return MATRIX[i][j];
    }

    // Is this a damage-consistent substitution?
    static bool is_damage_sub(char from, char to, bool is_5prime) {
        if (is_5prime) {
            // CT damage patterns
            return (from == 'R' && (to == 'W' || to == 'C')) ||
                   (from == 'Q' && to == '*') ||
                   (from == 'H' && to == 'Y') ||
                   (from == 'P' && (to == 'S' || to == 'L')) ||
                   (from == 'T' && (to == 'I' || to == 'M')) ||
                   (from == 'A' && to == 'V') ||
                   (from == 'S' && (to == 'F' || to == 'L'));
        } else {
            // GA damage patterns
            return (from == 'E' && to == 'K') ||
                   (from == 'D' && to == 'N') ||
                   (from == 'A' && to == 'T') ||
                   (from == 'G' && (to == 'R' || to == 'S' || to == 'E' || to == 'D')) ||
                   (from == 'V' && (to == 'I' || to == 'M')) ||
                   (from == 'R' && (to == 'K' || to == 'Q' || to == 'H')) ||
                   (from == 'S' && to == 'N') ||
                   (from == 'C' && to == 'Y') ||
                   (from == 'W' && to == '*');
        }
    }

    // Convert VTML score to probability that mismatch is damage
    // Low VTML = unusual sub = more likely damage
    // High VTML = common sub = could be natural variation
    static float score_to_damage_prob(int8_t vtml_score, float d_max) {
        if (vtml_score == INVALID) return 0.0f;

        // Logistic transform: P(damage) = d_max * sigmoid(-vtml_score / 4)
        // At VTML = -8: P ≈ d_max * 0.88 (high damage probability)
        // At VTML = 0:  P ≈ d_max * 0.50 (uncertain)
        // At VTML = +4: P ≈ d_max * 0.27 (lower probability)
        float sigmoid = 1.0f / (1.0f + std::exp(vtml_score / 4.0f));
        return d_max * sigmoid;
    }

private:
    // VTML20 matrix (from parsed file)
    static constexpr int8_t MATRIX[20][20] = {
        // A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y
        {  7,  -3,  -6,  -5,  -8,  -4,  -7,  -7,  -6,  -7,  -5,  -6,  -4,  -5,  -7,  -2,  -3,  -3,  -9,  -8}, // A
        { -3,  12, -14, -14, -13,  -7,  -6,  -5, -13, -12,  -4,  -8,  -9, -13,  -7,  -3,  -5,  -3, -15,  -4}, // C
        { -6, -14,   8,  -1, -16,  -6,  -4, -12,  -5, -15,  -9,  -1,  -6,  -4, -12,  -5,  -6,  -9, -10, -14}, // D
        { -5, -14,  -1,   7, -14,  -6,  -6, -10,  -2,  -8,  -8,  -5,  -6,  -1, -10,  -5,  -6,  -7, -16,  -7}, // E
        { -8, -13, -16, -14,   9, -11,  -5,  -5, -14,  -3,  -3, -10,  -9,  -8, -10,  -7,  -8,  -6,  -3,   0}, // F
        { -4,  -7,  -6,  -6, -11,   7,  -7, -15,  -7, -11, -10,  -5,  -8,  -8,  -7,  -4,  -8, -10,  -9, -10}, // G
        { -7,  -6,  -4,  -6,  -5,  -7,  10,  -9,  -5,  -7, -12,  -3,  -6,  -2,  -3,  -5,  -5,  -8,  -6,  -1}, // H
        { -7,  -5, -12, -10,  -5, -15,  -9,   7,  -9,  -2,  -2,  -9, -10,  -9,  -8,  -9,  -5,   1,  -6,  -8}, // I
        { -6, -13,  -5,  -2, -14,  -7,  -5,  -9,   7,  -8,  -5,  -3,  -6,  -2,   0,  -5,  -4,  -8,  -9,  -8}, // K
        { -7, -12, -15,  -8,  -3, -11,  -7,  -2,  -8,   6,   0,  -9,  -7,  -6,  -8,  -8,  -7,  -3,  -6,  -6}, // L
        { -5,  -4,  -9,  -8,  -3, -10, -12,  -2,  -5,   0,  10,  -7, -10,  -4,  -6,  -8,  -4,  -3, -13, -11}, // M
        { -6,  -8,  -1,  -5, -10,  -5,  -3,  -9,  -3,  -9,  -7,   8,  -8,  -4,  -5,  -2,  -4,  -9, -10,  -6}, // N
        { -4,  -9,  -6,  -6,  -9,  -8,  -6, -10,  -6,  -7, -10,  -8,   9,  -5,  -7,  -4,  -6,  -7,  -9, -15}, // P
        { -5, -13,  -4,  -1,  -8,  -8,  -2,  -9,  -2,  -6,  -4,  -4,  -5,   9,  -2,  -4,  -5,  -7, -15, -12}, // Q
        { -7,  -7, -12, -10, -10,  -7,  -3,  -8,   0,  -8,  -6,  -5,  -7,  -2,   8,  -6,  -6,  -9,  -8,  -7}, // R
        { -2,  -3,  -5,  -5,  -7,  -4,  -5,  -9,  -5,  -8,  -8,  -2,  -4,  -4,  -6,   7,  -1,  -8,  -8,  -6}, // S
        { -3,  -5,  -6,  -6,  -8,  -8,  -5,  -5,  -4,  -7,  -4,  -4,  -6,  -5,  -6,  -1,   8,  -4, -15,  -8}, // T
        { -3,  -3,  -9,  -7,  -6, -10,  -8,   1,  -8,  -3,  -3,  -9,  -7,  -7,  -9,  -8,  -4,   7, -13,  -8}, // V
        { -9, -15, -10, -16,  -3,  -9,  -6,  -6,  -9,  -6, -13, -10,  -9, -15,  -8,  -8, -15, -13,  12,  -3}, // W
        { -8,  -4, -14,  -7,   0, -10,  -1,  -8,  -8,  -6, -11,  -6, -15, -12,  -7,  -6,  -8,  -8,  -3,  10}, // Y
    };
};

// Result of damage-aware clustering
struct ClusterDamageResult {
    uint32_t cluster_id;
    uint32_t cluster_size;
    float damage_score;          // Aggregated damage evidence
    uint8_t ct_mismatches;       // CT-pattern mismatches at 5' end
    uint8_t ga_mismatches;       // GA-pattern mismatches at 3' end
    float avg_vtml_score;        // Average VTML score of mismatches (lower = more damage-like)
};

// VTML-aware consensus builder
class VTMLConsensusScorer {
public:
    static constexpr size_t K = 5;
    static constexpr size_t TERMINAL_WINDOW = 12;
    static constexpr int VTML_THRESHOLD = -8;  // Include all damage patterns

    struct PositionVote {
        char aa;
        uint16_t count;
        float vtml_weight;  // Sum of VTML scores supporting this AA
    };

    // Score a set of aligned proteins for damage evidence
    // offsets: AA offset for each protein to align to cluster coordinate system
    std::vector<ClusterDamageResult> score_clusters(
        const std::vector<std::string>& proteins,
        const std::vector<uint32_t>& cluster_ids,
        const std::vector<int32_t>& offsets,
        float d_max,
        float lambda)
    {
        // Group proteins by cluster
        std::unordered_map<uint32_t, std::vector<size_t>> clusters;
        for (size_t i = 0; i < proteins.size(); ++i) {
            clusters[cluster_ids[i]].push_back(i);
        }

        std::vector<ClusterDamageResult> results(proteins.size());

        for (const auto& [cid, members] : clusters) {
            if (members.size() < 2) {
                // Singleton - no consensus possible
                for (size_t idx : members) {
                    results[idx] = {cid, 1, 0.0f, 0, 0, 0.0f};
                }
                continue;
            }

            // Build position-wise consensus with proper alignment
            auto consensus = build_consensus_aligned(proteins, members, offsets);

            // Score each member against consensus
            for (size_t idx : members) {
                results[idx] = score_against_consensus_aligned(
                    proteins[idx], offsets[idx], consensus, cid, members.size(), d_max, lambda);
            }
        }

        return results;
    }

    // Overload for backwards compatibility (no offsets)
    std::vector<ClusterDamageResult> score_clusters(
        const std::vector<std::string>& proteins,
        const std::vector<uint32_t>& cluster_ids,
        float d_max,
        float lambda)
    {
        std::vector<int32_t> offsets(proteins.size(), 0);
        return score_clusters(proteins, cluster_ids, offsets, d_max, lambda);
    }

private:
    // Build consensus sequence from cluster members with proper alignment
    std::vector<std::vector<PositionVote>> build_consensus_aligned(
        const std::vector<std::string>& proteins,
        const std::vector<size_t>& members,
        const std::vector<int32_t>& offsets)
    {
        // Find max aligned position
        size_t max_aligned_pos = 0;
        for (size_t idx : members) {
            size_t end_pos = offsets[idx] + proteins[idx].length();
            max_aligned_pos = std::max(max_aligned_pos, end_pos);
        }

        // Count votes at each aligned position
        std::vector<std::unordered_map<char, uint16_t>> position_counts(max_aligned_pos);

        for (size_t idx : members) {
            const auto& prot = proteins[idx];
            int32_t offset = offsets[idx];
            for (size_t pos = 0; pos < prot.length(); ++pos) {
                size_t aligned_pos = offset + pos;
                if (aligned_pos < max_aligned_pos) {
                    position_counts[aligned_pos][prot[pos]]++;
                }
            }
        }

        // Convert to sorted vote lists
        std::vector<std::vector<PositionVote>> consensus(max_aligned_pos);
        for (size_t pos = 0; pos < max_aligned_pos; ++pos) {
            for (const auto& [aa, count] : position_counts[pos]) {
                consensus[pos].push_back({aa, count, 0.0f});
            }
            // Sort by count descending
            std::sort(consensus[pos].begin(), consensus[pos].end(),
                      [](const auto& a, const auto& b) { return a.count > b.count; });
        }

        return consensus;
    }

    // Legacy unaligned version
    std::vector<std::vector<PositionVote>> build_consensus(
        const std::vector<std::string>& proteins,
        const std::vector<size_t>& members)
    {
        std::vector<int32_t> dummy_offsets(proteins.size(), 0);
        return build_consensus_aligned(proteins, members, dummy_offsets);
    }

    // Score single protein against consensus with alignment offset
    ClusterDamageResult score_against_consensus_aligned(
        const std::string& protein,
        int32_t offset,
        const std::vector<std::vector<PositionVote>>& consensus,
        uint32_t cluster_id,
        uint32_t cluster_size,
        float d_max,
        float lambda)
    {
        ClusterDamageResult result = {cluster_id, cluster_size, 0.0f, 0, 0, 0.0f};

        float vtml_sum = 0.0f;
        int mismatch_count = 0;

        for (size_t pos = 0; pos < protein.length(); ++pos) {
            size_t aligned_pos = offset + pos;
            if (aligned_pos >= consensus.size()) continue;

            char my_aa = protein[pos];
            if (consensus[aligned_pos].empty()) continue;

            char consensus_aa = consensus[aligned_pos][0].aa;
            if (my_aa == consensus_aa) continue;

            // Mismatch - analyze it
            int8_t vtml = VTML20::score(consensus_aa, my_aa);
            if (vtml == VTML20::INVALID) continue;

            vtml_sum += vtml;
            mismatch_count++;

            // Check if damage-consistent (use protein-local position for terminal check)
            bool is_5prime = (pos < TERMINAL_WINDOW);
            bool is_3prime = (pos >= protein.length() - TERMINAL_WINDOW);

            if (is_5prime && VTML20::is_damage_sub(consensus_aa, my_aa, true)) {
                result.ct_mismatches++;
                // Weight by position (exponential decay from terminus)
                float pos_weight = std::exp(-lambda * static_cast<float>(pos));
                result.damage_score += VTML20::score_to_damage_prob(vtml, d_max) * pos_weight;
            }

            if (is_3prime && VTML20::is_damage_sub(consensus_aa, my_aa, false)) {
                result.ga_mismatches++;
                float dist_3p = static_cast<float>(protein.length() - 1 - pos);
                float pos_weight = std::exp(-lambda * dist_3p);
                result.damage_score += VTML20::score_to_damage_prob(vtml, d_max) * pos_weight;
            }
        }

        if (mismatch_count > 0) {
            result.avg_vtml_score = vtml_sum / mismatch_count;
        }

        return result;
    }

    // Legacy unaligned version
    ClusterDamageResult score_against_consensus(
        const std::string& protein,
        const std::vector<std::vector<PositionVote>>& consensus,
        uint32_t cluster_id,
        uint32_t cluster_size,
        float d_max,
        float lambda)
    {
        return score_against_consensus_aligned(protein, 0, consensus, cluster_id, cluster_size, d_max, lambda);
    }
};

} // namespace agp
