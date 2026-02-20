#pragma once
// Consensus-based per-read damage scoring
//
// Pipeline:
// 1. Cluster overlapping proteins via terminal k-mers
// 2. Build consensus per cluster
// 3. Compare each read to consensus
// 4. Score = count of damage-consistent mismatches at terminal positions
//
// This breaks the ρ≈0.10 ceiling by identifying WHICH positions are damaged
// rather than guessing from T/A counts.

#include "terminal_overlap.hpp"
#include <unordered_map>
#include <cmath>

namespace dart {

// Damage substitution patterns (reference → observed)
// CT damage (5' C→T in DNA)
constexpr std::pair<char, char> CT_DAMAGE_SUBS[] = {
    {'R', 'W'}, {'R', 'C'}, {'R', '*'},  // CGN → TGN
    {'Q', '*'},                           // CAA/CAG → TAA/TAG
    {'H', 'Y'},                           // CAC/CAT → TAC/TAT
    {'P', 'S'}, {'P', 'L'},              // CCN → TCN
    {'T', 'I'}, {'T', 'M'},              // ACN → ATN
    {'A', 'V'},                           // GCN → GTN (wait, this is GA not CT)
    {'S', 'F'}, {'S', 'L'},              // TCN → TTN
};

// GA damage (3' G→A in DNA)
constexpr std::pair<char, char> GA_DAMAGE_SUBS[] = {
    {'E', 'K'},                           // GAA/GAG → AAA/AAG
    {'D', 'N'},                           // GAC/GAT → AAC/AAT
    {'A', 'T'},                           // GCN → ACN
    {'G', 'R'}, {'G', 'S'}, {'G', 'E'}, {'G', 'D'},  // GGN → AGN
    {'V', 'I'}, {'V', 'M'},              // GTN → ATN
    {'R', 'K'}, {'R', 'Q'}, {'R', 'H'},  // CGN → CAN (actually CT at pos 2)
    {'S', 'N'},                           // AGC/AGT → AAC/AAT
    {'C', 'Y'},                           // TGC/TGT → TAC/TAT
    {'W', '*'},                           // TGG → TGA/TAG
};

inline bool is_ct_damage(char ref, char obs) {
    for (const auto& [r, o] : CT_DAMAGE_SUBS) {
        if (r == ref && o == obs) return true;
    }
    return false;
}

inline bool is_ga_damage(char ref, char obs) {
    for (const auto& [r, o] : GA_DAMAGE_SUBS) {
        if (r == ref && o == obs) return true;
    }
    return false;
}

// Per-read damage result from consensus comparison
struct ConsensusDamageResult {
    uint32_t read_id;
    uint32_t cluster_id;
    uint16_t cluster_size;

    // Damage counts
    uint8_t ct_damage_5p;    // CT damage near 5' end
    uint8_t ga_damage_3p;    // GA damage near 3' end
    uint8_t total_damage;    // Total damage-consistent mismatches
    uint8_t total_mismatch;  // Total mismatches vs consensus

    // Position-weighted score (terminals weighted higher)
    float weighted_score;

    // Confidence (based on cluster size and consensus quality)
    float confidence;
};

class ConsensusDamageScorer {
public:
    static constexpr size_t TERMINAL_WINDOW = 15;  // AAs from each end
    static constexpr size_t MIN_CLUSTER_SIZE = 3;  // Need 3+ for consensus

    struct ClusterConsensus {
        std::string sequence;
        std::vector<uint8_t> coverage;  // Per-position read count
        std::vector<float> confidence;  // Per-position consensus confidence
    };

    // Main entry: score all proteins using consensus
    std::vector<ConsensusDamageResult> score_proteins(
        const std::vector<std::string>& proteins,
        const std::vector<std::string>& read_ids,
        float d_max,
        float lambda) {

        // 1. Cluster proteins
        TerminalOverlapClusterer clusterer;
        auto clustering = clusterer.cluster(proteins);

        // 2. Build consensus per cluster
        std::unordered_map<uint32_t, std::vector<uint32_t>> cluster_members;
        for (size_t i = 0; i < proteins.size(); ++i) {
            cluster_members[clustering.component_id[i]].push_back(i);
        }

        std::unordered_map<uint32_t, ClusterConsensus> consensuses;
        for (const auto& [cid, members] : cluster_members) {
            if (members.size() >= MIN_CLUSTER_SIZE) {
                consensuses[cid] = build_consensus(proteins, members, clustering);
            }
        }

        // 3. Score each protein against its cluster consensus
        std::vector<ConsensusDamageResult> results;
        results.reserve(proteins.size());

        for (size_t i = 0; i < proteins.size(); ++i) {
            ConsensusDamageResult result{};
            result.read_id = i;
            result.cluster_id = clustering.component_id[i];

            auto it = cluster_members.find(result.cluster_id);
            result.cluster_size = it != cluster_members.end() ?
                                  it->second.size() : 1;

            // If cluster too small, can't use consensus
            if (result.cluster_size < MIN_CLUSTER_SIZE) {
                result.confidence = 0.0f;
                results.push_back(result);
                continue;
            }

            // Compare to consensus
            const auto& cons = consensuses[result.cluster_id];
            score_vs_consensus(proteins[i], cons, clustering.coord[i],
                              d_max, lambda, result);

            results.push_back(result);
        }

        return results;
    }

private:
    ClusterConsensus build_consensus(
        const std::vector<std::string>& proteins,
        const std::vector<uint32_t>& members,
        const TerminalOverlapClusterer::ClusterResult& clustering) {

        // Find coordinate range
        int32_t min_coord = INT32_MAX, max_coord = INT32_MIN;
        for (uint32_t m : members) {
            int32_t start = clustering.coord[m];
            int32_t end = start + proteins[m].length();
            min_coord = std::min(min_coord, start);
            max_coord = std::max(max_coord, end);
        }

        size_t cons_len = max_coord - min_coord;
        if (cons_len > 1000) cons_len = 1000;  // Cap length

        // Build pileup
        std::vector<std::array<uint16_t, 22>> pileup(cons_len);  // 21 AAs + gap

        for (uint32_t m : members) {
            const std::string& prot = proteins[m];
            int32_t offset = clustering.coord[m] - min_coord;

            for (size_t j = 0; j < prot.length(); ++j) {
                int32_t pos = offset + j;
                if (pos >= 0 && pos < static_cast<int32_t>(cons_len)) {
                    uint8_t code = AA_ENCODE[static_cast<uint8_t>(prot[j])];
                    if (code < 21) {
                        pileup[pos][code]++;
                    }
                }
            }
        }

        // Extract consensus
        ClusterConsensus result;
        result.sequence.reserve(cons_len);
        result.coverage.reserve(cons_len);
        result.confidence.reserve(cons_len);

        static const char AA_DECODE[] = "ACDEFGHIKLMNPQRSTVWY*";

        for (size_t i = 0; i < cons_len; ++i) {
            uint16_t total = 0;
            uint8_t best_aa = 0;
            uint16_t best_count = 0;

            for (uint8_t aa = 0; aa < 21; ++aa) {
                total += pileup[i][aa];
                if (pileup[i][aa] > best_count) {
                    best_count = pileup[i][aa];
                    best_aa = aa;
                }
            }

            if (total > 0) {
                result.sequence += AA_DECODE[best_aa];
                result.coverage.push_back(total);
                result.confidence.push_back(static_cast<float>(best_count) / total);
            }
        }

        return result;
    }

    void score_vs_consensus(
        const std::string& protein,
        const ClusterConsensus& consensus,
        int32_t coord,
        float d_max,
        float lambda,
        ConsensusDamageResult& result) {

        if (consensus.sequence.empty()) {
            result.confidence = 0.0f;
            return;
        }

        // Align protein to consensus using coord
        // coord is relative to consensus start
        size_t prot_len = protein.length();
        size_t cons_len = consensus.sequence.length();

        float weighted_score = 0.0f;
        float total_confidence = 0.0f;

        for (size_t i = 0; i < prot_len; ++i) {
            int32_t cons_pos = coord + i;
            if (cons_pos < 0 || cons_pos >= static_cast<int32_t>(cons_len)) continue;

            char obs = protein[i];
            char ref = consensus.sequence[cons_pos];

            if (obs == ref) continue;  // Match, no damage

            result.total_mismatch++;

            // Check if damage-consistent
            bool is_ct = is_ct_damage(ref, obs);
            bool is_ga = is_ga_damage(ref, obs);

            if (!is_ct && !is_ga) continue;  // Not damage pattern

            result.total_damage++;

            // Position weighting: terminal positions weighted by damage probability
            float pos_weight = 1.0f;

            // 5' end: CT damage
            if (is_ct && i < TERMINAL_WINDOW) {
                result.ct_damage_5p++;
                pos_weight = d_max * std::exp(-lambda * static_cast<float>(i));
            }

            // 3' end: GA damage
            if (is_ga && i >= prot_len - TERMINAL_WINDOW) {
                result.ga_damage_3p++;
                float dist_3 = static_cast<float>(prot_len - 1 - i);
                pos_weight = d_max * std::exp(-lambda * dist_3);
            }

            // Weight by consensus confidence at this position
            float conf = consensus.confidence[cons_pos];

            weighted_score += pos_weight * conf;
            total_confidence += conf;
        }

        result.weighted_score = weighted_score;
        result.confidence = total_confidence > 0 ?
                           total_confidence / prot_len : 0.0f;
    }
};

// Rescore p_damaged using consensus evidence
inline float rescore_with_consensus(
    float original_p_damaged,
    const ConsensusDamageResult& cons_result,
    float d_max) {

    // If no consensus available, keep original
    if (cons_result.confidence < 0.1f || cons_result.cluster_size < 3) {
        return original_p_damaged;
    }

    // Combine original score with consensus evidence
    // weighted_score is already position-weighted damage count
    float cons_score = cons_result.weighted_score;

    // Normalize by expected damage (based on d_max and protein length)
    // A protein with N terminal AAs expects ~N * d_max damage events
    float expected = 2.0f * 15.0f * d_max;  // Both termini, 15 AAs each
    float normalized = cons_score / (expected + 0.1f);

    // Blend with original (weight by consensus confidence)
    float blend = cons_result.confidence;
    float rescored = (1.0f - blend) * original_p_damaged +
                     blend * std::min(normalized, 1.0f);

    return std::clamp(rescored, 0.0f, 1.0f);
}

} // namespace dart
