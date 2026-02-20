#pragma once
// Consensus-based search protein X-masking
//
// Uses overlapping reads to identify likely damage positions and
// X-mask them in search proteins for better MMseqs2 alignment.
//
// IMPORTANT: We X-mask, NOT correct!
// - search_protein: X at damage positions (for MMseqs2 - X matches any AA)
// - protein: Original sequence preserved (for damage annotation later)
//
// This allows:
// 1. MMseqs2 to align through damaged positions (X is permissive)
// 2. dart damage-annotate to see original AAs and detect damage
//
// Strategy:
// 1. Cluster reads by genomic coordinates (same gene region)
// 2. Build consensus within each cluster
// 3. Identify positions where a read differs from consensus in damage-consistent way
// 4. X-mask those positions in search protein (preserve original in protein)

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>

namespace dart {

// Damage-consistent substitution patterns at protein level
struct DamagePattern {
    char consensus_aa;  // What consensus shows
    char observed_aa;   // What this read shows
    bool is_ct;         // CT damage (5' end)
    bool is_ga;         // GA damage (3' end)
};

// Check if substitution is damage-consistent
inline bool is_damage_consistent(char consensus, char observed, bool check_ct, bool check_ga) {
    if (check_ct) {
        // CT damage patterns (C→T in DNA → AA changes)
        if ((consensus == 'R' && (observed == 'W' || observed == 'C' || observed == '*')) ||
            (consensus == 'Q' && observed == '*') ||
            (consensus == 'H' && observed == 'Y') ||
            (consensus == 'P' && (observed == 'S' || observed == 'L')) ||
            (consensus == 'T' && (observed == 'I' || observed == 'M')) ||
            (consensus == 'A' && observed == 'V') ||
            (consensus == 'S' && (observed == 'F' || observed == 'L'))) {
            return true;
        }
    }
    if (check_ga) {
        // GA damage patterns (G→A in DNA → AA changes)
        if ((consensus == 'E' && observed == 'K') ||
            (consensus == 'D' && observed == 'N') ||
            (consensus == 'A' && observed == 'T') ||
            (consensus == 'G' && (observed == 'R' || observed == 'S' || observed == 'E' || observed == 'D')) ||
            (consensus == 'V' && (observed == 'I' || observed == 'M')) ||
            (consensus == 'R' && (observed == 'K' || observed == 'Q' || observed == 'H')) ||
            (consensus == 'S' && observed == 'N') ||
            (consensus == 'C' && observed == 'Y') ||
            (consensus == 'W' && observed == '*')) {
            return true;
        }
    }
    return false;
}

// Result of consensus-based correction
struct CorrectionResult {
    std::string corrected_search;  // Search protein with corrections
    size_t n_corrections;          // Number of positions corrected
    size_t n_ct_corrections;       // CT-pattern corrections
    size_t n_ga_corrections;       // GA-pattern corrections
    float confidence;              // Average confidence of corrections
};

// Position-level correction info
struct PositionCorrection {
    size_t pos;           // Position in protein
    char original;        // Original AA
    char corrected;       // Corrected AA (or 'X' for mask)
    float confidence;     // Confidence score
    bool is_ct;           // CT damage pattern
    bool is_ga;           // GA damage pattern
};

// Cluster of overlapping reads
struct ReadCluster {
    std::vector<size_t> member_indices;
    int32_t start;  // Cluster start coordinate
    int32_t end;    // Cluster end coordinate
};

// Consensus-based corrector
class ConsensusCorrector {
public:
    static constexpr size_t TERMINAL_WINDOW = 15;  // AAs to check at each end
    static constexpr size_t MIN_CLUSTER_SIZE = 3;  // Minimum for consensus
    static constexpr float MIN_CONSENSUS_FREQ = 0.6f;  // Minimum frequency for consensus AA
    static constexpr float MIN_DAMAGE_PROB = 0.1f;  // Minimum d_max to apply corrections

    struct ProteinInfo {
        std::string sequence;
        std::string search_protein;
        int32_t start;
        int32_t end;
        std::string gene_id;
    };

    // Correct search proteins using consensus from overlapping reads
    std::vector<CorrectionResult> correct_search_proteins(
        const std::vector<ProteinInfo>& proteins,
        float d_max,
        float lambda,
        bool use_x_mask = true)  // true = X-mask, false = correct to consensus
    {
        std::vector<CorrectionResult> results(proteins.size());

        // Initialize results with original search proteins
        for (size_t i = 0; i < proteins.size(); ++i) {
            results[i].corrected_search = proteins[i].search_protein;
            results[i].n_corrections = 0;
            results[i].n_ct_corrections = 0;
            results[i].n_ga_corrections = 0;
            results[i].confidence = 0.0f;
        }

        if (d_max < MIN_DAMAGE_PROB) {
            // No significant damage, don't apply corrections
            return results;
        }

        // Group by gene_id
        std::unordered_map<std::string, std::vector<size_t>> gene_groups;
        for (size_t i = 0; i < proteins.size(); ++i) {
            if (!proteins[i].gene_id.empty()) {
                gene_groups[proteins[i].gene_id].push_back(i);
            }
        }

        // Process each gene group
        for (const auto& [gene_id, members] : gene_groups) {
            // Build clusters of overlapping reads
            auto clusters = build_clusters(proteins, members);

            // Process each cluster
            for (const auto& cluster : clusters) {
                if (cluster.member_indices.size() < MIN_CLUSTER_SIZE) {
                    continue;
                }

                // Find corrections for each member
                for (size_t idx : cluster.member_indices) {
                    auto corrections = find_corrections(
                        proteins, idx, cluster, d_max, lambda);

                    // Apply corrections
                    apply_corrections(results[idx], corrections, use_x_mask);
                }
            }
        }

        return results;
    }

private:
    // Build clusters of coordinate-overlapping reads
    std::vector<ReadCluster> build_clusters(
        const std::vector<ProteinInfo>& proteins,
        const std::vector<size_t>& members)
    {
        if (members.empty()) return {};

        // Sort by start position
        std::vector<size_t> sorted = members;
        std::sort(sorted.begin(), sorted.end(),
            [&](size_t a, size_t b) { return proteins[a].start < proteins[b].start; });

        // Merge overlapping
        std::vector<ReadCluster> clusters;
        ReadCluster current;
        current.member_indices.push_back(sorted[0]);
        current.start = proteins[sorted[0]].start;
        current.end = proteins[sorted[0]].end;

        for (size_t i = 1; i < sorted.size(); ++i) {
            size_t idx = sorted[i];
            int32_t s = proteins[idx].start;
            int32_t e = proteins[idx].end;

            if (s <= current.end) {
                // Overlaps - extend cluster
                current.member_indices.push_back(idx);
                current.end = std::max(current.end, e);
            } else {
                // No overlap - start new cluster
                clusters.push_back(std::move(current));
                current = ReadCluster();
                current.member_indices.push_back(idx);
                current.start = s;
                current.end = e;
            }
        }
        clusters.push_back(std::move(current));

        return clusters;
    }

    // Find corrections for a protein based on cluster consensus
    std::vector<PositionCorrection> find_corrections(
        const std::vector<ProteinInfo>& proteins,
        size_t target_idx,
        const ReadCluster& cluster,
        float d_max,
        float lambda)
    {
        std::vector<PositionCorrection> corrections;
        const auto& target = proteins[target_idx];
        size_t prot_len = target.sequence.length();

        // For each terminal position, build consensus and check for damage
        for (size_t pos = 0; pos < prot_len; ++pos) {
            // Skip non-terminal positions
            bool is_5prime = (pos < TERMINAL_WINDOW);
            bool is_3prime = (pos >= prot_len - TERMINAL_WINDOW);
            if (!is_5prime && !is_3prime) continue;

            char target_aa = target.sequence[pos];

            // Compute position in genomic coordinates
            int32_t genomic_pos = target.start + static_cast<int32_t>(pos * 3);

            // Count AAs at this genomic position across cluster
            std::unordered_map<char, int> aa_counts;
            int total = 0;

            for (size_t other_idx : cluster.member_indices) {
                const auto& other = proteins[other_idx];

                // Find corresponding position in other protein
                int32_t offset = genomic_pos - other.start;
                if (offset < 0 || offset >= static_cast<int32_t>(other.sequence.length() * 3)) {
                    continue;
                }

                size_t other_pos = offset / 3;
                if (other_pos < other.sequence.length()) {
                    aa_counts[other.sequence[other_pos]]++;
                    total++;
                }
            }

            if (total < static_cast<int>(MIN_CLUSTER_SIZE)) continue;

            // Find consensus AA
            char consensus_aa = target_aa;
            int max_count = 0;
            for (const auto& [aa, count] : aa_counts) {
                if (count > max_count) {
                    max_count = count;
                    consensus_aa = aa;
                }
            }

            float consensus_freq = static_cast<float>(max_count) / total;
            if (consensus_freq < MIN_CONSENSUS_FREQ) continue;

            // Check if target differs from consensus in damage-consistent way
            if (target_aa != consensus_aa) {
                bool damage_ct = is_5prime && is_damage_consistent(consensus_aa, target_aa, true, false);
                bool damage_ga = is_3prime && is_damage_consistent(consensus_aa, target_aa, false, true);

                if (damage_ct || damage_ga) {
                    // Calculate position-weighted confidence
                    float pos_weight;
                    if (is_5prime) {
                        pos_weight = d_max * std::exp(-lambda * static_cast<float>(pos));
                    } else {
                        float dist_3p = static_cast<float>(prot_len - 1 - pos);
                        pos_weight = d_max * std::exp(-lambda * dist_3p);
                    }

                    float confidence = consensus_freq * pos_weight;

                    corrections.push_back({
                        pos,
                        target_aa,
                        consensus_aa,
                        confidence,
                        damage_ct,
                        damage_ga
                    });
                }
            }
        }

        return corrections;
    }

    // Apply corrections to result
    void apply_corrections(
        CorrectionResult& result,
        const std::vector<PositionCorrection>& corrections,
        bool use_x_mask)
    {
        if (corrections.empty()) return;

        float total_confidence = 0.0f;

        for (const auto& corr : corrections) {
            if (corr.pos < result.corrected_search.length()) {
                if (use_x_mask) {
                    result.corrected_search[corr.pos] = 'X';
                } else {
                    result.corrected_search[corr.pos] = corr.corrected;
                }

                result.n_corrections++;
                if (corr.is_ct) result.n_ct_corrections++;
                if (corr.is_ga) result.n_ga_corrections++;
                total_confidence += corr.confidence;
            }
        }

        if (result.n_corrections > 0) {
            result.confidence = total_confidence / result.n_corrections;
        }
    }
};

} // namespace dart
