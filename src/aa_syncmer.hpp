#pragma once
// Damage-tolerant AA syncmers for protein clustering
//
// Syncmers are k-mers where the minimum s-mer hash is at a specific position.
// This gives more even distribution than minimizers.
//
// For damage tolerance, we use canonical (damage-mapped) hashing:
// - W→R, Y→H, I→T, etc. before computing s-mer hash
// - Same syncmer selected even if damage changes one AA
//
// Workflow:
// 1. Extract syncmers from all proteins (damage-tolerant)
// 2. Cluster proteins by shared syncmers
// 3. Back-translate to nucleotides
// 4. Build NT consensus within clusters
// 5. Identify damage positions where read differs from consensus

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace agp {

// Canonical AA mapping (damaged → undamaged)
// This allows damage variants to hash to the same value
inline char canonical_aa(char c) {
    switch (c) {
        // CT damage products → original
        case 'W': case 'w': return 'R';  // R→W
        case 'Y': case 'y': return 'H';  // H→Y (note: also C→Y for GA)
        case 'I': case 'i': return 'T';  // T→I (note: also V→I for GA)
        case 'M': case 'm': return 'T';  // T→M (note: also V→M for GA)
        case 'F': case 'f': return 'S';  // S→F
        // GA damage products → original
        case 'K': case 'k': return 'E';  // E→K
        case 'N': case 'n': return 'D';  // D→N
        // Uppercase normalization
        case 'a': return 'A'; case 'c': return 'C'; case 'd': return 'D';
        case 'e': return 'E'; case 'g': return 'G'; case 'h': return 'H';
        case 'l': return 'L'; case 'p': return 'P'; case 'q': return 'Q';
        case 'r': return 'R'; case 's': return 'S'; case 't': return 'T';
        case 'v': return 'V';
        default: return (c >= 'A' && c <= 'Z') ? c : 'X';
    }
}

// Hash an s-mer using canonical AA mapping
inline uint32_t hash_smer_canonical(const char* seq, size_t s) {
    uint32_t h = 0;
    for (size_t i = 0; i < s; ++i) {
        char c = canonical_aa(seq[i]);
        if (c == 'X' || c == '*') return UINT32_MAX;  // Invalid
        h = h * 23 + (c - 'A');  // Simple polynomial hash
    }
    return h;
}

// Hash an s-mer directly (no canonicalization)
inline uint32_t hash_smer_direct(const char* seq, size_t s) {
    uint32_t h = 0;
    for (size_t i = 0; i < s; ++i) {
        char c = seq[i];
        if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        if (c < 'A' || c > 'Z' || c == '*') return UINT32_MAX;
        h = h * 23 + (c - 'A');
    }
    return h;
}

// Syncmer: k-mer where minimum s-mer is at position t
struct Syncmer {
    uint64_t hash;       // Hash of the full k-mer
    uint32_t protein_id; // Which protein this came from
    uint16_t position;   // Position in protein

    bool operator<(const Syncmer& o) const { return hash < o.hash; }
    bool operator==(const Syncmer& o) const { return hash == o.hash; }
};

// AA Syncmer extractor with damage tolerance
class AASyncmerExtractor {
public:
    // Tuned for ancient DNA: shorter windows, more overlap
    static constexpr size_t K = 5;   // k-mer size (syncmer window) - shorter for damaged proteins
    static constexpr size_t S = 3;   // s-mer size (for min selection)
    static constexpr size_t T = 0;   // Position where min s-mer must be (0 = open syncmer)

    // Extract syncmers from a protein sequence
    // Uses multi-hash: extracts syncmers for original AND damage variants
    std::vector<Syncmer> extract(const std::string& protein, uint32_t protein_id) const {
        std::vector<Syncmer> syncmers;
        std::unordered_set<uint64_t> seen_hashes;  // Deduplicate

        if (protein.length() < K) return syncmers;

        for (size_t pos = 0; pos <= protein.length() - K; ++pos) {
            const char* kmer = protein.c_str() + pos;

            // Generate original and damage variant k-mers
            std::vector<std::string> variants;
            variants.push_back(std::string(kmer, K));

            // Add single-AA damage variants
            bool is_5prime = (pos < 10);
            bool is_3prime = (pos >= protein.length() - 10);

            std::string base(kmer, K);
            for (size_t i = 0; i < K; ++i) {
                char orig = base[i];
                std::vector<char> subs;

                // CT damage variants (5' end)
                if (is_5prime) {
                    if (orig == 'R') subs = {'W'};
                    else if (orig == 'H') subs = {'Y'};
                    else if (orig == 'T') subs = {'I', 'M'};
                    else if (orig == 'A') subs = {'V'};
                    else if (orig == 'S') subs = {'F'};
                    // Reverse
                    else if (orig == 'W') subs = {'R'};
                    else if (orig == 'Y') subs = {'H'};
                    else if (orig == 'I' || orig == 'M') subs = {'T'};
                    else if (orig == 'V') subs = {'A'};
                    else if (orig == 'F') subs = {'S'};
                }
                // GA damage variants (3' end)
                if (is_3prime) {
                    if (orig == 'E') subs = {'K'};
                    else if (orig == 'D') subs = {'N'};
                    else if (orig == 'V') subs = {'I'};
                    else if (orig == 'G') subs = {'R', 'S'};
                    // Reverse
                    else if (orig == 'K') subs = {'E'};
                    else if (orig == 'N') subs = {'D'};
                    else if (orig == 'I') subs = {'V'};
                    else if (orig == 'R' || orig == 'S') subs.push_back('G');
                }

                for (char sub : subs) {
                    std::string var = base;
                    var[i] = sub;
                    variants.push_back(var);
                }
            }

            // Extract syncmers from each variant
            for (const auto& var : variants) {
                // Find position of minimum s-mer (using canonical hash)
                uint32_t min_hash = UINT32_MAX;
                size_t min_pos = 0;

                for (size_t i = 0; i <= K - S; ++i) {
                    uint32_t h = hash_smer_canonical(var.c_str() + i, S);
                    if (h < min_hash) {
                        min_hash = h;
                        min_pos = i;
                    }
                }

                // Check if this is a syncmer (min at first or last position)
                bool is_syncmer = (min_pos == T) || (min_pos == K - S);

                if (is_syncmer && min_hash != UINT32_MAX) {
                    // Compute full k-mer hash (canonical)
                    uint64_t kmer_hash = 0;
                    bool valid = true;
                    for (size_t i = 0; i < K && valid; ++i) {
                        char c = canonical_aa(var[i]);
                        if (c == 'X' || c == '*') valid = false;
                        else kmer_hash = kmer_hash * 23 + (c - 'A');
                    }

                    if (valid && seen_hashes.insert(kmer_hash).second) {
                        syncmers.push_back({kmer_hash, protein_id, static_cast<uint16_t>(pos)});
                    }
                }
            }
        }

        return syncmers;
    }

    // Extract syncmers from multiple proteins
    std::vector<Syncmer> extract_all(const std::vector<std::string>& proteins) const {
        std::vector<Syncmer> all_syncmers;

        for (uint32_t i = 0; i < proteins.size(); ++i) {
            auto syncmers = extract(proteins[i], i);
            all_syncmers.insert(all_syncmers.end(), syncmers.begin(), syncmers.end());
        }

        return all_syncmers;
    }
};

// Cluster proteins by shared syncmers using Union-Find
class SyncmerClusterer {
public:
    static constexpr size_t MIN_SHARED_SYNCMERS = 2;  // Minimum shared to cluster
    static constexpr size_t MAX_BUCKET_SIZE = 100;    // Skip overcrowded buckets

    struct ClusterResult {
        std::vector<uint32_t> cluster_id;  // Cluster ID for each protein
        uint32_t num_clusters;
    };

    ClusterResult cluster(const std::vector<std::string>& proteins) {
        size_t n = proteins.size();
        if (n == 0) return {{}, 0};

        // Extract all syncmers
        AASyncmerExtractor extractor;
        auto syncmers = extractor.extract_all(proteins);

        // Sort by hash to group matching syncmers
        std::sort(syncmers.begin(), syncmers.end());

        // Initialize Union-Find
        std::vector<uint32_t> parent(n);
        std::vector<uint32_t> rank(n, 0);
        for (size_t i = 0; i < n; ++i) parent[i] = i;

        auto find = [&](uint32_t x) -> uint32_t {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        };

        auto unite = [&](uint32_t x, uint32_t y) {
            uint32_t rx = find(x), ry = find(y);
            if (rx == ry) return;
            if (rank[rx] < rank[ry]) std::swap(rx, ry);
            parent[ry] = rx;
            if (rank[rx] == rank[ry]) rank[rx]++;
        };

        // Process buckets (proteins sharing same syncmer)
        size_t bucket_start = 0;
        while (bucket_start < syncmers.size()) {
            size_t bucket_end = bucket_start + 1;
            while (bucket_end < syncmers.size() &&
                   syncmers[bucket_end].hash == syncmers[bucket_start].hash) {
                ++bucket_end;
            }

            size_t bucket_size = bucket_end - bucket_start;

            if (bucket_size >= 2 && bucket_size <= MAX_BUCKET_SIZE) {
                // Unite all proteins in this bucket
                for (size_t i = bucket_start; i < bucket_end; ++i) {
                    for (size_t j = i + 1; j < bucket_end; ++j) {
                        if (syncmers[i].protein_id != syncmers[j].protein_id) {
                            unite(syncmers[i].protein_id, syncmers[j].protein_id);
                        }
                    }
                }
            }

            bucket_start = bucket_end;
        }

        // Assign cluster IDs
        ClusterResult result;
        result.cluster_id.resize(n);
        std::unordered_map<uint32_t, uint32_t> root_to_cluster;
        uint32_t next_cluster = 0;

        for (uint32_t i = 0; i < n; ++i) {
            uint32_t root = find(i);
            auto it = root_to_cluster.find(root);
            if (it == root_to_cluster.end()) {
                root_to_cluster[root] = next_cluster++;
            }
            result.cluster_id[i] = root_to_cluster[root];
        }

        result.num_clusters = next_cluster;
        return result;
    }
};

// Back-translate protein to possible codons
inline std::vector<std::string> back_translate_aa(char aa) {
    // Simplified: return most common codon(s) for each AA
    static const std::unordered_map<char, std::vector<std::string>> CODON_TABLE = {
        {'A', {"GCT", "GCC", "GCA", "GCG"}},
        {'C', {"TGT", "TGC"}},
        {'D', {"GAT", "GAC"}},
        {'E', {"GAA", "GAG"}},
        {'F', {"TTT", "TTC"}},
        {'G', {"GGT", "GGC", "GGA", "GGG"}},
        {'H', {"CAT", "CAC"}},
        {'I', {"ATT", "ATC", "ATA"}},
        {'K', {"AAA", "AAG"}},
        {'L', {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}},
        {'M', {"ATG"}},
        {'N', {"AAT", "AAC"}},
        {'P', {"CCT", "CCC", "CCA", "CCG"}},
        {'Q', {"CAA", "CAG"}},
        {'R', {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}},
        {'S', {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}},
        {'T', {"ACT", "ACC", "ACA", "ACG"}},
        {'V', {"GTT", "GTC", "GTA", "GTG"}},
        {'W', {"TGG"}},
        {'Y', {"TAT", "TAC"}},
        {'*', {"TAA", "TAG", "TGA"}},
    };

    auto it = CODON_TABLE.find(aa);
    if (it != CODON_TABLE.end()) return it->second;
    return {"NNN"};
}

} // namespace agp
