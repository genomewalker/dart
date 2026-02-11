#pragma once
// Terminal-seeded protein overlap clustering for damage detection
//
// Key ideas from Codex brainstorm:
// 1. Index only terminal k-mers (first/last ~12 AA)
// 2. Radix sort by seed hash for O(n) grouping
// 3. Diagonal consistency: require 2+ seeds on same diagonal
// 4. Union-Find with potentials: track coordinates without graph
// 5. Streaming pileup: consensus without storing edges
//
// Target: 10M proteins in <60s

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <immintrin.h>

namespace agp {

// 5-bit AA encoding: 20 AAs + X/stop → fits in 5 bits
// Pack 12 AAs into 64 bits (60 bits used)
constexpr uint8_t AA_ENCODE[256] = {
    // Default to 31 (invalid)
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31, 0,31, 1, 2, 3, 4, 5, 6, 7,31, 8, 9,10,11,31,  // @ABCDEFGHIJKLMNO
    12,13,14,15,16,31,17,18,20,19,31,31,31,31,31,31,  // PQRSTUVWXYZ
    31, 0,31, 1, 2, 3, 4, 5, 6, 7,31, 8, 9,10,11,31,  // lowercase
    12,13,14,15,16,31,17,18,20,19,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
    31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
};

// Seed: terminal k-mer with position info
struct Seed {
    uint64_t hash;       // k-mer hash (k=8 → 40 bits needed)
    uint32_t read_id;    // which protein
    uint16_t pos;        // position in protein (0 = prefix, len-k = suffix)
    uint8_t  is_suffix;  // 0=prefix seed, 1=suffix seed
    uint8_t  reserved;

    bool operator<(const Seed& o) const { return hash < o.hash; }
};

// Union-Find with potentials (offset tracking)
// Each element stores offset relative to its parent
class DSUWithOffset {
public:
    explicit DSUWithOffset(size_t n) : parent(n), rank(n, 0), offset(n, 0) {
        for (size_t i = 0; i < n; ++i) parent[i] = i;
    }

    // Find root and compute cumulative offset from element to root
    std::pair<uint32_t, int32_t> find(uint32_t x) {
        if (parent[x] == x) return {x, 0};
        auto [root, parent_offset] = find(parent[x]);
        // Path compression with offset update
        offset[x] += parent_offset;
        parent[x] = root;
        return {root, offset[x]};
    }

    // Unite two elements with known offset difference: coord[y] - coord[x] = delta
    bool unite(uint32_t x, uint32_t y, int32_t delta) {
        auto [rx, ox] = find(x);
        auto [ry, oy] = find(y);
        if (rx == ry) {
            // Check consistency: existing offset should match
            return (oy - ox) == delta;
        }
        // Union by rank
        if (rank[rx] < rank[ry]) {
            parent[rx] = ry;
            offset[rx] = oy - ox - delta;  // so that find(x) gives correct offset
        } else if (rank[rx] > rank[ry]) {
            parent[ry] = rx;
            offset[ry] = ox - oy + delta;
        } else {
            parent[ry] = rx;
            offset[ry] = ox - oy + delta;
            ++rank[rx];
        }
        return true;
    }

    // Get global coordinate for element (offset from component's origin)
    int32_t coord(uint32_t x) {
        auto [root, off] = find(x);
        return off;
    }

private:
    std::vector<uint32_t> parent;
    std::vector<uint8_t>  rank;
    std::vector<int32_t>  offset;
};

// Hash a k-mer (k=8 typical, uses 40 bits)
inline uint64_t hash_kmer(const char* seq, size_t k) {
    uint64_t h = 0;
    for (size_t i = 0; i < k; ++i) {
        uint8_t code = AA_ENCODE[static_cast<uint8_t>(seq[i])];
        if (code >= 21) return UINT64_MAX;  // Invalid AA
        h = h * 21 + code;
    }
    return h;
}

// Terminal overlap clusterer
class TerminalOverlapClusterer {
public:
    static constexpr size_t K = 6;              // k-mer size (balance sensitivity/specificity)
    static constexpr size_t TERMINAL_WINDOW = 12; // index first/last N AAs
    static constexpr size_t MIN_SEEDS = 2;      // require 2+ seeds on diagonal
    static constexpr size_t MAX_BUCKET_SIZE = 50;   // skip high-freq k-mers
    static constexpr int MAX_MISMATCH = 3;      // for SIMD verification

    struct ClusterResult {
        std::vector<uint32_t> component_id;  // component for each protein
        std::vector<int32_t>  coord;         // coordinate within component
        size_t num_components;
    };

    // Main entry point
    ClusterResult cluster(const std::vector<std::string>& proteins) {
        const size_t n = proteins.size();
        if (n == 0) return {{}, {}, 0};

        // 1. Extract terminal seeds
        std::vector<Seed> seeds;
        seeds.reserve(n * 4);  // ~2 prefix + 2 suffix seeds per protein

        for (uint32_t i = 0; i < n; ++i) {
            extract_seeds(proteins[i], i, seeds);
        }

        // 2. Radix sort by hash (O(n) for 64-bit keys)
        std::sort(seeds.begin(), seeds.end());  // TODO: radix sort

        // 3. Initialize DSU
        DSUWithOffset dsu(n);

        // 4. Process each hash bucket
        size_t bucket_start = 0;
        while (bucket_start < seeds.size()) {
            // Find bucket end
            size_t bucket_end = bucket_start + 1;
            while (bucket_end < seeds.size() &&
                   seeds[bucket_end].hash == seeds[bucket_start].hash) {
                ++bucket_end;
            }

            size_t bucket_size = bucket_end - bucket_start;

            // Skip high-frequency k-mers (likely repeats/low-complexity)
            if (bucket_size > 1 && bucket_size <= MAX_BUCKET_SIZE) {
                process_bucket(seeds, bucket_start, bucket_end, proteins, dsu);
            }

            bucket_start = bucket_end;
        }

        // 5. Extract results
        ClusterResult result;
        result.component_id.resize(n);
        result.coord.resize(n);

        std::unordered_map<uint32_t, uint32_t> root_to_component;
        uint32_t next_component = 0;

        for (uint32_t i = 0; i < n; ++i) {
            auto [root, off] = dsu.find(i);
            auto it = root_to_component.find(root);
            if (it == root_to_component.end()) {
                root_to_component[root] = next_component++;
            }
            result.component_id[i] = root_to_component[root];
            result.coord[i] = off;
        }

        result.num_components = next_component;
        return result;
    }

private:
    void extract_seeds(const std::string& prot, uint32_t read_id,
                       std::vector<Seed>& seeds) {
        if (prot.length() < K) return;

        // Prefix seeds (first TERMINAL_WINDOW - K + 1 positions)
        size_t prefix_end = std::min(TERMINAL_WINDOW, prot.length()) - K + 1;
        for (size_t pos = 0; pos < prefix_end; ++pos) {
            uint64_t h = hash_kmer(prot.c_str() + pos, K);
            if (h != UINT64_MAX) {
                seeds.push_back({h, read_id, static_cast<uint16_t>(pos), 0, 0});
            }
        }

        // Suffix seeds (last TERMINAL_WINDOW positions)
        size_t suffix_start = prot.length() > TERMINAL_WINDOW ?
                              prot.length() - TERMINAL_WINDOW : 0;
        for (size_t pos = std::max(suffix_start, prefix_end);
             pos <= prot.length() - K; ++pos) {
            uint64_t h = hash_kmer(prot.c_str() + pos, K);
            if (h != UINT64_MAX) {
                seeds.push_back({h, read_id, static_cast<uint16_t>(pos), 1, 0});
            }
        }
    }

    void process_bucket(const std::vector<Seed>& seeds,
                        size_t start, size_t end,
                        const std::vector<std::string>& proteins,
                        DSUWithOffset& dsu) {
        // For each pair in bucket, check diagonal consistency
        // Diagonal = pos_j - pos_i (implied overlap offset)

        // Build diagonal histogram for each read pair
        // (simplified: just check pairs directly for small buckets)

        for (size_t i = start; i < end; ++i) {
            for (size_t j = i + 1; j < end; ++j) {
                const Seed& si = seeds[i];
                const Seed& sj = seeds[j];

                if (si.read_id == sj.read_id) continue;

                // Compute implied offset: coord[sj] - coord[si]
                int32_t delta = static_cast<int32_t>(sj.pos) -
                               static_cast<int32_t>(si.pos);

                // Quick SIMD verification of overlap region
                if (verify_overlap(proteins[si.read_id], proteins[sj.read_id],
                                  delta)) {
                    dsu.unite(si.read_id, sj.read_id, delta);
                }
            }
        }
    }

    // SIMD-accelerated overlap verification
    bool verify_overlap(const std::string& a, const std::string& b,
                        int32_t delta) {
        // Overlap region: positions where both proteins have data
        int32_t a_start = std::max(0, -delta);
        int32_t b_start = std::max(0, delta);
        int32_t a_end = std::min(static_cast<int32_t>(a.length()),
                                 static_cast<int32_t>(b.length()) - delta);

        if (a_end <= a_start) return false;

        size_t overlap_len = a_end - a_start;
        if (overlap_len < K) return false;  // Overlap too short

        // Count mismatches (simple loop, TODO: AVX-512)
        int mismatches = 0;
        for (size_t i = 0; i < overlap_len && mismatches <= MAX_MISMATCH; ++i) {
            if (a[a_start + i] != b[b_start + i]) ++mismatches;
        }

        return mismatches <= MAX_MISMATCH;
    }
};

// Pileup for consensus building
class ConsensusPileup {
public:
    static constexpr size_t MAX_COORD = 500;  // Max contig length

    struct PositionCounts {
        std::array<uint16_t, 21> aa_counts{};  // 20 AAs + stop
        uint16_t total = 0;

        void add(char aa) {
            uint8_t code = AA_ENCODE[static_cast<uint8_t>(aa)];
            if (code < 21) {
                ++aa_counts[code];
                ++total;
            }
        }

        char consensus() const {
            static const char AA_DECODE[] = "ACDEFGHIKLMNPQRSTVWY*";
            uint8_t best = 0;
            for (uint8_t i = 1; i < 21; ++i) {
                if (aa_counts[i] > aa_counts[best]) best = i;
            }
            return aa_counts[best] > 0 ? AA_DECODE[best] : 'X';
        }
    };

    void add_protein(const std::string& prot, int32_t coord) {
        // Shift coord to positive range
        int32_t shifted = coord + static_cast<int32_t>(MAX_COORD / 2);
        if (shifted < 0) return;

        for (size_t i = 0; i < prot.length(); ++i) {
            size_t pos = shifted + i;
            if (pos < MAX_COORD) {
                if (pos >= positions.size()) {
                    positions.resize(pos + 1);
                }
                positions[pos].add(prot[i]);
            }
        }
    }

    std::string get_consensus() const {
        std::string result;
        result.reserve(positions.size());
        for (const auto& p : positions) {
            if (p.total > 0) result += p.consensus();
        }
        return result;
    }

    // Damage substitution table: reference → damaged
    // CT damage (5' C→T): R→W, R→C, R→*, Q→*, H→Y, P→S, P→L, T→I, T→M, A→V, S→F, S→L
    // GA damage (3' G→A): E→K, D→N, A→T, G→R, G→S, V→I, V→M, R→K, R→Q, R→H, G→E, G→D, S→N, C→Y, W→*
    static constexpr bool is_damage_sub(char ref, char obs) {
        // CT damage (5' end): C→T in DNA causes these AA changes
        if (ref == 'R' && (obs == 'W' || obs == 'C' || obs == '*')) return true;
        if (ref == 'Q' && obs == '*') return true;
        if (ref == 'H' && obs == 'Y') return true;
        if (ref == 'P' && (obs == 'S' || obs == 'L')) return true;
        if (ref == 'T' && (obs == 'I' || obs == 'M')) return true;
        if (ref == 'A' && obs == 'V') return true;
        if (ref == 'S' && (obs == 'F' || obs == 'L')) return true;

        // GA damage (3' end): G→A in DNA causes these AA changes
        if (ref == 'E' && obs == 'K') return true;
        if (ref == 'D' && obs == 'N') return true;
        if (ref == 'A' && obs == 'T') return true;
        if (ref == 'G' && (obs == 'R' || obs == 'S' || obs == 'E' || obs == 'D')) return true;
        if (ref == 'V' && (obs == 'I' || obs == 'M')) return true;
        if (ref == 'R' && (obs == 'K' || obs == 'Q' || obs == 'H')) return true;
        if (ref == 'S' && obs == 'N') return true;
        if (ref == 'C' && obs == 'Y') return true;
        if (ref == 'W' && obs == '*') return true;

        return false;
    }

    // Get damage candidates: positions where minority variants match damage patterns
    struct DamageCandidate {
        size_t pos;
        char consensus_aa;
        char variant_aa;
        uint16_t consensus_count;
        uint16_t variant_count;
        bool is_ct;  // true=CT damage (5'), false=GA damage (3')
    };

    std::vector<DamageCandidate> find_damage_candidates() const {
        std::vector<DamageCandidate> candidates;
        static const char AA_DECODE[] = "ACDEFGHIKLMNPQRSTVWY*";

        for (size_t i = 0; i < positions.size(); ++i) {
            const auto& p = positions[i];
            if (p.total < 3) continue;  // Need coverage

            char cons = p.consensus();
            uint8_t cons_code = AA_ENCODE[static_cast<uint8_t>(cons)];

            // Look for minority variants that match damage patterns
            for (uint8_t v = 0; v < 21; ++v) {
                if (v == cons_code) continue;
                if (p.aa_counts[v] == 0) continue;

                char var_aa = AA_DECODE[v];

                // Only keep if it's a damage-consistent substitution
                if (!is_damage_sub(cons, var_aa)) continue;

                // Minority check: consensus should be at least 2x variant
                if (p.aa_counts[cons_code] >= 2 * p.aa_counts[v]) {
                    // Determine if CT or GA based on substitution type
                    bool is_ct = (cons == 'R' && (var_aa == 'W' || var_aa == 'C' || var_aa == '*')) ||
                                 (cons == 'Q' && var_aa == '*') ||
                                 (cons == 'H' && var_aa == 'Y') ||
                                 (cons == 'P') || (cons == 'T') ||
                                 (cons == 'A' && var_aa == 'V') ||
                                 (cons == 'S' && (var_aa == 'F' || var_aa == 'L'));

                    candidates.push_back({
                        i, cons, var_aa,
                        p.aa_counts[cons_code], p.aa_counts[v],
                        is_ct
                    });
                }
            }
        }

        return candidates;
    }

private:
    std::vector<PositionCounts> positions;
};

} // namespace agp
