#pragma once
// Damage-aware protein k-mer indexing
//
// Key insight: instead of exact k-mer matching, we:
// 1. Canonicalize k-mers by mapping damaged→undamaged
// 2. This way "RQHPTV" and "WQYPTV" (damaged) hash to same bucket
// 3. Overlaps survive even with damage mutations
//
// Damage patterns (at protein level):
// CT damage (5'): R→{W,C,*}, Q→*, H→Y, P→{S,L}, T→{I,M}, A→V, S→{F,L}
// GA damage (3'): E→K, D→N, A→T, G→{R,S,E,D}, V→{I,M}, R→{K,Q,H}, S→N, C→Y, W→*

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>

namespace agp {

// Canonical AA: maps damaged variants back to undamaged form
// For example: W→R (because R→W is damage), Y→H (because H→Y is damage)
constexpr char CANONICAL_AA[256] = {
    // Initialize all to identity first, then override damaged forms
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    // @ABCDEFGHIJKLMNOPQRSTUVWXYZ
    0,'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
    'P','Q','R','S','T','U','V','W','X','Y','Z',0,0,0,0,0,
    // Lowercase
    0,'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
    'P','Q','R','S','T','U','V','W','X','Y','Z',0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

// Build damage canonicalization table
// Maps each AA to its "canonical" (undamaged) form
inline std::array<char, 256> build_canonical_table() {
    std::array<char, 256> table;
    // Identity for all
    for (int i = 0; i < 256; ++i) {
        table[i] = static_cast<char>(i);
    }

    // CT damage: damaged → undamaged
    // R→W: W is damaged form of R
    table['W'] = 'R';  table['w'] = 'R';
    // R→C: C could be damaged R (but C is also natural, so skip)
    // R→*: * is damaged form of R or Q
    table['*'] = 'Q';  // Map stop to Q (most common source)

    // H→Y: Y is damaged form of H
    table['Y'] = 'H';  table['y'] = 'H';

    // P→S: S could be damaged P (but S is natural, ambiguous)
    // P→L: L could be damaged P (but L is natural, ambiguous)

    // T→I: I could be damaged T
    table['I'] = 'T';  table['i'] = 'T';
    // T→M: M could be damaged T
    table['M'] = 'T';  table['m'] = 'T';

    // A→V: V could be damaged A (CT) - but also V→I (GA), skip

    // S→F: F could be damaged S
    table['F'] = 'S';  table['f'] = 'S';
    // S→L: L could be damaged S (ambiguous with P→L)

    // GA damage: damaged → undamaged
    // E→K: K could be damaged E
    table['K'] = 'E';  table['k'] = 'E';

    // D→N: N could be damaged D
    table['N'] = 'D';  table['n'] = 'D';

    // A→T: T could be damaged A (but T→I is CT damage, complex)

    // G→R: R could be damaged G (but R→W is CT, complex)
    // G→S: S could be damaged G
    // G→E: E could be damaged G (but E→K is GA)
    // G→D: D could be damaged G (but D→N is GA)

    // V→I: already mapped I→T above (choose one)
    // V→M: already mapped M→T above

    // C→Y: already mapped Y→H above (ambiguous)

    // For simplicity, use a reduced set of unambiguous mappings:
    // W→R, Y→H, I→T, M→T, F→S, K→E, N→D, *→Q

    return table;
}

// Global canonical table
inline const std::array<char, 256>& get_canonical_table() {
    static const auto table = build_canonical_table();
    return table;
}

// Canonicalize a k-mer: map all damaged AAs to undamaged form
inline std::string canonicalize_kmer(const std::string& kmer) {
    const auto& table = get_canonical_table();
    std::string result;
    result.reserve(kmer.size());
    for (char c : kmer) {
        char canonical = table[static_cast<uint8_t>(c)];
        result += (canonical != 0) ? canonical : c;
    }
    return result;
}

// Hash a canonical k-mer
inline uint64_t hash_canonical_kmer(const char* seq, size_t k) {
    const auto& table = get_canonical_table();

    // Simple rolling hash with canonical form
    uint64_t h = 0;
    for (size_t i = 0; i < k; ++i) {
        char c = seq[i];
        char canonical = table[static_cast<uint8_t>(c)];
        if (canonical == 0) canonical = c;

        // Map to 0-20
        uint8_t code;
        if (canonical >= 'A' && canonical <= 'Z') {
            code = canonical - 'A';
        } else {
            return UINT64_MAX;  // Invalid
        }

        h = h * 26 + code;
    }
    return h;
}

// Hash a k-mer directly (no canonicalization)
inline uint64_t hash_kmer_direct(const char* seq, size_t k) {
    uint64_t h = 0;
    for (size_t i = 0; i < k; ++i) {
        char c = seq[i];
        uint8_t code;
        if (c >= 'A' && c <= 'Z') {
            code = c - 'A';
        } else if (c >= 'a' && c <= 'z') {
            code = c - 'a';
        } else if (c == '*') {
            code = 26;  // Stop codon
        } else {
            return UINT64_MAX;
        }
        h = h * 27 + code;
    }
    return h;
}

// VTML20 score thresholds for filtering damage variants
// High-confidence substitutions (VTML20 >= -4) are always included
// Medium-confidence (-8 to -4) are included by default but can be filtered
// Low-confidence (< -8) are excluded
enum class DamageConfidence { HIGH = -4, MEDIUM = -8, LOW = -16 };

// Generate all hashes for a k-mer including damage variants
// This is the MULTI-HASH approach: index original + all single-damage variants
// BIDIRECTIONAL: generates both undamaged→damaged AND damaged→undamaged variants
// vtml_threshold: only include substitutions with VTML20 score >= this value (-6 recommended)
inline std::vector<uint64_t> hash_kmer_with_variants(const char* seq, size_t k, bool is_5prime,
                                                      int vtml_threshold = -8) {
    std::vector<uint64_t> hashes;

    // Original hash
    uint64_t orig = hash_kmer_direct(seq, k);
    if (orig != UINT64_MAX) {
        hashes.push_back(orig);
    }

    std::string kmer(seq, k);

    // Damage substitution with VTML20 scores: {from, to, score}
    // Stop codons (*) use score 0 to always be included (biologically critical)
    struct DamageSub { char from; char to; int vtml_score; };

    // CT damage (C→T at 5' end)
    static const DamageSub CT_DAMAGE_FWD[] = {
        {'R', 'W', -8}, {'R', 'C', -7},   // R→W/C
        {'Q', '*', 0},                     // Q→* (stop codon - always include)
        {'H', 'Y', -1},                    // H→Y (HIGH)
        {'P', 'S', -4}, {'P', 'L', -7},    // P→S (HIGH), P→L
        {'T', 'I', -5}, {'T', 'M', -4},    // T→I/M
        {'A', 'V', -3},                    // A→V (HIGH)
        {'S', 'F', -7}, {'S', 'L', -8},    // S→F/L
    };

    // GA damage (G→A at 3' end)
    static const DamageSub GA_DAMAGE_FWD[] = {
        {'E', 'K', -2},                    // E→K (HIGH)
        {'D', 'N', -1},                    // D→N (HIGH)
        {'A', 'T', -3},                    // A→T (HIGH)
        {'G', 'R', -7}, {'G', 'S', -4}, {'G', 'E', -6}, {'G', 'D', -6},
        {'V', 'I', 1}, {'V', 'M', -3},     // V→I (HIGHEST!), V→M (HIGH)
        {'R', 'K', 0}, {'R', 'Q', -2}, {'R', 'H', -3},  // R→K/Q/H (HIGH)
        {'S', 'N', -2},                    // S→N (HIGH)
        {'C', 'Y', -4},                    // C→Y (HIGH)
        {'W', '*', 0},                     // W→* (stop codon - always include)
    };

    // Reverse maps with VTML scores (same score - symmetric matrix)
    static const DamageSub CT_DAMAGE_REV[] = {
        {'W', 'R', -8}, {'C', 'R', -7}, {'Y', 'H', -1},
        {'S', 'P', -4}, {'L', 'P', -7}, {'L', 'S', -8},
        {'I', 'T', -5}, {'M', 'T', -4},
        {'V', 'A', -3}, {'F', 'S', -7},
        {'*', 'Q', 0}, {'*', 'R', 0},       // * can come from Q or R
    };
    static const DamageSub GA_DAMAGE_REV[] = {
        {'K', 'E', -2}, {'K', 'R', 0},
        {'N', 'D', -1}, {'N', 'S', -2},
        {'T', 'A', -3},
        {'R', 'G', -7}, {'S', 'G', -4}, {'E', 'G', -6}, {'D', 'G', -6},
        {'I', 'V', 1}, {'M', 'V', -3},
        {'Q', 'R', -2}, {'H', 'R', -3},
        {'Y', 'C', -4},
        {'*', 'W', 0},                      // * can come from W
    };

    // Generate variants with VTML score filtering
    auto generate_variants = [&](const DamageSub* damage_map, size_t n) {
        for (size_t pos = 0; pos < k; ++pos) {
            char orig_aa = kmer[pos];
            for (size_t i = 0; i < n; ++i) {
                const auto& sub = damage_map[i];
                if (orig_aa == sub.from && sub.vtml_score >= vtml_threshold) {
                    std::string variant = kmer;
                    variant[pos] = sub.to;
                    uint64_t h = hash_kmer_direct(variant.c_str(), k);
                    if (h != UINT64_MAX) {
                        hashes.push_back(h);
                    }
                }
            }
        }
    };

    if (is_5prime) {
        generate_variants(CT_DAMAGE_FWD, sizeof(CT_DAMAGE_FWD)/sizeof(CT_DAMAGE_FWD[0]));
        generate_variants(CT_DAMAGE_REV, sizeof(CT_DAMAGE_REV)/sizeof(CT_DAMAGE_REV[0]));
    } else {
        generate_variants(GA_DAMAGE_FWD, sizeof(GA_DAMAGE_FWD)/sizeof(GA_DAMAGE_FWD[0]));
        generate_variants(GA_DAMAGE_REV, sizeof(GA_DAMAGE_REV)/sizeof(GA_DAMAGE_REV[0]));
    }

    return hashes;
}

// Generate all damage variants of a k-mer (for validation/debugging)
inline std::vector<std::string> generate_damage_variants(const std::string& kmer) {
    std::vector<std::string> variants;
    variants.push_back(kmer);  // Original

    // Damage substitution map: original → possible damaged forms
    static const std::pair<char, std::vector<char>> DAMAGE_MAP[] = {
        {'R', {'W', 'C'}},      // CT damage
        {'Q', {'*'}},           // CT damage → stop
        {'H', {'Y'}},           // CT damage
        {'P', {'S', 'L'}},      // CT damage
        {'T', {'I', 'M'}},      // CT damage
        {'A', {'V', 'T'}},      // CT (V) or GA (T)
        {'S', {'F', 'L', 'N'}}, // CT (F,L) or GA (N)
        {'E', {'K'}},           // GA damage
        {'D', {'N'}},           // GA damage
        {'G', {'R', 'S', 'E', 'D'}}, // GA damage
        {'V', {'I', 'M'}},      // GA damage
        {'C', {'Y'}},           // GA damage
        {'W', {'*'}},           // GA damage → stop
    };

    // Generate single-mutation variants
    for (size_t pos = 0; pos < kmer.size(); ++pos) {
        char orig = kmer[pos];
        for (const auto& [from, tos] : DAMAGE_MAP) {
            if (orig == from) {
                for (char to : tos) {
                    std::string variant = kmer;
                    variant[pos] = to;
                    variants.push_back(variant);
                }
            }
        }
    }

    return variants;
}

// Damage-aware seed structure
struct DamageAwareSeed {
    uint64_t canonical_hash;  // Hash of canonicalized k-mer
    uint32_t read_id;
    uint16_t pos;             // Position in protein
    uint8_t  is_suffix;       // 0=prefix, 1=suffix
    uint8_t  n_damage_pos;    // Number of positions that could be damaged

    bool operator<(const DamageAwareSeed& o) const {
        return canonical_hash < o.canonical_hash;
    }
};

// Count how many positions in k-mer are "damageable"
inline uint8_t count_damageable(const char* seq, size_t k, bool is_5prime) {
    uint8_t count = 0;

    // CT damage at 5' end affects: R, Q, H, P, T, A, S
    // GA damage at 3' end affects: E, D, A, G, V, R, S, C, W

    if (is_5prime) {
        for (size_t i = 0; i < k; ++i) {
            char c = seq[i];
            if (c == 'R' || c == 'Q' || c == 'H' || c == 'P' ||
                c == 'T' || c == 'A' || c == 'S') {
                ++count;
            }
        }
    } else {
        for (size_t i = 0; i < k; ++i) {
            char c = seq[i];
            if (c == 'E' || c == 'D' || c == 'A' || c == 'G' ||
                c == 'V' || c == 'R' || c == 'S' || c == 'C' || c == 'W') {
                ++count;
            }
        }
    }

    return count;
}

// Damage-aware overlap clusterer
class DamageAwareClusterer {
public:
    static constexpr size_t K = 5;               // Shorter k for better survival
    static constexpr size_t TERMINAL_WINDOW = 10;
    static constexpr size_t MAX_BUCKET = 200;

    struct ClusterResult {
        std::vector<uint32_t> component_id;
        std::vector<int32_t> coord;
        size_t num_components;
    };

    ClusterResult cluster(const std::vector<std::string>& proteins) {
        const size_t n = proteins.size();
        if (n == 0) return {{}, {}, 0};

        // 1. Extract damage-aware seeds
        std::vector<DamageAwareSeed> seeds;
        seeds.reserve(n * 4);

        for (uint32_t i = 0; i < n; ++i) {
            extract_seeds(proteins[i], i, seeds);
        }

        // 2. Sort by canonical hash
        std::sort(seeds.begin(), seeds.end());

        // 3. Initialize Union-Find
        std::vector<uint32_t> parent(n);
        std::vector<int32_t> offset(n, 0);
        for (size_t i = 0; i < n; ++i) parent[i] = i;

        auto find = [&](uint32_t x) -> uint32_t {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];  // Path compression
                x = parent[x];
            }
            return x;
        };

        auto unite = [&](uint32_t x, uint32_t y, int32_t delta) {
            uint32_t rx = find(x), ry = find(y);
            if (rx != ry) {
                parent[rx] = ry;
            }
        };

        // 4. Process buckets
        size_t bucket_start = 0;
        while (bucket_start < seeds.size()) {
            size_t bucket_end = bucket_start + 1;
            while (bucket_end < seeds.size() &&
                   seeds[bucket_end].canonical_hash == seeds[bucket_start].canonical_hash) {
                ++bucket_end;
            }

            size_t bucket_size = bucket_end - bucket_start;
            if (bucket_size > 1 && bucket_size <= MAX_BUCKET) {
                // Unite all in bucket
                for (size_t i = bucket_start; i < bucket_end; ++i) {
                    for (size_t j = i + 1; j < bucket_end; ++j) {
                        if (seeds[i].read_id != seeds[j].read_id) {
                            int32_t delta = seeds[j].pos - seeds[i].pos;
                            unite(seeds[i].read_id, seeds[j].read_id, delta);
                        }
                    }
                }
            }

            bucket_start = bucket_end;
        }

        // 5. Extract results
        ClusterResult result;
        result.component_id.resize(n);
        result.coord.resize(n);

        std::unordered_map<uint32_t, uint32_t> root_to_cid;
        uint32_t next_cid = 0;

        for (uint32_t i = 0; i < n; ++i) {
            uint32_t root = find(i);
            auto it = root_to_cid.find(root);
            if (it == root_to_cid.end()) {
                root_to_cid[root] = next_cid++;
            }
            result.component_id[i] = root_to_cid[root];
            result.coord[i] = 0;  // Simplified - no offset tracking
        }

        result.num_components = next_cid;
        return result;
    }

private:
    void extract_seeds(const std::string& prot, uint32_t read_id,
                       std::vector<DamageAwareSeed>& seeds) {
        if (prot.length() < K) return;

        // Prefix seeds (5' end - CT damage)
        size_t prefix_end = std::min(TERMINAL_WINDOW, prot.length()) - K + 1;
        for (size_t pos = 0; pos < prefix_end; ++pos) {
            uint64_t h = hash_canonical_kmer(prot.c_str() + pos, K);
            if (h != UINT64_MAX) {
                uint8_t n_dmg = count_damageable(prot.c_str() + pos, K, true);
                seeds.push_back({h, read_id, static_cast<uint16_t>(pos), 0, n_dmg});
            }
        }

        // Suffix seeds (3' end - GA damage)
        size_t suffix_start = prot.length() > TERMINAL_WINDOW ?
                              prot.length() - TERMINAL_WINDOW : 0;
        for (size_t pos = std::max(suffix_start, prefix_end);
             pos <= prot.length() - K; ++pos) {
            uint64_t h = hash_canonical_kmer(prot.c_str() + pos, K);
            if (h != UINT64_MAX) {
                uint8_t n_dmg = count_damageable(prot.c_str() + pos, K, false);
                seeds.push_back({h, read_id, static_cast<uint16_t>(pos), 1, n_dmg});
            }
        }
    }
};

} // namespace agp
