#pragma once

/**
 * @file damage_aa_table.hpp
 * @brief Complete mapping of damage-induced amino acid changes
 *
 * C→T damage (5' on forward, 3' on reverse):
 *   - Deamination of cytosine to uracil (read as thymine)
 *   - Most impactful at codon position 1 (creates stops)
 *
 * G→A damage (3' on forward, 5' on reverse):
 *   - Deamination of guanine on complementary strand
 *   - TGG→TGA (W→*) is the only G→A stop-creating change
 */

#include <array>
#include <cstdint>

namespace dart {

// Amino acid encoding for compact storage
// Standard 1-letter codes + * for stop
constexpr int aa_to_idx(char aa) {
    switch (aa) {
        case 'A': return 0;  case 'C': return 1;  case 'D': return 2;
        case 'E': return 3;  case 'F': return 4;  case 'G': return 5;
        case 'H': return 6;  case 'I': return 7;  case 'K': return 8;
        case 'L': return 9;  case 'M': return 10; case 'N': return 11;
        case 'P': return 12; case 'Q': return 13; case 'R': return 14;
        case 'S': return 15; case 'T': return 16; case 'V': return 17;
        case 'W': return 18; case 'Y': return 19; case '*': return 20;
        default: return -1;
    }
}

constexpr char idx_to_aa(int idx) {
    constexpr char aas[] = "ACDEFGHIKLMNPQRSTVWY*";
    return (idx >= 0 && idx <= 20) ? aas[idx] : 'X';
}

/**
 * Damage-induced amino acid substitution
 */
struct DamageSubstitution {
    char original_codon[4];    // Original codon (before damage)
    char damaged_codon[4];     // Damaged codon (after C→T or G→A)
    char original_aa;          // Original amino acid
    char damaged_aa;           // Damaged amino acid
    int codon_position;        // Position of change (1, 2, or 3)
    bool creates_stop;         // True if damage creates premature stop
    bool is_synonymous;        // True if AA doesn't change
    float diagnostic_weight;   // Higher = more diagnostic of damage
};

// ============================================================================
// C→T DAMAGE SUBSTITUTIONS (5' end on forward strand)
// ============================================================================

// Position 1 C→T changes (most impactful - includes stop codons)
constexpr DamageSubstitution CT_POS1_SUBS[] = {
    // STOP-CREATING (highest diagnostic value)
    {"CAA", "TAA", 'Q', '*', 1, true,  false, 3.0f},  // Gln → Stop
    {"CAG", "TAG", 'Q', '*', 1, true,  false, 3.0f},  // Gln → Stop
    {"CGA", "TGA", 'R', '*', 1, true,  false, 3.0f},  // Arg → Stop

    // Non-synonymous changes
    {"CAT", "TAT", 'H', 'Y', 1, false, false, 1.5f},  // His → Tyr
    {"CAC", "TAC", 'H', 'Y', 1, false, false, 1.5f},  // His → Tyr
    {"CGT", "TGT", 'R', 'C', 1, false, false, 1.5f},  // Arg → Cys
    {"CGC", "TGC", 'R', 'C', 1, false, false, 1.5f},  // Arg → Cys
    {"CGG", "TGG", 'R', 'W', 1, false, false, 1.5f},  // Arg → Trp
    {"CCA", "TCA", 'P', 'S', 1, false, false, 1.0f},  // Pro → Ser
    {"CCT", "TCT", 'P', 'S', 1, false, false, 1.0f},  // Pro → Ser
    {"CCC", "TCC", 'P', 'S', 1, false, false, 1.0f},  // Pro → Ser
    {"CCG", "TCG", 'P', 'S', 1, false, false, 1.0f},  // Pro → Ser
    {"CTT", "TTT", 'L', 'F', 1, false, false, 1.0f},  // Leu → Phe
    {"CTC", "TTC", 'L', 'F', 1, false, false, 1.0f},  // Leu → Phe

    // Synonymous changes (lower diagnostic value)
    {"CTA", "TTA", 'L', 'L', 1, false, true,  0.3f},  // Leu → Leu
    {"CTG", "TTG", 'L', 'L', 1, false, true,  0.3f},  // Leu → Leu
};
constexpr size_t CT_POS1_COUNT = sizeof(CT_POS1_SUBS) / sizeof(CT_POS1_SUBS[0]);

// Position 2 C→T changes
constexpr DamageSubstitution CT_POS2_SUBS[] = {
    // T→I changes (Thr → Ile)
    {"ACT", "ATT", 'T', 'I', 2, false, false, 1.0f},
    {"ACC", "ATC", 'T', 'I', 2, false, false, 1.0f},
    {"ACA", "ATA", 'T', 'I', 2, false, false, 1.0f},
    {"ACG", "ATG", 'T', 'M', 2, false, false, 1.2f},  // Thr → Met (start codon context)

    // A→V changes (Ala → Val)
    {"GCT", "GTT", 'A', 'V', 2, false, false, 1.0f},
    {"GCC", "GTC", 'A', 'V', 2, false, false, 1.0f},
    {"GCA", "GTA", 'A', 'V', 2, false, false, 1.0f},
    {"GCG", "GTG", 'A', 'V', 2, false, false, 1.0f},

    // S→F changes (Ser → Phe)
    {"TCT", "TTT", 'S', 'F', 2, false, false, 1.0f},
    {"TCC", "TTC", 'S', 'F', 2, false, false, 1.0f},

    // S→L changes (Ser → Leu)
    {"TCA", "TTA", 'S', 'L', 2, false, false, 1.0f},
    {"TCG", "TTG", 'S', 'L', 2, false, false, 1.0f},

    // P→L changes (Pro → Leu)
    {"CCT", "CTT", 'P', 'L', 2, false, false, 1.0f},
    {"CCC", "CTC", 'P', 'L', 2, false, false, 1.0f},
    {"CCA", "CTA", 'P', 'L', 2, false, false, 1.0f},
    {"CCG", "CTG", 'P', 'L', 2, false, false, 1.0f},
};
constexpr size_t CT_POS2_COUNT = sizeof(CT_POS2_SUBS) / sizeof(CT_POS2_SUBS[0]);

// Position 3 C→T changes (mostly synonymous - wobble position)
constexpr DamageSubstitution CT_POS3_SUBS[] = {
    // Non-synonymous at position 3 (rare but exists)
    {"TGC", "TGT", 'C', 'C', 3, false, true,  0.2f},  // Cys → Cys (syn)
    {"AGC", "AGT", 'S', 'S', 3, false, true,  0.2f},  // Ser → Ser (syn)
    {"AAC", "AAT", 'N', 'N', 3, false, true,  0.2f},  // Asn → Asn (syn)
    {"GAC", "GAT", 'D', 'D', 3, false, true,  0.2f},  // Asp → Asp (syn)
    {"TTC", "TTT", 'F', 'F', 3, false, true,  0.2f},  // Phe → Phe (syn)
    {"TAC", "TAT", 'Y', 'Y', 3, false, true,  0.2f},  // Tyr → Tyr (syn)
    {"CAC", "CAT", 'H', 'H', 3, false, true,  0.2f},  // His → His (syn)
    {"ATC", "ATT", 'I', 'I', 3, false, true,  0.2f},  // Ile → Ile (syn)
    // Most position 3 changes are synonymous due to wobble
};
constexpr size_t CT_POS3_COUNT = sizeof(CT_POS3_SUBS) / sizeof(CT_POS3_SUBS[0]);

// ============================================================================
// G→A DAMAGE SUBSTITUTIONS (3' end on forward strand)
// ============================================================================

// Position 1 G→A changes
constexpr DamageSubstitution GA_POS1_SUBS[] = {
    // E→K changes (Glu → Lys)
    {"GAA", "AAA", 'E', 'K', 1, false, false, 1.5f},
    {"GAG", "AAG", 'E', 'K', 1, false, false, 1.5f},

    // D→N changes (Asp → Asn)
    {"GAT", "AAT", 'D', 'N', 1, false, false, 1.5f},
    {"GAC", "AAC", 'D', 'N', 1, false, false, 1.5f},

    // G→S changes (Gly → Ser)
    {"GGT", "AGT", 'G', 'S', 1, false, false, 1.2f},
    {"GGC", "AGC", 'G', 'S', 1, false, false, 1.2f},

    // G→R changes (Gly → Arg)
    {"GGA", "AGA", 'G', 'R', 1, false, false, 1.2f},
    {"GGG", "AGG", 'G', 'R', 1, false, false, 1.2f},

    // V→I changes (Val → Ile)
    {"GTT", "ATT", 'V', 'I', 1, false, false, 1.0f},
    {"GTC", "ATC", 'V', 'I', 1, false, false, 1.0f},
    {"GTA", "ATA", 'V', 'I', 1, false, false, 1.0f},
    {"GTG", "ATG", 'V', 'M', 1, false, false, 1.2f},  // Val → Met

    // A→T changes (Ala → Thr)
    {"GCT", "ACT", 'A', 'T', 1, false, false, 1.0f},
    {"GCC", "ACC", 'A', 'T', 1, false, false, 1.0f},
    {"GCA", "ACA", 'A', 'T', 1, false, false, 1.0f},
    {"GCG", "ACG", 'A', 'T', 1, false, false, 1.0f},
};
constexpr size_t GA_POS1_COUNT = sizeof(GA_POS1_SUBS) / sizeof(GA_POS1_SUBS[0]);

// Position 2 G→A changes
constexpr DamageSubstitution GA_POS2_SUBS[] = {
    // R→K changes (Arg → Lys)
    {"AGA", "AAA", 'R', 'K', 2, false, false, 1.0f},
    {"AGG", "AAG", 'R', 'K', 2, false, false, 1.0f},

    // S→N changes (Ser → Asn)
    {"AGT", "AAT", 'S', 'N', 2, false, false, 1.0f},
    {"AGC", "AAC", 'S', 'N', 2, false, false, 1.0f},

    // R→Q changes (Arg → Gln)
    {"CGT", "CAT", 'R', 'H', 2, false, false, 1.0f},  // Actually R→H
    {"CGC", "CAC", 'R', 'H', 2, false, false, 1.0f},
    {"CGA", "CAA", 'R', 'Q', 2, false, false, 1.0f},
    {"CGG", "CAG", 'R', 'Q', 2, false, false, 1.0f},

    // G→D changes (Gly → Asp)
    {"GGT", "GAT", 'G', 'D', 2, false, false, 1.2f},
    {"GGC", "GAC", 'G', 'D', 2, false, false, 1.2f},

    // G→E changes (Gly → Glu)
    {"GGA", "GAA", 'G', 'E', 2, false, false, 1.2f},
    {"GGG", "GAG", 'G', 'E', 2, false, false, 1.2f},
};
constexpr size_t GA_POS2_COUNT = sizeof(GA_POS2_SUBS) / sizeof(GA_POS2_SUBS[0]);

// Position 3 G→A changes (includes the important W→* stop)
constexpr DamageSubstitution GA_POS3_SUBS[] = {
    // STOP-CREATING (high diagnostic value)
    {"TGG", "TGA", 'W', '*', 3, true, false, 3.0f},  // Trp → Stop (ONLY G→A stop!)

    // Most position 3 G→A are synonymous
    {"CTG", "CTA", 'L', 'L', 3, false, true, 0.2f},
    {"TTG", "TTA", 'L', 'L', 3, false, true, 0.2f},
    {"ATG", "ATA", 'M', 'I', 3, false, false, 0.8f},  // Met → Ile (start codon loss!)
    {"GTG", "GTA", 'V', 'V', 3, false, true, 0.2f},
};
constexpr size_t GA_POS3_COUNT = sizeof(GA_POS3_SUBS) / sizeof(GA_POS3_SUBS[0]);

// ============================================================================
// LOOKUP FUNCTIONS
// ============================================================================

/**
 * Check if an amino acid change is consistent with C→T damage
 * Returns diagnostic weight (0 = not damage-consistent, >0 = damage-consistent)
 */
inline float is_ct_damage_change(char original_aa, char damaged_aa) {
    // Check position 1 changes
    for (size_t i = 0; i < CT_POS1_COUNT; ++i) {
        if (CT_POS1_SUBS[i].original_aa == original_aa &&
            CT_POS1_SUBS[i].damaged_aa == damaged_aa) {
            return CT_POS1_SUBS[i].diagnostic_weight;
        }
    }
    // Check position 2 changes
    for (size_t i = 0; i < CT_POS2_COUNT; ++i) {
        if (CT_POS2_SUBS[i].original_aa == original_aa &&
            CT_POS2_SUBS[i].damaged_aa == damaged_aa) {
            return CT_POS2_SUBS[i].diagnostic_weight;
        }
    }
    return 0.0f;
}

/**
 * Check if an amino acid change is consistent with G→A damage
 */
inline float is_ga_damage_change(char original_aa, char damaged_aa) {
    // Check position 1 changes
    for (size_t i = 0; i < GA_POS1_COUNT; ++i) {
        if (GA_POS1_SUBS[i].original_aa == original_aa &&
            GA_POS1_SUBS[i].damaged_aa == damaged_aa) {
            return GA_POS1_SUBS[i].diagnostic_weight;
        }
    }
    // Check position 2 changes
    for (size_t i = 0; i < GA_POS2_COUNT; ++i) {
        if (GA_POS2_SUBS[i].original_aa == original_aa &&
            GA_POS2_SUBS[i].damaged_aa == damaged_aa) {
            return GA_POS2_SUBS[i].diagnostic_weight;
        }
    }
    // Check position 3 changes
    for (size_t i = 0; i < GA_POS3_COUNT; ++i) {
        if (GA_POS3_SUBS[i].original_aa == original_aa &&
            GA_POS3_SUBS[i].damaged_aa == damaged_aa) {
            return GA_POS3_SUBS[i].diagnostic_weight;
        }
    }
    return 0.0f;
}

/**
 * Get the original amino acid if this could be damage-induced
 * Returns 'X' if not a recognized damage pattern
 */
inline char get_original_aa_ct(char damaged_aa, int codon_pos) {
    const DamageSubstitution* subs = nullptr;
    size_t count = 0;

    switch (codon_pos) {
        case 1: subs = CT_POS1_SUBS; count = CT_POS1_COUNT; break;
        case 2: subs = CT_POS2_SUBS; count = CT_POS2_COUNT; break;
        case 3: subs = CT_POS3_SUBS; count = CT_POS3_COUNT; break;
        default: return 'X';
    }

    for (size_t i = 0; i < count; ++i) {
        if (subs[i].damaged_aa == damaged_aa) {
            return subs[i].original_aa;
        }
    }
    return 'X';
}

inline char get_original_aa_ga(char damaged_aa, int codon_pos) {
    const DamageSubstitution* subs = nullptr;
    size_t count = 0;

    switch (codon_pos) {
        case 1: subs = GA_POS1_SUBS; count = GA_POS1_COUNT; break;
        case 2: subs = GA_POS2_SUBS; count = GA_POS2_COUNT; break;
        case 3: subs = GA_POS3_SUBS; count = GA_POS3_COUNT; break;
        default: return 'X';
    }

    for (size_t i = 0; i < count; ++i) {
        if (subs[i].damaged_aa == damaged_aa) {
            return subs[i].original_aa;
        }
    }
    return 'X';
}

// ============================================================================
// HEXAMER/DICODON INTEGRATION
// ============================================================================

/**
 * Compute hexamer code from 6-base sequence
 */
inline uint32_t hexamer_code(const char* seq) {
    auto base_to_int = [](char c) -> int {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    };
    uint32_t code = 0;
    for (int i = 0; i < 6; i++) {
        int b = base_to_int(seq[i]);
        if (b < 0) return UINT32_MAX;
        code = (code << 2) | static_cast<uint32_t>(b);
    }
    return code;
}

/**
 * Apply C→T damage to a hexamer at a specific position
 * Returns the damaged hexamer code, or UINT32_MAX if position has no C
 */
inline uint32_t apply_ct_damage_to_hexamer(uint32_t code, int pos) {
    // Extract base at position (0-5)
    int shift = (5 - pos) * 2;
    int base = (code >> shift) & 3;

    // C = 1, T = 3
    if (base != 1) return UINT32_MAX;  // Not a C, can't damage

    // Change C to T
    uint32_t mask = ~(3 << shift);
    return (code & mask) | (3 << shift);
}

/**
 * Apply G→A damage to a hexamer at a specific position
 */
inline uint32_t apply_ga_damage_to_hexamer(uint32_t code, int pos) {
    int shift = (5 - pos) * 2;
    int base = (code >> shift) & 3;

    // G = 2, A = 0
    if (base != 2) return UINT32_MAX;  // Not a G, can't damage

    uint32_t mask = ~(3 << shift);
    return (code & mask) | (0 << shift);
}

/**
 * Count C positions in a hexamer (potential C→T damage sites)
 */
inline int count_c_in_hexamer(uint32_t code) {
    int count = 0;
    for (int i = 0; i < 6; i++) {
        if ((code & 3) == 1) count++;  // C = 1
        code >>= 2;
    }
    return count;
}

/**
 * Count G positions in a hexamer (potential G→A damage sites)
 */
inline int count_g_in_hexamer(uint32_t code) {
    int count = 0;
    for (int i = 0; i < 6; i++) {
        if ((code & 3) == 2) count++;  // G = 2
        code >>= 2;
    }
    return count;
}

/**
 * Compute P(observed hexamer | damaged) for damage-aware scoring
 *
 * This marginalizes over all possible source hexamers that could produce
 * the observed hexamer via damage:
 *
 * P(obs | damaged) = Σ P(source) × P(source → obs | damage_rate)
 *
 * @param obs_code Observed hexamer code
 * @param damage_rate Position-specific damage rate (e.g., 0.28 at terminal)
 * @param is_ct True for C→T damage (5'), false for G→A damage (3')
 * @param hexamer_freqs Array of 4096 hexamer frequencies
 */
inline float compute_damaged_hexamer_prob(
    uint32_t obs_code,
    float damage_rate,
    bool is_ct,
    const float* hexamer_freqs)
{
    if (obs_code >= 4096) return 0.0f;

    float total_prob = 0.0f;

    // The observed hexamer could come from:
    // 1. Itself (no damage occurred at any position)
    // 2. A source hexamer where damage occurred at one or more positions

    // For simplicity, consider single-site damage (most common case)
    // Multi-site damage is rare (d^2 probability)

    // Case 1: No damage - observed is source
    float no_damage_prob = 1.0f;
    for (int i = 0; i < 6; i++) {
        int shift = (5 - i) * 2;
        int base = (obs_code >> shift) & 3;
        if (is_ct && base == 1) {  // C present - could have been damaged but wasn't
            no_damage_prob *= (1.0f - damage_rate);
        } else if (!is_ct && base == 2) {  // G present
            no_damage_prob *= (1.0f - damage_rate);
        }
    }
    total_prob += hexamer_freqs[obs_code] * no_damage_prob;

    // Case 2: Single-site damage
    for (int pos = 0; pos < 6; pos++) {
        uint32_t source_code;
        if (is_ct) {
            // Observed T at pos could be from C→T
            int shift = (5 - pos) * 2;
            int base = (obs_code >> shift) & 3;
            if (base != 3) continue;  // Not T, can't be from C→T

            // Source had C at this position
            uint32_t mask = ~(3 << shift);
            source_code = (obs_code & mask) | (1 << shift);
        } else {
            // Observed A at pos could be from G→A
            int shift = (5 - pos) * 2;
            int base = (obs_code >> shift) & 3;
            if (base != 0) continue;  // Not A, can't be from G→A

            source_code = (obs_code & ~(3 << shift)) | (2 << shift);
        }

        if (source_code < 4096) {
            // P(source) × P(damage at this site) × P(no damage at other sites)
            float transition_prob = damage_rate;
            for (int i = 0; i < 6; i++) {
                if (i == pos) continue;
                int shift = (5 - i) * 2;
                int base = (obs_code >> shift) & 3;
                if (is_ct && base == 1) {
                    transition_prob *= (1.0f - damage_rate);
                } else if (!is_ct && base == 2) {
                    transition_prob *= (1.0f - damage_rate);
                }
            }
            total_prob += hexamer_freqs[source_code] * transition_prob;
        }
    }

    return total_prob;
}

/**
 * Compute damage log-likelihood ratio for a hexamer
 *
 * LLR = log[ P(hexamer | damaged) / P(hexamer | undamaged) ]
 *
 * Positive = evidence for damage
 * Negative = evidence against damage
 */
inline float compute_hexamer_damage_llr(
    uint32_t obs_code,
    float damage_rate,
    bool is_ct,
    const float* hexamer_freqs)
{
    float p_damaged = compute_damaged_hexamer_prob(obs_code, damage_rate, is_ct, hexamer_freqs);
    float p_undamaged = hexamer_freqs[obs_code];

    if (p_undamaged < 1e-10f) return 0.0f;
    if (p_damaged < 1e-10f) return -10.0f;

    return std::log(p_damaged / p_undamaged);
}

} // namespace dart
