#pragma once

/**
 * @file damage_hexamer.hpp
 * @brief Damage-aware hexamer scoring for per-read damage probability
 *
 * Uses hexamer tables to compute likelihood ratio:
 * P(observed hexamer | damaged) / P(observed hexamer | undamaged)
 *
 * Key insight: damage (C→T at 5', G→A at 3') creates predictable shifts in
 * hexamer frequencies. By modeling which "source" hexamers could produce
 * the observed hexamer via damage, we can compute a damage likelihood.
 */

#include "agp/hexamer_tables.hpp"
#include <cmath>
#include <array>
#include <algorithm>

namespace agp {

// Nucleotide encoding: A=0, C=1, G=2, T=3
constexpr int base_to_idx(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

constexpr char idx_to_base(int idx) {
    constexpr char bases[] = "ACGT";
    return (idx >= 0 && idx < 4) ? bases[idx] : 'N';
}

// Decode hexamer code to 6-base sequence
inline void decode_hexamer(uint32_t code, char* out) {
    for (int i = 5; i >= 0; --i) {
        out[i] = idx_to_base(code & 3);
        code >>= 2;
    }
    out[6] = '\0';
}

/**
 * Compute P(observed hexamer | damaged) for 5' terminal position
 *
 * At 5' end, C→T damage occurs. For each position in the hexamer:
 * - If we observe T, it could be: (1) original T, or (2) C→T damage
 * - If we observe C, it must be: original C that wasn't damaged
 *
 * P(obs_hex | damaged) = sum over all source hexamers:
 *   P(source_hex) * P(source → obs via damage transitions)
 */
inline float compute_damaged_hexamer_prob_5prime(
    uint32_t obs_code,
    float d_max,
    float lambda,
    int hex_start_pos,  // Position of hexamer start from 5' end (0-indexed)
    Domain domain = Domain::GTDB) {

    if (obs_code >= 4096) return 0.0f;

    char obs_hex[7];
    decode_hexamer(obs_code, obs_hex);

    // For each position in the hexamer, compute transition probabilities
    // Position in read = hex_start_pos + i
    float total_prob = 0.0f;

    // Enumerate all possible source hexamers (up to 2^6 = 64 if all positions could be C→T)
    // Actually, we only need to consider positions where obs has T (could be from C)

    // Count T positions in observed hexamer
    int t_positions[6];
    int n_t = 0;
    for (int i = 0; i < 6; ++i) {
        if (obs_hex[i] == 'T') {
            t_positions[n_t++] = i;
        }
    }

    // Enumerate all 2^n_t combinations of "was this T originally C?"
    int n_combinations = 1 << n_t;

    for (int mask = 0; mask < n_combinations; ++mask) {
        // Build source hexamer
        char src_hex[7];
        for (int i = 0; i < 6; ++i) src_hex[i] = obs_hex[i];
        src_hex[6] = '\0';

        // For each T position, decide if it was originally C
        float transition_prob = 1.0f;
        for (int j = 0; j < n_t; ++j) {
            int pos_in_hex = t_positions[j];
            int pos_in_read = hex_start_pos + pos_in_hex;
            float d_pos = d_max * std::exp(-lambda * pos_in_read);

            if (mask & (1 << j)) {
                // This T was originally C (damaged)
                src_hex[pos_in_hex] = 'C';
                transition_prob *= d_pos;  // P(C→T)
            }
            // Else: T stays T with probability 1.0 (no change to transition_prob)
        }

        // For C positions in observed, they must be undamaged C
        for (int i = 0; i < 6; ++i) {
            if (obs_hex[i] == 'C') {
                int pos_in_read = hex_start_pos + i;
                float d_pos = d_max * std::exp(-lambda * pos_in_read);
                transition_prob *= (1.0f - d_pos);  // P(C stays C)
            }
        }

        // Get source hexamer frequency
        uint32_t src_code = encode_hexamer(src_hex);
        if (src_code == UINT32_MAX) continue;

        float src_freq = get_hexamer_freq(src_code, domain);
        total_prob += src_freq * transition_prob;
    }

    return total_prob;
}

/**
 * Compute P(observed hexamer | damaged) for 3' terminal position
 *
 * At 3' end, G→A damage occurs (on forward strand).
 */
inline float compute_damaged_hexamer_prob_3prime(
    uint32_t obs_code,
    float d_max,
    float lambda,
    int hex_dist_from_3prime,  // Distance of hexamer END from 3' end
    Domain domain = Domain::GTDB) {

    if (obs_code >= 4096) return 0.0f;

    char obs_hex[7];
    decode_hexamer(obs_code, obs_hex);

    float total_prob = 0.0f;

    // Count A positions (could be from G→A)
    int a_positions[6];
    int n_a = 0;
    for (int i = 0; i < 6; ++i) {
        if (obs_hex[i] == 'A') {
            a_positions[n_a++] = i;
        }
    }

    int n_combinations = 1 << n_a;

    for (int mask = 0; mask < n_combinations; ++mask) {
        char src_hex[7];
        for (int i = 0; i < 6; ++i) src_hex[i] = obs_hex[i];
        src_hex[6] = '\0';

        float transition_prob = 1.0f;
        for (int j = 0; j < n_a; ++j) {
            int pos_in_hex = a_positions[j];
            // Distance from 3' end: hex is at hex_dist_from_3prime, position within is (5 - pos_in_hex)
            int dist_from_3prime = hex_dist_from_3prime + (5 - pos_in_hex);
            float d_pos = d_max * std::exp(-lambda * dist_from_3prime);

            if (mask & (1 << j)) {
                // This A was originally G (damaged)
                src_hex[pos_in_hex] = 'G';
                transition_prob *= d_pos;
            }
        }

        // For G positions, they must be undamaged
        for (int i = 0; i < 6; ++i) {
            if (obs_hex[i] == 'G') {
                int dist_from_3prime = hex_dist_from_3prime + (5 - i);
                float d_pos = d_max * std::exp(-lambda * dist_from_3prime);
                transition_prob *= (1.0f - d_pos);
            }
        }

        uint32_t src_code = encode_hexamer(src_hex);
        if (src_code == UINT32_MAX) continue;

        float src_freq = get_hexamer_freq(src_code, domain);
        total_prob += src_freq * transition_prob;
    }

    return total_prob;
}

/**
 * Compute damage log-likelihood ratio for terminal hexamers
 *
 * LLR = log[ P(observed | damaged) / P(observed | undamaged) ]
 *
 * Positive LLR = evidence for damage
 * Negative LLR = evidence against damage
 */
inline float compute_damage_hexamer_llr(
    const std::string& seq,
    float d_max,
    float lambda,
    int terminal_window = 15,
    Domain domain = Domain::GTDB) {

    if (seq.length() < 12) return 0.0f;
    if (d_max < 0.01f) return 0.0f;

    float total_llr = 0.0f;
    int hex_count = 0;

    const size_t len = seq.length();

    // 5' terminal hexamers (C→T damage)
    for (size_t i = 0; i + 5 < len && i < static_cast<size_t>(terminal_window); i += 3) {
        uint32_t code = encode_hexamer(seq.c_str() + i);
        if (code == UINT32_MAX) continue;

        float p_undamaged = get_hexamer_freq(code, domain);
        float p_damaged = compute_damaged_hexamer_prob_5prime(code, d_max, lambda, i, domain);

        if (p_undamaged > 1e-10f && p_damaged > 1e-10f) {
            total_llr += std::log(p_damaged / p_undamaged);
            hex_count++;
        }
    }

    // 3' terminal hexamers (G→A damage)
    for (size_t dist = 0; dist + 5 < len && dist < static_cast<size_t>(terminal_window); dist += 3) {
        size_t pos = len - 6 - dist;
        if (pos >= len) continue;

        uint32_t code = encode_hexamer(seq.c_str() + pos);
        if (code == UINT32_MAX) continue;

        float p_undamaged = get_hexamer_freq(code, domain);
        float p_damaged = compute_damaged_hexamer_prob_3prime(code, d_max, lambda, dist, domain);

        if (p_undamaged > 1e-10f && p_damaged > 1e-10f) {
            total_llr += std::log(p_damaged / p_undamaged);
            hex_count++;
        }
    }

    return total_llr;
}

/**
 * Damage-consistent amino acid changes
 *
 * These AA substitutions can be caused by C→T or G→A damage:
 */
struct DamageAAChange {
    char from_aa;
    char to_aa;
    float weight;  // How diagnostic this change is (higher = more specific to damage)
};

// C→T induced changes (5' terminal)
constexpr DamageAAChange CT_AA_CHANGES[] = {
    {'Q', '*', 2.0f},  // CAA→TAA, CAG→TAG (Gln → Stop) - VERY diagnostic
    {'R', '*', 2.0f},  // CGA→TGA (Arg → Stop)
    {'W', '*', 2.0f},  // TGG→TGA (Trp → Stop) - actually G→A at pos 2
    {'S', 'L', 1.0f},  // TCA→TTA, TCG→TTG (Ser → Leu)
    {'P', 'L', 1.0f},  // CCA→CTA, CCG→CTG (Pro → Leu)
    {'T', 'I', 1.0f},  // ACA→ATA (Thr → Ile)
    {'T', 'M', 1.0f},  // ACG→ATG (Thr → Met)
    {'A', 'V', 1.0f},  // GCA→GTA, GCG→GTG (Ala → Val)
    {'S', 'F', 0.8f},  // TCC→TTC (Ser → Phe)
    {'P', 'S', 0.8f},  // CCT→TCT, CCC→TCC (Pro → Ser)
    {'H', 'Y', 0.8f},  // CAT→TAT, CAC→TAC (His → Tyr)
    {'R', 'C', 0.8f},  // CGT→TGT, CGC→TGC (Arg → Cys)
};

// G→A induced changes (3' terminal)
constexpr DamageAAChange GA_AA_CHANGES[] = {
    {'W', '*', 2.0f},  // TGG→TGA (Trp → Stop)
    {'G', 'D', 1.0f},  // GGT→GAT, GGC→GAC (Gly → Asp)
    {'G', 'E', 1.0f},  // GGA→GAA, GGG→GAG (Gly → Glu)
    {'G', 'S', 0.8f},  // GGT→AGT, GGC→AGC (Gly → Ser)
    {'A', 'T', 0.8f},  // GCT→ACT, GCC→ACC (Ala → Thr)
    {'E', 'K', 0.8f},  // GAA→AAA, GAG→AAG (Glu → Lys)
    {'D', 'N', 0.8f},  // GAT→AAT, GAC→AAC (Asp → Asn)
    {'R', 'K', 0.8f},  // AGA→AAA, AGG→AAG (Arg → Lys)
    {'S', 'N', 0.8f},  // AGT→AAT, AGC→AAC (Ser → Asn)
    {'V', 'I', 0.5f},  // GTT→ATT, GTC→ATC (Val → Ile)
    {'V', 'M', 0.5f},  // GTG→ATG (Val → Met)
};

/**
 * Score amino acid composition for damage evidence
 *
 * At 5' terminal: look for T-containing codons (could be from C→T)
 * At 3' terminal: look for A-containing codons (could be from G→A)
 */
inline float compute_damage_aa_llr(
    const std::string& protein,
    int terminal_aa = 5,
    bool is_5prime = true) {

    if (protein.length() < static_cast<size_t>(terminal_aa)) return 0.0f;

    float llr = 0.0f;

    // Amino acids enriched by damage
    // 5' (C→T): L, F, I, M, V, Y (T-containing codons more likely)
    // 3' (G→A): D, E, N, K (A-containing codons more likely)

    if (is_5prime) {
        // Check first terminal_aa amino acids
        for (int i = 0; i < terminal_aa && i < static_cast<int>(protein.length()); ++i) {
            char aa = protein[i];
            float weight = std::exp(-0.3f * i);  // Exponential decay from terminal

            // Damage-enriched AAs at 5' (from T-rich codons)
            switch (aa) {
                case 'L': llr += 0.15f * weight; break;  // TTA, TTG, CTN
                case 'F': llr += 0.20f * weight; break;  // TTT, TTC
                case 'I': llr += 0.15f * weight; break;  // ATT, ATC, ATA
                case 'M': llr += 0.10f * weight; break;  // ATG
                case 'V': llr += 0.10f * weight; break;  // GTN (from GCN via damage)
                case 'Y': llr += 0.15f * weight; break;  // TAT, TAC
                case '*': llr += 0.50f * weight; break;  // Stop codon - strong signal
                // Damage-depleted AAs (C-containing codons less likely after damage)
                case 'P': llr -= 0.10f * weight; break;  // CCN
                case 'H': llr -= 0.05f * weight; break;  // CAT, CAC
                case 'Q': llr -= 0.05f * weight; break;  // CAA, CAG
            }
        }
    } else {
        // Check last terminal_aa amino acids
        int start = std::max(0, static_cast<int>(protein.length()) - terminal_aa);
        for (int i = protein.length() - 1; i >= start; --i) {
            char aa = protein[i];
            float weight = std::exp(-0.3f * (protein.length() - 1 - i));

            // Damage-enriched AAs at 3' (from A-rich codons)
            switch (aa) {
                case 'D': llr += 0.15f * weight; break;  // GAT, GAC
                case 'E': llr += 0.15f * weight; break;  // GAA, GAG
                case 'N': llr += 0.15f * weight; break;  // AAT, AAC
                case 'K': llr += 0.15f * weight; break;  // AAA, AAG
                case 'I': llr += 0.10f * weight; break;  // ATT, ATC, ATA
                case '*': llr += 0.50f * weight; break;  // Stop codon
                // Damage-depleted (G-containing codons)
                case 'G': llr -= 0.10f * weight; break;  // GGN
                case 'A': llr -= 0.05f * weight; break;  // GCN
            }
        }
    }

    return llr;
}

} // namespace agp
