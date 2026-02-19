#pragma once

/**
 * @file forbidden_patterns.hpp
 * @brief Unorthodox frame scoring for frame discrimination
 *
 * This file contains multiple UNORTHODOX approaches for discriminating
 * between reading frames when traditional methods (stop codons) fail.
 *
 * KEY INSIGHT: ~40% of short reads have 0 stop codons in ALL 6 frames.
 * In these cases, we need alternative signals.
 *
 * APPROACHES IMPLEMENTED:
 * 1. Forbidden patterns: Stop-containing hexamers at codon boundaries
 * 2. Nucleotide periodicity: Position-specific nucleotide frequencies
 * 3. Rare codon avoidance: Penalize frames with many rare codons
 * 4. Dinucleotide bias: CpG depletion and other biases
 * 5. RNY rule: Purine-aNy-pYrimidine pattern in coding
 * 6. Amino acid composition: Penalize frames with unusual AA frequencies
 * 7. Low-complexity penalty: Penalize repetitive sequences
 */

#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include "agp/hexamer_tables.hpp"

namespace agp {
namespace forbidden {

// Lookup table for nucleotide to index conversion (A=0, C=1, G=2, T=3, other=-1)
// Handles both upper and lowercase efficiently
alignas(64) constexpr int8_t NT_LOOKUP[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 0-15
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 16-31
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 32-47
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 48-63
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1, // 64-79: A=65, C=67, G=71
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 80-95: T=84
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1, // 96-111: a=97, c=99, g=103
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 112-127: t=116
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 128-143
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 144-159
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 160-175
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 176-191
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 192-207
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 208-223
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 224-239
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1  // 240-255
};

// Rare codon flags (indexed by codon: b0*16 + b1*4 + b2)
// AGG=0*16+2*4+2=10, AGA=0*16+2*4+0=8, CGA=1*16+2*4+0=24, CGG=1*16+2*4+2=26, CTA=1*16+3*4+0=28, ATA=0*16+3*4+0=12
alignas(64) constexpr uint8_t RARE_CODON[64] = {
    0,0,0,0, 0,0,0,0, 1,0,1,0, 1,0,0,0,  // 0-15: AGA(8), AGG(10), ATA(12)
    0,0,0,0, 0,0,0,0, 1,0,1,0, 1,0,0,0,  // 16-31: CGA(24), CGG(26), CTA(28)
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  // 32-47
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0   // 48-63
};

//==============================================================================
// 1. FORBIDDEN PATTERNS (Stop-containing hexamers)
//==============================================================================

/**
 * Check if hexamer starts with a stop codon (TAA, TAG, TGA)
 */
inline bool is_forbidden_start(const char* hex) {
    if (hex[0] == 'T' || hex[0] == 't') {
        char c1 = (hex[1] >= 'a') ? (hex[1] - 32) : hex[1];
        char c2 = (hex[2] >= 'a') ? (hex[2] - 32) : hex[2];
        if (c1 == 'A' && (c2 == 'A' || c2 == 'G')) return true;  // TAA, TAG
        if (c1 == 'G' && c2 == 'A') return true;  // TGA
    }
    return false;
}

/**
 * Check if hexamer has stop codon at position 3-5
 */
inline bool is_forbidden_internal(const char* hex) {
    if (hex[3] == 'T' || hex[3] == 't') {
        char c1 = (hex[4] >= 'a') ? (hex[4] - 32) : hex[4];
        char c2 = (hex[5] >= 'a') ? (hex[5] - 32) : hex[5];
        if (c1 == 'A' && (c2 == 'A' || c2 == 'G')) return true;
        if (c1 == 'G' && c2 == 'A') return true;
    }
    return false;
}

/**
 * Count forbidden patterns (stop-containing hexamers at codon boundaries)
 */
inline int count_forbidden_patterns(const char* seq, size_t len, int frame) {
    if (len < 6 || frame < 0 || frame > 2) return 0;
    
    int forbidden_count = 0;
    for (size_t i = frame; i + 5 < len; i += 3) {
        if (is_forbidden_start(seq + i)) forbidden_count++;
        if (is_forbidden_internal(seq + i)) forbidden_count++;
    }
    return forbidden_count;
}

/**
 * Calculate forbidden pattern score (0 = many forbidden, 1 = none)
 */
inline float calculate_forbidden_score(const char* seq, size_t len, int frame) {
    int forbidden = count_forbidden_patterns(seq, len, frame);
    if (forbidden == 0) return 1.0f;
    return std::exp(-3.0f * forbidden);  // Softer penalty
}

//==============================================================================
// 2. NUCLEOTIDE PERIODICITY (Position-specific frequencies)
//==============================================================================

/**
 * Calculate nucleotide periodicity score based on codon position preferences.
 *
 * In coding regions, nucleotides show 3-periodic patterns:
 * - Position 0 (1st): High purine (A+G ~62%), especially G
 * - Position 1 (2nd): More balanced
 * - Position 2 (wobble): High GC (~65%), especially C
 */
inline float calculate_periodicity_score(const char* seq, size_t len, int frame) {
    if (len < 9 || frame < 0 || frame > 2) return 0.5f;

    int nt_counts[3][4] = {{0}};  // [codon_pos][A=0,C=1,G=2,T=3]
    int total[3] = {0};

    // Rolling position counter - avoids expensive modulo per iteration
    int pos = (3 - (frame % 3)) % 3;  // Starting codon position

    for (size_t i = 0; i < len; i++) {
        int8_t nt = NT_LOOKUP[static_cast<uint8_t>(seq[i])];
        if (nt >= 0) {
            nt_counts[pos][nt]++;
            total[pos]++;
        }
        if (++pos == 3) pos = 0;
    }

    if (total[0] < 3 || total[1] < 3 || total[2] < 3) return 0.5f;

    float purine_0 = static_cast<float>(nt_counts[0][0] + nt_counts[0][2]) / total[0];
    float purine_1 = static_cast<float>(nt_counts[1][0] + nt_counts[1][2]) / total[1];
    float gc_2 = static_cast<float>(nt_counts[2][1] + nt_counts[2][2]) / total[2];
    float g_0 = static_cast<float>(nt_counts[0][2]) / total[0];
    float c_2 = static_cast<float>(nt_counts[2][1]) / total[2];

    float score = 0.0f;
    score += (purine_0 - 0.5f) * 0.4f;
    score += (purine_0 - purine_1) * 0.3f;
    score += (gc_2 - 0.5f) * 0.4f;
    score += (c_2 - g_0) * 0.2f;

    return std::max(0.0f, std::min(1.0f, 0.5f + score));
}

//==============================================================================
// 3. RARE CODON AVOIDANCE
//==============================================================================

/**
 * Count rare codons in a frame using numeric lookup.
 *
 * Rare codons in bacteria: AGG, AGA, CGA, CGG (Arg), CTA (Leu), ATA (Ile)
 */
inline int count_rare_codons(const char* seq, size_t len, int frame) {
    if (len < 3 + static_cast<size_t>(frame)) return 0;

    int count = 0;
    for (size_t i = frame; i + 2 < len; i += 3) {
        int8_t b0 = NT_LOOKUP[static_cast<uint8_t>(seq[i])];
        int8_t b1 = NT_LOOKUP[static_cast<uint8_t>(seq[i+1])];
        int8_t b2 = NT_LOOKUP[static_cast<uint8_t>(seq[i+2])];

        if (b0 >= 0 && b1 >= 0 && b2 >= 0) {
            int idx = (b0 << 4) | (b1 << 2) | b2;
            count += RARE_CODON[idx];
        }
    }
    return count;
}

/**
 * Calculate rare codon penalty score (1 = no rare, 0 = many rare)
 */
inline float calculate_rare_codon_score(const char* seq, size_t len, int frame) {
    int rare = count_rare_codons(seq, len, frame);
    if (rare == 0) return 1.0f;
    int num_codons = (len - frame) / 3;
    float rare_frac = static_cast<float>(rare) / std::max(1, num_codons);
    // Expect ~2% rare codons in random sequence, ~0.5% in coding
    return std::exp(-10.0f * rare_frac);
}

//==============================================================================
// 4. DINUCLEOTIDE BIAS (CpG depletion, etc.)
//==============================================================================

/**
 * Calculate dinucleotide bias score.
 * 
 * Coding regions have specific dinucleotide patterns:
 * - CpG is depleted in most organisms (methylation)
 * - TA and TT at certain positions are avoided (ribosome stalling)
 */
inline float calculate_dinucleotide_score(const char* seq, size_t len, int frame) {
    if (len < 6) return 0.5f;
    
    // Count dinucleotides at codon boundaries (pos 2-3, spanning wobble to next codon)
    int cg_at_boundary = 0;  // CG at positions 2-0 (wobble to next first)
    int ta_at_boundary = 0;  // TA at positions 2-0
    int boundary_count = 0;
    
    for (size_t i = frame + 2; i + 1 < len; i += 3) {
        char c1 = (seq[i] >= 'a') ? (seq[i] - 32) : seq[i];
        char c2 = (seq[i+1] >= 'a') ? (seq[i+1] - 32) : seq[i+1];
        
        boundary_count++;
        if (c1 == 'C' && c2 == 'G') cg_at_boundary++;
        if (c1 == 'T' && c2 == 'A') ta_at_boundary++;
    }
    
    if (boundary_count < 3) return 0.5f;
    
    float cg_freq = static_cast<float>(cg_at_boundary) / boundary_count;
    float ta_freq = static_cast<float>(ta_at_boundary) / boundary_count;
    
    // CpG expected ~1% in coding (depleted from ~6% in random)
    // TA expected ~4% in coding
    float score = 0.5f;
    
    // Penalize high CpG (wrong frame may have more CpG at "boundaries")
    score -= cg_freq * 0.5f;
    
    // Slight penalty for high TA
    score -= ta_freq * 0.1f;
    
    return std::max(0.0f, std::min(1.0f, score + 0.1f));
}

//==============================================================================
// 5. RNY RULE (Purine-aNy-pYrimidine pattern)
//==============================================================================

/**
 * Calculate RNY rule compliance.
 * 
 * The RNY rule: codons often have pattern R-N-Y where:
 * - R = purine (A or G) at position 0
 * - N = any nucleotide at position 1
 * - Y = pyrimidine (C or T) at position 2
 * 
 * This is a weak signal but consistent across many organisms.
 */
inline float calculate_rny_score(const char* seq, size_t len, int frame) {
    if (len < 3 + static_cast<size_t>(frame)) return 0.5f;
    
    int rny_count = 0;
    int codon_count = 0;
    
    for (size_t i = frame; i + 2 < len; i += 3) {
        char c0 = (seq[i] >= 'a') ? (seq[i] - 32) : seq[i];
        char c2 = (seq[i+2] >= 'a') ? (seq[i+2] - 32) : seq[i+2];
        
        codon_count++;
        
        // Check RNY pattern
        bool is_R = (c0 == 'A' || c0 == 'G');
        bool is_Y = (c2 == 'C' || c2 == 'T');
        
        if (is_R && is_Y) rny_count++;
    }
    
    if (codon_count < 3) return 0.5f;
    
    float rny_frac = static_cast<float>(rny_count) / codon_count;
    
    // Expected ~25% random (0.5 * 0.5), ~35-40% in coding
    // Score higher for above-random RNY compliance
    return 0.5f + (rny_frac - 0.25f) * 1.5f;
}

//==============================================================================
// 6. AMINO ACID COMPOSITION PATTERNS
//==============================================================================

/**
 * Quick codon to amino acid (returns single char, 'X' for unknown, '*' for stop)
 */
inline char quick_translate(char c0, char c1, char c2) {
    // Uppercase
    c0 = (c0 >= 'a') ? (c0 - 32) : c0;
    c1 = (c1 >= 'a') ? (c1 - 32) : c1;
    c2 = (c2 >= 'a') ? (c2 - 32) : c2;
    
    if (c0 == 'T') {
        if (c1 == 'T') return (c2 == 'T' || c2 == 'C') ? 'F' : 'L';
        if (c1 == 'C') return 'S';
        if (c1 == 'A') return (c2 == 'T' || c2 == 'C') ? 'Y' : '*';
        if (c1 == 'G') return (c2 == 'T' || c2 == 'C') ? 'C' : ((c2 == 'A') ? '*' : 'W');
    }
    if (c0 == 'C') {
        if (c1 == 'T') return 'L';
        if (c1 == 'C') return 'P';
        if (c1 == 'A') return (c2 == 'T' || c2 == 'C') ? 'H' : 'Q';
        if (c1 == 'G') return 'R';
    }
    if (c0 == 'A') {
        if (c1 == 'T') return (c2 == 'G') ? 'M' : 'I';
        if (c1 == 'C') return 'T';
        if (c1 == 'A') return (c2 == 'T' || c2 == 'C') ? 'N' : 'K';
        if (c1 == 'G') return (c2 == 'T' || c2 == 'C') ? 'S' : 'R';
    }
    if (c0 == 'G') {
        if (c1 == 'T') return 'V';
        if (c1 == 'C') return 'A';
        if (c1 == 'A') return (c2 == 'T' || c2 == 'C') ? 'D' : 'E';
        if (c1 == 'G') return 'G';
    }
    return 'X';
}

/**
 * Calculate amino acid composition score.
 * 
 * Some amino acids are common (L, A, G, V, E) and others rare (W, C, M).
 * Wrong frames produce unusual compositions.
 */
inline float calculate_aa_composition_score(const char* seq, size_t len, int frame) {
    if (len < 3 + static_cast<size_t>(frame)) return 0.5f;
    
    // Count amino acids and check unusual composition patterns
    int total = 0;
    int common = 0;  // L, A, G, V, E, S, I, T (>5% each in most proteomes)
    int rare = 0;    // W, C, M, H (typically <3% each)
    int stop = 0;
    
    for (size_t i = frame; i + 2 < len; i += 3) {
        char aa = quick_translate(seq[i], seq[i+1], seq[i+2]);
        total++;
        
        switch (aa) {
            case 'L': case 'A': case 'G': case 'V': case 'E':
            case 'S': case 'I': case 'T': case 'R': case 'K':
                common++; break;
            case 'W': case 'C': case 'M': case 'H':
                rare++; break;
            case '*':
                stop++; break;
        }
    }
    
    if (total < 3) return 0.5f;
    
    float common_frac = static_cast<float>(common) / total;
    float rare_frac = static_cast<float>(rare) / total;
    float stop_frac = static_cast<float>(stop) / total;
    
    // Expected: common ~70%, rare ~8%, stops ~0%
    float score = 0.5f;
    score += (common_frac - 0.5f) * 0.3f;
    score -= rare_frac * 0.5f;
    score -= stop_frac * 2.0f;  // Strong penalty for internal stops
    
    return std::max(0.0f, std::min(1.0f, score + 0.2f));
}

//==============================================================================
// 7. LOW-COMPLEXITY PENALTY
//==============================================================================

/**
 * Calculate low-complexity penalty.
 * 
 * Repetitive sequences (AAAA, TTTT, etc.) are often non-coding or artifacts.
 * Penalize frames that have these at codon boundaries.
 */
inline float calculate_complexity_score(const char* seq, size_t len, int frame) {
    if (len < 9) return 1.0f;
    
    int homopolymer_codons = 0;
    int total_codons = 0;
    
    for (size_t i = frame; i + 2 < len; i += 3) {
        total_codons++;
        char c0 = (seq[i] >= 'a') ? (seq[i] - 32) : seq[i];
        char c1 = (seq[i+1] >= 'a') ? (seq[i+1] - 32) : seq[i+1];
        char c2 = (seq[i+2] >= 'a') ? (seq[i+2] - 32) : seq[i+2];
        
        // Check for homopolymer codon (AAA, TTT, CCC, GGG)
        if (c0 == c1 && c1 == c2) {
            homopolymer_codons++;
        }
    }
    
    if (total_codons < 3) return 1.0f;
    
    float homo_frac = static_cast<float>(homopolymer_codons) / total_codons;
    
    // Expected: ~1.5% random (4/64), penalize if higher
    return std::exp(-20.0f * homo_frac);
}

//==============================================================================
// 8. HEXAMER RARITY SCORE (use hexamer table to find unusual patterns)
//==============================================================================

/**
 * Calculate hexamer rarity score using the hexamer frequency table.
 * 
 * Wrong frames will have more rare (low-frequency) hexamers than correct frames.
 */
inline float calculate_hexamer_rarity_score(const char* seq, size_t len, int frame) {
    if (len < static_cast<size_t>(frame + 6)) return 0.5f;
    
    float log_prob_sum = 0.0f;
    int hexamer_count = 0;
    
    char hexbuf[7];
    hexbuf[6] = '\0';
    
    for (size_t i = frame; i + 5 < len; i += 3) {
        for (int j = 0; j < 6; j++) {
            char c = seq[i + j];
            hexbuf[j] = (c >= 'a') ? (c - 32) : c;
        }
        
        // Skip if any N
        bool has_n = false;
        for (int j = 0; j < 6; j++) {
            if (hexbuf[j] == 'N') { has_n = true; break; }
        }
        if (has_n) continue;
        
        uint32_t code = encode_hexamer(hexbuf);
        if (code == UINT32_MAX) continue;
        
        float freq = get_hexamer_freq(code);
        if (freq < 1e-8f) freq = 1e-8f;
        
        log_prob_sum += std::log(freq);
        hexamer_count++;
    }
    
    if (hexamer_count == 0) return 0.5f;
    
    // Average log probability per hexamer
    float avg_log_prob = log_prob_sum / hexamer_count;
    
    // Typical range: -8 (rare) to -5 (common)
    // Map to 0-1: -8 -> 0.2, -5 -> 0.8
    float normalized = (avg_log_prob + 8.0f) / 4.0f;
    return std::max(0.0f, std::min(1.0f, 0.3f + normalized * 0.5f));
}

//==============================================================================
// COMBINED SCORING
//==============================================================================

/**
 * Structure to hold all frame discrimination scores
 */
struct FrameScores {
    float forbidden;      // Stop-containing hexamers
    float periodicity;    // Nucleotide position bias
    float rare_codon;     // Rare codon avoidance
    float dinucleotide;   // CpG depletion etc.
    float rny;            // RNY rule compliance
    float aa_comp;        // Amino acid composition
    float complexity;     // Low-complexity penalty
    float hexamer_rarity; // Hexamer frequency
    
    float combined;       // Weighted combination
};

/**
 * Calculate all frame discrimination scores
 */
inline FrameScores calculate_all_scores(const char* seq, size_t len, int frame) {
    FrameScores s;
    s.forbidden = calculate_forbidden_score(seq, len, frame);
    s.periodicity = calculate_periodicity_score(seq, len, frame);
    s.rare_codon = calculate_rare_codon_score(seq, len, frame);
    s.dinucleotide = calculate_dinucleotide_score(seq, len, frame);
    s.rny = calculate_rny_score(seq, len, frame);
    s.aa_comp = calculate_aa_composition_score(seq, len, frame);
    s.complexity = calculate_complexity_score(seq, len, frame);
    s.hexamer_rarity = calculate_hexamer_rarity_score(seq, len, frame);
    
    // Weighted combination (tuned for maximum discrimination)
    // Higher weights for more reliable signals
    s.combined = 
        0.20f * s.forbidden +      // Stop patterns are highly discriminative
        0.15f * s.periodicity +    // Position-specific nucleotide bias
        0.10f * s.rare_codon +     // Rare codon avoidance
        0.08f * s.dinucleotide +   // Dinucleotide bias
        0.12f * s.rny +            // RNY rule
        0.15f * s.aa_comp +        // Amino acid composition
        0.05f * s.complexity +     // Complexity penalty
        0.15f * s.hexamer_rarity;  // Hexamer rarity
    
    return s;
}

/**
 * Calculate combined unorthodox frame score (single value, 0-1)
 */
inline float calculate_unorthodox_frame_score(const char* seq, size_t len, int frame) {
    return calculate_all_scores(seq, len, frame).combined;
}

/**
 * Select best frame using all unorthodox signals
 * Returns: (best_frame, confidence)
 */
inline std::pair<int, float> select_best_frame_unorthodox(const char* seq, size_t len) {
    float scores[3];
    scores[0] = calculate_unorthodox_frame_score(seq, len, 0);
    scores[1] = calculate_unorthodox_frame_score(seq, len, 1);
    scores[2] = calculate_unorthodox_frame_score(seq, len, 2);
    
    // Best frame has highest score
    int best = 0;
    if (scores[1] > scores[best]) best = 1;
    if (scores[2] > scores[best]) best = 2;
    
    // Confidence: difference from second best
    float sorted[3] = {scores[0], scores[1], scores[2]};
    std::sort(sorted, sorted + 3);
    
    float confidence = sorted[2] - sorted[1];  // Best minus second best
    
    return {best, confidence};
}

} // namespace forbidden
} // namespace agp
