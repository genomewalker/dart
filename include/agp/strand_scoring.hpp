#pragma once

/**
 * @file strand_scoring.hpp
 * @brief Strand discrimination using boundary-aware hexamer LLR tables
 *
 * Uses position-specific hexamer frequencies trained on GTDB bacteria/archaea:
 * - START region (first 30bp): ATG and start codon context
 * - STOP region (last 30bp): Stop codon context
 * - INTERIOR region: General coding patterns
 *
 * For strand discrimination, compute LLR for both forward and RC readings.
 * The reading with higher LLR is more likely the coding direction.
 */

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include "gtdb_strand_hexamer.hpp"

namespace agp {

// Strand LLR tables are now included from gtdb_strand_hexamer.hpp
// (STRAND_LLR_START, STRAND_LLR_STOP, STRAND_LLR_INTERIOR, STRAND_LLR_COMBINED)

namespace strand {

// Base encoding table (A=0, C=1, G=2, T=3, other=-1)
inline int8_t encode_base(char c) {
    static const int8_t base_map[256] = {
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
    };
    return base_map[(uint8_t)c];
}

// RC base table
inline char rc_base(char c) {
    static const char rc_table[256] = {
        'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
        'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
        'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
        'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
        'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
        'N','N','N','N','A','N','N','N','N','N','N','N','N','N','N','N',
        'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
        'N','N','N','N','A','N','N','N','N','N','N','N','N','N','N','N'
    };
    return rc_table[(uint8_t)c];
}

/**
 * Compute strand log-likelihood ratio for a sequence.
 *
 * Uses position-specific tables:
 * - First 30bp: START region (ATG context)
 * - Last 30bp: STOP region (stop codon context)
 * - Middle: INTERIOR region (general coding)
 *
 * @param seq DNA sequence (uppercase A/C/G/T)
 * @param len Length of sequence
 * @return Strand LLR (positive = forward bias, negative = reverse bias)
 */
inline float compute_llr(const char* seq, size_t len) {
    if (len < 6) return 0.0f;

    float total_llr = 0.0f;
    constexpr size_t BOUNDARY = 30;

    for (size_t i = 0; i + 5 < len; i++) {
        // Encode hexamer
        uint32_t code = 0;
        bool valid = true;
        for (int j = 0; j < 6; j++) {
            int8_t b = encode_base(seq[i + j]);
            if (b < 0) { valid = false; break; }
            code = (code << 2) | b;
        }
        if (!valid) continue;

        // Select table based on position
        float llr;
        if (i < BOUNDARY) {
            llr = STRAND_LLR_START[code];
        } else if (i + 6 > len - BOUNDARY) {
            llr = STRAND_LLR_STOP[code];
        } else {
            llr = STRAND_LLR_INTERIOR[code];
        }

        total_llr += llr;
    }

    return total_llr;
}

/**
 * Compute strand LLR using ONLY interior hexamers.
 * Skips first and last SKIP_TERMINAL positions to avoid damage bias.
 * For ancient DNA, this is more reliable than using terminal regions.
 */
inline float compute_llr_interior_only(const char* seq, size_t len, size_t skip_terminal = 15) {
    if (len < skip_terminal * 2 + 6) return 0.0f;  // Too short

    float total_llr = 0.0f;
    size_t start = skip_terminal;
    size_t end = len - skip_terminal;

    for (size_t i = start; i + 5 < end; i++) {
        uint32_t code = 0;
        bool valid = true;
        for (int j = 0; j < 6; j++) {
            int8_t b = encode_base(seq[i + j]);
            if (b < 0) { valid = false; break; }
            code = (code << 2) | b;
        }
        if (!valid) continue;

        total_llr += STRAND_LLR_INTERIOR[code];
    }

    return total_llr;
}

/**
 * Compute strand LLR using only combined (position-independent) table.
 * Faster but less accurate for fragments near gene boundaries.
 */
inline float compute_llr_combined(const char* seq, size_t len) {
    if (len < 6) return 0.0f;

    float total_llr = 0.0f;

    for (size_t i = 0; i + 5 < len; i++) {
        uint32_t code = 0;
        bool valid = true;
        for (int j = 0; j < 6; j++) {
            int8_t b = encode_base(seq[i + j]);
            if (b < 0) { valid = false; break; }
            code = (code << 2) | b;
        }
        if (!valid) continue;

        total_llr += STRAND_LLR_COMBINED[code];
    }

    return total_llr;
}

/**
 * Compute strand confidence from forward and reverse LLRs.
 *
 * @param llr_fwd LLR computed on forward reading
 * @param llr_rev LLR computed on reverse complement
 * @return Confidence in [0.5, 1.0] - 0.5 = no discrimination, 1.0 = perfect
 */
inline float llr_to_confidence(float llr_fwd, float llr_rev) {
    float llr_diff = llr_fwd - llr_rev;
    // Sigmoid-like transformation
    // At llr_diff = 10, confidence ≈ 0.95
    // At llr_diff = 5, confidence ≈ 0.85
    // At llr_diff = 0, confidence = 0.5
    float conf = 0.5f + 0.5f * std::tanh(llr_diff / 10.0f);
    return conf;
}

/**
 * Compute strand LLR for forward and reverse complement readings.
 *
 * @param seq DNA sequence
 * @param len Sequence length
 * @param llr_fwd [out] LLR for forward reading
 * @param llr_rev [out] LLR for reverse complement reading
 * @param interior_only If true, skip terminal 15bp to avoid damage bias
 * @return true if sufficient data for discrimination
 */
inline bool compute_both_llrs(const char* seq, size_t len,
                              float& llr_fwd, float& llr_rev,
                              bool interior_only = false) {
    size_t min_len = interior_only ? 45 : 30;  // Need enough interior
    if (len < min_len) {
        llr_fwd = llr_rev = 0.0f;
        return false;
    }

    // Build reverse complement (stack allocate for small sequences)
    char rc_stack[1024];
    char* rc_buf = (len <= sizeof(rc_stack)) ? rc_stack : new char[len + 1];

    for (size_t i = 0; i < len; i++) {
        rc_buf[len - 1 - i] = rc_base(seq[i]);
    }
    rc_buf[len] = '\0';

    if (interior_only) {
        // Use only interior hexamers (skip damaged terminals)
        llr_fwd = compute_llr_interior_only(seq, len, 15);
        llr_rev = compute_llr_interior_only(rc_buf, len, 15);
    } else {
        // Use all hexamers with position-specific tables
        llr_fwd = compute_llr(seq, len);
        llr_rev = compute_llr(rc_buf, len);
    }

    if (len > sizeof(rc_stack)) {
        delete[] rc_buf;
    }

    return true;
}

/**
 * Predict strand orientation for a sequence.
 *
 * @param seq DNA sequence
 * @param len Sequence length
 * @param confidence [out] Confidence in prediction (0.5-1.0)
 * @return true if forward strand predicted, false if reverse
 */
inline bool predict_strand(const char* seq, size_t len, float& confidence) {
    float llr_fwd, llr_rev;
    if (!compute_both_llrs(seq, len, llr_fwd, llr_rev)) {
        confidence = 0.5f;
        return true;  // Default to forward
    }

    confidence = llr_to_confidence(llr_fwd, llr_rev);
    return llr_fwd >= llr_rev;
}

/**
 * Get strand discrimination score (positive = forward, negative = reverse).
 *
 * @param seq DNA sequence
 * @param len Sequence length
 * @return LLR difference (fwd - rev). Magnitude indicates confidence.
 */
inline float get_strand_score(const char* seq, size_t len) {
    float llr_fwd, llr_rev;
    if (!compute_both_llrs(seq, len, llr_fwd, llr_rev)) {
        return 0.0f;
    }
    return llr_fwd - llr_rev;
}

} // namespace strand
} // namespace agp
