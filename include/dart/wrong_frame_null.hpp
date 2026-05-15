#pragma once
#include "dart/gtdb_hexamer_table.hpp"
#include <array>
#include <cmath>

// Wrong-frame null model for contrastive ORF scoring.
//
// Coding sequence read at a +1 or +2 frame shift produces a hexamer distribution
// that is NOT random — it is structured by the underlying codon usage. Scoring
// candidate frames against this null (instead of a uniform null) makes the score
// contrastive: a wrong frame that happens to look coding by coincidence scores near
// 0 rather than high.
//
// WRONG_FRAME_LLR[hex] = log P(hex | coding) - log P(hex | +1/+2 shift null)
//   Positive: hex is more likely in a true coding frame than in a shifted one.
//   Negative: hex is more common in frameshifted sequence.
//
// The null is derived from GTDB_HEXAMER_FREQ via a first-order codon Markov model:
//   F[c0][c1]  = GTDB_HEXAMER_FREQ[encode(c0, c1)]  (dicodon frequency)
//   M[c1]      = sum_{c0} F[c0][c1]                 (codon marginal)
//   T[c1][c2]  = F[c1][c2] / M[c1]                  (codon transition)
//   prob(c0,c1,c2) = F[c0][c1] * T[c1][c2]
//   plus1[hex(s[1..6])] += prob   (+1 shift hexamer)
//   plus2[hex(s[2..7])] += prob   (+2 shift hexamer)
//   null = 0.5 * plus1 + 0.5 * plus2

namespace dart {

inline std::array<float, 4096> compute_wrong_frame_null() {
    std::array<float, 4096> plus1{}, plus2{}, M{};

    for (int c1 = 0; c1 < 64; c1++)
        for (int c0 = 0; c0 < 64; c0++)
            M[c1] += GTDB_HEXAMER_FREQ[c0 * 64 + c1];

    for (int c0 = 0; c0 < 64; c0++) {
        for (int c1 = 0; c1 < 64; c1++) {
            float f01 = GTDB_HEXAMER_FREQ[c0 * 64 + c1];
            if (f01 < 1e-12f) continue;
            for (int c2 = 0; c2 < 64; c2++) {
                float f12 = GTDB_HEXAMER_FREQ[c1 * 64 + c2];
                if (f12 < 1e-12f || M[c1] < 1e-12f) continue;
                float prob = f01 * (f12 / M[c1]);
                if (prob < 1e-15f) continue;

                // Decode each codon to 3 bases (2 bits each: A=0,C=1,G=2,T=3)
                int b[9];
                b[0]=(c0>>4)&3; b[1]=(c0>>2)&3; b[2]=c0&3;
                b[3]=(c1>>4)&3; b[4]=(c1>>2)&3; b[5]=c1&3;
                b[6]=(c2>>4)&3; b[7]=(c2>>2)&3; b[8]=c2&3;

                // +1 shift: bases 1..6 → 12-bit code
                int h1 = b[1]*1024 + b[2]*256 + b[3]*64 + b[4]*16 + b[5]*4 + b[6];
                plus1[h1] += prob;

                // +2 shift: bases 2..7 → 12-bit code
                int h2 = b[2]*1024 + b[3]*256 + b[4]*64 + b[5]*16 + b[6]*4 + b[7];
                plus2[h2] += prob;
            }
        }
    }

    std::array<float, 4096> null_dist{};
    for (int h = 0; h < 4096; h++)
        null_dist[h] = 0.5f * plus1[h] + 0.5f * plus2[h];
    return null_dist;
}

inline const std::array<float, 4096>& get_wrong_frame_llr() {
    static const std::array<float, 4096> llr = []() {
        auto null = compute_wrong_frame_null();
        std::array<float, 4096> r{};
        for (int h = 0; h < 4096; h++)
            r[h] = std::log(GTDB_HEXAMER_FREQ[h] + 1e-10f)
                 - std::log(null[h]              + 1e-10f);
        return r;
    }();
    return llr;
}

}  // namespace dart
