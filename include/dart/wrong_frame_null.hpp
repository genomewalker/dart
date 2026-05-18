#pragma once
#include "dart/gtdb_hexamer_table.hpp"
#include <array>
#include <cmath>

// Wrong-frame null model for contrastive ORF scoring.
//
// Scores each hexamer as log P(hex|coding) - log P(hex|wrong-frame null).
// The null covers all 5 wrong-frame competitors:
//   +1, +2 forward frameshifts of coding sequence
//   RC frame 0, 1, 2 — reverse-complement of coding sequence read in each phase
//
// Oracle tests showed RC competitors are the dominant source of top-1 errors on
// damage-stop reads (+25% accuracy gain from removing RC competitors entirely).
// Including RC distributions in the null de-scores RC frames directly.
//
// The null is derived from GTDB_HEXAMER_FREQ via a first-order codon Markov model.

namespace dart {

// RC complement lookup: A↔T (0↔3), C↔G (1↔2)
inline constexpr int rc_base(int b) { return 3 - b; }

// Encode 6-base hexamer from individual base values (2-bit each)
inline int encode6(int b0, int b1, int b2, int b3, int b4, int b5) {
    return b0*1024 + b1*256 + b2*64 + b3*16 + b4*4 + b5;
}

inline std::array<float, 4096> compute_wrong_frame_null() {
    std::array<float, 4096> plus1{}, plus2{}, rev0{}, rev1{}, rev2{}, M{};

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

                // Decode 9 bases of the codon triple (A=0,C=1,G=2,T=3)
                int b[9];
                b[0]=(c0>>4)&3; b[1]=(c0>>2)&3; b[2]=c0&3;
                b[3]=(c1>>4)&3; b[4]=(c1>>2)&3; b[5]=c1&3;
                b[6]=(c2>>4)&3; b[7]=(c2>>2)&3; b[8]=c2&3;

                // Forward frameshifts: bases 1..6 and 2..7
                plus1[encode6(b[1],b[2],b[3],b[4],b[5],b[6])] += prob;
                plus2[encode6(b[2],b[3],b[4],b[5],b[6],b[7])] += prob;

                // Reverse-complement of 9-base sequence: bases read 8..3, 7..2, 6..1
                // RC base r[i] = complement(b[8-i])
                int r[9];
                for (int i = 0; i < 9; i++) r[i] = rc_base(b[8-i]);

                // RC frames 0, 1, 2 read from the RC sequence
                rev0[encode6(r[0],r[1],r[2],r[3],r[4],r[5])] += prob;
                rev1[encode6(r[1],r[2],r[3],r[4],r[5],r[6])] += prob;
                rev2[encode6(r[2],r[3],r[4],r[5],r[6],r[7])] += prob;
            }
        }
    }

    // Equal-weight 5-component mixture over all wrong-frame competitors
    std::array<float, 4096> null_dist{};
    for (int h = 0; h < 4096; h++)
        null_dist[h] = 0.2f * (plus1[h] + plus2[h] + rev0[h] + rev1[h] + rev2[h]);
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

// RC-strand null: for a hexamer h observed in an RC-oriented frame, the underlying
// coding sequence had hexamer rc(h). So the contrastive score is wf_llr[rc(h)].
// rc_complement(h): reverse the 6 bases and complement each (A↔T=0↔3, C↔G=1↔2).
inline int rc_complement_code(int h) {
    int b0 = (h >> 10) & 3, b1 = (h >> 8) & 3, b2 = (h >> 6) & 3;
    int b3 = (h >>  4) & 3, b4 = (h >> 2) & 3, b5 =  h       & 3;
    return encode6(3-b5, 3-b4, 3-b3, 3-b2, 3-b1, 3-b0);
}

inline const std::array<float, 4096>& get_rc_frame_llr() {
    static const std::array<float, 4096> rc_llr = []() {
        const auto& fwd = get_wrong_frame_llr();
        std::array<float, 4096> r{};
        for (int h = 0; h < 4096; h++)
            r[h] = fwd[rc_complement_code(h)];
        return r;
    }();
    return rc_llr;
}

}  // namespace dart
