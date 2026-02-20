/**
 * Codon and amino acid frequency tables
 *
 * Encoding: T=0, C=1, A=2, G=3
 * Index = base1*16 + base2*4 + base3
 */

#include "dart/codon_tables.hpp"
#include <array>

namespace dart {

// Convert base to index: T=0, C=1, A=2, G=3
static inline int base_to_idx(char c) {
    switch (fast_upper(c)) {
        case 'T': return 0;
        case 'C': return 1;
        case 'A': return 2;
        case 'G': return 3;
        default: return -1;
    }
}

// Convert codon to array index (0-63)
static inline int codon_to_idx(char c1, char c2, char c3) {
    int i1 = base_to_idx(c1);
    int i2 = base_to_idx(c2);
    int i3 = base_to_idx(c3);
    if (i1 < 0 || i2 < 0 || i3 < 0) return -1;
    return i1 * 16 + i2 * 4 + i3;
}

// Static codon frequency array (E. coli RSCU)
// Order: TTT=0, TTC=1, TTA=2, TTG=3, TCT=4, ... GGG=63
static const std::array<float, 64> CODON_FREQ = []() {
    std::array<float, 64> arr{};
    // Fill with defaults
    arr.fill(0.1f);

    // TTx (0-3): Phe, Phe, Leu, Leu
    arr[0] = 0.57f;  // TTT
    arr[1] = 0.43f;  // TTC
    arr[2] = 0.11f;  // TTA
    arr[3] = 0.11f;  // TTG

    // TCx (4-7): Ser
    arr[4] = 0.17f;  // TCT
    arr[5] = 0.15f;  // TCC
    arr[6] = 0.14f;  // TCA
    arr[7] = 0.14f;  // TCG

    // TAx (8-11): Tyr, Tyr, Stop, Stop
    arr[8] = 0.59f;  // TAT
    arr[9] = 0.41f;  // TAC
    arr[10] = 0.61f; // TAA (stop)
    arr[11] = 0.09f; // TAG (stop)

    // TGx (12-15): Cys, Cys, Stop, Trp
    arr[12] = 0.46f; // TGT
    arr[13] = 0.54f; // TGC
    arr[14] = 0.30f; // TGA (stop)
    arr[15] = 1.00f; // TGG

    // CTx (16-19): Leu
    arr[16] = 0.10f; // CTT
    arr[17] = 0.10f; // CTC
    arr[18] = 0.03f; // CTA
    arr[19] = 0.55f; // CTG

    // CCx (20-23): Pro
    arr[20] = 0.18f; // CCT
    arr[21] = 0.13f; // CCC
    arr[22] = 0.20f; // CCA
    arr[23] = 0.49f; // CCG

    // CAx (24-27): His, His, Gln, Gln
    arr[24] = 0.57f; // CAT
    arr[25] = 0.43f; // CAC
    arr[26] = 0.34f; // CAA
    arr[27] = 0.66f; // CAG

    // CGx (28-31): Arg
    arr[28] = 0.36f; // CGT
    arr[29] = 0.36f; // CGC
    arr[30] = 0.07f; // CGA
    arr[31] = 0.11f; // CGG

    // ATx (32-35): Ile, Ile, Ile, Met
    arr[32] = 0.49f; // ATT
    arr[33] = 0.39f; // ATC
    arr[34] = 0.11f; // ATA
    arr[35] = 1.00f; // ATG

    // ACx (36-39): Thr
    arr[36] = 0.19f; // ACT
    arr[37] = 0.40f; // ACC
    arr[38] = 0.17f; // ACA
    arr[39] = 0.25f; // ACG

    // AAx (40-43): Asn, Asn, Lys, Lys
    arr[40] = 0.49f; // AAT
    arr[41] = 0.51f; // AAC
    arr[42] = 0.74f; // AAA
    arr[43] = 0.26f; // AAG

    // AGx (44-47): Ser, Ser, Arg, Arg
    arr[44] = 0.16f; // AGT
    arr[45] = 0.25f; // AGC
    arr[46] = 0.07f; // AGA
    arr[47] = 0.04f; // AGG

    // GTx (48-51): Val
    arr[48] = 0.28f; // GTT
    arr[49] = 0.20f; // GTC
    arr[50] = 0.17f; // GTA
    arr[51] = 0.35f; // GTG

    // GCx (52-55): Ala
    arr[52] = 0.18f; // GCT
    arr[53] = 0.26f; // GCC
    arr[54] = 0.23f; // GCA
    arr[55] = 0.33f; // GCG

    // GAx (56-59): Asp, Asp, Glu, Glu
    arr[56] = 0.63f; // GAT
    arr[57] = 0.37f; // GAC
    arr[58] = 0.68f; // GAA
    arr[59] = 0.32f; // GAG

    // GGx (60-63): Gly
    arr[60] = 0.35f; // GGT
    arr[61] = 0.37f; // GGC
    arr[62] = 0.13f; // GGA
    arr[63] = 0.15f; // GGG

    return arr;
}();

// Codon to amino acid translation table
static const std::array<char, 64> CODON_TO_AA = []() {
    std::array<char, 64> arr{};

    // TTx: Phe, Phe, Leu, Leu
    arr[0] = 'F'; arr[1] = 'F'; arr[2] = 'L'; arr[3] = 'L';
    // TCx: Ser
    arr[4] = 'S'; arr[5] = 'S'; arr[6] = 'S'; arr[7] = 'S';
    // TAx: Tyr, Tyr, Stop, Stop
    arr[8] = 'Y'; arr[9] = 'Y'; arr[10] = '*'; arr[11] = '*';
    // TGx: Cys, Cys, Stop, Trp
    arr[12] = 'C'; arr[13] = 'C'; arr[14] = '*'; arr[15] = 'W';

    // CTx: Leu
    arr[16] = 'L'; arr[17] = 'L'; arr[18] = 'L'; arr[19] = 'L';
    // CCx: Pro
    arr[20] = 'P'; arr[21] = 'P'; arr[22] = 'P'; arr[23] = 'P';
    // CAx: His, His, Gln, Gln
    arr[24] = 'H'; arr[25] = 'H'; arr[26] = 'Q'; arr[27] = 'Q';
    // CGx: Arg
    arr[28] = 'R'; arr[29] = 'R'; arr[30] = 'R'; arr[31] = 'R';

    // ATx: Ile, Ile, Ile, Met
    arr[32] = 'I'; arr[33] = 'I'; arr[34] = 'I'; arr[35] = 'M';
    // ACx: Thr
    arr[36] = 'T'; arr[37] = 'T'; arr[38] = 'T'; arr[39] = 'T';
    // AAx: Asn, Asn, Lys, Lys
    arr[40] = 'N'; arr[41] = 'N'; arr[42] = 'K'; arr[43] = 'K';
    // AGx: Ser, Ser, Arg, Arg
    arr[44] = 'S'; arr[45] = 'S'; arr[46] = 'R'; arr[47] = 'R';

    // GTx: Val
    arr[48] = 'V'; arr[49] = 'V'; arr[50] = 'V'; arr[51] = 'V';
    // GCx: Ala
    arr[52] = 'A'; arr[53] = 'A'; arr[54] = 'A'; arr[55] = 'A';
    // GAx: Asp, Asp, Glu, Glu
    arr[56] = 'D'; arr[57] = 'D'; arr[58] = 'E'; arr[59] = 'E';
    // GGx: Gly
    arr[60] = 'G'; arr[61] = 'G'; arr[62] = 'G'; arr[63] = 'G';

    return arr;
}();

// Amino acid frequencies (indexed by char - 'A')
static const std::array<float, 26> AA_FREQ = []() {
    std::array<float, 26> arr{};
    arr.fill(0.01f);  // Default for unknown

    arr['A' - 'A'] = 0.082569f;
    arr['C' - 'A'] = 0.013888f;
    arr['D' - 'A'] = 0.054628f;
    arr['E' - 'A'] = 0.067171f;
    arr['F' - 'A'] = 0.038694f;
    arr['G' - 'A'] = 0.070735f;
    arr['H' - 'A'] = 0.022788f;
    arr['I' - 'A'] = 0.059073f;
    arr['K' - 'A'] = 0.057985f;
    arr['L' - 'A'] = 0.096498f;
    arr['M' - 'A'] = 0.024120f;
    arr['N' - 'A'] = 0.040637f;
    arr['P' - 'A'] = 0.047481f;
    arr['Q' - 'A'] = 0.039324f;
    arr['R' - 'A'] = 0.055287f;
    arr['S' - 'A'] = 0.066602f;
    arr['T' - 'A'] = 0.053658f;
    arr['V' - 'A'] = 0.068565f;
    arr['W' - 'A'] = 0.011053f;
    arr['Y' - 'A'] = 0.029246f;

    return arr;
}();

// Fast codon frequency lookup
float get_codon_freq_fast(char c1, char c2, char c3) {
    int idx = codon_to_idx(c1, c2, c3);
    if (idx < 0) return 0.1f;
    return CODON_FREQ[idx];
}

// Fast amino acid translation
char translate_codon_fast(char c1, char c2, char c3) {
    int idx = codon_to_idx(c1, c2, c3);
    if (idx < 0) return 'X';
    return CODON_TO_AA[idx];
}

// Fast AA frequency lookup
float get_aa_freq_fast(char aa) {
    if (aa == '*') return 0.0001f;
    if (aa < 'A' || aa > 'Z') return 0.01f;
    return AA_FREQ[aa - 'A'];
}

} // namespace dart
