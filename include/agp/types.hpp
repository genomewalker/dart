#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <cctype>

namespace agp {

// Basic sequence types
using Sequence = std::string;
using Position = uint32_t;
using Score = float;

// Nucleotide encoding
enum class Nucleotide : uint8_t {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4  // Unknown
};

// Convert char to nucleotide
inline Nucleotide char_to_nt(char c) {
    switch(c) {
        case 'A': case 'a': return Nucleotide::A;
        case 'C': case 'c': return Nucleotide::C;
        case 'G': case 'g': return Nucleotide::G;
        case 'T': case 't': return Nucleotide::T;
        default: return Nucleotide::N;
    }
}

// Convert nucleotide to char
inline char nt_to_char(Nucleotide nt) {
    switch(nt) {
        case Nucleotide::A: return 'A';
        case Nucleotide::C: return 'C';
        case Nucleotide::G: return 'G';
        case Nucleotide::T: return 'T';
        default: return 'N';
    }
}

// Quality score utilities
/**
 * Convert Phred quality score (ASCII character) to error probability
 * Phred+33 encoding: Q = -10 * log10(P)
 * P = 10^(-Q/10)
 */
inline float phred_to_error_prob(char qual_char) {
    int phred = static_cast<int>(qual_char) - 33;
    if (phred < 0) phred = 0;
    if (phred > 60) phred = 60;
    return std::pow(10.0f, -phred / 10.0f);
}

/**
 * Convert Phred quality to correct base probability
 */
inline float phred_to_correct_prob(char qual_char) {
    return 1.0f - phred_to_error_prob(qual_char);
}

/**
 * Quality-weighted nucleotide probabilities
 * Given observed nucleotide and quality, compute P(true_nt | obs_nt, quality)
 */
struct QualityWeightedNucleotide {
    float prob[4];  // Probability for each nucleotide A, C, G, T

    QualityWeightedNucleotide(Nucleotide observed, char quality) {
        float error_prob = phred_to_error_prob(quality);
        float correct_prob = 1.0f - error_prob;

        // Distribute error probability uniformly among other 3 bases
        float error_per_base = error_prob / 3.0f;

        for (int i = 0; i < 4; ++i) {
            if (i == static_cast<int>(observed)) {
                prob[i] = correct_prob;
            } else {
                prob[i] = error_per_base;
            }
        }
    }

    float get(Nucleotide nt) const {
        if (nt == Nucleotide::N) return 0.25f;  // Uniform for N
        return prob[static_cast<int>(nt)];
    }
};

// HMM states
enum class State : uint8_t {
    // Forward strand gene states (6-periodic for codon positions)
    GENE_FWD_0 = 0,
    GENE_FWD_1 = 1,
    GENE_FWD_2 = 2,
    GENE_FWD_3 = 3,
    GENE_FWD_4 = 4,
    GENE_FWD_5 = 5,

    // Reverse strand gene states
    GENE_REV_0 = 6,
    GENE_REV_1 = 7,
    GENE_REV_2 = 8,
    GENE_REV_3 = 9,
    GENE_REV_4 = 10,
    GENE_REV_5 = 11,

    // Start codon states (forward)
    START_FWD_0 = 12,
    START_FWD_1 = 13,
    START_FWD_2 = 14,

    // Start codon states (reverse)
    START_REV_0 = 15,
    START_REV_1 = 16,
    START_REV_2 = 17,

    // Stop codon states (forward)
    STOP_FWD_0 = 18,
    STOP_FWD_1 = 19,
    STOP_FWD_2 = 20,

    // Stop codon states (reverse)
    STOP_REV_0 = 21,
    STOP_REV_1 = 22,
    STOP_REV_2 = 23,

    // Intergenic
    INTERGENIC = 24,

    NUM_STATES = 25
};

constexpr size_t NUM_STATES = 25;
constexpr size_t NUM_NUCLEOTIDES = 4;

// Gene prediction result
// Cache-aligned for optimal performance with 1B+ reads
struct alignas(64) Gene {
    Position start;          // Start position (0-based)
    Position end;            // End position (0-based, exclusive)
    bool is_forward;         // Strand orientation
    bool is_fragment;        // True if fragment (no start/stop required)
    int frame = 0;           // Reading frame (0, 1, or 2)
    Score score;             // Prediction score
    Score ancient_prob;      // Probability of being ancient DNA
    Score coding_prob;       // Probability of being coding (for fragments)
    Score frame_score;       // Frame selection confidence score
    Score damage_score = 0;  // Observed damage magnitude (0-6 scale: T's at 5' + A's at 3')
    Score strand_conf = 0.5f; // Strand confidence (0.5 = uncertain, 1.0 = confident)
    std::string sequence;    // Predicted gene sequence (potentially damaged)
    std::string protein;     // Translated protein sequence (from damaged DNA)
    std::string corrected_sequence;   // Damage-corrected DNA sequence
    std::string corrected_protein;    // Damage-corrected protein sequence
    size_t dna_corrections = 0;       // Number of DNA corrections made
    size_t aa_corrections = 0;        // Number of amino acid changes from correction
    char correction_pattern = 'N';    // 'F'=forward strand, 'R'=reverse strand, 'N'=none
};

// Damage profile for a sequence
struct DamageProfile {
    std::vector<float> ct_prob_5prime;  // C->T probability from 5' end
    std::vector<float> ct_prob_3prime;  // C->T probability from 3' end
    std::vector<float> ga_prob_5prime;  // G->A probability from 5' end (rev complement)
    std::vector<float> ga_prob_3prime;  // G->A probability from 3' end (rev complement)
    float lambda_5prime = 0.1f;         // Decay parameter for 5' end
    float lambda_3prime = 0.1f;         // Decay parameter for 3' end
    float delta_max = 0.4f;             // Maximum damage rate
    float delta_background = 0.01f;     // Background damage rate

    // Calculate damage probability at position
    float get_ct_damage(Position pos, Position seq_len) const;
    float get_ga_damage(Position pos, Position seq_len) const;
};

// Damage transition matrix (4x4 for ACGT)
using DamageMatrix = std::array<std::array<float, 4>, 4>;

// Codon translation table (standard genetic code)
struct CodonTable {
    static inline char translate_codon(char c1, char c2, char c3) {
        c1 = std::toupper(c1);
        c2 = std::toupper(c2);
        c3 = std::toupper(c3);

        // TTx codons
        if (c1 == 'T' && c2 == 'T') {
            if (c3 == 'T' || c3 == 'C') return 'F';  // Phe
            return 'L';  // Leu
        }
        // TCx codons -> Ser
        if (c1 == 'T' && c2 == 'C') return 'S';
        // TAx codons
        if (c1 == 'T' && c2 == 'A') {
            if (c3 == 'T' || c3 == 'C') return 'Y';  // Tyr
            return '*';  // Stop (TAA, TAG)
        }
        // TGx codons
        if (c1 == 'T' && c2 == 'G') {
            if (c3 == 'T' || c3 == 'C') return 'C';  // Cys
            if (c3 == 'A') return '*';  // Stop (TGA)
            return 'W';  // Trp (TGG)
        }
        // CTx codons -> Leu
        if (c1 == 'C' && c2 == 'T') return 'L';
        // CCx codons -> Pro
        if (c1 == 'C' && c2 == 'C') return 'P';
        // CAx codons
        if (c1 == 'C' && c2 == 'A') {
            if (c3 == 'T' || c3 == 'C') return 'H';  // His
            return 'Q';  // Gln
        }
        // CGx codons -> Arg
        if (c1 == 'C' && c2 == 'G') return 'R';
        // ATx codons
        if (c1 == 'A' && c2 == 'T') {
            if (c3 == 'G') return 'M';  // Met
            return 'I';  // Ile
        }
        // ACx codons -> Thr
        if (c1 == 'A' && c2 == 'C') return 'T';
        // AAx codons
        if (c1 == 'A' && c2 == 'A') {
            if (c3 == 'T' || c3 == 'C') return 'N';  // Asn
            return 'K';  // Lys
        }
        // AGx codons
        if (c1 == 'A' && c2 == 'G') {
            if (c3 == 'T' || c3 == 'C') return 'S';  // Ser
            return 'R';  // Arg
        }
        // GTx codons -> Val
        if (c1 == 'G' && c2 == 'T') return 'V';
        // GCx codons -> Ala
        if (c1 == 'G' && c2 == 'C') return 'A';
        // GAx codons
        if (c1 == 'G' && c2 == 'A') {
            if (c3 == 'T' || c3 == 'C') return 'D';  // Asp
            return 'E';  // Glu
        }
        // GGx codons -> Gly
        if (c1 == 'G' && c2 == 'G') return 'G';

        return 'X';  // Unknown
    }

    static inline bool is_stop_codon(char c1, char c2, char c3) {
        c1 = std::toupper(c1);
        c2 = std::toupper(c2);
        c3 = std::toupper(c3);
        return (c1 == 'T' && c2 == 'A' && (c3 == 'A' || c3 == 'G')) ||
               (c1 == 'T' && c2 == 'G' && c3 == 'A');
    }

    static inline bool is_start_codon(char c1, char c2, char c3) {
        c1 = std::toupper(c1);
        c2 = std::toupper(c2);
        c3 = std::toupper(c3);
        // ATG is the standard start, also GTG and TTG in bacteria
        return (c1 == 'A' && c2 == 'T' && c3 == 'G') ||
               (c1 == 'G' && c2 == 'T' && c3 == 'G') ||
               (c1 == 'T' && c2 == 'T' && c3 == 'G');
    }
};

// Training data structure
struct TrainingData {
    std::vector<std::string> gene_sequences;
    std::vector<std::string> intergenic_sequences;
};

} // namespace agp
