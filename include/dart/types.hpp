#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <cctype>

namespace dart {

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
 * Normalize raw coding score to [0,1] range using calibrated sigmoid.
 * Raw scores typically range [0.75, 1.08] with mean ~0.86.
 * Sigmoid parameters calibrated so:
 *   - score 0.75 -> ~0.11 (low confidence)
 *   - score 0.85 -> ~0.50 (median)
 *   - score 0.95 -> ~0.88 (high confidence)
 *   - score 1.05 -> ~0.99 (very high confidence)
 *
 * This normalization is applied at output only (not during frame selection)
 * to preserve ranking while producing intuitive [0,1] scores for FDR.
 */
inline float normalize_coding_score(float raw_score) {
    constexpr float center = 0.85f;  // Median of typical score distribution
    constexpr float k = 22.0f;       // Steepness (calibrated for 0.2 score spread)
    return 1.0f / (1.0f + std::exp(-k * (raw_score - center)));
}

/**
 * Calculate confidence score from protein length and optional signals.
 *
 * Calibrated against benchmark data showing cross-species error rates:
 *   - <15aa:  18.9% errors → ~0.35 confidence
 *   - 15-19aa: 12.3% errors → ~0.45 confidence
 *   - 20-29aa: 8.6% errors → ~0.60 confidence
 *   - 30-49aa: 5.9% errors → ~0.75 confidence
 *   - 50+aa:  ~3% errors → ~0.90 confidence
 *
 * Returns value in [0,1] range suitable for user-facing filtering thresholds.
 */
inline float calculate_confidence(size_t protein_length, float coding_prob = 0.0f) {
    // Length-based confidence: sigmoid centered at 25aa
    // Maps: 10aa→0.27, 15aa→0.38, 20aa→0.50, 30aa→0.73, 40aa→0.88
    constexpr float length_center = 25.0f;
    constexpr float length_k = 0.12f;  // Steepness
    float length_conf = 1.0f / (1.0f + std::exp(-length_k * (protein_length - length_center)));

    // If coding probability is available, incorporate it
    if (coding_prob > 0.0f) {
        // Geometric mean of length and coding confidence
        return std::sqrt(length_conf * coding_prob);
    }

    return length_conf;
}

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
    Score damage_signal;     // Terminal damage signal strength (calibrated 0-1, ~P(terminal damage))
    Score coding_prob;       // Probability of being coding (for fragments)
    Score frame_score;       // Frame selection confidence score
    Score damage_score = 0;  // Observed damage magnitude (0-100% scale: terminal pattern match)
    Score p_read_damaged = 0; // P(read has damage) - codon-aware Bayesian inference (0-1)
    Score damage_log_lr = 0;  // Log-likelihood ratio for damage (for aggregation across reads)
    Score damage_info = 0;    // Informativeness: sum of pC*D at terminals (weight for aggregation)
    Score strand_conf = 0.5f; // Strand confidence (0.5 = uncertain, 1.0 = confident)
    std::string sequence;    // Predicted gene sequence (potentially damaged)
    std::string protein;     // Translated protein sequence (from damaged DNA, stops as *)
    std::string search_protein;       // Search-optimized: terminal damage stops as X (for MMseqs2)
    std::string corrected_sequence;   // Damage-corrected DNA sequence
    std::string corrected_protein;    // Observed protein with X markers (for debugging)
    size_t dna_corrections = 0;       // Number of DNA corrections made
    size_t aa_corrections = 0;        // Number of amino acid changes from correction
    size_t stop_fixes = 0;            // Number of stop codons corrected (damage-induced)
    size_t stop_restorations = 0;     // Number of stop codons restored to sense codons
    Score corrected_coding_prob = 0;  // Coding probability after correction (0 = not computed)
    char correction_pattern = 'N';    // 'F'=forward strand, 'R'=reverse strand, 'N'=none
    Score p_correct = 0.0f;           // Calibrated probability of correct prediction (0-1) [DEPRECATED]
    Score damage_pctile = 0.0f;       // Percentile rank of damage_signal within sample (0-100)
    size_t internal_stops = 0;        // Number of internal stop codons in predicted protein
    Score frame_pmax = 0.0f;          // Best frame posterior (max of 6-way)
    Score frame_margin = 0.0f;        // Frame confidence margin (p_max - p_second)
    Score confidence = 0.0f;          // Combined confidence score (0-1) for filtering

    // ORF enumeration mode fields
    size_t orf_rank = 0;              // Rank in sorted length order (0 = longest)
    size_t passed_stops = 0;          // Number of damage-probable stops passed over (marked as X)
    size_t x_count = 0;               // Number of X amino acids (ambiguous positions)
    Score orf_length_score = 0.0f;    // Length-based score for ORF ranking
    Score per_read_damage_evidence = 0.0f; // Sum of p_damage_stop over convertible terminal stops

    // Frameshift detection fields
    bool has_frameshift = false;      // True if frameshift was detected in this read
    size_t frameshift_region = 0;     // Region index (0 = first region)
    size_t n_frameshift_regions = 1;  // Total number of regions from this read
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
        c1 = std::toupper(static_cast<unsigned char>(c1));
        c2 = std::toupper(static_cast<unsigned char>(c2));
        c3 = std::toupper(static_cast<unsigned char>(c3));

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
        c1 = std::toupper(static_cast<unsigned char>(c1));
        c2 = std::toupper(static_cast<unsigned char>(c2));
        c3 = std::toupper(static_cast<unsigned char>(c3));
        return (c1 == 'T' && c2 == 'A' && (c3 == 'A' || c3 == 'G')) ||
               (c1 == 'T' && c2 == 'G' && c3 == 'A');
    }

    static inline bool is_start_codon(char c1, char c2, char c3) {
        c1 = std::toupper(static_cast<unsigned char>(c1));
        c2 = std::toupper(static_cast<unsigned char>(c2));
        c3 = std::toupper(static_cast<unsigned char>(c3));
        // ATG is the standard start, also GTG and TTG in bacteria
        return (c1 == 'A' && c2 == 'T' && c3 == 'G') ||
               (c1 == 'G' && c2 == 'T' && c3 == 'G') ||
               (c1 == 'T' && c2 == 'T' && c3 == 'G');
    }
};

} // namespace dart
