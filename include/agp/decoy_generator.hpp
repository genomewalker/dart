// Decoy sequence generator for FDR calibration
// Generates frame-preserving decoy sequences for target-decoy FDR estimation

#pragma once

#include <string>
#include <vector>
#include <random>
#include <unordered_map>
#include <array>
#include <algorithm>
#include <cctype>

namespace agp {

// Synonymous codon tables for each amino acid
// Key: amino acid, Value: vector of codons encoding that amino acid
inline const std::unordered_map<char, std::vector<std::string>>& get_synonymous_codons() {
    static const std::unordered_map<char, std::vector<std::string>> table = {
        {'A', {"GCT", "GCC", "GCA", "GCG"}},
        {'R', {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}},
        {'N', {"AAT", "AAC"}},
        {'D', {"GAT", "GAC"}},
        {'C', {"TGT", "TGC"}},
        {'Q', {"CAA", "CAG"}},
        {'E', {"GAA", "GAG"}},
        {'G', {"GGT", "GGC", "GGA", "GGG"}},
        {'H', {"CAT", "CAC"}},
        {'I', {"ATT", "ATC", "ATA"}},
        {'L', {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}},
        {'K', {"AAA", "AAG"}},
        {'M', {"ATG"}},  // Only one codon
        {'F', {"TTT", "TTC"}},
        {'P', {"CCT", "CCC", "CCA", "CCG"}},
        {'S', {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}},
        {'T', {"ACT", "ACC", "ACA", "ACG"}},
        {'W', {"TGG"}},  // Only one codon
        {'Y', {"TAT", "TAC"}},
        {'V', {"GTT", "GTC", "GTA", "GTG"}},
        {'*', {"TAA", "TAG", "TGA"}}  // Stop codons
    };
    return table;
}

// Codon to amino acid lookup
inline char codon_to_aa(const std::string& codon) {
    static const std::unordered_map<std::string, char> table = {
        {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
        {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
        {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
        {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
        {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
        {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
        {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
        {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };
    auto it = table.find(codon);
    return (it != table.end()) ? it->second : 'X';
}

enum class DecoyStrategy {
    SYNONYMOUS_SHUFFLE,   // Shuffle synonymous codons (preserves AA sequence)
    CODON_REVERSE,        // Reverse codon order
    DINUCLEOTIDE_SHUFFLE, // Shuffle preserving dinucleotide frequencies
    RANDOM_FRAME          // Random sequence with matched length
};

class DecoyGenerator {
public:
    explicit DecoyGenerator(uint64_t seed = 42) : rng_(seed) {}

    // Generate a single decoy from target sequence
    std::string generate(const std::string& target, DecoyStrategy strategy) {
        switch (strategy) {
            case DecoyStrategy::SYNONYMOUS_SHUFFLE:
                return synonymous_shuffle(target);
            case DecoyStrategy::CODON_REVERSE:
                return codon_reverse(target);
            case DecoyStrategy::DINUCLEOTIDE_SHUFFLE:
                return dinucleotide_shuffle(target);
            case DecoyStrategy::RANDOM_FRAME:
                return random_frame(target);
            default:
                return synonymous_shuffle(target);
        }
    }

    // Generate multiple decoys per target
    std::vector<std::string> generate_batch(
        const std::vector<std::string>& targets,
        DecoyStrategy strategy,
        size_t decoys_per_target = 1
    ) {
        std::vector<std::string> decoys;
        decoys.reserve(targets.size() * decoys_per_target);

        for (const auto& target : targets) {
            for (size_t i = 0; i < decoys_per_target; ++i) {
                decoys.push_back(generate(target, strategy));
            }
        }
        return decoys;
    }

    // Set seed for reproducibility
    void set_seed(uint64_t seed) {
        rng_.seed(seed);
    }

private:
    std::mt19937_64 rng_;

    // Synonymous codon shuffle: preserve AA sequence, randomize codon choice
    std::string synonymous_shuffle(const std::string& seq) {
        if (seq.length() < 3) return seq;

        std::string result;
        result.reserve(seq.length());

        const auto& syn_table = get_synonymous_codons();

        // Process each codon
        for (size_t i = 0; i + 2 < seq.length(); i += 3) {
            std::string codon = seq.substr(i, 3);
            // Convert to uppercase
            for (char& c : codon) c = std::toupper(static_cast<unsigned char>(c));

            char aa = codon_to_aa(codon);

            if (aa == 'X') {
                // Unknown codon, keep as-is
                result += codon;
            } else {
                auto it = syn_table.find(aa);
                if (it != syn_table.end() && it->second.size() > 1) {
                    // Randomly select a synonymous codon
                    std::uniform_int_distribution<size_t> dist(0, it->second.size() - 1);
                    result += it->second[dist(rng_)];
                } else {
                    // Only one codon for this AA (M, W) or stop
                    result += codon;
                }
            }
        }

        // Handle remaining bases (if not multiple of 3)
        size_t remainder = seq.length() % 3;
        if (remainder > 0) {
            result += seq.substr(seq.length() - remainder);
        }

        return result;
    }

    // Codon reverse: reverse the order of codons (not nucleotides)
    std::string codon_reverse(const std::string& seq) {
        if (seq.length() < 3) return seq;

        std::vector<std::string> codons;
        for (size_t i = 0; i + 2 < seq.length(); i += 3) {
            codons.push_back(seq.substr(i, 3));
        }

        std::string result;
        result.reserve(seq.length());

        // Reverse codon order
        for (auto it = codons.rbegin(); it != codons.rend(); ++it) {
            result += *it;
        }

        // Handle remainder
        size_t remainder = seq.length() % 3;
        if (remainder > 0) {
            result += seq.substr(seq.length() - remainder);
        }

        return result;
    }

    // Dinucleotide shuffle: preserve dinucleotide frequencies
    std::string dinucleotide_shuffle(const std::string& seq) {
        if (seq.length() < 4) return seq;

        // Build dinucleotide transition graph
        std::unordered_map<char, std::vector<char>> transitions;
        for (size_t i = 0; i + 1 < seq.length(); ++i) {
            char from = std::toupper(static_cast<unsigned char>(seq[i]));
            char to = std::toupper(static_cast<unsigned char>(seq[i + 1]));
            if (from == 'A' || from == 'C' || from == 'G' || from == 'T') {
                transitions[from].push_back(to);
            }
        }

        // Shuffle each transition list
        for (auto& [base, nexts] : transitions) {
            std::shuffle(nexts.begin(), nexts.end(), rng_);
        }

        // Generate shuffled sequence
        std::string result;
        result.reserve(seq.length());
        result += std::toupper(static_cast<unsigned char>(seq[0]));

        std::unordered_map<char, size_t> indices = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};

        for (size_t i = 1; i < seq.length(); ++i) {
            char prev = result.back();
            if (transitions.count(prev) && indices[prev] < transitions[prev].size()) {
                result += transitions[prev][indices[prev]++];
            } else {
                // Fallback: random base
                static const char bases[] = "ACGT";
                std::uniform_int_distribution<int> dist(0, 3);
                result += bases[dist(rng_)];
            }
        }

        return result;
    }

    // Random frame: generate random sequence with same length
    std::string random_frame(const std::string& seq) {
        static const char bases[] = "ACGT";
        std::uniform_int_distribution<int> dist(0, 3);

        std::string result;
        result.reserve(seq.length());

        for (size_t i = 0; i < seq.length(); ++i) {
            result += bases[dist(rng_)];
        }

        return result;
    }
};

}  // namespace agp
