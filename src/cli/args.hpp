#ifndef AGP_CLI_ARGS_HPP
#define AGP_CLI_ARGS_HPP

#include <string>
#include <cstddef>
#include <cstdint>

namespace agp {
namespace cli {

// Library type for ancient DNA analysis
enum class LibraryType {
    UNKNOWN,          // Auto-detect from damage pattern
    DOUBLE_STRANDED,  // Standard double-stranded library
    SINGLE_STRANDED   // Single-stranded library prep
};

struct Options {
    // Input/Output
    std::string input_file;
    std::string output_file = "predictions.gff";
    std::string fasta_nt;
    std::string fasta_nt_corrected;  // Damage-corrected nucleotide output
    std::string fasta_aa;            // Protein output (observed, stops as *)
    std::string fasta_aa_masked;     // Search-optimized: terminal damage stops masked as X
    std::string summary_file;        // Summary statistics output file (JSON format)
    std::string damage_index;        // Binary damage index (.agd) for post-mapping annotation
    std::string bdamage_file;        // External metaDMG bdamage file for damage profile

    // Taxonomic domain
    std::string domain_name = "gtdb";  // Default to bacteria/archaea for ancient DNA

    // Basic options
    size_t min_length = 30;
    int num_threads = 0;
    bool verbose = false;
    bool use_damage = true;
    bool aggregate_damage = true;    // Two-pass damage detection
    bool damage_only = false;        // Run Pass 1 only, skip gene prediction
    LibraryType forced_library_type = LibraryType::UNKNOWN;

    // ORF enumeration parameters
    size_t orf_min_aa = 10;          // Minimum ORF length in amino acids
    bool adaptive_orf = false;       // Adaptive mode: output ORFs within score threshold
};

// Convert CLI LibraryType to string for display
inline const char* library_type_to_string(LibraryType lt) {
    switch (lt) {
        case LibraryType::DOUBLE_STRANDED: return "double-stranded";
        case LibraryType::SINGLE_STRANDED: return "single-stranded";
        default: return "auto";
    }
}

// Print version string to stdout
void print_version();

// Print usage/help to stdout
void print_usage(const char* program_name);

// Parse command-line arguments into Options struct
Options parse_args(int argc, char* argv[]);

}  // namespace cli
}  // namespace agp

#endif  // AGP_CLI_ARGS_HPP
