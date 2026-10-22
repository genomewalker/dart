#ifndef DART_CLI_ARGS_HPP
#define DART_CLI_ARGS_HPP

#include <string>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace dart {
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
    std::string fasta_aa_corrected;  // Protein output (Bayesian-corrected ancestor AAs at damage stops)
    std::string fasta_aa_masked;     // Search-optimized: terminal damage stops masked as X
    int decoy_every = 0;             // --decoy-every N with --fasta-aa: reversed proteins of every N-th read -> <aa>.decoy.<ext>; 0 (default) = off
    std::string summary_file;        // Summary statistics output file (JSON format)
    std::string damage_index;        // Binary damage index (.agd) for post-mapping annotation

    // Library damage profile sidecar (single-engine: fqdup-fit DamageSummary).
    // When set, predict applies it as the authoritative damage gate instead of
    // re-deriving the profile from the reads it annotates. Required by default:
    // self-profiling is disabled unless allow_self_profile is set (standalone use).
    std::string library_profile;     // path to a taph DamageSummary sidecar (.tdmg)
    bool allow_self_profile = false; // opt-in escape hatch for the self-profile fallback

    // Self-profile input decoupling (§8.1): damage is profiled on this FASTQ, ORF
    // prediction always runs on -i. Empty = profile -i. The pipeline profiles the
    // MERGED pre-dedup reads (high count -> DETECTED) while predicting on the dedup set.
    std::string profile_input;       // FASTQ to profile damage on; empty = input_file
    size_t      profile_sample = 0;  // reads sampled for self-profile adapter pre-scan; 0 = all

    // Taxonomic domain
    std::string domain_name = "meta";  // Default to mixed-community hexamer tables for ancient metagenomes

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
    bool adaptive_orf = true;        // Adaptive mode: output ORFs within score threshold
};

// Exception used by parse_args to signal an early/controlled CLI exit.
class ParseArgsExit : public std::runtime_error {
public:
    ParseArgsExit(int exit_code, const std::string& message = "")
        : std::runtime_error(message), exit_code_(exit_code) {}

    int exit_code() const noexcept { return exit_code_; }

private:
    int exit_code_ = 1;
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

// Decoy peptide file next to --fasta-aa: ".decoy" inserted before the extension so
// compression detection still works (proteins.faa.gz -> proteins.decoy.faa.gz).
std::string decoy_fasta_path(const std::string& fasta_aa);

}  // namespace cli
}  // namespace dart

#endif // DART_CLI_ARGS_HPP
