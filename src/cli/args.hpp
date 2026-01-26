#ifndef AGP_CLI_ARGS_HPP
#define AGP_CLI_ARGS_HPP

#include <string>
#include <cstddef>
#include <cstdint>

namespace agp {
namespace cli {

// Library type for ancient DNA analysis
// Matches SampleDamageProfile::LibraryType but defined here to avoid circular deps
enum class LibraryType {
    UNKNOWN,          // Auto-detect from damage pattern
    DOUBLE_STRANDED,  // Standard double-stranded library
    SINGLE_STRANDED   // Single-stranded library prep
};

struct Options {
    std::string input_file;
    std::string output_file = "predictions.gff";
    std::string fasta_nt;
    std::string fasta_nt_corrected;  // Damage-corrected nucleotide output
    std::string fasta_aa;
    std::string fasta_corrected;  // Damage-corrected protein output for benchmarking
    std::string summary_file;     // Summary statistics output file (JSON format)
    std::string domain_name = "gtdb";  // Default to bacteria/archaea for ancient DNA
    std::string bdamage_file;     // External metaDMG bdamage file for damage profile
    bool use_damage = true;
    size_t min_length = 30;
    float min_coding_prob = 0.5f;
    int num_threads = 0;
    bool verbose = false;
    bool dual_strand = true;       // Default ON for short aDNA
    bool best_strand_only = false; // Default OFF: adaptive-orf mode decides output
    bool aggregate_damage = true;  // Default ON for better damage detection
    bool all_frames = false;       // Output all 6 reading frames per read
    bool metagenome_mode = false;  // Use weighted ensemble scoring across all domains
    bool use_bayesian = false;     // Use Bayesian frame selection (hexamer LLR posteriors)
    bool iterative_damage = false; // Iterative damage refinement for low-signal samples
    float iterative_threshold = 0.5f;  // Min coding_prob for damage re-estimation
    float strand_weight = 1.0f;    // Weight for damage_signal in strand selection (0 = coding only, 1 = full bias)
    LibraryType forced_library_type = LibraryType::UNKNOWN;  // UNKNOWN = auto-detect
    float orf_confidence_threshold = 0.3f;   // Adaptive ORF output: score gap threshold (default ON)
                                             // When >= 0: output 2 ORFs if gap < threshold, 1 ORF if gap >= threshold
                                             // Default 0.3 achieves 87% accuracy with 1.47 ORFs/read average
    float min_confidence = 0.7f;             // Minimum calibrated confidence to output prediction (0 = output all)
                                             // Uses isotonic calibration: raw_score -> P(correct)
                                             // Recommended thresholds: 0.80 (high), 0.90 (strict)
    bool compute_fdr = false;                // Enable decoy-based FDR estimation
    size_t fdr_decoys = 1;                   // Number of decoys per target sequence
    float fdr_threshold = -1.0f;             // FDR threshold for filtering (-1 = no filtering, report only)
    uint64_t fdr_seed = 42;                  // Random seed for decoy generation (reproducibility)
    std::string fdr_strategy = "dinucleotide"; // Decoy strategy: synonymous, dinucleotide, random

    // Damage percentile filtering (contaminant removal)
    float min_damage_pctile = 0.0f;  // Minimum damage percentile to output (0 = no filtering)
                                      // Reads below this percentile are potential contaminants
                                      // Recommended: 10-20% for contaminated samples
    bool output_damage_pctile = true; // Include damage_pctile in output (default: yes)
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
// Exits with code 0 for --help/--version
// Exits with code 1 for errors (missing input, unknown options)
Options parse_args(int argc, char* argv[]);

}  // namespace cli
}  // namespace agp

#endif  // AGP_CLI_ARGS_HPP
