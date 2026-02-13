#include "args.hpp"
#include "agp/version.h"
#include <iostream>
#include <algorithm>

namespace agp {
namespace cli {

void print_version() {
    std::cout << "agp " << AGP_VERSION << "\n";
}

void print_usage(const char* program_name) {
    std::cout << "Ancient Gene Predictor v" << AGP_VERSION << "\n\n";
    std::cout << "Usage: " << program_name << " -i <input> -o <output> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -i, --input <file>       Input FASTA/FASTQ file (or .gz)\n";
    std::cout << "  -o, --output <file>      Output GFF3 file (default: predictions.gff)\n";
    std::cout << "  --fasta-nt <file>        Output nucleotide FASTA\n";
    std::cout << "  --fasta-nt-corrected <file> Output damage-corrected nucleotide FASTA\n";
    std::cout << "  --fasta-aa <file>        Output amino acid FASTA (observed, stops as *)\n";
    std::cout << "  --fasta-aa-masked <file> Search-optimized: terminal damage stops as X\n";
    std::cout << "  --summary <file>         Output sample statistics (JSON format)\n";
    std::cout << "  --damage-index <file>    Output binary damage index (.agd) for damage-annotate\n";
    std::cout << "  --min-length <int>       Minimum sequence length (default: 30)\n";
    std::cout << "  -t, --threads <int>      Number of threads (default: auto)\n";
    std::cout << "  --no-damage              Disable damage detection\n";
    std::cout << "  --no-aggregate           Disable two-pass damage aggregation\n";
    std::cout << "  --damage-only            Run damage detection only, skip gene prediction\n";
    std::cout << "  --domain <name>          Taxonomic domain (default: gtdb)\n";
    std::cout << "                           Options: gtdb, fungi, plant, protozoa, invertebrate, viral\n";
    std::cout << "  --library-type <type>    Force library type: ds, ss, or auto (default)\n";
    std::cout << "\nORF Parameters:\n";
    std::cout << "  --orf-min-aa <int>       Minimum ORF length in amino acids (default: 10)\n";
    std::cout << "  --adaptive               Adaptive ORF selection (score-based threshold)\n";
    std::cout << "\n";
    std::cout << "  -v, --verbose            Verbose output\n";
    std::cout << "  -V, --version            Show version and exit\n";
    std::cout << "  -h, --help               Show this help message\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  " << program_name << " -i reads.fq.gz -o predictions.gff --fasta-aa proteins.faa\n";
    std::cout << "  " << program_name << " -i reads.fq.gz -o predictions.gff --no-damage\n";
}

Options parse_args(int argc, char* argv[]) {
    Options opts;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            exit(0);
        } else if (arg == "-V" || arg == "--version") {
            print_version();
            exit(0);
        } else if (arg == "-i" || arg == "--input") {
            if (i + 1 < argc) opts.input_file = argv[++i];
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) opts.output_file = argv[++i];
        } else if (arg == "--fasta-nt") {
            if (i + 1 < argc) opts.fasta_nt = argv[++i];
        } else if (arg == "--fasta-nt-corrected") {
            if (i + 1 < argc) opts.fasta_nt_corrected = argv[++i];
        } else if (arg == "--fasta-aa") {
            if (i + 1 < argc) opts.fasta_aa = argv[++i];
        } else if (arg == "--fasta-aa-masked") {
            if (i + 1 < argc) opts.fasta_aa_masked = argv[++i];
        } else if (arg == "--summary") {
            if (i + 1 < argc) opts.summary_file = argv[++i];
        } else if (arg == "--damage-index") {
            if (i + 1 < argc) opts.damage_index = argv[++i];
        } else if (arg == "--no-damage") {
            opts.use_damage = false;
        } else if (arg == "--min-length") {
            if (i + 1 < argc) opts.min_length = std::stoul(argv[++i]);
        } else if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) opts.num_threads = std::stoi(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            opts.verbose = true;
        } else if (arg == "--no-aggregate") {
            opts.aggregate_damage = false;
        } else if (arg == "--damage-only") {
            opts.damage_only = true;
        } else if (arg == "--library-type") {
            if (i + 1 < argc) {
                std::string lib_type = argv[++i];
                if (lib_type == "ds" || lib_type == "double-stranded") {
                    opts.forced_library_type = LibraryType::DOUBLE_STRANDED;
                } else if (lib_type == "ss" || lib_type == "single-stranded") {
                    opts.forced_library_type = LibraryType::SINGLE_STRANDED;
                } else if (lib_type == "auto") {
                    opts.forced_library_type = LibraryType::UNKNOWN;
                } else {
                    std::cerr << "Error: Unknown library type '" << lib_type << "'\n";
                    exit(1);
                }
            }
        } else if (arg == "--domain") {
            if (i + 1 < argc) opts.domain_name = argv[++i];
        } else if (arg == "--orf-min-aa") {
            if (i + 1 < argc) {
                opts.orf_min_aa = std::stoul(argv[++i]);
                if (opts.orf_min_aa < 1) {
                    std::cerr << "Error: --orf-min-aa must be >= 1\n";
                    exit(1);
                }
            }
        } else if (arg == "--adaptive") {
            opts.adaptive_orf = true;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            exit(1);
        }
    }

    if (opts.input_file.empty()) {
        std::cerr << "Error: No input file specified\n\n";
        print_usage(argv[0]);
        exit(1);
    }

    return opts;
}

}  // namespace cli
}  // namespace agp
