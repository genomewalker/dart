#include "args.hpp"
#include "dart/version.h"
#include <iostream>
#include <algorithm>
#include <string>

namespace dart {
namespace cli {

void print_version() {
    std::cout << "dart " << DART_VERSION << "\n";
}

void print_usage(const char* program_name) {
    std::cout << "DART v" << DART_VERSION << "\n\n";
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

        auto require_value = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) {
                throw ParseArgsExit(1, "Error: Missing value for " + flag);
            }
            return argv[++i];
        };

        auto parse_size = [&](const std::string& flag, const std::string& value) -> size_t {
            try {
                size_t idx = 0;
                size_t parsed = std::stoull(value, &idx);
                if (idx != value.size()) {
                    throw ParseArgsExit(1, "Error: Invalid integer for " + flag + ": " + value);
                }
                return parsed;
            } catch (const ParseArgsExit&) {
                throw;
            } catch (...) {
                throw ParseArgsExit(1, "Error: Invalid integer for " + flag + ": " + value);
            }
        };

        auto parse_int = [&](const std::string& flag, const std::string& value) -> int {
            try {
                size_t idx = 0;
                int parsed = std::stoi(value, &idx);
                if (idx != value.size()) {
                    throw ParseArgsExit(1, "Error: Invalid integer for " + flag + ": " + value);
                }
                return parsed;
            } catch (const ParseArgsExit&) {
                throw;
            } catch (...) {
                throw ParseArgsExit(1, "Error: Invalid integer for " + flag + ": " + value);
            }
        };

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            throw ParseArgsExit(0);
        } else if (arg == "-V" || arg == "--version") {
            print_version();
            throw ParseArgsExit(0);
        } else if (arg == "-i" || arg == "--input") {
            opts.input_file = require_value(arg);
        } else if (arg == "-o" || arg == "--output") {
            opts.output_file = require_value(arg);
        } else if (arg == "--fasta-nt") {
            opts.fasta_nt = require_value(arg);
        } else if (arg == "--fasta-nt-corrected") {
            opts.fasta_nt_corrected = require_value(arg);
        } else if (arg == "--fasta-aa") {
            opts.fasta_aa = require_value(arg);
        } else if (arg == "--fasta-aa-masked") {
            opts.fasta_aa_masked = require_value(arg);
        } else if (arg == "--summary") {
            opts.summary_file = require_value(arg);
        } else if (arg == "--damage-index") {
            opts.damage_index = require_value(arg);
        } else if (arg == "--no-damage") {
            opts.use_damage = false;
        } else if (arg == "--min-length") {
            opts.min_length = parse_size(arg, require_value(arg));
            if (opts.min_length < 1) {
                throw ParseArgsExit(1, "Error: --min-length must be >= 1");
            }
        } else if (arg == "-t" || arg == "--threads") {
            opts.num_threads = parse_int(arg, require_value(arg));
            if (opts.num_threads < 1) {
                throw ParseArgsExit(1, "Error: --threads must be >= 1");
            }
        } else if (arg == "-v" || arg == "--verbose") {
            opts.verbose = true;
        } else if (arg == "--no-aggregate") {
            opts.aggregate_damage = false;
        } else if (arg == "--damage-only") {
            opts.damage_only = true;
        } else if (arg == "--library-type") {
            std::string lib_type = require_value(arg);
            if (lib_type == "ds" || lib_type == "double-stranded") {
                opts.forced_library_type = LibraryType::DOUBLE_STRANDED;
            } else if (lib_type == "ss" || lib_type == "single-stranded") {
                opts.forced_library_type = LibraryType::SINGLE_STRANDED;
            } else if (lib_type == "auto") {
                opts.forced_library_type = LibraryType::UNKNOWN;
            } else {
                throw ParseArgsExit(1, "Error: Unknown library type '" + lib_type + "'");
            }
        } else if (arg == "--domain") {
            opts.domain_name = require_value(arg);
        } else if (arg == "--orf-min-aa") {
            opts.orf_min_aa = parse_size(arg, require_value(arg));
            if (opts.orf_min_aa < 1) {
                throw ParseArgsExit(1, "Error: --orf-min-aa must be >= 1");
            }
        } else if (arg == "--adaptive") {
            opts.adaptive_orf = true;
        } else {
            throw ParseArgsExit(1, "Error: Unknown option: " + arg);
        }
    }

    if (opts.input_file.empty()) {
        throw ParseArgsExit(1, "Error: No input file specified");
    }

    return opts;
}

}  // namespace cli
}  // namespace dart
