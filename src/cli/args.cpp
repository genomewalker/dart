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
    std::cout << "  --fasta-nt <file>        Output nucleotide FASTA (original)\n";
    std::cout << "  --fasta-nt-corrected <file> Output damage-corrected nucleotide FASTA\n";
    std::cout << "  --fasta-aa <file>        Output amino acid FASTA (original)\n";
    std::cout << "  --fasta-corrected <file> Output damage-corrected amino acid FASTA\n";
    std::cout << "  --summary <file>         Output sample statistics (JSON format)\n";
    std::cout << "  --min-length <int>       Minimum sequence length (default: 30)\n";
    std::cout << "  --min-coding-prob <f>    Minimum coding probability (default: 0.5)\n";
    std::cout << "  --threads <int>          Number of threads (default: auto)\n";
    std::cout << "  --no-damage              Disable damage detection\n";
    std::cout << "  --single-strand          Score only forward strand\n";
    std::cout << "  --both-strands           Output predictions for both strands\n";
    std::cout << "  --no-aggregate           Disable two-pass damage aggregation\n";
    std::cout << "  --all-frames             Output all 6 reading frames per read\n";
    std::cout << "  --bayesian               Use Bayesian frame selection (hexamer posterior)\n";
    std::cout << "                           Emits 0-2 frames based on posterior confidence\n";
    std::cout << "  --domain <name>          Taxonomic domain for hexamer scoring (default: gtdb)\n";
    std::cout << "                           Options: gtdb, fungi, plant, protozoa, invertebrate,\n";
    std::cout << "                                    viral, mammal, vertebrate, meta\n";
    std::cout << "                           'meta' enables ensemble scoring + iterative damage\n";
    std::cout << "  --iterative-damage       Enable iterative damage refinement (4-pass mode)\n";
    std::cout << "                           Re-estimates damage using only high-scoring genes\n";
    std::cout << "                           Recommended for metagenomes with mixed damage signals\n";
    std::cout << "  --iterative-threshold <f> Min coding probability for damage re-estimation (default: 0.5)\n";
    std::cout << "  --strand-weight <f>      Weight for damage_signal in best-strand selection (default: 1.0)\n";
    std::cout << "                           0.0 = pure coding score, 1.0 = full damage bias\n";
    std::cout << "  --library-type <type>    Force library type: ds (double-stranded), ss (single-stranded),\n";
    std::cout << "                           or auto (default, auto-detect from damage pattern)\n";
    std::cout << "  --bdamage <file>         Import damage profile from metaDMG bdamage file\n";
    std::cout << "                           Skips internal damage detection, uses metaDMG estimates\n";
    std::cout << "                           Generate with: metaDMG-cpp getdamage -r 0 <bam>\n";
    std::cout << "  --no-adaptive-orf        Disable adaptive ORF output (always output best strand)\n";
    std::cout << "                           Default: adaptive ON (1 ORF when confident, 2 when uncertain)\n";
    std::cout << "  --orf-confidence <f>     Score gap threshold for adaptive ORF output\n";
    std::cout << "                           Lower = more single-ORF output, higher = more dual-ORF\n";
    std::cout << "                           0 = always 1 ORF (strict), large = always 2 ORFs\n";
    std::cout << "                           Default: 0.3 (outputs 2 ORFs when scores differ by <0.3)\n";
    std::cout << "  --min-confidence <f>     Minimum calibrated confidence P(correct) for output\n";
    std::cout << "                           0.70 = default (balanced for ancient DNA)\n";
    std::cout << "                           0.80 = ~75% coverage, ~81% accuracy\n";
    std::cout << "                           0.85 = ~30% coverage, ~85% accuracy\n";
    std::cout << "                           0.90 = ~13% coverage, ~87% accuracy\n";
    std::cout << "  --confidence <mode>      Preset confidence filtering mode:\n";
    std::cout << "                           standard = no filtering\n";
    std::cout << "                           high = --min-confidence 0.80\n";
    std::cout << "                           strict = --min-confidence 0.85\n";
    std::cout << "  --compute-fdr            Enable decoy-based FDR estimation\n";
    std::cout << "                           Generates synonymous codon shuffles as decoys\n";
    std::cout << "                           Reports FDR at various thresholds in --summary\n";
    std::cout << "  --fdr-decoys <int>       Number of decoys per target (default: 1)\n";
    std::cout << "  --fdr-threshold <f>      Filter predictions by FDR (e.g., 0.05 for 5% FDR)\n";
    std::cout << "                           Requires --compute-fdr; outputs only predictions\n";
    std::cout << "                           passing the empirically-determined threshold\n";
    std::cout << "  --fdr-seed <int>         Random seed for decoy generation (default: 42)\n";
    std::cout << "  --fdr-strategy <name>    Decoy generation strategy (default: dinucleotide)\n";
    std::cout << "                           synonymous = shuffle synonymous codons (preserves AA)\n";
    std::cout << "                           dinucleotide = preserve dinucleotide frequencies\n";
    std::cout << "                           random = random sequence (same length)\n";
    std::cout << "  --min-damage-pctile <f>  Minimum damage percentile to output (0 = no filtering)\n";
    std::cout << "                           Reads below this percentile are potential contaminants\n";
    std::cout << "                           Recommended: 10-20 for contaminated samples\n";
    std::cout << "  --no-damage-pctile       Do not include damage_pctile in output\n";
    std::cout << "  -v, --verbose            Verbose output\n";
    std::cout << "  -V, --version            Show version and exit\n";
    std::cout << "  -h, --help               Show this help message\n";
    std::cout << "\n";
    std::cout << "Defaults:\n";
    std::cout << "  - Adaptive ORF output: 1 prediction when confident, 2 when uncertain\n";
    std::cout << "  - Confidence filtering: min_confidence=0.7, min_coding_prob=0.5\n";
    std::cout << "  - Aggregate damage detection (two-pass for sample-level statistics)\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  # Standard usage (best strand + aggregate damage)\n";
    std::cout << "  " << program_name << " -i reads.fastq -o predictions.gff\n\n";
    std::cout << "  # Output both strands per read\n";
    std::cout << "  " << program_name << " -i reads.fastq -o predictions.gff --both-strands\n\n";
    std::cout << "  # Undamaged DNA (skip damage detection)\n";
    std::cout << "  " << program_name << " -i reads.fastq -o predictions.gff --no-damage\n\n";
    std::cout << "  # Iterative damage refinement for mixed metagenomes\n";
    std::cout << "  " << program_name << " -i reads.fastq -o predictions.gff --iterative-damage\n\n";
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
        } else if (arg == "--fasta-corrected") {
            if (i + 1 < argc) opts.fasta_corrected = argv[++i];
        } else if (arg == "--summary") {
            if (i + 1 < argc) opts.summary_file = argv[++i];
        } else if (arg == "--no-damage") {
            opts.use_damage = false;
        } else if (arg == "--min-coding-prob") {
            if (i + 1 < argc) opts.min_coding_prob = std::stof(argv[++i]);
        } else if (arg == "--min-length") {
            if (i + 1 < argc) opts.min_length = std::stoul(argv[++i]);
        } else if (arg == "--threads") {
            if (i + 1 < argc) opts.num_threads = std::stoi(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            opts.verbose = true;
        } else if (arg == "--single-strand") {
            opts.dual_strand = false;
        } else if (arg == "--both-strands") {
            opts.best_strand_only = false;
        } else if (arg == "--no-aggregate") {
            opts.aggregate_damage = false;
        } else if (arg == "--aggregate-damage") {
            opts.aggregate_damage = true;
        } else if (arg == "--all-frames") {
            opts.all_frames = true;
            opts.best_strand_only = false;  // All-frames implies outputting all strands
        } else if (arg == "--bayesian") {
            opts.use_bayesian = true;
        } else if (arg == "--iterative-damage") {
            opts.iterative_damage = true;
        } else if (arg == "--iterative-threshold") {
            if (i + 1 < argc) opts.iterative_threshold = std::stof(argv[++i]);
        } else if (arg == "--strand-weight") {
            if (i + 1 < argc) opts.strand_weight = std::clamp(std::stof(argv[++i]), 0.0f, 1.0f);
        } else if (arg == "--library-type") {
            if (i + 1 < argc) {
                std::string lib_type = argv[++i];
                if (lib_type == "ds" || lib_type == "double" || lib_type == "double-stranded") {
                    opts.forced_library_type = LibraryType::DOUBLE_STRANDED;
                } else if (lib_type == "ss" || lib_type == "single" || lib_type == "single-stranded") {
                    opts.forced_library_type = LibraryType::SINGLE_STRANDED;
                } else if (lib_type == "auto") {
                    opts.forced_library_type = LibraryType::UNKNOWN;
                } else {
                    std::cerr << "Error: Unknown library type '" << lib_type << "'\n";
                    std::cerr << "Valid options: ds (double-stranded), ss (single-stranded), auto\n";
                    exit(1);
                }
            }
        } else if (arg == "--domain") {
            if (i + 1 < argc) {
                opts.domain_name = argv[++i];
                // "meta" enables weighted ensemble scoring across all domains
                // and automatically enables iterative damage refinement
                if (opts.domain_name == "meta" || opts.domain_name == "metagenome" || opts.domain_name == "all") {
                    opts.metagenome_mode = true;
                    opts.iterative_damage = true;  // Auto-enable for metagenomes
                }
            }
        } else if (arg == "--bdamage") {
            if (i + 1 < argc) opts.bdamage_file = argv[++i];
        } else if (arg == "--orf-confidence" || arg == "--orf-confidence-threshold") {
            if (i + 1 < argc) {
                opts.orf_confidence_threshold = std::stof(argv[++i]);
                // When adaptive mode is enabled, override best_strand_only
                // The actual decision is made per-read based on score gap
                if (opts.orf_confidence_threshold >= 0.0f) {
                    opts.best_strand_only = false;  // Allow both strands to be considered
                }
            }
        } else if (arg == "--no-adaptive-orf") {
            // Disable adaptive ORF output - always output only best strand
            opts.orf_confidence_threshold = -1.0f;
            opts.best_strand_only = true;
        } else if (arg == "--adaptive-orf") {
            // Explicit enable (already default, but keep for backwards compat)
            opts.orf_confidence_threshold = 0.3f;
            opts.best_strand_only = false;
        } else if (arg == "--min-confidence") {
            if (i + 1 < argc) {
                opts.min_confidence = std::stof(argv[++i]);
                if (opts.min_confidence < 0.0f || opts.min_confidence > 1.0f) {
                    std::cerr << "Error: --min-confidence must be between 0 and 1\n";
                    exit(1);
                }
            }
        } else if (arg == "--confidence") {
            if (i + 1 < argc) {
                std::string mode = argv[++i];
                if (mode == "standard" || mode == "none") {
                    opts.min_confidence = 0.0f;
                } else if (mode == "high") {
                    opts.min_confidence = 0.80f;
                } else if (mode == "strict") {
                    opts.min_confidence = 0.85f;
                } else {
                    std::cerr << "Error: Unknown confidence mode '" << mode << "'\n";
                    std::cerr << "Valid options: standard, high, strict\n";
                    exit(1);
                }
            }
        } else if (arg == "--compute-fdr") {
            opts.compute_fdr = true;
        } else if (arg == "--fdr-decoys") {
            if (i + 1 < argc) {
                opts.fdr_decoys = std::stoul(argv[++i]);
                if (opts.fdr_decoys == 0) {
                    std::cerr << "Error: --fdr-decoys must be >= 1\n";
                    exit(1);
                }
            }
        } else if (arg == "--fdr-threshold") {
            if (i + 1 < argc) {
                opts.fdr_threshold = std::stof(argv[++i]);
                if (opts.fdr_threshold < 0.0f || opts.fdr_threshold > 1.0f) {
                    std::cerr << "Error: --fdr-threshold must be between 0 and 1\n";
                    exit(1);
                }
                opts.compute_fdr = true;  // Implicitly enable FDR computation
            }
        } else if (arg == "--fdr-seed") {
            if (i + 1 < argc) {
                opts.fdr_seed = std::stoull(argv[++i]);
            }
        } else if (arg == "--fdr-strategy") {
            if (i + 1 < argc) {
                opts.fdr_strategy = argv[++i];
                if (opts.fdr_strategy != "synonymous" &&
                    opts.fdr_strategy != "dinucleotide" &&
                    opts.fdr_strategy != "random") {
                    std::cerr << "Error: Unknown FDR strategy '" << opts.fdr_strategy << "'\n";
                    std::cerr << "Valid options: synonymous, dinucleotide, random\n";
                    exit(1);
                }
            }
        } else if (arg == "--min-damage-pctile") {
            if (i + 1 < argc) opts.min_damage_pctile = std::stof(argv[++i]);
        } else if (arg == "--no-damage-pctile") {
            opts.output_damage_pctile = false;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            std::cerr << "All arguments must be named (e.g., -i input.fa -o output.gff)\n";
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
