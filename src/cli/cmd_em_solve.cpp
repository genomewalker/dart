// agp em-solve: EM on binary index
//
// Runs SQUAREM-accelerated EM on pre-indexed alignment data.
// All data is mmap'd - no parsing overhead.

#include "subcommand.hpp"
#include "agp/version.h"
#include "agp/em_index.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <chrono>
#include <iomanip>

namespace agp {
namespace cli {

int cmd_em_solve(int argc, char* argv[]) {
    std::string index_file;
    std::string output_file;
    std::string summary_file;
    double lambda_b = 3.0;
    uint32_t max_iters = 100;
    double tol = 1e-4;
    bool use_damage = true;
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "--index") == 0 || strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
            index_file = argv[++i];
        } else if ((strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--summary") == 0 && i + 1 < argc) {
            summary_file = argv[++i];
        } else if (strcmp(argv[i], "--lambda") == 0 && i + 1 < argc) {
            lambda_b = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--max-iters") == 0 && i + 1 < argc) {
            max_iters = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) {
            tol = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--no-damage") == 0) {
            use_damage = false;
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            std::cerr << "Usage: agp em-solve --index <input.emi> -o <output.tsv> [options]\n\n";
            std::cerr << "Run SQUAREM-accelerated EM on pre-indexed alignment data.\n\n";
            std::cerr << "Required:\n";
            std::cerr << "  --index, -i <file>     Input binary index (.emi from em-prepare)\n";
            std::cerr << "  --output, -o <file>    Output TSV with read assignments\n\n";
            std::cerr << "Options:\n";
            std::cerr << "  --summary <file>       Write per-reference summary TSV\n";
            std::cerr << "  --lambda <float>       Bit score temperature (default: 3.0)\n";
            std::cerr << "  --max-iters <int>      Maximum EM iterations (default: 100)\n";
            std::cerr << "  --tol <float>          Convergence tolerance (default: 1e-4)\n";
            std::cerr << "  --no-damage            Disable damage-aware EM\n";
            std::cerr << "  -v, --verbose          Verbose output\n";
            std::cerr << "  --help, -h             Show this help\n";
            return 0;
        }
    }

    if (index_file.empty() || output_file.empty()) {
        std::cerr << "Error: --index and --output are required\n";
        std::cerr << "Use --help for usage.\n";
        return 1;
    }

    auto t_start = std::chrono::steady_clock::now();

    // Load binary index (mmap'd)
    EMIndexReader index(index_file);

    if (verbose) {
        std::cerr << "Loaded index:\n";
        std::cerr << "  Reads: " << index.num_reads() << "\n";
        std::cerr << "  References: " << index.num_refs() << "\n";
        std::cerr << "  Alignments: " << index.num_alignments() << "\n";
        std::cerr << "  d_max: " << index.d_max() << "\n";
        std::cerr << "  lambda: " << index.lambda() << "\n";
    }

    auto t_load = std::chrono::steady_clock::now();

    // Run SQUAREM EM
    EMResult result = em_solve_indexed(index, lambda_b, max_iters, tol, use_damage);

    auto t_em = std::chrono::steady_clock::now();

    if (verbose) {
        auto load_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_load - t_start).count();
        auto em_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_em - t_load).count();
        std::cerr << "EM completed:\n";
        std::cerr << "  Iterations: " << result.iterations << "\n";
        std::cerr << "  Log-likelihood: " << std::fixed << std::setprecision(2) << result.log_likelihood << "\n";
        std::cerr << "  Pi (ancient): " << std::setprecision(4) << result.pi << "\n";
        std::cerr << "  Load time: " << load_ms << " ms\n";
        std::cerr << "  EM time: " << em_ms << " ms\n";
        std::cerr << "  Throughput: " << (index.num_alignments() * 1000 / (em_ms + 1)) << " aln/sec\n";
    }

    // Write read assignments
    std::ofstream out(output_file);
    if (!out) {
        std::cerr << "Error: Cannot open output file: " << output_file << "\n";
        return 1;
    }

    out << "read_id\tbest_ref\tgamma\tgamma_ancient\n";

    for (uint32_t r = 0; r < index.num_reads(); ++r) {
        uint64_t start = index.offset(r);
        uint64_t end = index.offset(r + 1);

        // Find best assignment
        double best_gamma = -1.0;
        double best_gamma_ancient = 0.0;
        uint32_t best_ref = 0;

        for (uint64_t j = start; j < end; ++j) {
            if (result.gamma[j] > best_gamma) {
                best_gamma = result.gamma[j];
                best_ref = index.alignments()[j].ref_idx;
                if (use_damage && !result.gamma_ancient.empty()) {
                    best_gamma_ancient = result.gamma_ancient[j];
                }
            }
        }

        out << index.read_name(r) << '\t'
            << index.ref_name(best_ref) << '\t'
            << std::fixed << std::setprecision(4) << best_gamma << '\t'
            << std::setprecision(4) << best_gamma_ancient << '\n';
    }

    out.close();

    // Write per-reference summary
    if (!summary_file.empty()) {
        std::ofstream sum(summary_file);
        if (!sum) {
            std::cerr << "Error: Cannot open summary file: " << summary_file << "\n";
            return 1;
        }

        // Accumulate per-reference stats
        std::vector<double> ref_n_eff(index.num_refs(), 0.0);
        std::vector<double> ref_n_ancient(index.num_refs(), 0.0);
        std::vector<uint32_t> ref_n_reads(index.num_refs(), 0);

        for (uint64_t j = 0; j < index.num_alignments(); ++j) {
            uint32_t ref_idx = index.alignments()[j].ref_idx;
            ref_n_eff[ref_idx] += result.gamma[j];
            if (use_damage && !result.gamma_ancient.empty()) {
                ref_n_ancient[ref_idx] += result.gamma_ancient[j];
            }
        }

        // Count reads per ref (best assignment)
        for (uint32_t r = 0; r < index.num_reads(); ++r) {
            uint64_t start = index.offset(r);
            uint64_t end = index.offset(r + 1);
            double best_gamma = -1.0;
            uint32_t best_ref = 0;
            for (uint64_t j = start; j < end; ++j) {
                if (result.gamma[j] > best_gamma) {
                    best_gamma = result.gamma[j];
                    best_ref = index.alignments()[j].ref_idx;
                }
            }
            ref_n_reads[best_ref]++;
        }

        sum << "ref_id\tn_reads\tn_effective\tn_ancient\tn_modern\tenrichment\n";

        for (uint32_t t = 0; t < index.num_refs(); ++t) {
            double n_eff = ref_n_eff[t];
            double n_anc = ref_n_ancient[t];
            double n_mod = n_eff - n_anc;

            // Enrichment: log2(ancient/modern) - log2(pi/(1-pi))
            double enrichment = 0.0;
            if (n_anc > 0.01 && n_mod > 0.01 && result.pi > 0.01 && result.pi < 0.99) {
                double prior_ratio = result.pi / (1.0 - result.pi);
                double obs_ratio = n_anc / n_mod;
                enrichment = std::log2(obs_ratio) - std::log2(prior_ratio);
            }

            if (ref_n_reads[t] > 0) {  // Only output refs with reads
                sum << index.ref_name(t) << '\t'
                    << ref_n_reads[t] << '\t'
                    << std::fixed << std::setprecision(2) << n_eff << '\t'
                    << std::setprecision(2) << n_anc << '\t'
                    << std::setprecision(2) << n_mod << '\t'
                    << std::setprecision(4) << enrichment << '\n';
            }
        }

        sum.close();

        if (verbose) {
            std::cerr << "Summary written: " << summary_file << "\n";
        }
    }

    auto t_end = std::chrono::steady_clock::now();
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    if (verbose) {
        std::cerr << "Total time: " << total_ms << " ms\n";
    }

    return 0;
}

namespace {
    struct EMSolveRegistrar {
        EMSolveRegistrar() {
            SubcommandRegistry::instance().register_command(
                "em-solve",
                "Run EM on binary index",
                cmd_em_solve, 51);
        }
    };
    static EMSolveRegistrar registrar;
}

}  // namespace cli
}  // namespace agp
