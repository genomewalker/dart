// agp em-solve2: Run EM on columnar index with parallel row groups
//
// Features:
// - Parallel row group processing
// - Predicate pushdown for filtering
// - SQUAREM acceleration

#include "subcommand.hpp"
#include "agp/version.h"
#include "agp/columnar_index.hpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <numeric>

namespace agp {
namespace cli {

int cmd_em_solve2(int argc, char* argv[]) {
    std::string index_file;
    std::string output_file;
    double lambda_b = 3.0;
    uint32_t max_iters = 100;
    double tol = 1e-4;
    double alpha_prior = 1.0;  // Dirichlet prior (prevents rich-get-richer)
    bool use_damage = true;
    bool verbose = false;

    // Filter options
    float min_bit_score = 0.0f;
    float max_evalue = 1e10f;
    float min_identity = 0.0f;

    for (int i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "--index") == 0 || strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
            index_file = argv[++i];
        } else if ((strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--lambda") == 0 && i + 1 < argc) {
            lambda_b = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--max-iters") == 0 && i + 1 < argc) {
            max_iters = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) {
            tol = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--alpha") == 0 && i + 1 < argc) {
            alpha_prior = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--no-damage") == 0) {
            use_damage = false;
        } else if (strcmp(argv[i], "--min-bit-score") == 0 && i + 1 < argc) {
            min_bit_score = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--max-evalue") == 0 && i + 1 < argc) {
            max_evalue = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--min-identity") == 0 && i + 1 < argc) {
            min_identity = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            std::cerr << "Usage: agp em-solve2 --index <file.emi2> -o <output.tsv> [options]\n\n";
            std::cerr << "Run SQUAREM-accelerated EM on columnar index.\n\n";
            std::cerr << "Features:\n";
            std::cerr << "  - Parallel row group processing\n";
            std::cerr << "  - Predicate pushdown (skip chunks by statistics)\n";
            std::cerr << "  - SQUAREM acceleration (2-3x fewer iterations)\n\n";
            std::cerr << "Required:\n";
            std::cerr << "  --index, -i <file>     Columnar index (.emi2)\n";
            std::cerr << "  --output, -o <file>    Output TSV with assignments\n\n";
            std::cerr << "EM Options:\n";
            std::cerr << "  --lambda <float>       Score temperature (default: 3.0)\n";
            std::cerr << "  --max-iters <int>      Maximum iterations (default: 100)\n";
            std::cerr << "  --tol <float>          Convergence tolerance (default: 1e-4)\n";
            std::cerr << "  --alpha <float>        Dirichlet prior (default: 1.0, prevents rich-get-richer)\n";
            std::cerr << "  --no-damage            Disable damage-aware EM\n\n";
            std::cerr << "Filter Options (predicate pushdown):\n";
            std::cerr << "  --min-bit-score <f>    Minimum bit score\n";
            std::cerr << "  --max-evalue <f>       Maximum e-value\n";
            std::cerr << "  --min-identity <f>     Minimum identity (0-1)\n\n";
            std::cerr << "Other:\n";
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

    // Open columnar index
    ColumnarIndexReader reader(index_file);
    if (!reader.is_valid()) {
        std::cerr << "Error: Cannot open columnar index: " << index_file << "\n";
        return 1;
    }

    if (verbose) {
        std::cerr << "Loaded columnar index:\n";
        std::cerr << "  Alignments: " << reader.num_alignments() << "\n";
        std::cerr << "  Reads: " << reader.num_reads() << "\n";
        std::cerr << "  References: " << reader.num_refs() << "\n";
        std::cerr << "  Row groups: " << reader.num_row_groups() << "\n";
        std::cerr << "  d_max: " << reader.d_max() << "\n";
        std::cerr << "  lambda: " << reader.lambda() << "\n";
    }

    // Add filters for predicate pushdown
    if (min_bit_score > 0.0f) {
        FilterPredicate fp;
        fp.column = ColumnID::BIT_SCORE;
        fp.op = FilterPredicate::Op::GE;
        fp.value_f = min_bit_score;
        reader.add_filter(fp);
        if (verbose) {
            std::cerr << "Filter: bit_score >= " << min_bit_score << "\n";
        }
    }

    if (max_evalue < 1e9f) {
        FilterPredicate fp;
        fp.column = ColumnID::EVALUE_LOG10;
        fp.op = FilterPredicate::Op::LE;
        fp.value_f = std::log10(max_evalue);
        reader.add_filter(fp);
        if (verbose) {
            std::cerr << "Filter: evalue <= " << max_evalue << "\n";
        }
    }

    if (min_identity > 0.0f) {
        FilterPredicate fp;
        fp.column = ColumnID::IDENTITY_Q;
        fp.op = FilterPredicate::Op::GE;
        fp.value_u = static_cast<uint32_t>(min_identity * 65535.0f);
        reader.add_filter(fp);
        if (verbose) {
            std::cerr << "Filter: identity >= " << min_identity << "\n";
        }
    }

    // Run EM
    auto t_em_start = std::chrono::steady_clock::now();
    auto result = em_solve_columnar(reader, lambda_b, max_iters, tol, use_damage, alpha_prior);
    auto t_em_end = std::chrono::steady_clock::now();

    if (verbose) {
        auto em_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_em_end - t_em_start).count();
        std::cerr << "EM completed:\n";
        std::cerr << "  Iterations: " << result.iterations << "\n";
        std::cerr << "  Time: " << em_ms << " ms\n";
        std::cerr << "  Log-likelihood: " << result.log_likelihood << "\n";
        if (use_damage) {
            std::cerr << "  Ancient fraction (pi): " << result.pi << "\n";
        }
    }

    // Output results
    std::ofstream out(output_file);
    if (!out) {
        std::cerr << "Error: Cannot create output file: " << output_file << "\n";
        return 1;
    }

    // Sort references by weight (descending)
    std::vector<uint32_t> sorted_refs(reader.num_refs());
    std::iota(sorted_refs.begin(), sorted_refs.end(), 0);
    std::sort(sorted_refs.begin(), sorted_refs.end(), [&](uint32_t a, uint32_t b) {
        return result.weights[a] > result.weights[b];
    });

    // Write per-reference summary
    out << "# EM results\n";
    out << "# Iterations: " << result.iterations << "\n";
    out << "# Log-likelihood: " << result.log_likelihood << "\n";
    out << "ref_name\tweight\teffective_count\n";

    double total_reads = static_cast<double>(reader.num_reads());
    for (uint32_t ref_idx : sorted_refs) {
        double weight = result.weights[ref_idx];
        if (weight < 1e-8) continue;  // Skip negligible weights

        out << reader.ref_name(ref_idx) << "\t"
            << std::setprecision(6) << weight << "\t"
            << std::setprecision(1) << std::fixed << (weight * total_reads) << "\n";
    }

    out.close();

    auto t_end = std::chrono::steady_clock::now();
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    if (verbose) {
        std::cerr << "Written: " << output_file << "\n";
        std::cerr << "Total time: " << total_ms << " ms\n";
    }

    return 0;
}

namespace {
    struct EMSolve2Registrar {
        EMSolve2Registrar() {
            SubcommandRegistry::instance().register_command(
                "em-solve2",
                "Run EM on columnar index with parallel processing",
                cmd_em_solve2, 52);
        }
    };
    static EMSolve2Registrar registrar;
}

}  // namespace cli
}  // namespace agp
