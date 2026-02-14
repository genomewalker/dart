// agp hits2emi: Convert MMseqs2 hits TSV to columnar EM index (.emi)
//
// DuckDB-style streaming architecture:
// - Hybrid mode: mmap for small files, streaming for large
// - Spilling dictionary for billion-row support
// - Row groups for parallel EM processing
// - Lossless storage including qaln/taln

#include "subcommand.hpp"
#include "agp/version.h"
#include "agp/fast_columnar_writer.hpp"

#include <iostream>
#include <cstring>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace agp {
namespace cli {

int cmd_hits2emi(int argc, char* argv[]) {
    std::string tsv_file;
    std::string output_file;
    float d_max = 0.0f;
    float lambda = 0.3f;
    bool verbose = false;
    bool fast = false;
    size_t memory_gb = 16;  // Default 16GB
    bool force_streaming = false;
    int threads = 0;  // 0 = use all available

    for (int i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "--tsv") == 0 || strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
            tsv_file = argv[++i];
        } else if ((strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--d-max") == 0 && i + 1 < argc) {
            d_max = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda") == 0 && i + 1 < argc) {
            lambda = std::stof(argv[++i]);
        } else if ((strcmp(argv[i], "--memory") == 0 || strcmp(argv[i], "-m") == 0) && i + 1 < argc) {
            memory_gb = std::stoul(argv[++i]);
        } else if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0) && i + 1 < argc) {
            threads = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--streaming") == 0) {
            force_streaming = true;
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--fast") == 0) {
            fast = true;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            std::cerr << "Usage: agp hits2emi -i <hits.tsv[.gz]> -o <output.emi> [options]\n\n";
            std::cerr << "Convert MMseqs2 hits to columnar EM index (DuckDB-style streaming).\n\n";
            std::cerr << "Features:\n";
            std::cerr << "  - Hybrid mode: mmap for small files, streaming for large\n";
            std::cerr << "  - Auto-detects gzip compression (.gz files)\n";
            std::cerr << "  - Row groups (64K) for parallel EM processing\n";
            std::cerr << "  - Spilling dictionary for billion-row support\n";
            std::cerr << "  - Lossless storage (all TSV fields + qaln/taln)\n\n";
            std::cerr << "Required:\n";
            std::cerr << "  --tsv, -i <file>       Input MMseqs2 hits (16-column, .tsv or .tsv.gz)\n";
            std::cerr << "  --output, -o <file>    Output columnar index (.emi)\n\n";
            std::cerr << "Options:\n";
            std::cerr << "  --threads, -t <N>      Number of threads (default: all available)\n";
            std::cerr << "  --memory, -m <GB>      Memory limit in GB (default: 16)\n";
            std::cerr << "                         Files < limit use fast mmap, larger use streaming\n";
            std::cerr << "  --streaming            Force streaming mode (bounded memory)\n";
            std::cerr << "  --d-max <float>        Sample damage level (0-1), stored in index\n";
            std::cerr << "  --lambda <float>       Damage decay parameter (default: 0.3)\n";
            std::cerr << "  --fast                 Skip alignment strings (faster, EM doesn't need them)\n";
            std::cerr << "  -v, --verbose          Verbose output\n";
            std::cerr << "  --help, -h             Show this help\n";
            return 0;
        }
    }

    if (tsv_file.empty() || output_file.empty()) {
        std::cerr << "Error: --tsv and --output are required\n";
        std::cerr << "Use --help for usage.\n";
        return 1;
    }

    auto t_start = std::chrono::steady_clock::now();

    try {
        // Set thread count
#ifdef _OPENMP
        if (threads > 0) {
            omp_set_num_threads(threads);
        }
#endif

        FastColumnarWriter::Config config;
        config.d_max = d_max;
        config.lambda = lambda;
        config.verbose = verbose;
        config.skip_alignments = fast;
        config.memory_limit = memory_gb * 1024ULL * 1024ULL * 1024ULL;
        config.force_streaming = force_streaming;

        uint64_t n_rows = FastColumnarWriter::convert(tsv_file, output_file, config);

        auto t_end = std::chrono::steady_clock::now();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

        if (verbose) {
            std::cerr << "Written: " << output_file << "\n";
            std::cerr << "Total time: " << total_ms << " ms\n";
            std::cerr << "Throughput: " << (n_rows * 1000 / (total_ms + 1)) << " rows/sec\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

namespace {
    struct Hits2EmiRegistrar {
        Hits2EmiRegistrar() {
            SubcommandRegistry::instance().register_command(
                "hits2emi",
                "Convert MMseqs2 hits to columnar EM index (.emi)",
                cmd_hits2emi, 51);
        }
    };
    static Hits2EmiRegistrar registrar;
}

}  // namespace cli
}  // namespace agp
