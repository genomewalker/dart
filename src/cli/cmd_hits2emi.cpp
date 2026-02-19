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
#include "agp/log_utils.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <cctype>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace agp {
namespace cli {

namespace {

static uint64_t detect_available_memory_bytes() {
    // Linux: prefer MemAvailable from /proc/meminfo for a realistic budget.
    {
        std::ifstream meminfo("/proc/meminfo");
        if (meminfo.is_open()) {
            std::string key;
            uint64_t value_kb = 0;
            std::string unit;
            while (meminfo >> key >> value_kb >> unit) {
                if (key == "MemAvailable:") {
                    return value_kb * 1024ULL;
                }
            }
        }
    }

    // Fallback: total physical RAM.
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGESIZE);
    if (pages > 0 && page_size > 0) {
        return static_cast<uint64_t>(pages) * static_cast<uint64_t>(page_size);
    }

    // Last-resort fallback.
    return 16ULL * 1024ULL * 1024ULL * 1024ULL;
}

static uint64_t default_memory_limit_bytes() {
    uint64_t available = detect_available_memory_bytes();
    uint64_t limit = available / 10ULL;  // 10%
    constexpr uint64_t MIN_LIMIT = 256ULL * 1024ULL * 1024ULL;  // 256 MB floor
    if (limit < MIN_LIMIT) limit = MIN_LIMIT;
    return limit;
}

static uint64_t parse_memory_arg(const std::string& raw) {
    if (raw.empty()) {
        throw std::runtime_error("Invalid --memory value: empty");
    }

    std::string s = raw;
    char unit = '\0';

    // Support optional trailing 'B' (e.g., 4GB, 512MB).
    if (!s.empty() && (s.back() == 'B' || s.back() == 'b')) {
        s.pop_back();
    }

    if (!s.empty() && std::isalpha(static_cast<unsigned char>(s.back()))) {
        unit = static_cast<char>(std::toupper(static_cast<unsigned char>(s.back())));
        s.pop_back();
    }

    if (s.empty()) {
        throw std::runtime_error("Invalid --memory value: " + raw);
    }

    char* end = nullptr;
    double value = std::strtod(s.c_str(), &end);
    if (!end || *end != '\0' || !std::isfinite(value) || value <= 0.0) {
        throw std::runtime_error("Invalid --memory value: " + raw);
    }

    uint64_t multiplier = 0;
    if (unit == 'K') {
        multiplier = 1024ULL;
    } else if (unit == 'M') {
        multiplier = 1024ULL * 1024ULL;
    } else if (unit == 'G') {
        multiplier = 1024ULL * 1024ULL * 1024ULL;
    } else if (unit == 'T') {
        multiplier = 1024ULL * 1024ULL * 1024ULL * 1024ULL;
    } else if (unit == '\0') {
        throw std::runtime_error("Invalid --memory value: " + raw +
                                 " (explicit suffix required: K/M/G/T, optional B)");
    } else {
        throw std::runtime_error("Invalid --memory suffix in value: " + raw +
                                 " (use K/M/G/T, optional B)");
    }

    long double bytes_ld = static_cast<long double>(value) *
                           static_cast<long double>(multiplier);
    if (bytes_ld > static_cast<long double>(std::numeric_limits<uint64_t>::max())) {
        throw std::runtime_error("Invalid --memory value too large: " + raw);
    }
    return static_cast<uint64_t>(bytes_ld);
}

}  // namespace

int cmd_hits2emi(int argc, char* argv[]) {
    std::string tsv_file;
    std::string output_file;
    float d_max = 0.0f;
    float lambda = 0.3f;
    bool d_max_set_by_user = false;
    bool lambda_set_by_user = false;
    bool verbose = false;
    uint64_t memory_limit_bytes = default_memory_limit_bytes();  // Default: 10% available RAM
    bool force_streaming = false;
    int threads = 0;  // 0 = use all available
    std::string damage_index_file;

    for (int i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "--tsv") == 0 || strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
            tsv_file = argv[++i];
        } else if ((strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--d-max") == 0 && i + 1 < argc) {
            d_max = std::stof(argv[++i]);
            d_max_set_by_user = true;
        } else if (strcmp(argv[i], "--lambda") == 0 && i + 1 < argc) {
            lambda = std::stof(argv[++i]);
            lambda_set_by_user = true;
        } else if ((strcmp(argv[i], "--memory") == 0 || strcmp(argv[i], "-m") == 0) && i + 1 < argc) {
            memory_limit_bytes = parse_memory_arg(argv[++i]);
        } else if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0) && i + 1 < argc) {
            threads = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--damage-index") == 0 && i + 1 < argc) {
            damage_index_file = argv[++i];
        } else if (strcmp(argv[i], "--streaming") == 0) {
            force_streaming = true;
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            std::cout << "Usage: agp hits2emi -i <hits.tsv[.gz]> -o <output.emi> [options]\n\n";
            std::cout << "Convert MMseqs2 hits to columnar EM index (DuckDB-style streaming).\n\n";
            std::cout << "Features:\n";
            std::cout << "  - Hybrid mode for plain files: mmap for small, streaming for large\n";
            std::cout << "  - Auto-detects gzip compression (.gz files, defaults to streaming)\n";
            std::cout << "  - Row groups (64K) for parallel EM processing\n";
            std::cout << "  - Spilling dictionary for billion-row support\n";
            std::cout << "  - Lossless storage (all TSV fields + qaln/taln)\n\n";
            std::cout << "Required:\n";
            std::cout << "  --tsv, -i <file>       Input MMseqs2 hits (16-column, .tsv or .tsv.gz)\n";
            std::cout << "  --output, -o <file>    Output columnar index (.emi)\n\n";
            std::cout << "Options:\n";
            std::cout << "  --threads, -t <N>      Number of threads (default: all available)\n";
            std::cout << "  --memory, -m <SIZE>    Memory limit (default: 10% of available RAM)\n";
            std::cout << "                         Suffixes: K, M, G, T (optional B; e.g. 512M, 4GB)\n";
            std::cout << "                         Plain files < limit use fast mmap, larger use streaming\n";
            std::cout << "  --streaming            Force streaming mode (bounded memory)\n";
            std::cout << "  --d-max <float>        Sample damage level (0-1), stored in index\n";
            std::cout << "  --lambda <float>       Damage decay parameter (default: 0.3)\n";
            std::cout << "  --damage-index <file>  .agd from predict; embeds per-read p_damaged\n";
            std::cout << "                         and alignment-level damage likelihood columns\n";
            std::cout << "  -v, --verbose          Verbose output\n";
            std::cout << "  --help, -h             Show this help\n";
            return 0;
        } else if (argv[i][0] == '-') {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            std::cerr << "Run 'agp hits2emi --help' for usage.\n";
            return 1;
        } else {
            std::cerr << "Unexpected positional argument: " << argv[i] << "\n";
            std::cerr << "Run 'agp hits2emi --help' for usage.\n";
            return 1;
        }
    }

    if (tsv_file.empty() || output_file.empty()) {
        std::cerr << "Error: --tsv and --output are required\n";
        std::cerr << "Run 'agp hits2emi --help' for usage.\n";
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
        config.d_max_set_by_user = d_max_set_by_user;
        config.lambda_set_by_user = lambda_set_by_user;
        config.verbose = verbose;
        config.skip_alignments = false;  // Alignment strings required for downstream annotate-damage.
        config.memory_limit = memory_limit_bytes;
        config.force_streaming = force_streaming;
        config.damage_index_path = damage_index_file;

        uint64_t n_rows = FastColumnarWriter::convert(tsv_file, output_file, config);

        auto t_end = std::chrono::steady_clock::now();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

        if (verbose) {
            std::cerr << "Written: " << output_file << "\n";
            std::cerr << "Total time: " << agp::log_utils::format_duration_ms(total_ms) << "\n";
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
