// dart damage-profile: Terminal nucleotide damage profile
//
// Computes T/(T+C) at 5' and A/(A+G) at 3' positions for reads
// that mapped to proteins. Outputs per-protein or aggregate profiles.

#include "subcommand.hpp"
#include "dart/version.h"
#include "dart/log_utils.hpp"
#include "dart/sequence_io.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <zlib.h>

namespace dart {
namespace cli {

static constexpr int MAX_POSITIONS = 25;

struct NucleotideCounts {
    std::array<uint64_t, MAX_POSITIONS> T_5prime{};
    std::array<uint64_t, MAX_POSITIONS> C_5prime{};
    std::array<uint64_t, MAX_POSITIONS> A_3prime{};
    std::array<uint64_t, MAX_POSITIONS> G_3prime{};
    uint64_t n_reads = 0;

    void add(const NucleotideCounts& other) {
        for (int i = 0; i < MAX_POSITIONS; ++i) {
            T_5prime[i] += other.T_5prime[i];
            C_5prime[i] += other.C_5prime[i];
            A_3prime[i] += other.A_3prime[i];
            G_3prime[i] += other.G_3prime[i];
        }
        n_reads += other.n_reads;
    }
};

// Extract read name from query_id (remove _+_0, _-_1 etc suffix)
static std::string extract_read_name(const std::string& query_id) {
    // Pattern: readname_[+-]_[012]
    size_t pos = query_id.rfind('_');
    if (pos != std::string::npos && pos > 0) {
        size_t pos2 = query_id.rfind('_', pos - 1);
        if (pos2 != std::string::npos) {
            char strand = query_id[pos2 + 1];
            if (strand == '+' || strand == '-') {
                return query_id.substr(0, pos2);
            }
        }
    }
    return query_id;
}

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n\n"
              << "Compute terminal nucleotide damage profiles for protein-mapped reads.\n\n"
              << "Required:\n"
              << "  -i, --input FILE       Input FASTQ file (gzipped supported)\n"
              << "  --hits FILE            MMseqs2 hits file (query_id in first column)\n"
              << "  -o, --output FILE      Output file (.tsv.gz)\n\n"
              << "Options:\n"
              << "  --aggregate            Output aggregate profile instead of per-protein\n"
              << "  --min-reads N          Min reads per protein (default: 10)\n"
              << "  --max-pos N            Max positions to track (default: 20)\n"
              << "  -h, --help             Show this help\n\n"
              << "Output columns (per-protein mode):\n"
              << "  target, n_reads, TC_ratio_0..N, AG_ratio_0..N\n\n"
              << "Output columns (aggregate mode):\n"
              << "  position, T_5prime, C_5prime, TC_ratio, A_3prime, G_3prime, AG_ratio\n";
}

int cmd_damage_profile(int argc, char* argv[]) {
    auto run_start = std::chrono::steady_clock::now();
    std::string input_file;
    std::string hits_file;
    std::string output_file;
    bool aggregate = false;
    int min_reads = 10;
    int max_pos = 20;

    // Parse arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if ((strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) && i + 1 < argc) {
            input_file = argv[++i];
        } else if (strcmp(argv[i], "--hits") == 0 && i + 1 < argc) {
            hits_file = argv[++i];
        } else if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--aggregate") == 0) {
            aggregate = true;
        } else if (strcmp(argv[i], "--min-reads") == 0 && i + 1 < argc) {
            min_reads = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--max-pos") == 0 && i + 1 < argc) {
            max_pos = std::stoi(argv[++i]);
        } else if (argv[i][0] == '-') {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            std::cerr << "Run 'dart damage-profile --help' for usage.\n";
            return 1;
        } else {
            std::cerr << "Unexpected positional argument: " << argv[i] << "\n";
            std::cerr << "Run 'dart damage-profile --help' for usage.\n";
            return 1;
        }
    }

    if (input_file.empty() || hits_file.empty() || output_file.empty()) {
        std::cerr << "Error: --input, --hits, and --output are required.\n";
        std::cerr << "Run 'dart damage-profile --help' for usage.\n";
        return 1;
    }

    if (max_pos > MAX_POSITIONS) max_pos = MAX_POSITIONS;

    // Step 1: Parse hits file to get read->target mapping
    auto step1_start = std::chrono::steady_clock::now();
    std::cerr << "Reading hits file..." << std::endl;
    std::unordered_map<std::string, std::string> read_to_target;

    std::ifstream hits(hits_file);
    if (!hits) {
        std::cerr << "Error: Cannot open hits file: " << hits_file << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(hits, line)) {
        if (line.empty()) continue;

        // Find first two tab-separated fields
        size_t tab1 = line.find('\t');
        if (tab1 == std::string::npos) continue;
        size_t tab2 = line.find('\t', tab1 + 1);

        std::string query_id = line.substr(0, tab1);
        std::string target = line.substr(tab1 + 1, tab2 - tab1 - 1);

        std::string read_name = extract_read_name(query_id);

        // Keep first hit per read (best hit)
        if (read_to_target.find(read_name) == read_to_target.end()) {
            read_to_target[read_name] = target;
        }
    }
    hits.close();
    auto step1_end = std::chrono::steady_clock::now();

    std::cerr << "  Found " << read_to_target.size() << " read-protein mappings" << std::endl;
    std::cerr << "  Step 1 runtime: " << dart::log_utils::format_elapsed(step1_start, step1_end) << std::endl;

    // Step 2: Read FASTQ and compute profiles
    auto step2_start = std::chrono::steady_clock::now();
    std::cerr << "Processing reads..." << std::endl;

    std::unordered_map<std::string, NucleotideCounts> target_counts;
    NucleotideCounts total_counts;
    uint64_t reads_processed = 0;
    uint64_t reads_matched = 0;

    // Open input (gzip or plain)
    gzFile gz = gzopen(input_file.c_str(), "r");
    if (!gz) {
        std::cerr << "Error: Cannot open input file: " << input_file << std::endl;
        return 1;
    }

    char buffer[65536];
    std::string header, seq, plus_line, qual;

    while (gzgets(gz, buffer, sizeof(buffer))) {
        header = buffer;
        if (header.empty() || header[0] != '@') continue;

        // Read sequence
        if (!gzgets(gz, buffer, sizeof(buffer))) break;
        seq = buffer;
        if (!seq.empty() && seq.back() == '\n') seq.pop_back();

        // Skip + and qual
        gzgets(gz, buffer, sizeof(buffer));
        gzgets(gz, buffer, sizeof(buffer));

        reads_processed++;

        // Extract read name
        size_t space = header.find(' ');
        std::string read_name = header.substr(1, space - 1);
        if (read_name.back() == '\n') read_name.pop_back();

        auto it = read_to_target.find(read_name);
        if (it == read_to_target.end()) continue;

        reads_matched++;
        const std::string& target = it->second;

        // Count nucleotides
        NucleotideCounts counts;
        counts.n_reads = 1;
        size_t seq_len = seq.size();

        for (int pos = 0; pos < max_pos && pos < (int)seq_len; ++pos) {
            char base = std::toupper(static_cast<unsigned char>(seq[pos]));
            if (base == 'T') counts.T_5prime[pos]++;
            else if (base == 'C') counts.C_5prime[pos]++;

            char base3 = std::toupper(static_cast<unsigned char>(seq[seq_len - 1 - pos]));
            if (base3 == 'A') counts.A_3prime[pos]++;
            else if (base3 == 'G') counts.G_3prime[pos]++;
        }

        target_counts[target].add(counts);
        total_counts.add(counts);

        if (reads_processed % 1000000 == 0) {
            std::cerr << "  Processed " << reads_processed / 1000000 << "M reads, "
                      << reads_matched << " matched" << std::endl;
        }
    }

    gzclose(gz);
    auto step2_end = std::chrono::steady_clock::now();

    std::cerr << "  Total: " << reads_processed << " reads, "
              << reads_matched << " matched to " << target_counts.size() << " proteins" << std::endl;
    std::cerr << "  Step 2 runtime: " << dart::log_utils::format_elapsed(step2_start, step2_end) << std::endl;

    // Step 3: Write output
    auto step3_start = std::chrono::steady_clock::now();
    std::cerr << "Writing output..." << std::endl;

    gzFile out = gzopen(output_file.c_str(), "wb");
    if (!out) {
        std::cerr << "Error: Cannot open output file: " << output_file << std::endl;
        return 1;
    }

    std::ostringstream oss;

    if (aggregate) {
        // Aggregate mode: one row per position
        oss << "position\tT_5prime\tC_5prime\tTC_ratio\tA_3prime\tG_3prime\tAG_ratio\n";

        for (int pos = 0; pos < max_pos; ++pos) {
            uint64_t T = total_counts.T_5prime[pos];
            uint64_t C = total_counts.C_5prime[pos];
            uint64_t A = total_counts.A_3prime[pos];
            uint64_t G = total_counts.G_3prime[pos];

            double tc_ratio = (T + C > 0) ? (double)T / (T + C) : 0.5;
            double ag_ratio = (A + G > 0) ? (double)A / (A + G) : 0.5;

            oss << pos << "\t" << T << "\t" << C << "\t"
                << std::fixed << std::setprecision(4) << tc_ratio << "\t"
                << A << "\t" << G << "\t" << ag_ratio << "\n";
        }
    } else {
        // Per-protein mode: one row per protein
        oss << "target\tn_reads";
        for (int pos = 0; pos < max_pos; ++pos) {
            oss << "\tTC_" << pos;
        }
        for (int pos = 0; pos < max_pos; ++pos) {
            oss << "\tAG_" << pos;
        }
        oss << "\n";

        // Sort by read count descending
        std::vector<std::pair<std::string, NucleotideCounts>> sorted_targets(
            target_counts.begin(), target_counts.end());
        std::sort(sorted_targets.begin(), sorted_targets.end(),
                  [](const auto& a, const auto& b) { return a.second.n_reads > b.second.n_reads; });

        for (const auto& [target, counts] : sorted_targets) {
            if ((int)counts.n_reads < min_reads) continue;

            oss << target << "\t" << counts.n_reads;

            // TC ratios
            for (int pos = 0; pos < max_pos; ++pos) {
                uint64_t T = counts.T_5prime[pos];
                uint64_t C = counts.C_5prime[pos];
                double ratio = (T + C > 0) ? (double)T / (T + C) : 0.5;
                oss << "\t" << std::fixed << std::setprecision(4) << ratio;
            }

            // AG ratios
            for (int pos = 0; pos < max_pos; ++pos) {
                uint64_t A = counts.A_3prime[pos];
                uint64_t G = counts.G_3prime[pos];
                double ratio = (A + G > 0) ? (double)A / (A + G) : 0.5;
                oss << "\t" << std::fixed << std::setprecision(4) << ratio;
            }

            oss << "\n";
        }
    }

    std::string output_str = oss.str();
    gzwrite(out, output_str.c_str(), output_str.size());
    gzclose(out);
    auto step3_end = std::chrono::steady_clock::now();
    auto run_end = std::chrono::steady_clock::now();

    std::cerr << "Done. Output: " << output_file << std::endl;
    std::cerr << "  Step 3 runtime: " << dart::log_utils::format_elapsed(step3_start, step3_end) << std::endl;
    std::cerr << "  Total runtime: " << dart::log_utils::format_elapsed(run_start, run_end) << std::endl;

    return 0;
}

// Register subcommand
static struct DamageProfileRegistrar {
    DamageProfileRegistrar() {
        SubcommandRegistry::instance().register_command(
            "damage-profile",
            "Compute terminal nucleotide damage profiles",
            cmd_damage_profile,
            4  // Order after damage-annotate
        );
    }
} damage_profile_registrar;

}  // namespace cli
}  // namespace dart
