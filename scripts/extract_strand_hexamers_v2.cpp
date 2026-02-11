// Strand hexamer extractor v2 - produces ANTISYMMETRIC tables
// Key insight: For strand discrimination, we need LLR(X) = -LLR(RC(X))
//
// Approach:
// 1. Count hexamers in forward strand genes (5'→3' direction)
// 2. For each hexamer X, compute:
//    - count_fwd(X) = occurrences in forward reading
//    - count_fwd(RC(X)) = occurrences of RC in forward reading
// 3. Antisymmetric LLR: LLR(X) = log(count_fwd(X) / count_fwd(RC(X)))
//    This guarantees LLR(X) = -LLR(RC(X))

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <chrono>

namespace fs = std::filesystem;

inline int encode_base(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

inline int encode_hexamer(const char* s) {
    int code = 0;
    for (int i = 0; i < 6; i++) {
        int b = encode_base(s[i]);
        if (b < 0) return -1;
        code = (code << 2) | b;
    }
    return code;
}

inline int rc_hexamer_code(int code) {
    int rc = 0;
    for (int i = 0; i < 6; i++) {
        int base = (code >> (2 * (5 - i))) & 3;
        int rc_base = 3 - base;  // A↔T (0↔3), C↔G (1↔2)
        rc = (rc << 2) | rc_base;
    }
    return rc;
}

std::string decode_hexamer(int code) {
    static const char bases[] = "ACGT";
    char hex[7];
    for (int i = 5; i >= 0; i--) {
        hex[i] = bases[code & 3];
        code >>= 2;
    }
    hex[6] = '\0';
    return std::string(hex);
}

// Thread-local counters
struct ThreadCounters {
    uint64_t coding[4096] = {0};  // Forward strand coding only
    uint64_t total = 0;
};

uint64_t g_coding[4096] = {0};
uint64_t g_total = 0;

std::atomic<uint64_t> files_processed{0};
std::atomic<uint64_t> genes_processed{0};

bool parse_prodigal_header(const std::string& header,
                           std::string& contig,
                           int64_t& start, int64_t& end, int& strand) {
    size_t hash1 = header.find('#');
    if (hash1 == std::string::npos) return false;
    size_t hash2 = header.find('#', hash1 + 1);
    if (hash2 == std::string::npos) return false;
    size_t hash3 = header.find('#', hash2 + 1);
    if (hash3 == std::string::npos) return false;
    size_t hash4 = header.find('#', hash3 + 1);
    if (hash4 == std::string::npos) return false;

    contig = header.substr(1, hash1 - 1);
    while (!contig.empty() && contig.back() == ' ') contig.pop_back();
    size_t last_underscore = contig.rfind('_');
    if (last_underscore != std::string::npos) {
        contig = contig.substr(0, last_underscore);
    }

    try {
        start = std::stoll(header.substr(hash1 + 1, hash2 - hash1 - 1));
        end = std::stoll(header.substr(hash2 + 1, hash3 - hash2 - 1));
        strand = std::stoi(header.substr(hash3 + 1, hash4 - hash3 - 1));
    } catch (...) {
        return false;
    }
    return true;
}

void process_file(const std::string& path, ThreadCounters& tc) {
    std::string cmd;
    if (path.ends_with(".gz")) {
        cmd = "pigz -dc '" + path + "' 2>/dev/null || zcat '" + path + "' 2>/dev/null";
    } else {
        cmd = "cat '" + path + "'";
    }

    FILE* fp = popen(cmd.c_str(), "r");
    if (!fp) return;

    char line[65536];
    std::string current_seq;
    std::string current_contig;
    int64_t current_start = 0, current_end = 0;
    int current_strand = 0;
    bool in_gene = false;

    auto process_gene = [&]() {
        if (current_seq.length() < 6) return;

        // Extract hexamers from coding sequence (forward reading)
        for (size_t i = 0; i + 5 < current_seq.length(); i++) {
            int code = encode_hexamer(current_seq.c_str() + i);
            if (code >= 0) {
                tc.coding[code]++;
                tc.total++;
            }
        }
        genes_processed++;
    };

    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
            line[--len] = '\0';
        }

        if (line[0] == '>') {
            if (in_gene && !current_seq.empty()) {
                process_gene();
            }
            std::string header(line);
            in_gene = parse_prodigal_header(header, current_contig,
                                            current_start, current_end, current_strand);
            current_seq.clear();
        } else if (in_gene) {
            current_seq += line;
        }
    }

    if (in_gene && !current_seq.empty()) {
        process_gene();
    }

    pclose(fp);
    files_processed++;
}

void worker_thread(const std::vector<std::string>& files,
                   size_t start_idx, size_t end_idx,
                   ThreadCounters& tc) {
    for (size_t i = start_idx; i < end_idx; i++) {
        process_file(files[i], tc);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s --cds <dir> --output <prefix> [--threads N]\n", argv[0]);
        fprintf(stderr, "\nGenerates ANTISYMMETRIC strand LLR table where LLR(X) = -LLR(RC(X))\n");
        return 1;
    }

    std::string cds_dir, output_prefix;
    int num_threads = std::thread::hardware_concurrency();

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--cds") == 0 && i + 1 < argc) {
            cds_dir = argv[++i];
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        }
    }

    std::vector<std::string> files;
    for (const auto& entry : fs::recursive_directory_iterator(cds_dir)) {
        if (entry.is_regular_file()) {
            std::string p = entry.path().string();
            if (p.ends_with(".fna") || p.ends_with(".fna.gz") ||
                p.ends_with(".fasta") || p.ends_with(".fasta.gz")) {
                files.push_back(p);
            }
        }
    }

    fprintf(stderr, "Found %zu files, using %d threads\n", files.size(), num_threads);
    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<ThreadCounters> thread_counters(num_threads);
    std::vector<std::thread> threads;

    size_t files_per_thread = (files.size() + num_threads - 1) / num_threads;
    for (int t = 0; t < num_threads; t++) {
        size_t start = t * files_per_thread;
        size_t end = std::min(start + files_per_thread, files.size());
        if (start < files.size()) {
            threads.emplace_back(worker_thread, std::cref(files), start, end,
                                 std::ref(thread_counters[t]));
        }
    }

    std::atomic<bool> done{false};
    std::thread progress_thread([&]() {
        while (!done) {
            std::this_thread::sleep_for(std::chrono::seconds(5));
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(now - t0).count();
            fprintf(stderr, "\r[%.0fs] Files: %lu/%zu, Genes: %lu",
                    elapsed, files_processed.load(), files.size(), genes_processed.load());
        }
    });

    for (auto& t : threads) t.join();
    done = true;
    progress_thread.join();

    auto t1 = std::chrono::high_resolution_clock::now();
    fprintf(stderr, "\nExtraction: %.1fs\n", std::chrono::duration<double>(t1 - t0).count());

    // Merge counters
    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i < 4096; i++) {
            g_coding[i] += thread_counters[t].coding[i];
        }
        g_total += thread_counters[t].total;
    }

    fprintf(stderr, "Total hexamers: %lu\n", g_total);

    // Compute ANTISYMMETRIC strand LLR
    // LLR(X) = log2(count(X) / count(RC(X)))
    // This guarantees LLR(X) = -LLR(RC(X))
    double strand_llr[4096];
    for (int i = 0; i < 4096; i++) {
        int rc = rc_hexamer_code(i);
        // Add pseudocount for smoothing
        double count_x = g_coding[i] + 1.0;
        double count_rc = g_coding[rc] + 1.0;
        strand_llr[i] = log2(count_x / count_rc);
    }

    // Verify antisymmetry
    double max_error = 0;
    for (int i = 0; i < 4096; i++) {
        int rc = rc_hexamer_code(i);
        double error = std::abs(strand_llr[i] + strand_llr[rc]);
        if (error > max_error) max_error = error;
    }
    fprintf(stderr, "Antisymmetry check: max |LLR(X) + LLR(RC(X))| = %.6f\n", max_error);

    // Write C++ header
    std::string header_path = output_prefix + "_strand_llr_antisym.hpp";
    std::ofstream hdr(header_path);
    hdr << "// Auto-generated ANTISYMMETRIC strand LLR table\n";
    hdr << "// Property: LLR(X) = -LLR(RC(X)) for all hexamers\n";
    hdr << "// Total coding hexamers: " << g_total << "\n";
    hdr << "// Usage: sum over read hexamers; positive = forward strand\n";
    hdr << "#pragma once\n\n";
    hdr << "namespace agp {\n\n";

    hdr << "inline constexpr float STRAND_LLR_ANTISYM[4096] = {\n";
    for (int i = 0; i < 4096; i++) {
        if (i % 8 == 0) hdr << "    ";
        hdr << strand_llr[i] << "f";
        if (i < 4095) hdr << ",";
        if (i % 8 == 7) hdr << "\n";
    }
    hdr << "};\n\n";

    hdr << "} // namespace agp\n";
    hdr.close();

    // Print statistics
    std::vector<std::pair<double, int>> sorted_llr;
    for (int i = 0; i < 4096; i++) {
        sorted_llr.emplace_back(strand_llr[i], i);
    }
    std::sort(sorted_llr.begin(), sorted_llr.end());

    fprintf(stderr, "\nTop 20 forward-biased hexamers:\n");
    for (int i = 4095; i >= 4076; i--) {
        int code = sorted_llr[i].second;
        int rc = rc_hexamer_code(code);
        fprintf(stderr, "  %s: LLR=%+.4f (count=%lu, rc_count=%lu)\n",
                decode_hexamer(code).c_str(), sorted_llr[i].first,
                g_coding[code], g_coding[rc]);
    }

    fprintf(stderr, "\nTop 20 reverse-biased hexamers (should be RC of above):\n");
    for (int i = 0; i < 20; i++) {
        int code = sorted_llr[i].second;
        int rc = rc_hexamer_code(code);
        fprintf(stderr, "  %s: LLR=%+.4f (count=%lu, rc_count=%lu)\n",
                decode_hexamer(code).c_str(), sorted_llr[i].first,
                g_coding[code], g_coding[rc]);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    fprintf(stderr, "\nTotal time: %.1fs\n", std::chrono::duration<double>(t2 - t0).count());
    fprintf(stderr, "Output: %s\n", header_path.c_str());

    return 0;
}
