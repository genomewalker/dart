// Strand hexamer extractor
// Uses DNA 6-mers (4096) instead of AA 6-mers (85M)
// Per-thread counters, no atomics, merge at end
// ~10 minutes for full GTDB

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

// Encoding: A=0, C=1, G=2, T=3
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

inline char rc_base(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

// Thread-local counters (no atomics needed)
struct ThreadCounters {
    uint64_t coding_fwd[4096] = {0};    // Forward strand coding
    uint64_t coding_rev[4096] = {0};    // Reverse strand coding (RC of fwd gene)
    uint64_t intergenic[4096] = {0};    // Between genes
    uint64_t total_coding_fwd = 0;
    uint64_t total_coding_rev = 0;
    uint64_t total_intergenic = 0;
};

// Global counters for merging
uint64_t g_coding_fwd[4096] = {0};
uint64_t g_coding_rev[4096] = {0};
uint64_t g_intergenic[4096] = {0};
uint64_t g_total_fwd = 0;
uint64_t g_total_rev = 0;
uint64_t g_total_inter = 0;

std::atomic<uint64_t> files_processed{0};
std::atomic<uint64_t> genes_processed{0};

// Parse Prodigal header: >contig_geneN # start # end # strand # ...
bool parse_prodigal_header(const std::string& header,
                           std::string& contig,
                           int64_t& start, int64_t& end, int& strand) {
    // Format: >NAME # START # END # STRAND # ATTRS
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

    // Extract last part before _N for contig name
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
    // Use zcat/pigz for decompression
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

    // Track genes per contig for intergenic extraction
    std::string last_contig;
    int64_t last_end = 0;

    auto process_gene = [&]() {
        if (current_seq.length() < 6) return;

        // Extract hexamers from coding sequence
        for (size_t i = 0; i + 5 < current_seq.length(); i++) {
            int code = encode_hexamer(current_seq.c_str() + i);
            if (code >= 0) {
                tc.coding_fwd[code]++;
                tc.total_coding_fwd++;
            }
        }

        // Extract hexamers from reverse complement (RC of forward = what reverse strand looks like)
        std::string rc_seq(current_seq.length(), 'N');
        for (size_t i = 0; i < current_seq.length(); i++) {
            rc_seq[current_seq.length() - 1 - i] = rc_base(current_seq[i]);
        }
        for (size_t i = 0; i + 5 < rc_seq.length(); i++) {
            int code = encode_hexamer(rc_seq.c_str() + i);
            if (code >= 0) {
                tc.coding_rev[code]++;
                tc.total_coding_rev++;
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
            // Process previous gene
            if (in_gene && !current_seq.empty()) {
                process_gene();
            }

            // Parse new header
            std::string header(line);
            in_gene = parse_prodigal_header(header, current_contig,
                                            current_start, current_end, current_strand);
            current_seq.clear();
        } else if (in_gene) {
            current_seq += line;
        }
    }

    // Process last gene
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

int main(int argc, char* argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s --cds <dir> --output <prefix> [--threads N]\n", argv[0]);
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

    // Collect files
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

    // Create per-thread counters
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

    // Progress thread
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
    double extract_time = std::chrono::duration<double>(t1 - t0).count();
    fprintf(stderr, "\nExtraction: %.1fs\n", extract_time);

    // Merge counters
    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i < 4096; i++) {
            g_coding_fwd[i] += thread_counters[t].coding_fwd[i];
            g_coding_rev[i] += thread_counters[t].coding_rev[i];
            g_intergenic[i] += thread_counters[t].intergenic[i];
        }
        g_total_fwd += thread_counters[t].total_coding_fwd;
        g_total_rev += thread_counters[t].total_coding_rev;
        g_total_inter += thread_counters[t].total_intergenic;
    }

    fprintf(stderr, "Total coding fwd hexamers: %lu\n", g_total_fwd);
    fprintf(stderr, "Total coding rev hexamers: %lu\n", g_total_rev);

    // Compute frequencies and strand LLR
    double fwd_freq[4096], rev_freq[4096], strand_llr[4096];
    for (int i = 0; i < 4096; i++) {
        fwd_freq[i] = (g_coding_fwd[i] + 1.0) / (g_total_fwd + 4096.0);  // Laplace smoothing
        rev_freq[i] = (g_coding_rev[i] + 1.0) / (g_total_rev + 4096.0);
        strand_llr[i] = log2(fwd_freq[i] / rev_freq[i]);
    }

    // Write C++ header
    std::string header_path = output_prefix + "_strand_hexamer.hpp";
    std::ofstream hdr(header_path);
    hdr << "// Auto-generated strand hexamer table\n";
    hdr << "// Forward coding: " << g_total_fwd << " hexamers\n";
    hdr << "// Reverse (RC): " << g_total_rev << " hexamers\n";
    hdr << "#pragma once\n\n";
    hdr << "namespace agp {\n\n";

    // Forward frequencies
    hdr << "inline constexpr float STRAND_FWD_FREQ[4096] = {\n";
    for (int i = 0; i < 4096; i++) {
        if (i % 8 == 0) hdr << "    ";
        hdr << fwd_freq[i];
        if (i < 4095) hdr << ",";
        if (i % 8 == 7) hdr << "\n";
    }
    hdr << "};\n\n";

    // Reverse frequencies
    hdr << "inline constexpr float STRAND_REV_FREQ[4096] = {\n";
    for (int i = 0; i < 4096; i++) {
        if (i % 8 == 0) hdr << "    ";
        hdr << rev_freq[i];
        if (i < 4095) hdr << ",";
        if (i % 8 == 7) hdr << "\n";
    }
    hdr << "};\n\n";

    // Strand LLR (log2(P(hex|fwd)/P(hex|rev)))
    hdr << "// Positive = forward bias, Negative = reverse bias\n";
    hdr << "inline constexpr float STRAND_LLR[4096] = {\n";
    for (int i = 0; i < 4096; i++) {
        if (i % 8 == 0) hdr << "    ";
        hdr << strand_llr[i] << "f";
        if (i < 4095) hdr << ",";
        if (i % 8 == 7) hdr << "\n";
    }
    hdr << "};\n\n";

    hdr << "} // namespace agp\n";
    hdr.close();

    // Print top discriminative hexamers
    std::vector<std::pair<double, int>> sorted_llr;
    for (int i = 0; i < 4096; i++) {
        sorted_llr.emplace_back(strand_llr[i], i);
    }
    std::sort(sorted_llr.begin(), sorted_llr.end());

    fprintf(stderr, "\nTop 20 forward-biased hexamers:\n");
    for (int i = 4095; i >= 4076; i--) {
        int code = sorted_llr[i].second;
        fprintf(stderr, "  %s: LLR=%.4f (fwd=%.2e, rev=%.2e)\n",
                decode_hexamer(code).c_str(), sorted_llr[i].first,
                fwd_freq[code], rev_freq[code]);
    }

    fprintf(stderr, "\nTop 20 reverse-biased hexamers:\n");
    for (int i = 0; i < 20; i++) {
        int code = sorted_llr[i].second;
        fprintf(stderr, "  %s: LLR=%.4f (fwd=%.2e, rev=%.2e)\n",
                decode_hexamer(code).c_str(), sorted_llr[i].first,
                fwd_freq[code], rev_freq[code]);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    fprintf(stderr, "\nTotal time: %.1fs\n", std::chrono::duration<double>(t2 - t0).count());
    fprintf(stderr, "Output: %s\n", header_path.c_str());

    return 0;
}
