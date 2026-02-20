// Boundary-aware strand hexamer extractor
// Extracts hexamers from:
//   - Start region (first 30bp, includes ATG)
//   - Stop region (last 30bp, includes stop codon)
//   - Interior (middle of gene, > 30bp from ends)
// Also extracts RC versions for strand discrimination
//
// Key insight: Start patterns are at 5' end of coding, stop patterns at 3' end
// When read in reverse strand, patterns are flipped (RC(stop)...RC(start))

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

constexpr int BOUNDARY_SIZE = 30;  // bp from start/stop to sample

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

inline std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.length(), 'N');
    for (size_t i = 0; i < seq.length(); i++) {
        rc[seq.length() - 1 - i] = rc_base(seq[i]);
    }
    return rc;
}

struct ThreadCounters {
    // Forward strand patterns (coding direction)
    uint64_t start_fwd[4096] = {0};     // Near ATG (5' end of gene)
    uint64_t stop_fwd[4096] = {0};      // Near stop codon (3' end of gene)
    uint64_t interior_fwd[4096] = {0};  // Middle of gene

    // Reverse complement patterns (what we see reading wrong strand)
    uint64_t start_rc[4096] = {0};      // RC of start region
    uint64_t stop_rc[4096] = {0};       // RC of stop region
    uint64_t interior_rc[4096] = {0};   // RC of interior

    uint64_t total_start_fwd = 0, total_stop_fwd = 0, total_interior_fwd = 0;
    uint64_t total_start_rc = 0, total_stop_rc = 0, total_interior_rc = 0;
};

std::atomic<uint64_t> files_processed{0};
std::atomic<uint64_t> genes_processed{0};

void extract_hexamers(const std::string& seq, uint64_t* counts, uint64_t& total) {
    if (seq.length() < 6) return;
    for (size_t i = 0; i + 5 < seq.length(); i++) {
        int code = encode_hexamer(seq.c_str() + i);
        if (code >= 0) {
            counts[code]++;
            total++;
        }
    }
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
    bool in_seq = false;

    auto process_gene = [&]() {
        if (current_seq.length() < 60) return;  // Need at least 60bp for meaningful regions

        // Extract regions
        std::string start_region = current_seq.substr(0, BOUNDARY_SIZE);
        std::string stop_region = current_seq.substr(current_seq.length() - BOUNDARY_SIZE);

        // Interior: middle, avoiding boundaries
        size_t interior_start = BOUNDARY_SIZE;
        size_t interior_end = current_seq.length() - BOUNDARY_SIZE;
        std::string interior_region;
        if (interior_end > interior_start) {
            interior_region = current_seq.substr(interior_start, interior_end - interior_start);
        }

        // Forward strand patterns
        extract_hexamers(start_region, tc.start_fwd, tc.total_start_fwd);
        extract_hexamers(stop_region, tc.stop_fwd, tc.total_stop_fwd);
        if (!interior_region.empty()) {
            extract_hexamers(interior_region, tc.interior_fwd, tc.total_interior_fwd);
        }

        // Reverse complement patterns (what we'd see reading reverse strand)
        std::string start_rc = reverse_complement(start_region);
        std::string stop_rc = reverse_complement(stop_region);
        extract_hexamers(start_rc, tc.start_rc, tc.total_start_rc);
        extract_hexamers(stop_rc, tc.stop_rc, tc.total_stop_rc);
        if (!interior_region.empty()) {
            std::string interior_rc = reverse_complement(interior_region);
            extract_hexamers(interior_rc, tc.interior_rc, tc.total_interior_rc);
        }

        genes_processed++;
    };

    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
            line[--len] = '\0';
        }

        if (line[0] == '>') {
            if (in_seq && !current_seq.empty()) {
                process_gene();
            }
            in_seq = true;
            current_seq.clear();
        } else if (in_seq) {
            current_seq += line;
        }
    }

    if (in_seq && !current_seq.empty()) {
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
    uint64_t g_start_fwd[4096] = {0}, g_stop_fwd[4096] = {0}, g_interior_fwd[4096] = {0};
    uint64_t g_start_rc[4096] = {0}, g_stop_rc[4096] = {0}, g_interior_rc[4096] = {0};
    uint64_t g_total_start_fwd = 0, g_total_stop_fwd = 0, g_total_interior_fwd = 0;
    uint64_t g_total_start_rc = 0, g_total_stop_rc = 0, g_total_interior_rc = 0;

    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i < 4096; i++) {
            g_start_fwd[i] += thread_counters[t].start_fwd[i];
            g_stop_fwd[i] += thread_counters[t].stop_fwd[i];
            g_interior_fwd[i] += thread_counters[t].interior_fwd[i];
            g_start_rc[i] += thread_counters[t].start_rc[i];
            g_stop_rc[i] += thread_counters[t].stop_rc[i];
            g_interior_rc[i] += thread_counters[t].interior_rc[i];
        }
        g_total_start_fwd += thread_counters[t].total_start_fwd;
        g_total_stop_fwd += thread_counters[t].total_stop_fwd;
        g_total_interior_fwd += thread_counters[t].total_interior_fwd;
        g_total_start_rc += thread_counters[t].total_start_rc;
        g_total_stop_rc += thread_counters[t].total_stop_rc;
        g_total_interior_rc += thread_counters[t].total_interior_rc;
    }

    fprintf(stderr, "Total hexamers:\n");
    fprintf(stderr, "  Start (fwd): %lu, Stop (fwd): %lu, Interior (fwd): %lu\n",
            g_total_start_fwd, g_total_stop_fwd, g_total_interior_fwd);
    fprintf(stderr, "  Start (rc): %lu, Stop (rc): %lu, Interior (rc): %lu\n",
            g_total_start_rc, g_total_stop_rc, g_total_interior_rc);

    // Compute frequencies with Laplace smoothing
    auto compute_freq = [](uint64_t* counts, uint64_t total, double* freq) {
        for (int i = 0; i < 4096; i++) {
            freq[i] = (counts[i] + 1.0) / (total + 4096.0);
        }
    };

    double freq_start_fwd[4096], freq_stop_fwd[4096], freq_interior_fwd[4096];
    double freq_start_rc[4096], freq_stop_rc[4096], freq_interior_rc[4096];

    compute_freq(g_start_fwd, g_total_start_fwd, freq_start_fwd);
    compute_freq(g_stop_fwd, g_total_stop_fwd, freq_stop_fwd);
    compute_freq(g_interior_fwd, g_total_interior_fwd, freq_interior_fwd);
    compute_freq(g_start_rc, g_total_start_rc, freq_start_rc);
    compute_freq(g_stop_rc, g_total_stop_rc, freq_stop_rc);
    compute_freq(g_interior_rc, g_total_interior_rc, freq_interior_rc);

    // Compute strand LLR: log2(P(hex|fwd) / P(hex|rc))
    // For strand discrimination: fwd = start_fwd + interior_fwd + stop_fwd
    //                            rc  = start_rc + interior_rc + stop_rc
    // But we want position-specific signals too

    double strand_llr_start[4096], strand_llr_stop[4096], strand_llr_interior[4096];
    for (int i = 0; i < 4096; i++) {
        strand_llr_start[i] = log2(freq_start_fwd[i] / freq_start_rc[i]);
        strand_llr_stop[i] = log2(freq_stop_fwd[i] / freq_stop_rc[i]);
        strand_llr_interior[i] = log2(freq_interior_fwd[i] / freq_interior_rc[i]);
    }

    // Write C++ header
    std::string header_path = output_prefix + "_strand_hexamer.hpp";
    std::ofstream hdr(header_path);
    hdr << "// Auto-generated boundary-aware strand hexamer table\n";
    hdr << "// Start region (near ATG): " << g_total_start_fwd << " hexamers\n";
    hdr << "// Stop region (near stop): " << g_total_stop_fwd << " hexamers\n";
    hdr << "// Interior region: " << g_total_interior_fwd << " hexamers\n";
    hdr << "#pragma once\n\n";
    hdr << "namespace dart {\n\n";

    // Write all frequency tables
    auto write_table = [&](const char* name, double* freq) {
        hdr << "inline constexpr float " << name << "[4096] = {\n";
        for (int i = 0; i < 4096; i++) {
            if (i % 8 == 0) hdr << "    ";
            hdr << freq[i];
            if (i < 4095) hdr << ",";
            if (i % 8 == 7) hdr << "\n";
        }
        hdr << "};\n\n";
    };

    write_table("STRAND_START_FWD_FREQ", freq_start_fwd);
    write_table("STRAND_START_RC_FREQ", freq_start_rc);
    write_table("STRAND_STOP_FWD_FREQ", freq_stop_fwd);
    write_table("STRAND_STOP_RC_FREQ", freq_stop_rc);
    write_table("STRAND_INTERIOR_FWD_FREQ", freq_interior_fwd);
    write_table("STRAND_INTERIOR_RC_FREQ", freq_interior_rc);

    // Write LLR tables (use fixed format to avoid "0f" parsing issue)
    auto write_llr_table = [&](const char* name, double* llr) {
        hdr << "// Positive = forward bias, Negative = reverse bias\n";
        hdr << "inline constexpr float " << name << "[4096] = {\n";
        for (int i = 0; i < 4096; i++) {
            if (i % 8 == 0) hdr << "    ";
            char buf[32];
            snprintf(buf, sizeof(buf), "%.6ff", llr[i]);
            hdr << buf;
            if (i < 4095) hdr << ",";
            if (i % 8 == 7) hdr << "\n";
        }
        hdr << "};\n\n";
    };

    write_llr_table("STRAND_LLR_START", strand_llr_start);
    write_llr_table("STRAND_LLR_STOP", strand_llr_stop);
    write_llr_table("STRAND_LLR_INTERIOR", strand_llr_interior);

    // Combined LLR (weighted average)
    hdr << "// Combined strand LLR (start*0.3 + interior*0.4 + stop*0.3)\n";
    hdr << "inline constexpr float STRAND_LLR_COMBINED[4096] = {\n";
    for (int i = 0; i < 4096; i++) {
        if (i % 8 == 0) hdr << "    ";
        double combined = 0.3 * strand_llr_start[i] + 0.4 * strand_llr_interior[i] + 0.3 * strand_llr_stop[i];
        char buf[32];
        snprintf(buf, sizeof(buf), "%.6ff", combined);
        hdr << buf;
        if (i < 4095) hdr << ",";
        if (i % 8 == 7) hdr << "\n";
    }
    hdr << "};\n\n";

    hdr << "} // namespace dart\n";
    hdr.close();

    // Print top discriminative hexamers for each region
    auto print_top = [&](const char* name, double* llr) {
        std::vector<std::pair<double, int>> sorted;
        for (int i = 0; i < 4096; i++) {
            sorted.emplace_back(llr[i], i);
        }
        std::sort(sorted.begin(), sorted.end());

        fprintf(stderr, "\nTop 10 forward-biased %s hexamers:\n", name);
        for (int i = 4095; i >= 4086; i--) {
            int code = sorted[i].second;
            fprintf(stderr, "  %s: LLR=%.4f\n", decode_hexamer(code).c_str(), sorted[i].first);
        }

        fprintf(stderr, "Top 10 reverse-biased %s hexamers:\n", name);
        for (int i = 0; i < 10; i++) {
            int code = sorted[i].second;
            fprintf(stderr, "  %s: LLR=%.4f\n", decode_hexamer(code).c_str(), sorted[i].first);
        }
    };

    print_top("START", strand_llr_start);
    print_top("STOP", strand_llr_stop);
    print_top("INTERIOR", strand_llr_interior);

    auto t2 = std::chrono::high_resolution_clock::now();
    fprintf(stderr, "\nTotal time: %.1fs\n", std::chrono::duration<double>(t2 - t0).count());
    fprintf(stderr, "Output: %s\n", header_path.c_str());

    return 0;
}
