// Extract coding vs intergenic hexamers for strand discrimination
// Uses CDS coordinates to identify intergenic regions in genomes
//
// Key insight: Intergenic regions have different hexamer patterns than coding
// The transition zones (coding→intergenic, intergenic→coding) are asymmetric

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <atomic>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <chrono>

namespace fs = std::filesystem;

struct GeneCoord {
    int64_t start;
    int64_t end;
    int strand;  // 1 = forward, -1 = reverse
};

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

struct ThreadCounters {
    uint64_t coding[4096] = {0};
    uint64_t intergenic[4096] = {0};
    uint64_t boundary_5p[4096] = {0};  // 30bp upstream of gene start
    uint64_t boundary_3p[4096] = {0};  // 30bp downstream of gene end
    uint64_t total_coding = 0;
    uint64_t total_intergenic = 0;
    uint64_t total_5p = 0;
    uint64_t total_3p = 0;
};

std::atomic<uint64_t> files_processed{0};
std::atomic<uint64_t> genes_processed{0};
std::atomic<uint64_t> intergenic_regions{0};

// Parse Prodigal header: >contig_geneN # start # end # strand # ...
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

void process_genome_pair(const std::string& cds_path, const std::string& genome_path,
                         ThreadCounters& tc) {
    // Load gene coordinates from CDS file
    std::map<std::string, std::vector<GeneCoord>> genes_by_contig;

    std::string cmd = "pigz -dc '" + cds_path + "' 2>/dev/null || zcat '" + cds_path + "' 2>/dev/null";
    FILE* fp = popen(cmd.c_str(), "r");
    if (!fp) return;

    char line[65536];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] != '>') continue;

        std::string contig;
        int64_t start, end;
        int strand;
        if (parse_prodigal_header(line, contig, start, end, strand)) {
            genes_by_contig[contig].push_back({start, end, strand});
            genes_processed++;
        }
    }
    pclose(fp);

    // Sort genes by position for each contig
    for (auto& [contig, genes] : genes_by_contig) {
        std::sort(genes.begin(), genes.end(),
                  [](const GeneCoord& a, const GeneCoord& b) { return a.start < b.start; });
    }

    // Load genome and extract hexamers
    cmd = "pigz -dc '" + genome_path + "' 2>/dev/null || zcat '" + genome_path + "' 2>/dev/null";
    fp = popen(cmd.c_str(), "r");
    if (!fp) return;

    std::string current_contig;
    std::string current_seq;

    auto process_contig = [&]() {
        if (current_seq.empty()) return;

        auto it = genes_by_contig.find(current_contig);
        if (it == genes_by_contig.end()) return;

        const auto& genes = it->second;
        if (genes.empty()) return;

        // Mark regions as coding or intergenic
        std::vector<bool> is_coding(current_seq.length(), false);
        for (const auto& g : genes) {
            int64_t s = (g.start - 1 > 0) ? (g.start - 1) : 0;  // 1-based to 0-based
            int64_t e = ((int64_t)current_seq.length() < g.end) ? (int64_t)current_seq.length() : g.end;
            for (int64_t i = s; i < e; i++) {
                is_coding[i] = true;
            }
        }

        // Extract hexamers from intergenic regions
        for (size_t i = 0; i + 5 < current_seq.length(); i++) {
            // Check if this hexamer is fully in intergenic
            bool all_intergenic = true;
            for (size_t j = i; j < i + 6; j++) {
                if (is_coding[j]) {
                    all_intergenic = false;
                    break;
                }
            }

            if (all_intergenic) {
                int code = encode_hexamer(current_seq.c_str() + i);
                if (code >= 0) {
                    tc.intergenic[code]++;
                    tc.total_intergenic++;
                }
            }
        }

        intergenic_regions++;

        // Extract boundary hexamers (30bp around gene boundaries)
        for (const auto& g : genes) {
            int64_t gene_start = g.start - 1;  // 0-based
            int64_t gene_end = g.end;

            // 5' boundary: 30bp upstream of gene start
            if (gene_start >= 30) {
                std::string boundary = current_seq.substr(gene_start - 30, 30);
                extract_hexamers(boundary, tc.boundary_5p, tc.total_5p);
            }

            // 3' boundary: 30bp downstream of gene end
            if (gene_end + 30 <= (int64_t)current_seq.length()) {
                std::string boundary = current_seq.substr(gene_end, 30);
                extract_hexamers(boundary, tc.boundary_3p, tc.total_3p);
            }
        }
    };

    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
            line[--len] = '\0';
        }

        if (line[0] == '>') {
            process_contig();
            current_contig = std::string(line + 1);
            // Extract just the contig ID (first word)
            size_t space = current_contig.find(' ');
            if (space != std::string::npos) {
                current_contig = current_contig.substr(0, space);
            }
            current_seq.clear();
        } else {
            current_seq += line;
        }
    }
    process_contig();

    pclose(fp);
    files_processed++;
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
    if (argc < 7) {
        fprintf(stderr, "Usage: %s --cds <dir> --genomes <dir> --output <prefix> [--threads N]\n", argv[0]);
        return 1;
    }

    std::string cds_dir, genome_dir, output_prefix;
    int num_threads = std::thread::hardware_concurrency();

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--cds") == 0 && i + 1 < argc) {
            cds_dir = argv[++i];
        } else if (strcmp(argv[i], "--genomes") == 0 && i + 1 < argc) {
            genome_dir = argv[++i];
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        }
    }

    // Find matching CDS and genome files
    std::vector<std::pair<std::string, std::string>> file_pairs;

    for (const auto& entry : fs::directory_iterator(cds_dir)) {
        if (!entry.is_regular_file()) continue;
        std::string cds_path = entry.path().string();

        // Extract genome ID from CDS filename (e.g., RS_GCF_039595165.1_protein.fna.gz -> GCF_039595165.1)
        std::string filename = entry.path().filename().string();
        size_t underscore = filename.find('_');
        if (underscore == std::string::npos) continue;

        std::string genome_id = filename.substr(underscore + 1);
        size_t protein = genome_id.find("_protein");
        if (protein != std::string::npos) {
            genome_id = genome_id.substr(0, protein);
        }

        // Find matching genome file
        std::string genome_path = genome_dir + "/" + genome_id + "_genomic.fna.gz";
        if (fs::exists(genome_path)) {
            file_pairs.emplace_back(cds_path, genome_path);
        }
    }

    fprintf(stderr, "Found %zu genome/CDS pairs, using %d threads\n", file_pairs.size(), num_threads);

    if (file_pairs.empty()) {
        fprintf(stderr, "No matching files found!\n");
        return 1;
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    // Process in parallel
    std::vector<ThreadCounters> thread_counters(num_threads);
    std::vector<std::thread> threads;

    size_t pairs_per_thread = (file_pairs.size() + num_threads - 1) / num_threads;

    for (int t = 0; t < num_threads; t++) {
        size_t start = t * pairs_per_thread;
        size_t end = std::min(start + pairs_per_thread, file_pairs.size());

        threads.emplace_back([&file_pairs, &thread_counters, start, end, t]() {
            for (size_t i = start; i < end; i++) {
                process_genome_pair(file_pairs[i].first, file_pairs[i].second,
                                    thread_counters[t]);
            }
        });
    }

    // Progress thread
    std::atomic<bool> done{false};
    std::thread progress_thread([&]() {
        while (!done) {
            std::this_thread::sleep_for(std::chrono::seconds(10));
            fprintf(stderr, "\r[%lus] Files: %lu, Genes: %lu, Intergenic: %lu",
                    (unsigned long)std::chrono::duration_cast<std::chrono::seconds>(
                        std::chrono::high_resolution_clock::now() - t0).count(),
                    files_processed.load(), genes_processed.load(), intergenic_regions.load());
        }
    });

    for (auto& t : threads) t.join();
    done = true;
    progress_thread.join();

    // Merge counters
    uint64_t g_coding[4096] = {0}, g_intergenic[4096] = {0};
    uint64_t g_boundary_5p[4096] = {0}, g_boundary_3p[4096] = {0};
    uint64_t g_total_coding = 0, g_total_intergenic = 0;
    uint64_t g_total_5p = 0, g_total_3p = 0;

    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i < 4096; i++) {
            g_intergenic[i] += thread_counters[t].intergenic[i];
            g_boundary_5p[i] += thread_counters[t].boundary_5p[i];
            g_boundary_3p[i] += thread_counters[t].boundary_3p[i];
        }
        g_total_intergenic += thread_counters[t].total_intergenic;
        g_total_5p += thread_counters[t].total_5p;
        g_total_3p += thread_counters[t].total_3p;
    }

    fprintf(stderr, "\nTotal intergenic hexamers: %lu\n", g_total_intergenic);
    fprintf(stderr, "Total 5' boundary hexamers: %lu\n", g_total_5p);
    fprintf(stderr, "Total 3' boundary hexamers: %lu\n", g_total_3p);

    // Load coding hexamers from existing table
    // For now, output just the intergenic tables

    // Write header
    std::string header_path = output_prefix + "_intergenic_hexamer.hpp";
    std::ofstream hdr(header_path);
    hdr << "// Auto-generated intergenic hexamer table\n";
    hdr << "// Intergenic: " << g_total_intergenic << " hexamers\n";
    hdr << "// 5' boundary (upstream of ATG): " << g_total_5p << " hexamers\n";
    hdr << "// 3' boundary (downstream of stop): " << g_total_3p << " hexamers\n";
    hdr << "#pragma once\n\n";
    hdr << "namespace agp {\n\n";

    // Compute frequencies
    auto write_freq = [&](const char* name, uint64_t* counts, uint64_t total) {
        hdr << "inline constexpr float " << name << "[4096] = {\n";
        for (int i = 0; i < 4096; i++) {
            if (i % 8 == 0) hdr << "    ";
            double freq = (counts[i] + 1.0) / (total + 4096.0);
            char buf[32];
            snprintf(buf, sizeof(buf), "%.6ef", freq);
            hdr << buf;
            if (i < 4095) hdr << ",";
            if (i % 8 == 7) hdr << "\n";
        }
        hdr << "};\n\n";
    };

    write_freq("INTERGENIC_FREQ", g_intergenic, g_total_intergenic);
    write_freq("BOUNDARY_5P_FREQ", g_boundary_5p, g_total_5p);
    write_freq("BOUNDARY_3P_FREQ", g_boundary_3p, g_total_3p);

    hdr << "} // namespace agp\n";
    hdr.close();

    fprintf(stderr, "Output: %s\n", header_path.c_str());

    return 0;
}
