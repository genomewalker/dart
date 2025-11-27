/*
 * Extract hexamer frequencies from RefSeq CDS sequences
 * Outputs C++ header file with constexpr array
 *
 * Usage: ./extract_refseq_hexamers <directory> <domain_name> [num_threads]
 * Example: ./extract_refseq_hexamers /path/to/refseq/fungi FUNGI 32
 *
 * Outputs: <domain_name>_hexamer_table.hpp
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <dirent.h>
#include <zlib.h>
#include <cstring>
#include <atomic>
#include <omp.h>
#include <iomanip>
#include <cctype>

// Encode hexamer as integer (0-4095) for fast counting
inline uint32_t encode_hexamer(const char* seq) {
    static const int8_t base_map[256] = {
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,  // A=0, C=1, G=2
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  // T=3
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,  // a=0, c=1, g=2
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1   // t=3
    };

    uint32_t code = 0;
    for (int i = 0; i < 6; i++) {
        int8_t b = base_map[(uint8_t)seq[i]];
        if (b < 0) return UINT32_MAX;  // Invalid base (N, etc.)
        code = (code << 2) | b;
    }
    return code;
}

std::string decode_hexamer(uint32_t code) {
    static const char bases[] = "ACGT";
    std::string hex(6, 'N');
    for (int i = 5; i >= 0; i--) {
        hex[i] = bases[code & 3];
        code >>= 2;
    }
    return hex;
}

// Process a single gzipped FASTA file
void process_file(const std::string& filepath, uint64_t* local_counts,
                  uint64_t& local_hexamers, uint64_t& local_sequences) {
    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) {
        std::cerr << "Warning: Cannot open " << filepath << std::endl;
        return;
    }

    char buffer[65536];
    std::string current_seq;

    while (gzgets(file, buffer, sizeof(buffer)) != nullptr) {
        size_t len = strlen(buffer);
        if (len > 0 && buffer[len-1] == '\n') buffer[--len] = '\0';
        if (len > 0 && buffer[len-1] == '\r') buffer[--len] = '\0';

        if (buffer[0] == '>') {
            // Process previous sequence
            if (current_seq.length() >= 6) {
                // Count hexamers every 3bp (codon-aligned for CDS)
                for (size_t i = 0; i + 6 <= current_seq.length(); i += 3) {
                    uint32_t code = encode_hexamer(current_seq.c_str() + i);
                    if (code != UINT32_MAX) {
                        local_counts[code]++;
                        local_hexamers++;
                    }
                }
                local_sequences++;
            }
            current_seq.clear();
        } else {
            current_seq += buffer;
        }
    }

    // Process last sequence
    if (current_seq.length() >= 6) {
        for (size_t i = 0; i + 6 <= current_seq.length(); i += 3) {
            uint32_t code = encode_hexamer(current_seq.c_str() + i);
            if (code != UINT32_MAX) {
                local_counts[code]++;
                local_hexamers++;
            }
        }
        local_sequences++;
    }

    gzclose(file);
}

// Find all matching files in directory
std::vector<std::string> find_fasta_files(const std::string& dir_path) {
    std::vector<std::string> files;

    DIR* dir = opendir(dir_path.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << dir_path << std::endl;
        return files;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string name = entry->d_name;
        // Match .fna.gz and .fasta.gz files
        if ((name.length() > 7 && name.substr(name.length() - 7) == ".fna.gz") ||
            (name.length() > 9 && name.substr(name.length() - 9) == ".fasta.gz")) {
            files.push_back(dir_path + "/" + name);
        }
    }
    closedir(dir);

    std::sort(files.begin(), files.end());
    return files;
}

// Convert domain name to uppercase
std::string to_upper(const std::string& s) {
    std::string result = s;
    for (char& c : result) c = std::toupper(c);
    return result;
}

// Convert domain name to lowercase
std::string to_lower(const std::string& s) {
    std::string result = s;
    for (char& c : result) c = std::tolower(c);
    return result;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <directory> <domain_name> [num_threads]\n";
        std::cerr << "\nExamples:\n";
        std::cerr << "  " << argv[0] << " /path/to/refseq/fungi FUNGI 32\n";
        std::cerr << "  " << argv[0] << " /path/to/refseq/viral VIRAL\n";
        std::cerr << "\nOutput: <domain_name>_hexamer_table.hpp\n";
        return 1;
    }

    std::string dir_path = argv[1];
    std::string domain_name = argv[2];
    std::string domain_upper = to_upper(domain_name);
    std::string domain_lower = to_lower(domain_name);

    int num_threads = (argc > 3) ? std::atoi(argv[3]) : omp_get_max_threads();
    omp_set_num_threads(num_threads);

    std::cerr << "Domain: " << domain_name << std::endl;
    std::cerr << "Using " << num_threads << " threads\n";

    // Find all FASTA files
    std::vector<std::string> files = find_fasta_files(dir_path);
    if (files.empty()) {
        std::cerr << "Error: No .fna.gz or .fasta.gz files found in " << dir_path << std::endl;
        return 1;
    }
    std::cerr << "Found " << files.size() << " files to process\n";

    // Global counts
    uint64_t global_counts[4096] = {0};
    std::atomic<uint64_t> total_files{0};
    std::atomic<uint64_t> total_sequences{0};
    std::atomic<uint64_t> total_hexamers{0};

    // Parallel processing
    #pragma omp parallel
    {
        // Thread-local counts
        uint64_t local_counts[4096] = {0};
        uint64_t local_sequences = 0;
        uint64_t local_hexamers = 0;

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < files.size(); i++) {
            process_file(files[i], local_counts, local_hexamers, local_sequences);

            uint64_t done = ++total_files;
            #pragma omp critical
            std::cerr << "Processed " << done << "/" << files.size()
                      << " files (" << files[i].substr(files[i].rfind('/') + 1) << ")\r" << std::flush;
        }

        // Merge local counts to global
        #pragma omp critical
        {
            for (int i = 0; i < 4096; i++) {
                global_counts[i] += local_counts[i];
            }
            total_sequences += local_sequences;
            total_hexamers += local_hexamers;
        }
    }

    std::cerr << "\n\nProcessed " << total_files.load() << " files, "
              << total_sequences.load() << " sequences, "
              << total_hexamers.load() << " hexamers\n";

    // Calculate frequencies
    uint64_t total_hex = total_hexamers.load();
    if (total_hex == 0) {
        std::cerr << "Error: No hexamers extracted!\n";
        return 1;
    }

    // Output C++ header file
    std::string output_path = domain_lower + "_hexamer_table.hpp";
    std::ofstream out(output_path);

    out << "#pragma once\n\n";
    out << "// " << domain_upper << " hexamer frequency table\n";
    out << "// Generated from RefSeq " << domain_lower << " CDS sequences\n";
    out << "// Total sequences: " << total_sequences.load() << "\n";
    out << "// Total hexamers: " << total_hexamers.load() << "\n";
    out << "// Files processed: " << total_files.load() << "\n";
    out << "//\n";
    out << "// Usage: float freq = " << domain_upper << "_HEXAMER_FREQ[encode_hexamer(\"ATGATG\")];\n\n";

    out << "#include <cstdint>\n\n";
    out << "namespace agp {\n\n";

    // Write frequency array
    out << "constexpr float " << domain_upper << "_HEXAMER_FREQ[4096] = {\n";

    for (int i = 0; i < 4096; i += 8) {
        out << "    ";
        for (int j = 0; j < 8; j++) {
            int idx = i + j;
            double freq = static_cast<double>(global_counts[idx]) / total_hex;
            out << std::scientific << std::setprecision(6) << freq << "f";
            if (idx < 4095) out << ", ";
        }
        out << "\n";
    }
    out << "};\n\n";

    out << "} // namespace agp\n";

    std::cerr << "Wrote output to " << output_path << "\n";

    // Also output JSON for reference
    std::string json_path = domain_lower + "_hexamers.json";
    std::ofstream json_out(json_path);

    // Sort hexamers by frequency for JSON
    std::vector<std::pair<uint32_t, uint64_t>> sorted_hexamers;
    for (uint32_t i = 0; i < 4096; i++) {
        sorted_hexamers.push_back({i, global_counts[i]});
    }
    std::sort(sorted_hexamers.begin(), sorted_hexamers.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    json_out << "{\n";
    json_out << "  \"domain\": \"" << domain_name << "\",\n";
    json_out << "  \"total_files\": " << total_files.load() << ",\n";
    json_out << "  \"total_sequences\": " << total_sequences.load() << ",\n";
    json_out << "  \"total_hexamers\": " << total_hexamers.load() << ",\n\n";
    json_out << "  \"top_20_hexamers\": {\n";

    for (size_t i = 0; i < std::min(size_t(20), sorted_hexamers.size()); i++) {
        std::string hex = decode_hexamer(sorted_hexamers[i].first);
        double freq = static_cast<double>(sorted_hexamers[i].second) / total_hex;
        json_out << "    \"" << hex << "\": " << std::scientific << freq;
        if (i < 19) json_out << ",";
        json_out << "\n";
    }
    json_out << "  }\n}\n";

    std::cerr << "Wrote JSON to " << json_path << "\n";

    // Print top 20 hexamers
    std::cerr << "\n=== TOP 20 HEXAMERS (" << domain_upper << ") ===\n";
    for (size_t i = 0; i < std::min(size_t(20), sorted_hexamers.size()); i++) {
        std::string hex = decode_hexamer(sorted_hexamers[i].first);
        double freq = static_cast<double>(sorted_hexamers[i].second) / total_hex;
        std::cerr << hex << ": " << std::fixed << std::setprecision(6) << freq
                  << " (" << sorted_hexamers[i].second << ")\n";
    }

    return 0;
}
