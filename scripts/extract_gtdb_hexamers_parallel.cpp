/*
 * Extract hexamer frequencies from GTDB nucleotide CDS sequences
 * PARALLEL version using OpenMP
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <dirent.h>
#include <zlib.h>
#include <cstring>
#include <atomic>
#include <omp.h>

// Encode hexamer as integer for faster counting
inline uint32_t encode_hexamer(const char* seq) {
    static const int8_t base_map[256] = {
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
    };

    uint32_t code = 0;
    for (int i = 0; i < 6; i++) {
        int8_t b = base_map[(uint8_t)seq[i]];
        if (b < 0) return UINT32_MAX;
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

// Process a single file, return local counts
void process_file(const std::string& filepath, uint64_t* local_counts,
                  uint64_t& local_hexamers, uint64_t& local_sequences) {
    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) return;

    char buffer[65536];
    std::string current_seq;

    while (gzgets(file, buffer, sizeof(buffer)) != nullptr) {
        size_t len = strlen(buffer);
        if (len > 0 && buffer[len-1] == '\n') buffer[--len] = '\0';

        if (buffer[0] == '>') {
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
            current_seq.clear();
        } else {
            current_seq += buffer;
        }
    }

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

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <directory> <output.json> [num_threads]\n";
        return 1;
    }

    std::string dir_path = argv[1];
    std::string output_path = argv[2];
    int num_threads = (argc > 3) ? std::atoi(argv[3]) : omp_get_max_threads();

    omp_set_num_threads(num_threads);
    std::cerr << "Using " << num_threads << " threads\n";

    // List files
    DIR* dir = opendir(dir_path.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory\n";
        return 1;
    }

    std::vector<std::string> files;
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string name = entry->d_name;
        if (name.length() > 7 && name.substr(name.length() - 7) == ".fna.gz") {
            files.push_back(dir_path + "/" + name);
        }
    }
    closedir(dir);

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

        #pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < files.size(); i++) {
            process_file(files[i], local_counts, local_hexamers, local_sequences);

            uint64_t done = ++total_files;
            if (done % 5000 == 0) {
                #pragma omp critical
                std::cerr << "Processed " << done << "/" << files.size() << " files\r" << std::flush;
            }
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

    std::cerr << "\nProcessed " << total_files.load() << " files, "
              << total_sequences.load() << " sequences, "
              << total_hexamers.load() << " hexamers\n";

    // Sort and output
    std::vector<std::pair<uint32_t, uint64_t>> sorted_hexamers;
    for (uint32_t i = 0; i < 4096; i++) {
        if (global_counts[i] > 0) {
            sorted_hexamers.push_back({i, global_counts[i]});
        }
    }
    std::sort(sorted_hexamers.begin(), sorted_hexamers.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    std::ofstream out(output_path);
    out << "{\n";
    out << "  \"total_files\": " << total_files.load() << ",\n";
    out << "  \"total_sequences\": " << total_sequences.load() << ",\n";
    out << "  \"total_hexamers\": " << total_hexamers.load() << ",\n\n";
    out << "  \"hexamer_frequencies\": {\n";

    uint64_t total_hex = total_hexamers.load();
    for (size_t i = 0; i < sorted_hexamers.size(); i++) {
        std::string hex = decode_hexamer(sorted_hexamers[i].first);
        double freq = static_cast<double>(sorted_hexamers[i].second) / total_hex;
        out << "    \"" << hex << "\": " << freq;
        if (i < sorted_hexamers.size() - 1) out << ",";
        out << "\n";
    }
    out << "  },\n\n";

    out << "  \"hexamer_counts\": {\n";
    for (size_t i = 0; i < sorted_hexamers.size(); i++) {
        std::string hex = decode_hexamer(sorted_hexamers[i].first);
        out << "    \"" << hex << "\": " << sorted_hexamers[i].second;
        if (i < sorted_hexamers.size() - 1) out << ",";
        out << "\n";
    }
    out << "  }\n";
    out << "}\n";

    std::cerr << "Wrote output to " << output_path << "\n";

    // Print top 20
    std::cerr << "\n=== TOP 20 HEXAMERS ===\n";
    for (size_t i = 0; i < std::min(size_t(20), sorted_hexamers.size()); i++) {
        std::string hex = decode_hexamer(sorted_hexamers[i].first);
        double freq = static_cast<double>(sorted_hexamers[i].second) / total_hex;
        std::cerr << hex << ": " << freq << " (" << sorted_hexamers[i].second << ")\n";
    }

    return 0;
}
