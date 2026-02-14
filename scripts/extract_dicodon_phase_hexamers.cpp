/*
 * Extract DICODON PHASE hexamer frequencies from CDS sequences
 * 
 * Unlike regular hexamers that only capture "coding vs non-coding",
 * dicodon phase hexamers capture FRAME INFORMATION by tracking
 * hexamers that start at codon position 0, 1, or 2.
 * 
 * Frame discrimination logic:
 * - A hexamer at true frame 0 starts at codon boundary (pos 0)
 * - A hexamer at "frame 1" starts at codon position 1
 * - A hexamer at "frame 2" starts at codon position 2
 * 
 * For frame selection:
 * - If we're testing frame X, extract hexamers at positions starting with X mod 3
 * - Compare P(hexamer | frame=0) vs P(hexamer | frame≠0)
 * 
 * This script extracts:
 * 1. PHASE0_HEXAMER_FREQ[4096] - hexamers at true codon boundaries (frame 0)
 * 2. PHASE1_HEXAMER_FREQ[4096] - hexamers offset by 1 from codon boundary  
 * 3. PHASE2_HEXAMER_FREQ[4096] - hexamers offset by 2 from codon boundary
 * 4. Codon position nucleotide frequencies for wobble analysis
 * 
 * Usage:
 *   ./extract_dicodon_phase_hexamers --domain gtdb <cds_dir> [options]
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdint>
#include <dirent.h>
#include <zlib.h>
#include <cstring>
#include <atomic>
#include <omp.h>
#include <iomanip>
#include <cctype>
#include <cmath>
#include <ctime>

// Constants
constexpr size_t CACHE_LINE_SIZE = 64;
constexpr size_t MIN_CDS_LENGTH = 30;  // Minimum length for processing
constexpr size_t READ_BUFFER_SIZE = 1 << 16;

// Hexamer encoding
alignas(CACHE_LINE_SIZE) static const int8_t BASE_MAP[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
};

__attribute__((always_inline))
inline uint32_t encode_hexamer(const char* seq) {
    int8_t b0 = BASE_MAP[(uint8_t)seq[0]];
    int8_t b1 = BASE_MAP[(uint8_t)seq[1]];
    int8_t b2 = BASE_MAP[(uint8_t)seq[2]];
    int8_t b3 = BASE_MAP[(uint8_t)seq[3]];
    int8_t b4 = BASE_MAP[(uint8_t)seq[4]];
    int8_t b5 = BASE_MAP[(uint8_t)seq[5]];
    if ((b0 | b1 | b2 | b3 | b4 | b5) < 0) return UINT32_MAX;
    return (b0 << 10) | (b1 << 8) | (b2 << 6) | (b3 << 4) | (b4 << 2) | b5;
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

// CDS validation
inline bool has_valid_start(const char* seq, size_t len) {
    if (len < 3) return false;
    uint8_t c0 = seq[0] & 0xDF;
    uint8_t c1 = seq[1] & 0xDF;
    uint8_t c2 = seq[2] & 0xDF;
    return (c1 == 'T' && c2 == 'G' && (c0 == 'A' || c0 == 'G' || c0 == 'T'));
}

inline bool has_valid_stop(const char* seq, size_t len) {
    if (len < 3) return false;
    const char* end = seq + len - 3;
    uint8_t c0 = end[0] & 0xDF;
    uint8_t c1 = end[1] & 0xDF;
    uint8_t c2 = end[2] & 0xDF;
    return (c0 == 'T' && ((c1 == 'A' && (c2 == 'A' || c2 == 'G')) || (c1 == 'G' && c2 == 'A')));
}

// Thread-local state
struct ThreadState {
    // Dicodon phase hexamer frequencies
    // Phase 0: hexamers starting at codon boundary (positions 0, 3, 6, ...)
    // Phase 1: hexamers starting at codon position 1 (positions 1, 4, 7, ...)
    // Phase 2: hexamers starting at codon position 2 (positions 2, 5, 8, ...)
    std::array<uint64_t, 4096> phase0_counts{};
    std::array<uint64_t, 4096> phase1_counts{};
    std::array<uint64_t, 4096> phase2_counts{};
    uint64_t phase0_total = 0;
    uint64_t phase1_total = 0;
    uint64_t phase2_total = 0;
    
    // Codon position nucleotide frequencies
    // Position 0: first codon position (non-synonymous)
    // Position 1: second codon position (non-synonymous)
    // Position 2: third codon position (wobble/synonymous)
    std::array<uint64_t, 4> pos0_nt{};  // A, C, G, T counts at codon pos 0
    std::array<uint64_t, 4> pos1_nt{};  // A, C, G, T counts at codon pos 1
    std::array<uint64_t, 4> pos2_nt{};  // A, C, G, T counts at codon pos 2 (wobble)
    
    uint64_t sequences = 0;
    uint64_t valid_cds = 0;
    
    std::string seq_buffer;
    
    ThreadState() {
        seq_buffer.reserve(100000);
    }
};

// File processing
void process_file(const std::string& filepath, ThreadState& state) {
    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) return;
    
    gzbuffer(file, READ_BUFFER_SIZE);
    
    char buffer[READ_BUFFER_SIZE];
    std::string& current_seq = state.seq_buffer;
    current_seq.clear();
    
    auto process_sequence = [&]() {
        const size_t len = current_seq.length();
        if (len < MIN_CDS_LENGTH) return;
        
        state.sequences++;
        
        const char* seq = current_seq.c_str();
        
        // Validate CDS
        if (!has_valid_start(seq, len)) return;
        if (!has_valid_stop(seq, len)) return;
        if (len % 3 != 0) return;
        
        state.valid_cds++;
        
        // Extract dicodon phase hexamers
        // We slide a 6bp window and track which codon phase it starts at
        for (size_t i = 0; i + 5 < len; i++) {
            uint32_t code = encode_hexamer(seq + i);
            if (code == UINT32_MAX) continue;
            
            // Determine which phase this hexamer is in
            size_t phase = i % 3;
            
            switch (phase) {
                case 0:
                    state.phase0_counts[code]++;
                    state.phase0_total++;
                    break;
                case 1:
                    state.phase1_counts[code]++;
                    state.phase1_total++;
                    break;
                case 2:
                    state.phase2_counts[code]++;
                    state.phase2_total++;
                    break;
            }
        }
        
        // Extract codon position nucleotide frequencies
        for (size_t i = 0; i + 2 < len; i += 3) {
            for (int pos = 0; pos < 3; pos++) {
                int8_t base = BASE_MAP[(uint8_t)seq[i + pos]];
                if (base >= 0) {
                    switch (pos) {
                        case 0: state.pos0_nt[base]++; break;
                        case 1: state.pos1_nt[base]++; break;
                        case 2: state.pos2_nt[base]++; break;
                    }
                }
            }
        }
    };
    
    while (gzgets(file, buffer, sizeof(buffer)) != nullptr) {
        size_t slen = strlen(buffer);
        while (slen > 0 && (buffer[slen-1] == '\n' || buffer[slen-1] == '\r')) slen--;
        buffer[slen] = '\0';
        
        if (buffer[0] == '>') {
            process_sequence();
            current_seq.clear();
        } else {
            current_seq.append(buffer, slen);
        }
    }
    
    process_sequence();
    gzclose(file);
}

// File discovery
void find_fasta_files(const std::string& dir_path, std::vector<std::string>& files) {
    DIR* dir = opendir(dir_path.c_str());
    if (!dir) return;
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string name = entry->d_name;
        if (name == "." || name == "..") continue;
        
        std::string full_path = dir_path + "/" + name;
        
        DIR* subdir = opendir(full_path.c_str());
        if (subdir) {
            closedir(subdir);
            find_fasta_files(full_path, files);
        } else {
            size_t nlen = name.length();
            if ((nlen > 7 && name.compare(nlen - 7, 7, ".fna.gz") == 0) ||
                (nlen > 9 && name.compare(nlen - 9, 9, ".fasta.gz") == 0) ||
                (nlen > 6 && name.compare(nlen - 6, 6, ".fa.gz") == 0)) {
                files.push_back(full_path);
            }
        }
    }
    closedir(dir);
}

// Output generation
void write_dicodon_phase_header(const std::string& output_dir, const std::string& domain,
                                 const std::array<uint64_t, 4096>& phase0,
                                 const std::array<uint64_t, 4096>& phase1,
                                 const std::array<uint64_t, 4096>& phase2,
                                 uint64_t total0, uint64_t total1, uint64_t total2,
                                 const std::array<uint64_t, 4>& pos0_nt,
                                 const std::array<uint64_t, 4>& pos1_nt,
                                 const std::array<uint64_t, 4>& pos2_nt,
                                 uint64_t valid_cds, uint64_t files_processed) {
    
    std::string upper_domain = domain;
    std::transform(upper_domain.begin(), upper_domain.end(), upper_domain.begin(), ::toupper);
    
    std::string filename = output_dir + "/" + domain + "_dicodon_phase.hpp";
    std::ofstream out(filename);
    
    // Get current date
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char date_buf[32];
    strftime(date_buf, sizeof(date_buf), "%Y-%m-%d", timeinfo);
    
    out << "#pragma once\n\n";
    out << "/**\n";
    out << " * @file " << domain << "_dicodon_phase.hpp\n";
    out << " * @brief Dicodon phase hexamer frequencies for " << upper_domain << "\n";
    out << " *\n";
    out << " * DICODON PHASE HEXAMERS for FRAME DISCRIMINATION\n";
    out << " *\n";
    out << " * Unlike regular hexamers that only capture 'coding vs non-coding',\n";
    out << " * dicodon phase hexamers capture FRAME INFORMATION by tracking\n";
    out << " * hexamers that start at codon position 0, 1, or 2.\n";
    out << " *\n";
    out << " * Phase 0: hexamers at true codon boundaries (frame 0)\n";
    out << " * Phase 1: hexamers offset by 1 from codon boundary\n";
    out << " * Phase 2: hexamers offset by 2 from codon boundary\n";
    out << " *\n";
    out << " * For frame selection:\n";
    out << " * - Score(frame) = sum of log(P(hex | phase=frame) / P(hex | uniform))\n";
    out << " * - Higher score = more likely to be correct frame\n";
    out << " *\n";
    out << " * Generated: " << date_buf << "\n";
    out << " * Files processed: " << files_processed << "\n";
    out << " * Valid CDS: " << valid_cds << "\n";
    out << " * Phase 0 hexamers: " << total0 << "\n";
    out << " * Phase 1 hexamers: " << total1 << "\n";
    out << " * Phase 2 hexamers: " << total2 << "\n";
    out << " */\n\n";
    
    out << "#include <cstdint>\n";
    out << "#include <cmath>\n\n";
    out << "namespace agp {\n\n";
    
    // Helper to write frequency array
    auto write_freq_array = [&](const std::string& name, const std::array<uint64_t, 4096>& counts, uint64_t total) {
        out << "constexpr float " << name << "[4096] = {\n";
        for (int i = 0; i < 4096; i += 8) {
            out << "    ";
            for (int j = 0; j < 8; j++) {
                int idx = i + j;
                // Laplace smoothing
                double freq = (static_cast<double>(counts[idx]) + 0.1) / (total + 0.1 * 4096);
                out << std::scientific << std::setprecision(6) << freq << "f";
                if (idx < 4095) out << ", ";
            }
            out << "\n";
        }
        out << "};\n\n";
    };
    
    write_freq_array(upper_domain + "_PHASE0_HEXAMER_FREQ", phase0, total0);
    write_freq_array(upper_domain + "_PHASE1_HEXAMER_FREQ", phase1, total1);
    write_freq_array(upper_domain + "_PHASE2_HEXAMER_FREQ", phase2, total2);
    
    // Write codon position nucleotide frequencies
    out << "// Codon position nucleotide frequencies\n";
    out << "// Order: A=0, C=1, G=2, T=3\n";
    
    auto write_nt_array = [&](const std::string& name, const std::array<uint64_t, 4>& counts) {
        uint64_t total = counts[0] + counts[1] + counts[2] + counts[3];
        out << "constexpr float " << name << "[4] = {\n    ";
        for (int i = 0; i < 4; i++) {
            double freq = static_cast<double>(counts[i]) / total;
            out << std::fixed << std::setprecision(6) << freq << "f";
            if (i < 3) out << ", ";
        }
        out << "\n};\n\n";
    };
    
    write_nt_array(upper_domain + "_CODON_POS0_NT_FREQ", pos0_nt);
    write_nt_array(upper_domain + "_CODON_POS1_NT_FREQ", pos1_nt);
    write_nt_array(upper_domain + "_CODON_POS2_NT_FREQ", pos2_nt);
    
    // Compute and write GC content by codon position
    double gc0 = static_cast<double>(pos0_nt[1] + pos0_nt[2]) / (pos0_nt[0] + pos0_nt[1] + pos0_nt[2] + pos0_nt[3]);
    double gc1 = static_cast<double>(pos1_nt[1] + pos1_nt[2]) / (pos1_nt[0] + pos1_nt[1] + pos1_nt[2] + pos1_nt[3]);
    double gc2 = static_cast<double>(pos2_nt[1] + pos2_nt[2]) / (pos2_nt[0] + pos2_nt[1] + pos2_nt[2] + pos2_nt[3]);
    
    out << "// GC content by codon position\n";
    out << "constexpr float " << upper_domain << "_GC1 = " << std::fixed << std::setprecision(6) << gc0 << "f;\n";
    out << "constexpr float " << upper_domain << "_GC2 = " << std::fixed << std::setprecision(6) << gc1 << "f;\n";
    out << "constexpr float " << upper_domain << "_GC3 = " << std::fixed << std::setprecision(6) << gc2 << "f;  // Wobble position\n\n";
    
    // Write scoring function
    out << "/**\n";
    out << " * Calculate frame score using dicodon phase hexamers\n";
    out << " * \n";
    out << " * @param seq DNA sequence\n";
    out << " * @param frame Reading frame (0, 1, or 2)\n";
    out << " * @return Log-likelihood ratio score (higher = more likely this frame)\n";
    out << " */\n";
    out << "inline float calculate_" << domain << "_frame_score(const char* seq, size_t len, int frame) {\n";
    out << "    if (len < 6 || frame < 0 || frame > 2) return 0.0f;\n";
    out << "    \n";
    out << "    float llr_sum = 0.0f;\n";
    out << "    int count = 0;\n";
    out << "    constexpr float UNIFORM = 1.0f / 4096.0f;\n";
    out << "    \n";
    out << "    extern uint32_t encode_hexamer(const char*);\n";
    out << "    \n";
    out << "    // Score each hexamer position\n";
    out << "    for (size_t i = 0; i + 5 < len; i++) {\n";
    out << "        uint32_t code = encode_hexamer(seq + i);\n";
    out << "        if (code >= 4096) continue;\n";
    out << "        \n";
    out << "        // What phase would this position be in if we're in 'frame'?\n";
    out << "        // phase = (i - frame + 300) % 3\n";
    out << "        int phase = (static_cast<int>(i) - frame + 300) % 3;\n";
    out << "        \n";
    out << "        float freq;\n";
    out << "        switch (phase) {\n";
    out << "            case 0: freq = " << upper_domain << "_PHASE0_HEXAMER_FREQ[code]; break;\n";
    out << "            case 1: freq = " << upper_domain << "_PHASE1_HEXAMER_FREQ[code]; break;\n";
    out << "            case 2: freq = " << upper_domain << "_PHASE2_HEXAMER_FREQ[code]; break;\n";
    out << "            default: freq = UNIFORM;\n";
    out << "        }\n";
    out << "        \n";
    out << "        if (freq > 1e-10f) {\n";
    out << "            llr_sum += std::log(freq / UNIFORM);\n";
    out << "            count++;\n";
    out << "        }\n";
    out << "    }\n";
    out << "    \n";
    out << "    return count > 0 ? llr_sum / count : 0.0f;\n";
    out << "}\n\n";
    
    // Write combined frame selection function
    out << "/**\n";
    out << " * Select best frame using dicodon phase analysis\n";
    out << " * \n";
    out << " * @param seq DNA sequence\n";
    out << " * @param len Sequence length\n";
    out << " * @return Best frame (0, 1, or 2) and score difference from second best\n";
    out << " */\n";
    out << "inline std::pair<int, float> select_" << domain << "_best_frame(const char* seq, size_t len) {\n";
    out << "    float scores[3];\n";
    out << "    scores[0] = calculate_" << domain << "_frame_score(seq, len, 0);\n";
    out << "    scores[1] = calculate_" << domain << "_frame_score(seq, len, 1);\n";
    out << "    scores[2] = calculate_" << domain << "_frame_score(seq, len, 2);\n";
    out << "    \n";
    out << "    int best = 0;\n";
    out << "    if (scores[1] > scores[best]) best = 1;\n";
    out << "    if (scores[2] > scores[best]) best = 2;\n";
    out << "    \n";
    out << "    // Find second best for confidence calculation\n";
    out << "    float second_best = (best == 0) ? scores[1] : scores[0];\n";
    out << "    for (int i = 0; i < 3; i++) {\n";
    out << "        if (i != best && scores[i] > second_best) second_best = scores[i];\n";
    out << "    }\n";
    out << "    \n";
    out << "    return {best, scores[best] - second_best};\n";
    out << "}\n\n";
    
    out << "} // namespace agp\n";
    
    std::cerr << "  Wrote: " << filename << "\n";
}

// Main
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Dicodon Phase Hexamer Extractor for AGP\n\n";
        std::cerr << "Usage: " << argv[0] << " --domain <domain> <directory1> [directory2 ...] [options]\n\n";
        std::cerr << "Domains: gtdb, fungi, plant, protozoa, invertebrate, viral,\n";
        std::cerr << "         vertebrate_mammalian, vertebrate_other\n\n";
        std::cerr << "Options:\n";
        std::cerr << "  --threads N       Number of threads (default: all)\n";
        std::cerr << "  --output-dir DIR  Output directory (default: include/agp)\n\n";
        std::cerr << "Output: {domain}_dicodon_phase.hpp\n";
        return 1;
    }
    
    // Parse arguments
    std::string domain;
    std::vector<std::string> directories;
    int num_threads = omp_get_max_threads();
    std::string output_dir = "include/agp";
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--domain" && i + 1 < argc) {
            domain = argv[++i];
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::atoi(argv[++i]);
        } else if (arg == "--output-dir" && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg[0] != '-') {
            directories.push_back(arg);
        }
    }
    
    if (domain.empty()) {
        std::cerr << "Error: --domain is required\n";
        return 1;
    }
    
    if (directories.empty()) {
        std::cerr << "Error: At least one CDS directory is required\n";
        return 1;
    }
    
    omp_set_num_threads(num_threads);
    
    std::cerr << "Dicodon phase hexamer extractor\n";
    std::cerr << "Domain: " << domain << "\n";
    std::cerr << "Threads: " << num_threads << "\n";
    std::cerr << "Output dir: " << output_dir << "\n";
    
    // Find all files
    std::vector<std::string> files;
    for (const auto& dir : directories) {
        std::cerr << "Scanning: " << dir << "\n";
        find_fasta_files(dir, files);
    }
    
    if (files.empty()) {
        std::cerr << "Error: No CDS files found\n";
        return 1;
    }
    
    std::sort(files.begin(), files.end());
    std::cerr << "Files to process: " << files.size() << "\n";
    
    // Global counts
    std::array<uint64_t, 4096> global_phase0{};
    std::array<uint64_t, 4096> global_phase1{};
    std::array<uint64_t, 4096> global_phase2{};
    std::atomic<uint64_t> total_phase0{0};
    std::atomic<uint64_t> total_phase1{0};
    std::atomic<uint64_t> total_phase2{0};
    
    std::array<uint64_t, 4> global_pos0_nt{};
    std::array<uint64_t, 4> global_pos1_nt{};
    std::array<uint64_t, 4> global_pos2_nt{};
    
    std::atomic<uint64_t> total_files{0};
    std::atomic<uint64_t> total_sequences{0};
    std::atomic<uint64_t> total_valid_cds{0};
    
    // Parallel processing
    #pragma omp parallel
    {
        ThreadState state;
        
        #pragma omp for schedule(dynamic, 4)
        for (size_t i = 0; i < files.size(); i++) {
            process_file(files[i], state);
            
            uint64_t done = ++total_files;
            if (done % 500 == 0 || done == files.size()) {
                #pragma omp critical
                std::cerr << "Progress: " << done << "/" << files.size() << "\r" << std::flush;
            }
        }
        
        // Merge thread-local to global
        #pragma omp critical
        {
            for (int j = 0; j < 4096; j++) {
                global_phase0[j] += state.phase0_counts[j];
                global_phase1[j] += state.phase1_counts[j];
                global_phase2[j] += state.phase2_counts[j];
            }
            total_phase0 += state.phase0_total;
            total_phase1 += state.phase1_total;
            total_phase2 += state.phase2_total;
            
            for (int j = 0; j < 4; j++) {
                global_pos0_nt[j] += state.pos0_nt[j];
                global_pos1_nt[j] += state.pos1_nt[j];
                global_pos2_nt[j] += state.pos2_nt[j];
            }
            
            total_sequences += state.sequences;
            total_valid_cds += state.valid_cds;
        }
    }
    
    std::cerr << "\n\nResults\n";
    std::cerr << "Files processed: " << total_files.load() << "\n";
    std::cerr << "Total sequences: " << total_sequences.load() << "\n";
    std::cerr << "Valid CDS: " << total_valid_cds.load() << "\n";
    std::cerr << "Phase 0 hexamers: " << total_phase0.load() << "\n";
    std::cerr << "Phase 1 hexamers: " << total_phase1.load() << "\n";
    std::cerr << "Phase 2 hexamers: " << total_phase2.load() << "\n";
    
    // Show codon position stats
    std::cerr << "\nCodon-position nucleotide frequencies\n";
    const char* bases = "ACGT";
    for (int pos = 0; pos < 3; pos++) {
        const auto& counts = (pos == 0) ? global_pos0_nt : (pos == 1) ? global_pos1_nt : global_pos2_nt;
        uint64_t total = counts[0] + counts[1] + counts[2] + counts[3];
        std::cerr << "Position " << (pos + 1) << ": ";
        for (int b = 0; b < 4; b++) {
            std::cerr << bases[b] << "=" << std::fixed << std::setprecision(1) 
                      << (100.0 * counts[b] / total) << "% ";
        }
        double gc = 100.0 * (counts[1] + counts[2]) / total;
        std::cerr << "(GC=" << std::fixed << std::setprecision(1) << gc << "%)\n";
    }
    
    if (total_phase0.load() == 0) {
        std::cerr << "Error: No hexamers extracted!\n";
        return 1;
    }
    
    // Write output
    std::cerr << "\nWriting output\n";
    
    write_dicodon_phase_header(output_dir, domain,
                                global_phase0, global_phase1, global_phase2,
                                total_phase0.load(), total_phase1.load(), total_phase2.load(),
                                global_pos0_nt, global_pos1_nt, global_pos2_nt,
                                total_valid_cds.load(), total_files.load());
    
    std::cerr << "\nDone!\n";
    
    return 0;
}
