/*
 * Unified Hexamer Table Extractor for AGP
 *
 * Extracts ALL hexamer tables needed by AGP in a single pass:
 * 1. Overall hexamer frequencies (for coding potential scoring)
 * 2. Positional hexamer frequencies (START/INTERNAL/END regions)
 * 3. Damage likelihood ratios (5' C→T and 3' G→A patterns)
 *
 * Supports all 8 domains:
 * - gtdb (bacteria/archaea)
 * - fungi, plant, protozoa, invertebrate, viral
 * - vertebrate_mammalian, vertebrate_other
 *
 * Usage:
 *   ./extract_unified_hexamers --domain gtdb <cds_dir1> [cds_dir2 ...] [options]
 *   ./extract_unified_hexamers --domain meta <all_cds_dirs...> [options]
 *
 * Outputs:
 *   include/dart/{domain}_hexamer_table.hpp         - Overall frequencies
 *   include/dart/{domain}_positional_hexamer.hpp    - Position-specific (START/INTERNAL/END)
 *   include/dart/{domain}_damage_likelihood.hpp     - Damage LLR tables
 *
 * Or in meta mode:
 *   include/dart/meta_hexamer_tables.hpp            - Combined weighted tables
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
#include <map>
#include <new>
#include <ctime>

// Compile-time constants
constexpr size_t CACHE_LINE_SIZE = 64;
constexpr size_t START_REGION_SIZE = 30;   // First 30bp for START hexamers
constexpr size_t END_REGION_SIZE = 30;     // Last 30bp for END hexamers
constexpr size_t MIN_CDS_LENGTH = 90;      // Need distinct regions
constexpr size_t READ_BUFFER_SIZE = 1 << 16;  // 64KB

// Hexamer encoding
alignas(CACHE_LINE_SIZE) static const int8_t BASE_MAP[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,  // A=0, C=1, G=2
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  // T=3
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,  // a=0, c=1, g=2
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  // t=3
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

// Get complement of a base
inline char complement(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        default: return 'N';
    }
}

// CDS validation
__attribute__((always_inline))
inline bool has_valid_start(const char* seq, size_t len) {
    if (len < 3) return false;
    uint8_t c0 = seq[0] & 0xDF;
    uint8_t c1 = seq[1] & 0xDF;
    uint8_t c2 = seq[2] & 0xDF;
    return (c1 == 'T' && c2 == 'G' && (c0 == 'A' || c0 == 'G' || c0 == 'T'));
}

__attribute__((always_inline))
inline bool has_valid_stop(const char* seq, size_t len) {
    if (len < 3) return false;
    const char* end = seq + len - 3;
    uint8_t c0 = end[0] & 0xDF;
    uint8_t c1 = end[1] & 0xDF;
    uint8_t c2 = end[2] & 0xDF;
    return (c0 == 'T' && ((c1 == 'A' && (c2 == 'A' || c2 == 'G')) || (c1 == 'G' && c2 == 'A')));
}

// Aligned count arrays
struct alignas(CACHE_LINE_SIZE) AlignedCounts {
    uint64_t counts[4096];
    uint64_t total;

    AlignedCounts() : total(0) {
        std::memset(counts, 0, sizeof(counts));
    }

    void merge_into(uint64_t* dest, std::atomic<uint64_t>& dest_total) const {
        for (int i = 0; i < 4096; i++) {
            dest[i] += counts[i];
        }
        dest_total += total;
    }
};

// Reverse complement
inline std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.length(), 'N');
    for (size_t i = 0; i < seq.length(); i++) {
        rc[seq.length() - 1 - i] = complement(seq[i]);
    }
    return rc;
}

// Thread-local state
struct ThreadState {
    // Overall hexamer frequencies
    AlignedCounts overall;

    // Positional hexamer frequencies
    AlignedCounts start_region;
    AlignedCounts internal_region;
    AlignedCounts end_region;

    // Sense/Antisense hexamer frequencies for strand discrimination
    AlignedCounts sense_hexamers;      // Hexamers from CDS (sense strand)
    AlignedCounts antisense_hexamers;  // Hexamers from RC of CDS (antisense strand)

    // Terminal nucleotide counts for damage LLR
    // 5' end: count C and T at first 15 positions
    // 3' end: count G and A at last 15 positions
    std::array<uint64_t, 15> t_count_5prime{};  // T at position i from 5' end
    std::array<uint64_t, 15> c_count_5prime{};  // C at position i from 5' end
    std::array<uint64_t, 15> a_count_3prime{};  // A at position i from 3' end
    std::array<uint64_t, 15> g_count_3prime{};  // G at position i from 3' end

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

    alignas(CACHE_LINE_SIZE) char buffer[READ_BUFFER_SIZE];
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

        // 1) Overall hexamer frequencies
        // Scan all hexamers at codon boundaries (frame 0)
        for (size_t i = 0; i + 5 < len; i += 3) {
            uint32_t code = encode_hexamer(seq + i);
            if (code != UINT32_MAX) {
                state.overall.counts[code]++;
                state.overall.total++;
            }
        }

        // 2) Positional hexamer frequencies
        // START region: first 30bp
        const char* ptr = seq;
        const char* start_end = seq + std::min(START_REGION_SIZE, len) - 5;
        while (ptr < start_end) {
            uint32_t code = encode_hexamer(ptr);
            if (code != UINT32_MAX) {
                state.start_region.counts[code]++;
                state.start_region.total++;
            }
            ptr += 3;
        }

        // END region: last 30bp before stop
        const size_t end_start = len - END_REGION_SIZE - 3;
        if (end_start > START_REGION_SIZE) {
            ptr = seq + end_start;
            const char* end_end = seq + len - 3 - 5;
            while (ptr < end_end) {
                uint32_t code = encode_hexamer(ptr);
                if (code != UINT32_MAX) {
                    state.end_region.counts[code]++;
                    state.end_region.total++;
                }
                ptr += 3;
            }
        }

        // INTERNAL region: between start and end
        const size_t internal_start = START_REGION_SIZE;
        const size_t internal_end = len - END_REGION_SIZE - 3;
        if (internal_end > internal_start) {
            ptr = seq + internal_start;
            const char* int_end = seq + internal_end - 5;
            while (ptr < int_end) {
                uint32_t code = encode_hexamer(ptr);
                if (code != UINT32_MAX) {
                    state.internal_region.counts[code]++;
                    state.internal_region.total++;
                }
                ptr += 3;
            }
        }

        // 3) Terminal nucleotide counts for damage LLR
        // 5' end: T and C counts at first 15 positions
        for (size_t i = 0; i < 15 && i < len; i++) {
            char c = seq[i] & 0xDF;  // Uppercase
            if (c == 'T') state.t_count_5prime[i]++;
            else if (c == 'C') state.c_count_5prime[i]++;
        }

        // 3' end: A and G counts at last 15 positions (before stop)
        size_t cds_end = len - 3;  // Exclude stop codon
        for (size_t i = 0; i < 15 && i < cds_end; i++) {
            size_t pos = cds_end - 1 - i;  // Position from end
            char c = seq[pos] & 0xDF;
            if (c == 'A') state.a_count_3prime[i]++;
            else if (c == 'G') state.g_count_3prime[i]++;
        }

        // 4) Sense/antisense hexamer frequencies
        // For strand discrimination: compare hexamers from sense vs antisense strands
        // Sense strand: the CDS as provided (coding strand)
        // Antisense strand: reverse complement of the CDS
        //
        // IMPORTANT: Extract at CODON BOUNDARIES (every 3rd position, frame 0)
        // This ensures hexamers capture the dicodon patterns that distinguish
        // sense from antisense strands. Frame-aligned hexamers have stronger
        // discriminative power because they reflect codon usage patterns.

        // Extract hexamers from sense strand at codon boundaries (frame 0)
        for (size_t i = 0; i + 5 < len; i += 3) {
            uint32_t code = encode_hexamer(seq + i);
            if (code != UINT32_MAX) {
                state.sense_hexamers.counts[code]++;
                state.sense_hexamers.total++;
            }
        }

        // Extract hexamers from antisense strand at codon boundaries
        // Note: For antisense, we also use frame 0 relative to the RC sequence
        // This means we're sampling at positions that would be frame 0 if this
        // were the coding strand (which it isn't, hence antisense patterns)
        std::string antisense = reverse_complement(current_seq);
        const char* anti_seq = antisense.c_str();
        for (size_t i = 0; i + 5 < len; i += 3) {
            uint32_t code = encode_hexamer(anti_seq + i);
            if (code != UINT32_MAX) {
                state.antisense_hexamers.counts[code]++;
                state.antisense_hexamers.total++;
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

        size_t nlen = name.length();
        
        // Check if it's a FASTA file first (most common case)
        if ((nlen > 7 && name.compare(nlen - 7, 7, ".fna.gz") == 0) ||
            (nlen > 9 && name.compare(nlen - 9, 9, ".fasta.gz") == 0) ||
            (nlen > 6 && name.compare(nlen - 6, 6, ".fa.gz") == 0)) {
            files.push_back(dir_path + "/" + name);
        } else if (entry->d_type == DT_DIR) {
            // Only recurse into directories
            find_fasta_files(dir_path + "/" + name, files);
        } else if (entry->d_type == DT_UNKNOWN) {
            // Fallback for filesystems that don't support d_type
            std::string full_path = dir_path + "/" + name;
            DIR* subdir = opendir(full_path.c_str());
            if (subdir) {
                closedir(subdir);
                find_fasta_files(full_path, files);
            }
        }
    }
    closedir(dir);
}

// Output generation
void write_frequency_array(std::ofstream& out, const std::string& name,
                           const uint64_t* counts, uint64_t total,
                           bool add_pseudocount = false) {
    out << "constexpr float " << name << "[4096] = {\n";

    for (int i = 0; i < 4096; i += 8) {
        out << "    ";
        for (int j = 0; j < 8; j++) {
            int idx = i + j;
            double freq;
            if (add_pseudocount) {
                // Add pseudocount for Laplace smoothing
                freq = (static_cast<double>(counts[idx]) + 0.1) / (total + 0.1 * 4096);
            } else {
                freq = (total > 0) ? static_cast<double>(counts[idx]) / total : 0.0;
            }
            out << std::scientific << std::setprecision(6) << freq << "f";
            if (idx < 4095) out << ", ";
        }
        out << "\n";
    }
    out << "};\n\n";
}

void write_damage_llr_array(std::ofstream& out, const std::string& name,
                            const uint64_t* counts, uint64_t total,
                            bool is_5prime) {
    // For damage LLR calculation:
    // 5' end: T hexamers are damage products, C hexamers are precursors
    //         LLR = log2(P(C-precursor) / P(T-form)) for T-starting hexamers
    // 3' end: A hexamers are damage products, G hexamers are precursors
    //         LLR = log2(P(G-precursor) / P(A-form)) for A-ending hexamers

    std::array<float, 4096> llr{};

    for (uint32_t code = 0; code < 4096; code++) {
        std::string hex = decode_hexamer(code);

        if (is_5prime) {
            // Check if hexamer starts with T (potential C→T damage product)
            if (hex[0] == 'T') {
                // Get C-precursor frequency
                std::string precursor = hex;
                precursor[0] = 'C';
                uint32_t prec_code = encode_hexamer(precursor.c_str());

                if (prec_code < 4096 && counts[code] > 0) {
                    double t_freq = static_cast<double>(counts[code]) / total;
                    double c_freq = static_cast<double>(counts[prec_code]) / total;

                    // LLR = log2(C_freq / T_freq)
                    // Positive if C is more common → consistent with damage
                    if (t_freq > 1e-9) {
                        double ratio = c_freq / t_freq;
                        if (ratio > 1.5) {  // Only informative if ratio > 1.5
                            llr[code] = std::min(12.0f, static_cast<float>(std::log2(ratio)));
                        }
                    }
                }
            }
        } else {
            // Check if hexamer ends with A (potential G→A damage product)
            if (hex[5] == 'A') {
                // Get G-precursor frequency
                std::string precursor = hex;
                precursor[5] = 'G';
                uint32_t prec_code = encode_hexamer(precursor.c_str());

                if (prec_code < 4096 && counts[code] > 0) {
                    double a_freq = static_cast<double>(counts[code]) / total;
                    double g_freq = static_cast<double>(counts[prec_code]) / total;

                    if (a_freq > 1e-9) {
                        double ratio = g_freq / a_freq;
                        if (ratio > 1.5) {
                            llr[code] = std::min(12.0f, static_cast<float>(std::log2(ratio)));
                        }
                    }
                }
            }
        }
    }

    out << "static constexpr float " << name << "[4096] = {\n";
    for (int i = 0; i < 4096; i += 16) {
        out << "    ";
        for (int j = 0; j < 16; j++) {
            int idx = i + j;
            out << std::fixed << std::setprecision(4) << llr[idx] << "f";
            if (idx < 4095) out << ", ";
        }
        out << "\n";
    }
    out << "};\n\n";
}

// Main output functions
void write_hexamer_table(const std::string& output_dir, const std::string& domain,
                         const uint64_t* counts, uint64_t total,
                         uint64_t files_processed, uint64_t valid_cds) {
    std::string upper_domain = domain;
    std::transform(upper_domain.begin(), upper_domain.end(), upper_domain.begin(), ::toupper);

    std::string filename = output_dir + "/" + domain + "_hexamer_table.hpp";
    std::ofstream out(filename);

    out << "#pragma once\n\n";
    out << "// " << upper_domain << " hexamer frequency table (AUTO-GENERATED)\n";
    out << "// Files processed: " << files_processed << "\n";
    out << "// Valid CDS: " << valid_cds << "\n";
    out << "// Total hexamers: " << total << "\n\n";
    out << "#include <cstdint>\n\n";
    out << "namespace dart {\n\n";
    out << "// Note: Use encode_hexamer() from hexamer_tables.hpp\n\n";

    write_frequency_array(out, upper_domain + "_HEXAMER_FREQ", counts, total, true);

    out << "} // namespace dart\n";

    std::cerr << "  Wrote: " << filename << "\n";
}

void write_positional_table(const std::string& output_dir, const std::string& domain,
                            const uint64_t* start_counts, uint64_t start_total,
                            const uint64_t* internal_counts, uint64_t internal_total,
                            const uint64_t* end_counts, uint64_t end_total,
                            uint64_t valid_cds) {
    std::string upper_domain = domain;
    std::transform(upper_domain.begin(), upper_domain.end(), upper_domain.begin(), ::toupper);

    std::string filename = output_dir + "/" + domain + "_positional_hexamer.hpp";
    std::ofstream out(filename);

    out << "#pragma once\n\n";
    out << "// " << upper_domain << " positional hexamer tables (AUTO-GENERATED)\n";
    out << "// Valid CDS: " << valid_cds << "\n";
    out << "// START hexamers: " << start_total << "\n";
    out << "// INTERNAL hexamers: " << internal_total << "\n";
    out << "// END hexamers: " << end_total << "\n\n";
    out << "#include <cstdint>\n\n";
    out << "namespace dart {\n\n";

    write_frequency_array(out, upper_domain + "_START_HEXAMER_FREQ", start_counts, start_total);
    write_frequency_array(out, upper_domain + "_INTERNAL_HEXAMER_FREQ", internal_counts, internal_total);
    write_frequency_array(out, upper_domain + "_END_HEXAMER_FREQ", end_counts, end_total);

    // Combined weighted table (20% start + 60% internal + 20% end)
    out << "// Combined: 20% start + 60% internal + 20% end\n";
    out << "constexpr float " << upper_domain << "_POSITIONAL_HEXAMER_FREQ[4096] = {\n";

    for (int i = 0; i < 4096; i += 8) {
        out << "    ";
        for (int j = 0; j < 8; j++) {
            int idx = i + j;
            double fs = (start_total > 0) ? static_cast<double>(start_counts[idx]) / start_total : 0.0;
            double fi = (internal_total > 0) ? static_cast<double>(internal_counts[idx]) / internal_total : 0.0;
            double fe = (end_total > 0) ? static_cast<double>(end_counts[idx]) / end_total : 0.0;
            double combined = 0.2 * fs + 0.6 * fi + 0.2 * fe;
            out << std::scientific << std::setprecision(6) << combined << "f";
            if (idx < 4095) out << ", ";
        }
        out << "\n";
    }
    out << "};\n\n";

    out << "} // namespace dart\n";

    std::cerr << "  Wrote: " << filename << "\n";
}

void write_damage_likelihood_table(const std::string& output_dir, const std::string& domain,
                                   const uint64_t* counts, uint64_t total,
                                   uint64_t valid_cds) {
    std::string upper_domain = domain;
    std::transform(upper_domain.begin(), upper_domain.end(), upper_domain.begin(), ::toupper);

    std::string filename = output_dir + "/" + domain + "_damage_likelihood.hpp";
    std::ofstream out(filename);

    out << "#pragma once\n\n";
    out << "// " << upper_domain << " damage likelihood ratios (AUTO-GENERATED)\n";
    out << "// Based on hexamer frequency analysis of " << valid_cds << " CDS sequences\n";
    out << "//\n";
    out << "// 5' C→T damage: log2(P(C-precursor) / P(T-form)) for T-starting hexamers\n";
    out << "// 3' G→A damage: log2(P(G-precursor) / P(A-form)) for A-ending hexamers\n";
    out << "//\n";
    out << "// Positive values indicate damage-consistent patterns\n\n";
    out << "#include <cstdint>\n\n";
    out << "namespace dart {\n\n";

    write_damage_llr_array(out, upper_domain + "_DAMAGE_LLR_5PRIME", counts, total, true);
    write_damage_llr_array(out, upper_domain + "_DAMAGE_LLR_3PRIME", counts, total, false);

    // Helper function for computing damage score
    out << "// Compute damage score for a sequence using " << upper_domain << " likelihood ratios\n";
    out << "inline float compute_" << domain << "_damage_score(const char* seq, size_t len) {\n";
    out << "    if (len < 6) return 0.0f;\n";
    out << "    \n";
    out << "    extern uint32_t encode_hexamer(const char*);\n";
    out << "    \n";
    out << "    float score_5prime = 0.0f;\n";
    out << "    float score_3prime = 0.0f;\n";
    out << "    \n";
    out << "    // Check 5' end (first 3 hexamers, position-weighted)\n";
    out << "    for (int i = 0; i < 3 && i + 6 <= (int)len; i++) {\n";
    out << "        uint32_t code = encode_hexamer(seq + i);\n";
    out << "        if (code < 4096) {\n";
    out << "            float weight = 1.0f / (1 << i);  // 1.0, 0.5, 0.25\n";
    out << "            score_5prime += " << upper_domain << "_DAMAGE_LLR_5PRIME[code] * weight;\n";
    out << "        }\n";
    out << "    }\n";
    out << "    \n";
    out << "    // Check 3' end (last 3 hexamers, position-weighted)\n";
    out << "    for (int i = 0; i < 3 && i + 6 <= (int)len; i++) {\n";
    out << "        size_t pos = len - 6 - i;\n";
    out << "        uint32_t code = encode_hexamer(seq + pos);\n";
    out << "        if (code < 4096) {\n";
    out << "            float weight = 1.0f / (1 << i);\n";
    out << "            score_3prime += " << upper_domain << "_DAMAGE_LLR_3PRIME[code] * weight;\n";
    out << "        }\n";
    out << "    }\n";
    out << "    \n";
    out << "    return score_5prime + score_3prime;\n";
    out << "}\n\n";

    out << "} // namespace dart\n";

    std::cerr << "  Wrote: " << filename << "\n";
}

// Strand hexamer table output
void write_strand_hexamer_table(const std::string& output_dir, const std::string& domain,
                                const uint64_t* sense_counts, uint64_t sense_total,
                                const uint64_t* antisense_counts, uint64_t antisense_total,
                                uint64_t valid_cds) {
    std::string upper_domain = domain;
    std::transform(upper_domain.begin(), upper_domain.end(), upper_domain.begin(), ::toupper);

    std::string filename = output_dir + "/" + domain + "_strand_hexamer.hpp";
    std::ofstream out(filename);

    out << "#pragma once\n\n";
    out << "// " << upper_domain << " strand-discriminative hexamer tables (AUTO-GENERATED)\n";
    out << "// Based on " << valid_cds << " CDS sequences\n";
    out << "//\n";
    out << "// Sense hexamers: extracted from coding strand (as annotated in CDS)\n";
    out << "// Antisense hexamers: extracted from reverse complement of CDS\n";
    out << "//\n";
    out << "// Use log-likelihood ratios (LLR) to determine strand orientation:\n";
    out << "//   LLR > 0: more likely sense (forward) strand\n";
    out << "//   LLR < 0: more likely antisense (reverse) strand\n\n";
    out << "#include <cstdint>\n";
    out << "#include <cmath>\n\n";
    out << "namespace dart {\n";
    out << "namespace strand {\n\n";

    // Write sense hexamer frequencies
    out << "// Sense strand hexamer frequencies (from CDS)\n";
    out << "static constexpr float " << upper_domain << "_SENSE_FREQ[4096] = {\n";
    for (int i = 0; i < 4096; i += 8) {
        out << "    ";
        for (int j = 0; j < 8; j++) {
            int idx = i + j;
            // Use Laplace smoothing with pseudocount
            double freq = (static_cast<double>(sense_counts[idx]) + 0.1) / (sense_total + 0.1 * 4096);
            out << std::scientific << std::setprecision(6) << freq << "f";
            if (idx < 4095) out << ", ";
        }
        out << "\n";
    }
    out << "};\n\n";

    // Write antisense hexamer frequencies
    out << "// Antisense strand hexamer frequencies (from RC of CDS)\n";
    out << "static constexpr float " << upper_domain << "_ANTISENSE_FREQ[4096] = {\n";
    for (int i = 0; i < 4096; i += 8) {
        out << "    ";
        for (int j = 0; j < 8; j++) {
            int idx = i + j;
            double freq = (static_cast<double>(antisense_counts[idx]) + 0.1) / (antisense_total + 0.1 * 4096);
            out << std::scientific << std::setprecision(6) << freq << "f";
            if (idx < 4095) out << ", ";
        }
        out << "\n";
    }
    out << "};\n\n";

    // Pre-compute and write log-likelihood ratios
    out << "// Pre-computed log-likelihood ratios: log2(P_sense / P_antisense)\n";
    out << "// Positive = more likely sense strand, Negative = more likely antisense\n";
    out << "static constexpr float " << upper_domain << "_STRAND_LLR[4096] = {\n";
    for (int i = 0; i < 4096; i += 8) {
        out << "    ";
        for (int j = 0; j < 8; j++) {
            int idx = i + j;
            // Use Laplace smoothing
            double sense_freq = (static_cast<double>(sense_counts[idx]) + 0.1) / (sense_total + 0.1 * 4096);
            double antisense_freq = (static_cast<double>(antisense_counts[idx]) + 0.1) / (antisense_total + 0.1 * 4096);
            double llr = std::log2(sense_freq / antisense_freq);
            // Clamp to reasonable range
            llr = std::max(-10.0, std::min(10.0, llr));
            out << std::fixed << std::setprecision(4) << llr << "f";
            if (idx < 4095) out << ", ";
        }
        out << "\n";
    }
    out << "};\n\n";

    // Write inline function to calculate strand LLR for a sequence
    out << "/**\n";
    out << " * Calculate strand log-likelihood ratio for a sequence.\n";
    out << " * \n";
    out << " * @param seq The DNA sequence\n";
    out << " * @param len Length of the sequence\n";
    out << " * @return Positive if more likely sense strand, negative if antisense\n";
    out << " */\n";
    out << "inline float calculate_" << domain << "_strand_llr(const char* seq, size_t len) {\n";
    out << "    if (len < 6) return 0.0f;\n";
    out << "    \n";
    out << "    // Hexamer encoding function\n";
    out << "    auto encode_hex = [](const char* s) -> uint32_t {\n";
    out << "        static const int8_t base_map[256] = {\n";
    out << "            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n";
    out << "            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n";
    out << "            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n";
    out << "            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n";
    out << "            -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,\n";
    out << "            -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n";
    out << "            -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,\n";
    out << "            -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1\n";
    out << "        };\n";
    out << "        uint32_t code = 0;\n";
    out << "        for (int i = 0; i < 6; i++) {\n";
    out << "            int8_t b = base_map[(uint8_t)s[i]];\n";
    out << "            if (b < 0) return UINT32_MAX;\n";
    out << "            code = (code << 2) | b;\n";
    out << "        }\n";
    out << "        return code;\n";
    out << "    };\n";
    out << "    \n";
    out << "    float total_llr = 0.0f;\n";
    out << "    int count = 0;\n";
    out << "    \n";
    out << "    // Sum LLR for all hexamers in the sequence\n";
    out << "    for (size_t i = 0; i + 5 < len; i++) {\n";
    out << "        uint32_t code = encode_hex(seq + i);\n";
    out << "        if (code < 4096) {\n";
    out << "            total_llr += " << upper_domain << "_STRAND_LLR[code];\n";
    out << "            count++;\n";
    out << "        }\n";
    out << "    }\n";
    out << "    \n";
    out << "    // Return average LLR (normalized by sequence length)\n";
    out << "    return (count > 0) ? total_llr / count : 0.0f;\n";
    out << "}\n\n";

    out << "} // namespace strand\n";
    out << "} // namespace dart\n";

    std::cerr << "  Wrote: " << filename << "\n";
}

// Manifest output
void write_manifest(const std::string& output_dir, const std::string& domain,
                    uint64_t files_processed, uint64_t total_sequences,
                    uint64_t valid_cds, uint64_t overall_hexamers,
                    uint64_t start_hexamers, uint64_t internal_hexamers,
                    uint64_t end_hexamers, uint64_t sense_hexamers, uint64_t antisense_hexamers,
                    const std::vector<std::string>& source_dirs) {
    std::string filename = output_dir + "/" + domain + "_hexamer_manifest.json";
    std::ofstream out(filename);

    // Get current date
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char date_buf[32];
    strftime(date_buf, sizeof(date_buf), "%Y-%m-%d", timeinfo);

    out << "{\n";
    out << "  \"domain\": \"" << domain << "\",\n";
    out << "  \"generated\": \"" << date_buf << "\",\n";
    out << "  \"source_directories\": [";
    for (size_t i = 0; i < source_dirs.size(); i++) {
        out << "\n    \"" << source_dirs[i] << "\"";
        if (i < source_dirs.size() - 1) out << ",";
    }
    out << "\n  ],\n";
    out << "  \"files_processed\": " << files_processed << ",\n";
    out << "  \"total_sequences\": " << total_sequences << ",\n";
    out << "  \"valid_cds\": " << valid_cds << ",\n";
    out << "  \"hexamer_counts\": {\n";
    out << "    \"overall\": " << overall_hexamers << ",\n";
    out << "    \"start\": " << start_hexamers << ",\n";
    out << "    \"internal\": " << internal_hexamers << ",\n";
    out << "    \"end\": " << end_hexamers << ",\n";
    out << "    \"sense\": " << sense_hexamers << ",\n";
    out << "    \"antisense\": " << antisense_hexamers << "\n";
    out << "  },\n";
    out << "  \"output_files\": [\n";
    out << "    \"" << domain << "_hexamer_table.hpp\",\n";
    out << "    \"" << domain << "_positional_hexamer.hpp\",\n";
    out << "    \"" << domain << "_damage_likelihood.hpp\",\n";
    out << "    \"" << domain << "_strand_hexamer.hpp\"\n";
    out << "  ]\n";
    out << "}\n";

    std::cerr << "  Wrote: " << filename << "\n";
}

// Main
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Unified Hexamer Table Extractor for AGP\n\n";
        std::cerr << "Usage: " << argv[0] << " --domain <domain> <directory1> [directory2 ...] [options]\n\n";
        std::cerr << "Domains:\n";
        std::cerr << "  gtdb                   - Bacteria/Archaea (GTDB)\n";
        std::cerr << "  fungi                  - Fungi\n";
        std::cerr << "  plant                  - Plants\n";
        std::cerr << "  protozoa               - Protozoa\n";
        std::cerr << "  invertebrate           - Invertebrates\n";
        std::cerr << "  viral                  - Viruses\n";
        std::cerr << "  vertebrate_mammalian   - Vertebrate mammals\n";
        std::cerr << "  vertebrate_other       - Other vertebrates\n";
        std::cerr << "  meta                   - Generate combined meta-domain tables\n\n";
        std::cerr << "Options:\n";
        std::cerr << "  --threads N            - Number of threads (default: all)\n";
        std::cerr << "  --output-dir DIR       - Output directory (default: include/agp)\n\n";
        std::cerr << "Outputs:\n";
        std::cerr << "  {domain}_hexamer_table.hpp        - Overall hexamer frequencies\n";
        std::cerr << "  {domain}_positional_hexamer.hpp   - Position-specific frequencies\n";
        std::cerr << "  {domain}_damage_likelihood.hpp    - Damage LLR tables\n";
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

    std::cerr << "Unified hexamer extractor\n";
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
    alignas(CACHE_LINE_SIZE) uint64_t global_overall[4096] = {0};
    alignas(CACHE_LINE_SIZE) uint64_t global_start[4096] = {0};
    alignas(CACHE_LINE_SIZE) uint64_t global_internal[4096] = {0};
    alignas(CACHE_LINE_SIZE) uint64_t global_end[4096] = {0};
    alignas(CACHE_LINE_SIZE) uint64_t global_sense[4096] = {0};
    alignas(CACHE_LINE_SIZE) uint64_t global_antisense[4096] = {0};

    std::atomic<uint64_t> total_files{0};
    std::atomic<uint64_t> total_sequences{0};
    std::atomic<uint64_t> total_valid_cds{0};
    std::atomic<uint64_t> total_overall_hex{0};
    std::atomic<uint64_t> total_start_hex{0};
    std::atomic<uint64_t> total_internal_hex{0};
    std::atomic<uint64_t> total_end_hex{0};
    std::atomic<uint64_t> total_sense_hex{0};
    std::atomic<uint64_t> total_antisense_hex{0};

    // Global terminal nucleotide counts
    std::array<std::atomic<uint64_t>, 15> global_t_5prime{};
    std::array<std::atomic<uint64_t>, 15> global_c_5prime{};
    std::array<std::atomic<uint64_t>, 15> global_a_3prime{};
    std::array<std::atomic<uint64_t>, 15> global_g_3prime{};

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
            state.overall.merge_into(global_overall, total_overall_hex);
            state.start_region.merge_into(global_start, total_start_hex);
            state.internal_region.merge_into(global_internal, total_internal_hex);
            state.end_region.merge_into(global_end, total_end_hex);
            state.sense_hexamers.merge_into(global_sense, total_sense_hex);
            state.antisense_hexamers.merge_into(global_antisense, total_antisense_hex);
            total_sequences += state.sequences;
            total_valid_cds += state.valid_cds;

            for (int j = 0; j < 15; j++) {
                global_t_5prime[j] += state.t_count_5prime[j];
                global_c_5prime[j] += state.c_count_5prime[j];
                global_a_3prime[j] += state.a_count_3prime[j];
                global_g_3prime[j] += state.g_count_3prime[j];
            }
        }
    }

    std::cerr << "\n\nResults\n";
    std::cerr << "Files processed: " << total_files.load() << "\n";
    std::cerr << "Total sequences: " << total_sequences.load() << "\n";
    std::cerr << "Valid CDS: " << total_valid_cds.load() << "\n";
    std::cerr << "Overall hexamers: " << total_overall_hex.load() << "\n";
    std::cerr << "START hexamers: " << total_start_hex.load() << "\n";
    std::cerr << "INTERNAL hexamers: " << total_internal_hex.load() << "\n";
    std::cerr << "END hexamers: " << total_end_hex.load() << "\n";
    std::cerr << "SENSE hexamers: " << total_sense_hex.load() << "\n";
    std::cerr << "ANTISENSE hexamers: " << total_antisense_hex.load() << "\n";

    if (total_overall_hex.load() == 0) {
        std::cerr << "Error: No hexamers extracted!\n";
        return 1;
    }

    // Show terminal damage bias (should be ~equal for undamaged samples)
    std::cerr << "\nTerminal nucleotide bias\n";
    std::cerr << "5' end T/(T+C) by position:\n  ";
    for (int i = 0; i < 10; i++) {
        uint64_t t = global_t_5prime[i].load();
        uint64_t c = global_c_5prime[i].load();
        if (t + c > 0) {
            std::cerr << std::fixed << std::setprecision(3)
                      << (100.0 * t / (t + c)) << "% ";
        }
    }
    std::cerr << "\n3' end A/(A+G) by position:\n  ";
    for (int i = 0; i < 10; i++) {
        uint64_t a = global_a_3prime[i].load();
        uint64_t g = global_g_3prime[i].load();
        if (a + g > 0) {
            std::cerr << std::fixed << std::setprecision(3)
                      << (100.0 * a / (a + g)) << "% ";
        }
    }
    std::cerr << "\n";

    // Write output files
    std::cerr << "\nWriting output files\n";

    write_hexamer_table(output_dir, domain,
                        global_overall, total_overall_hex.load(),
                        total_files.load(), total_valid_cds.load());

    write_positional_table(output_dir, domain,
                           global_start, total_start_hex.load(),
                           global_internal, total_internal_hex.load(),
                           global_end, total_end_hex.load(),
                           total_valid_cds.load());

    write_damage_likelihood_table(output_dir, domain,
                                  global_overall, total_overall_hex.load(),
                                  total_valid_cds.load());

    write_strand_hexamer_table(output_dir, domain,
                               global_sense, total_sense_hex.load(),
                               global_antisense, total_antisense_hex.load(),
                               total_valid_cds.load());

    write_manifest(output_dir, domain,
                   total_files.load(), total_sequences.load(),
                   total_valid_cds.load(), total_overall_hex.load(),
                   total_start_hex.load(), total_internal_hex.load(),
                   total_end_hex.load(), total_sense_hex.load(), total_antisense_hex.load(),
                   directories);

    std::cerr << "\nDone!\n";

    return 0;
}
