/**
 * Extract codon bigram (transition) and dipeptide probabilities from CDS sequences.
 *
 * For each adjacent codon pair in CDS regions, count occurrences.
 * Output: 64x64 codon transition matrix + 20x20 dipeptide matrix
 *
 * This enables damage correction by comparing:
 *   P(observed_codon | context) vs P(undamaged_codon | context)
 *
 * Usage: extract_codon_bigrams --domain <domain> <cds_dirs...> --output-dir <dir>
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <dirent.h>
#include <zlib.h>
#include <atomic>
#include <omp.h>
#include <iomanip>

// ============== CONSTANTS ==============
constexpr size_t CACHE_LINE_SIZE = 64;
constexpr size_t MIN_CDS_LENGTH = 30;
constexpr size_t READ_BUFFER_SIZE = 1 << 16;

// ============== BASE ENCODING ==============
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

inline int encode_codon(const char* s) {
    int8_t b0 = BASE_MAP[(uint8_t)s[0]];
    int8_t b1 = BASE_MAP[(uint8_t)s[1]];
    int8_t b2 = BASE_MAP[(uint8_t)s[2]];
    if ((b0 | b1 | b2) < 0) return -1;
    return (b0 << 4) | (b1 << 2) | b2;  // 0-63
}

inline std::string decode_codon(int code) {
    static const char bases[] = "ACGT";
    std::string s(3, ' ');
    s[0] = bases[(code >> 4) & 3];
    s[1] = bases[(code >> 2) & 3];
    s[2] = bases[code & 3];
    return s;
}

// Standard genetic code
inline char translate_codon(int code) {
    static const char GENETIC_CODE[] =
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
    if (code < 0 || code >= 64) return 'X';
    return GENETIC_CODE[code];
}

// AA encoding: ACDEFGHIKLMNPQRSTVWY = 0-19
inline int aa_to_idx(char aa) {
    static const int8_t AA_MAP[256] = {
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 1, 2, 3, 4, 5, 6,-1, 7, 8, 9,10,11,-1,  // A,C,D,E,F,G,H,I,K,L,M,N
        12,13,14,15,16,-1,17,18,-1,19,-1,-1,-1,-1,-1,-1,  // P,Q,R,S,T,V,W,Y
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    };
    return AA_MAP[(uint8_t)aa];
}

inline char idx_to_aa(int idx) {
    static const char* AA_ORDER = "ACDEFGHIKLMNPQRSTVWY";
    return (idx >= 0 && idx < 20) ? AA_ORDER[idx] : '?';
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

// ============== THREAD-LOCAL STATE ==============
struct ThreadState {
    std::array<std::array<uint64_t, 64>, 64> codon_bigram{};  // [prev][curr]
    std::array<uint64_t, 64> codon_unigram{};
    std::array<std::array<uint64_t, 20>, 20> dipeptide{};  // [prev][curr]
    std::array<uint64_t, 20> aa_unigram{};

    uint64_t total_bigrams = 0;
    uint64_t total_codons = 0;
    uint64_t total_dipeptides = 0;
    uint64_t total_aa = 0;
    uint64_t sequences = 0;
    uint64_t valid_cds = 0;

    std::string seq_buffer;

    ThreadState() {
        seq_buffer.reserve(100000);
        for (auto& row : codon_bigram) row.fill(0);
        codon_unigram.fill(0);
        for (auto& row : dipeptide) row.fill(0);
        aa_unigram.fill(0);
    }
};

// ============== GLOBAL AGGREGATED STATE ==============
std::array<std::array<std::atomic<uint64_t>, 64>, 64> g_codon_bigram;
std::array<std::atomic<uint64_t>, 64> g_codon_unigram;
std::array<std::array<std::atomic<uint64_t>, 20>, 20> g_dipeptide;
std::array<std::atomic<uint64_t>, 20> g_aa_unigram;
std::atomic<uint64_t> g_total_bigrams{0};
std::atomic<uint64_t> g_total_codons{0};
std::atomic<uint64_t> g_total_dipeptides{0};
std::atomic<uint64_t> g_total_aa{0};
std::atomic<uint64_t> g_sequences{0};
std::atomic<uint64_t> g_valid_cds{0};

void aggregate_state(const ThreadState& state) {
    for (int i = 0; i < 64; i++) {
        g_codon_unigram[i] += state.codon_unigram[i];
        for (int j = 0; j < 64; j++) {
            g_codon_bigram[i][j] += state.codon_bigram[i][j];
        }
    }
    for (int i = 0; i < 20; i++) {
        g_aa_unigram[i] += state.aa_unigram[i];
        for (int j = 0; j < 20; j++) {
            g_dipeptide[i][j] += state.dipeptide[i][j];
        }
    }
    g_total_bigrams += state.total_bigrams;
    g_total_codons += state.total_codons;
    g_total_dipeptides += state.total_dipeptides;
    g_total_aa += state.total_aa;
    g_sequences += state.sequences;
    g_valid_cds += state.valid_cds;
}

// ============== FILE PROCESSING ==============
void process_file(const std::string& filepath, ThreadState& state) {
    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) return;

    gzbuffer(file, READ_BUFFER_SIZE);

    char buffer[READ_BUFFER_SIZE];
    std::string& seq = state.seq_buffer;
    seq.clear();

    auto process_sequence = [&]() {
        const size_t len = seq.length();
        if (len < MIN_CDS_LENGTH) return;

        state.sequences++;

        const char* s = seq.c_str();

        // Validate CDS
        if (!has_valid_start(s, len)) return;
        if (!has_valid_stop(s, len)) return;
        if (len % 3 != 0) return;

        state.valid_cds++;

        // Process codons
        int prev_codon = -1;
        int prev_aa = -1;

        for (size_t i = 0; i + 2 < len; i += 3) {
            int curr_codon = encode_codon(s + i);
            if (curr_codon < 0) {
                prev_codon = -1;
                prev_aa = -1;
                continue;
            }

            char aa = translate_codon(curr_codon);
            int curr_aa = aa_to_idx(aa);

            // Codon unigram
            state.codon_unigram[curr_codon]++;
            state.total_codons++;

            // Codon bigram
            if (prev_codon >= 0) {
                state.codon_bigram[prev_codon][curr_codon]++;
                state.total_bigrams++;
            }

            // AA processing (skip stops)
            if (curr_aa >= 0) {
                state.aa_unigram[curr_aa]++;
                state.total_aa++;

                // Dipeptide
                if (prev_aa >= 0) {
                    state.dipeptide[prev_aa][curr_aa]++;
                    state.total_dipeptides++;
                }
                prev_aa = curr_aa;
            } else {
                prev_aa = -1;
            }

            prev_codon = curr_codon;
        }
    };

    while (gzgets(file, buffer, sizeof(buffer)) != nullptr) {
        size_t slen = strlen(buffer);
        while (slen > 0 && (buffer[slen-1] == '\n' || buffer[slen-1] == '\r')) slen--;
        buffer[slen] = '\0';

        if (buffer[0] == '>') {
            process_sequence();
            seq.clear();
        } else {
            seq.append(buffer, slen);
        }
    }

    process_sequence();
    gzclose(file);
}

// ============== FILE DISCOVERY ==============
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

// ============== OUTPUT ==============
void output_codon_bigram_table(const std::string& prefix, const std::string& domain) {
    std::string filename = prefix + "/" + domain + "_codon_bigram.hpp";
    std::ofstream out(filename);

    std::string domain_upper = domain;
    for (char& c : domain_upper) c = std::toupper(c);

    out << "#pragma once\n";
    out << "// Auto-generated codon bigram LLR for domain: " << domain << "\n";
    out << "// Total codons: " << g_total_codons << ", bigrams: " << g_total_bigrams << "\n\n";
    out << "#include <array>\n\n";
    out << "namespace dart {\n";
    out << "namespace codon_bigram {\n\n";

    out << "inline constexpr std::array<std::array<float, 64>, 64> " << domain_upper << "_BIGRAM_LLR = {{\n";

    for (int prev = 0; prev < 64; prev++) {
        out << "    {{ // " << decode_codon(prev) << " (" << translate_codon(prev) << ")\n        ";

        uint64_t prev_total = 0;
        for (int curr = 0; curr < 64; curr++) {
            prev_total += g_codon_bigram[prev][curr].load();
        }

        for (int curr = 0; curr < 64; curr++) {
            double p_cond = (prev_total > 0) ?
                (double)g_codon_bigram[prev][curr].load() / prev_total : 0.0;
            double p_uni = (g_total_codons > 0) ?
                (double)g_codon_unigram[curr].load() / g_total_codons : 1.0/64;

            double llr = 0.0;
            if (p_cond > 0 && p_uni > 0) {
                llr = std::log2(p_cond / p_uni);
            } else if (p_cond == 0) {
                llr = -10.0;
            }

            out << std::fixed << std::setprecision(3) << llr << "f";
            if (curr < 63) out << ",";
            if ((curr + 1) % 8 == 0 && curr < 63) out << "\n        ";
        }
        out << "\n    }}";
        if (prev < 63) out << ",";
        out << "\n";
    }

    out << "}};\n\n";
    out << "} // namespace codon_bigram\n";
    out << "} // namespace dart\n";

    std::cerr << "Wrote " << filename << "\n";
}

void output_dipeptide_table(const std::string& prefix, const std::string& domain) {
    std::string filename = prefix + "/" + domain + "_dipeptide.hpp";
    std::ofstream out(filename);

    std::string domain_upper = domain;
    for (char& c : domain_upper) c = std::toupper(c);

    out << "#pragma once\n";
    out << "// Auto-generated dipeptide LLR for domain: " << domain << "\n";
    out << "// Total AA: " << g_total_aa << ", dipeptides: " << g_total_dipeptides << "\n\n";
    out << "#include <array>\n\n";
    out << "namespace dart {\n";
    out << "namespace dipeptide {\n\n";

    // AA order reference
    out << "// AA order: ACDEFGHIKLMNPQRSTVWY (indices 0-19)\n\n";

    // AA frequencies
    out << "inline constexpr std::array<float, 20> " << domain_upper << "_AA_FREQ = {{\n    ";
    for (int i = 0; i < 20; i++) {
        double freq = (g_total_aa > 0) ? (double)g_aa_unigram[i].load() / g_total_aa : 0.05;
        out << std::fixed << std::setprecision(5) << freq << "f";
        if (i < 19) out << ", ";
        if ((i + 1) % 5 == 0 && i < 19) out << "\n    ";
    }
    out << "\n}};\n\n";

    // Dipeptide LLR
    out << "inline constexpr std::array<std::array<float, 20>, 20> " << domain_upper << "_DIPEPTIDE_LLR = {{\n";

    for (int aa1 = 0; aa1 < 20; aa1++) {
        out << "    {{ // " << idx_to_aa(aa1) << " followed by:\n        ";

        for (int aa2 = 0; aa2 < 20; aa2++) {
            double p_joint = (g_total_dipeptides > 0) ?
                (double)g_dipeptide[aa1][aa2].load() / g_total_dipeptides : 0.0;
            double p_ind = (g_total_aa > 0) ?
                ((double)g_aa_unigram[aa1].load() / g_total_aa) *
                ((double)g_aa_unigram[aa2].load() / g_total_aa) : 1.0/400;

            double llr = 0.0;
            if (p_joint > 0 && p_ind > 0) {
                llr = std::log2(p_joint / p_ind);
            } else if (p_joint == 0) {
                llr = -10.0;
            }

            out << std::fixed << std::setprecision(3) << llr << "f";
            if (aa2 < 19) out << ", ";
            if ((aa2 + 1) % 5 == 0 && aa2 < 19) out << "\n        ";
        }
        out << "\n    }}";
        if (aa1 < 19) out << ",";
        out << "\n";
    }

    out << "}};\n\n";
    out << "} // namespace dipeptide\n";
    out << "} // namespace dart\n";

    std::cerr << "Wrote " << filename << "\n";
}

void print_stats() {
    std::cerr << "\n=== Extraction Statistics ===\n";
    std::cerr << "Sequences: " << g_sequences << "\n";
    std::cerr << "Valid CDS: " << g_valid_cds << "\n";
    std::cerr << "Total codons: " << g_total_codons << "\n";
    std::cerr << "Total bigrams: " << g_total_bigrams << "\n";
    std::cerr << "Total AA: " << g_total_aa << "\n";
    std::cerr << "Total dipeptides: " << g_total_dipeptides << "\n";

    // Find rare dipeptides
    std::cerr << "\n=== Rarest dipeptides (potential damage markers) ===\n";
    std::vector<std::tuple<double, int, int>> rare;

    for (int aa1 = 0; aa1 < 20; aa1++) {
        for (int aa2 = 0; aa2 < 20; aa2++) {
            double p_joint = (g_total_dipeptides > 0) ?
                (double)g_dipeptide[aa1][aa2].load() / g_total_dipeptides : 0.0;
            double p_ind = (g_total_aa > 0) ?
                ((double)g_aa_unigram[aa1].load() / g_total_aa) *
                ((double)g_aa_unigram[aa2].load() / g_total_aa) : 1.0/400;

            if (p_ind > 0) {
                double llr = (p_joint > 0) ? std::log2(p_joint / p_ind) : -10.0;
                rare.push_back({llr, aa1, aa2});
            }
        }
    }

    std::sort(rare.begin(), rare.end());

    for (int i = 0; i < std::min(15, (int)rare.size()); i++) {
        auto [llr, aa1, aa2] = rare[i];
        std::cerr << "  " << idx_to_aa(aa1) << idx_to_aa(aa2)
                  << " LLR=" << std::fixed << std::setprecision(2) << llr
                  << " (count=" << g_dipeptide[aa1][aa2].load() << ")\n";
    }
}

int main(int argc, char* argv[]) {
    std::string domain;
    std::string output_dir = ".";
    std::vector<std::string> input_dirs;
    int threads = 16;

    // Parse arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--domain" && i + 1 < argc) {
            domain = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg == "--threads" && i + 1 < argc) {
            threads = std::stoi(argv[++i]);
        } else if (arg[0] != '-') {
            input_dirs.push_back(arg);
        }
    }

    if (domain.empty() || input_dirs.empty()) {
        std::cerr << "Usage: " << argv[0] << " --domain <domain> <cds_dirs...>\n";
        std::cerr << "Options:\n";
        std::cerr << "  --output-dir <dir>  Output directory (default: .)\n";
        std::cerr << "  --threads <n>       Number of threads (default: 16)\n";
        return 1;
    }

    // Initialize global state
    for (auto& row : g_codon_bigram) for (auto& v : row) v = 0;
    for (auto& v : g_codon_unigram) v = 0;
    for (auto& row : g_dipeptide) for (auto& v : row) v = 0;
    for (auto& v : g_aa_unigram) v = 0;

    // Find files
    std::vector<std::string> files;
    for (const auto& dir : input_dirs) {
        find_fasta_files(dir, files);
    }

    std::cerr << "Found " << files.size() << " files\n";

    // Process in parallel
    omp_set_num_threads(threads);
    std::atomic<size_t> file_idx{0};
    std::atomic<size_t> files_done{0};

    #pragma omp parallel
    {
        ThreadState state;

        while (true) {
            size_t idx = file_idx.fetch_add(1);
            if (idx >= files.size()) break;

            process_file(files[idx], state);

            size_t done = ++files_done;
            if (done % 1000 == 0) {
                std::cerr << "Processed " << done << "/" << files.size()
                          << " files (" << state.total_codons << " codons)...\n";
            }
        }

        #pragma omp critical
        {
            aggregate_state(state);
        }
    }

    // Output
    print_stats();
    output_codon_bigram_table(output_dir, domain);
    output_dipeptide_table(output_dir, domain);

    return 0;
}
