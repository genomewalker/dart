/**
 * Extract damage-correction context patterns from GTDB CDS sequences.
 *
 * Extracts three types of patterns for ancient DNA damage correction:
 * 1. N-terminal AA positional frequencies (positions 1-10 after start Met)
 * 2. Codon-level bigram context (CGG vs TGG preceding codons)
 * 3. Rare W-containing tripeptides (forbidden patterns)
 *
 * These patterns enable high-confidence correction of C→T damage:
 * - CGG (Arg) → TGG (Trp) at 5' end
 * - CGC/CGT (Arg) → TGC/TGT (Cys) at 5' end
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <dirent.h>
#include <zlib.h>
#include <atomic>
#include <omp.h>
#include <iomanip>

constexpr size_t READ_BUFFER_SIZE = 1 << 16;
constexpr size_t MIN_CDS_LENGTH = 30;
constexpr int N_TERMINAL_POSITIONS = 10;  // Track first 10 AAs after Met
constexpr int NUM_DISTANCE_BINS = 5;      // Distance bins: 0-14, 15-29, 30-49, 50-99, 100+

// Base encoding
alignas(64) static const int8_t BASE_MAP[256] = {
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

// Genetic code
static const char* GENETIC_CODE =
    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

inline int encode_codon(const char* s) {
    int8_t b0 = BASE_MAP[(uint8_t)s[0]];
    int8_t b1 = BASE_MAP[(uint8_t)s[1]];
    int8_t b2 = BASE_MAP[(uint8_t)s[2]];
    if ((b0 | b1 | b2) < 0) return -1;
    return (b0 << 4) | (b1 << 2) | b2;
}

inline char translate_codon(int code) {
    if (code < 0 || code >= 64) return 'X';
    return GENETIC_CODE[code];
}

// Decode codon index to string
inline std::string decode_codon(int c) {
    static const char* BASES = "ACGT";
    std::string s(3, ' ');
    s[0] = BASES[(c >> 4) & 3];
    s[1] = BASES[(c >> 2) & 3];
    s[2] = BASES[c & 3];
    return s;
}

// AA encoding: standard order ACDEFGHIKLMNPQRSTVWY (alphabetical except C before D)
// Index: A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9, M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19
inline int aa_to_idx(char aa) {
    switch (aa & 0xDF) {  // Case-insensitive
        case 'A': return 0;
        case 'C': return 1;
        case 'D': return 2;
        case 'E': return 3;
        case 'F': return 4;
        case 'G': return 5;
        case 'H': return 6;
        case 'I': return 7;
        case 'K': return 8;
        case 'L': return 9;
        case 'M': return 10;
        case 'N': return 11;
        case 'P': return 12;
        case 'Q': return 13;
        case 'R': return 14;
        case 'S': return 15;
        case 'T': return 16;
        case 'V': return 17;
        case 'W': return 18;
        case 'Y': return 19;
        case '*': return 20;
        default: return -1;
    }
}

inline char idx_to_aa(int idx) {
    static const char* AAS = "ACDEFGHIKLMNPQRSTVWY*";
    return (idx >= 0 && idx < 21) ? AAS[idx] : 'X';
}

// CDS validation
inline bool has_valid_start(const char* seq, size_t len) {
    if (len < 3) return false;
    uint8_t c0 = seq[0] & 0xDF, c1 = seq[1] & 0xDF, c2 = seq[2] & 0xDF;
    return (c1 == 'T' && c2 == 'G' && (c0 == 'A' || c0 == 'G' || c0 == 'T'));
}

inline bool has_valid_stop(const char* seq, size_t len) {
    if (len < 3) return false;
    const char* end = seq + len - 3;
    uint8_t c0 = end[0] & 0xDF, c1 = end[1] & 0xDF, c2 = end[2] & 0xDF;
    return (c0 == 'T' && ((c1 == 'A' && (c2 == 'A' || c2 == 'G')) || (c1 == 'G' && c2 == 'A')));
}

// Codon indices for damage-relevant codons
constexpr int CGG = 0b100110;  // 38 - Arg (damage-susceptible)
constexpr int TGG = 0b110110;  // 54 - Trp (damage product)
constexpr int CGC = 0b100101;  // 37 - Arg (damage-susceptible)
constexpr int CGT = 0b100111;  // 39 - Arg (damage-susceptible)
constexpr int TGC = 0b110101;  // 53 - Cys (damage product)
constexpr int TGT = 0b110111;  // 55 - Cys (damage product)
constexpr int CGA = 0b100100;  // 36 - Arg (→ stop TGA)
constexpr int AGA = 0b001100;  // 12 - Arg
constexpr int AGG = 0b001110;  // 14 - Arg

// Distance bin helper (for window-based generalization)
inline int get_distance_bin(int distance_from_5prime) {
    if (distance_from_5prime < 15) return 0;      // Damage zone
    if (distance_from_5prime < 30) return 1;      // Near damage zone
    if (distance_from_5prime < 50) return 2;      // Intermediate
    if (distance_from_5prime < 100) return 3;     // Far
    return 4;                                      // Very far (baseline)
}

// Thread state
struct ThreadState {
    // N-terminal positional AA frequencies (positions 1-10 after Met)
    std::array<std::array<uint64_t, 21>, N_TERMINAL_POSITIONS> nterminal_aa{};

    // Distance-binned AA frequencies (for sliding window generalization)
    std::array<std::array<uint64_t, 21>, NUM_DISTANCE_BINS> distance_binned_aa{};

    // Global AA frequencies for comparison
    std::array<uint64_t, 21> global_aa{};

    // Codon bigram: preceding codon counts for damage-relevant codons
    // Index by [target_codon][preceding_codon]
    std::array<std::array<uint64_t, 64>, 64> codon_bigram{};

    // Tripeptide counts (21^3 = 9261)
    std::array<uint64_t, 21*21*21> tripeptide{};

    // W-containing tripeptide patterns specifically
    // XWX patterns where W is in middle position
    std::array<uint64_t, 21*21> xwx_pattern{};  // [X1][X3] where middle is W
    // XRX patterns for comparison
    std::array<uint64_t, 21*21> xrx_pattern{};

    uint64_t total_aa = 0;
    uint64_t total_tripeptides = 0;
    uint64_t valid_cds = 0;

    std::string seq_buffer;

    ThreadState() {
        seq_buffer.reserve(100000);
        for (auto& row : nterminal_aa) row.fill(0);
        for (auto& row : distance_binned_aa) row.fill(0);
        global_aa.fill(0);
        for (auto& row : codon_bigram) row.fill(0);
        tripeptide.fill(0);
        xwx_pattern.fill(0);
        xrx_pattern.fill(0);
    }
};

// Global aggregated state
std::array<std::array<std::atomic<uint64_t>, 21>, N_TERMINAL_POSITIONS> g_nterminal_aa;
std::array<std::array<std::atomic<uint64_t>, 21>, NUM_DISTANCE_BINS> g_distance_binned_aa;
std::array<std::atomic<uint64_t>, 21> g_global_aa;
std::array<std::array<std::atomic<uint64_t>, 64>, 64> g_codon_bigram;
std::array<std::atomic<uint64_t>, 21*21*21> g_tripeptide;
std::array<std::atomic<uint64_t>, 21*21> g_xwx_pattern;
std::array<std::atomic<uint64_t>, 21*21> g_xrx_pattern;
std::atomic<uint64_t> g_total_aa{0};
std::atomic<uint64_t> g_total_tripeptides{0};
std::atomic<uint64_t> g_valid_cds{0};

void aggregate_state(const ThreadState& s) {
    for (int p = 0; p < N_TERMINAL_POSITIONS; p++) {
        for (int a = 0; a < 21; a++) {
            g_nterminal_aa[p][a] += s.nterminal_aa[p][a];
        }
    }
    for (int b = 0; b < NUM_DISTANCE_BINS; b++) {
        for (int a = 0; a < 21; a++) {
            g_distance_binned_aa[b][a] += s.distance_binned_aa[b][a];
        }
    }
    for (int a = 0; a < 21; a++) {
        g_global_aa[a] += s.global_aa[a];
    }
    for (int i = 0; i < 64; i++) {
        for (int j = 0; j < 64; j++) {
            g_codon_bigram[i][j] += s.codon_bigram[i][j];
        }
    }
    for (size_t i = 0; i < s.tripeptide.size(); i++) {
        g_tripeptide[i] += s.tripeptide[i];
    }
    for (size_t i = 0; i < s.xwx_pattern.size(); i++) {
        g_xwx_pattern[i] += s.xwx_pattern[i];
        g_xrx_pattern[i] += s.xrx_pattern[i];
    }
    g_total_aa += s.total_aa;
    g_total_tripeptides += s.total_tripeptides;
    g_valid_cds += s.valid_cds;
}

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

        const char* s = seq.c_str();
        if (!has_valid_start(s, len)) return;
        if (!has_valid_stop(s, len)) return;
        if (len % 3 != 0) return;

        state.valid_cds++;

        size_t num_codons = len / 3;
        std::vector<int> codons(num_codons);
        std::vector<int> aas(num_codons);

        // Translate all codons
        for (size_t i = 0; i < num_codons; i++) {
            int codon = encode_codon(s + i * 3);
            codons[i] = codon;
            if (codon >= 0) {
                char aa = translate_codon(codon);
                aas[i] = aa_to_idx(aa);
            } else {
                aas[i] = -1;
            }
        }

        // Extract patterns
        for (size_t i = 0; i < num_codons; i++) {
            if (aas[i] < 0 || aas[i] >= 20) continue;  // Skip stops/unknown

            // N-terminal positional frequencies (positions 1-10 after Met)
            // Position 0 is always Met, so we track positions 1+
            if (i > 0 && i <= N_TERMINAL_POSITIONS) {
                state.nterminal_aa[i-1][aas[i]]++;
            }

            // Distance-binned AA frequencies (for window-based generalization)
            // Distance from 5' is measured in codons (= AA positions)
            int bin = get_distance_bin(i);
            state.distance_binned_aa[bin][aas[i]]++;

            // Global AA frequencies
            state.global_aa[aas[i]]++;
            state.total_aa++;

            // Codon bigram (preceding codon context)
            if (i > 0 && codons[i] >= 0 && codons[i-1] >= 0) {
                state.codon_bigram[codons[i]][codons[i-1]]++;
            }

            // Tripeptides
            if (i >= 2 && aas[i-2] >= 0 && aas[i-1] >= 0 && aas[i-2] < 21 && aas[i-1] < 21 && aas[i] < 21) {
                int idx = aas[i-2] * 21 * 21 + aas[i-1] * 21 + aas[i];
                state.tripeptide[idx]++;
                state.total_tripeptides++;

                // Track XWX patterns (W=19 in our encoding, but need to check actual char)
                char middle_aa = idx_to_aa(aas[i-1]);
                if (middle_aa == 'W') {
                    int xwx_idx = aas[i-2] * 21 + aas[i];
                    state.xwx_pattern[xwx_idx]++;
                } else if (middle_aa == 'R') {
                    int xrx_idx = aas[i-2] * 21 + aas[i];
                    state.xrx_pattern[xrx_idx]++;
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
            seq.clear();
        } else {
            seq.append(buffer, slen);
        }
    }

    process_sequence();
    gzclose(file);
}

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
            if ((nlen > 7 && name.compare(nlen-7, 7, ".fna.gz") == 0) ||
                (nlen > 9 && name.compare(nlen-9, 9, ".fasta.gz") == 0)) {
                files.push_back(full_path);
            }
        }
    }
    closedir(dir);
}

void print_nterminal_bias() {
    std::cerr << "\n=== N-TERMINAL POSITIONAL AA FREQUENCIES ===\n";
    std::cerr << "Position-specific frequencies vs global (for damage-relevant AAs)\n\n";

    double total_global = g_total_aa.load();
    if (total_global == 0) return;

    // Focus on W, R, C (damage-relevant)
    std::vector<std::pair<char, int>> focus_aas = {{'W', -1}, {'R', -1}, {'C', -1}};
    for (auto& [aa, idx] : focus_aas) {
        idx = aa_to_idx(aa);
    }

    std::cerr << "Pos\t";
    for (auto& [aa, idx] : focus_aas) {
        std::cerr << aa << "_freq\t" << aa << "_ratio\t";
    }
    std::cerr << "\n";

    for (int p = 0; p < N_TERMINAL_POSITIONS; p++) {
        uint64_t pos_total = 0;
        for (int a = 0; a < 20; a++) {
            pos_total += g_nterminal_aa[p][a].load();
        }
        if (pos_total == 0) continue;

        std::cerr << (p+1) << "\t";
        for (auto& [aa, idx] : focus_aas) {
            double pos_freq = (double)g_nterminal_aa[p][idx].load() / pos_total;
            double global_freq = (double)g_global_aa[idx].load() / total_global;
            double ratio = (global_freq > 0) ? pos_freq / global_freq : 0;
            std::cerr << std::fixed << std::setprecision(4) << pos_freq << "\t"
                      << std::setprecision(2) << ratio << "\t";
        }
        std::cerr << "\n";
    }

    // Output as C++ code
    std::cout << "\n// N-terminal AA frequency ratios (position-specific / global)\n";
    std::cout << "// Position 1 = first AA after Met\n";
    std::cout << "static constexpr float NTERMINAL_FREQ_RATIO[" << N_TERMINAL_POSITIONS << "][20] = {\n";
    for (int p = 0; p < N_TERMINAL_POSITIONS; p++) {
        uint64_t pos_total = 0;
        for (int a = 0; a < 20; a++) {
            pos_total += g_nterminal_aa[p][a].load();
        }

        std::cout << "    {";
        for (int a = 0; a < 20; a++) {
            double pos_freq = (pos_total > 0) ? (double)g_nterminal_aa[p][a].load() / pos_total : 0;
            double global_freq = (double)g_global_aa[a].load() / total_global;
            double ratio = (global_freq > 0) ? pos_freq / global_freq : 1.0;
            std::cout << std::fixed << std::setprecision(3) << ratio;
            if (a < 19) std::cout << ", ";
        }
        std::cout << "}";
        if (p < N_TERMINAL_POSITIONS - 1) std::cout << ",";
        std::cout << "  // pos " << (p+1) << "\n";
    }
    std::cout << "};\n";
}

void print_distance_binned_ratios() {
    std::cerr << "\n=== DISTANCE-BINNED AA FREQUENCIES ===\n";
    std::cerr << "For window-based generalization (independent of absolute position)\n\n";

    static const char* BIN_NAMES[] = {"0-14", "15-29", "30-49", "50-99", "100+"};

    double total_global = g_total_aa.load();
    if (total_global == 0) return;

    // Focus on W, R, C
    int w_idx = aa_to_idx('W');
    int r_idx = aa_to_idx('R');
    int c_idx = aa_to_idx('C');

    std::cerr << "Bin\tW_ratio\tR_ratio\tC_ratio\n";
    for (int b = 0; b < NUM_DISTANCE_BINS; b++) {
        uint64_t bin_total = 0;
        for (int a = 0; a < 20; a++) {
            bin_total += g_distance_binned_aa[b][a].load();
        }
        if (bin_total == 0) continue;

        double w_freq = (double)g_distance_binned_aa[b][w_idx].load() / bin_total;
        double r_freq = (double)g_distance_binned_aa[b][r_idx].load() / bin_total;
        double c_freq = (double)g_distance_binned_aa[b][c_idx].load() / bin_total;

        double w_global = (double)g_global_aa[w_idx].load() / total_global;
        double r_global = (double)g_global_aa[r_idx].load() / total_global;
        double c_global = (double)g_global_aa[c_idx].load() / total_global;

        std::cerr << BIN_NAMES[b] << "\t"
                  << std::fixed << std::setprecision(2) << w_freq / w_global << "\t"
                  << r_freq / r_global << "\t"
                  << c_freq / c_global << "\n";
    }

    // Output as C++ code
    std::cout << "\n// Distance-binned AA frequency ratios (bin_freq / global_freq)\n";
    std::cout << "// Bins: 0-14 (damage zone), 15-29, 30-49, 50-99, 100+ (baseline)\n";
    std::cout << "static constexpr float DISTANCE_BINNED_FREQ_RATIO[" << NUM_DISTANCE_BINS << "][20] = {\n";
    for (int b = 0; b < NUM_DISTANCE_BINS; b++) {
        uint64_t bin_total = 0;
        for (int a = 0; a < 20; a++) {
            bin_total += g_distance_binned_aa[b][a].load();
        }

        std::cout << "    {";
        for (int a = 0; a < 20; a++) {
            double bin_freq = (bin_total > 0) ? (double)g_distance_binned_aa[b][a].load() / bin_total : 0;
            double global_freq = (double)g_global_aa[a].load() / total_global;
            double ratio = (global_freq > 0) ? bin_freq / global_freq : 1.0;
            std::cout << std::fixed << std::setprecision(3) << ratio;
            if (a < 19) std::cout << ", ";
        }
        std::cout << "}";
        if (b < NUM_DISTANCE_BINS - 1) std::cout << ",";
        std::cout << "  // " << BIN_NAMES[b] << "\n";
    }
    std::cout << "};\n";
}

void print_codon_bigram_context() {
    std::cerr << "\n=== CODON BIGRAM CONTEXT FOR CGG vs TGG ===\n";
    std::cerr << "What codons precede CGG (Arg) vs TGG (Trp)?\n\n";

    // Get totals for normalization
    uint64_t cgg_total = 0, tgg_total = 0;
    for (int j = 0; j < 64; j++) {
        cgg_total += g_codon_bigram[CGG][j].load();
        tgg_total += g_codon_bigram[TGG][j].load();
    }

    if (cgg_total == 0 || tgg_total == 0) return;

    // Compute LLR for each preceding codon: log2(P(codon|CGG) / P(codon|TGG))
    std::vector<std::tuple<double, int, double, double>> bigram_llrs;  // LLR, codon, P_CGG, P_TGG

    for (int j = 0; j < 64; j++) {
        double p_cgg = (double)g_codon_bigram[CGG][j].load() / cgg_total;
        double p_tgg = (double)g_codon_bigram[TGG][j].load() / tgg_total;

        double llr = 0;
        if (p_cgg > 0 && p_tgg > 0) {
            llr = std::log2(p_cgg / p_tgg);
        } else if (p_cgg > 0) {
            llr = 5.0;  // Cap
        } else if (p_tgg > 0) {
            llr = -5.0;
        }

        bigram_llrs.push_back({llr, j, p_cgg, p_tgg});
    }

    std::sort(bigram_llrs.rbegin(), bigram_llrs.rend());

    std::cerr << "Top 15 codons that favor CGG over TGG (strong R→W damage signal):\n";
    std::cerr << "Codon\tAA\tLLR\tP(|CGG)\tP(|TGG)\tOdds\n";
    for (int i = 0; i < 15 && i < (int)bigram_llrs.size(); i++) {
        auto [llr, codon, p_cgg, p_tgg] = bigram_llrs[i];
        std::cerr << decode_codon(codon) << "\t" << translate_codon(codon) << "\t"
                  << std::fixed << std::setprecision(2) << llr << "\t"
                  << std::setprecision(4) << p_cgg << "\t" << p_tgg << "\t"
                  << std::setprecision(1) << std::pow(2.0, llr) << "x\n";
    }

    // Output as C++ lookup table
    std::cout << "\n// Codon bigram LLR: log2(P(codon precedes CGG) / P(codon precedes TGG))\n";
    std::cout << "// Positive = preceding codon favors original being CGG (Arg)\n";
    std::cout << "static constexpr float CODON_BIGRAM_LLR_CGG_TGG[64] = {\n";
    for (int j = 0; j < 64; j++) {
        double p_cgg = (cgg_total > 0) ? (double)g_codon_bigram[CGG][j].load() / cgg_total : 0;
        double p_tgg = (tgg_total > 0) ? (double)g_codon_bigram[TGG][j].load() / tgg_total : 0;

        double llr = 0;
        if (p_cgg > 1e-9 && p_tgg > 1e-9) {
            llr = std::log2(p_cgg / p_tgg);
            llr = std::max(-3.0, std::min(3.0, llr));  // Clamp
        }

        std::cout << std::fixed << std::setprecision(3) << llr;
        if (j < 63) std::cout << ", ";
        if ((j + 1) % 8 == 0) std::cout << "  // " << decode_codon(j-7) << "-" << decode_codon(j) << "\n";
    }
    std::cout << "};\n";

    // Same for CGC/CGT vs TGC/TGT (Arg→Cys)
    std::cerr << "\n=== CODON BIGRAM CONTEXT FOR CG[CT] vs TG[CT] ===\n";

    uint64_t cgc_total = 0, cgc_cgt_total = 0, tgc_total = 0, tgc_tgt_total = 0;
    for (int j = 0; j < 64; j++) {
        cgc_total += g_codon_bigram[CGC][j].load();
        cgc_cgt_total += g_codon_bigram[CGC][j].load() + g_codon_bigram[CGT][j].load();
        tgc_total += g_codon_bigram[TGC][j].load();
        tgc_tgt_total += g_codon_bigram[TGC][j].load() + g_codon_bigram[TGT][j].load();
    }

    if (cgc_cgt_total > 0 && tgc_tgt_total > 0) {
        std::cout << "\n// Codon bigram LLR for CGC/CGT vs TGC/TGT (Arg→Cys)\n";
        std::cout << "static constexpr float CODON_BIGRAM_LLR_CGX_TGX[64] = {\n";
        for (int j = 0; j < 64; j++) {
            double p_cgx = (double)(g_codon_bigram[CGC][j].load() + g_codon_bigram[CGT][j].load()) / cgc_cgt_total;
            double p_tgx = (double)(g_codon_bigram[TGC][j].load() + g_codon_bigram[TGT][j].load()) / tgc_tgt_total;

            double llr = 0;
            if (p_cgx > 1e-9 && p_tgx > 1e-9) {
                llr = std::log2(p_cgx / p_tgx);
                llr = std::max(-3.0, std::min(3.0, llr));
            }

            std::cout << std::fixed << std::setprecision(3) << llr;
            if (j < 63) std::cout << ", ";
            if ((j + 1) % 8 == 0) std::cout << "\n";
        }
        std::cout << "};\n";
    }
}

void print_forbidden_w_tripeptides() {
    std::cerr << "\n=== FORBIDDEN W-CONTAINING TRIPEPTIDES ===\n";
    std::cerr << "XWX patterns that are extremely rare (likely damage artifacts)\n\n";

    double total_tripeptides = g_total_tripeptides.load();
    if (total_tripeptides == 0) return;

    // Calculate expected frequencies under independence
    double total_aa = g_total_aa.load();

    // Collect XWX vs XRX ratios
    std::vector<std::tuple<double, int, int, uint64_t, uint64_t>> xwx_llrs;
    // LLR, first_aa, third_aa, xwx_count, xrx_count

    int w_idx = aa_to_idx('W');
    int r_idx = aa_to_idx('R');

    for (int a1 = 0; a1 < 20; a1++) {
        for (int a3 = 0; a3 < 20; a3++) {
            uint64_t xwx_count = g_xwx_pattern[a1 * 21 + a3].load();
            uint64_t xrx_count = g_xrx_pattern[a1 * 21 + a3].load();

            // LLR = log2(P(XWX observed) / P(XWX expected by independence))
            // Expected = P(a1) * P(W) * P(a3) * total
            double p1 = (double)g_global_aa[a1].load() / total_aa;
            double pw = (double)g_global_aa[w_idx].load() / total_aa;
            double p3 = (double)g_global_aa[a3].load() / total_aa;
            double expected_xwx = p1 * pw * p3 * total_tripeptides;

            double llr = -10.0;
            if (xwx_count > 0 && expected_xwx > 0) {
                llr = std::log2((double)xwx_count / expected_xwx);
            }

            xwx_llrs.push_back({llr, a1, a3, xwx_count, xrx_count});
        }
    }

    std::sort(xwx_llrs.begin(), xwx_llrs.end());

    std::cerr << "Top 30 most underrepresented XWX patterns:\n";
    std::cerr << "Pattern\tLLR\tCount\tXRX_count\tRatio\n";
    for (int i = 0; i < 30 && i < (int)xwx_llrs.size(); i++) {
        auto [llr, a1, a3, xwx, xrx] = xwx_llrs[i];
        char c1 = idx_to_aa(a1);
        char c3 = idx_to_aa(a3);
        double ratio = (xwx > 0) ? (double)xrx / xwx : 999;
        std::cerr << c1 << "W" << c3 << "\t"
                  << std::fixed << std::setprecision(2) << llr << "\t"
                  << xwx << "\t" << xrx << "\t"
                  << std::setprecision(1) << ratio << "x\n";
    }

    // Output forbidden patterns as C++ code
    std::cout << "\n// Forbidden W tripeptide patterns (XWX with LLR < -2)\n";
    std::cout << "// If we see these patterns with W at damage position, very likely was R\n";
    std::cout << "struct ForbiddenWPattern {\n";
    std::cout << "    char aa_before;\n";
    std::cout << "    char aa_after;\n";
    std::cout << "    float llr;        // log2(observed/expected)\n";
    std::cout << "    float r_ratio;    // How much more common XRX is than XWX\n";
    std::cout << "};\n\n";
    std::cout << "static constexpr ForbiddenWPattern FORBIDDEN_W_PATTERNS[] = {\n";

    // Sort by R/W ratio instead of LLR (more directly useful for correction)
    std::sort(xwx_llrs.begin(), xwx_llrs.end(), [](const auto& a, const auto& b) {
        auto [llr_a, a1_a, a3_a, xwx_a, xrx_a] = a;
        auto [llr_b, a1_b, a3_b, xwx_b, xrx_b] = b;
        double ratio_a = (xwx_a > 0) ? (double)xrx_a / xwx_a : 0;
        double ratio_b = (xwx_b > 0) ? (double)xrx_b / xwx_b : 0;
        return ratio_a > ratio_b;  // Descending by ratio
    });

    int count = 0;
    for (auto& [llr, a1, a3, xwx, xrx] : xwx_llrs) {
        double ratio = (xwx > 0) ? (double)xrx / xwx : 0;
        if (ratio < 5.0) break;  // Only include patterns with R:W ratio >= 5
        if (count >= 100) break;  // Limit output

        char c1 = idx_to_aa(a1);
        char c3 = idx_to_aa(a3);

        std::cout << "    {'" << c1 << "', '" << c3 << "', "
                  << std::fixed << std::setprecision(2) << llr << "f, "
                  << std::setprecision(1) << ratio << "f},\n";
        count++;
    }
    std::cout << "};\n";
    std::cout << "static constexpr size_t NUM_FORBIDDEN_W_PATTERNS = " << count << ";\n";
}

int main(int argc, char* argv[]) {
    std::string domain = "gtdb";
    std::vector<std::string> input_dirs;
    int threads = 16;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--domain" && i + 1 < argc) {
            domain = argv[++i];
        } else if (arg == "--threads" && i + 1 < argc) {
            threads = std::stoi(argv[++i]);
        } else if (arg[0] != '-') {
            input_dirs.push_back(arg);
        }
    }

    if (input_dirs.empty()) {
        std::cerr << "Usage: " << argv[0] << " [--domain name] [--threads N] <cds_dirs...>\n";
        std::cerr << "\nExtracts damage-correction context patterns from GTDB CDS sequences.\n";
        return 1;
    }

    // Initialize atomics
    for (auto& row : g_nterminal_aa) for (auto& v : row) v = 0;
    for (auto& row : g_distance_binned_aa) for (auto& v : row) v = 0;
    for (auto& v : g_global_aa) v = 0;
    for (auto& row : g_codon_bigram) for (auto& v : row) v = 0;
    for (auto& v : g_tripeptide) v = 0;
    for (auto& v : g_xwx_pattern) v = 0;
    for (auto& v : g_xrx_pattern) v = 0;

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
            if (done % 5000 == 0) {
                std::cerr << "Processed " << done << "/" << files.size() << " files\n";
            }
        }

        #pragma omp critical
        { aggregate_state(state); }
    }

    std::cerr << "\n=== STATISTICS ===\n";
    std::cerr << "Valid CDS: " << g_valid_cds << "\n";
    std::cerr << "Total AA: " << g_total_aa << "\n";
    std::cerr << "Total tripeptides: " << g_total_tripeptides << "\n";

    // Output header
    std::cout << "/**\n";
    std::cout << " * Damage-correction context tables extracted from GTDB.\n";
    std::cout << " * Domain: " << domain << "\n";
    std::cout << " * Valid CDS: " << g_valid_cds << "\n";
    std::cout << " * Total AA: " << g_total_aa << "\n";
    std::cout << " * Auto-generated by extract_damage_context\n";
    std::cout << " */\n\n";
    std::cout << "#pragma once\n\n";
    std::cout << "#include <cstddef>\n\n";
    std::cout << "namespace dart {\n\n";

    print_nterminal_bias();
    print_distance_binned_ratios();
    print_codon_bigram_context();
    print_forbidden_w_tripeptides();

    std::cout << "\n} // namespace dart\n";

    return 0;
}
