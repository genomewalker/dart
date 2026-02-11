/**
 * Extract tripeptide frequencies and rare patterns from CDS sequences.
 *
 * Tripeptides provide stronger signal than dipeptides because:
 * - 8000 combinations vs 400 (more sparse, more extreme values)
 * - Context-dependent: captures 3-AA motifs that may be forbidden
 *
 * Also extracts:
 * - Codon context around damage-susceptible codons
 * - Position-specific codon usage (start/internal/end)
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

// AA encoding using standard order ARNDCQEGHILKMFPSTWYV
inline int aa_to_idx(char aa) {
    static const int8_t map[256] = {
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 4, 3, 6,13, 7, 8,-1,11,12,10, 9,14,-1,  // A=0,C=4,D=3,E=6,F=13,G=7,H=8,I=9,K=11,L=12,M=10,N=14
        15,16, 1,17,18,-1,19,19,-1,20,-1,-1,-1,-1,-1,-1,  // P=15,Q=16,R=1,S=17,T=18,V=19,W=19,Y=20
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
    return map[(uint8_t)aa];
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

// Thread state
struct ThreadState {
    // Tripeptide counts: 21^3 = 9261 (including stop)
    std::array<uint64_t, 21*21*21> tripeptide{};
    std::array<uint64_t, 21> aa_count{};

    // Codon context: for each codon, what codons precede/follow it
    std::array<std::array<uint64_t, 64>, 64> codon_context_before{};
    std::array<std::array<uint64_t, 64>, 64> codon_context_after{};

    // Position-specific codon usage
    std::array<uint64_t, 64> codon_at_start{};    // Codon 0 (after start)
    std::array<uint64_t, 64> codon_at_end{};      // Codon -2 (before stop)
    std::array<uint64_t, 64> codon_internal{};    // All other codons

    uint64_t total_tripeptides = 0;
    uint64_t total_aa = 0;
    uint64_t valid_cds = 0;

    std::string seq_buffer;

    ThreadState() {
        seq_buffer.reserve(100000);
        tripeptide.fill(0);
        aa_count.fill(0);
        codon_at_start.fill(0);
        codon_at_end.fill(0);
        codon_internal.fill(0);
        for (auto& row : codon_context_before) row.fill(0);
        for (auto& row : codon_context_after) row.fill(0);
    }
};

// Global aggregated state
std::array<std::atomic<uint64_t>, 21*21*21> g_tripeptide;
std::array<std::atomic<uint64_t>, 21> g_aa_count;
std::array<std::array<std::atomic<uint64_t>, 64>, 64> g_codon_context_before;
std::array<std::array<std::atomic<uint64_t>, 64>, 64> g_codon_context_after;
std::array<std::atomic<uint64_t>, 64> g_codon_at_start;
std::array<std::atomic<uint64_t>, 64> g_codon_at_end;
std::array<std::atomic<uint64_t>, 64> g_codon_internal;
std::atomic<uint64_t> g_total_tripeptides{0};
std::atomic<uint64_t> g_total_aa{0};
std::atomic<uint64_t> g_valid_cds{0};

void aggregate_state(const ThreadState& s) {
    for (size_t i = 0; i < s.tripeptide.size(); i++) {
        g_tripeptide[i] += s.tripeptide[i];
    }
    for (int i = 0; i < 21; i++) {
        g_aa_count[i] += s.aa_count[i];
    }
    for (int i = 0; i < 64; i++) {
        g_codon_at_start[i] += s.codon_at_start[i];
        g_codon_at_end[i] += s.codon_at_end[i];
        g_codon_internal[i] += s.codon_internal[i];
        for (int j = 0; j < 64; j++) {
            g_codon_context_before[i][j] += s.codon_context_before[i][j];
            g_codon_context_after[i][j] += s.codon_context_after[i][j];
        }
    }
    g_total_tripeptides += s.total_tripeptides;
    g_total_aa += s.total_aa;
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
                if (aas[i] < 0 && aa != '*') aas[i] = 20;  // Unknown -> index 20
            } else {
                aas[i] = -1;
            }
        }

        // Extract patterns
        for (size_t i = 0; i < num_codons; i++) {
            if (codons[i] < 0) continue;

            // Position-specific codon usage
            if (i == 1) {  // First codon after start
                state.codon_at_start[codons[i]]++;
            } else if (i == num_codons - 2) {  // Last codon before stop
                state.codon_at_end[codons[i]]++;
            } else if (i > 1 && i < num_codons - 2) {
                state.codon_internal[codons[i]]++;
            }

            // Codon context
            if (i > 0 && codons[i-1] >= 0) {
                state.codon_context_before[codons[i]][codons[i-1]]++;
            }
            if (i < num_codons - 1 && codons[i+1] >= 0) {
                state.codon_context_after[codons[i]][codons[i+1]]++;
            }

            // AA counts (skip stops)
            if (aas[i] >= 0 && aas[i] < 20) {
                state.aa_count[aas[i]]++;
                state.total_aa++;
            }

            // Tripeptides (skip if any unknown)
            if (i >= 2 && aas[i] >= 0 && aas[i-1] >= 0 && aas[i-2] >= 0 &&
                aas[i] < 21 && aas[i-1] < 21 && aas[i-2] < 21) {
                int idx = aas[i-2] * 21 * 21 + aas[i-1] * 21 + aas[i];
                state.tripeptide[idx]++;
                state.total_tripeptides++;
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

void print_rare_tripeptides() {
    // AA order for output
    static const char* AA_ORDER = "ARNDCQEGHILKMFPSTWYV*";

    std::cerr << "\n=== RAREST TRIPEPTIDES ===\n";

    // Collect all tripeptides with their LLR
    std::vector<std::tuple<double, int, int, int, uint64_t>> rare;

    double total = g_total_tripeptides.load();
    if (total == 0) return;

    for (int a1 = 0; a1 < 20; a1++) {
        for (int a2 = 0; a2 < 20; a2++) {
            for (int a3 = 0; a3 < 20; a3++) {
                int idx = a1 * 21 * 21 + a2 * 21 + a3;
                uint64_t count = g_tripeptide[idx].load();

                // Expected by independence
                double p1 = (double)g_aa_count[a1].load() / g_total_aa.load();
                double p2 = (double)g_aa_count[a2].load() / g_total_aa.load();
                double p3 = (double)g_aa_count[a3].load() / g_total_aa.load();
                double expected = p1 * p2 * p3 * total;

                double llr = -10.0;
                if (count > 0 && expected > 0) {
                    llr = std::log2((double)count / expected);
                }

                rare.push_back({llr, a1, a2, a3, count});
            }
        }
    }

    std::sort(rare.begin(), rare.end());

    std::cerr << "Top 30 underrepresented tripeptides:\n";
    for (int i = 0; i < 30 && i < (int)rare.size(); i++) {
        auto [llr, a1, a2, a3, count] = rare[i];
        std::cerr << "  " << AA_ORDER[a1] << AA_ORDER[a2] << AA_ORDER[a3]
                  << " LLR=" << std::fixed << std::setprecision(2) << llr
                  << " (count=" << count << ")\n";
    }
}

void print_damage_codon_context() {
    static const char* BASES = "ACGT";

    auto decode_codon = [](int c) {
        std::string s(3, ' ');
        s[0] = BASES[(c >> 4) & 3];
        s[1] = BASES[(c >> 2) & 3];
        s[2] = BASES[c & 3];
        return s;
    };

    // Damage-susceptible codons (contain C that can become T)
    std::vector<std::pair<int, std::string>> damage_codons = {
        {0b010010, "CGG→TGG (R→W)"},  // CGG = 22
        {0b010001, "CGC→TGC (R→C)"},  // CGC = 21
        {0b010011, "CGT→TGT (R→C)"},  // CGT = 23
        {0b010000, "CGA→TGA (R→*)"},  // CGA = 20
        {0b000100, "CAA→TAA (Q→*)"},  // CAA = 4
        {0b000110, "CAG→TAG (Q→*)"},  // CAG = 6
    };

    std::cerr << "\n=== CODON CONTEXT FOR DAMAGE-SUSCEPTIBLE CODONS ===\n";

    for (auto& [codon_idx, desc] : damage_codons) {
        std::cerr << "\n" << desc << " - preceding codons:\n";

        uint64_t total_before = 0;
        for (int j = 0; j < 64; j++) {
            total_before += g_codon_context_before[codon_idx][j].load();
        }

        if (total_before == 0) continue;

        // Find most common preceding codons
        std::vector<std::pair<double, int>> precs;
        for (int j = 0; j < 64; j++) {
            double freq = (double)g_codon_context_before[codon_idx][j].load() / total_before;
            precs.push_back({freq, j});
        }
        std::sort(precs.rbegin(), precs.rend());

        for (int i = 0; i < 5; i++) {
            std::cerr << "  " << decode_codon(precs[i].second)
                      << " (" << translate_codon(precs[i].second) << "): "
                      << std::fixed << std::setprecision(1) << precs[i].first * 100 << "%\n";
        }
    }
}

int main(int argc, char* argv[]) {
    std::string domain;
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

    if (domain.empty() || input_dirs.empty()) {
        std::cerr << "Usage: " << argv[0] << " --domain <domain> <cds_dirs...>\n";
        return 1;
    }

    // Initialize
    for (auto& v : g_tripeptide) v = 0;
    for (auto& v : g_aa_count) v = 0;
    for (auto& v : g_codon_at_start) v = 0;
    for (auto& v : g_codon_at_end) v = 0;
    for (auto& v : g_codon_internal) v = 0;
    for (auto& row : g_codon_context_before) for (auto& v : row) v = 0;
    for (auto& row : g_codon_context_after) for (auto& v : row) v = 0;

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

            if (++files_done % 1000 == 0) {
                std::cerr << "Processed " << files_done << "/" << files.size() << "\n";
            }
        }

        #pragma omp critical
        { aggregate_state(state); }
    }

    std::cerr << "\n=== STATISTICS ===\n";
    std::cerr << "Valid CDS: " << g_valid_cds << "\n";
    std::cerr << "Total AA: " << g_total_aa << "\n";
    std::cerr << "Total tripeptides: " << g_total_tripeptides << "\n";

    print_rare_tripeptides();
    print_damage_codon_context();

    return 0;
}
