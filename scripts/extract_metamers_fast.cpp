/*
 * AA 6-mer extractor using direct array indexing
 *
 * Key optimization: 6-mer AA = 21^6 = 85,766,121 possible values
 * Use direct array indexing (O(1)) instead of hash map (O(1) avg with collisions)
 *
 * Compile:
 *   g++ -O3 -march=native -fopenmp -o extract_metamers_fast extract_metamers_fast.cpp -lz
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <atomic>
#include <dirent.h>
#include <zlib.h>
#include <omp.h>
#include <sys/stat.h>
#include <chrono>

// ============== CONSTANTS ==============
constexpr size_t KMER_AA = 6;
constexpr size_t KMER_NT = 18;  // 6 codons = 18 nt
constexpr size_t BOUNDARY_CODONS = 10;  // codons near start/stop
constexpr size_t MIN_GENE_CODONS = 30;
constexpr size_t AA_SPACE = 85766121;  // 21^6
constexpr size_t READ_BUFFER = 1 << 18;

// AA encoding: 0-19 = standard AAs, 20 = stop
alignas(64) static const int8_t CODON_TO_AA[64] = {
    11, 2, 11, 2, 16, 16, 16, 16,  // AAA-ACT
    1, 15, 1, 15, 9, 9, 12, 9,     // AGA-ATT
    5, 8, 5, 8, 14, 14, 14, 14,    // CAA-CCT
    1, 1, 1, 1, 10, 10, 10, 10,    // CGA-CTT
    6, 3, 6, 3, 0, 0, 0, 0,        // GAA-GCT
    7, 7, 7, 7, 19, 19, 19, 19,    // GGA-GTT
    20, 18, 20, 18, 15, 15, 15, 15, // TAA-TCT
    20, 4, 17, 4, 10, 13, 10, 13   // TGA-TTT
};

alignas(64) static const int8_t NT_MAP[256] = {
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

// Powers of 21 for encoding
constexpr uint32_t POW21[7] = {1, 21, 441, 9261, 194481, 4084101, 85766121};

// Encode 6 AAs to index (0 to 85766120)
inline int64_t encode_aa6(const int8_t* aas) {
    int64_t code = 0;
    for (int i = 0; i < 6; i++) {
        if (aas[i] < 0 || aas[i] > 20) return -1;
        code += aas[i] * POW21[i];
    }
    return code;
}

// Translate codon to AA
inline int8_t translate_codon(const char* seq) {
    int8_t b0 = NT_MAP[(uint8_t)seq[0]];
    int8_t b1 = NT_MAP[(uint8_t)seq[1]];
    int8_t b2 = NT_MAP[(uint8_t)seq[2]];
    if ((b0 | b1 | b2) < 0) return -1;
    return CODON_TO_AA[(b0 << 4) | (b1 << 2) | b2];
}

// Global counters using atomic arrays
struct Counters {
    std::vector<std::atomic<uint32_t>> coding;
    std::vector<std::atomic<uint32_t>> boundary_5p;
    std::vector<std::atomic<uint32_t>> boundary_3p;
    std::vector<std::atomic<uint32_t>> rc_coding;  // Reverse complement

    std::atomic<uint64_t> total_coding{0};
    std::atomic<uint64_t> total_5p{0};
    std::atomic<uint64_t> total_3p{0};
    std::atomic<uint64_t> total_rc{0};
    std::atomic<uint64_t> genes{0};

    Counters() : coding(AA_SPACE), boundary_5p(AA_SPACE),
                 boundary_3p(AA_SPACE), rc_coding(AA_SPACE) {
        for (size_t i = 0; i < AA_SPACE; i++) {
            coding[i] = 0;
            boundary_5p[i] = 0;
            boundary_3p[i] = 0;
            rc_coding[i] = 0;
        }
    }
};

// Reverse complement a sequence
inline char complement_base(char c) {
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 't'; case 't': return 'a';
        case 'c': return 'g'; case 'g': return 'c';
        default: return 'N';
    }
}

void reverse_complement(const std::string& seq, std::string& rc) {
    rc.resize(seq.size());
    for (size_t i = 0; i < seq.size(); i++) {
        rc[seq.size() - 1 - i] = complement_base(seq[i]);
    }
}

// Process a single gene sequence
void process_gene(const std::string& seq, Counters& cnt) {
    size_t len = seq.size();
    size_t n_codons = len / 3;
    if (n_codons < MIN_GENE_CODONS) return;

    cnt.genes++;

    const char* s = seq.c_str();
    int8_t aas[6];

    // 1. Extract CODING 6-mers (interior, skip boundaries)
    size_t start_codon = BOUNDARY_CODONS;
    size_t end_codon = n_codons - BOUNDARY_CODONS - 6;

    for (size_t c = start_codon; c < end_codon; c++) {
        bool valid = true;
        for (int i = 0; i < 6 && valid; i++) {
            aas[i] = translate_codon(s + (c + i) * 3);
            if (aas[i] < 0 || aas[i] == 20) valid = false;  // Skip stops
        }
        if (valid) {
            int64_t code = encode_aa6(aas);
            if (code >= 0) {
                cnt.coding[code]++;
                cnt.total_coding++;
            }
        }
    }

    // 2. Extract 5' BOUNDARY 6-mers (first BOUNDARY_CODONS codons)
    for (size_t c = 0; c < BOUNDARY_CODONS && c + 6 <= n_codons; c++) {
        bool valid = true;
        for (int i = 0; i < 6 && valid; i++) {
            aas[i] = translate_codon(s + (c + i) * 3);
            if (aas[i] < 0 || aas[i] == 20) valid = false;
        }
        if (valid) {
            int64_t code = encode_aa6(aas);
            if (code >= 0) {
                cnt.boundary_5p[code]++;
                cnt.total_5p++;
            }
        }
    }

    // 3. Extract 3' BOUNDARY 6-mers (last BOUNDARY_CODONS codons)
    size_t start_3p = (n_codons > BOUNDARY_CODONS + 6) ? n_codons - BOUNDARY_CODONS - 6 : 0;
    for (size_t c = start_3p; c + 6 <= n_codons; c++) {
        bool valid = true;
        for (int i = 0; i < 6 && valid; i++) {
            aas[i] = translate_codon(s + (c + i) * 3);
            if (aas[i] < 0 || aas[i] == 20) valid = false;
        }
        if (valid) {
            int64_t code = encode_aa6(aas);
            if (code >= 0) {
                cnt.boundary_3p[code]++;
                cnt.total_3p++;
            }
        }
    }

    // 4. Extract RC 6-mers (for strand discrimination training)
    std::string rc;
    reverse_complement(seq, rc);
    const char* r = rc.c_str();

    // Sample RC at same density as coding
    for (size_t c = start_codon; c < end_codon; c++) {
        bool valid = true;
        for (int i = 0; i < 6 && valid; i++) {
            aas[i] = translate_codon(r + (c + i) * 3);
            if (aas[i] < 0) valid = false;  // Allow stops in RC
        }
        if (valid) {
            int64_t code = encode_aa6(aas);
            if (code >= 0) {
                cnt.rc_coding[code]++;
                cnt.total_rc++;
            }
        }
    }
}

// Process a CDS file
void process_file(const std::string& path, Counters& cnt) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) return;

    gzbuffer(f, READ_BUFFER);
    char buf[READ_BUFFER];
    std::string seq;
    seq.reserve(50000);

    while (gzgets(f, buf, sizeof(buf))) {
        size_t len = strlen(buf);
        while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) len--;

        if (buf[0] == '>') {
            if (!seq.empty()) process_gene(seq, cnt);
            seq.clear();
        } else {
            seq.append(buf, len);
        }
    }
    if (!seq.empty()) process_gene(seq, cnt);

    gzclose(f);
}

// Find files
void find_files(const std::string& dir, std::vector<std::string>& files) {
    DIR* d = opendir(dir.c_str());
    if (!d) return;

    struct dirent* e;
    while ((e = readdir(d))) {
        std::string name = e->d_name;
        if (name == "." || name == "..") continue;

        std::string path = dir + "/" + name;
        struct stat st;
        if (stat(path.c_str(), &st) == 0) {
            if (S_ISDIR(st.st_mode)) {
                find_files(path, files);
            } else if (name.find(".fna.gz") != std::string::npos ||
                       name.find(".fa.gz") != std::string::npos) {
                files.push_back(path);
            }
        }
    }
    closedir(d);
}

// Decode AA 6-mer to string
std::string decode_aa6(int64_t code) {
    static const char AA_CHARS[] = "ARNDCQEGHILKMFPSTWYV*";
    std::string s(6, ' ');
    for (int i = 0; i < 6; i++) {
        s[i] = AA_CHARS[code % 21];
        code /= 21;
    }
    return s;
}

int main(int argc, char* argv[]) {
    std::string cds_dir;
    std::string output = "metamers";
    int threads = omp_get_max_threads();

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--cds" && i+1 < argc) cds_dir = argv[++i];
        else if (a == "--output" && i+1 < argc) output = argv[++i];
        else if (a == "--threads" && i+1 < argc) threads = std::stoi(argv[++i]);
        else if (a == "--help") {
            std::cerr << "Usage: " << argv[0] << " --cds <dir> [--output <prefix>] [--threads <n>]\n";
            return 0;
        }
    }

    if (cds_dir.empty()) {
        std::cerr << "Error: --cds required\n";
        return 1;
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    std::cerr << "=== Fast AA 6-mer Extractor ===\n";
    std::cerr << "Threads: " << threads << "\n";
    std::cerr << "Allocating " << (AA_SPACE * 4 * 4 / 1e9) << " GB for counters...\n";

    Counters cnt;

    std::vector<std::string> files;
    find_files(cds_dir, files);
    std::cerr << "Found " << files.size() << " files\n";

    std::atomic<size_t> done{0};
    omp_set_num_threads(threads);

    #pragma omp parallel for schedule(dynamic, 10)
    for (size_t i = 0; i < files.size(); i++) {
        process_file(files[i], cnt);
        size_t d = ++done;
        if (d % 1000 == 0 || d == files.size()) {
            #pragma omp critical
            std::cerr << "\rProcessed " << d << "/" << files.size()
                      << " (" << cnt.genes << " genes)" << std::flush;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    std::cerr << "\n\n=== Results ===\n";
    std::cerr << "Time: " << secs << " seconds\n";
    std::cerr << "Genes: " << cnt.genes << "\n";
    std::cerr << "Coding 6-mers: " << cnt.total_coding << "\n";
    std::cerr << "5' boundary: " << cnt.total_5p << "\n";
    std::cerr << "3' boundary: " << cnt.total_3p << "\n";
    std::cerr << "RC 6-mers: " << cnt.total_rc << "\n";

    // Write output - only non-zero entries
    std::cerr << "\nWriting tables...\n";

    auto write_table = [&](const std::string& name, const std::vector<std::atomic<uint32_t>>& arr, uint64_t total) {
        std::ofstream out(output + "_" + name + ".tsv");
        out << "# " << name << " AA 6-mer frequencies\n";
        out << "# Total: " << total << "\n";

        std::vector<std::pair<int64_t, uint32_t>> entries;
        for (size_t i = 0; i < AA_SPACE; i++) {
            uint32_t c = arr[i];
            if (c > 0) entries.push_back({i, c});
        }

        std::sort(entries.begin(), entries.end(),
                  [](auto& a, auto& b) { return a.second > b.second; });

        size_t n = std::min(entries.size(), (size_t)500000);
        for (size_t i = 0; i < n; i++) {
            out << entries[i].first << "\t" << decode_aa6(entries[i].first)
                << "\t" << entries[i].second << "\t"
                << (double)entries[i].second / total << "\n";
        }
        std::cerr << "  " << name << ": " << entries.size() << " unique, wrote top " << n << "\n";
    };

    write_table("coding", cnt.coding, cnt.total_coding);
    write_table("boundary_5p", cnt.boundary_5p, cnt.total_5p);
    write_table("boundary_3p", cnt.boundary_3p, cnt.total_3p);
    write_table("rc", cnt.rc_coding, cnt.total_rc);

    // Write C++ header with log-likelihood ratios for strand scoring
    std::ofstream hpp(output + "_strand.hpp");
    hpp << R"(#pragma once
// Strand-discriminative AA 6-mer scoring
// log(P(6mer|coding) / P(6mer|RC))

#include <cstdint>
#include <cmath>
#include <array>

namespace agp {

// AA 6-mer to index: code = sum(aa[i] * 21^i) for i=0..5
// AA encoding: A=0,R=1,N=2,D=3,C=4,Q=5,E=6,G=7,H=8,I=9,L=10,K=11,M=12,F=13,P=14,S=15,T=16,W=17,Y=18,V=19,*=20

constexpr size_t AA6_SPACE = 85766121;

inline int8_t translate_codon_fast(const char* s) {
    static const int8_t NT[256] = {['A']=0,['a']=0,['C']=1,['c']=1,['G']=2,['g']=2,['T']=3,['t']=3};
    static const int8_t AA[64] = {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,7,7,7,7,19,19,19,19,20,18,20,18,15,15,15,15,20,4,17,4,10,13,10,13};
    int c = (NT[(uint8_t)s[0]] << 4) | (NT[(uint8_t)s[1]] << 2) | NT[(uint8_t)s[2]];
    return AA[c];
}

constexpr uint32_t POW21[6] = {1, 21, 441, 9261, 194481, 4084101};

inline int64_t encode_aa6_fast(const char* seq) {
    int64_t code = 0;
    for (int i = 0; i < 6; i++) {
        int8_t aa = translate_codon_fast(seq + i*3);
        if (aa < 0 || aa > 20) return -1;
        code += aa * POW21[i];
    }
    return code;
}

)";

    // Compute log-likelihood ratios
    std::vector<float> strand_llr(AA_SPACE, 0.0f);
    double pseudo = 1.0;  // Pseudocount
    for (size_t i = 0; i < AA_SPACE; i++) {
        double p_coding = (cnt.coding[i] + pseudo) / (cnt.total_coding + AA_SPACE * pseudo);
        double p_rc = (cnt.rc_coding[i] + pseudo) / (cnt.total_rc + AA_SPACE * pseudo);
        strand_llr[i] = std::log(p_coding / p_rc);
    }

    // Write top discriminative 6-mers
    std::vector<std::pair<int64_t, float>> sorted_llr;
    for (size_t i = 0; i < AA_SPACE; i++) {
        if (cnt.coding[i] >= 10 || cnt.rc_coding[i] >= 10) {
            sorted_llr.push_back({i, strand_llr[i]});
        }
    }
    std::sort(sorted_llr.begin(), sorted_llr.end(),
              [](auto& a, auto& b) { return std::abs(a.second) > std::abs(b.second); });

    hpp << "// Top strand-discriminative 6-mers (by |LLR|)\n";
    hpp << "// Positive = prefer coding, Negative = prefer RC\n";
    hpp << "inline const std::array<std::pair<int64_t, float>, 100000>& get_strand_llr() {\n";
    hpp << "    static const std::array<std::pair<int64_t, float>, 100000> llr = {{\n";

    size_t n_llr = std::min(sorted_llr.size(), (size_t)100000);
    for (size_t i = 0; i < n_llr; i++) {
        hpp << "        {" << sorted_llr[i].first << ", " << sorted_llr[i].second << "f},\n";
    }
    hpp << "    }};\n    return llr;\n}\n\n";

    hpp << "} // namespace agp\n";

    std::cerr << "Done! Files: " << output << "_*.tsv, " << output << "_strand.hpp\n";

    return 0;
}
