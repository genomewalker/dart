/*
 * ULTRA-FAST 6-mer AA extractor using direct array indexing
 *
 * Key optimizations:
 * 1. 6-mer AA = 21^6 = 85,766,121 values → direct array (340MB per category)
 * 2. Atomic increments - no locks, no merging
 * 3. Minimal memory allocations
 *
 * Compile:
 *   g++ -O3 -march=native -fopenmp -o extract_aa6_ultra extract_aa6_ultra.cpp -lz
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <dirent.h>
#include <zlib.h>
#include <omp.h>
#include <sys/stat.h>
#include <chrono>
#include <algorithm>

constexpr uint32_t AA6_SIZE = 85766121;  // 21^6
constexpr uint32_t POW21[7] = {1, 21, 441, 9261, 194481, 4084101, 85766121};
constexpr size_t MIN_CODONS = 20;
constexpr size_t BOUNDARY = 8;  // codons

// Codon to AA
alignas(64) static const int8_t CODON_AA[64] = {
    11, 2, 11, 2, 16, 16, 16, 16, 1, 15, 1, 15, 9, 9, 12, 9,
    5, 8, 5, 8, 14, 14, 14, 14, 1, 1, 1, 1, 10, 10, 10, 10,
    6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7, 19, 19, 19, 19,
    20, 18, 20, 18, 15, 15, 15, 15, 20, 4, 17, 4, 10, 13, 10, 13
};

// NT to 2-bit (A=0,C=1,G=2,T=3)
alignas(64) static const int8_t NT[256] = {
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

// RC base
inline char rc_base(char c) {
    switch(c) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        default: return 'N';
    }
}

// Translate codon to AA index (0-20, -1 for invalid)
inline int8_t translate(const char* s) {
    int b0 = NT[(uint8_t)s[0]], b1 = NT[(uint8_t)s[1]], b2 = NT[(uint8_t)s[2]];
    if ((b0 | b1 | b2) < 0) return -1;
    return CODON_AA[(b0 << 4) | (b1 << 2) | b2];
}

// Encode 6-mer AA to index
inline uint32_t encode6(const int8_t* aa) {
    return aa[0] + aa[1]*21 + aa[2]*441 + aa[3]*9261 + aa[4]*194481 + aa[5]*4084101;
}

// Decode to string
std::string decode6(uint32_t code) {
    static const char AA[] = "ARNDCQEGHILKMFPSTWYV*";
    std::string s(6, ' ');
    for (int i = 0; i < 6; i++) { s[i] = AA[code % 21]; code /= 21; }
    return s;
}

// Global atomic counters
struct Counters {
    std::atomic<uint32_t>* coding;
    std::atomic<uint32_t>* boundary_5p;
    std::atomic<uint32_t>* boundary_3p;
    std::atomic<uint32_t>* rc;

    std::atomic<uint64_t> total_coding{0};
    std::atomic<uint64_t> total_5p{0};
    std::atomic<uint64_t> total_3p{0};
    std::atomic<uint64_t> total_rc{0};
    std::atomic<uint64_t> genes{0};

    Counters() {
        coding = new std::atomic<uint32_t>[AA6_SIZE]();
        boundary_5p = new std::atomic<uint32_t>[AA6_SIZE]();
        boundary_3p = new std::atomic<uint32_t>[AA6_SIZE]();
        rc = new std::atomic<uint32_t>[AA6_SIZE]();
    }

    ~Counters() {
        delete[] coding;
        delete[] boundary_5p;
        delete[] boundary_3p;
        delete[] rc;
    }
};

void process_gene(const char* seq, size_t len, Counters& c) {
    size_t n = len / 3;
    if (n < MIN_CODONS) return;
    c.genes++;

    int8_t aa[6];
    auto extract = [&](size_t start, size_t end, std::atomic<uint32_t>* arr, std::atomic<uint64_t>& total) {
        for (size_t i = start; i + 6 <= end; i++) {
            bool ok = true;
            for (int j = 0; j < 6 && ok; j++) {
                aa[j] = translate(seq + (i+j)*3);
                if (aa[j] < 0 || aa[j] == 20) ok = false;  // Skip invalid/stop
            }
            if (ok) {
                uint32_t code = encode6(aa);
                arr[code].fetch_add(1, std::memory_order_relaxed);
                total.fetch_add(1, std::memory_order_relaxed);
            }
        }
    };

    // Coding (interior)
    if (n > 2*BOUNDARY + 6) extract(BOUNDARY, n - BOUNDARY, c.coding, c.total_coding);

    // 5' boundary
    extract(0, std::min(n, BOUNDARY + 6), c.boundary_5p, c.total_5p);

    // 3' boundary
    size_t s3 = (n > BOUNDARY + 6) ? n - BOUNDARY - 6 : 0;
    extract(s3, n, c.boundary_3p, c.total_3p);

    // RC
    thread_local std::string rc_buf;
    rc_buf.resize(len);
    for (size_t i = 0; i < len; i++) rc_buf[len-1-i] = RC[(uint8_t)seq[i]];

    if (n > 2*BOUNDARY + 6) {
        const char* r = rc_buf.c_str();
        for (size_t i = BOUNDARY; i + 6 <= n - BOUNDARY; i++) {
            bool ok = true;
            for (int j = 0; j < 6 && ok; j++) {
                aa[j] = translate(r + (i+j)*3);
                if (aa[j] < 0) ok = false;  // Allow stops in RC
            }
            if (ok) {
                uint32_t code = encode6(aa);
                c.rc[code].fetch_add(1, std::memory_order_relaxed);
                c.total_rc.fetch_add(1, std::memory_order_relaxed);
            }
        }
    }
}

void process_file(const std::string& path, Counters& c) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) return;
    gzbuffer(f, 1 << 18);

    char buf[1 << 18];
    thread_local std::string seq;
    seq.clear();
    seq.reserve(50000);

    while (gzgets(f, buf, sizeof(buf))) {
        size_t len = strlen(buf);
        while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) len--;
        if (buf[0] == '>') {
            if (!seq.empty()) process_gene(seq.c_str(), seq.size(), c);
            seq.clear();
        } else {
            seq.append(buf, len);
        }
    }
    if (!seq.empty()) process_gene(seq.c_str(), seq.size(), c);
    gzclose(f);
}

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
            if (S_ISDIR(st.st_mode)) find_files(path, files);
            else if (name.find(".fna") != std::string::npos || name.find(".fa") != std::string::npos)
                files.push_back(path);
        }
    }
    closedir(d);
}

int main(int argc, char* argv[]) {
    std::string cds_dir, output = "aa6";
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
    if (cds_dir.empty()) { std::cerr << "Error: --cds required\n"; return 1; }

    auto t0 = std::chrono::high_resolution_clock::now();

    std::cerr << "=== 6-mer AA Extractor ===\n";
    std::cerr << "Allocating " << (AA6_SIZE * 4 * 4 / 1e9) << " GB for counters...\n";

    Counters cnt;

    std::vector<std::string> files;
    find_files(cds_dir, files);
    std::cerr << "Found " << files.size() << " files\n";
    std::cerr << "Using " << threads << " threads\n";

    std::atomic<size_t> done{0};
    omp_set_num_threads(threads);

    #pragma omp parallel for schedule(dynamic, 10)
    for (size_t i = 0; i < files.size(); i++) {
        process_file(files[i], cnt);
        size_t d = ++done;
        if (d % 500 == 0 || d == files.size()) {
            #pragma omp critical
            std::cerr << "\r" << d << "/" << files.size() << " files, " << cnt.genes << " genes" << std::flush;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    std::cerr << "\n\n=== Results ===\n";
    std::cerr << "Time: " << secs << "s (" << (files.size() / secs) << " files/sec)\n";
    std::cerr << "Genes: " << cnt.genes << "\n";
    std::cerr << "Coding: " << cnt.total_coding << "\n";
    std::cerr << "5' boundary: " << cnt.total_5p << "\n";
    std::cerr << "3' boundary: " << cnt.total_3p << "\n";
    std::cerr << "RC: " << cnt.total_rc << "\n";

    // Count non-zero
    size_t nz_coding = 0, nz_rc = 0;
    for (size_t i = 0; i < AA6_SIZE; i++) {
        if (cnt.coding[i] > 0) nz_coding++;
        if (cnt.rc[i] > 0) nz_rc++;
    }
    std::cerr << "Unique coding: " << nz_coding << " / " << AA6_SIZE << "\n";
    std::cerr << "Unique RC: " << nz_rc << "\n";

    // Write TSV (top 500k)
    std::cerr << "\nWriting output...\n";

    auto write_tsv = [](const std::string& path, std::atomic<uint32_t>* arr, uint64_t total) {
        std::vector<std::pair<uint32_t, uint32_t>> entries;
        entries.reserve(AA6_SIZE / 10);
        for (uint32_t i = 0; i < AA6_SIZE; i++) {
            uint32_t c = arr[i].load(std::memory_order_relaxed);
            if (c > 0) entries.push_back({i, c});
        }
        std::partial_sort(entries.begin(),
                          entries.begin() + std::min(entries.size(), (size_t)500000),
                          entries.end(),
                          [](auto& a, auto& b) { return a.second > b.second; });

        std::ofstream out(path);
        out << "# Total: " << total << "\n# Unique: " << entries.size() << "\n";
        size_t n = std::min(entries.size(), (size_t)500000);
        for (size_t i = 0; i < n; i++) {
            out << entries[i].first << "\t" << decode6(entries[i].first)
                << "\t" << entries[i].second << "\t"
                << (double)entries[i].second / total << "\n";
        }
        std::cerr << "  " << path << ": " << n << " entries\n";
    };

    write_tsv(output + "_coding.tsv", cnt.coding, cnt.total_coding);
    write_tsv(output + "_5p.tsv", cnt.boundary_5p, cnt.total_5p);
    write_tsv(output + "_3p.tsv", cnt.boundary_3p, cnt.total_3p);
    write_tsv(output + "_rc.tsv", cnt.rc, cnt.total_rc);

    // Write C++ header with strand LLR
    std::cerr << "Computing strand LLR...\n";
    std::ofstream hpp(output + "_strand.hpp");
    hpp << "#pragma once\n// 6-mer AA strand discrimination\n\n";
    hpp << "#include <cstdint>\n#include <array>\n#include <cmath>\n\n";
    hpp << "namespace agp {\n\n";
    hpp << "constexpr uint32_t AA6_SIZE = " << AA6_SIZE << ";\n";
    hpp << "constexpr uint32_t POW21[6] = {1, 21, 441, 9261, 194481, 4084101};\n\n";

    // Compute LLR array (only store significant ones)
    std::vector<std::pair<uint32_t, float>> llr_list;
    double pseudo = 1.0;
    double tc = cnt.total_coding + AA6_SIZE * pseudo;
    double tr = cnt.total_rc + AA6_SIZE * pseudo;

    for (uint32_t i = 0; i < AA6_SIZE; i++) {
        uint32_t cc = cnt.coding[i], cr = cnt.rc[i];
        if (cc > 5 || cr > 5) {
            double pc = (cc + pseudo) / tc;
            double pr = (cr + pseudo) / tr;
            float llr = std::log(pc / pr);
            if (std::abs(llr) > 0.2f) llr_list.push_back({i, llr});
        }
    }

    std::sort(llr_list.begin(), llr_list.end(),
              [](auto& a, auto& b) { return std::abs(a.second) > std::abs(b.second); });

    size_t n_llr = std::min(llr_list.size(), (size_t)100000);
    hpp << "// Top " << n_llr << " strand-discriminative 6-mers\n";
    hpp << "inline const std::array<std::pair<uint32_t, float>, " << n_llr << ">& get_strand_llr() {\n";
    hpp << "    static const std::array<std::pair<uint32_t, float>, " << n_llr << "> llr = {{\n";
    for (size_t i = 0; i < n_llr; i++) {
        hpp << "        {" << llr_list[i].first << ", " << llr_list[i].second << "f},\n";
    }
    hpp << "    }};\n    return llr;\n}\n\n";
    hpp << "} // namespace agp\n";

    std::cerr << "Done! " << n_llr << " strand-discriminative 6-mers\n";

    return 0;
}
