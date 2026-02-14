/*
 * Metamer extractor using Metabuli-style encoding
 *
 * Encoding: 64-bit = (8 AAs × 5 bits) << 24 | (8 codons × 3 bits)
 * Storage: Sorted vectors, merge per thread
 *
 * Compile:
 *   g++ -O3 -march=native -fopenmp -o extract_metamers_v2 extract_metamers_v2.cpp -lz
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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
#include <unordered_map>

constexpr size_t KMER_LEN = 8;
constexpr size_t MIN_GENE_CODONS = 30;
constexpr size_t BOUNDARY_CODONS = 10;
constexpr size_t READ_BUF = 1 << 18;

// Codon to AA (0-20, 20=stop, -1=invalid)
alignas(64) static const int8_t CODON_AA[64] = {
    11, 2, 11, 2, 16, 16, 16, 16, 1, 15, 1, 15, 9, 9, 12, 9,
    5, 8, 5, 8, 14, 14, 14, 14, 1, 1, 1, 1, 10, 10, 10, 10,
    6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7, 19, 19, 19, 19,
    20, 18, 20, 18, 15, 15, 15, 15, 20, 4, 17, 4, 10, 13, 10, 13
};

// Codon to variant index (0-5 within each AA)
alignas(64) static int8_t CODON_VAR[64];

// NT to 2-bit
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

void init_codon_var() {
    int8_t cnt[21] = {0};
    for (int c = 0; c < 64; c++) {
        int8_t aa = CODON_AA[c];
        CODON_VAR[c] = cnt[aa]++;
    }
}

// Encode 8 codons to metamer
inline uint64_t encode_metamer(const char* seq) {
    uint64_t aa_part = 0, dna_part = 0;
    for (int i = 0; i < 8; i++) {
        int b0 = NT[(uint8_t)seq[i*3]];
        int b1 = NT[(uint8_t)seq[i*3+1]];
        int b2 = NT[(uint8_t)seq[i*3+2]];
        if ((b0 | b1 | b2) < 0) return UINT64_MAX;
        int codon = (b0 << 4) | (b1 << 2) | b2;
        int8_t aa = CODON_AA[codon];
        if (aa < 0) return UINT64_MAX;
        int8_t var = CODON_VAR[codon];
        aa_part = (aa_part << 5) | aa;
        dna_part = (dna_part << 3) | var;
    }
    return (aa_part << 24) | dna_part;
}

// Reverse complement
inline char rc_base(char c) {
    switch(c) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        default: return 'N';
    }
}

void reverse_complement(const std::string& s, std::string& rc) {
    rc.resize(s.size());
    for (size_t i = 0; i < s.size(); i++)
        rc[s.size()-1-i] = rc_base(s[i]);
}

// Decode metamer for output
std::string decode_aa(uint64_t m) {
    static const char AA[] = "ARNDCQEGHILKMFPSTWYV*";
    uint64_t aa_part = m >> 24;
    std::string s(8, ' ');
    for (int i = 7; i >= 0; i--) {
        s[i] = AA[aa_part & 0x1F];
        aa_part >>= 5;
    }
    return s;
}

// Thread-local collection
struct ThreadData {
    std::vector<uint64_t> coding;
    std::vector<uint64_t> boundary_5p;
    std::vector<uint64_t> boundary_3p;
    std::vector<uint64_t> rc_coding;

    void reserve(size_t n) {
        coding.reserve(n);
        boundary_5p.reserve(n/10);
        boundary_3p.reserve(n/10);
        rc_coding.reserve(n);
    }
};

std::atomic<uint64_t> genes_done{0};

void process_gene(const std::string& seq, ThreadData& td) {
    size_t len = seq.size();
    size_t n_codons = len / 3;
    if (n_codons < MIN_GENE_CODONS) return;

    genes_done++;
    const char* s = seq.c_str();

    // 1. CODING metamers (interior)
    size_t start = BOUNDARY_CODONS;
    size_t end = n_codons - BOUNDARY_CODONS - 8;
    for (size_t c = start; c < end; c++) {
        uint64_t m = encode_metamer(s + c*3);
        if (m != UINT64_MAX && (m >> 24) != 0x108421084210ULL) // Skip all-stop
            td.coding.push_back(m);
    }

    // 2. 5' BOUNDARY
    for (size_t c = 0; c < BOUNDARY_CODONS && c+8 <= n_codons; c++) {
        uint64_t m = encode_metamer(s + c*3);
        if (m != UINT64_MAX) td.boundary_5p.push_back(m);
    }

    // 3. 3' BOUNDARY
    size_t s3 = (n_codons > BOUNDARY_CODONS + 8) ? n_codons - BOUNDARY_CODONS - 8 : 0;
    for (size_t c = s3; c + 8 <= n_codons; c++) {
        uint64_t m = encode_metamer(s + c*3);
        if (m != UINT64_MAX) td.boundary_3p.push_back(m);
    }

    // 4. RC metamers (for strand discrimination)
    std::string rc;
    reverse_complement(seq, rc);
    const char* r = rc.c_str();
    for (size_t c = start; c < end; c++) {
        uint64_t m = encode_metamer(r + c*3);
        if (m != UINT64_MAX) td.rc_coding.push_back(m);
    }
}

void process_file(const std::string& path, ThreadData& td) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) return;
    gzbuffer(f, READ_BUF);

    char buf[READ_BUF];
    std::string seq;
    seq.reserve(50000);

    while (gzgets(f, buf, sizeof(buf))) {
        size_t len = strlen(buf);
        while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) len--;
        if (buf[0] == '>') {
            if (!seq.empty()) process_gene(seq, td);
            seq.clear();
        } else {
            seq.append(buf, len);
        }
    }
    if (!seq.empty()) process_gene(seq, td);
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
            else if (name.find(".fna.gz") != std::string::npos ||
                     name.find(".fa.gz") != std::string::npos)
                files.push_back(path);
        }
    }
    closedir(d);
}

// Count frequencies from sorted vector
std::unordered_map<uint64_t, uint32_t> count_freq(std::vector<uint64_t>& v) {
    std::sort(v.begin(), v.end());
    std::unordered_map<uint64_t, uint32_t> freq;
    freq.reserve(v.size() / 10);

    uint64_t prev = UINT64_MAX;
    uint32_t cnt = 0;
    for (uint64_t m : v) {
        if (m == prev) cnt++;
        else {
            if (prev != UINT64_MAX && cnt > 0) freq[prev] = cnt;
            prev = m;
            cnt = 1;
        }
    }
    if (prev != UINT64_MAX && cnt > 0) freq[prev] = cnt;
    return freq;
}

int main(int argc, char* argv[]) {
    std::string cds_dir, output = "metamers";
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

    init_codon_var();
    auto t0 = std::chrono::high_resolution_clock::now();

    std::cerr << "=== Metabuli-style Metamer Extractor ===\n";
    std::cerr << "Threads: " << threads << "\n";

    std::vector<std::string> files;
    find_files(cds_dir, files);
    std::cerr << "Found " << files.size() << " files\n";

    // Thread-local data
    std::vector<ThreadData> tdata(threads);
    for (auto& td : tdata) td.reserve(10000000);

    std::atomic<size_t> done{0};
    omp_set_num_threads(threads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < files.size(); i++) {
            process_file(files[i], tdata[tid]);
            size_t d = ++done;
            if (d % 500 == 0 || d == files.size()) {
                #pragma omp critical
                std::cerr << "\rProcessed " << d << "/" << files.size()
                          << " (" << genes_done << " genes)" << std::flush;
            }
        }
    }

    std::cerr << "\n\nMerging thread data...\n";

    // Merge all thread vectors
    std::vector<uint64_t> all_coding, all_5p, all_3p, all_rc;
    for (auto& td : tdata) {
        all_coding.insert(all_coding.end(), td.coding.begin(), td.coding.end());
        all_5p.insert(all_5p.end(), td.boundary_5p.begin(), td.boundary_5p.end());
        all_3p.insert(all_3p.end(), td.boundary_3p.begin(), td.boundary_3p.end());
        all_rc.insert(all_rc.end(), td.rc_coding.begin(), td.rc_coding.end());
        td.coding.clear(); td.coding.shrink_to_fit();
        td.boundary_5p.clear(); td.boundary_5p.shrink_to_fit();
        td.boundary_3p.clear(); td.boundary_3p.shrink_to_fit();
        td.rc_coding.clear(); td.rc_coding.shrink_to_fit();
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    std::cerr << "\n=== Results ===\n";
    std::cerr << "Time: " << secs << " s\n";
    std::cerr << "Genes: " << genes_done << "\n";
    std::cerr << "Coding: " << all_coding.size() << "\n";
    std::cerr << "5' boundary: " << all_5p.size() << "\n";
    std::cerr << "3' boundary: " << all_3p.size() << "\n";
    std::cerr << "RC: " << all_rc.size() << "\n";

    // Count frequencies
    std::cerr << "\nCounting frequencies...\n";
    auto freq_coding = count_freq(all_coding);
    auto freq_5p = count_freq(all_5p);
    auto freq_3p = count_freq(all_3p);
    auto freq_rc = count_freq(all_rc);

    std::cerr << "Unique coding: " << freq_coding.size() << "\n";
    std::cerr << "Unique 5': " << freq_5p.size() << "\n";
    std::cerr << "Unique 3': " << freq_3p.size() << "\n";
    std::cerr << "Unique RC: " << freq_rc.size() << "\n";

    // Write TSV output
    auto write_tsv = [](const std::string& path, const std::unordered_map<uint64_t, uint32_t>& freq, uint64_t total) {
        std::vector<std::pair<uint64_t, uint32_t>> sorted(freq.begin(), freq.end());
        std::sort(sorted.begin(), sorted.end(), [](auto& a, auto& b) { return a.second > b.second; });

        std::ofstream out(path);
        out << "# Total: " << total << "\n# Unique: " << freq.size() << "\n";
        size_t n = std::min(sorted.size(), (size_t)500000);
        for (size_t i = 0; i < n; i++) {
            out << sorted[i].first << "\t" << decode_aa(sorted[i].first)
                << "\t" << sorted[i].second << "\t"
                << (double)sorted[i].second / total << "\n";
        }
    };

    std::cerr << "Writing output...\n";
    write_tsv(output + "_coding.tsv", freq_coding, all_coding.size());
    write_tsv(output + "_5p.tsv", freq_5p, all_5p.size());
    write_tsv(output + "_3p.tsv", freq_3p, all_3p.size());
    write_tsv(output + "_rc.tsv", freq_rc, all_rc.size());

    // Write C++ header with strand LLR
    std::ofstream hpp(output + "_strand.hpp");
    hpp << "#pragma once\n// Strand-discriminative 8-mer metamer scoring\n\n";
    hpp << "#include <cstdint>\n#include <unordered_map>\n#include <cmath>\n\n";
    hpp << "namespace agp {\n\n";

    // Compute LLR
    double pseudo = 0.5;
    double total_coding = all_coding.size();
    double total_rc = all_rc.size();

    std::vector<std::pair<uint64_t, float>> llr_entries;
    for (auto& [m, c] : freq_coding) {
        double p_c = (c + pseudo) / (total_coding + pseudo * freq_coding.size());
        double rc_cnt = freq_rc.count(m) ? freq_rc[m] : 0;
        double p_r = (rc_cnt + pseudo) / (total_rc + pseudo * freq_rc.size());
        float llr = std::log(p_c / p_r);
        if (std::abs(llr) > 0.1f) llr_entries.push_back({m, llr});
    }
    std::sort(llr_entries.begin(), llr_entries.end(),
              [](auto& a, auto& b) { return std::abs(a.second) > std::abs(b.second); });

    size_t n_llr = std::min(llr_entries.size(), (size_t)100000);
    hpp << "inline const std::unordered_map<uint64_t, float>& get_strand_llr() {\n";
    hpp << "    static const std::unordered_map<uint64_t, float> llr = {\n";
    for (size_t i = 0; i < n_llr; i++) {
        hpp << "        {" << llr_entries[i].first << "ULL, " << llr_entries[i].second << "f},\n";
    }
    hpp << "    };\n    return llr;\n}\n\n} // namespace agp\n";

    std::cerr << "Done! Wrote " << n_llr << " strand-discriminative metamers\n";

    return 0;
}
