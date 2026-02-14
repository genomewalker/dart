/*
 * Metamer extractor for strand-discriminative training
 *
 * Metamers = 8-mer amino acids + codon encoding (like Metabuli)
 * This captures both AA conservation AND strand-specific codon usage.
 *
 * Categories extracted:
 * 1. CODING - interior of genes (intra-gene)
 * 2. BOUNDARY_5P - near start codon (gene start region)
 * 3. BOUNDARY_3P - near stop codon (gene end region)
 * 4. INTERGENIC - between genes (NT k-mers only)
 *
 * Output: frequency tables for strand-discriminative scoring
 *
 * Compile:
 *   g++ -O3 -march=native -fopenmp -o extract_metamers extract_metamers.cpp -lz -lpthread
 *
 * Usage:
 *   ./extract_metamers --cds <protein_fna_dir> --genomes <genome_dir> --output <prefix>
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <atomic>
#include <thread>
#include <mutex>
#include <dirent.h>
#include <zlib.h>
#include <omp.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

// Constants
constexpr size_t KMER_AA = 6;           // 6-mer amino acids (21^6 = 85M possible)
constexpr size_t KMER_NT = 18;          // 18-mer nucleotides (6 codons)
constexpr size_t BOUNDARY_WINDOW = 24;  // bases around start/stop (8 codons)
constexpr size_t MIN_GENE_LEN = 60;     // minimum gene length (20 codons)
constexpr size_t INTERGENIC_KMER = 12;  // 12-mer for intergenic regions
constexpr size_t READ_BUFFER = 1 << 18; // 256KB read buffer

// Encoding tables

// Nucleotide to 2-bit: A=0, C=1, G=2, T=3
alignas(64) static const int8_t NT_MAP[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1, // A=0, C=1, G=2
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // T=3
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1, // a=0, c=1, g=2
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // t=3
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
};

// Codon to amino acid (0-20, where 20 = stop, 21 = invalid)
// Standard genetic code
alignas(64) static const uint8_t CODON_TO_AA[64] = {
    // AAA=K, AAC=N, AAG=K, AAT=N, ACA=T, ACC=T, ACG=T, ACT=T
    11, 2, 11, 2, 16, 16, 16, 16,
    // AGA=R, AGC=S, AGG=R, AGT=S, ATA=I, ATC=I, ATG=M, ATT=I
    1, 15, 1, 15, 9, 9, 12, 9,
    // CAA=Q, CAC=H, CAG=Q, CAT=H, CCA=P, CCC=P, CCG=P, CCT=P
    5, 8, 5, 8, 14, 14, 14, 14,
    // CGA=R, CGC=R, CGG=R, CGT=R, CTA=L, CTC=L, CTG=L, CTT=L
    1, 1, 1, 1, 10, 10, 10, 10,
    // GAA=E, GAC=D, GAG=E, GAT=D, GCA=A, GCC=A, GCG=A, GCT=A
    6, 3, 6, 3, 0, 0, 0, 0,
    // GGA=G, GGC=G, GGG=G, GGT=G, GTA=V, GTC=V, GTG=V, GTT=V
    7, 7, 7, 7, 19, 19, 19, 19,
    // TAA=*, TAC=Y, TAG=*, TAT=Y, TCA=S, TCC=S, TCG=S, TCT=S
    20, 18, 20, 18, 15, 15, 15, 15,
    // TGA=*, TGC=C, TGG=W, TGT=C, TTA=L, TTC=F, TTG=L, TTT=F
    20, 4, 17, 4, 10, 13, 10, 13
};

// Amino acid to number of synonymous codons (for codon variant encoding)
alignas(64) static const uint8_t AA_CODON_COUNT[21] = {
    4, // A - Ala (GCN)
    6, // R - Arg (CGN, AGR)
    2, // N - Asn (AAY)
    2, // D - Asp (GAY)
    2, // C - Cys (TGY)
    2, // Q - Gln (CAR)
    2, // E - Glu (GAR)
    4, // G - Gly (GGN)
    2, // H - His (CAY)
    3, // I - Ile (ATH)
    6, // L - Leu (CTN, TTR)
    2, // K - Lys (AAR)
    1, // M - Met (ATG)
    2, // F - Phe (TTY)
    4, // P - Pro (CCN)
    6, // S - Ser (TCN, AGY)
    4, // T - Thr (ACN)
    1, // W - Trp (TGG)
    2, // Y - Tyr (TAY)
    4, // V - Val (GTN)
    3  // * - Stop (TAR, TGA)
};

// Codon to variant index within amino acid
alignas(64) static uint8_t CODON_TO_VARIANT[64];

// Initialize codon variant table
void init_codon_variants() {
    uint8_t aa_variant_count[22] = {0};
    for (int c = 0; c < 64; c++) {
        uint8_t aa = CODON_TO_AA[c];
        CODON_TO_VARIANT[c] = aa_variant_count[aa]++;
    }
}

// Metamer encoding

// Metamer: 8 amino acids + codon variants
// Encoding:
//   - AA part: 8 × 5 bits = 40 bits (21 values per position)
//   - Codon part: use separate 64-bit for codon variants
struct Metamer {
    uint64_t aa_code;     // 8 AAs encoded in 40 bits
    uint64_t codon_code;  // 8 codon variants (3 bits each = 24 bits)

    bool operator==(const Metamer& other) const {
        return aa_code == other.aa_code && codon_code == other.codon_code;
    }
};

struct MetamerHash {
    size_t operator()(const Metamer& m) const {
        return m.aa_code ^ (m.codon_code * 0x9e3779b97f4a7c15ULL);
    }
};

// Encode codon (3 nucleotides) to 6-bit code
inline int encode_codon(const char* seq) {
    int8_t b0 = NT_MAP[(uint8_t)seq[0]];
    int8_t b1 = NT_MAP[(uint8_t)seq[1]];
    int8_t b2 = NT_MAP[(uint8_t)seq[2]];
    if ((b0 | b1 | b2) < 0) return -1;
    return (b0 << 4) | (b1 << 2) | b2;
}

// Encode 8 codons (24 nt) to metamer
inline bool encode_metamer(const char* seq, Metamer& m) {
    m.aa_code = 0;
    m.codon_code = 0;

    for (int i = 0; i < 8; i++) {
        int codon = encode_codon(seq + i * 3);
        if (codon < 0) return false;

        uint8_t aa = CODON_TO_AA[codon];
        if (aa > 20) return false;  // Invalid

        uint8_t variant = CODON_TO_VARIANT[codon];

        m.aa_code |= (uint64_t)aa << (i * 5);
        m.codon_code |= (uint64_t)variant << (i * 3);
    }
    return true;
}

// Encode nucleotide k-mer (for intergenic)
inline uint64_t encode_nt_kmer(const char* seq, int k) {
    uint64_t code = 0;
    for (int i = 0; i < k; i++) {
        int8_t b = NT_MAP[(uint8_t)seq[i]];
        if (b < 0) return UINT64_MAX;
        code = (code << 2) | b;
    }
    return code;
}

// Thread-safe counters

struct MetamerCounts {
    std::unordered_map<Metamer, uint64_t, MetamerHash> coding;
    std::unordered_map<Metamer, uint64_t, MetamerHash> boundary_5p;
    std::unordered_map<Metamer, uint64_t, MetamerHash> boundary_3p;
    std::vector<std::atomic<uint64_t>> intergenic_12mer;  // 4^12 = 16M (heap allocated)

    std::atomic<uint64_t> total_coding{0};
    std::atomic<uint64_t> total_5p{0};
    std::atomic<uint64_t> total_3p{0};
    std::atomic<uint64_t> total_intergenic{0};
    std::atomic<uint64_t> genes_processed{0};

    std::mutex mtx;  // For hash map updates

    MetamerCounts() : intergenic_12mer(1ULL << 24) {
        for (auto& a : intergenic_12mer) a = 0;
    }

    void add_coding(const Metamer& m) {
        std::lock_guard<std::mutex> lock(mtx);
        coding[m]++;
        total_coding++;
    }

    void add_5p(const Metamer& m) {
        std::lock_guard<std::mutex> lock(mtx);
        boundary_5p[m]++;
        total_5p++;
    }

    void add_3p(const Metamer& m) {
        std::lock_guard<std::mutex> lock(mtx);
        boundary_3p[m]++;
        total_3p++;
    }

    void add_intergenic(uint64_t code) {
        if (code < intergenic_12mer.size()) {
            intergenic_12mer[code]++;
            total_intergenic++;
        }
    }
};

// Thread-local batch for reducing lock contention
struct ThreadBatch {
    std::unordered_map<Metamer, uint64_t, MetamerHash> coding;
    std::unordered_map<Metamer, uint64_t, MetamerHash> boundary_5p;
    std::unordered_map<Metamer, uint64_t, MetamerHash> boundary_3p;
    // Skip intergenic regions and focus on coding/boundary classes

    uint64_t total_coding = 0;
    uint64_t total_5p = 0;
    uint64_t total_3p = 0;

    ThreadBatch() {}

    void flush(MetamerCounts& global) {
        {
            std::lock_guard<std::mutex> lock(global.mtx);
            for (auto& [m, c] : coding) global.coding[m] += c;
            for (auto& [m, c] : boundary_5p) global.boundary_5p[m] += c;
            for (auto& [m, c] : boundary_3p) global.boundary_3p[m] += c;
        }
        global.total_coding += total_coding;
        global.total_5p += total_5p;
        global.total_3p += total_3p;

        // Clear batch
        coding.clear();
        boundary_5p.clear();
        boundary_3p.clear();
        total_coding = total_5p = total_3p = 0;
    }
};

// Gene parsing

struct Gene {
    std::string contig;
    size_t start;
    size_t end;
    int strand;  // 1 = forward, -1 = reverse
    std::string sequence;
};

// Parse Prodigal header: >contig_id # start # end # strand # ...
bool parse_prodigal_header(const std::string& header, Gene& gene) {
    // Format: >BA000021.3_1 # 182 # 2068 # 1 # ID=...
    size_t pos = 0;

    // Skip '>'
    if (header[0] == '>') pos = 1;

    // Get contig_gene
    size_t space = header.find(' ', pos);
    if (space == std::string::npos) return false;
    gene.contig = header.substr(pos, space - pos);

    // Parse # start # end # strand
    int field = 0;
    size_t i = space;
    while (i < header.size() && field < 3) {
        if (header[i] == '#') {
            i++;
            while (i < header.size() && header[i] == ' ') i++;
            size_t num_start = i;
            while (i < header.size() && header[i] != ' ' && header[i] != '#') i++;
            std::string num = header.substr(num_start, i - num_start);

            switch (field) {
                case 0: gene.start = std::stoull(num); break;
                case 1: gene.end = std::stoull(num); break;
                case 2: gene.strand = std::stoi(num); break;
            }
            field++;
        } else {
            i++;
        }
    }

    return field >= 3;
}

// File processing

void process_cds_file(const std::string& filepath, ThreadBatch& batch, MetamerCounts& global) {
    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) return;

    gzbuffer(file, READ_BUFFER);

    char buffer[READ_BUFFER];
    Gene current_gene;
    std::string current_seq;
    current_seq.reserve(100000);

    auto process_gene = [&]() {
        if (current_seq.size() < MIN_GENE_LEN) return;

        const char* seq = current_seq.c_str();
        size_t len = current_seq.size();

        // Ensure length is multiple of 3
        len = (len / 3) * 3;
        if (len < KMER_NT) return;

        global.genes_processed++;

        // Extract CODING metamers (interior, skip first and last 30bp)
        size_t coding_start = BOUNDARY_WINDOW;
        size_t coding_end = len - BOUNDARY_WINDOW - KMER_NT;

        if (coding_end > coding_start) {
            for (size_t i = coding_start; i <= coding_end; i += 3) {
                Metamer m;
                if (encode_metamer(seq + i, m)) {
                    batch.coding[m]++;
                    batch.total_coding++;
                }
            }
        }

        // Extract 5' BOUNDARY metamers (first 30bp region)
        for (size_t i = 0; i + KMER_NT <= BOUNDARY_WINDOW + KMER_NT && i + KMER_NT <= len; i += 3) {
            Metamer m;
            if (encode_metamer(seq + i, m)) {
                batch.boundary_5p[m]++;
                batch.total_5p++;
            }
        }

        // Extract 3' BOUNDARY metamers (last 30bp region)
        size_t start_3p = (len > BOUNDARY_WINDOW + KMER_NT) ? len - BOUNDARY_WINDOW - KMER_NT : 0;
        start_3p = (start_3p / 3) * 3;  // Align to codon
        for (size_t i = start_3p; i + KMER_NT <= len; i += 3) {
            Metamer m;
            if (encode_metamer(seq + i, m)) {
                batch.boundary_3p[m]++;
                batch.total_3p++;
            }
        }
    };

    while (gzgets(file, buffer, sizeof(buffer)) != nullptr) {
        size_t slen = strlen(buffer);
        while (slen > 0 && (buffer[slen-1] == '\n' || buffer[slen-1] == '\r')) slen--;
        buffer[slen] = '\0';

        if (buffer[0] == '>') {
            process_gene();
            current_seq.clear();
            parse_prodigal_header(buffer, current_gene);
        } else {
            current_seq.append(buffer, slen);
        }
    }
    process_gene();

    gzclose(file);

    // Flush batch periodically
    if (batch.total_coding > 1000000) {
        batch.flush(global);
    }
}

// Find all files in directory
void find_files(const std::string& dir, std::vector<std::string>& files, const std::string& suffix = ".fna.gz") {
    DIR* d = opendir(dir.c_str());
    if (!d) return;

    struct dirent* entry;
    while ((entry = readdir(d)) != nullptr) {
        std::string name = entry->d_name;
        if (name == "." || name == "..") continue;

        std::string path = dir + "/" + name;

        struct stat st;
        if (stat(path.c_str(), &st) == 0) {
            if (S_ISDIR(st.st_mode)) {
                find_files(path, files, suffix);
            } else if (name.size() >= suffix.size() &&
                       name.substr(name.size() - suffix.size()) == suffix) {
                files.push_back(path);
            }
        }
    }
    closedir(d);
}

// Output

void write_metamer_table(const std::string& filename,
                         const std::unordered_map<Metamer, uint64_t, MetamerHash>& counts,
                         uint64_t total, const std::string& name,
                         uint64_t min_count = 5, size_t max_entries = 500000) {
    std::ofstream out(filename);
    out << "# " << name << " metamer frequencies\n";
    out << "# Total: " << total << "\n";
    out << "# Unique: " << counts.size() << "\n";
    out << "# Min count filter: " << min_count << "\n";
    out << "# Format: aa_code codon_code count frequency\n";

    // Filter and collect
    std::vector<std::pair<Metamer, uint64_t>> filtered;
    filtered.reserve(std::min(counts.size(), max_entries));
    for (const auto& [m, c] : counts) {
        if (c >= min_count) filtered.push_back({m, c});
    }

    // Partial sort to get top entries
    size_t n = std::min(filtered.size(), max_entries);
    std::partial_sort(filtered.begin(), filtered.begin() + n, filtered.end(),
                      [](const auto& a, const auto& b) { return a.second > b.second; });

    for (size_t i = 0; i < n; i++) {
        double freq = (double)filtered[i].second / total;
        out << filtered[i].first.aa_code << "\t" << filtered[i].first.codon_code
            << "\t" << filtered[i].second << "\t" << freq << "\n";
    }
    std::cerr << "  Wrote " << n << " " << name << " metamers (filtered from " << counts.size() << ")\n";
}

void write_cpp_header(const std::string& filename,
                      const std::unordered_map<Metamer, uint64_t, MetamerHash>& coding,
                      const std::unordered_map<Metamer, uint64_t, MetamerHash>& boundary_5p,
                      const std::unordered_map<Metamer, uint64_t, MetamerHash>& boundary_3p,
                      uint64_t total_coding, uint64_t total_5p, uint64_t total_3p) {
    std::ofstream out(filename);

    out << R"(#pragma once
/**
 * @file metamer_tables.hpp
 * @brief Metamer frequency tables for strand-discriminative scoring
 *
 * Generated by extract_metamers from GTDB data.
 *
 * Metamers encode 8 amino acids + codon variants:
 * - aa_code: 8 × 5 bits = 40 bits
 * - codon_code: 8 × 3 bits = 24 bits
 *
 * Categories:
 * - CODING: interior of genes
 * - BOUNDARY_5P: near start codon (first 30bp)
 * - BOUNDARY_3P: near stop codon (last 30bp)
 */

#include <cstdint>
#include <unordered_map>

namespace agp {

struct Metamer {
    uint64_t aa_code;
    uint64_t codon_code;

    bool operator==(const Metamer& other) const {
        return aa_code == other.aa_code && codon_code == other.codon_code;
    }
};

struct MetamerHash {
    size_t operator()(const Metamer& m) const {
        return m.aa_code ^ (m.codon_code * 0x9e3779b97f4a7c15ULL);
    }
};

)";

    // Write top metamers as static arrays for fast lookup
    auto write_top_metamers = [&](const std::string& name,
                                   const std::unordered_map<Metamer, uint64_t, MetamerHash>& counts,
                                   uint64_t total, size_t max_entries = 100000) {
        std::vector<std::pair<Metamer, uint64_t>> sorted(counts.begin(), counts.end());
        std::sort(sorted.begin(), sorted.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });

        size_t n = std::min(sorted.size(), max_entries);

        out << "// " << name << " metamers (top " << n << " of " << counts.size() << ")\n";
        out << "inline const std::unordered_map<Metamer, float, MetamerHash>& get_" << name << "_freq() {\n";
        out << "    static std::unordered_map<Metamer, float, MetamerHash> table = {\n";

        for (size_t i = 0; i < n; i++) {
            double freq = (double)sorted[i].second / total;
            out << "        {{" << sorted[i].first.aa_code << "ULL, "
                << sorted[i].first.codon_code << "ULL}, " << freq << "f},\n";
        }

        out << "    };\n";
        out << "    return table;\n";
        out << "}\n\n";
    };

    write_top_metamers("coding", coding, total_coding);
    write_top_metamers("boundary_5p", boundary_5p, total_5p);
    write_top_metamers("boundary_3p", boundary_3p, total_3p);

    out << "} // namespace agp\n";
}

// Main

int main(int argc, char* argv[]) {
    std::string cds_dir;
    std::string output_prefix = "metamers";
    int num_threads = omp_get_max_threads();

    // Parse arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--cds" && i + 1 < argc) {
            cds_dir = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            std::cerr << "Usage: " << argv[0] << " --cds <protein_fna_dir> [--output <prefix>] [--threads <n>]\n";
            return 0;
        }
    }

    if (cds_dir.empty()) {
        std::cerr << "Error: --cds directory required\n";
        std::cerr << "Usage: " << argv[0] << " --cds <protein_fna_dir> [--output <prefix>] [--threads <n>]\n";
        return 1;
    }

    init_codon_variants();

    std::cerr << "Metamer extractor for strand-discriminative training\n";
    std::cerr << "CDS directory: " << cds_dir << "\n";
    std::cerr << "Output prefix: " << output_prefix << "\n";
    std::cerr << "Threads: " << num_threads << "\n\n";

    // Find all CDS files
    std::vector<std::string> files;
    find_files(cds_dir, files, "_protein.fna.gz");
    if (files.empty()) {
        find_files(cds_dir, files, ".fna.gz");
    }
    if (files.empty()) {
        find_files(cds_dir, files, ".fa.gz");
    }

    std::cerr << "Found " << files.size() << " files\n";

    if (files.empty()) {
        std::cerr << "Error: No files found in " << cds_dir << "\n";
        return 1;
    }

    // Process files in parallel
    MetamerCounts global;
    std::atomic<size_t> files_done{0};

    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        ThreadBatch batch;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < files.size(); i++) {
            process_cds_file(files[i], batch, global);

            size_t done = ++files_done;
            if (done % 1000 == 0 || done == files.size()) {
                #pragma omp critical
                std::cerr << "\rProcessed " << done << "/" << files.size()
                          << " files, " << global.genes_processed << " genes" << std::flush;
            }
        }

        batch.flush(global);
    }

    std::cerr << "\n\nResults\n";
    std::cerr << "Genes processed: " << global.genes_processed << "\n";
    std::cerr << "Coding metamers: " << global.total_coding << " (" << global.coding.size() << " unique)\n";
    std::cerr << "5' boundary metamers: " << global.total_5p << " (" << global.boundary_5p.size() << " unique)\n";
    std::cerr << "3' boundary metamers: " << global.total_3p << " (" << global.boundary_3p.size() << " unique)\n";

    // Write outputs
    std::cerr << "\nWriting output files...\n";

    write_metamer_table(output_prefix + "_coding.tsv", global.coding, global.total_coding, "Coding");
    write_metamer_table(output_prefix + "_boundary_5p.tsv", global.boundary_5p, global.total_5p, "5' Boundary");
    write_metamer_table(output_prefix + "_boundary_3p.tsv", global.boundary_3p, global.total_3p, "3' Boundary");

    write_cpp_header(output_prefix + "_tables.hpp",
                     global.coding, global.boundary_5p, global.boundary_3p,
                     global.total_coding, global.total_5p, global.total_3p);

    std::cerr << "Done! Output files:\n";
    std::cerr << "  " << output_prefix << "_coding.tsv\n";
    std::cerr << "  " << output_prefix << "_boundary_5p.tsv\n";
    std::cerr << "  " << output_prefix << "_boundary_3p.tsv\n";
    std::cerr << "  " << output_prefix << "_tables.hpp\n";

    return 0;
}
