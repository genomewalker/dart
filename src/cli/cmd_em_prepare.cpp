// agp em-prepare: Convert TSV alignments to binary index for EM
//
// Two-phase workflow:
//   1. em-prepare: TSV → binary .emi (slow, once)
//   2. em-solve: mmap .emi → run EM (fast, repeatable)

#include "subcommand.hpp"
#include "agp/version.h"
#include "agp/em_index.hpp"
#include "agp/fast_tsv.hpp"
#include "agp/damage_index_reader.hpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>
#include <string>
#include <chrono>

namespace agp {
namespace cli {

int cmd_em_prepare(int argc, char* argv[]) {
    std::string tsv_file;
    std::string output_file;
    std::string agd_file;
    float d_max = 0.0f;
    float lambda = 0.3f;
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "--tsv") == 0 || strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
            tsv_file = argv[++i];
        } else if ((strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--damage-index") == 0 && i + 1 < argc) {
            agd_file = argv[++i];
        } else if (strcmp(argv[i], "--d-max") == 0 && i + 1 < argc) {
            d_max = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda") == 0 && i + 1 < argc) {
            lambda = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            std::cerr << "Usage: agp em-prepare --tsv <hits.tsv> -o <output.emi> [options]\n\n";
            std::cerr << "Convert MMseqs2 TSV alignments to binary index for EM.\n\n";
            std::cerr << "Required:\n";
            std::cerr << "  --tsv, -i <file>       Input MMseqs2 hits (16-column format)\n";
            std::cerr << "  --output, -o <file>    Output binary index (.emi)\n\n";
            std::cerr << "Options:\n";
            std::cerr << "  --damage-index <file>  AGD index for damage scores (from agp predict)\n";
            std::cerr << "  --d-max <float>        Sample damage level (0-1), stored in index\n";
            std::cerr << "  --lambda <float>       Damage decay parameter (default: 0.3)\n";
            std::cerr << "  -v, --verbose          Verbose output\n";
            std::cerr << "  --help, -h             Show this help\n";
            return 0;
        }
    }

    if (tsv_file.empty() || output_file.empty()) {
        std::cerr << "Error: --tsv and --output are required\n";
        std::cerr << "Use --help for usage.\n";
        return 1;
    }

    auto t_start = std::chrono::steady_clock::now();

    // Load AGD index if provided
    std::unique_ptr<DamageIndexReader> agd;
    if (!agd_file.empty()) {
        agd = std::make_unique<DamageIndexReader>(agd_file);
        if (!agd->is_valid()) {
            std::cerr << "Error: Cannot load damage index: " << agd_file << "\n";
            return 1;
        }
        if (verbose) {
            std::cerr << "Loaded AGD index: " << agd->record_count() << " entries\n";
        }
    }

    // Memory-map TSV file
    MappedTSV tsv(tsv_file.c_str());
    if (!tsv.valid()) {
        std::cerr << "Error: Cannot mmap TSV file: " << tsv_file << "\n";
        return 1;
    }
    if (verbose) {
        std::cerr << "Mapped TSV: " << (tsv.size() / (1024*1024)) << " MB\n";
    }

    // Build name -> index mappings
    std::unordered_map<std::string, uint32_t> read_to_idx;
    std::unordered_map<std::string, uint32_t> ref_to_idx;
    std::unordered_map<uint32_t, uint32_t> ref_lengths;

    // First pass: collect unique names and max ref lengths
    struct TempAlignment {
        uint32_t read_idx;
        uint32_t ref_idx;
        float bit_score;
        float damage_score;
        float evalue_log10;
        uint16_t identity_q;
        uint16_t aln_len;
        uint16_t tstart;
        uint16_t tend;
        uint16_t tlen;
        uint16_t qlen;
    };
    std::vector<TempAlignment> temp_alns;
    temp_alns.reserve(tsv.size() / 100);

    size_t n_rows = parallel_parse_tsv(tsv, [&](const TSVLine& line) {
        if (line.n_fields < 12) return;

        // Parse fields: query, target, fident, alnlen, mismatch, gapopen,
        //               qstart, qend, tstart, tend, evalue, bits, [qlen, tlen]
        std::string query(line.fields[0]);
        std::string target(line.fields[1]);
        float fident = fast_stof(line.fields[2]);
        uint32_t alnlen = fast_stou(line.fields[3]);
        uint32_t tstart = fast_stou(line.fields[8]);
        uint32_t tend = fast_stou(line.fields[9]);
        float evalue = fast_stof(line.fields[10]);
        float bits = fast_stof(line.fields[11]);
        uint32_t qlen = (line.n_fields >= 13) ? fast_stou(line.fields[12]) : 0;
        uint32_t tlen = (line.n_fields >= 14) ? fast_stou(line.fields[13]) : 0;

        // Thread-safe index assignment
        static std::mutex mu;
        std::lock_guard<std::mutex> lock(mu);

        // Get or create read index
        auto [it_r, inserted_r] = read_to_idx.emplace(query, static_cast<uint32_t>(read_to_idx.size()));
        uint32_t read_idx = it_r->second;

        // Get or create ref index
        auto [it_t, inserted_t] = ref_to_idx.emplace(target, static_cast<uint32_t>(ref_to_idx.size()));
        uint32_t ref_idx = it_t->second;

        // Track ref length
        if (tlen > 0) {
            auto& stored = ref_lengths[ref_idx];
            if (tlen > stored) stored = tlen;
        }

        // Lookup damage score from AGD
        float damage_score = 0.0f;
        if (agd) {
            auto rec = agd->find(query);
            if (rec) {
                // p_damaged_q is quantized to 0-255
                damage_score = static_cast<float>(rec->p_damaged_q) / 255.0f;
            }
        }

        TempAlignment ta;
        ta.read_idx = read_idx;
        ta.ref_idx = ref_idx;
        ta.bit_score = bits;
        ta.damage_score = damage_score;
        ta.evalue_log10 = (evalue > 0) ? std::log10(evalue) : -300.0f;
        ta.identity_q = static_cast<uint16_t>(std::min(fident, 1.0f) * 65535.0f);
        ta.aln_len = static_cast<uint16_t>(std::min(alnlen, 65535u));
        ta.tstart = static_cast<uint16_t>(std::min(tstart, 65535u));
        ta.tend = static_cast<uint16_t>(std::min(tend, 65535u));
        ta.tlen = static_cast<uint16_t>(std::min(tlen, 65535u));
        ta.qlen = static_cast<uint16_t>(std::min(qlen, 65535u));
        temp_alns.push_back(ta);
    }, true, 12);  // skip header, min 12 fields

    if (verbose) {
        auto t_parse = std::chrono::steady_clock::now();
        auto parse_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_parse - t_start).count();
        std::cerr << "Parsed TSV: " << n_rows << " rows, "
                  << read_to_idx.size() << " reads, "
                  << ref_to_idx.size() << " refs in "
                  << parse_ms << " ms\n";
    }

    // Sort by (read_idx, bit_score DESC) for CSR layout
    // CRITICAL: Within each read, best hit (highest bit_score) must come first
    // so that damage-annotate's "first alignment per read" is the best hit
    std::sort(temp_alns.begin(), temp_alns.end(),
        [](const TempAlignment& a, const TempAlignment& b) {
            if (a.read_idx != b.read_idx) return a.read_idx < b.read_idx;
            return a.bit_score > b.bit_score;  // Descending by score
        });

    // Build EMI index
    EMIndexWriter writer(output_file);

    // Add references
    std::vector<std::string> ref_names(ref_to_idx.size());
    for (const auto& [name, idx] : ref_to_idx) {
        ref_names[idx] = name;
    }
    for (const auto& name : ref_names) {
        writer.add_ref(name);
    }

    // Group alignments by read and add to index
    std::vector<std::string> read_names(read_to_idx.size());
    for (const auto& [name, idx] : read_to_idx) {
        read_names[idx] = name;
    }

    uint32_t current_read = UINT32_MAX;
    std::vector<EMIAlignment> read_alns;
    read_alns.reserve(16);

    for (const auto& ta : temp_alns) {
        if (ta.read_idx != current_read) {
            // Flush previous read
            if (current_read != UINT32_MAX && !read_alns.empty()) {
                writer.add_read(read_names[current_read], read_alns);
            }
            current_read = ta.read_idx;
            read_alns.clear();
        }

        EMIAlignment ea;
        ea.ref_idx = ta.ref_idx;
        ea.bit_score = ta.bit_score;
        ea.damage_score = ta.damage_score;
        ea.evalue_log10 = ta.evalue_log10;
        ea.identity_q = ta.identity_q;
        ea.aln_len = ta.aln_len;
        ea.tstart = ta.tstart;
        ea.tend = ta.tend;
        ea.tlen = ta.tlen;
        ea.qlen = ta.qlen;
        read_alns.push_back(ea);
    }
    // Flush last read
    if (current_read != UINT32_MAX && !read_alns.empty()) {
        writer.add_read(read_names[current_read], read_alns);
    }

    // Finalize and write
    writer.finalize(d_max, lambda);

    auto t_end = std::chrono::steady_clock::now();
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    if (verbose) {
        std::cerr << "Written: " << output_file << "\n";
        std::cerr << "Total time: " << total_ms << " ms\n";
        std::cerr << "Throughput: " << (n_rows * 1000 / (total_ms + 1)) << " rows/sec\n";
    }

    return 0;
}

namespace {
    struct EMPrepareRegistrar {
        EMPrepareRegistrar() {
            SubcommandRegistry::instance().register_command(
                "em-prepare",
                "Convert TSV alignments to binary index for fast EM",
                cmd_em_prepare, 50);
        }
    };
    static EMPrepareRegistrar registrar;
}

}  // namespace cli
}  // namespace agp
