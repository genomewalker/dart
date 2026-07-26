// Regression: a DART ORF id must identify exactly one predicted ORF.
//
// The bug this guards: src/main.cpp named every ORF "<read>_<strand>_<frame>"
// and discarded gene_num, so two ORFs predicted in the same frame of the same
// read shared an id. Measured on production, four libraries: 0.48-0.77% of ORF
// records carried a name already used (894,474 of 137,771,742 on LV7008956626).
// Downstream, damage-annotate interns query ids by name, so the two ORFs'
// alignments were fused into one dense id and became inseparable.
//
// Two halves:
//   1. the id grammar parses, for legacy 3-field ids and current 4-field ids,
//      including a two-digit ordinal (the case the old backward digit scan
//      silently mis-parsed);
//   2. `dart predict` on a read carrying two same-frame ORFs emits two distinct
//      FASTA ids. This half FAILS against the pre-fix code.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "dart/damage_index.hpp"

namespace {

int failures = 0;

bool expect(bool cond, const std::string& what) {
    if (!cond) { std::cerr << "  FAIL: " << what << "\n"; ++failures; }
    return cond;
}

void check_id_grammar() {
    std::cout << "Running ORF id grammar checks...\n";

    struct Case {
        const char* id;
        bool        ok;
        const char* read;
        char        strand;
        int         frame;
        int         ordinal;
    };
    // ordinal -1 means "absent" (legacy id).
    const Case cases[] = {
        {"readA_+_0",              true,  "readA",              '+', 0, -1},
        {"readA_-_2",              true,  "readA",              '-', 2, -1},
        {"readA_+_0_0",            true,  "readA",              '+', 0,  0},
        {"readA_-_1_7",            true,  "readA",              '-', 1,  7},
        // Two-digit ordinal. The pre-fix backward digit scan bailed here and
        // silently returned the hardcoded defaults {'+', 0}.
        {"readA_+_2_10",           true,  "readA",              '+', 2, 10},
        {"readA_-_0_123",          true,  "readA",              '-', 0, 123},
        // Illumina-style name: colons, dots and digits before the marker.
        {"A00706:851:H5FH3DSX7:1:1110:25000:13000_-_1_4", true,
         "A00706:851:H5FH3DSX7:1:1110:25000:13000", '-', 1, 4},
        // Not ORF ids: must be rejected, and strip must leave them alone.
        {"plain_read_name",        false, "plain_read_name",    '+', 0, -1},
        {"readA_+_",               false, "readA_+_",           '+', 0, -1},
        {"readA_+_0_",             false, "readA_+_0_",         '+', 0, -1},
        {"readA_+_0_x",            false, "readA_+_0_x",        '+', 0, -1},
        {"readA_*_0",              false, "readA_*_0",          '+', 0, -1},
        {"_+_0",                   false, "_+_0",               '+', 0, -1},
        {"",                       false, "",                   '+', 0, -1},
    };

    for (const Case& c : cases) {
        dart::DartOrfId p = dart::parse_dart_orf_id(c.id);
        std::string tag = std::string("'") + c.id + "'";
        if (!expect(p.ok == c.ok, tag + " ok=" + (p.ok ? "true" : "false") +
                                  ", expected " + (c.ok ? "true" : "false"))) continue;
        expect(p.read == c.read, tag + " read='" + std::string(p.read) + "', expected '" + c.read + "'");
        expect(p.strand == c.strand, tag + " strand");
        expect(p.frame == c.frame, tag + " frame=" + std::to_string(p.frame) +
                                   ", expected " + std::to_string(c.frame));
        expect(p.ordinal == c.ordinal, tag + " ordinal=" + std::to_string(p.ordinal) +
                                       ", expected " + std::to_string(c.ordinal));
        // strip_dart_suffix must agree with the parse on every input, including
        // the rejects: the .agd read key and the parsed read name are the same
        // thing, and a divergence between them is the silent-miss failure mode.
        expect(dart::strip_dart_suffix(c.id) == c.read,
               tag + " strip='" + std::string(dart::strip_dart_suffix(c.id)) + "'");
    }

    // Every ordinal of one read+frame must strip to the same read key, or the
    // damage index splits one read across N entries.
    for (int n = 0; n < 20; ++n) {
        std::string id = "readZ_+_1_" + std::to_string(n);
        expect(dart::strip_dart_suffix(id) == "readZ", "strip('" + id + "')");
    }
}

// Four 150 nt reads, each of which `dart predict` splits into two ORFs within
// one (strand, frame). They were not hand-designed: they are the first hits from
// a scan of 200,000 uniform-random reads, 34,658 of which (17.3%) produced a
// colliding id under the pre-fix naming. The predictor's default adaptive mode
// (frame_selector.cpp:1106) emits best-per-strand plus a margin-gated second
// fragment, and that second fragment may sit in the frame already emitted --
// e.g. read 0 yields "_+_0" at 1..69 and again at 73..150, two stop-separated
// fragments of forward frame 0.
const char* const kTwoOrfReads[] = {
    "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGT"
    "AAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAAC",
    "AGTTGAAATTCTGTTTTGTTTGTTCTATGCGGTCCTGTCGGAACAATACGTTTCTGTTTTACCTAACCCATTGCG"
    "TCAGTATTATACGCTGGCGGCGAATTTCGGCGCCCAGCAGCACTTCTATATCAGCCTCCTATACAATGACTACAC",
    "AAATTGAAGAGCTGCCACCGAATCTCTTGCTACCCGAGTGACTACGGCATTAGAGGCTCTTAGACTCGTTTACGG"
    "GCCACGGGTACGGAACGTACATAACACTCACGCTATATACCCAGTTAGAACTCCTGGGAACAGTCCAAATAGCCC",
    "GCCATACCTGTTAAGCACTGGGACTGAGCTGGTTGACACCGCCGCATAGTTTTAGTTTGCGCGGGGCAAATACGA"
    "GCCACTGTGACGTTCTGCATATTAGCCCACCGTATTCAAAAGAGCCCCTTTGGACACGCGTCGCTAAAATCTGAT",
};

struct GffGene { std::string seq; char strand; int frame; };

void check_predict_emits_distinct_ids(const char* dart_bin) {
    std::cout << "Running predict ORF-id uniqueness regression...\n";

    char tmpl[] = "/projects/caeg/scratch/user/tmp/dart_orfid_XXXXXX";
    const char* dir = mkdtemp(tmpl);
    if (!expect(dir != nullptr, "mkdtemp")) return;

    std::string reads = std::string(dir) + "/reads.fq";
    std::string gff   = std::string(dir) + "/pred.gff";
    std::string faa   = std::string(dir) + "/pred.faa";
    {
        // Replicated so the self-profile pass has a population to work with; the
        // ORF call is per read, so replication does not manufacture the collision.
        std::ofstream fq(reads);
        for (int rep = 0; rep < 50; ++rep)
            for (int k = 0; k < 4; ++k) {
                std::string seq = kTwoOrfReads[k];
                fq << "@synthread" << k << "_" << rep << "\n" << seq << "\n+\n"
                   << std::string(seq.size(), 'I') << "\n";
            }
    }

    std::string cmd = std::string(dart_bin) + " predict -i " + reads +
                      " -o " + gff + " --fasta-aa-masked " + faa +
                      " --orf-min-aa 10 --min-length 30 --allow-self-profile --threads 2 > " +
                      std::string(dir) + "/log 2>&1";
    int rc = std::system(cmd.c_str());
    if (!expect(rc == 0, "dart predict exited " + std::to_string(rc) + " (cmd: " + cmd + ")")) return;

    std::vector<GffGene> genes;
    {
        std::ifstream in(gff);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::vector<std::string> f;
            std::stringstream ss(line);
            std::string t;
            while (std::getline(ss, t, '\t')) f.push_back(t);
            if (f.size() < 8) continue;
            genes.push_back({f[0], f[6].empty() ? '?' : f[6][0], std::atoi(f[7].c_str())});
        }
    }

    std::vector<std::string> fasta_ids;
    {
        std::ifstream in(faa);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] != '>') continue;
            std::string h = line.substr(1);
            size_t sp = h.find_first_of(" \t");
            fasta_ids.push_back(sp == std::string::npos ? h : h.substr(0, sp));
        }
    }

    // The gate must be able to fire. If the fixture stops producing two ORFs in
    // one frame of one read, this test would pass vacuously -- so fail loudly
    // instead and make whoever changed the predictor supply a new fixture.
    std::map<std::string, int> per_frame;
    for (const GffGene& g : genes)
        per_frame[g.seq + "\t" + std::string(1, g.strand) + "\t" + std::to_string(g.frame)]++;
    int max_per_frame = 0;
    for (const auto& kv : per_frame) max_per_frame = std::max(max_per_frame, kv.second);
    if (!expect(max_per_frame >= 2,
                "FIXTURE NO LONGER EXERCISES THE BUG: no (read,strand,frame) has >=2 predicted "
                "genes (max " + std::to_string(max_per_frame) + " over " +
                std::to_string(genes.size()) + " genes). This test cannot fire; replace the "
                "fixture, do not delete the test.")) return;

    expect(fasta_ids.size() == genes.size(),
           "FASTA records (" + std::to_string(fasta_ids.size()) + ") != GFF genes (" +
           std::to_string(genes.size()) + ")");

    std::set<std::string> uniq(fasta_ids.begin(), fasta_ids.end());
    if (!expect(uniq.size() == fasta_ids.size(),
                "ORF ids collide: " + std::to_string(fasta_ids.size()) + " records but only " +
                std::to_string(uniq.size()) + " distinct ids -- two ORFs share a name")) {
        std::map<std::string, int> seen;
        for (const std::string& id : fasta_ids) seen[id]++;
        int shown = 0;
        for (const auto& kv : seen)
            if (kv.second > 1 && shown++ < 5)
                std::cerr << "    duplicate id: " << kv.first << " x" << kv.second << "\n";
    }

    // Ids must remain parseable and must all reduce to the read they came from.
    for (const std::string& id : fasta_ids) {
        dart::DartOrfId p = dart::parse_dart_orf_id(id);
        if (!expect(p.ok, "emitted id '" + id + "' does not parse")) continue;
        expect(p.read.rfind("synthread", 0) == 0,
               "emitted id '" + id + "' strips to '" + std::string(p.read) + "'");
    }

    std::cout << "  genes=" << genes.size() << " fasta_ids=" << fasta_ids.size()
              << " distinct=" << uniq.size() << " max_genes_per_frame=" << max_per_frame << "\n";
}

} // namespace

int main(int argc, char** argv) {
    check_id_grammar();
    if (argc > 1) check_predict_emits_distinct_ids(argv[1]);
    else { std::cerr << "usage: " << argv[0] << " <path-to-dart>\n"; return 2; }

    if (failures) {
        std::cout << "ORF ID UNIQUENESS: " << failures << " FAILURE(S)\n";
        return 1;
    }
    std::cout << "ORF ID UNIQUENESS: ALL PASS\n";
    return 0;
}
