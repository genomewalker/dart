// emi-merge end-to-end golden: two EMI v7 inputs over overlapping reads and overlapping targets,
// merged by `dart emi-merge`, checked by `dart emi-verify`, and compared row by row against a
// hand-built expected union. The expected table follows the contract documented at the top of
// src/cli/cmd_emi_merge.cpp:
//   - reads are unioned BY NAME (across inputs and for duplicate names within one input);
//   - refs are unioned BY NAME, first-encounter order (input 0 dense order, then new names);
//   - duplicate (read,target) pairs keep the best bit_score; ties keep the gather-order row
//     (input 0 before input 1, lower dense id before higher);
//   - survivors are ordered per read by (bit desc, global ref asc);
//   - per-read annot comes from the first occurrence (input 0, lowest dense id);
//   - every other payload field passes through unchanged.
// Two negatives: inputs with different damage params are refused, and emi-verify rejects a file
// that carries a duplicate (read,target) pair.
//
// Usage: dart-emi-merge-golden <dart binary>

#include "dart/emi_v7_io.hpp"

#include <sys/wait.h>

#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

using namespace dart::emiv7;

namespace {

int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s\n", msg); failures++; } } while (0)

std::string make_tmpdir() {
    const char* base = std::getenv("TMPDIR");
    std::string tmpl = std::string(base && *base ? base : "/scratch/tmp") + "/dart_emi_merge_XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    if (!mkdtemp(buf.data())) { printf("mkdtemp failed\n"); std::exit(2); }
    return std::string(buf.data());
}

int run(const std::string& cmd, const std::string& log) {
    const int rc = std::system((cmd + " > " + log + " 2>&1").c_str());
    return WIFEXITED(rc) ? WEXITSTATUS(rc) : -1;
}

// One alignment with every payload field derived from `tag`, so a wrong survivor is caught on
// any column, not just bit_score. qaln carries a readable label.
Row mk(uint32_t read, uint32_t ref, float bit, const char* label, uint16_t tag) {
    Row r;
    r.read_idx = read; r.ref_idx = ref; r.bit_score = bit;
    r.damage_score = tag / 1000.0f;
    r.evalue_log10 = -1.0f - tag;
    r.identity_q = uint16_t(1000 + tag); r.aln_len = uint16_t(20 + tag);
    r.qstart = tag; r.qend = uint16_t(tag + 10);
    r.tstart = uint16_t(2 * tag); r.tend = uint16_t(2 * tag + 10);
    r.tlen = uint16_t(300 + tag); r.qlen = uint16_t(30 + tag);
    r.dmg_k = tag; r.dmg_m = uint16_t(tag + 1);
    r.dmg_ll_a = -0.5f * tag; r.dmg_ll_m = -0.25f * tag;
    r.qaln = label; r.taln = std::string(label) + "_t";
    return r;
}

EmiReadAnnot annot(uint64_t id_hash) {
    EmiReadAnnot a{}; a.id_hash = id_hash; a.seq_len = uint16_t(id_hash % 200); return a;
}

bool same_row(const Row& a, const Row& b) {
    return a.bit_score == b.bit_score && a.damage_score == b.damage_score &&
           a.evalue_log10 == b.evalue_log10 && a.identity_q == b.identity_q &&
           a.aln_len == b.aln_len && a.qstart == b.qstart && a.qend == b.qend &&
           a.tstart == b.tstart && a.tend == b.tend && a.tlen == b.tlen && a.qlen == b.qlen &&
           a.dmg_k == b.dmg_k && a.dmg_m == b.dmg_m && a.dmg_ll_a == b.dmg_ll_a &&
           a.dmg_ll_m == b.dmg_ll_m && a.qaln == b.qaln && a.taln == b.taln;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) { printf("usage: %s <dart binary>\n", argv[0]); return 2; }
    const std::string dart = argv[1];
    const std::string dir = make_tmpdir();
    const float d_max = 0.3f, lambda = 0.2f;

    // ---- input A: dense reads 0..4 are r0 r1 r2 r3 r3 (r3 duplicated), refs t0..t3 ----
    const std::vector<std::string> a_reads = {"r0", "r1", "r2", "r3", "r3"};
    const std::vector<std::string> a_refs  = {"t0", "t1", "t2", "t3"};
    std::vector<EmiReadAnnot> a_annots;
    for (uint64_t i = 0; i < a_reads.size(); ++i) a_annots.push_back(annot(100 + i));
    std::vector<Row> a_rows = {
        mk(0, 0, 100.f, "A_r0_t0", 1),    // ties with B r0/t0 at 100: input 0 wins
        mk(0, 1,  90.f, "A_r0_t1", 2),
        mk(1, 1,  80.f, "A_r1_t1", 3),
        mk(1, 2,  80.f, "A_r1_t2", 4),    // loses to B r1/t2 at 85
        mk(2, 2,  70.f, "A_r2_t2", 5),    // two HSPs to one target within A: 70 survives
        mk(2, 2,  60.f, "A_r2_t2_hsp2", 6),
        mk(3, 0,  50.f, "A_r3_t0", 7),    // r3 (dense 3) and r3 (dense 4) are one read
        mk(4, 3,  55.f, "A_r3_t3", 8),
        mk(4, 0,  45.f, "A_r3_t0_dup", 9),  // loses to dense-3 r3/t0 at 50
    };
    // ---- input B: reads r5 r1 r0 r6 (different order), refs t2 t4 t0 t5 (t2,t0 overlap A) ----
    const std::vector<std::string> b_reads = {"r5", "r1", "r0", "r6"};
    const std::vector<std::string> b_refs  = {"t2", "t4", "t0", "t5"};
    std::vector<EmiReadAnnot> b_annots;
    for (uint64_t i = 0; i < b_reads.size(); ++i) b_annots.push_back(annot(200 + i));
    std::vector<Row> b_rows = {
        mk(0, 1,  30.f, "B_r5_t4", 11),
        mk(1, 0,  85.f, "B_r1_t2", 12),   // beats A r1/t2
        mk(1, 1,  20.f, "B_r1_t4", 13),
        mk(2, 2, 100.f, "B_r0_t0", 14),   // tie, loses to A
        mk(2, 3,  95.f, "B_r0_t5", 15),
        mk(3, 2,  10.f, "B_r6_t0", 16),
    };

    const std::string a_path = dir + "/a.emi", b_path = dir + "/b.emi", m_path = dir + "/merged.emi";
    Writer::write(a_path, a_rows, a_annots, a_reads, a_refs, EMI_V7_BLOCK_TARGET_ROWS, d_max, lambda);
    Writer::write(b_path, b_rows, b_annots, b_reads, b_refs, EMI_V7_BLOCK_TARGET_ROWS, d_max, lambda);

    // Precondition: Writer keeps both HSPs of A r2/t2 (dedup is emi-merge's job, not the writer's).
    { Reader a(a_path); CHECK(a.rows_for_read(2).size() == 2, "writer keeps duplicate (read,target) rows"); }

    CHECK(run(dart + " emi-merge " + m_path + " 2 " + a_path + " " + b_path, dir + "/merge.log") == 0, "emi-merge exit 0");
    CHECK(run(dart + " emi-verify " + m_path + " 2", dir + "/verify.log") == 0, "emi-verify exit 0 on merged");

    // ---- hand-built expected union: read name -> rows in (bit desc, ref asc) order ----
    const std::vector<std::string> exp_reads = {"r0", "r1", "r2", "r3", "r5", "r6"};
    const std::vector<std::string> exp_refs  = {"t0", "t1", "t2", "t3", "t4", "t5"};
    std::map<std::string, std::vector<std::pair<std::string, Row>>> expect;
    expect["r0"] = {{"t0", a_rows[0]}, {"t5", b_rows[4]}, {"t1", a_rows[1]}};
    expect["r1"] = {{"t2", b_rows[1]}, {"t1", a_rows[2]}, {"t4", b_rows[2]}};
    expect["r2"] = {{"t2", a_rows[4]}};
    expect["r3"] = {{"t3", a_rows[7]}, {"t0", a_rows[6]}};
    expect["r5"] = {{"t4", b_rows[0]}};
    expect["r6"] = {{"t0", b_rows[5]}};
    const std::map<std::string, uint64_t> exp_annot = {
        {"r0", 100}, {"r1", 101}, {"r2", 102}, {"r3", 103}, {"r5", 200}, {"r6", 203}};

    Reader m(m_path);
    CHECK(m.num_reads() == exp_reads.size(), "merged num_reads == 6");
    CHECK(m.num_refs() == exp_refs.size(), "merged num_refs == 6");
    CHECK(m.num_alignments() == 11, "merged num_alignments == 11");
    CHECK(m.d_max() == d_max && m.lambda() == lambda, "damage params carried");
    for (uint32_t r = 0; r < exp_reads.size() && r < m.num_reads(); ++r)
        CHECK(m.read_name(r) == exp_reads[r], "read dict first-encounter order");
    for (uint32_t f = 0; f < exp_refs.size() && f < m.num_refs(); ++f)
        CHECK(m.ref_name(f) == exp_refs[f], "ref dict first-encounter order");

    uint64_t seen = 0;
    for (uint32_t r = 0; r < m.num_reads(); ++r) {
        const std::string name = m.read_name(r);
        auto it = expect.find(name);
        if (it == expect.end()) { CHECK(false, "unexpected read name"); continue; }
        const auto& exp = it->second;
        auto got = m.rows_for_read(r);
        CHECK(m.degree(r) == exp.size(), "degree == survivors");
        CHECK(got.size() == exp.size(), "row count == survivors");
        CHECK(m.annot(r).id_hash == exp_annot.at(name), "annot from first occurrence");
        for (size_t k = 0; k < got.size() && k < exp.size(); ++k) {
            CHECK(got[k].read_idx == r, "read_idx remapped");
            CHECK(m.ref_name(got[k].ref_idx) == exp[k].first, "survivor target and order");
            CHECK(same_row(got[k], exp[k].second), "survivor payload byte-for-byte");
        }
        seen += got.size();
    }
    CHECK(seen == 11, "scanned rows == 11");

    // ---- negative: damage params must match across inputs ----
    {
        const std::string c_path = dir + "/c.emi";
        Writer::write(c_path, b_rows, b_annots, b_reads, b_refs, EMI_V7_BLOCK_TARGET_ROWS, 0.5f, lambda);
        CHECK(run(dart + " emi-merge " + dir + "/bad.emi 2 " + a_path + " " + c_path, dir + "/merge_bad.log") == 1,
              "emi-merge refuses differing d_max");
    }
    // ---- negative: emi-verify rejects a duplicate (read,target) pair (input A has one) ----
    CHECK(run(dart + " emi-verify " + a_path + " 2", dir + "/verify_a.log") == 1, "emi-verify fails on dup pair");

    if (failures == 0) { printf("EMI-MERGE GOLDEN: ALL PASS (%s)\n", dir.c_str()); return 0; }
    printf("EMI-MERGE GOLDEN: %d FAILURES (%s)\n", failures, dir.c_str());
    return 1;
}
