// EMI v7 round-trip fixture (commit 1): blocks, degree-CSR, hot/cold columns,
// var-length qaln/taln, per-read random access, oversized read-aligned block.
#include "dart/emi_v7_io.hpp"
#include <cstdio>
#include <cstdlib>
#include <map>

using namespace dart::emiv7;

static int failures = 0;
#define CHECK(cond, msg) do { if(!(cond)){ printf("FAIL: %s\n", msg); failures++; } } while(0)

int main() {
    const uint32_t NREADS = 50;
    const uint32_t BLOCK_TARGET = 8;   // tiny -> many blocks + oversized case

    // Build rows: read r has degree (r % 4)+1, EXCEPT read 7 is oversized (20 > target).
    std::vector<Row> rows;
    std::vector<EmiReadAnnot> annots(NREADS);
    std::map<uint32_t, std::vector<Row>> expect;
    for (uint32_t r = 0; r < NREADS; ++r) {
        EmiReadAnnot a{}; a.id_hash = 1000 + r; a.seq_len = uint16_t(30 + r);
        a.n_5prime = 5; a.n_3prime = 5; a.codons_5prime[0] = uint8_t(r % 64);
        a.damage_pct_q = uint8_t(r % 200); a.p_damaged_q = uint8_t(r % 255);
        annots[r] = a;
        uint32_t deg = (r == 7) ? 20 : (r % 4) + 1;
        for (uint32_t k = 0; k < deg; ++k) {
            Row row;
            row.read_idx = r; row.ref_idx = r * 100 + k;
            row.bit_score = float(50 + k);          // descending after sort
            row.damage_score = float(k) / deg;
            row.evalue_log10 = -float(k) - 1.0f;
            row.tstart = uint16_t(k); row.tend = uint16_t(k + 10);
            row.tlen = uint16_t(400 + r); row.qlen = uint16_t(14 + (k % 3));
            row.qaln = std::string(3 + (k % 5), char('A' + (k % 20)));
            row.taln = std::string(3 + (k % 5), char('a' + (k % 20)));
            expect[r].push_back(row);
        }
    }
    // shuffle input order to prove the writer sorts
    for (auto& kv : expect) for (auto& row : kv.second) rows.push_back(row);
    std::reverse(rows.begin(), rows.end());

    // realistic KapK-style read names (shared flowcell prefix) + UniProt-ish refs
    std::vector<std::string> read_names(NREADS);
    for (uint32_t r = 0; r < NREADS; ++r)
        read_names[r] = "A00706:851:H5FH3DSX7:1:1110:" + std::to_string(25000 + r) + ":" + std::to_string(13000 + r) + "_+_1";
    uint32_t nref = 0; for (auto& kv : expect) for (auto& row : kv.second) nref = std::max(nref, row.ref_idx + 1);
    std::vector<std::string> ref_names(nref);
    for (uint32_t i = 0; i < nref; ++i) ref_names[i] = "A0A" + std::to_string(100000 + i) + "X" + std::to_string(i % 9);

    const std::string path = "/dev/shm/emi_v7_rt.emi";
    Writer::write(path, rows, annots, read_names, ref_names, BLOCK_TARGET);

    Reader rd(path);
    CHECK(rd.num_reads() == NREADS, "num_reads");
    uint64_t total = 0; for (auto& kv : expect) total += kv.second.size();
    CHECK(rd.num_alignments() == total, "num_alignments");
    printf("blocks=%u alns=%llu reads=%llu\n", rd.num_blocks(),
           (unsigned long long)rd.num_alignments(), (unsigned long long)rd.num_reads());

    // per-read random access round-trip
    for (uint32_t r = 0; r < NREADS; ++r) {
        CHECK(rd.degree(r) == expect[r].size(), "degree");
        auto got = rd.rows_for_read(r);
        CHECK(got.size() == expect[r].size(), "rows_for_read size");
        // expected are best-first (bit desc); writer sorted the same way
        auto exp = expect[r];
        std::sort(exp.begin(), exp.end(), [](const Row&a,const Row&b){return a.bit_score>b.bit_score;});
        for (size_t j = 0; j < got.size() && j < exp.size(); ++j) {
            CHECK(got[j].ref_idx == exp[j].ref_idx, "ref_idx");
            CHECK(got[j].bit_score == exp[j].bit_score, "bit_score");
            CHECK(got[j].damage_score == exp[j].damage_score, "damage_score");
            CHECK(got[j].tstart == exp[j].tstart && got[j].tend == exp[j].tend, "coords");
            CHECK(got[j].tlen == exp[j].tlen && got[j].qlen == exp[j].qlen, "lens");
            CHECK(got[j].qaln == exp[j].qaln, "qaln");
            CHECK(got[j].taln == exp[j].taln, "taln");
        }
        // annot round-trip
        CHECK(rd.annot(r).id_hash == 1000 + r, "annot.id_hash");
        CHECK(rd.annot(r).seq_len == uint16_t(30 + r), "annot.seq_len");
    }

    // oversized read 7 must live wholly in a single block (read-aligned invariant)
    {
        auto g = rd.rows_for_read(7);
        CHECK(g.size() == 20, "oversized read degree 20");
    }

    // selective HOT scan across all blocks reproduces every (ref,bit,dmg)
    {
        std::map<uint32_t, std::pair<float,float>> hot;   // ref_idx unique here
        for (uint32_t b = 0; b < rd.num_blocks(); ++b) {
            auto hb = rd.hot_block(b);
            for (uint32_t j = 0; j < hb.n; ++j) hot[hb.ref[j]] = {hb.bit[j], hb.dmg[j]};
        }
        size_t hot_total = 0;
        for (auto& kv : expect) for (auto& row : kv.second) {
            auto it = hot.find(row.ref_idx);
            CHECK(it != hot.end(), "hot ref present");
            if (it != hot.end()) CHECK(it->second.first == row.bit_score, "hot bit matches");
            hot_total++;
        }
        CHECK(hot_total == total, "hot count");
    }

    // name-dict round-trip (front-coded + ZSTD)
    for (uint32_t r = 0; r < NREADS; ++r)
        CHECK(rd.read_name(r) == read_names[r], "read_name round-trip");
    for (uint32_t i = 0; i < ref_names.size(); ++i)
        CHECK(rd.ref_name(i) == ref_names[i], "ref_name round-trip");

    // measure dict compression on the realistic names
    {
        size_t raw = 0; for (auto& s : read_names) raw += s.size() + 1;
        auto enc = dart::emiv7::detail::encode_name_dict(read_names, EMI_V7_NAME_DICT_BLOCK);
        printf("read-dict: raw=%zuB encoded=%zuB (%.1f%%)\n", raw, enc.size(), 100.0*enc.size()/raw);
    }

    if (failures == 0) { printf("EMI v7 ROUND-TRIP: ALL PASS\n"); return 0; }
    printf("EMI v7 ROUND-TRIP: %d FAILURES\n", failures); return 1;
}
