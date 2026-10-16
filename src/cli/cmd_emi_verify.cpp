// emi-verify: structural + invariant check for a v7 .emi. Confirms the file is internally
// consistent (degree CSR == per-read row count, totals match the header, footer magic) and that
// the joint-EM invariant holds: within every read the alignment rows carry DISTINCT targets (no
// duplicate (read,target) pairs) and are in the canonical per-read order (bit_score desc, ref
// asc). Used to gate emi-merge output before the EM runs. Exits non-zero on any violation.
#include "subcommand.hpp"
#include "dart/emi_v7_io.hpp"
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <atomic>
#include <string>
#include <thread>
#include <vector>
#include <unordered_set>
#include <algorithm>

namespace dart {
namespace cli {

using namespace dart::emiv7;

int cmd_emi_verify(int argc, char* argv[]) {
    if (argc >= 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        printf("usage: emi-verify <in.emi> [threads=8]\n\n"
               "Check an EMI v7: degree CSR, totals vs header, footer magic, no duplicate\n"
               "(read, target) pairs, canonical per-read order. Exits non-zero on any violation.\n");
        return 0;
    }
    if (argc < 2) { printf("usage: emi-verify <in.emi> [threads]\n"); return 2; }
    int nthreads = argc >= 3 ? atoi(argv[2]) : 8; if (nthreads < 1) nthreads = 1;

    Reader in(argv[1]);
    const uint32_t nreads = (uint32_t)in.num_reads();
    const uint64_t nrefs = in.num_refs();
    const uint64_t nal = in.num_alignments();

    std::atomic<uint64_t> total_rows{0};
    std::atomic<uint64_t> dup_pairs{0};     // (read,target) collisions within a read
    std::atomic<uint64_t> order_viol{0};    // rows not in (bit desc, ref asc) order
    std::atomic<uint64_t> deg_mismatch{0};  // degree(read) != rows_for_read().size()
    std::atomic<uint64_t> ref_oob{0};       // ref_idx >= num_refs
    std::atomic<uint32_t> next{0};
    const uint32_t chunk = 4096;

    auto worker = [&] {
        std::unordered_set<uint32_t> seen;
        uint64_t t = 0, dp = 0, ov = 0, dm = 0, ob = 0;
        for (;;) {
            uint32_t s = next.fetch_add(chunk);
            if (s >= nreads) break;
            uint32_t e = std::min(nreads, s + chunk);
            for (uint32_t r = s; r < e; ++r) {
                std::vector<Row> rows = in.rows_for_read(r);
                if (rows.size() != in.degree(r)) dm++;
                t += rows.size();
                seen.clear();
                for (size_t k = 0; k < rows.size(); ++k) {
                    if (rows[k].ref_idx >= nrefs) ob++;
                    if (!seen.insert(rows[k].ref_idx).second) dp++;
                    if (k > 0) {
                        const Row& a = rows[k - 1]; const Row& b = rows[k];
                        bool ok = (a.bit_score > b.bit_score) ||
                                  (a.bit_score == b.bit_score && a.ref_idx <= b.ref_idx);
                        if (!ok) ov++;
                    }
                }
            }
        }
        total_rows += t; dup_pairs += dp; order_viol += ov; deg_mismatch += dm; ref_oob += ob;
    };
    std::vector<std::thread> pool;
    for (int i = 0; i < nthreads; ++i) pool.emplace_back(worker);
    for (auto& th : pool) th.join();

    const uint64_t tr = total_rows.load();
    printf("emi-verify %s\n", argv[1]);
    printf("  reads=%u refs=%llu alignments(header)=%llu rows(scanned)=%llu\n",
           nreads, (unsigned long long)nrefs, (unsigned long long)nal, (unsigned long long)tr);
    printf("  dup (read,target) pairs = %llu\n", (unsigned long long)dup_pairs.load());
    printf("  per-read order violations = %llu\n", (unsigned long long)order_viol.load());
    printf("  degree/rowcount mismatch  = %llu\n", (unsigned long long)deg_mismatch.load());
    printf("  ref_idx out of range      = %llu\n", (unsigned long long)ref_oob.load());
    bool ok = (tr == nal) && dup_pairs.load() == 0 && order_viol.load() == 0 &&
              deg_mismatch.load() == 0 && ref_oob.load() == 0;
    if (tr != nal) printf("  FAIL rows(scanned)=%llu != header alignments=%llu\n",
                          (unsigned long long)tr, (unsigned long long)nal);
    printf("VERIFY: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

namespace {
    struct EmiVerifyRegistrar {
        EmiVerifyRegistrar() {
            SubcommandRegistry::instance().register_command(
                "emi-verify",
                "Check a v7 .emi: no duplicate (read,target) pairs, per-read order, count integrity",
                cmd_emi_verify, 54);
        }
    };
    static EmiVerifyRegistrar registrar;
}

}  // namespace cli
}  // namespace dart
