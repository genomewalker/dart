// emi-merge memory gate: peak RSS of `dart emi-merge` on a synthetic two-input union of a few
// million rows must stay under a fixed bound. emi-merge is streaming (pass2 holds at most
// 3*threads+4 encoded blocks in its reorder ring), so its peak must not grow with the row count;
// this gate catches a regression that makes rows resident (e.g. gathering every block before
// writing). The thread count is pinned because the ring, and so the peak, scales with it.
//
// Peak RSS is the kernel's hiwater_rss counter: sampled as VmHWM from /proc/<pid>/status while
// the child runs, and taken authoritatively from wait4()'s ru_maxrss when it exits (the /proc
// entry is gone by then; both report the same counter in kB). The gate uses the larger of the two.
//
// The merged row count is also checked analytically and the output passes `dart emi-verify`,
// so the multi-block pass2 path (reorder ring, post-pass2 degree section, degree escape for a
// read with >= 65535 targets) is exercised on real sizes, which the golden fixture is too small
// to reach (one merged block).
//
// Usage: dart-emi-merge-memory-gate <dart binary>

#include "dart/emi_v7_io.hpp"

#include <sys/resource.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using namespace dart::emiv7;

namespace {

// Synthetic layout: input A has NA reads, input B has NB reads of which OVERLAP share A's names.
// Read i of A hits targets (i*7 + k) % NREF for k in [0,DEG); the same read in B hits
// (i*7 + k + DEG/2) % NREF, so DEG/2 (read,target) pairs are duplicated across inputs and the
// merged row count is NA*DEG + NB*DEG - OVERLAP*DEG/2 (+ the escape read's extra targets).
constexpr uint32_t NA = 500000, NB = 500000, OVERLAP = 250000, DEG = 4, NREF = 100000;
constexpr uint32_t ESCAPE_DEG = 70000;      // one A read with >= EMI_V7_DEGREE_ESCAPE targets
constexpr int THREADS = 4;
// Measured in a Release build with THREADS=4: peak_rss_kb=240704 (235 MB), the
// same from the /proc VmHWM samples and from ru_maxrss. Bound = 1.5x that, rounded to 360 MB.
constexpr long BOUND_KB = 360L * 1024L;

int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s\n", msg); failures++; } } while (0)

std::string make_tmpdir() {
    const char* base = std::getenv("TMPDIR");
    std::string tmpl = std::string(base && *base ? base : "/tmp") + "/dart_emi_merge_mem_XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    if (!mkdtemp(buf.data())) { printf("mkdtemp failed\n"); std::exit(2); }
    return std::string(buf.data());
}

std::string read_name(uint32_t i) { return "A00706:851:H5FH3DSX7:1:1110:" + std::to_string(i) + "_+_1"; }

Row mk(uint32_t read, uint32_t ref, float bit) {
    Row r; r.read_idx = read; r.ref_idx = ref; r.bit_score = bit; r.damage_score = 0.1f;
    r.evalue_log10 = -5.0f; r.tlen = 300; r.qlen = 30; r.aln_len = 30;
    r.qaln = "ACDEFGHI"; r.taln = "ACDEFGHI";
    return r;
}

void write_input(const std::string& path, uint32_t nreads, uint32_t first_name, uint32_t koff, bool escape) {
    std::vector<Row> rows; rows.reserve((size_t)nreads * DEG + (escape ? ESCAPE_DEG : 0));
    std::vector<std::string> names(nreads);
    for (uint32_t r = 0; r < nreads; ++r) {
        const uint32_t i = first_name + r;
        names[r] = read_name(i);
        for (uint32_t k = 0; k < DEG; ++k) rows.push_back(mk(r, (uint32_t)(((uint64_t)i * 7 + k + koff) % NREF), 100.f - k));
    }
    if (escape) for (uint32_t k = 0; k < ESCAPE_DEG; ++k) rows.push_back(mk(nreads - 1, k, 1.0f + k * 1e-3f));
    std::vector<std::string> refs(NREF);
    for (uint32_t f = 0; f < NREF; ++f) refs[f] = "A0A" + std::to_string(100000 + f);
    std::vector<EmiReadAnnot> annots(nreads);
    Writer::write(path, std::move(rows), annots, names, refs, EMI_V7_BLOCK_TARGET_ROWS, 0.3f, 0.2f);
}

long vmhwm_kb(pid_t pid) {
    std::ifstream st("/proc/" + std::to_string(pid) + "/status");
    std::string line;
    while (std::getline(st, line))
        if (line.rfind("VmHWM:", 0) == 0) return std::atol(line.c_str() + 6);
    return -1;
}

// fork/exec argv, sample VmHWM while it runs, return exit status; peak_kb = max(sampled, ru_maxrss).
int run_measured(const std::vector<std::string>& argv, const std::string& log, long& peak_kb) {
    fflush(stdout);   // a buffered line must not be inherited (and re-flushed) by the child
    pid_t pid = fork();
    if (pid == 0) {
        if (!freopen(log.c_str(), "w", stdout)) _exit(126);
        dup2(fileno(stdout), fileno(stderr));
        std::vector<char*> av; for (const auto& s : argv) av.push_back(const_cast<char*>(s.c_str())); av.push_back(nullptr);
        execv(av[0], av.data()); _exit(127);
    }
    long sampled = 0; int status = 0; struct rusage ru{};
    const timespec tick{0, 2000000};   // 2 ms
    for (;;) {
        pid_t w = wait4(pid, &status, WNOHANG, &ru);
        if (w == pid) break;
        if (w < 0) { peak_kb = -1; return -1; }
        long v = vmhwm_kb(pid); if (v > sampled) sampled = v;
        nanosleep(&tick, nullptr);
    }
    peak_kb = ru.ru_maxrss > sampled ? ru.ru_maxrss : sampled;
    printf("vmhwm_sampled_kb=%ld ru_maxrss_kb=%ld\n", sampled, (long)ru.ru_maxrss);
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) { printf("usage: %s <dart binary>\n", argv[0]); return 2; }
    const std::string dart = argv[1];
    const std::string dir = make_tmpdir();
    const std::string a_path = dir + "/a.emi", b_path = dir + "/b.emi", m_path = dir + "/merged.emi";

    // Generate in a child so the writer's resident rows never count toward anything measured here.
    {
        fflush(stdout);
        pid_t pid = fork();
        if (pid == 0) {
            write_input(a_path, NA, 0, 0, true);
            write_input(b_path, NB, NA - OVERLAP, DEG / 2, false);
            _exit(0);
        }
        int st = 0; waitpid(pid, &st, 0);
        if (!WIFEXITED(st) || WEXITSTATUS(st) != 0) { printf("FAIL: input generation\n"); return 1; }
    }
    const uint64_t expected_rows = (uint64_t)NA * DEG + (uint64_t)NB * DEG - (uint64_t)OVERLAP * (DEG / 2) + ESCAPE_DEG;

    long peak_kb = 0;
    const int rc = run_measured({dart, "emi-merge", m_path, std::to_string(THREADS), a_path, b_path}, dir + "/merge.log", peak_kb);
    CHECK(rc == 0, "emi-merge exit 0");
    printf("peak_rss_kb=%ld bound_kb=%ld input_rows=%llu merged_rows_expected=%llu\n", peak_kb, BOUND_KB,
           (unsigned long long)((uint64_t)NA * DEG + (uint64_t)NB * DEG + ESCAPE_DEG), (unsigned long long)expected_rows);
    CHECK(peak_kb > 0, "peak RSS measured");
    CHECK(peak_kb <= BOUND_KB, "emi-merge peak RSS within bound");

    {
        Reader m(m_path);
        CHECK(m.num_alignments() == expected_rows, "merged row count == analytic union");
        CHECK(m.num_reads() == NA + NB - OVERLAP, "merged read count");
        CHECK(m.num_refs() == NREF, "merged ref count");
        CHECK(m.num_blocks() > 1, "merge produced several blocks");
        // read NA-1 is in both inputs (DEG + DEG/2 distinct regular targets) and carries the escape targets
        CHECK(m.degree(NA - 1) == DEG + DEG / 2 + ESCAPE_DEG, "escape-degree read kept all distinct targets");
    }
    long verify_kb = 0;
    CHECK(run_measured({dart, "emi-verify", m_path, std::to_string(THREADS)}, dir + "/verify.log", verify_kb) == 0,
          "emi-verify exit 0 on stress merge");

    if (failures == 0) {
        std::remove(a_path.c_str()); std::remove(b_path.c_str()); std::remove(m_path.c_str());
        printf("EMI-MERGE MEMORY GATE: ALL PASS (%s)\n", dir.c_str()); return 0;
    }
    printf("EMI-MERGE MEMORY GATE: %d FAILURES (%s)\n", failures, dir.c_str());
    return 1;
}
