/**
 * Regression tests for `dart profile`, which runs the libtaph library profiler
 * (the same one as `fqdup profile`) and writes its JSON and .tdmg.
 *
 * Fixtures:
 *   DS damaged  — tutorial reads
 *   SS damaged  — dart simulate --library-type ss
 *   Modern      — random reads (uniform composition, no damage signal)
 *
 * Usage: test-profile-regressions <path-to-dart> <path-to-tutorial-reads.fq.gz>
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace {

std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out += c;
    }
    out += "'";
    return out;
}

bool run_command(const std::string& cmd, std::string& out, std::string& err_msg) {
    std::array<char, 256> buf{};
    FILE* pipe = popen((cmd + " 2>&1").c_str(), "r");
    if (!pipe) { err_msg = "popen failed"; return false; }
    std::ostringstream oss;
    while (fgets(buf.data(), buf.size(), pipe)) oss << buf.data();
    int rc = pclose(pipe);
    out = oss.str();
    if (WIFEXITED(rc) && WEXITSTATUS(rc) == 0) return true;
    err_msg = "exit " + std::to_string(WEXITSTATUS(rc)) + ": " + out;
    return false;
}

bool expect(bool cond, const std::string& msg) {
    if (!cond) std::cerr << "FAILED: " << msg << "\n";
    return cond;
}

std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

// Number after "key": at the first occurrence of `key` past `anchor`; NaN if absent.
double number_after(const std::string& json, const std::string& anchor, const std::string& key) {
    size_t p = json.find(anchor);
    if (p == std::string::npos) return std::nan("");
    p = json.find("\"" + key + "\":", p);
    if (p == std::string::npos) return std::nan("");
    return std::strtod(json.c_str() + p + key.size() + 3, nullptr);
}

bool has(const std::string& json, const std::string& s) { return json.find(s) != std::string::npos; }

void write_random_fastq(const std::string& path, int n, int len, uint64_t seed) {
    std::ofstream out(path);
    uint64_t rng = seed;
    for (int i = 0; i < n; ++i) {
        out << "@rnd_" << i << "\n";
        for (int j = 0; j < len; ++j) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            out << "ACGT"[(rng >> 33) & 3];
        }
        out << "\n+\n" << std::string(len, 'I') << "\n";
    }
}

// Runs `dart profile` on `fq` with `threads` workers; the JSON lands in `json`.
bool profile(const std::string& dart, const std::string& fq, const std::string& out_prefix,
             int threads, std::string& json) {
    std::string out, err;
    const std::string cmd = shell_quote(dart) + " profile -i " + shell_quote(fq) +
        " -t " + std::to_string(threads) + " -o " + shell_quote(out_prefix + ".json") +
        " --damage-summary " + shell_quote(out_prefix + ".tdmg");
    if (!run_command(cmd, out, err)) { std::cerr << "FAILED: " << cmd << ": " << err << "\n"; return false; }
    json = slurp(out_prefix + ".json");
    return true;
}

const std::string kDeamVerdict = "\"verdict\": {\n      \"deamination\":";

bool test_ds_damaged(const std::string& dart, const std::string& fq, const std::string& tmp) {
    std::string json;
    if (!profile(dart, fq, tmp + "/ds", 4, json)) return false;
    bool ok = true;
    ok = expect(has(json, "\"library_type\": \"double-stranded\""), "tutorial reads are double-stranded") && ok;
    const double d5 = number_after(json, kDeamVerdict, "d_max_5prime");
    ok = expect(d5 > 0.05, "tutorial 5' d_max > 0.05 (got " + std::to_string(d5) + ")") && ok;
    ok = expect(!slurp(tmp + "/ds.tdmg").empty(), ".tdmg written") && ok;
    if (ok) std::cout << "[PASS] ds_damaged\n";
    return ok;
}

bool test_ss_damaged(const std::string& dart, const std::string& tmp) {
    const std::string fq = tmp + "/ss.fq";
    std::string out, err;
    if (!run_command(shell_quote(dart) + " simulate -n 20000 --library-type ss --nv 0.05 --lambda 0.25"
                     " --delta-s 0.90 --delta-d 0.02 --seed 42 -o " + shell_quote(fq), out, err)) {
        std::cerr << "FAILED: simulate: " << err << "\n";
        return false;
    }
    std::string json;
    if (!profile(dart, fq, tmp + "/ss", 4, json)) return false;
    bool ok = true;
    ok = expect(has(json, "\"library_type\": \"single-stranded\""), "simulated ss library is single-stranded") && ok;
    const double d5 = number_after(json, kDeamVerdict, "d_max_5prime");
    ok = expect(d5 > 0.05, "ss 5' d_max > 0.05 (got " + std::to_string(d5) + ")") && ok;
    if (ok) std::cout << "[PASS] ss_damaged\n";
    return ok;
}

bool test_modern(const std::string& dart, const std::string& tmp) {
    const std::string fq = tmp + "/modern.fq";
    write_random_fastq(fq, 20000, 80, 7);
    std::string json;
    if (!profile(dart, fq, tmp + "/modern", 4, json)) return false;
    bool ok = true;
    const double d5 = number_after(json, kDeamVerdict, "d_max_5prime");
    ok = expect(d5 < 0.05, "random reads 5' d_max < 0.05 (got " + std::to_string(d5) + ")") && ok;
    ok = expect(!has(json, "\"state\": \"DETECTED\""), "random reads are not DETECTED") && ok;
    if (ok) std::cout << "[PASS] modern\n";
    return ok;
}

// Output does not depend on the thread count, and the JSON carries the sections and
// the libtaph SHA the profile came from.
bool test_schema_and_threads(const std::string& dart, const std::string& fq, const std::string& tmp) {
    std::string j1, j8, version, err;
    if (!profile(dart, fq, tmp + "/t1", 1, j1) || !profile(dart, fq, tmp + "/t8", 8, j8)) return false;
    bool ok = true;
    ok = expect(j1 == j8, "JSON identical with 1 and 8 threads") && ok;
    ok = expect(slurp(tmp + "/t1.tdmg") == slurp(tmp + "/t8.tdmg"), ".tdmg identical with 1 and 8 threads") && ok;
    for (const char* section : {"\"identity\"", "\"qc\"", "\"library_call\"", "\"channels\"",
                                "\"characterization\"", "\"diagnostics\""})
        ok = expect(has(j1, section), std::string("section ") + section) && ok;
    if (run_command(shell_quote(dart) + " --version", version, err)) {
        const size_t p = version.find("libtaph ");
        const std::string sha = p == std::string::npos ? "" : version.substr(p + 8, 40);
        ok = expect(!sha.empty() && has(j1, "\"libtaph_git_sha\": \"" + sha + "\""),
                    "identity.libtaph_git_sha matches dart --version") && ok;
    }
    if (ok) std::cout << "[PASS] schema_and_threads\n";
    return ok;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <path-to-dart> <path-to-tutorial-reads.fq.gz>\n";
        return 2;
    }
    const std::string dart = argv[1], tutorial = argv[2];
    const char* base = std::getenv("TMPDIR");
    std::string pattern = std::string(base && *base ? base : ".") + "/dart_profile_regress_XXXXXX";
    std::vector<char> tmpl(pattern.begin(), pattern.end());
    tmpl.push_back('\0');
    if (!mkdtemp(tmpl.data())) { std::cerr << "Failed to create temp dir\n"; return 2; }
    const std::string tmp = tmpl.data();

    bool ok = true;
    ok = test_ds_damaged(dart, tutorial, tmp) && ok;
    ok = test_ss_damaged(dart, tmp) && ok;
    ok = test_modern(dart, tmp) && ok;
    ok = test_schema_and_threads(dart, tutorial, tmp) && ok;
    if (!ok) return 1;
    std::cout << "All profile regressions passed.\n";
    return 0;
}
