/**
 * Regression tests for `dart profile` JSON output.
 *
 * Fixture strategy:
 *   DS damaged  — tutorial reads (real ancient DNA, no ORF composition bias)
 *   SS damaged  — dart simulate --library-type ss (verifies SS classification path)
 *   Modern      — inline LCG random reads (uniform composition, no damage signal)
 *
 * dart simulate embeds ORFs which create strong codon-phase composition bias at
 * read termini, saturating pos0 and confusing library-type detection. It is
 * only usable here for SS reads where the classification is robust to that bias.
 *
 * Usage: test-profile-regressions <path-to-dart> <path-to-tutorial-reads.fq.gz>
 */

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

namespace {

static std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out += c;
    }
    out += "'";
    return out;
}

static bool run_command(const std::string& cmd, std::string& out, std::string& err_msg) {
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

static bool expect(bool cond, const std::string& msg) {
    if (!cond) { std::cerr << "FAILED: " << msg << "\n"; return false; }
    return true;
}

// ── JSON field extractors ────────────────────────────────────────────────────

static bool get_json_double(const std::string& json, const std::string& key, double& val) {
    std::string pat = "\"" + key + "\": ";
    size_t p = json.find(pat);
    if (p == std::string::npos) return false;
    try { val = std::stod(json.substr(p + pat.size())); return true; }
    catch (...) { return false; }
}

static bool get_json_bool(const std::string& json, const std::string& key, bool& val) {
    std::string pat = "\"" + key + "\": ";
    size_t p = json.find(pat);
    if (p == std::string::npos) return false;
    p += pat.size();
    if (json.substr(p, 4) == "true")  { val = true;  return true; }
    if (json.substr(p, 5) == "false") { val = false; return true; }
    return false;
}

static bool get_json_string(const std::string& json, const std::string& key, std::string& val) {
    std::string pat = "\"" + key + "\": \"";
    size_t p = json.find(pat);
    if (p == std::string::npos) return false;
    p += pat.size();
    size_t e = json.find('"', p);
    if (e == std::string::npos) return false;
    val = json.substr(p, e - p);
    return true;
}

// ── Modern fixture: inline random-sequence FASTQ ────────────────────────────
// Generates reads with uniform base composition (no ORF, no damage signal).
// LCG parameters from Knuth MMIX.

static bool write_random_fastq(const std::string& path, int n, int len, uint64_t seed, std::string& err) {
    std::ofstream out(path);
    if (!out) { err = "cannot open " + path; return false; }
    uint64_t rng = seed;
    auto next_base = [&]() -> char {
        rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
        return "ACGT"[(rng >> 33) & 3];
    };
    for (int i = 0; i < n; ++i) {
        out << "@rnd_" << i << "\n";
        for (int j = 0; j < len; ++j) out << next_base();
        out << "\n+\n";
        out << std::string(len, 'I') << "\n";
    }
    if (!out) { err = "write error"; return false; }
    return true;
}

static bool run_profile(const std::string& dart_bin, const std::string& fq,
                        std::string& json, std::string& err) {
    return run_command(shell_quote(dart_bin) + " profile " + shell_quote(fq), json, err);
}

// ── Test cases ───────────────────────────────────────────────────────────────

// DS damaged: tutorial reads are real DS ancient DNA — no pos0 artifact, clean
// terminal composition, hex_end_asymmetry reliably computed.
static bool test_profile_ds_damaged(const std::string& dart_bin,
                                    const std::string& tutorial_reads) {
    std::cout << "Running profile_ds_damaged...\n";
    std::string json, err;
    if (!run_profile(dart_bin, tutorial_reads, json, err)) {
        std::cerr << "[FAIL] ds_damaged: profile: " << err << "\n";
        return false;
    }

    bool ok = true;
    double d_max; bool damage_validated, d5_flat, d5_blunting; std::string lib_type;

    ok = expect(get_json_double(json, "d_max", d_max),             "d_max present")    && ok;
    ok = expect(d_max > 5.0, "d_max > 5% (got " + std::to_string(d_max) + ")")         && ok;
    ok = expect(get_json_bool(json, "damage_validated", damage_validated),
                "damage_validated present")                                              && ok;
    ok = expect(damage_validated,                                  "damage_validated")  && ok;
    ok = expect(get_json_string(json, "library_type", lib_type),  "library_type present") && ok;
    ok = expect(lib_type == "double-stranded",                    "library_type DS")   && ok;

    // v1.2.0 flatness/blunting (real ancient DS reads should show real decay, not flat)
    ok = expect(json.find("\"d5_profile_flat\"") != std::string::npos,
                "d5_profile_flat key present")                                          && ok;
    ok = expect(json.find("\"d3_profile_flat\"") != std::string::npos,
                "d3_profile_flat key present")                                          && ok;
    ok = expect(json.find("\"d5_blunting_suspected\"") != std::string::npos,
                "d5_blunting_suspected key present")                                    && ok;
    ok = expect(json.find("\"d5_max_rate_pos0_4\"") != std::string::npos,
                "d5_max_rate_pos0_4 key present")                                       && ok;
    ok = expect(json.find("\"d3_max_rate_pos0_4\"") != std::string::npos,
                "d3_max_rate_pos0_4 key present")                                       && ok;
    if (get_json_bool(json, "d5_profile_flat", d5_flat))
        ok = expect(!d5_flat,    "d5_profile_flat false for damaged DS")                && ok;
    if (get_json_bool(json, "d5_blunting_suspected", d5_blunting))
        ok = expect(!d5_blunting,"d5_blunting_suspected false for normal DS")           && ok;

    // Tutorial reads have no pos0 artifact → hex_end_asymmetry must be present, not null
    ok = expect(json.find("\"position_0_artifact_5prime\": false") != std::string::npos,
                "no pos0 artifact for tutorial reads")                                  && ok;
    ok = expect(json.find("\"hex_end_asymmetry\": {") != std::string::npos,
                "hex_end_asymmetry object present (not null)")                          && ok;
    ok = expect(json.find("\"rc_excess_jsd\"") != std::string::npos,
                "rc_excess_jsd present")                                                && ok;

    if (ok) std::cout << "[PASS] profile_ds_damaged\n";
    return ok;
}

// SS damaged: simulate reads reliably classify as SS and validate damage
// despite ORF-induced pos0 artifact (hex_end_asymmetry will be null).
static bool test_profile_ss_damaged(const std::string& dart_bin, const std::string& tmpdir) {
    std::cout << "Running profile_ss_damaged...\n";
    std::string fq = tmpdir + "/ss_damaged.fq";
    std::string cmd = shell_quote(dart_bin)
        + " simulate -n 20000 --library-type ss"
        + " --nv 0.05 --lambda 0.25 --delta-s 0.90 --delta-d 0.02"
        + " --seed 42 -o " + shell_quote(fq);
    std::string out, err;
    if (!run_command(cmd, out, err)) {
        std::cerr << "[FAIL] ss_damaged: simulate: " << err << "\n";
        return false;
    }
    std::string json;
    if (!run_profile(dart_bin, fq, json, err)) {
        std::cerr << "[FAIL] ss_damaged: profile: " << err << "\n";
        return false;
    }

    bool ok = true;
    double d_max; bool damage_validated, d5_flat; std::string lib_type;

    ok = expect(get_json_double(json, "d_max", d_max),            "d_max present")     && ok;
    ok = expect(d_max > 5.0, "d_max > 5% (got " + std::to_string(d_max) + ")")         && ok;
    ok = expect(get_json_bool(json, "damage_validated", damage_validated),
                "damage_validated present")                                              && ok;
    ok = expect(damage_validated,                                 "damage_validated")   && ok;
    ok = expect(get_json_string(json, "library_type", lib_type), "library_type present") && ok;
    ok = expect(lib_type == "single-stranded",                   "library_type SS")    && ok;

    ok = expect(json.find("\"d5_profile_flat\"") != std::string::npos,
                "d5_profile_flat key present")                                          && ok;
    ok = expect(json.find("\"d5_blunting_suspected\"") != std::string::npos,
                "d5_blunting_suspected key present")                                    && ok;
    if (get_json_bool(json, "d5_profile_flat", d5_flat))
        ok = expect(!d5_flat, "d5_profile_flat false for damaged SS")                   && ok;

    if (ok) std::cout << "[PASS] profile_ss_damaged\n";
    return ok;
}

// Modern: uniform-composition random reads have no terminal enrichment.
static bool test_profile_modern(const std::string& dart_bin, const std::string& tmpdir) {
    std::cout << "Running profile_modern...\n";
    std::string fq = tmpdir + "/modern.fq";
    std::string err;
    if (!write_random_fastq(fq, 20000, 80, 0xdeadbeefULL, err)) {
        std::cerr << "[FAIL] modern: generate reads: " << err << "\n";
        return false;
    }
    std::string json;
    if (!run_profile(dart_bin, fq, json, err)) {
        std::cerr << "[FAIL] modern: profile: " << err << "\n";
        return false;
    }

    bool ok = true;
    double d_max; bool damage_validated;

    ok = expect(get_json_double(json, "d_max", d_max),             "d_max present")    && ok;
    ok = expect(d_max < 15.0, "d_max < 15% for undamaged (got " + std::to_string(d_max) + ")") && ok;
    ok = expect(get_json_bool(json, "damage_validated", damage_validated),
                "damage_validated present")                                              && ok;
    ok = expect(!damage_validated,                                 "damage_validated false") && ok;

    // Flatness fields must be present regardless of sample type
    ok = expect(json.find("\"d5_profile_flat\"") != std::string::npos,
                "d5_profile_flat key present")                                          && ok;
    ok = expect(json.find("\"d3_profile_flat\"") != std::string::npos,
                "d3_profile_flat key present")                                          && ok;
    ok = expect(json.find("\"d5_max_rate_pos0_4\"") != std::string::npos,
                "d5_max_rate_pos0_4 key present")                                       && ok;

    double rate;
    if (get_json_double(json, "d5_max_rate_pos0_4", rate))
        ok = expect(rate < 0.10, "d5_max_rate_pos0_4 low for modern (got "
                    + std::to_string(rate) + ")")                                       && ok;

    if (ok) std::cout << "[PASS] profile_modern\n";
    return ok;
}

// Schema completeness: verify every v1.2.0 field is in the JSON output.
// Uses tutorial reads (rerun — no shared state with the DS test).
static bool test_profile_json_schema(const std::string& dart_bin,
                                     const std::string& tutorial_reads) {
    std::cout << "Running profile_json_schema...\n";
    std::string json, err;
    if (!run_profile(dart_bin, tutorial_reads, json, err)) {
        std::cerr << "[FAIL] json_schema: " << err << "\n";
        return false;
    }

    bool ok = true;
    // Core fields (pre-v1.2)
    for (const char* key : {
        "d_max", "d_max_5prime", "d_max_3prime",
        "lambda_5prime", "lambda_3prime",
        "channel_b_valid", "channel_b_llr", "channel_b3_llr",
        "damage_validated", "position_0_artifact_5prime",
        "terminal_shift_5prime", "inverted_pattern_5prime", "inverted_pattern_3prime",
        "d_max_from_channel_b", "channel_c_valid", "ox_d_max",
        "ox_damage_detected", "depurination_detected", "library_type"
    }) {
        ok = expect(json.find(std::string("\"") + key + "\"") != std::string::npos,
                    std::string("field present: ") + key) && ok;
    }
    // v1.2.0 flatness/blunting
    for (const char* key : {
        "d5_profile_flat", "d3_profile_flat", "d5_blunting_suspected",
        "d5_max_rate_pos0_4", "d3_max_rate_pos0_4"
    }) {
        ok = expect(json.find(std::string("\"") + key + "\"") != std::string::npos,
                    std::string("v1.2 field present: ") + key) && ok;
    }
    ok = expect(json.find("\"hex_end_asymmetry\"") != std::string::npos,
                "v1.2 field present: hex_end_asymmetry") && ok;

    if (ok) std::cout << "[PASS] profile_json_schema\n";
    return ok;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <path-to-dart> <path-to-tutorial-reads.fq.gz>\n";
        return 2;
    }
    const std::string dart_bin     = argv[1];
    const std::string tutorial_fq  = argv[2];

    char tmp_template[] = "/tmp/dart_profile_regress_XXXXXX";
    char* tmp = mkdtemp(tmp_template);
    if (!tmp) { std::cerr << "Failed to create temp dir\n"; return 2; }
    const std::string tmpdir = tmp;

    bool ok = true;
    ok = test_profile_ds_damaged(dart_bin, tutorial_fq)    && ok;
    ok = test_profile_ss_damaged(dart_bin, tmpdir)         && ok;
    ok = test_profile_modern(dart_bin, tmpdir)             && ok;
    ok = test_profile_json_schema(dart_bin, tutorial_fq)   && ok;

    if (!ok) return 1;
    std::cout << "All profile regressions passed.\n";
    return 0;
}
