#include "dart/damage_index.hpp"
#include "dart/fdr_estimator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <sys/wait.h>
#include <unistd.h>

namespace {

struct Table {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
};

static std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out += c;
    }
    out += "'";
    return out;
}

static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (start <= line.size()) {
        size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
}

static bool run_command_capture(const std::string& cmd, std::string& out, std::string& err) {
    std::array<char, 256> buf{};
    std::string full_cmd = cmd + " 2>&1";
    FILE* pipe = popen(full_cmd.c_str(), "r");
    if (!pipe) {
        err = "popen failed";
        return false;
    }
    while (fgets(buf.data(), static_cast<int>(buf.size()), pipe)) {
        out += buf.data();
    }
    int rc = pclose(pipe);
    if (rc == -1) {
        err = "pclose failed";
        return false;
    }
    if (!WIFEXITED(rc) || WEXITSTATUS(rc) != 0) {
        std::ostringstream oss;
        oss << "command failed (" << rc << "): " << cmd << "\n" << out;
        err = oss.str();
        return false;
    }
    return true;
}

static bool convert_hits_to_emi(const std::string& dart_bin,
                                const std::string& hits_path,
                                const std::string& emi_path,
                                std::string& err) {
    std::string out;
    std::string cmd = shell_quote(dart_bin) + " hits2emi -i " + shell_quote(hits_path) +
                      " -o " + shell_quote(emi_path) + " --threads 2";
    return run_command_capture(cmd, out, err);
}

static bool parse_table(const std::string& text, Table& t, std::string& err) {
    std::istringstream iss(text);
    std::string line;
    bool have_header = false;
    while (std::getline(iss, line)) {
        // Captured output interleaves stderr log lines (e.g. the libtaph gate warning); only
        // tab-separated lines belong to the table.
        if (line.find('\t') == std::string::npos) continue;
        if (!have_header) {
            t.header = split_tabs(line);
            have_header = true;
            continue;
        }
        t.rows.push_back(split_tabs(line));
    }
    if (!have_header) {
        err = "missing header";
        return false;
    }
    return true;
}

static std::unordered_map<std::string, std::string>
row_to_map(const Table& t, const std::vector<std::string>& row) {
    std::unordered_map<std::string, std::string> m;
    const size_t n = std::min(t.header.size(), row.size());
    for (size_t i = 0; i < n; ++i) {
        m[t.header[i]] = row[i];
    }
    return m;
}

static bool write_file(const std::string& path, const std::string& content, std::string& err) {
    std::ofstream out(path);
    if (!out.is_open()) {
        err = "cannot open " + path;
        return false;
    }
    out << content;
    if (!out) {
        err = "write failed for " + path;
        return false;
    }
    return true;
}

static bool read_file(const std::string& path, std::string& content, std::string& err) {
    std::ifstream in(path);
    if (!in.is_open()) {
        err = "cannot open " + path;
        return false;
    }
    std::ostringstream oss;
    oss << in.rdbuf();
    if (!in.good() && !in.eof()) {
        err = "read failed for " + path;
        return false;
    }
    content = oss.str();
    return true;
}


static bool get_row_for_query(
    const Table& t,
    const std::string& query_id,
    std::unordered_map<std::string, std::string>& out) {
    for (const auto& row : t.rows) {
        auto m = row_to_map(t, row);
        auto it = m.find("query_id");
        if (it != m.end() && it->second == query_id) {
            out = std::move(m);
            return true;
        }
    }
    return false;
}

static bool expect(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "FAILED: " << msg << "\n";
        return false;
    }
    return true;
}

static bool test_em_target_gamma_consistency(const std::string& dart_bin, const std::string& tmpdir) {
    std::cout << "Running EM target/gamma consistency regression...\n";
    std::ostringstream hits;
    hits << "readX_+_0\ttargetA\t0.90\t4\t1\t0\t1\t4\t1\t4\t1e-20\t100\t4\t100\tWCDE\tRCDE\n";
    hits << "readX_+_0\ttargetB\t0.90\t4\t1\t0\t1\t4\t1\t4\t1e-20\t99\t4\t100\tWCDE\tRCDE\n";
    for (int i = 1; i <= 200; ++i) {
        hits << "readB" << i << "_+_0\ttargetB\t0.90\t4\t0\t0\t1\t4\t1\t4\t1e-20\t100\t4\t100\tACDE\tACDE\n";
    }
    std::string hits_path = tmpdir + "/em_hits.tsv";
    std::string err;
    if (!write_file(hits_path, hits.str(), err)) {
        std::cerr << err << "\n";
        return false;
    }
    const std::string emi_path = tmpdir + "/em_hits.emi2";
    if (!convert_hits_to_emi(dart_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string out;
    std::string cmd = shell_quote(dart_bin) + " damage-annotate --no-null --emi " + shell_quote(emi_path) + " --em";
    if (!run_command_capture(cmd, out, err)) {
        std::cerr << err << "\n";
        return false;
    }

    Table t;
    if (!parse_table(out, t, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::unordered_map<std::string, std::string> row;
    if (!get_row_for_query(t, "readX_+_0", row)) {
        std::cerr << "readX_+_0 not found in output\n";
        return false;
    }

    double gamma = std::stod(row["gamma"]);
    int best_hit = std::stoi(row["best_hit"]);

    if (!expect(best_hit == 0, "best_hit should be 0 when EM prefers a different target")) return false;
    if (!expect(gamma > 0.0, "gamma should be positive for reported target")) return false;
    if (!expect(gamma < 0.9, "gamma should reflect reported target, not EM-best alternative")) return false;
    return true;
}


static bool test_multimap_summary_invariance(const std::string& dart_bin, const std::string& tmpdir) {
    std::cout << "Running multimap summary invariance regression...\n";
    std::ostringstream hits;
    hits << "shared_+_0\tprotA\t0.833\t9\t1\t0\t1\t6\t1\t6\t1e-20\t100\t6\t100\tWAAAAA\tRAAAAA\n";
    hits << "shared_+_0\tprotB\t0.833\t9\t1\t0\t1\t6\t1\t6\t1e-20\t99\t6\t100\tWAAAAA\tRAAAAA\n";
    for (int i = 0; i < 40; ++i)
        hits << "anchor" << i << "_+_0\tprotB\t1\t6\t0\t0\t1\t6\t1\t6\t1e-20\t100\t6\t100\tRAAAAA\tRAAAAA\n";
    for (int i = 0; i < 5; ++i)
        hits << "anchorA" << i << "_+_0\tprotA\t1\t6\t0\t0\t1\t6\t1\t6\t1e-20\t100\t6\t100\tRAAAAA\tRAAAAA\n";
    std::string err;
    const std::string emi = tmpdir + "/multimap.emi";
    if (!write_file(tmpdir + "/multimap.tsv", hits.str(), err) ||
        !convert_hits_to_emi(dart_bin, tmpdir + "/multimap.tsv", emi, err)) {
        std::cerr << err << "\n";
        return false;
    }
    std::array<std::string, 2> expected;
    for (int mode = 0; mode < 4; ++mode) {
        const std::string prefix = tmpdir + "/multimap_" + std::to_string(mode);
        std::string out;
        const std::string cmd = shell_quote(dart_bin) + " damage-annotate --no-null -i " + shell_quote(emi) +
            " --d-max 0.3 --min-reads 1 --min-breadth 0 --min-depth 0"
            " --min-positional-score 0 --min-terminal-ratio 0 -t " + ((mode & 1) ? "4" : "1") +
            ((mode & 2) ? " --emit-multimap --multimap-k 1" : "") +
            " -o " + shell_quote(prefix + ".reads") +
            " --protein-summary " + shell_quote(prefix + ".protein") +
            " --gene-summary " + shell_quote(prefix + ".gene");
        if (!run_command_capture(cmd, out, err)) {
            std::cerr << err << "\n";
            return false;
        }
        int table = 0;
        for (const auto& suffix : {".protein", ".gene"}) {
            std::string text;
            if (!read_file(prefix + suffix, text, err)) return false;
            if (mode == 0) expected[table] = text;
            else if (!expect(text == expected[table],
                             "multimap/thread mode " + std::to_string(mode) + suffix + " must match best-hit summary"))
                return false;
            ++table;
        }
    }
    return true;
}

static bool test_protein_mapping_filters(const std::string& dart_bin, const std::string& tmpdir) {
    std::cout << "Running protein mapping filter regression...\n";
    const std::string hits_path = tmpdir + "/mapping_hits.tsv";
    const std::string qaln = "AAAAAAAAAAAAAAAAAAAF";
    const std::string taln = "AAAAAAAAAAAAAAAAAAAE";

    std::ostringstream hits;
    // Protein with broad starts + terminal support (should pass).
    hits << "good_r1_+_0\tprot_good\t0.95\t20\t1\t0\t1\t20\t1\t20\t1e-20\t120\t20\t100\t" << qaln << "\t" << taln << "\n";
    hits << "good_r2_+_0\tprot_good\t0.95\t20\t1\t0\t1\t20\t80\t99\t1e-20\t120\t20\t100\t" << qaln << "\t" << taln << "\n";
    hits << "good_r3_+_0\tprot_good\t0.95\t20\t1\t0\t1\t20\t40\t59\t1e-20\t120\t20\t100\t" << qaln << "\t" << taln << "\n";
    // Protein with collapsed middle-only starts (should fail).
    hits << "bad_r1_+_0\tprot_bad\t0.95\t20\t1\t0\t1\t20\t40\t59\t1e-20\t120\t20\t100\t" << qaln << "\t" << taln << "\n";
    hits << "bad_r2_+_0\tprot_bad\t0.95\t20\t1\t0\t1\t20\t40\t59\t1e-20\t120\t20\t100\t" << qaln << "\t" << taln << "\n";
    hits << "bad_r3_+_0\tprot_bad\t0.95\t20\t1\t0\t1\t20\t40\t59\t1e-20\t120\t20\t100\t" << qaln << "\t" << taln << "\n";

    std::string err;
    if (!write_file(hits_path, hits.str(), err)) {
        std::cerr << err << "\n";
        return false;
    }
    const std::string emi_path = tmpdir + "/mapping_hits.emi2";
    if (!convert_hits_to_emi(dart_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    const std::string gene_summary = tmpdir + "/gene_summary.tsv";
    const std::string protein_summary = tmpdir + "/protein_summary.tsv";
    std::string out;
    std::string cmd = shell_quote(dart_bin) + " damage-annotate --no-null --emi " + shell_quote(emi_path) +
                      " --gene-summary " + shell_quote(gene_summary) +
                      " --protein-summary " + shell_quote(protein_summary) +
                      " --min-reads 1 --min-breadth 0 --min-depth 0 " +
                      " --min-positional-score 0.4 --min-terminal-ratio 0.5";
    if (!run_command_capture(cmd, out, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string gene_txt, protein_txt;
    if (!read_file(gene_summary, gene_txt, err) ||
        !read_file(protein_summary, protein_txt, err)) {
        std::cerr << err << "\n";
        return false;
    }

    Table tg, tp;
    if (!parse_table(gene_txt, tg, err) || !parse_table(protein_txt, tp, err)) {
        std::cerr << err << "\n";
        return false;
    }

    bool gene_has_good = false, gene_has_bad = false;
    for (const auto& row : tg.rows) {
        auto m = row_to_map(tg, row);
        if (m["target_id"] == "prot_good") gene_has_good = true;
        if (m["target_id"] == "prot_bad") gene_has_bad = true;
    }
    if (!expect(gene_has_good, "gene summary should keep prot_good")) return false;
    if (!expect(!gene_has_bad, "gene summary should filter prot_bad by mapping thresholds")) return false;

    bool good_pass = false, bad_fail = false;
    for (const auto& row : tp.rows) {
        auto m = row_to_map(tp, row);
        if (m["protein_id"] == "prot_good" && m["pass_mapping_filter"] == "1") good_pass = true;
        if (m["protein_id"] == "prot_bad" && m["pass_mapping_filter"] == "0") bad_fail = true;
    }
    if (!expect(good_pass, "protein summary should mark prot_good as pass_mapping_filter=1")) return false;
    if (!expect(bad_fail, "protein summary should mark prot_bad as pass_mapping_filter=0")) return false;

    return true;
}


static bool test_refine_damage_zero_mismatch_denominator(const std::string& dart_bin, const std::string& tmpdir) {
    std::cout << "Running --refine-damage zero-mismatch denominator regression...\n";

    // Create synthetic data:
    //   5 reads with damage-consistent R->W at AA position 0 (5' terminus)
    //   5 reads with no mismatches but same susceptible R in target at position 0
    // All reads map to the same target.
    //
    // With d_max=1.0, lambda=0.3, decay at dist=0 is exp(-0.3*0)=1.0.
    // Correct MLE (denominator includes all reads):
    //   sum_gamma_k=5, sum_gamma_decay=10 -> d_max_hat = 0.5
    // Biased estimate (denominator from mismatch reads only):
    //   sum_gamma_k=5, sum_gamma_decay=5  -> d_max_hat = 1.0
    //
    // This regression catches upward bias from excluding zero-mismatch reads.

    const std::string hits_path = tmpdir + "/refine_hits.tsv";
    std::ostringstream hits;
    for (int i = 1; i <= 5; ++i) {
        // R->W damage-consistent substitution at position 0 (dist_5=0)
        hits << "hitR_" << i << "_+_0\tprotR\t0.833\t6\t0\t0\t1\t6\t1\t6\t1e-10\t60\t6\t10\t"
             << "WAAAAA\tRAAAAA\n";
    }
    for (int i = 1; i <= 5; ++i) {
        // No mismatches, but R at target position 0 is a susceptible site
        hits << "noHit_" << i << "_+_0\tprotR\t1.00\t6\t0\t0\t1\t6\t1\t6\t1e-10\t60\t6\t10\t"
             << "RAAAAA\tRAAAAA\n";
    }

    std::string err;
    if (!write_file(hits_path, hits.str(), err)) {
        std::cerr << err << "\n";
        return false;
    }
    const std::string emi_path = tmpdir + "/refine_hits.emi2";
    if (!convert_hits_to_emi(dart_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string out;
    // Use d_max=1.0 so any downward refit proves zero-mismatch reads enter the denominator.
    std::string cmd = shell_quote(dart_bin) + " damage-annotate --no-null --emi " + shell_quote(emi_path) +
                      " --d-max 1.0 --lambda 0.3 --refine-damage -v";
    if (!run_command_capture(cmd, out, err)) {
        std::cerr << err << "\n";
        return false;
    }

    const std::string marker = "[refine-damage]";
    size_t pos = out.find(marker);
    if (!expect(pos != std::string::npos, "--refine-damage should emit a [refine-damage] verbose line")) {
        return false;
    }

    // Parse d_max_new from "d_max: <old> -> <new> ..."
    size_t arrow = out.find("-> ", pos);
    if (!expect(arrow != std::string::npos, "[refine-damage] line should contain '->'")) return false;
    double d_max_new = std::stod(out.substr(arrow + 3));

    // Unbiased estimate: 5 hits / 10 decays = 0.5
    // Biased estimate (old behaviour): 5 hits / 5 decays = 1.0
    // We allow generous headroom for EM weighting effects.
    if (!expect(d_max_new < 0.8,
                "refine-damage d_max should be < 0.8 (zero-mismatch reads in denominator)")) {
        std::cerr << "  got d_max_new=" << d_max_new << " (expected ~0.5)\n";
        return false;
    }
    if (!expect(d_max_new > 0.0, "refine-damage d_max should be > 0")) return false;
    return true;
}

static bool test_ancient_threshold_cli_flag(const std::string& dart_bin, const std::string& tmpdir) {
    std::cout << "Running --ancient-threshold CLI flag regression...\n";
    const std::string hits_path = tmpdir + "/athresh_hits.tsv";
    // One R->W damage-consistent substitution at the terminal position.
    // Site-only posterior (no AGD p_read) lands ~0.38 — above modern_threshold (0.25)
    // but below the default ancient_threshold (0.60), so classification is Undetermined.
    // With --ancient-threshold 0.20 the same posterior triggers AncientConfident.
    std::ostringstream hits;
    hits << "athresh_r1_+_0\tatarget\t0.90\t6\t1\t0\t1\t6\t1\t6\t1e-10\t40\t6\t100\tWAAAAA\tRAAAAA\n";
    std::string err;
    if (!write_file(hits_path, hits.str(), err)) {
        std::cerr << err << "\n";
        return false;
    }
    const std::string emi_path = tmpdir + "/athresh_hits.emi2";
    if (!convert_hits_to_emi(dart_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    // Default threshold (0.60): site-only posterior ~0.38 -> Undetermined -> is_damaged=0
    std::string out_default, out_low;
    std::string cmd_default = shell_quote(dart_bin) + " damage-annotate --no-null --emi " + shell_quote(emi_path) +
                              " --d-max 0.3";
    if (!run_command_capture(cmd_default, out_default, err)) {
        std::cerr << err << "\n";
        return false;
    }

    // Lower ancient-threshold to 0.20 (and modern-threshold to 0.10 to satisfy ordering constraint):
    // same posterior should become damage_class="ancient"
    std::string cmd_low = shell_quote(dart_bin) + " damage-annotate --no-null --emi " + shell_quote(emi_path) +
                          " --d-max 0.3 --ancient-threshold 0.20 --modern-threshold 0.10";
    if (!run_command_capture(cmd_low, out_low, err)) {
        std::cerr << err << "\n";
        return false;
    }

    Table t_default, t_low;
    if (!parse_table(out_default, t_default, err) || !parse_table(out_low, t_low, err)) {
        std::cerr << err << "\n";
        return false;
    }
    std::unordered_map<std::string, std::string> r_default, r_low;
    if (!get_row_for_query(t_default, "athresh_r1_+_0", r_default) ||
        !get_row_for_query(t_low, "athresh_r1_+_0", r_low)) {
        std::cerr << "athresh_r1_+_0 missing in output\n";
        return false;
    }

    const double posterior = std::stod(r_default["posterior"]);
    if (!expect(posterior > 0.20 && posterior < 0.60,
                "default run: posterior should be in (0.20, 0.60)")) {
        std::cerr << "  got posterior=" << posterior << "\n";
        return false;
    }
    // With default ancient_threshold=0.60, posterior ~0.38 -> damage_class="uncertain"
    if (!expect(r_default["damage_class"] == "uncertain",
                "default threshold (0.60): damage_class should be 'uncertain'")) return false;

    // With --ancient-threshold 0.20, same posterior (>0.20) -> damage_class="ancient"
    if (!expect(r_low["damage_class"] == "ancient",
                "--ancient-threshold 0.20: damage_class should be 'ancient'")) return false;
    return true;
}

static bool test_fdr_scan_non_monotone() {
    std::cout << "Running FDR scan non-monotone regression (unit test)...\n";

    // Construct a non-monotone FDR scenario:
    //   - 3 targets at high scores (0.9, 0.8, 0.7)
    //   - 1 decoy with score 0.65 (spike causing FDR > target right at 0.65)
    //   - 100 targets at score 0.01 (dilute the decoy at low thresholds)
    //
    // At target_fdr=0.05, step=0.01:
    //   threshold=0.70: 3 targets, 0 decoys -> FDR=0 (best so far)
    //   threshold=0.65: 3 targets, 1 decoy  -> FDR=0.25 (exceeds target -> old code breaks)
    //   threshold=0.01: 103 targets, 1 decoy -> FDR~0.01 (below target -> correct best)
    //
    // Old code (with break) returns 0.70.
    // Fixed code (no break) returns 0.01.

    std::vector<float> target_scores = {0.9f, 0.8f, 0.7f};
    for (int i = 0; i < 100; ++i) target_scores.push_back(0.01f);
    std::vector<float> decoy_scores = {0.65f};

    const float threshold = dart::FDREstimator::find_threshold_for_fdr(
        0.05f, target_scores, decoy_scores, 0.0f, -1.0f, 0.01f);

    // Fixed code should find a threshold <= 0.01 (the low-score targets satisfy FDR<=0.05)
    if (!expect(threshold <= 0.02f,
                "FDR scan should find threshold ~0.01 after decoy spike (no early break)")) {
        std::cerr << "  got threshold=" << threshold << " (expected <=0.02)\n";
        return false;
    }
    // Verify FDR is actually satisfied at the returned threshold
    const float fdr_at_result = dart::FDREstimator::estimate_fdr(
        threshold, target_scores, decoy_scores);
    if (!expect(fdr_at_result <= 0.05f,
                "returned threshold should satisfy target FDR=0.05")) {
        std::cerr << "  fdr=" << fdr_at_result << " at threshold=" << threshold << "\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path-to-dart>\n";
        return 2;
    }
    const std::string dart_bin = argv[1];

    const char* base = std::getenv("TMPDIR");
    std::string pattern = std::string(base && *base ? base : ".") + "/dart_regress_XXXXXX";
    std::vector<char> tmp_template(pattern.begin(), pattern.end());
    tmp_template.push_back('\0');
    char* tmp = mkdtemp(tmp_template.data());
    if (!tmp) {
        std::cerr << "Failed to create temp dir\n";
        return 2;
    }
    const std::string tmpdir = tmp;
    // The EM tape needs an O_DIRECT-capable local dir; keep a caller's DART_EM_TAPE_DIR.
    setenv("DART_EM_TAPE_DIR", (tmpdir + "/emtape").c_str(), 0);

    bool ok = true;
    if (argc == 3 && std::string(argv[2]) == "--multimap-only")
        return test_multimap_summary_invariance(dart_bin, tmpdir) ? 0 : 1;
    ok = test_multimap_summary_invariance(dart_bin, tmpdir) && ok;
    ok = test_em_target_gamma_consistency(dart_bin, tmpdir) && ok;
    ok = test_protein_mapping_filters(dart_bin, tmpdir) && ok;
    ok = test_refine_damage_zero_mismatch_denominator(dart_bin, tmpdir) && ok;
    ok = test_ancient_threshold_cli_flag(dart_bin, tmpdir) && ok;
    ok = test_fdr_scan_non_monotone() && ok;

    if (!ok) return 1;
    std::cout << "All damage-annotate regressions passed.\n";
    return 0;
}
