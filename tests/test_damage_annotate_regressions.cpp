#include "agp/damage_index.hpp"

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

static bool convert_hits_to_emi(const std::string& agp_bin,
                                const std::string& hits_path,
                                const std::string& emi_path,
                                std::string& err) {
    std::string out;
    std::string cmd = shell_quote(agp_bin) + " hits2emi -i " + shell_quote(hits_path) +
                      " -o " + shell_quote(emi_path) + " --threads 2";
    return run_command_capture(cmd, out, err);
}

static bool parse_table(const std::string& text, Table& t, std::string& err) {
    std::istringstream iss(text);
    std::string line;
    bool have_header = false;
    while (std::getline(iss, line)) {
        if (line.empty()) continue;
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

static bool write_single_record_agd(
    const std::string& path,
    const std::string& read_id_with_suffix,
    float p_read,
    std::string& err) {
    agp::AgdHeader header{};
    header.magic = agp::AGD_MAGIC;
    header.version = agp::AGD_VERSION;
    header.num_records = 1;
    header.num_buckets = 8;
    header.d_max = 0.30f;
    header.lambda = 0.30f;
    header.library_type = 2;
    header.damage_validated = 1;
    header.damage_artifact = 0;
    header.channel_b_valid = 1;
    header.stop_decay_llr = 1000.0f;
    header.terminal_shift = 0.05f;

    std::vector<agp::AgdBucket> buckets(header.num_buckets);
    for (auto& b : buckets) b.record_offset = 0xFFFFFFFFu;

    agp::AgdRecord rec{};
    rec.id_hash = agp::fnv1a_hash(agp::strip_agp_suffix(read_id_with_suffix));
    rec.seq_len = 120;
    rec.frame_strand = agp::encode_frame_strand(0, false);
    rec.damage_pct_q = agp::quantize_damage_pct(10.0f);
    rec.p_damaged_q = agp::quantize_probability(p_read);
    rec.n_5prime = 0;
    rec.n_3prime = 0;
    std::memset(rec.codons_5prime, 255, sizeof(rec.codons_5prime));
    std::memset(rec.codons_3prime, 255, sizeof(rec.codons_3prime));

    const uint64_t bucket_idx = rec.id_hash % header.num_buckets;
    buckets[bucket_idx].record_offset = 0;

    const uint32_t chain = 0xFFFFFFFFu;

    std::ofstream out(path, std::ios::binary);
    if (!out.is_open()) {
        err = "cannot open " + path;
        return false;
    }
    out.write(reinterpret_cast<const char*>(&header), sizeof(header));
    out.write(reinterpret_cast<const char*>(buckets.data()),
              static_cast<std::streamsize>(buckets.size() * sizeof(agp::AgdBucket)));
    out.write(reinterpret_cast<const char*>(&rec), sizeof(rec));
    out.write(reinterpret_cast<const char*>(&chain), sizeof(chain));
    if (!out) {
        err = "failed writing AGD file";
        return false;
    }
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

static bool approx_zero(double x, double eps = 1e-6) {
    return std::fabs(x) <= eps;
}

static bool expect(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "FAILED: " << msg << "\n";
        return false;
    }
    return true;
}

static bool test_em_target_gamma_consistency(const std::string& agp_bin, const std::string& tmpdir) {
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
    if (!convert_hits_to_emi(agp_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string out;
    std::string cmd = shell_quote(agp_bin) + " damage-annotate --emi " + shell_quote(emi_path) + " --em";
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

static bool test_max_dist_recomputes_metrics(const std::string& agp_bin, const std::string& tmpdir) {
    std::cout << "Running --max-dist metric recompute regression...\n";
    const std::string hits_path = tmpdir + "/maxdist_hits.tsv";
    const std::string taln = "AAAAAAAAAARAAAAAAAAA";
    const std::string qaln = "AAAAAAAAAAWAAAAAAAAA";
    std::ostringstream hits;
    hits << "readM_+_0\ttargetM\t0.90\t20\t1\t0\t1\t20\t1\t20\t1e-30\t120\t20\t100\t"
         << qaln << "\t" << taln << "\n";
    std::string err;
    if (!write_file(hits_path, hits.str(), err)) {
        std::cerr << err << "\n";
        return false;
    }
    const std::string emi_path = tmpdir + "/maxdist_hits.emi2";
    if (!convert_hits_to_emi(agp_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string out1, out2;
    std::string cmd1 = shell_quote(agp_bin) + " damage-annotate --emi " + shell_quote(emi_path) +
                       " --d-max 0.3 --lambda 0.3";
    if (!run_command_capture(cmd1, out1, err)) {
        std::cerr << err << "\n";
        return false;
    }
    std::string cmd2 = shell_quote(agp_bin) + " damage-annotate --emi " + shell_quote(emi_path) +
                       " --d-max 0.3 --lambda 0.3 --max-dist 0";
    if (!run_command_capture(cmd2, out2, err)) {
        std::cerr << err << "\n";
        return false;
    }

    Table t1, t2;
    if (!parse_table(out1, t1, err) || !parse_table(out2, t2, err)) {
        std::cerr << err << "\n";
        return false;
    }
    std::unordered_map<std::string, std::string> r1, r2;
    if (!get_row_for_query(t1, "readM_+_0", r1) || !get_row_for_query(t2, "readM_+_0", r2)) {
        std::cerr << "readM_+_0 missing in output\n";
        return false;
    }

    if (!expect(r1["damage_consistent"] == "1", "unfiltered run should retain one damage site")) return false;
    if (!expect(r2["damage_consistent"] == "0", "max-dist run should remove distant site")) return false;

    double p1 = std::stod(r1["p_damaged"]);
    double p2 = std::stod(r2["p_damaged"]);
    double maxp2 = std::stod(r2["max_p_damage"]);
    if (!expect(p1 > 0.0, "unfiltered p_damaged should be > 0")) return false;
    if (!expect(approx_zero(p2), "filtered p_damaged should be recomputed to 0")) return false;
    if (!expect(approx_zero(maxp2), "filtered max_p_damage should be recomputed to 0")) return false;
    return true;
}

static bool test_threshold_is_damaged(const std::string& agp_bin, const std::string& tmpdir) {
    std::cout << "Running threshold/is_damaged regression...\n";
    const std::string hits_path = tmpdir + "/threshold_hits.tsv";
    std::ostringstream hits;
    // Use a damage-consistent substitution so posterior falls in a threshold-sensitive range.
    hits << "readT_+_0\ttargetT\t0.90\t4\t1\t0\t1\t4\t1\t4\t1e-10\t40\t4\t100\tWCDE\tRCDE\n";
    std::string err;
    if (!write_file(hits_path, hits.str(), err)) {
        std::cerr << err << "\n";
        return false;
    }
    const std::string emi_path = tmpdir + "/threshold_hits.emi2";
    if (!convert_hits_to_emi(agp_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    const std::string agd_path = tmpdir + "/threshold.agd";
    if (!write_single_record_agd(agd_path, "readT_+_0", 0.6f, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string out_lo, out_hi;
    std::string cmd_hi = shell_quote(agp_bin) + " damage-annotate --emi " + shell_quote(emi_path) +
                         " --damage-index " + shell_quote(agd_path) +
                         " --d-max 0.3 --no-identity --threshold 1.0";
    if (!run_command_capture(cmd_hi, out_hi, err)) {
        std::cerr << err << "\n";
        return false;
    }
    std::string cmd_lo = shell_quote(agp_bin) + " damage-annotate --emi " + shell_quote(emi_path) +
                         " --damage-index " + shell_quote(agd_path) +
                         " --d-max 0.3 --no-identity --threshold 0.0";
    if (!run_command_capture(cmd_lo, out_lo, err)) {
        std::cerr << err << "\n";
        return false;
    }

    Table thi, tlo;
    if (!parse_table(out_hi, thi, err) || !parse_table(out_lo, tlo, err)) {
        std::cerr << err << "\n";
        return false;
    }
    std::unordered_map<std::string, std::string> rhi, rlo;
    if (!get_row_for_query(thi, "readT_+_0", rhi) || !get_row_for_query(tlo, "readT_+_0", rlo)) {
        std::cerr << "readT_+_0 missing in output\n";
        return false;
    }

    if (!expect(rhi["posterior"] != "NA", "posterior should be informative with AGD p_read")) return false;
    if (!expect(rhi["damage_informative"] == "1", "damage_informative should be 1 with AGD evidence")) return false;

    const double posterior = std::stod(rhi["posterior"]);
    if (!expect(posterior > 0.0 && posterior < 1.0, "posterior should be strictly between 0 and 1")) return false;

    if (!expect(rhi["is_damaged"] == "0", "threshold 1.0 should mark read as not damaged")) return false;
    if (!expect(rlo["is_damaged"] == "1", "threshold 0.0 should mark read as damaged")) return false;
    return true;
}

static bool test_protein_mapping_filters(const std::string& agp_bin, const std::string& tmpdir) {
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
    if (!convert_hits_to_emi(agp_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    const std::string gene_summary = tmpdir + "/gene_summary.tsv";
    const std::string protein_summary = tmpdir + "/protein_summary.tsv";
    const std::string protein_filtered = tmpdir + "/protein_filtered.tsv";
    const std::string combined_output = tmpdir + "/combined_output.tsv";
    std::string out;
    std::string cmd = shell_quote(agp_bin) + " damage-annotate --emi " + shell_quote(emi_path) +
                      " --gene-summary " + shell_quote(gene_summary) +
                      " --protein-summary " + shell_quote(protein_summary) +
                      " --protein-filtered " + shell_quote(protein_filtered) +
                      " --combined-output " + shell_quote(combined_output) +
                      " --min-reads 1 --min-breadth 0 --min-depth 0 " +
                      " --min-positional-score 0.4 --min-terminal-ratio 0.5";
    if (!run_command_capture(cmd, out, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string gene_txt, protein_txt, protein_filt_txt, combined_txt;
    if (!read_file(gene_summary, gene_txt, err) ||
        !read_file(protein_summary, protein_txt, err) ||
        !read_file(protein_filtered, protein_filt_txt, err) ||
        !read_file(combined_output, combined_txt, err)) {
        std::cerr << err << "\n";
        return false;
    }

    Table tg, tp, tpf, tc;
    if (!parse_table(gene_txt, tg, err) || !parse_table(protein_txt, tp, err) ||
        !parse_table(protein_filt_txt, tpf, err) || !parse_table(combined_txt, tc, err)) {
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

    if (!expect(tpf.rows.size() == 1, "protein-filtered should contain exactly one passing protein")) return false;
    auto only = row_to_map(tpf, tpf.rows.front());
    if (!expect(only["protein_id"] == "prot_good", "protein-filtered should contain prot_good")) return false;

    bool combined_good_pass = false, combined_bad_fail = false;
    for (const auto& row : tc.rows) {
        auto m = row_to_map(tc, row);
        if (m["target_id"] == "prot_good" && m["pass_mapping_filter"] == "1") combined_good_pass = true;
        if (m["target_id"] == "prot_bad" && m["pass_mapping_filter"] == "0") combined_bad_fail = true;
    }
    if (!expect(combined_good_pass, "combined output should mark prot_good as pass_mapping_filter=1")) return false;
    if (!expect(combined_bad_fail, "combined output should mark prot_bad as pass_mapping_filter=0")) return false;
    return true;
}

static bool test_blast8_unique_export(const std::string& agp_bin, const std::string& tmpdir) {
    std::cout << "Running BLAST8 unique export regression...\n";
    const std::string hits_path = tmpdir + "/blast8_unique_hits.tsv";
    std::ostringstream hits;
    hits << "readU_+_0\ttargetU\t0.95\t6\t0\t0\t1\t6\t1\t6\t1e-20\t80\t6\t120\tACDEFG\tACDEFG\n";
    hits << "readM_+_0\ttargetA\t0.95\t6\t0\t0\t1\t6\t1\t6\t1e-20\t80\t6\t120\tACDEFG\tACDEFG\n";
    hits << "readM_+_0\ttargetB\t0.90\t6\t1\t0\t1\t6\t1\t6\t1e-10\t60\t6\t120\tACDEFG\tACDEFA\n";

    std::string err;
    if (!write_file(hits_path, hits.str(), err)) {
        std::cerr << err << "\n";
        return false;
    }
    const std::string emi_path = tmpdir + "/blast8_unique_hits.emi2";
    if (!convert_hits_to_emi(agp_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    const std::string blast8_unique = tmpdir + "/blast8_unique.tsv";
    std::string out;
    std::string cmd = shell_quote(agp_bin) + " damage-annotate --emi " + shell_quote(emi_path) +
                      " --blast8-unique " + shell_quote(blast8_unique);
    if (!run_command_capture(cmd, out, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string blast_txt;
    if (!read_file(blast8_unique, blast_txt, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::istringstream iss(blast_txt);
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(iss, line)) {
        if (!line.empty()) lines.push_back(line);
    }

    if (!expect(lines.size() == 1, "blast8-unique should contain only one uniquely-mapped read")) return false;
    if (!expect(lines[0].find("readU_+_0\ttargetU\t") == 0, "blast8-unique should retain readU->targetU")) return false;
    if (!expect(lines[0].find("readM_+_0") == std::string::npos, "blast8-unique should exclude multi-mapped readM")) return false;
    return true;
}

static bool test_analysis_prefix_bundle(const std::string& agp_bin, const std::string& tmpdir) {
    std::cout << "Running analysis-prefix bundle regression...\n";
    const std::string hits_path = tmpdir + "/bundle_hits.tsv";
    std::ostringstream hits;
    hits << "r1_+_0\tprotA\t0.95\t6\t0\t0\t1\t6\t1\t6\t1e-20\t80\t6\t120\tACDEFG\tACDEFG\n";
    hits << "r2_+_0\tprotA\t0.95\t6\t1\t0\t1\t6\t10\t15\t1e-10\t70\t6\t120\tACDEFA\tACDEFG\n";
    std::string err;
    if (!write_file(hits_path, hits.str(), err)) {
        std::cerr << err << "\n";
        return false;
    }
    const std::string emi_path = tmpdir + "/bundle_hits.emi2";
    if (!convert_hits_to_emi(agp_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    const std::string map_path = tmpdir + "/bundle_map.tsv";
    if (!write_file(map_path, "protA\tK00001\n", err)) {
        std::cerr << err << "\n";
        return false;
    }

    const std::string prefix = tmpdir + "/bundle";
    std::string out;
    std::string cmd = shell_quote(agp_bin) + " damage-annotate --emi " + shell_quote(emi_path) +
                      " --map " + shell_quote(map_path) +
                      " --analysis-prefix " + shell_quote(prefix) +
                      " --min-reads 1 --min-breadth 0 --min-depth 0";
    if (!run_command_capture(cmd, out, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string reads_txt, proteins_txt, blast8_txt, categories_txt;
    if (!read_file(prefix + ".reads.tsv", reads_txt, err) ||
        !read_file(prefix + ".protein.tsv", proteins_txt, err) ||
        !read_file(prefix + ".blast8.unique.tsv", blast8_txt, err) ||
        !read_file(prefix + ".categories.tsv", categories_txt, err)) {
        std::cerr << err << "\n";
        return false;
    }

    if (!expect(!reads_txt.empty(), "bundle reads output should not be empty")) return false;
    if (!expect(proteins_txt.find("pass_mapping_filter\tfilter_fail") != std::string::npos,
                "bundle proteins output should include pass_mapping_filter and filter_fail")) return false;
    if (!expect(categories_txt.find("n_damaged_genes\tdamaged_gene_frac\tdamage_enrichment") != std::string::npos,
                "bundle categories output should include damage enrichment columns")) return false;
    if (!expect(!blast8_txt.empty(), "bundle blast8 unique output should not be empty")) return false;
    return true;
}

static bool test_analysis_explicit_named_outputs(const std::string& agp_bin, const std::string& tmpdir) {
    std::cout << "Running explicit named analysis outputs regression...\n";
    const std::string hits_path = tmpdir + "/named_hits.tsv";
    std::ostringstream hits;
    hits << "n1_+_0\tprotN\t0.95\t6\t0\t0\t1\t6\t1\t6\t1e-20\t80\t6\t120\tACDEFG\tACDEFG\n";
    std::string err;
    if (!write_file(hits_path, hits.str(), err)) {
        std::cerr << err << "\n";
        return false;
    }
    const std::string emi_path = tmpdir + "/named_hits.emi2";
    if (!convert_hits_to_emi(agp_bin, hits_path, emi_path, err)) {
        std::cerr << err << "\n";
        return false;
    }

    const std::string map_path = tmpdir + "/named_map.tsv";
    if (!write_file(map_path, "protN\tK99999\n", err)) {
        std::cerr << err << "\n";
        return false;
    }

    const std::string proteins_out = tmpdir + "/my.proteins.custom.tsv";
    const std::string blast8_out = tmpdir + "/my.unique.custom.blast8";
    const std::string categories_out = tmpdir + "/my.categories.custom.tsv";
    std::string out;
    std::string cmd = shell_quote(agp_bin) + " damage-annotate --emi " + shell_quote(emi_path) +
                      " --map " + shell_quote(map_path) +
                      " --analysis-proteins " + shell_quote(proteins_out) +
                      " --analysis-blast8 " + shell_quote(blast8_out) +
                      " --analysis-categories " + shell_quote(categories_out) +
                      " --min-reads 1 --min-breadth 0 --min-depth 0";
    if (!run_command_capture(cmd, out, err)) {
        std::cerr << err << "\n";
        return false;
    }

    std::string ptxt, btxt, ctxt;
    if (!read_file(proteins_out, ptxt, err) ||
        !read_file(blast8_out, btxt, err) ||
        !read_file(categories_out, ctxt, err)) {
        std::cerr << err << "\n";
        return false;
    }
    if (!expect(!ptxt.empty(), "explicit proteins output should exist and be non-empty")) return false;
    if (!expect(!btxt.empty(), "explicit blast8 output should exist and be non-empty")) return false;
    if (!expect(!ctxt.empty(), "explicit categories output should exist and be non-empty")) return false;
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path-to-agp>\n";
        return 2;
    }
    const std::string agp_bin = argv[1];

    char tmp_template[] = "/tmp/agp_regress_XXXXXX";
    char* tmp = mkdtemp(tmp_template);
    if (!tmp) {
        std::cerr << "Failed to create temp dir\n";
        return 2;
    }
    const std::string tmpdir = tmp;

    bool ok = true;
    ok = test_em_target_gamma_consistency(agp_bin, tmpdir) && ok;
    ok = test_max_dist_recomputes_metrics(agp_bin, tmpdir) && ok;
    ok = test_threshold_is_damaged(agp_bin, tmpdir) && ok;
    ok = test_protein_mapping_filters(agp_bin, tmpdir) && ok;
    ok = test_blast8_unique_export(agp_bin, tmpdir) && ok;
    ok = test_analysis_prefix_bundle(agp_bin, tmpdir) && ok;
    ok = test_analysis_explicit_named_outputs(agp_bin, tmpdir) && ok;

    if (!ok) return 1;
    std::cout << "All damage-annotate regressions passed.\n";
    return 0;
}
