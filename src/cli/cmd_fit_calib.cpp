// dart fit-calib — fit the coherent-LR calibration model (f_true + f_0) from a
// LABELED synthetic search.
//
// Inputs:
//   --hits   mmseqs m8 (query target pident alnlen ... evalue bits)
//   --truth  the ground-truth TSV written by `dart simulate --source-faa`
//            (read_id cohort source_protein source_aa_len frag_nt_start
//             frag_nt_end strand is_damaged ct_positions ga_positions)
//
// f_true is fit from TRUE-SOURCE hits (ancient+modern reads whose alignment target
// equals the read's source protein) — the correct labeled distribution, NOT the
// winner-truncated observed target-hit distribution. f_0 is fit from the per-query
// BEST bit of SPURIOUS reads (sourced from held-out proteins that are absent from
// the search DB), which is exactly the null/decoy score distribution. Both are
// binned by aligned length so they describe the same statistic on the same bins.

#include "subcommand.hpp"
#include "dart/em_reassign.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace dart {
namespace cli {
namespace {

struct TruthRow { std::string cohort; std::string source; };

// dart predict emits ORF headers as "<read_id>_<strand>_<frame>" (strand ∈ {+,-},
// frame a small integer). The production AA search is keyed on those ORF ids, so
// map an m8 query back to its read_id: strip a trailing "_<+|->_<digits>" when
// present, else return the query unchanged (read-keyed searches pass through).
static std::string base_read_id(const std::string& q) {
    const size_t p2 = q.rfind('_');
    if (p2 == std::string::npos || p2 == 0 || p2 + 1 >= q.size()) return q;
    for (size_t i = p2 + 1; i < q.size(); ++i)
        if (!std::isdigit(static_cast<unsigned char>(q[i]))) return q;  // frame not numeric
    const size_t p1 = q.rfind('_', p2 - 1);
    if (p1 == std::string::npos) return q;
    const std::string strand = q.substr(p1 + 1, p2 - p1 - 1);
    if (strand != "+" && strand != "-") return q;
    return q.substr(0, p1);
}

static std::unordered_map<std::string, TruthRow> load_truth(const std::string& path) {
    std::unordered_map<std::string, TruthRow> t;
    std::ifstream f(path);
    if (!f) { std::cerr << "fit-calib: cannot open --truth " << path << "\n"; return t; }
    std::string line;
    bool first = true;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (first) { first = false; if (line.rfind("read_id", 0) == 0) continue; }
        std::istringstream ss(line);
        std::string rid, cohort, source;
        if (!std::getline(ss, rid, '\t')) continue;
        if (!std::getline(ss, cohort, '\t')) continue;
        if (!std::getline(ss, source, '\t')) continue;
        t.emplace(rid, TruthRow{cohort, source});
    }
    return t;
}

static void print_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " fit-calib --hits HITS.m8 --truth TRUTH.tsv --out CALIB.model\n\n"
        "Fit the coherent-LR calibration model (f_true + f_0) from a labeled synthetic\n"
        "search. f_true = Gaussian(bit) of true-source hits per length bin; f_0 = Gumbel\n"
        "of per-query best spurious bits per length bin.\n\n"
        "  --hits FILE      mmseqs m8 (query target pident alnlen ... evalue bits)\n"
        "  --truth FILE     ground-truth TSV from `dart simulate --source-faa`\n"
        "  --out FILE       output calibration model\n"
        "  --min-true N     min true hits/bin before pooled fallback (default 100)\n"
        "  --min-null N     min spurious queries/bin before pooled fallback (default 200)\n"
        "  --bits-col N     1-based m8 column holding the bit score (default 12)\n"
        "  --alnlen-col N   1-based m8 column holding the aligned length (default 4)\n"
        "  --nu FLOAT       Student-t dof for the v2 length-continuous f_true (default 6; <=0 → Gaussian)\n"
        "  --legacy-bins    keep only the v1 discrete-bin model (no length-continuous fit)\n";
}

} // namespace

int cmd_fit_calib(int argc, char* argv[]) {
    std::string hits_file, truth_file, out_file;
    uint32_t min_true = 100, min_null = 200;
    int bits_col = 12, alnlen_col = 4;   // mmseqs default m8 layout (1-based)
    double nu = 6.0;          // Student-t dof for f_true (<=0 → Gaussian)
    bool legacy_bins = false; // keep only the v1 discrete-bin model

    for (int i = 1; i < argc; ++i) {
        auto need = [&](const char* f) { if (i + 1 >= argc) { std::cerr << "fit-calib: " << f << " needs an argument\n"; return false; } return true; };
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(argv[0]); return 0; }
        else if (!strcmp(argv[i], "--hits"))       { if (!need(argv[i])) return 1; hits_file = argv[++i]; }
        else if (!strcmp(argv[i], "--truth"))      { if (!need(argv[i])) return 1; truth_file = argv[++i]; }
        else if (!strcmp(argv[i], "--out"))        { if (!need(argv[i])) return 1; out_file = argv[++i]; }
        else if (!strcmp(argv[i], "--min-true"))   { if (!need(argv[i])) return 1; min_true = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (!strcmp(argv[i], "--min-null"))   { if (!need(argv[i])) return 1; min_null = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (!strcmp(argv[i], "--bits-col"))   { if (!need(argv[i])) return 1; bits_col = std::stoi(argv[++i]); }
        else if (!strcmp(argv[i], "--alnlen-col")) { if (!need(argv[i])) return 1; alnlen_col = std::stoi(argv[++i]); }
        else if (!strcmp(argv[i], "--nu"))          { if (!need(argv[i])) return 1; nu = std::stod(argv[++i]); }
        else if (!strcmp(argv[i], "--legacy-bins")) { legacy_bins = true; }
        else { std::cerr << "fit-calib: unknown option '" << argv[i] << "'\n"; print_usage(argv[0]); return 1; }
    }
    if (hits_file.empty() || truth_file.empty() || out_file.empty()) {
        print_usage(argv[0]); return 1;
    }

    auto truth = load_truth(truth_file);
    if (truth.empty()) { std::cerr << "fit-calib: no truth rows\n"; return 1; }

    std::ifstream hf(hits_file);
    if (!hf) { std::cerr << "fit-calib: cannot open --hits " << hits_file << "\n"; return 1; }

    // f_true: per-bin true-source hit bits. f_0: per-query best spurious (bit,bin).
    std::vector<std::vector<double>> true_bits(kNumNullBins);
    std::unordered_map<std::string, std::pair<double, uint16_t>> null_best;  // query → (best bit, alnlen)
    // Raw (aln_len, bit) pairs for the v2 length-continuous fit (sol short-fragment fix).
    std::vector<std::pair<double, double>> true_LB;
    double censor = std::numeric_limits<double>::infinity();
    size_t n_true = 0, n_null_rows = 0;

    std::string line;
    std::vector<std::string> col;
    while (std::getline(hf, line)) {
        if (line.empty()) continue;
        col.clear();
        std::istringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, '\t')) col.push_back(tok);
        if (static_cast<int>(col.size()) < std::max(bits_col, alnlen_col)) continue;

        const std::string& query = col[0];
        const std::string& target = col[1];
        // Exact read-keyed match first; else fold the ORF id back to its read_id.
        auto it = truth.find(query);
        std::string read_key = query;
        if (it == truth.end()) {
            read_key = base_read_id(query);
            it = truth.find(read_key);
            if (it == truth.end()) continue;
        }
        const TruthRow& tr = it->second;

        double bits = 0.0; uint16_t alnlen = 0;
        try {
            bits = std::stod(col[static_cast<size_t>(bits_col - 1)]);
            alnlen = static_cast<uint16_t>(std::stoul(col[static_cast<size_t>(alnlen_col - 1)]));
        } catch (...) { continue; }
        int b = null_length_bin(alnlen);

        if ((tr.cohort == "ancient" || tr.cohort == "modern") && target == tr.source) {
            true_bits[static_cast<size_t>(b)].push_back(bits);
            true_LB.emplace_back(static_cast<double>(alnlen), bits);
            censor = std::min(censor, bits);
            ++n_true;
        } else if (tr.cohort == "spurious") {
            // Null is per-READ (the EM's null column takes the max over a read's
            // alignments), so track the best spurious bit per read_id, not per ORF.
            auto& best = null_best[read_key];
            if (best.second == 0 || bits > best.first) best = {bits, alnlen};
            ++n_null_rows;
        }
    }

    if (n_true == 0) { std::cerr << "fit-calib: no labeled true-source hits found\n"; return 1; }

    std::vector<std::vector<double>> null_bits(kNumNullBins);
    std::vector<std::pair<double, double>> null_LB;
    for (const auto& kv : null_best) {
        int b = null_length_bin(kv.second.second);
        null_bits[static_cast<size_t>(b)].push_back(kv.second.first);
        null_LB.emplace_back(static_cast<double>(kv.second.second), kv.second.first);
    }

    CalibModel m;
    // v1 discrete sections written for provenance regardless of the active path.
    m.true_bins = fit_true_model(true_bits, min_true);
    m.null = fit_null_model(null_bits, min_null);
    m.censor = std::isfinite(censor) ? static_cast<float>(censor) : 0.0f;

    // v2 length-continuous fit is the default active model (sol short-fragment fix);
    // --legacy-bins keeps the v1 discrete-bin path for the regression baseline.
    if (!legacy_bins) {
        if (!fit_calib_continuous(true_LB, null_LB, nu, m)) {
            std::cerr << "fit-calib: too few true hits for the length-continuous fit "
                         "(need >=10); falling back to v1 discrete bins\n";
        }
    }

    if (!write_calib_model(out_file, m)) {
        std::cerr << "fit-calib: cannot write --out " << out_file << "\n"; return 1;
    }

    size_t n_null_q = null_best.size();
    std::cerr << "fit-calib: f_true from " << n_true << " true-source hits, f_0 from "
              << n_null_q << " spurious queries (" << n_null_rows << " rows), censor="
              << m.censor << " -> " << out_file << "\n";
    if (m.cont) {
        std::cerr << "  [v2 length-continuous] f_true mu(L)=" << m.ft_mu_a << "+" << m.ft_mu_b
                  << "*L  sigma(L)=exp(" << m.ft_ls_a << "+" << m.ft_ls_b << "*L)  nu=" << m.ft_nu << "\n"
                  << "                        f_0 K-A loc(L)=" << m.f0_loc_a << "+" << m.f0_loc_b
                  << "*L  beta=" << m.f0_beta << "\n";
        // Print the implied mu_true / loc_null at a few lengths to sanity-check separation.
        for (int L : {10, 13, 17, 25, 40}) {
            const double mut = m.ft_mu_a + m.ft_mu_b * L;
            const double loc = m.f0_loc_a + m.f0_loc_b * L;
            std::cerr << "                        L=" << L << ": mu_true=" << mut
                      << "  loc_null=" << loc << "  gap=" << (mut - loc) << "\n";
        }
    }
    for (int b = 0; b < kNumNullBins; ++b) {
        const auto& t = m.true_bins[static_cast<size_t>(b)];
        std::cerr << "  bin " << b << " [" << kNullBinEdges[b] << "," << kNullBinEdges[b + 1] << ")"
                  << " f_true N=" << t.n << " mu=" << t.mu0 << " sd=" << t.sd0
                  << " | f_0 mu=" << m.null.mu[static_cast<size_t>(b)]
                  << " beta=" << m.null.beta[static_cast<size_t>(b)] << "\n";
    }
    return 0;
}

namespace {
    struct FitCalibRegistrar {
        FitCalibRegistrar() {
            SubcommandRegistry::instance().register_command(
                "fit-calib",
                "Fit the coherent-LR calibration model (f_true + f_0) from labeled hits",
                cmd_fit_calib, 31);
        }
    } fit_calib_registrar;
}

}  // namespace cli
}  // namespace dart
