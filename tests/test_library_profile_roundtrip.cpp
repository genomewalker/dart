// Contract: the TDM1 library damage summary that fqdup writes is the summary DART applies.
//
// fqdup writes the sidecar with taph::write_damage_summary; `dart predict --library-profile`
// reads it with taph::read_damage_summary and stamps it into the .agd header that
// damage-annotate scores against. Both sides must build the same pinned libtaph.
//
//   1. write -> read through the pinned libtaph preserves every field;
//   2. `dart predict --library-profile` stamps those values into the .agd header;
//   3. if tests/data/tdm1/fqdup_sample.tdm1 exists (a sidecar written by a real fqdup run),
//      it parses, matches fqdup_sample.expected.tsv, and predict accepts it.

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "dart/damage_index_reader.hpp"
#include "taph/damage_summary.hpp"

namespace {

int failures = 0;

bool expect(bool cond, const std::string& what) {
    if (!cond) { std::cerr << "  FAIL: " << what << "\n"; ++failures; }
    return cond;
}

bool near(double a, double b, double tol = 1e-6) { return std::fabs(a - b) <= tol; }

// Every serialized field, by the name fqdup's expected.tsv uses.
const std::vector<std::pair<std::string, std::function<double(const taph::DamageSummary&)>>> kFields = {
    {"d_max_combined",        [](const taph::DamageSummary& s) { return s.d_max_combined; }},
    {"d_max_5prime",          [](const taph::DamageSummary& s) { return s.d_max_5prime; }},
    {"d_max_3prime",          [](const taph::DamageSummary& s) { return s.d_max_3prime; }},
    {"lambda_5prime",         [](const taph::DamageSummary& s) { return s.lambda_5prime; }},
    {"lambda_3prime",         [](const taph::DamageSummary& s) { return s.lambda_3prime; }},
    {"fit_baseline_5prime",   [](const taph::DamageSummary& s) { return s.fit_baseline_5prime; }},
    {"fit_baseline_3prime",   [](const taph::DamageSummary& s) { return s.fit_baseline_3prime; }},
    {"terminal_shift_5prime", [](const taph::DamageSummary& s) { return s.terminal_shift_5prime; }},
    {"stop_decay_llr_5prime", [](const taph::DamageSummary& s) { return s.stop_decay_llr_5prime; }},
    {"d_max_cooc_point",      [](const taph::DamageSummary& s) { return s.d_max_cooc_point; }},
    {"pi_point",              [](const taph::DamageSummary& s) { return s.pi_point; }},
    {"pi_shape_lambda",       [](const taph::DamageSummary& s) { return s.pi_shape_lambda; }},
    {"pi_shape_baseline",     [](const taph::DamageSummary& s) { return s.pi_shape_baseline; }},
    {"pi_shape_3p_lambda",    [](const taph::DamageSummary& s) { return s.pi_shape_3p_lambda; }},
    {"pi_shape_3p_baseline",  [](const taph::DamageSummary& s) { return s.pi_shape_3p_baseline; }},
    {"pi_damaged_point",      [](const taph::DamageSummary& s) { return s.pi_damaged_point; }},
    {"library_type",          [](const taph::DamageSummary& s) { return s.library_type; }},
    {"damage_validated",      [](const taph::DamageSummary& s) { return s.damage_validated; }},
    {"damage_artifact",       [](const taph::DamageSummary& s) { return s.damage_artifact; }},
    {"channel_b_valid",       [](const taph::DamageSummary& s) { return s.channel_b_valid; }},
    {"pi_shape_fitted",       [](const taph::DamageSummary& s) { return s.pi_shape_fitted; }},
    {"pi_shape_3p_fitted",    [](const taph::DamageSummary& s) { return s.pi_shape_3p_fitted; }},
    {"pi_state",              [](const taph::DamageSummary& s) { return s.pi_state; }},
};

std::string make_tmpdir() {
    const char* base = std::getenv("TMPDIR");
    std::string tmpl = std::string(base && *base ? base : "/tmp") + "/dart_tdm1_XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    const char* dir = mkdtemp(buf.data());
    return dir ? std::string(dir) : std::string();
}

// Deterministic random reads: predict only needs a population to run on; the damage
// parameters under test come from the library profile, not from these reads.
void write_reads(const std::string& path) {
    std::ofstream fq(path);
    uint64_t x = 0x9E3779B97F4A7C15ULL;
    const char nt[] = "ACGT";
    for (int r = 0; r < 400; ++r) {
        std::string seq(90, 'A');
        for (char& c : seq) { x ^= x << 13; x ^= x >> 7; x ^= x << 17; c = nt[x & 3]; }
        fq << "@r" << r << "\n" << seq << "\n+\n" << std::string(seq.size(), 'I') << "\n";
    }
}

// Run predict with the sidecar and check the .agd header carries it.
void check_predict_applies(const std::string& dart, const std::string& dir,
                           const std::string& tdm1, const taph::DamageSummary& s,
                           const std::string& label) {
    const std::string reads = dir + "/reads.fq";
    const std::string agd = dir + "/" + label + ".agd";
    const std::string cmd = dart + " predict -i " + reads + " -o " + dir + "/" + label + ".gff" +
                            " --damage-index " + agd + " --library-profile " + tdm1 +
                            " --threads 2 > " + dir + "/" + label + ".log 2>&1";
    if (!expect(std::system(cmd.c_str()) == 0, label + ": dart predict failed (see " + dir + ")")) return;

    dart::DamageIndexReader rd(agd, false);
    if (!expect(rd.is_valid(), label + ": .agd unreadable")) return;
    expect(near(rd.d_max(), s.d_max_combined), label + ": agd d_max != d_max_combined");
    expect(near(rd.lambda(), (s.lambda_5prime + s.lambda_3prime) / 2.0f), label + ": agd lambda");
    expect(rd.library_type() == s.library_type, label + ": agd library_type");
    expect(rd.pi_state() == s.pi_state, label + ": agd pi_state");
    expect(near(rd.pi_point(), s.pi_point), label + ": agd pi_point");
    if (s.pi_damaged_point >= 0.0f)
        expect(near(rd.damage_prior_point(), s.pi_damaged_point), label + ": agd damage_prior_point");
}

void check_synthetic(const std::string& dart, const std::string& dir) {
    std::cout << "Running TDM1 write/read + predict --library-profile round trip...\n";
    taph::DamageSummary s;
    s.d_max_combined = 0.271f;   s.d_max_5prime = 0.268f;     s.d_max_3prime = 0.143f;
    s.lambda_5prime = 0.41f;     s.lambda_3prime = 0.37f;
    s.fit_baseline_5prime = 0.012f; s.fit_baseline_3prime = 0.009f;
    s.terminal_shift_5prime = 0.118f; s.stop_decay_llr_5prime = 42.5f;
    s.d_max_cooc_point = 0.255;  s.pi_point = 0.31;
    s.pi_shape_lambda = 0.44;    s.pi_shape_baseline = 0.011;
    s.pi_shape_3p_lambda = 0.39; s.pi_shape_3p_baseline = 0.008;
    s.pi_damaged_point = 0.42f;
    s.library_type = 2; s.damage_validated = 1; s.damage_artifact = 0; s.channel_b_valid = 1;
    s.pi_shape_fitted = 1; s.pi_shape_3p_fitted = 1; s.pi_state = 0;  // DETECTED

    const std::string path = dir + "/synthetic.tdm1";
    if (!expect(taph::write_damage_summary(path, s), "write_damage_summary")) return;
    auto r = taph::read_damage_summary(path);
    if (!expect(r.has_value(), "read_damage_summary of a just-written file")) return;
    for (const auto& [name, get] : kFields)
        expect(get(*r) == get(s), "round trip changed " + name);

    check_predict_applies(dart, dir, path, s, "synthetic");
}

void check_fqdup_fixture(const std::string& dart, const std::string& dir, const std::string& data) {
    const std::string tdm1 = data + "/tdm1/fqdup_sample.tdm1";
    const std::string exp = data + "/tdm1/fqdup_sample.expected.tsv";
    struct stat st{};
    if (stat(tdm1.c_str(), &st) != 0) {
        std::cout << "  no fqdup fixture at " << tdm1 << "; skipping the producer half\n";
        return;
    }
    std::cout << "Running fqdup-written TDM1 fixture...\n";
    auto s = taph::read_damage_summary(tdm1);
    if (!expect(s.has_value(), "fqdup fixture does not parse with the pinned libtaph")) return;

    std::map<std::string, double> want;
    std::ifstream in(exp);
    if (!expect(static_cast<bool>(in), "missing " + exp)) return;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        std::string k; double v;
        if (ls >> k >> v) want[k] = v;
    }
    expect(!want.empty(), "fixture expected.tsv has no fields");
    for (const auto& [name, get] : kFields) {
        auto it = want.find(name);
        if (it == want.end()) continue;
        expect(near(get(*s), it->second, 1e-5), "fqdup fixture " + name);
        want.erase(it);
    }
    for (const auto& kv : want)
        if (kv.first != "version") expect(false, "expected.tsv names an unknown field: " + kv.first);

    check_predict_applies(dart, dir, tdm1, *s, "fqdup");
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <path-to-dart> <tests/data dir>\n";
        return 2;
    }
    const std::string dir = make_tmpdir();
    if (dir.empty()) { std::cerr << "mkdtemp failed\n"; return 2; }
    write_reads(dir + "/reads.fq");

    check_synthetic(argv[1], dir);
    check_fqdup_fixture(argv[1], dir, argv[2]);

    if (failures) { std::cerr << failures << " failure(s); outputs kept in " << dir << "\n"; return 1; }
    std::filesystem::remove_all(dir);
    std::cout << "TDM1 library-profile contract holds.\n";
    return 0;
}
