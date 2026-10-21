// Contract: the decoy null model path, end to end on the tutorial data.
//
//   1. .null sidecar round trip: write -> read preserves every field, rewrite is byte-identical,
//      a bad magic, a truncated file and an unknown version are refused with a message;
//   2. `dart predict --fasta-aa X` writes the reversed proteins of every 10th input read to the
//      decoy file with the DECOY_ prefix, and nothing else; --decoy-every 0 writes no file;
//   3. `dart hits2emi` with and without --decoy-hits writes byte-identical EMIs, the sidecar
//      exists only with the flag, and `dart null` on the same decoy hits fits the same bins;
//   4. `dart damage-annotate` finds <emi>.null by itself and gives the same output as
//      --em-null-model; without a sidecar and without --no-null it refuses.
//
// Usage: dart-null-model-tests <dart binary> <examples/tutorial dir>

#include "dart/null_model.hpp"

#include <zlib.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

int failures = 0;

bool expect(bool cond, const std::string& what) {
    if (!cond) { std::cerr << "  FAIL: " << what << "\n"; ++failures; }
    return cond;
}

std::string make_tmpdir() {
    const char* base = std::getenv("TMPDIR");
    std::string tmpl = std::string(base && *base ? base : "/tmp") + "/dart_null_XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    if (!mkdtemp(buf.data())) { std::cerr << "mkdtemp failed\n"; std::exit(2); }
    return std::string(buf.data());
}

int run(const std::string& cmd, const std::string& log) {
    const std::string full = cmd + " > " + log + " 2>&1";
    const int rc = std::system(full.c_str());
    return WIFEXITED(rc) ? WEXITSTATUS(rc) : -1;
}

std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

bool exists(const std::string& path) { return std::ifstream(path).good(); }

std::vector<std::string> gz_lines(const std::string& path) {
    std::vector<std::string> out;
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) return out;
    std::vector<char> buf(1 << 20);
    std::string cur;
    int n;
    while ((n = gzread(f, buf.data(), static_cast<unsigned>(buf.size()))) > 0) {
        for (int i = 0; i < n; ++i) {
            if (buf[i] == '\n') { out.push_back(cur); cur.clear(); }
            else cur.push_back(buf[i]);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    gzclose(f);
    return out;
}

// FASTA: id (first token after '>') -> sequence, in file order.
std::vector<std::pair<std::string, std::string>> read_fasta(const std::string& path) {
    std::vector<std::pair<std::string, std::string>> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            const size_t sp = line.find(' ');
            out.emplace_back(line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1), "");
        } else if (!out.empty()) {
            out.back().second += line;
        }
    }
    return out;
}

// Protein ids are <read id>_<strand>_<frame>_<n>; strip the three trailing fields.
std::string read_of(const std::string& protein_id) {
    size_t p = protein_id.size();
    for (int k = 0; k < 3; ++k) { p = protein_id.rfind('_', p - 1); if (p == std::string::npos) return protein_id; }
    return protein_id.substr(0, p);
}

void test_roundtrip(const std::string& dir) {
    std::cout << "round trip... ";
    std::vector<uint32_t> n = {250, 300, 400, 10, 500, 600};
    std::vector<double> sum(6), sq(6);
    for (int b = 0; b < 6; ++b) { const double mean = 20.0 + 3.0 * b, sd = 2.0 + 0.25 * b; sum[b] = n[b] * mean; sq[b] = n[b] * (sd * sd + mean * mean); }
    dart::NullModel m = dart::fit_null_model_moments(n, sum, sq, 200);
    expect(m.mu.size() == 6 && m.n == n, "fit keeps bins and counts");
    expect(std::fabs(m.beta[0] - std::sqrt(6.0 * 4.0 * n[0] / (n[0] - 1.0)) / M_PI) < 1e-4, "beta from the n-1 sample sd");
    uint32_t pn = 0; double ps = 0, pss = 0;
    for (int b = 0; b < 6; ++b) { pn += n[b]; ps += sum[b]; pss += sq[b]; }
    const dart::NullModel pooled = dart::fit_null_model_moments({pn}, {ps}, {pss}, 200);
    expect(m.mu[3] == pooled.mu[0] && m.beta[3] == pooled.beta[0] && m.mu[3] != m.mu[0], "sparse bin n=10 takes the pooled fit");
    dart::NullProvenance prov;
    prov.dart_version = "test"; prov.created_utc = 1700000000; prov.n_decoy_queries = 2060; prov.n_decoy_hits = 9999;
    prov.db_label = "tutorial"; prov.inputs.push_back({"a.tsv", 123, 456}); prov.inputs.push_back({"b.tsv", 7, 8});
    const std::string path = dir + "/rt.null";
    std::string err;
    expect(dart::write_null_model(path, m, prov, &err), "write: " + err);
    dart::NullModel r; dart::NullProvenance rp;
    expect(dart::read_null_model(path, r, &rp, &err), "read: " + err);
    expect(r.mu == m.mu && r.beta == m.beta && r.n == m.n && r.sum == m.sum && r.sum_sq == m.sum_sq, "fields survive the round trip");
    expect(!rp.legacy_tsv && rp.dart_version == "test" && rp.created_utc == 1700000000 && rp.n_decoy_queries == 2060 &&
           rp.n_decoy_hits == 9999 && rp.db_label == "tutorial" && rp.inputs.size() == 2 && rp.inputs[1].path == "b.tsv" &&
           rp.inputs[0].size == 123 && rp.inputs[1].mtime == 8, "provenance survives the round trip");
    expect(dart::write_null_model(path + ".2", r, rp, &err) && slurp(path) == slurp(path + ".2"), "rewrite is byte-identical");
    std::string bytes = slurp(path);
    { std::string bad = bytes; bad[0] = 'X'; std::ofstream(dir + "/bad.null", std::ios::binary) << bad;
      dart::NullModel x; expect(!dart::read_null_model(dir + "/bad.null", x, nullptr, &err) && !err.empty(), "bad magic refused: " + err); }
    { std::ofstream(dir + "/trunc.null", std::ios::binary) << bytes.substr(0, bytes.size() - 9);
      dart::NullModel x; expect(!dart::read_null_model(dir + "/trunc.null", x, nullptr, &err) && err.find("truncated") != std::string::npos, "truncated refused: " + err); }
    { std::string v2 = bytes; v2[4] = 2; std::ofstream(dir + "/v2.null", std::ios::binary) << v2;
      dart::NullModel x; expect(!dart::read_null_model(dir + "/v2.null", x, nullptr, &err) && err.find("version") != std::string::npos, "version 2 refused: " + err); }
    // Pooling two sidecars equals one fit on the summed moments.
    dart::DecoyStats s; s.add(r); s.add(r);
    dart::NullModel p = dart::fit_null_model(s);
    std::vector<uint32_t> n2(6); std::vector<double> sum2(6), sq2(6);
    for (int b = 0; b < 6; ++b) { n2[b] = 2 * n[b]; sum2[b] = 2 * sum[b]; sq2[b] = 2 * sq[b]; }
    dart::NullModel d = dart::fit_null_model_moments(n2, sum2, sq2, 200);
    expect(p.mu == d.mu && p.beta == d.beta && p.n == d.n, "pooling sidecars equals pooling moments");
    std::cout << (failures ? "FAILED\n" : "ok\n");
}

void test_predict_decoys(const std::string& dart, const std::string& tut, const std::string& dir) {
    std::cout << "predict decoys... " << std::flush;
    const std::string faa = dir + "/p.faa";
    const std::string cmd = dart + " predict -i " + tut + "/reads.fq.gz -o " + dir + "/p.gff --fasta-aa " + faa + " --allow-self-profile -t 4";
    expect(run(cmd, dir + "/predict.log") == 0, "predict exits 0");
    const std::string decoy = dir + "/p.decoy.faa";
    expect(exists(decoy), "decoy file p.decoy.faa written");
    auto real = read_fasta(faa);
    auto dec = read_fasta(decoy);
    std::unordered_map<std::string, std::string> real_seq(real.begin(), real.end());
    // Read ordinal per read id from the FASTQ.
    std::unordered_map<std::string, size_t> ordinal;
    { auto lines = gz_lines(tut + "/reads.fq.gz"); size_t k = 0;
      for (size_t i = 0; i + 3 < lines.size(); i += 4) { const std::string& h = lines[i]; const size_t sp = h.find(' '); ordinal[h.substr(1, sp == std::string::npos ? std::string::npos : sp - 1)] = k++; } }
    expect(!real.empty() && !dec.empty(), "proteins and decoys present");
    size_t bad_prefix = 0, bad_seq = 0, bad_ordinal = 0;
    std::unordered_set<std::string> decoy_ids;
    for (const auto& kv : dec) {
        if (kv.first.rfind("DECOY_", 0) != 0) { ++bad_prefix; continue; }
        const std::string id = kv.first.substr(6);
        decoy_ids.insert(id);
        auto it = real_seq.find(id);
        if (it == real_seq.end() || std::string(it->second.rbegin(), it->second.rend()) != kv.second) ++bad_seq;
        auto o = ordinal.find(read_of(id));
        if (o == ordinal.end() || o->second % 10 != 0) ++bad_ordinal;
    }
    size_t missing = 0;
    for (const auto& kv : real) {
        auto o = ordinal.find(read_of(kv.first));
        if (o != ordinal.end() && o->second % 10 == 0 && !decoy_ids.count(kv.first)) ++missing;
    }
    expect(bad_prefix == 0, "every decoy id has the DECOY_ prefix");
    expect(bad_seq == 0, "every decoy is the reversed real protein (" + std::to_string(bad_seq) + " differ)");
    expect(bad_ordinal == 0, "every decoy comes from a read with ordinal % 10 == 0 (" + std::to_string(bad_ordinal) + " do not)");
    expect(missing == 0, "every protein of an ordinal % 10 == 0 read has a decoy (" + std::to_string(missing) + " missing)");
    expect(dec.size() == decoy_ids.size(), "decoy ids are unique");
    // --decoy-every 0: no decoy file.
    const std::string faa0 = dir + "/p0.faa";
    expect(run(dart + " predict -i " + tut + "/reads.fq.gz -o " + dir + "/p0.gff --fasta-aa " + faa0 + " --allow-self-profile --decoy-every 0 -t 4", dir + "/predict0.log") == 0, "predict --decoy-every 0 exits 0");
    expect(!exists(dir + "/p0.decoy.faa"), "--decoy-every 0 writes no decoy file");
    expect(slurp(faa0) == slurp(faa), "--decoy-every 0 leaves the protein file unchanged");
    std::cout << (failures ? "FAILED\n" : "ok\n");
}

void test_hits2emi_and_annotate(const std::string& dart, const std::string& tut, const std::string& dir) {
    std::cout << "hits2emi, dart null, damage-annotate... " << std::flush;
    auto lines = gz_lines(tut + "/expected/hits.tsv.gz");
    expect(lines.size() > 1000, "tutorial hits present");
    // Real hits as-is; decoys: every 7th query renamed DECOY_, bit score halved.
    std::ofstream real(dir + "/real.tsv"), decoy(dir + "/decoy.tsv");
    decoy << "query\ttarget\tfident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits\tqlen\ttlen\tqaln\ttaln\n";
    std::string last_q; size_t qi = 0; uint64_t decoy_rows = 0;
    for (const auto& l : lines) {
        if (l.empty() || l.rfind("query", 0) == 0) continue;
        real << l << "\n";
        const size_t t1 = l.find('\t');
        const std::string q = l.substr(0, t1);
        if (q != last_q) { last_q = q; ++qi; }
        if (qi % 7 != 0) continue;
        std::vector<std::string> f; std::stringstream ss(l); std::string x;
        while (std::getline(ss, x, '\t')) f.push_back(x);
        if (f.size() < 16) continue;
        f[0] = "DECOY_" + f[0];
        f[11] = std::to_string(std::strtod(f[11].c_str(), nullptr) * 0.5);
        for (size_t i = 0; i < f.size(); ++i) decoy << (i ? "\t" : "") << f[i];
        decoy << "\n";
        ++decoy_rows;
    }
    real.close(); decoy.close();
    const std::string agd = tut + "/expected/predictions.agd";
    const std::string a = dir + "/a.emi", b = dir + "/b.emi";
    expect(run(dart + " hits2emi -i " + dir + "/real.tsv -o " + a + " --damage-index " + agd + " -t 4", dir + "/h2e_a.log") == 0, "hits2emi without decoys exits 0");
    expect(run(dart + " hits2emi -i " + dir + "/real.tsv -o " + b + " --damage-index " + agd + " --decoy-hits " + dir + "/decoy.tsv --db-label tutorial -t 4", dir + "/h2e_b.log") == 0, "hits2emi --decoy-hits exits 0");
    expect(slurp(a) == slurp(b), "EMI is byte-identical with and without --decoy-hits");
    expect(!exists(a + ".null"), "no sidecar without --decoy-hits");
    expect(exists(b + ".null"), "sidecar b.emi.null written");
    dart::NullModel mb; dart::NullProvenance pb; std::string err;
    expect(dart::read_null_model(b + ".null", mb, &pb, &err), "sidecar reads: " + err);
    expect(pb.db_label == "tutorial" && pb.n_decoy_queries > 0 && pb.inputs.size() == 1, "sidecar provenance");
    expect(pb.n_decoy_hits == decoy_rows, "the header line is not counted as a decoy hit (" + std::to_string(pb.n_decoy_hits) + " vs " + std::to_string(decoy_rows) + ")");
    expect(run(dart + " null -o " + dir + "/n.null --db-label tutorial " + dir + "/decoy.tsv", dir + "/null.log") == 0, "dart null exits 0");
    dart::NullModel mn; dart::NullProvenance pn;
    expect(dart::read_null_model(dir + "/n.null", mn, &pn, &err), "dart null sidecar reads: " + err);
    expect(mn.mu == mb.mu && mn.beta == mb.beta && mn.n == mb.n && mn.sum == mb.sum, "dart null fits the same bins as hits2emi --decoy-hits");
    expect(run(dart + " damage-info " + b + ".null", dir + "/info.log") == 0 && slurp(dir + "/info.log").find("decoy queries") != std::string::npos, "damage-info prints the sidecar");
    // Annotate: auto-discovered sidecar == explicit --em-null-model.
    const std::string common = " --damage-index " + agd + " -t 4";
    expect(run(dart + " damage-annotate --emi " + b + " -o " + dir + "/ann1.tsv" + common, dir + "/ann1.log") == 0, "annotate with auto-discovered sidecar exits 0");
    expect(slurp(dir + "/ann1.log").find("sidecar next to the EMI") != std::string::npos, "annotate reports the sidecar source");
    expect(run(dart + " damage-annotate --emi " + b + " -o " + dir + "/ann2.tsv --em-null-model " + b + ".null" + common, dir + "/ann2.log") == 0, "annotate with --em-null-model exits 0");
    expect(slurp(dir + "/ann1.tsv") == slurp(dir + "/ann2.tsv"), "auto-discovered and explicit null give byte-identical output");
    expect(run(dart + " damage-annotate --emi " + a + " -o " + dir + "/ann3.tsv" + common, dir + "/ann3.log") != 0, "annotate without a null and without --no-null refuses");
    expect(slurp(dir + "/ann3.log").find("refusing EM annotation without a null model") != std::string::npos, "refusal message");
    std::cout << (failures ? "FAILED\n" : "ok\n");
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: " << argv[0] << " <dart> <examples/tutorial>\n"; return 2; }
    const std::string dart = argv[1], tut = argv[2];
    const std::string dir = make_tmpdir();
    test_roundtrip(dir);
    test_predict_decoys(dart, tut, dir);
    test_hits2emi_and_annotate(dart, tut, dir);
    if (failures) { std::cerr << failures << " failure(s); outputs kept in " << dir << "\n"; return 1; }
    std::system(("rm -rf " + dir).c_str());
    std::cout << "PASSED\n";
    return 0;
}
