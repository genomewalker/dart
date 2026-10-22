#include "dart/null_model.hpp"
#include "dart/fast_columnar_writer.hpp"   // ChunkedFileReader (plain or .gz, skips the header)

#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <sys/stat.h>

namespace dart {

namespace {

constexpr char kMagic[4]    = {'D', 'N', 'U', 'L'};
constexpr char kMagicEnd[4] = {'L', 'U', 'N', 'D'};
constexpr uint16_t kVersion = 2;   // 2 adds n_decoy_submitted after n_decoy_hits; 1 still reads
constexpr const char* kPrefix = "DECOY_";
constexpr size_t kPrefixLen = 6;

template <class T> void put(std::string& b, T v) { b.append(reinterpret_cast<const char*>(&v), sizeof(T)); }
void put_str(std::string& b, const std::string& s) {
    put<uint32_t>(b, static_cast<uint32_t>(s.size()));
    b.append(s);
}

struct Cursor {
    const char* p; const char* end; bool ok = true;
    template <class T> T get() {
        T v{};
        if (end - p < static_cast<ptrdiff_t>(sizeof(T))) { ok = false; return v; }
        std::memcpy(&v, p, sizeof(T)); p += sizeof(T); return v;
    }
    std::string get_str() {
        const uint32_t n = get<uint32_t>();
        if (!ok || end - p < static_cast<ptrdiff_t>(n)) { ok = false; return {}; }
        std::string s(p, n); p += n; return s;
    }
};

// Method-of-moments Gumbel from (n, sum, sum_sq): beta = sqrt(6)*sd/pi, mu = mean - gamma*beta,
// sample variance with n-1, as fit_null_model always did.
bool gumbel_from_moments(uint32_t n, double sum, double sum_sq, float& mu, float& beta) {
    if (n < 2) return false;
    const double mean = sum / n;
    const double var = (sum_sq - n * mean * mean) / (n - 1.0);
    if (!(var > 0.0)) return false;
    const double b = std::sqrt(6.0) * std::sqrt(var) / M_PI;
    beta = static_cast<float>(b);
    mu = static_cast<float>(mean - 0.5772156649015329 * b);
    return true;
}

struct QMax { float bits; uint16_t aln_len; };

void reduce_into(const std::unordered_map<std::string, QMax>& qmax, DecoyStats& stats) {
    for (const auto& kv : qmax) stats.add_query(kv.second.aln_len, kv.second.bits);
}

}  // namespace

void DecoyStats::add_query(uint16_t aln_len, double bits) {
    const int b = null_length_bin(aln_len);
    n[b] += 1; sum[b] += bits; sum_sq[b] += bits * bits; n_queries += 1;
}

void DecoyStats::add(const DecoyStats& o) {
    for (int b = 0; b < kNumNullBins; ++b) { n[b] += o.n[b]; sum[b] += o.sum[b]; sum_sq[b] += o.sum_sq[b]; }
    n_queries += o.n_queries; n_hits += o.n_hits;
}

void DecoyStats::add(const NullModel& m) {
    for (int b = 0; b < kNumNullBins && b < static_cast<int>(m.n.size()); ++b) {
        n[b] += m.n[b];
        if (b < static_cast<int>(m.sum.size())) { sum[b] += m.sum[b]; sum_sq[b] += m.sum_sq[b]; }
        n_queries += m.n[b];
    }
}

uint64_t accumulate_decoy_hits(const std::string& path, DecoyStats& stats,
                               NullProvenance& prov, uint64_t* skipped, bool* by_prefix) {
    std::unordered_map<std::string, QMax> prefixed, plain;
    uint64_t hits = 0;
    ChunkedFileReader reader(path);
    std::string field;
    bool first = true;
    while (true) {
        std::string_view line = reader.next_line();
        if (line.empty()) break;
        if (first) {
            first = false;
            if (line.size() >= 6 && line.compare(0, 6, "query\t") == 0) continue;   // MMseqs2 header
        }
        const char* p = line.data();
        const char* end = p + line.size();
        const char* f[13];
        const char* fp = p;
        bool complete = true;
        for (int i = 0; i < 12; ++i) {
            f[i] = fp;
            const char* tab = simd_find_tab(fp, end);
            if (tab >= end) { complete = false; break; }
            fp = tab + 1;
        }
        if (!complete) continue;
        f[12] = fp;
        const char* bits_end = simd_find_tab(fp, end);
        ++hits;
        const std::string_view q(f[0], static_cast<size_t>(f[1] - f[0] - 1));
        field.assign(f[3], static_cast<size_t>(f[4] - f[3] - 1));
        const uint32_t alnlen = static_cast<uint32_t>(std::strtoul(field.c_str(), nullptr, 10));
        field.assign(f[11], static_cast<size_t>(bits_end - f[11]));
        const float bits = std::strtof(field.c_str(), nullptr);
        auto& map = (q.size() >= kPrefixLen && q.compare(0, kPrefixLen, kPrefix) == 0) ? prefixed : plain;
        auto it = map.find(std::string(q));
        if (it == map.end()) {
            map.emplace(std::string(q), QMax{bits, static_cast<uint16_t>(std::min(alnlen, 65535u))});
        } else if (bits > it->second.bits) {
            it->second = QMax{bits, static_cast<uint16_t>(std::min(alnlen, 65535u))};
        }
    }
    const auto& take = prefixed.empty() ? plain : prefixed;
    if (skipped) *skipped = prefixed.empty() ? 0 : plain.size();
    if (by_prefix) *by_prefix = !prefixed.empty();
    reduce_into(take, stats);
    stats.n_hits += hits;
    prov.n_decoy_queries += take.size();
    prov.n_decoy_hits += hits;
    NullProvenance::Input in;
    in.path = path;
    struct stat st;
    if (stat(path.c_str(), &st) == 0) { in.size = static_cast<uint64_t>(st.st_size); in.mtime = static_cast<uint64_t>(st.st_mtime); }
    prov.inputs.push_back(in);
    return take.size();
}

NullModel fit_null_model(const DecoyStats& s, uint32_t min_n) {
    return fit_null_model_moments(s.n, s.sum, s.sum_sq, min_n);
}

bool write_null_model(const std::string& path, const NullModel& m,
                      const NullProvenance& prov, std::string* err) {
    const int nb = static_cast<int>(m.mu.size());
    if (nb != kNumNullBins || m.beta.size() != m.mu.size() || m.n.size() != m.mu.size()) {
        if (err) *err = "null model has " + std::to_string(nb) + " bins, expected " + std::to_string(kNumNullBins);
        return false;
    }
    std::string b;
    b.append(kMagic, 4);
    put<uint16_t>(b, kVersion);
    put<uint16_t>(b, 0);
    put<uint16_t>(b, static_cast<uint16_t>(nb));
    for (int i = 0; i <= nb; ++i) put<uint16_t>(b, kNullBinEdges[i]);
    for (int i = 0; i < nb; ++i) {
        put<double>(b, m.mu[i]);
        put<double>(b, m.beta[i]);
        put<uint32_t>(b, m.n[i]);
        put<double>(b, i < static_cast<int>(m.sum.size()) ? m.sum[i] : 0.0);
        put<double>(b, i < static_cast<int>(m.sum_sq.size()) ? m.sum_sq[i] : 0.0);
    }
    put_str(b, prov.dart_version);
    put<uint64_t>(b, prov.created_utc ? prov.created_utc : static_cast<uint64_t>(std::time(nullptr)));
    put<uint64_t>(b, prov.n_decoy_queries);
    put<uint64_t>(b, prov.n_decoy_hits);
    put<uint64_t>(b, prov.n_decoy_submitted);
    put_str(b, prov.db_label);
    put<uint16_t>(b, static_cast<uint16_t>(std::min<size_t>(prov.inputs.size(), 65535)));
    for (size_t i = 0; i < prov.inputs.size() && i < 65535; ++i) {
        put_str(b, prov.inputs[i].path);
        put<uint64_t>(b, prov.inputs[i].size);
        put<uint64_t>(b, prov.inputs[i].mtime);
    }
    b.append(kMagicEnd, 4);
    std::ofstream out(path, std::ios::binary);
    if (!out || !out.write(b.data(), static_cast<std::streamsize>(b.size()))) {
        if (err) *err = "cannot write " + path;
        return false;
    }
    return true;
}

bool is_null_sidecar(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    char magic[4];
    return in && in.read(magic, 4) && std::memcmp(magic, kMagic, 4) == 0;
}

static bool read_null_model_tsv(const std::string& path, NullModel& m) {
    std::ifstream in(path);
    if (!in) return false;
    m.mu.clear(); m.beta.clear(); m.n.clear(); m.sum.clear(); m.sum_sq.clear();
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (first) { first = false; if (line.rfind("bin", 0) == 0) continue; }
        std::istringstream ss(line);
        int bin; double lo, hi, mu, beta; long n = 0;
        if (!(ss >> bin >> lo >> hi >> mu >> beta)) continue;
        ss >> n;
        m.mu.push_back(static_cast<float>(mu));
        m.beta.push_back(static_cast<float>(beta));
        m.n.push_back(static_cast<uint32_t>(std::max(0L, n)));
    }
    return m.valid();
}

bool read_null_model(const std::string& path, NullModel& m, NullProvenance* prov, std::string* err) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { if (err) *err = "cannot open " + path; return false; }
    std::string buf((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (buf.size() < 4 || std::memcmp(buf.data(), kMagic, 4) != 0) {
        if (!read_null_model_tsv(path, m)) { if (err) *err = path + " is neither a .null sidecar nor a null TSV"; return false; }
        if (prov) { *prov = NullProvenance{}; prov->legacy_tsv = true; }
        return true;
    }
    if (buf.size() < 8 || std::memcmp(buf.data() + buf.size() - 4, kMagicEnd, 4) != 0) {
        if (err) *err = path + ": truncated or corrupt null sidecar (bad trailing magic)";
        return false;
    }
    Cursor c{buf.data() + 4, buf.data() + buf.size() - 4};
    const uint16_t version = c.get<uint16_t>();
    c.get<uint16_t>();   // flags
    if (version < 1 || version > kVersion) {
        if (err) *err = path + ": null sidecar version " + std::to_string(version) + ", this build reads versions 1-" + std::to_string(kVersion);
        return false;
    }
    const uint16_t nb = c.get<uint16_t>();
    if (nb != kNumNullBins) { if (err) *err = path + ": " + std::to_string(nb) + " bins, this build uses " + std::to_string(kNumNullBins); return false; }
    for (int i = 0; i <= nb; ++i) {
        const uint16_t e = c.get<uint16_t>();
        if (c.ok && e != kNullBinEdges[i]) { if (err) *err = path + ": bin edge " + std::to_string(e) + " at " + std::to_string(i) + " differs from this build's " + std::to_string(kNullBinEdges[i]); return false; }
    }
    m.mu.assign(nb, 0.0f); m.beta.assign(nb, 0.0f); m.n.assign(nb, 0u); m.sum.assign(nb, 0.0); m.sum_sq.assign(nb, 0.0);
    for (int i = 0; i < nb; ++i) {
        m.mu[i] = static_cast<float>(c.get<double>());
        m.beta[i] = static_cast<float>(c.get<double>());
        m.n[i] = c.get<uint32_t>();
        m.sum[i] = c.get<double>();
        m.sum_sq[i] = c.get<double>();
    }
    NullProvenance p;
    p.dart_version = c.get_str();
    p.created_utc = c.get<uint64_t>();
    p.n_decoy_queries = c.get<uint64_t>();
    p.n_decoy_hits = c.get<uint64_t>();
    if (version >= 2) p.n_decoy_submitted = c.get<uint64_t>();
    p.db_label = c.get_str();
    const uint16_t ni = c.get<uint16_t>();
    for (uint16_t i = 0; i < ni && c.ok; ++i) {
        NullProvenance::Input in;
        in.path = c.get_str(); in.size = c.get<uint64_t>(); in.mtime = c.get<uint64_t>();
        p.inputs.push_back(in);
    }
    if (!c.ok || c.p != c.end) { if (err) *err = path + ": null sidecar body does not match its header"; return false; }
    if (prov) *prov = p;
    return true;
}

NullFitCounts null_fit_counts(const NullModel& m, uint32_t min_n) {
    NullFitCounts c;
    for (size_t b = 0; b < m.beta.size(); ++b) {
        if (!(m.beta[b] > 0.0f)) ++c.unfit;
        else if (b < m.n.size() && m.n[b] >= min_n) ++c.own;
        else ++c.pooled;
    }
    return c;
}

std::string null_fit_summary(const NullModel& m) {
    const NullFitCounts c = null_fit_counts(m);
    return std::to_string(c.own) + " own, " + std::to_string(c.pooled) + " pooled, " +
           std::to_string(c.unfit) + " unfit of " + std::to_string(m.beta.size()) + " bins";
}

void print_null_model(std::ostream& os, const NullModel& m, const NullProvenance& prov) {
    os << "=== DART null model (.null) ===\n";
    if (prov.legacy_tsv) os << "format: deprecated TSV sidecar (lambda scale, no provenance)\n";
    else {
        os << "dart: " << prov.dart_version << "\n";
        char buf[32]; const std::time_t t = static_cast<std::time_t>(prov.created_utc);
        std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S UTC", std::gmtime(&t));
        os << "created: " << buf << "\n";
        os << "decoy queries: " << prov.n_decoy_queries << " (" << prov.n_decoy_hits << " hits)\n";
        os << "decoys submitted: " << (prov.n_decoy_submitted ? std::to_string(prov.n_decoy_submitted) : "unknown") << "\n";
        os << "db label: " << (prov.db_label.empty() ? "-" : prov.db_label) << "\n";
        for (const auto& in : prov.inputs) os << "input: " << in.path << " size=" << in.size << " mtime=" << in.mtime << "\n";
    }
    os << "fit: " << null_fit_summary(m) << "\n";
    os << "bin\tlo\thi\tn\tmu\tbeta\n";
    for (size_t b = 0; b < m.mu.size(); ++b)
        os << b << '\t' << kNullBinEdges[b] << '\t' << kNullBinEdges[b + 1] << '\t' << m.n[b] << '\t'
           << std::setprecision(6) << m.mu[b] << '\t' << m.beta[b] << '\n';
}

uint64_t count_fasta_records(const std::string& path) {
    ChunkedFileReader reader(path);
    uint64_t n = 0;
    while (true) {
        const std::string_view line = reader.next_line();
        if (line.data() == nullptr) break;   // a blank line is an empty view into the buffer, EOF is {}
        if (!line.empty() && line[0] == '>') ++n;
    }
    return n;
}

uint64_t null_decoy_denominator(const NullModel& m, const NullProvenance& prov, bool* submitted) {
    if (submitted) *submitted = prov.n_decoy_submitted > 0;
    if (prov.n_decoy_submitted > 0) return prov.n_decoy_submitted;
    if (prov.n_decoy_queries > 0) return prov.n_decoy_queries;
    uint64_t n = 0;
    for (uint32_t v : m.n) n += v;
    return n;
}

namespace {
bool chance_bin(const NullModel& m, uint64_t n_decoys, int b, double& w) {
    if (b < 0 || b >= static_cast<int>(m.beta.size()) || b >= static_cast<int>(m.n.size())) return false;
    if (!(m.beta[b] > 0.0f) || n_decoys == 0) return false;
    w = static_cast<double>(m.n[b]) / static_cast<double>(n_decoys);
    return true;
}
}  // namespace

double null_chance_p(const NullModel& m, double scale, uint64_t n_decoys, float bits, uint16_t aln_len) {
    const int b = null_length_bin(aln_len);
    double w = 0.0;
    if (!chance_bin(m, n_decoys, b, w)) return std::numeric_limits<double>::quiet_NaN();
    const double z = (scale * static_cast<double>(bits) - m.mu[b]) / m.beta[b];
    const double tail = (z <= -40.0) ? 1.0 : -std::expm1(-std::exp(-z));
    return w * tail;
}

double null_chance_bits_at(const NullModel& m, double scale, uint64_t n_decoys, int b, double p) {
    double w = 0.0;
    if (!chance_bin(m, n_decoys, b, w) || !(p < w) || !(p > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double z = -std::log(-std::log1p(-p / w));
    return (m.mu[b] + m.beta[b] * z) / scale;
}

}  // namespace dart
