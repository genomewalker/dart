/**
 * @file cmd_ref_stats.cpp
 * @brief Per-reference coverage/identity stats from a search-hits TSV (ONE sample/library).
 *
 * Emits per (sample, target): ref_length, n_reads, n_alns, breadth (covered/ref_length, 0-1),
 * depth_mean (aligned residues / covered), avg_identity (mean fident). breadth is a TRUE fraction
 * because the hits TSV carries tlen (col 14 of the LOCKED mmseqs FMT). Run once per library (no
 * cross-lib pooling): lean tlen-bit bitset for the covered union + running sums, streaming, low memory.
 *
 * Input FMT (LOCKED, 0-based fields used): target=1, fident=2, alnlen=3, tstart=8, tend=9, tlen=13.
 *   query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln
 */
#include "subcommand.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Acc {
    uint32_t tlen = 0;
    std::vector<uint64_t> bits;   // covered-position bitset (lazy, size tlen)
    double   sum_fident = 0.0;
    uint64_t sum_alnlen = 0;
    uint32_t n_reads = 0;         // distinct queries
    uint32_t n_alns  = 0;
    uint64_t last_q  = 0;         // 1-based ordinal of the last query counted (0 = none)

    void ensure(uint32_t L) {
        if (bits.empty() && L > 0) { tlen = L; bits.assign((static_cast<size_t>(L) + 63) / 64, 0ULL); }
    }
    void cover(uint32_t s0, uint32_t e0) {            // [s0, e0) 0-based, clamped to tlen
        if (e0 > tlen) e0 = tlen;
        if (s0 > tlen) s0 = tlen;
        for (uint32_t p = s0; p < e0; ++p) bits[p >> 6] |= (1ULL << (p & 63));
    }
    uint32_t covered() const {
        uint32_t c = 0;
        for (uint64_t w : bits) c += static_cast<uint32_t>(__builtin_popcountll(w));
        return c;
    }
};

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " ref-stats --hits <hits.tsv> --sample <ID> --out <stats.tsv> [--min-breadth F]\n\n"
              << "Per-reference coverage/identity stats for ONE sample (library). breadth is covered/ref_length (0-1).\n"
              << "Input: search-hits TSV in the LOCKED 16-col mmseqs FMT (tlen must be column 14).\n\n"
              << "Options:\n"
              << "  --hits FILE       Search-hits TSV for one library (required)\n"
              << "  --sample ID       Sample/library id, written as first column (required)\n"
              << "  --out FILE        Output stats TSV (required)\n"
              << "  --min-breadth F   Only emit targets with breadth >= F (default 0)\n"
              << "  -h, --help        Show this help\n";
}

// Return pointer to field `idx` (0-based) start and set its length; returns nullptr if missing.
inline const char* field(const char* p, const char* end, int idx, size_t& len) {
    int f = 0;
    const char* s = p;
    while (p < end) {
        if (*p == '\t' || *p == '\n') {
            if (f == idx) { len = static_cast<size_t>(p - s); return s; }
            ++f; s = p + 1;
            if (*p == '\n') break;
        }
        ++p;
    }
    if (f == idx) { len = static_cast<size_t>(p - s); return s; }
    return nullptr;
}

}  // namespace

int cmd_ref_stats(int argc, char* argv[]) {
    std::string hits, sample, out;
    double min_breadth = 0.0;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--hits") && i + 1 < argc) hits = argv[++i];
        else if (!std::strcmp(argv[i], "--sample") && i + 1 < argc) sample = argv[++i];
        else if (!std::strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
        else if (!std::strcmp(argv[i], "--min-breadth") && i + 1 < argc) min_breadth = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-h") || !std::strcmp(argv[i], "--help")) { print_usage(argv[0]); return 0; }
    }
    if (hits.empty() || sample.empty() || out.empty()) { print_usage(argv[0]); return 1; }

    std::ifstream in(hits, std::ios::binary);
    if (!in) { std::cerr << "ref-stats: cannot open " << hits << "\n"; return 1; }

    std::unordered_map<std::string, Acc> acc;
    acc.reserve(1u << 20);

    std::string line, prev_q;
    line.reserve(512);
    uint64_t nrows = 0;
    // ceiling: n_reads is exact only when each query's rows are contiguous (mmseqs convertalis
    // order); a query split across non-adjacent blocks is counted once per block.
    uint64_t q_ord = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const char* p = line.data();
        const char* end = p + line.size();
        size_t l_t, l_fi, l_al, l_ts, l_te, l_tl;
        const char* f_t  = field(p, end, 1, l_t);
        const char* f_fi = field(p, end, 2, l_fi);
        const char* f_al = field(p, end, 3, l_al);
        const char* f_ts = field(p, end, 8, l_ts);
        const char* f_te = field(p, end, 9, l_te);
        const char* f_tl = field(p, end, 13, l_tl);
        if (!f_t || !f_ts || !f_te || !f_tl) continue;   // malformed / short line
        uint32_t tstart = static_cast<uint32_t>(std::strtoul(f_ts, nullptr, 10));
        uint32_t tend   = static_cast<uint32_t>(std::strtoul(f_te, nullptr, 10));
        uint32_t tlen   = static_cast<uint32_t>(std::strtoul(f_tl, nullptr, 10));
        if (tlen == 0 || tend < tstart) continue;
        const size_t l_q = static_cast<size_t>(f_t - p) - 1;
        if (q_ord == 0 || prev_q.compare(0, std::string::npos, p, l_q) != 0) { prev_q.assign(p, l_q); ++q_ord; }
        double fident = f_fi ? std::strtod(f_fi, nullptr) : 0.0;
        uint32_t alnlen = f_al ? static_cast<uint32_t>(std::strtoul(f_al, nullptr, 10)) : (tend - tstart + 1);

        Acc& a = acc[std::string(f_t, l_t)];
        a.ensure(tlen);
        // tstart/tend are 1-based inclusive on the target -> [tstart-1, tend)
        uint32_t s0 = tstart > 0 ? tstart - 1 : 0;
        a.cover(s0, tend);
        a.sum_fident += fident;
        a.sum_alnlen += alnlen;
        if (a.last_q != q_ord) { a.last_q = q_ord; ++a.n_reads; }
        ++a.n_alns;
        ++nrows;
    }

    std::ofstream os(out, std::ios::binary);
    if (!os) { std::cerr << "ref-stats: cannot write " << out << "\n"; return 1; }
    os << "sample_id\tref_id\tref_length\tn_reads\tn_alns\tbreadth\tdepth_mean\tavg_identity\n";
    uint64_t emitted = 0;
    char buf[512];
    for (const auto& kv : acc) {
        const Acc& a = kv.second;
        if (a.tlen == 0) continue;
        uint32_t cov = a.covered();
        double breadth = static_cast<double>(cov) / static_cast<double>(a.tlen);
        if (breadth < min_breadth) continue;
        double depth = cov > 0 ? static_cast<double>(a.sum_alnlen) / static_cast<double>(cov) : 0.0;
        double avg_id = a.n_alns > 0 ? a.sum_fident / static_cast<double>(a.n_alns) : 0.0;
        int m = std::snprintf(buf, sizeof(buf), "%s\t%s\t%u\t%u\t%u\t%.4f\t%.3f\t%.4f\n",
                              sample.c_str(), kv.first.c_str(), a.tlen, a.n_reads, a.n_alns,
                              breadth, depth, avg_id);
        os.write(buf, m);
        ++emitted;
    }
    std::cerr << "ref-stats " << sample << ": rows=" << nrows << " targets=" << acc.size()
              << " emitted=" << emitted << " -> " << out << "\n";
    return 0;
}

namespace {
    struct RefStatsRegistrar {
        RefStatsRegistrar() {
            dart::cli::SubcommandRegistry::instance().register_command(
                "ref-stats",
                "Per-reference coverage/identity stats from a hits TSV (one sample)",
                cmd_ref_stats, 46);
        }
    };
    static RefStatsRegistrar registrar;
}
