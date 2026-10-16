/**
 * @file cmd_ref_stats.cpp
 * @brief Per-reference coverage/identity stats for ONE sample/library, from a hits TSV or an EMI v7.
 *
 * Emits per (sample, target): ref_length, n_reads, n_alns, breadth (covered/ref_length, 0-1),
 * depth_mean (aligned residues / covered), avg_identity (mean fident). breadth is a TRUE fraction
 * because both inputs carry tlen. Run once per library (no cross-lib pooling): lean tlen-bit bitset
 * for the covered union + running sums.
 *
 * --hits FMT (LOCKED, 0-based fields used): query=0, target=1, fident=2, alnlen=3, tstart=8, tend=9, tlen=13.
 *   query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln
 * --emi: hits_v7.emi columns ref_idx, identity_q (fident*65535), aln_len, tstart, tend, tlen.
 */
#include "subcommand.hpp"
#include "dart/emi_v7_io.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
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
    uint32_t covered(uint32_t lo, uint32_t hi) const {   // covered positions in [lo, hi)
        uint32_t c = 0;
        for (uint32_t p = lo; p < hi; ++p) c += (bits[p >> 6] >> (p & 63)) & 1ULL;
        return c;
    }
};

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " ref-stats (--hits <hits.tsv> | --emi <hits_v7.emi>) --sample <ID> --out <stats.tsv>\n"
              << "                 [--min-breadth F] [--trim] [--threads N]\n\n"
              << "Per-reference coverage/identity stats for ONE sample (library). breadth is covered/ref_length (0-1).\n\n"
              << "Options:\n"
              << "  --hits FILE       Search-hits TSV in the LOCKED 16-col mmseqs FMT (tlen = column 14)\n"
              << "  --emi FILE        EMI v7 hits file for one library (alternative to --hits)\n"
              << "  --sample ID       Sample/library id, written as first column (required)\n"
              << "  --out FILE        Output stats TSV (required)\n"
              << "  --min-breadth F   Only emit targets with breadth >= F (default 0)\n"
              << "  --trim            x-filter breadth: ignore ceil(avg alnlen)+10 positions at each end,\n"
              << "                    still divide by ref_length (default: untrimmed)\n"
              << "  --threads N       Threads for --emi block decode (default: all)\n"
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

// Writes one stats row; returns false when filtered out by min_breadth.
bool write_row(std::ofstream& os, const std::string& sample, const std::string& ref, const Acc& a,
               bool trim, double min_breadth) {
    uint32_t cov = a.covered();
    double breadth = static_cast<double>(cov) / static_cast<double>(a.tlen);
    if (trim) {   // x-filter trim_coverage_by_subject: skipped when the trim reaches tlen
        const uint64_t t = (a.sum_alnlen + a.n_alns - 1) / a.n_alns + 10;
        if (t < a.tlen) {
            const uint32_t lo = static_cast<uint32_t>(t), hi = a.tlen - lo;
            breadth = lo < hi ? static_cast<double>(a.covered(lo, hi)) / static_cast<double>(a.tlen) : 0.0;
        }
    }
    if (breadth < min_breadth) return false;
    double depth = cov > 0 ? static_cast<double>(a.sum_alnlen) / static_cast<double>(cov) : 0.0;
    double avg_id = a.n_alns > 0 ? a.sum_fident / static_cast<double>(a.n_alns) : 0.0;
    char buf[512];
    int m = std::snprintf(buf, sizeof(buf), "%s\t%s\t%u\t%u\t%u\t%.4f\t%.3f\t%.4f\n",
                          sample.c_str(), ref.c_str(), a.tlen, a.n_reads, a.n_alns,
                          breadth, depth, avg_id);
    os.write(buf, m);
    return true;
}

}  // namespace

int cmd_ref_stats(int argc, char* argv[]) {
    std::string hits, emi, sample, out;
    double min_breadth = 0.0;
    bool trim = false;
    int threads = 0;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--hits") && i + 1 < argc) hits = argv[++i];
        else if (!std::strcmp(argv[i], "--emi") && i + 1 < argc) emi = argv[++i];
        else if (!std::strcmp(argv[i], "--sample") && i + 1 < argc) sample = argv[++i];
        else if (!std::strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
        else if (!std::strcmp(argv[i], "--min-breadth") && i + 1 < argc) min_breadth = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--trim")) trim = true;
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc) threads = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-h") || !std::strcmp(argv[i], "--help")) { print_usage(argv[0]); return 0; }
    }
    if (hits.empty() == emi.empty() || sample.empty() || out.empty()) { print_usage(argv[0]); return 1; }

    std::ofstream os(out, std::ios::binary);
    if (!os) { std::cerr << "ref-stats: cannot write " << out << "\n"; return 1; }
    os << "sample_id\tref_id\tref_length\tn_reads\tn_alns\tbreadth\tdepth_mean\tavg_identity\n";
    uint64_t nrows = 0, ntargets = 0, emitted = 0;

    if (!emi.empty()) try {
        if (threads > 0) omp_set_num_threads(threads);
        dart::emiv7::Reader rd(emi);
        std::vector<Acc> acc(rd.num_refs());
        // ceiling: one lock per alignment row; contention grows with threads. upgrade: per-thread
        // ref-sharded buffers merged per block.
        std::array<std::mutex, 1 << 16> locks;
        std::vector<uint64_t> rows_per_block(rd.num_blocks(), 0);
        rd.parallel_scan([&](uint32_t b, uint32_t n, const uint32_t* read, const uint32_t* ref,
                             const float*, const float*, const float*,
                             const uint16_t*, const uint16_t*, const float*, const float*,
                             const uint16_t* iq, const uint16_t* al, const uint16_t*, const uint16_t*,
                             const uint16_t* ts, const uint16_t* te, const uint16_t* tl, const uint16_t*,
                             const uint16_t*, const uint16_t*) {
            // Blocks are read-aligned, so a read's rows are contiguous here: n_reads counts each
            // (read, ref) pair once without cross-block state.
            std::unordered_set<uint32_t> seen;
            for (uint32_t j = 0; j < n; ++j) {
                if (j == 0 || read[j] != read[j - 1]) seen.clear();
                if (tl[j] == 0 || te[j] < ts[j]) continue;
                const bool first = seen.insert(ref[j]).second;
                std::lock_guard<std::mutex> g(locks[ref[j] & 0xFFFF]);
                Acc& a = acc[ref[j]];
                a.ensure(tl[j]);
                a.cover(ts[j] > 0 ? ts[j] - 1u : 0u, te[j]);   // 1-based inclusive -> [s-1, e)
                a.sum_fident += iq[j] / 65535.0;
                // ceiling: hits2emi-stream stores aln_len as max(q span, t span)+1 (mmseqs_db.hpp:76),
                // not mmseqs alnLength, so gapped rows run short in depth_mean and the --trim length.
                // upgrade: use get_taln(b, j).size() (decodes the cold aln strings).
                a.sum_alnlen += al[j];
                a.n_reads += first;
                ++a.n_alns;
            }
            rows_per_block[b] = n;
        });
        for (uint64_t r : rows_per_block) nrows += r;
        for (uint32_t i = 0; i < acc.size(); ++i) {
            if (acc[i].tlen == 0) continue;
            ++ntargets;
            emitted += write_row(os, sample, rd.ref_name(i), acc[i], trim, min_breadth);
        }
    } catch (const std::exception& e) {
        std::cerr << "ref-stats: " << emi << ": " << e.what() << "\n";
        return 1;
    } else {
        std::ifstream in(hits, std::ios::binary);
        if (!in) { std::cerr << "ref-stats: cannot open " << hits << "\n"; return 1; }
        std::unordered_map<std::string, Acc> acc;
        acc.reserve(1u << 20);
        std::string line, prev_q;
        line.reserve(512);
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
        for (const auto& kv : acc) {
            if (kv.second.tlen == 0) continue;
            ++ntargets;
            emitted += write_row(os, sample, kv.first, kv.second, trim, min_breadth);
        }
    }

    std::cerr << "ref-stats " << sample << ": rows=" << nrows << " targets=" << ntargets
              << " emitted=" << emitted << " -> " << out << "\n";
    return 0;
}

namespace {
    struct RefStatsRegistrar {
        RefStatsRegistrar() {
            dart::cli::SubcommandRegistry::instance().register_command(
                "ref-stats",
                "Per-reference coverage/identity stats from a hits TSV or EMI v7 (one sample)",
                cmd_ref_stats, 46);
        }
    };
    static RefStatsRegistrar registrar;
}
