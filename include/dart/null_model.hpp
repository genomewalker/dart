// Empirical null model sidecar (.null): decoy-fitted Gumbel per aligned-length bin.
//
// predict writes reversed decoy peptides next to --fasta-aa; the user searches them with
// the real search's parameters; hits2emi --decoy-hits (or dart null) reduces the decoy
// hits to the best bit score per decoy query, bins by aligned length, fits a Gumbel per
// bin and writes the sidecar next to the EMI. damage-annotate finds <emi>.null by itself.
//
// The sidecar is on the raw bit scale; annotate multiplies mu and beta by its bit-to-nat
// lambda at load (Gumbel location and scale are linear in a positive multiplier).
#pragma once

#include "dart/em_reassign.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace dart {

struct NullProvenance {
    struct Input {
        std::string path;
        uint64_t size = 0;
        uint64_t mtime = 0;
    };
    std::string dart_version;
    uint64_t created_utc = 0;
    uint64_t n_decoy_queries = 0;
    uint64_t n_decoy_hits = 0;
    std::string db_label;
    std::vector<Input> inputs;
    bool legacy_tsv = false;   // set by read_null_model for the deprecated TSV sidecar; not written
};

// Per-bin moments of the best decoy hit per query, raw bit scale. Enough to refit and to
// pool sidecars exactly.
struct DecoyStats {
    std::vector<uint32_t> n      = std::vector<uint32_t>(kNumNullBins, 0u);
    std::vector<double>   sum    = std::vector<double>(kNumNullBins, 0.0);
    std::vector<double>   sum_sq = std::vector<double>(kNumNullBins, 0.0);
    uint64_t n_queries = 0;
    uint64_t n_hits = 0;

    void add_query(uint16_t aln_len, double bits);
    void add(const DecoyStats& o);
    void add(const NullModel& m);   // from a sidecar's stored moments
};

// One decoy hits file (16-column MMseqs2 TSV, plain or .gz). Rule: when the file has any
// query id starting with DECOY_, only those rows are decoys; otherwise every row is.
// A first line whose query field is "query" (an MMseqs2 header) is skipped.
// Returns the number of decoy queries taken from this file; *skipped gets the number of
// unprefixed queries ignored and *by_prefix whether the DECOY_ rule selected the rows.
// Appends the file to prov.inputs and adds its hits to prov.
uint64_t accumulate_decoy_hits(const std::string& path, DecoyStats& stats,
                               NullProvenance& prov, uint64_t* skipped = nullptr,
                               bool* by_prefix = nullptr);

NullModel fit_null_model(const DecoyStats& s, uint32_t min_n = 200);

// Binary sidecar: "DNUL" ... "LUND", little-endian, see wiki/File-Formats-and-Memory.md.
bool write_null_model(const std::string& path, const NullModel& m,
                      const NullProvenance& prov, std::string* err = nullptr);

// Dispatches on the leading magic: binary sidecar, else the deprecated TSV
// (bin lo hi mu beta n, already on the lambda scale; prov->legacy_tsv is set).
bool read_null_model(const std::string& path, NullModel& m,
                     NullProvenance* prov = nullptr, std::string* err = nullptr);

bool is_null_sidecar(const std::string& path);   // leading "DNUL"

inline std::string null_sidecar_path(const std::string& emi_path) { return emi_path + ".null"; }

// How each bin got its Gumbel: its own fit (n >= min_n), the fit pooled over every bin, or
// none (beta 0: fewer than two queries or zero variance in the pool). With no fitted bin
// annotate runs EM in memory with the null coordinate off (read_null_log_density returns
// -inf), which gives the output of --no-null without the tape.
struct NullFitCounts { int own = 0, pooled = 0, unfit = 0; };
NullFitCounts null_fit_counts(const NullModel& m, uint32_t min_n = 200);
inline int fitted_bins(const NullModel& m) { const NullFitCounts c = null_fit_counts(m); return c.own + c.pooled; }
std::string null_fit_summary(const NullModel& m);   // "k own, p pooled, u unfit of 6 bins"

void print_null_model(std::ostream& os, const NullModel& m, const NullProvenance& prov);

}  // namespace dart
