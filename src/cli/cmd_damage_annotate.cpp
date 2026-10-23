// dart damage-annotate: Post-mapping damage annotation
//
// Given EMI alignments (with qaln/taln), identifies amino acid positions
// where the query (observed protein) differs from the target (reference)
// in a damage-consistent way. C->T deamination (5' end) and G->A (3' end)
// cause specific amino acid substitution patterns that can be identified
// with ~90% precision when a reference protein is available.

#include "subcommand.hpp"
#include "dart/version.h"
#include "dart/null_model.hpp"
#include "dart/damage_subs.hpp"
#include "dart/damage_stats.hpp"
#include "dart/damage_index_reader.hpp"
#include "dart/bayesian_damage_score.hpp"
#include "dart/em_reassign.hpp"
#include "dart/fast_tsv.hpp"
#include "dart/theta_d_coordinator.hpp"
#include "dart/coverage_em.hpp"
#include "dart/columnar_index.hpp"
#include "dart/emi_v7_io.hpp"
#include "dart/log_utils.hpp"
#include "dart/mmap_array.hpp"
#include "dart/damage_index.hpp"      // fnv1a_hash for the .dcx calib_contract_id
#include "dart/dcx_external_sort.hpp"
#include "dart/dart_columnar.hpp"     // .dcx interchange writer (--dcx)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <numeric>
#ifdef _OPENMP
#include <parallel/algorithm>
#endif
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <string_view>
#ifdef __linux__
#include <malloc.h>
#endif
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <cstdio>
#include <cerrno>
#include <sstream>
#include <atomic>
#include <mutex>
#include <limits>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dart {
namespace cli {

// DamageSubstitution + DAMAGE_SUBS + NUM_DAMAGE_SUBS now live in the shared header
// dart/damage_subs.hpp (single source; also consumed by the damage-aware submat emitter).
// They are in namespace dart::cli, matching this translation unit.

// O(1) lookup table indexed by [target_aa][query_aa]
// Stores index+1 into DAMAGE_SUBS (0 = no match)
struct DamageLookup {
    std::array<std::array<uint8_t, 128>, 128> table;

    constexpr DamageLookup() : table{} {
        for (size_t i = 0; i < NUM_DAMAGE_SUBS; ++i) {
            auto t = static_cast<unsigned char>(DAMAGE_SUBS[i].target_aa);
            auto q = static_cast<unsigned char>(DAMAGE_SUBS[i].query_aa);
            if (t < 128 && q < 128) {
                table[t][q] = static_cast<uint8_t>(i + 1);
            }
        }
    }
};

static constexpr DamageLookup DAMAGE_LOOKUP{};

// For opportunity counting in Bayesian scoring: which target AAs can be produced
// by C->T or G->A damage in at least one codon context.
struct DamageTargetClassLookup {
    std::array<uint8_t, 128> class_mask;

    constexpr DamageTargetClassLookup() : class_mask{} {
        for (size_t i = 0; i < NUM_DAMAGE_SUBS; ++i) {
            auto t = static_cast<unsigned char>(DAMAGE_SUBS[i].target_aa);
            if (t >= 128) continue;
            if (DAMAGE_SUBS[i].damage_class == 'C') class_mask[t] |= 0x1;
            if (DAMAGE_SUBS[i].damage_class == 'G') class_mask[t] |= 0x2;
        }
    }
};

static constexpr DamageTargetClassLookup DAMAGE_TARGET_CLASSES{};

static const DamageSubstitution* find_damage_sub(char target_aa, char query_aa) {
    auto t = static_cast<unsigned char>(target_aa);
    auto q = static_cast<unsigned char>(query_aa);
    if (t >= 128 || q >= 128) return nullptr;
    uint8_t idx = DAMAGE_LOOKUP.table[t][q];
    return idx ? &DAMAGE_SUBS[idx - 1] : nullptr;
}

struct DamageSite {
    size_t query_pos;
    size_t target_pos;
    char target_aa;
    char query_aa;
    char damage_class;     // 'C' or 'G'
    bool high_confidence;
    size_t dist_5prime;    // nt distance from 5' end
    size_t dist_3prime;    // nt distance from 3' end
    float p_damage;        // positional damage probability
};

// One best-hit read x target alignment staged for the .dcx --dcx emit (one row
// per read = its argmax target; em_keep carried as a SOFT weight, not a filter,
// so em_keep=0 rows are emitted too — the consumer weights, never WHEREs).
// Grain A scalars + the concatenated per-covered-position Grain B columns,
// built by walking qaln/taln inside the annotate scan (the only place the aln
// strings are live). posterior is patched in after the Bayesian output loop
// computes it; the record is then sorted by (target_ord, tstart) and appended
// to the dcx::Batch. Coverage policy: emit ALL covered target positions (both
// match and mismatch), presence bit set per covered column, p_damage / dist /
// damage_class taken from the read's DamageSite at mismatch positions and left
// at 0/None elsewhere — faithful to the format (bit k => covered) and lets
// convergence compute depth while filtering mismatches itself (ref_aa<>obs_aa).
using DcxReadRec = dcx::Record;

struct ProteinDamageSummary {
    uint32_t read_idx = 0;
    // query_id and target_id omitted — reconstructed at output time via reader.read_name/ref_name
    float evalue;
    float bits;
    float fident;
    size_t alnlen;
    size_t total_mismatches;
    size_t damage_consistent;   // Terminal damage sites (p_damage >= 5%)
    size_t non_damage_mismatches;
    size_t ct_sites;
    size_t ga_sites;
    size_t high_conf_sites;
    float damage_fraction;
    float max_p_damage;         // Highest p_damage among terminal sites
    float p_damaged;            // Overall probability read is damaged
    std::vector<DamageSite> sites;  // Cleared+shrunk in callback after aggregation
    // Channel D: Oxidative damage tracking (G→T at DNA level)
    // Gly→Cys (GGN→TGN) and Gly→Trp indicate G→T at codon position 1
    uint32_t gly_to_cys = 0;    // Gly→Cys substitutions (oxidation signature)
    uint32_t cys_to_gly = 0;    // Cys→Gly substitutions (reverse for ratio)
    uint32_t gly_to_trp = 0;    // Gly→Trp substitutions (strong oxidation signal)
    uint32_t trp_to_gly = 0;    // Trp→Gly substitutions (reverse)
    uint32_t gly_total = 0;     // Total Gly positions in target (denominator)
    // Pre-aggregated site scalars (populated from sites before clearing)
    float info_sum = 0.0f;      // sum(site.p_damage) → ReadDamageObs.info
    uint32_t terminal_5 = 0;    // C-terminal 5' site count → pa.terminal_5
    uint32_t terminal_3 = 0;    // G-terminal 3' site count → pa.terminal_3
    // Posterior stored after score_params setup for ReadRefEntry building
    float posterior = 0.0f;
    uint8_t classification = 0; // dart::AncientClassification cast to uint8_t
    // Synonymous damage (from damage index, invisible at protein level)
    size_t syn_5prime = 0;      // Synonymous C→T at 5' terminus
    size_t syn_3prime = 0;      // Synonymous G→A at 3' terminus
    bool has_syn_data = false;  // True if damage index was used
    // Per-read p_damaged from DART predict (stored in AGD index)
    float p_read = 0.0f;        // Per-read damage probability from predict
    bool has_p_read = false;    // True if p_read was retrieved from index
    // AA-level Bayesian evidence terms (derived from alignment + site calls)
    uint32_t aa_m_opportunities = 0;
    uint32_t aa_k_hits = 0;
    float aa_sum_qA = 0.0f;
    float aa_sum_log_qA_hits = 0.0f;
    float aa_sum_exp_decay = 0.0f;
    size_t qstart_0;
    size_t qlen = 0;  // Total read length (for Bayesian scoring)
    size_t tstart_0 = 0;  // 0-based target start position
    size_t tlen = 0;      // Target (reference) length
    // Short protein filtering metrics
    float delta_bits = 0.0f;    // bits_best - bits_second
    float bpa = 0.0f;           // bits / alnlen (bits-per-aligned-AA)
    uint8_t length_bin = 0;     // 0=15-24, 1=25-39, 2=40-59, 3=60+
    float z_bpa = 0.0f;         // Length-normalized BPA z-score
};

// Compact per-read data extracted from ProteinDamageSummary for protein_agg construction.
// Uses pre-aggregated site scalars to avoid per-read DamageSite vector heap allocations.
struct ProteinAggInput {
    uint32_t read_idx;
    uint32_t ref_idx;
    uint16_t tlen;
    float p_damaged;            // annotate reference-BASED posterior (diagnostic only; NOT auth)
    float p_read = 0.0f;        // AGD reference-FREE per-read posterior (drives authentication)
    bool has_p_read = false;    // true if p_read came from the AGD index
    uint32_t damage_consistent;
    float info_sum;       // sum(site.p_damage) → ReadDamageObs.info
    uint32_t terminal_5;  // C-terminal 5' site count → pa.terminal_5
    uint32_t terminal_3;  // G-terminal 3' site count → pa.terminal_3
    float posterior;      // Bayesian posterior (p_read + AA evidence); 0 if unavailable
    uint32_t gly_to_cys;  // Gly→Cys substitutions (oxidation signal)
    uint32_t gly_total;   // Total Gly positions in target (denominator)
    uint32_t aa_m_opportunities;  // #2: C→T/G→A eligible opportunities → ProteinReadObs.n_eligible
    uint32_t aa_k_hits;           // #2: damage-consistent mismatches → ProteinReadObs.k_damaged
    float gamma = 1.0f;         // real EM mapping responsibility (weights per-read auth evidence)
    bool em_ambiguous = false;  // EM confidence gate: exclude from damage/auth evidence
};

// Compact per-read assignment data sorted by ref_idx for streaming gene/protein output.
// Replaces the gene_agg and protein_agg unordered_maps with a sort-then-stream pass.
struct ReadRefEntry {
    uint32_t ref_idx;
    uint32_t qstart_0;
    uint32_t tstart_0;
    float gamma;
    float fident;
    float bits;
    float bits_second;
    uint32_t aln_len_on_target;
    uint16_t tlen;
    uint8_t length_bin;
    bool has_damage_obs;
    bool auth_keep;
    float posterior;
    uint32_t aa_m_opportunities;  // BUG1: C→T/G→A opportunity count (carried for clean reads too); uint32 avoids truncation
    uint8_t classification;  // dart::AncientClassification cast to uint8_t
    uint8_t em_ambiguous;    // BUG1: EM-ambiguous → exclude clean read from per-opportunity evidence
};  // packed per-read annotation row

// --- Binomial provenance-LR helpers (per-protein ancient authentication) ---
// dart authenticates whether a target protein's supporting reads carry a
// genuinely ancient damage signal via a binomial likelihood-ratio test:
//   H1: damaged-read rate = pi  (the library LSD-joint pi_damaged)
//   H0: damaged-read rate = eps0 (modern/background damaged-read rate)
// The exact binomial upper tail P(K>=k | N, p) is I_p(k, N-k+1); the
// regularized incomplete beta is evaluated by Lentz's continued fraction (NR 6.4).
namespace {
double dart_betacf(double a, double b, double x) {
    const int MAXIT = 300; const double EPS = 3.0e-12; const double FPMIN = 1.0e-300;
    const double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap;
    if (std::fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d; double h = d;
    for (int m = 1; m <= MAXIT; ++m) {
        const double m2 = 2.0 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if (std::fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c; if (std::fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if (std::fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c; if (std::fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; const double del = d * c; h *= del;
        if (std::fabs(del - 1.0) < EPS) break;
    }
    return h;
}
double dart_betai(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double lbt = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b)
                     + a * std::log(x) + b * std::log(1.0 - x);
    const double bt = std::exp(lbt);
    return (x < (a + 1.0) / (a + b + 2.0))
        ? bt * dart_betacf(a, b, x) / a
        : 1.0 - bt * dart_betacf(b, a, 1.0 - x) / b;
}
// Exact binomial upper tail P(K >= k | N, p).
double binom_upper_tail(int k, int N, double p) {
    if (k <= 0) return 1.0;
    if (k > N)  return 0.0;
    return dart_betai(static_cast<double>(k), static_cast<double>(N - k + 1), p);
}
// Smallest N for which the one-sided level-`alpha` binomial test of H0:rate=eps0
// reaches >=`power` power against H1:rate=pi. Exact (reuses binom_upper_tail):
// find the critical count k_c(N) = min k with P(K>=k|N,eps0)<=alpha, then the
// smallest N whose power P(K>=k_c|N,pi)>=power. Returns nmax when pi<=eps0 (no
// ancient signal above background) or no N up to nmax suffices. Library-level
// constant (pi/eps0/alpha are uniform across the library) — computed once.
int min_n_for_power(double pi, double eps0, double alpha, double power, int nmax) {
    if (pi <= eps0) return nmax;
    for (int N = 1; N <= nmax; ++N) {
        int kc = N + 1;
        for (int k = 1; k <= N; ++k) {
            if (binom_upper_tail(k, N, eps0) <= alpha) { kc = k; break; }
        }
        if (kc <= N && binom_upper_tail(kc, N, pi) >= power) return N;
    }
    return nmax;
}
} // namespace

// Parse strand and frame from a DART ORF id: "readname_+_1" or "readname_+_1_3".
// Returns {strand, frame}, defaulting to {'+', 0} when the id is not an ORF id.
//
// This used to scan backwards from the end taking the last '0'-'2' digit as the
// frame. That answered correctly only while the frame was the final field: once
// an ORF ordinal followed it, the scan either absorbed the ordinal as the frame
// or bailed and returned the defaults, silently. parse_dart_orf_id anchors on
// the '_<strand>_' marker and reads forward, so field count cannot shift the
// answer -- it either parses or reports failure.
static std::pair<char, int> parse_strand_and_frame(const std::string& query_id) {
    dart::DartOrfId p = dart::parse_dart_orf_id(query_id);
    return {p.strand, p.frame};
}

// Legacy wrapper for code that only needs strand
static char parse_strand(const std::string& query_id) {
    return parse_strand_and_frame(query_id).first;
}

// Compute positional damage probability: d_max * exp(-lambda * dist)
static float positional_damage_prob(size_t dist_from_terminus, float d_max, float lambda) {
    if (d_max <= 0.0f) return 0.0f;
    return d_max * std::exp(-lambda * static_cast<float>(dist_from_terminus));
}

static inline float evalue_from_log10(float evalue_log10) {
    if (!std::isfinite(evalue_log10) || evalue_log10 <= -300.0f) return 0.0f;
    return std::pow(10.0f, evalue_log10);
}

// Annotate a single alignment
static ProteinDamageSummary annotate_alignment(
    const std::string& query_id,
    std::string_view qaln, std::string_view taln,
    size_t qstart, size_t tstart, size_t qlen,
    float evalue, float bits, float fident,
    float d_max, float lambda, float prior_ancient)
{
    constexpr float Q_BASELINE = 0.005f;
    // Fix 1: AA_SCALE derived from codon usage analysis
    // P(damage-consistent AA event | C→T nucleotide event) ≈ 0.30
    // This accounts for: codon structure, reading frame, which AA changes we count
    constexpr float AA_SCALE = 0.30f;
    constexpr float EPS = 1e-9f;

    char strand = parse_strand(query_id);

    std::vector<DamageSite> sites;
    size_t total_mismatches = 0;
    uint32_t aa_m_opportunities = 0;
    uint32_t aa_k_hits = 0;
    float aa_sum_qA = 0.0f;
    float aa_sum_log_qA_hits = 0.0f;

    float sum_exp_decay = 0.0f;  // Fix 2: data-driven w = sum_exp_decay / m

    // Channel D: Oxidative damage tracking (G→T at DNA level)
    // Gly→Cys (GGN→TGN) and Gly→Trp indicate G→T at codon position 1
    uint32_t gly_to_cys = 0, cys_to_gly = 0;
    uint32_t gly_to_trp = 0, trp_to_gly = 0;
    uint32_t gly_total = 0;

    size_t q_pos = qstart;  // 0-based
    size_t t_pos = tstart;

    for (size_t i = 0; i < qaln.size() && i < taln.size(); ++i) {
        char q_aa = qaln[i];
        char t_aa = taln[i];

        if (q_aa == '-') {
            t_pos++;
            continue;
        }
        if (t_aa == '-') {
            q_pos++;
            continue;
        }

        size_t dist_5 = 0;
        size_t dist_3 = 0;
        if (strand == '+') {
            dist_5 = q_pos * 3;
            dist_3 = (qlen > 0 && q_pos < qlen) ? (qlen - 1 - q_pos) * 3 : 0;
        } else {
            dist_5 = (qlen > 0 && q_pos < qlen) ? (qlen - 1 - q_pos) * 3 : 0;
            dist_3 = q_pos * 3;
        }

        const auto t_code = static_cast<unsigned char>(t_aa);
        const uint8_t class_mask = (t_code < 128) ? DAMAGE_TARGET_CLASSES.class_mask[t_code] : 0;
        if (class_mask != 0) {
            aa_m_opportunities++;

            // Position-weighted hazard for this opportunity
            const float decay_5 = (class_mask & 0x1) ? std::exp(-lambda * static_cast<float>(dist_5)) : 0.0f;
            const float decay_3 = (class_mask & 0x2) ? std::exp(-lambda * static_cast<float>(dist_3)) : 0.0f;
            const float decay = std::max(decay_5, decay_3);
            sum_exp_decay += decay;

            // Compute q1 = P(damage | ancient) using bounded combine
            // h = d_max * exp(-lambda * dist)
            // AA_SCALE is not applied here: position-weighted formula uses raw d_max
            // to compute P(damage at this position | ancient) directly.
            const float h = d_max * decay;
            const float q0 = Q_BASELINE;
            // Bounded combine: q1 = 1 - (1-q0)*(1-h) to handle large hazards
            const float q1 = 1.0f - (1.0f - q0) * (1.0f - h);
            const float q0_safe = std::clamp(q0, EPS, 1.0f - EPS);
            const float q1_safe = std::clamp(q1, EPS, 1.0f - EPS);

            // Check if this opportunity has a damage-consistent hit
            bool is_hit = false;
            if (q_aa != t_aa) {
                const DamageSubstitution* sub = find_damage_sub(t_aa, q_aa);
                if (sub) {
                    is_hit = true;
                }
            }

            aa_sum_qA += h;  // Legacy: sum of hazards
        }

        if (q_aa != t_aa) {
            total_mismatches++;

            const DamageSubstitution* sub = find_damage_sub(t_aa, q_aa);
            if (sub) {
                float p_dmg;
                if (sub->damage_class == 'C') {
                    p_dmg = positional_damage_prob(dist_5, d_max, lambda);
                } else {
                    p_dmg = positional_damage_prob(dist_3, d_max, lambda);
                }
                aa_k_hits++;
                aa_sum_log_qA_hits += std::log((p_dmg + Q_BASELINE) / Q_BASELINE);

                sites.push_back({
                    q_pos, t_pos,
                    t_aa, q_aa,
                    sub->damage_class,
                    sub->high_confidence,
                    dist_5, dist_3,
                    p_dmg
                });
            }
        }

        // Channel D: Track Gly substitutions for oxidative damage detection
        // G→T at DNA codon pos 1: GGN→TGN = Gly→Cys, GGN→TGG = Gly→Trp
        if (t_aa == 'G') {
            gly_total++;
            if (q_aa == 'C') gly_to_cys++;
            else if (q_aa == 'W') gly_to_trp++;
        } else if (t_aa == 'C' && q_aa == 'G') {
            cys_to_gly++;
        } else if (t_aa == 'W' && q_aa == 'G') {
            trp_to_gly++;
        }

        q_pos++;
        t_pos++;
    }

    // Count ALL damage-consistent sites (not just terminal)
    // Use positional weighting for probability, but don't filter sites out
    // This improves recall while maintaining precision via scoring

    // Bayesian LLR approach: logit(posterior) = logit(prior) + Σ LLR_i
    // LLR_i = log(P(damage-consistent mismatch | damaged) / P(... | not damaged))
    //
    // For a damage-consistent mismatch at position with positional damage rate p_damage:
    // - P(mismatch | damaged read) = p_damage (expected from damage model)
    // - P(mismatch | not damaged) = seq_error_rate (~0.01 for Illumina)
    // LLR = log(p_damage / seq_error_rate) — but capped to avoid infinities
    //
    // For interior sites (p_damage ≈ 0), use a small baseline rate reflecting
    // that some interior damage still occurs.
    constexpr float SEQ_ERROR_RATE = 0.01f;        // P(damage-like mismatch | not damaged)
    constexpr float MIN_DAMAGE_RATE = 0.01f;       // Interior sites: LLR=0 (no evidence)
    // Empirical prior from the Pass-1 EM ancient fraction (pi) instead of a fixed 0.5;
    // a fixed 0.5 mis-scales p_damaged across libraries of different damage levels.
    const float PRIOR_DAMAGED = std::clamp(prior_ancient, 0.02f, 0.98f);

    size_t ct = 0, ga = 0, hc = 0;
    float max_p_damage = 0.0f;
    float llr_sum = 0.0f;   // Log-likelihood ratio accumulator

    for (const auto& s : sites) {
        if (s.damage_class == 'C') ct++;
        else ga++;
        if (s.high_confidence) hc++;
        if (s.p_damage > max_p_damage) max_p_damage = s.p_damage;

        // Bayesian LLR: each damage-consistent site provides evidence
        // P(site | damaged) = max(p_damage, MIN_DAMAGE_RATE) - position-dependent
        // P(site | not damaged) = SEQ_ERROR_RATE - sequencing error baseline
        float p_site_given_damaged = std::max(s.p_damage, MIN_DAMAGE_RATE);
        float p_site_given_undamaged = SEQ_ERROR_RATE;

        // LLR for this site (capped to avoid extreme values)
        float llr = std::log(p_site_given_damaged / p_site_given_undamaged);
        llr = std::clamp(llr, -3.0f, 5.0f);  // Cap at ~0.05 to ~150 odds ratio
        llr_sum += llr;
    }

    // Bayesian posterior: logit(posterior) = logit(prior) + Σ LLR
    // With no sites: return prior (no evidence to update)
    float p_damaged = PRIOR_DAMAGED;
    if (!sites.empty()) {
        float logit_prior = std::log(PRIOR_DAMAGED / (1.0f - PRIOR_DAMAGED));
        float logit_posterior = logit_prior + llr_sum;
        // Clamp to avoid numerical issues
        logit_posterior = std::clamp(logit_posterior, -10.0f, 10.0f);
        p_damaged = 1.0f / (1.0f + std::exp(-logit_posterior));
    } else {
        p_damaged = PRIOR_DAMAGED;  // No AA-level evidence → return prior
    }

    // damage_consistent now includes ALL damage-consistent sites (not just terminal)
    // terminal_damage is the subset with high p_damage
    size_t all_damage_consistent = sites.size();

    ProteinDamageSummary out{};
    out.read_idx = 0;
    out.evalue = evalue;
    out.bits = bits;
    out.fident = fident;
    out.alnlen = qaln.size();
    out.total_mismatches = total_mismatches;
    out.damage_consistent = all_damage_consistent;
    out.non_damage_mismatches = total_mismatches - all_damage_consistent;
    out.ct_sites = ct;
    out.ga_sites = ga;
    out.high_conf_sites = hc;
    out.damage_fraction = total_mismatches > 0
        ? static_cast<float>(all_damage_consistent) / static_cast<float>(total_mismatches)
        : 0.0f;
    out.max_p_damage = max_p_damage;
    out.p_damaged = p_damaged;
    out.sites = std::move(sites);
    out.aa_m_opportunities = aa_m_opportunities;
    out.aa_k_hits = aa_k_hits;
    out.aa_sum_qA = aa_sum_qA;
    out.aa_sum_log_qA_hits = aa_sum_log_qA_hits;
    out.aa_sum_exp_decay = sum_exp_decay;
    // Channel D: Oxidative damage counts
    out.gly_to_cys = gly_to_cys;
    out.cys_to_gly = cys_to_gly;
    out.gly_to_trp = gly_to_trp;
    out.trp_to_gly = trp_to_gly;
    out.gly_total = gly_total;
    return out;
}

// Expected terminal/middle coverage ratio for a uniformly-sampled gene.
//
// Short reads (length R) systematically under-cover protein termini: position p
// can only be reached by reads starting in [max(0, p-R+1), min(p, L-R)], so
// terminal positions have fewer covering read starts than interior positions.
// This function computes the theoretical ratio for a real gene so that
// terminal_middle_ratio() can be normalised to ~1.0 regardless of R or L.
// Results are cached because (L, R) repeat frequently across proteins.
static float expected_terminal_ratio(uint32_t L, uint32_t R) {
    if (L < 4 || R == 0) return 1.0f;
    R = std::min(R, L);

    static std::unordered_map<uint64_t, float> cache;
    const uint64_t key = (static_cast<uint64_t>(L) << 16) | R;
    const auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    const uint32_t edge = std::max(1u, std::min(L / 10u, L / 2u));
    double term_sum = 0.0, mid_sum = 0.0;
    uint32_t term_n = 0, mid_n = 0;
    for (uint32_t p = 0; p < L; ++p) {
        const uint32_t s_min = (p + 1 >= R) ? (p + 1 - R) : 0u;
        const uint32_t s_max = (L >= R) ? std::min(p, L - R) : 0u;
        const uint32_t cnt = (s_max >= s_min) ? (s_max - s_min + 1) : 0u;
        if (p < edge || p >= L - edge) { term_sum += cnt; ++term_n; }
        else                           { mid_sum  += cnt; ++mid_n;  }
    }
    const float ratio = (mid_n > 0 && mid_sum > 0.0)
        ? static_cast<float>((term_sum / term_n) / (mid_sum / mid_n))
        : 1.0f;
    cache.emplace(key, ratio);
    return ratio;
}

// Per-protein aggregated damage summary
struct ProteinAggregate {
    std::string protein_id;
    uint32_t tlen = 0;            // Target length (AA)
    size_t n_reads = 0;           // Total assigned reads mapping to this protein
    double n_reads_eff = 0.0;     // EM-effective assigned reads (sum gamma)
    size_t n_damage_reads = 0;    // Assigned reads with mismatch evidence (abundance/coverage base)
    size_t n_auth_reads = 0;      // Reads with an AGD p_read (denominator for the columns below)
    size_t n_damaged = 0;         // AGD reads with p_read > 0.5
    size_t total_damage_sites = 0; // Sum of damage_consistent across reads
    size_t terminal_5 = 0;        // C→T sites at 5' end (first 3 codons)
    size_t terminal_3 = 0;        // G→A sites at 3' end (first 3 codons)
    float max_p_damaged = 0.0f;   // Maximum AGD p_read across reads
    float sum_p_damaged = 0.0f;   // Sum of AGD p_read (for averaging over n_auth_reads)

    // Per-protein authentication accumulators (Stage 2b). p_protein_damaged =
    // sigmoid(logit(prior) + Σ gamma_i·llr_i), llr_i = the PRIOR-FREE per-read damage LLR.
    // The prior is added ONCE at output (library level) — NOT once per read: at LOW_ABUNDANCE
    // read_ancient_posterior refuses, so the .agd p_read already IS sigmoid(damage LLR) with
    // no baked prior; logit(p_read) recovers that raw LLR. Subtracting logit(prior) per read
    // (the earlier bug) injected +logit(prior)·Σgamma, making the score climb with COVERAGE, not
    // damage, and never fall below the prior. A read at the baseline has llr≈0 → contributes 0
    // (spike at prior); an anti-damage read has llr<0 → pushes BELOW the prior.
    double auth_sum_g_llr = 0.0;   // Σ gamma_i · logit(p_read_i)  [= Σ gamma_i · llr_i]
    double auth_sum_g2    = 0.0;   // Σ gamma_i²  — effective-sample-size proxy for the band

    // Per-read observations for Beta-Binomial model
    std::vector<ReadDamageObs> read_obs;  // Collected for proper statistical inference

    // Short protein filtering aggregates
    float sum_delta_bits = 0.0f;
    float sum_bpa = 0.0f;
    std::array<uint32_t, 4> length_bin_counts{};  // per-bin read counts

    // Mapping-pattern aggregates for spurious-hit filtering
    std::unordered_set<uint32_t> unique_starts;
    std::unordered_set<uint64_t> unique_start_pairs;
    uint32_t min_start = UINT32_MAX;
    uint32_t max_start = 0;
    std::vector<float> coverage; // Per-position coverage depth (gamma-weighted)
    double left_cov = 0.0;   // Alignment overlap in left terminal window
    double right_cov = 0.0;  // Alignment overlap in right terminal window
    double middle_cov = 0.0; // Alignment overlap in middle region
    double sum_aln_len = 0.0;
    double sum_aln_len_sq = 0.0;

    void add_assignment(uint32_t qstart, uint32_t tstart, uint32_t aln_len_on_target,
                        float gamma) {
        ++n_reads;
        n_reads_eff += gamma;

        unique_starts.insert(tstart);
        unique_start_pairs.insert((static_cast<uint64_t>(qstart) << 32) |
                                  static_cast<uint64_t>(tstart));
        if (tstart < min_start) min_start = tstart;
        if (tstart > max_start) max_start = tstart;

        if (tlen > 0 && coverage.empty()) {
            coverage.assign(tlen, 0.0f);
        }
        if (!coverage.empty()) {
            const size_t start = std::min(static_cast<size_t>(tstart), coverage.size());
            const size_t end = std::min(start + static_cast<size_t>(aln_len_on_target),
                                        coverage.size());
            for (size_t p = start; p < end; ++p) {
                coverage[p] += gamma;
            }
        }
    }

    float breadth() const {
        if (coverage.empty()) return 0.0f;
        uint32_t covered = 0;
        for (float v : coverage) {
            if (v > 0.0f) ++covered;
        }
        return static_cast<float>(covered) / static_cast<float>(coverage.size());
    }

    float depth_mean() const {
        if (coverage.empty()) return 0.0f;
        double total = 0.0;
        for (float v : coverage) total += v;
        return static_cast<float>(total / static_cast<double>(coverage.size()));
    }

    size_t n_distinct_starts() const {
        return unique_start_pairs.size();
    }

    float start_dispersion() const {
        if (n_reads == 0) return 0.0f;
        return static_cast<float>(unique_start_pairs.size()) /
               static_cast<float>(n_reads);
    }

    float start_diversity() const {
        if (n_reads == 0 || tlen == 0) return 0.0f;
        const uint32_t max_possible = std::min<uint32_t>(static_cast<uint32_t>(n_reads), tlen);
        if (max_possible == 0) return 0.0f;
        return static_cast<float>(unique_starts.size()) / static_cast<float>(max_possible);
    }

    float start_span() const {
        if (tlen == 0 || min_start > max_start) return 0.0f;
        return static_cast<float>(max_start - min_start + 1) / static_cast<float>(tlen);
    }

    float positional_score() const {
        return std::sqrt(start_diversity() * start_span());
    }

    float terminal_middle_ratio() const {
        if (tlen == 0) return 0.0f;
        size_t edge = std::max<size_t>(1, static_cast<size_t>(tlen) / 10);
        edge = std::min(edge, static_cast<size_t>(tlen) / 2);
        if (edge == 0) return 0.0f;
        const size_t term_len = 2 * edge;
        const size_t mid_len = (static_cast<size_t>(tlen) > term_len)
            ? (static_cast<size_t>(tlen) - term_len)
            : 0;
        // Use min(left_mean, right_mean) instead of their average.
        // This catches one-sided domain hits where reads pile up on only one
        // terminus: an N-terminal domain hit has right_mean ≈ 0, so min ≈ 0
        // and rho_norm << 1, correctly flagging it.  For authentic full-length
        // hits, left_mean ≈ right_mean so min ≈ average and rho_norm ≈ 1.
        const double left_mean = left_cov / static_cast<double>(edge);
        const double right_mean = right_cov / static_cast<double>(edge);
        const double term_mean = std::min(left_mean, right_mean);
        if (mid_len == 0) {
            return term_mean > 0.0 ? std::numeric_limits<float>::infinity() : 0.0f;
        }
        const double mid_mean = middle_cov / static_cast<double>(mid_len);
        if (mid_mean <= 0.0) {
            return term_mean > 0.0 ? std::numeric_limits<float>::infinity() : 0.0f;
        }
        const float observed = static_cast<float>(term_mean / mid_mean);
        // Normalise by the expected ratio for a real gene given mean read length.
        // expected_terminal_ratio averages both edges; for symmetric authentic
        // genes, left_mean ≈ right_mean so min ≈ average and normalisation holds.
        if (n_reads > 0) {
            const uint32_t R = static_cast<uint32_t>(sum_aln_len / static_cast<double>(n_reads));
            const float expected = expected_terminal_ratio(tlen, R);
            if (expected > 0.01f) return observed / expected;
        }
        return observed;
    }

    float avg_aln_len() const {
        return n_reads > 0 ? static_cast<float>(sum_aln_len / static_cast<double>(n_reads)) : 0.0f;
    }

    float std_aln_len() const {
        if (n_reads <= 1) return 0.0f;
        const double mean = sum_aln_len / static_cast<double>(n_reads);
        const double var = std::max(0.0, (sum_aln_len_sq / static_cast<double>(n_reads)) - (mean * mean));
        return static_cast<float>(std::sqrt(var));
    }
};

// Wilson score confidence interval for binomial proportion
// Returns (lower, upper) bounds for the true proportion given observed successes/n
inline std::pair<float, float> wilson_ci(uint32_t successes, uint32_t n, float z = 1.96f) {
    if (n == 0) return {0.0f, 0.0f};
    float p = static_cast<float>(successes) / n;
    float z2 = z * z;
    float denom = 1.0f + z2 / n;
    float center = (p + z2 / (2.0f * n)) / denom;
    float margin = z * std::sqrt((p * (1.0f - p) + z2 / (4.0f * n)) / n) / denom;
    return {std::max(0.0f, center - margin), std::min(1.0f, center + margin)};
}

// Wilson score interval for effective (possibly fractional) counts.
inline std::pair<float, float> wilson_ci_effective(double successes, double n, float z = 1.96f) {
    if (n <= 0.0) return {0.0f, 0.0f};
    successes = std::clamp(successes, 0.0, n);
    const double p = successes / n;
    const double z2 = static_cast<double>(z) * static_cast<double>(z);
    const double denom = 1.0 + z2 / n;
    const double center = (p + z2 / (2.0 * n)) / denom;
    const double margin = (static_cast<double>(z) / denom) *
        std::sqrt((p * (1.0 - p) + z2 / (4.0 * n)) / n);
    return {
        static_cast<float>(std::max(0.0, center - margin)),
        static_cast<float>(std::min(1.0, center + margin))
    };
}

// In-place quantile via nth_element. q is in [0, 1].
inline float quantile_inplace(std::vector<float>& values, double q) {
    if (values.empty()) return 0.0f;
    const double qc = std::clamp(q, 0.0, 1.0);
    const size_t idx = static_cast<size_t>(qc * static_cast<double>(values.size() - 1));
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(idx), values.end());
    return values[idx];
}

// Per-gene accumulator for functional profiling
// Aggregates read-level damage annotations into gene-level summaries
// Supports both best-hit (gamma=1.0) and EM-weighted (gamma<1.0) modes
struct GeneSummaryAccumulator {
    std::string target_id;
    uint32_t tlen = 0;
    uint32_t n_reads = 0;             // Assigned read count
    double n_reads_eff = 0.0;         // Effective assigned reads (sum of gamma)
    uint32_t n_damage_reads = 0;      // Assigned reads with mismatch evidence
    double n_damage_reads_eff = 0.0;  // Effective mismatch-evidence reads (sum gamma)
    double n_ancient_conf = 0.0;      // Gamma-weighted ancient confident (damage-evidence reads)
    double n_ancient_likely = 0.0;    // Gamma-weighted ancient likely (damage-evidence reads)
    double n_undetermined = 0.0;      // Gamma-weighted undetermined (damage-evidence reads)
    double n_modern_conf = 0.0;       // Gamma-weighted modern confident (damage-evidence reads)
    double sum_posterior = 0.0;       // Gamma-weighted posterior over damage-evidence reads
    double sum_fident = 0.0;          // Gamma-weighted identity over assigned reads
    double sum_fident_sq = 0.0;       // Gamma-weighted identity^2 (for std)
    double sum_aln_len = 0.0;       // Sum of alignment lengths (for avg)
    double sum_aln_len_sq = 0.0;    // Sum of squared aln lengths (for std)
    std::vector<float> coverage;    // Per-position coverage depth (gamma-weighted)
    std::unordered_set<uint32_t> unique_starts;  // Distinct start positions
    uint32_t min_start = UINT32_MAX;
    uint32_t max_start = 0;

    // Add assigned read with EM weight (gamma). For best-hit mode, gamma=1.0.
    void add_assignment(float fident, size_t tstart, size_t aln_len_on_target,
                        float gamma = 1.0f) {
        ++n_reads;
        n_reads_eff += gamma;
        sum_fident    += gamma * fident;
        sum_fident_sq += gamma * fident * fident;
        sum_aln_len += gamma * aln_len_on_target;
        sum_aln_len_sq += gamma * aln_len_on_target * aln_len_on_target;

        // Track start position diversity
        uint32_t start32 = static_cast<uint32_t>(tstart);
        unique_starts.insert(start32);
        if (start32 < min_start) min_start = start32;
        if (start32 > max_start) max_start = start32;

        // Update per-position coverage with gamma weight
        if (tlen > 0 && coverage.empty()) {
            coverage.resize(tlen, 0.0f);
        }
        size_t end = std::min(tstart + aln_len_on_target, static_cast<size_t>(tlen));
        for (size_t p = tstart; p < end; ++p) {
            coverage[p] += gamma;  // Gamma-weighted coverage
        }
    }

    // Add damage-evidence observation for classification/posterior summaries.
    void add_damage_observation(AncientClassification cls, float posterior, float gamma = 1.0f) {
        ++n_damage_reads;
        n_damage_reads_eff += gamma;
        sum_posterior += gamma * posterior;

        switch (cls) {
            case AncientClassification::AncientConfident: n_ancient_conf += gamma; break;
            case AncientClassification::AncientLikely:    n_ancient_likely += gamma; break;
            case AncientClassification::Undetermined:     n_undetermined += gamma; break;
            case AncientClassification::ModernConfident:  n_modern_conf += gamma; break;
        }
    }

    float breadth() const {
        if (coverage.empty()) return 0.0f;
        uint32_t covered = 0;
        for (float v : coverage) {
            if (v > 0.0f) ++covered;
        }
        return static_cast<float>(covered) / static_cast<float>(coverage.size());
    }

    // Abundance estimate: gamma-weighted depth (sum of coverage / length)
    float depth_mean() const {
        if (coverage.empty()) return 0.0f;
        double total = 0.0;
        for (float v : coverage) total += v;
        return static_cast<float>(total / coverage.size());
    }

    // Coverage standard deviation
    float depth_std() const {
        if (coverage.empty()) return 0.0f;
        double mean = depth_mean();
        double sum_sq = 0.0;
        for (float v : coverage) {
            double diff = v - mean;
            sum_sq += diff * diff;
        }
        return static_cast<float>(std::sqrt(sum_sq / coverage.size()));
    }

    // Coverage evenness (coefficient of variation): std / mean
    // Lower = more even coverage, higher = more uneven
    float depth_evenness() const {
        float mean = depth_mean();
        return mean > 0 ? depth_std() / mean : 0.0f;
    }

    // Average alignment length
    float avg_aln_len() const {
        return n_reads_eff > 0 ? static_cast<float>(sum_aln_len / n_reads_eff) : 0.0f;
    }

    // Standard deviation of alignment length
    float std_aln_len() const {
        if (n_reads_eff <= 1) return 0.0f;
        double mean = sum_aln_len / n_reads_eff;
        double var = (sum_aln_len_sq / n_reads_eff) - (mean * mean);
        return static_cast<float>(std::sqrt(std::max(0.0, var)));
    }

    // Number of unique start positions
    uint32_t n_unique_starts() const {
        return static_cast<uint32_t>(unique_starts.size());
    }

    // Start diversity: n_unique_starts / min(n_reads, tlen)
    // High = reads start at different positions (authentic)
    // Low = all reads start at same position (spurious)
    float start_diversity() const {
        if (n_reads == 0 || tlen == 0) return 0.0f;
        uint32_t max_possible = std::min(n_reads, tlen);
        return static_cast<float>(unique_starts.size()) / static_cast<float>(max_possible);
    }

    // Start span: (max_start - min_start + 1) / tlen
    // Measures how spread out the start positions are across the gene
    float start_span() const {
        if (tlen == 0 || min_start > max_start) return 0.0f;
        return static_cast<float>(max_start - min_start + 1) / static_cast<float>(tlen);
    }

    // Positional score: sqrt(diversity * span) [0-1]
    // Combined metric for filtering spurious matches
    float positional_score() const {
        return std::sqrt(start_diversity() * start_span());
    }

    // Terminal/middle coverage ratio, normalised by the expected ratio for a
    // real gene given mean read alignment length.  Raw ratio depends on both
    // protein length and read length even for authentic genes (short reads
    // under-cover termini geometrically); normalisation makes the metric
    // comparable across proteins and datasets.  A value ~1.0 indicates
    // terminal coverage consistent with a real gene; <<1.0 indicates reads
    // avoiding both ends (typical interior domain / conserved motif hit).
    float terminal_middle_ratio() const {
        if (coverage.empty()) return 0.0f;
        const size_t n = coverage.size();
        size_t edge = std::max<size_t>(1, n / 10);
        edge = std::min(edge, n / 2);
        if (edge == 0) return 0.0f;

        double left_sum = 0.0, right_sum = 0.0;
        double mid_sum = 0.0;
        size_t mid_n = 0;
        for (size_t i = 0; i < n; ++i) {
            if (i < edge) {
                left_sum += coverage[i];
            } else if (i >= n - edge) {
                right_sum += coverage[i];
            } else {
                mid_sum += coverage[i];
                ++mid_n;
            }
        }

        // Use min(left_mean, right_mean) so that one-sided domain hits
        // (where only one terminus is covered) are correctly flagged.
        const double left_mean = left_sum / static_cast<double>(edge);
        const double right_mean = right_sum / static_cast<double>(edge);
        const double term_mean = std::min(left_mean, right_mean);
        if (mid_n == 0) {
            return term_mean > 0.0 ? std::numeric_limits<float>::infinity() : 0.0f;
        }
        const double mid_mean = mid_sum / static_cast<double>(mid_n);
        if (mid_mean <= 0.0) {
            return term_mean > 0.0 ? std::numeric_limits<float>::infinity() : 0.0f;
        }
        const float observed = static_cast<float>(term_mean / mid_mean);
        if (n_reads > 0) {
            const uint32_t R = static_cast<uint32_t>(sum_aln_len / static_cast<double>(n_reads));
            const float expected = expected_terminal_ratio(static_cast<uint32_t>(n), R);
            if (expected > 0.01f) return observed / expected;
        }
        return observed;
    }

    // Effective read count (sum of gamma weights)
    float eff_reads() const {
        return static_cast<float>(n_reads_eff);
    }

    float avg_posterior() const {
        return n_damage_reads_eff > 0 ? static_cast<float>(sum_posterior / n_damage_reads_eff) : 0.0f;
    }

    float avg_fident() const {
        return n_reads_eff > 0 ? static_cast<float>(sum_fident / n_reads_eff) : 0.0f;
    }

    float std_fident() const {
        if (n_reads_eff <= 1) return 0.0f;
        double mean = sum_fident / n_reads_eff;
        double var = (sum_fident_sq / n_reads_eff) - (mean * mean);
        return static_cast<float>(std::sqrt(std::max(0.0, var)));
    }

    // Ancient fraction based on effective counts
    float ancient_frac() const {
        double total = n_ancient_conf + n_ancient_likely + n_modern_conf;
        return total > 0 ? static_cast<float>((n_ancient_conf + n_ancient_likely) / total) : 0.0f;
    }
};

// Functional profiling accumulator
struct FunctionalAccumulator {
    std::string function_id;
    std::string db_type;  // e.g. "MAP" (generic mapping source)
    uint32_t n_genes = 0;
    uint32_t n_damaged_genes = 0;  // Genes with mean posterior >= threshold
    uint32_t n_reads = 0;
    float n_ancient = 0.0f;  // Can be fractional from EM
    float n_modern = 0.0f;
    float n_undetermined = 0.0f;
    double sum_posterior = 0.0;
    // Coverage stats (per-gene depth_mean aggregated across genes)
    double sum_coverage = 0.0;
    double sum_coverage_sq = 0.0;
    // Identity and read-length stats (reads-weighted across all reads in group)
    double sum_fident = 0.0;
    double sum_fident_sq = 0.0;
    double sum_aln_len = 0.0;
    double sum_aln_len_sq = 0.0;
    double sum_eff_reads = 0.0;  // denominator for reads-weighted stats

    void add_gene(uint32_t gene_reads, float gene_ancient, float gene_modern,
                  float gene_undetermined, float gene_mean_posterior,
                  bool gene_is_damaged,
                  float gene_depth_mean,
                  double gene_sum_fident, double gene_sum_fident_sq,
                  double gene_sum_aln_len, double gene_sum_aln_len_sq,
                  double gene_eff_reads) {
        ++n_genes;
        if (gene_is_damaged) ++n_damaged_genes;
        n_reads += gene_reads;
        n_ancient += gene_ancient;
        n_modern += gene_modern;
        n_undetermined += gene_undetermined;
        sum_posterior += gene_mean_posterior * gene_reads;
        sum_coverage    += gene_depth_mean;
        sum_coverage_sq += static_cast<double>(gene_depth_mean) * gene_depth_mean;
        sum_fident      += gene_sum_fident;
        sum_fident_sq   += gene_sum_fident_sq;
        sum_aln_len     += gene_sum_aln_len;
        sum_aln_len_sq  += gene_sum_aln_len_sq;
        sum_eff_reads   += gene_eff_reads;
    }

    float ancient_frac() const {
        float total = n_ancient + n_modern;
        return total > 0 ? n_ancient / total : 0.0f;
    }

    float mean_posterior() const {
        return n_reads > 0 ? static_cast<float>(sum_posterior / n_reads) : 0.0f;
    }

    float damaged_gene_frac() const {
        return n_genes > 0 ? static_cast<float>(n_damaged_genes) / static_cast<float>(n_genes) : 0.0f;
    }

    float coverage_mean() const {
        return n_genes > 0 ? static_cast<float>(sum_coverage / n_genes) : 0.0f;
    }

    float coverage_stdev() const {
        if (n_genes < 2) return 0.0f;
        double mean = sum_coverage / n_genes;
        double var = (sum_coverage_sq / n_genes) - (mean * mean);
        return static_cast<float>(std::sqrt(std::max(0.0, var)));
    }

    float avg_identity() const {
        return sum_eff_reads > 0 ? static_cast<float>(sum_fident / sum_eff_reads) : 0.0f;
    }

    float stdev_identity() const {
        if (sum_eff_reads <= 1) return 0.0f;
        double mean = sum_fident / sum_eff_reads;
        double var = (sum_fident_sq / sum_eff_reads) - (mean * mean);
        return static_cast<float>(std::sqrt(std::max(0.0, var)));
    }

    float avg_read_length() const {
        return sum_eff_reads > 0 ? static_cast<float>(sum_aln_len / sum_eff_reads) : 0.0f;
    }

    float stdev_read_length() const {
        if (sum_eff_reads <= 1) return 0.0f;
        double mean = sum_aln_len / sum_eff_reads;
        double var = (sum_aln_len_sq / sum_eff_reads) - (mean * mean);
        return static_cast<float>(std::sqrt(std::max(0.0, var)));
    }
};

// Load a 2-column mapping file: gene_id<TAB>function_id
static std::unordered_map<std::string, std::string> load_mapping_file(
    const std::string& path, bool verbose = false)
{
    std::unordered_map<std::string, std::string> mapping;
    std::ifstream in(path);
    if (!in.is_open()) {
        if (verbose) {
            std::cerr << "Warning: Cannot open mapping file: " << path << "\n";
        }
        return mapping;
    }

    std::string line;
    // Skip header if present
    if (std::getline(in, line)) {
        // Check if it's a header (starts with "gene" or "#")
        if (line.empty() || line[0] == '#' ||
            line.find("gene") == 0 || line.find("Gene") == 0) {
            // It's a header, skip it
        } else {
            // Not a header, process it
            size_t tab = line.find('\t');
            if (tab != std::string::npos) {
                mapping[line.substr(0, tab)] = line.substr(tab + 1);
            }
        }
    }

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t tab = line.find('\t');
        if (tab != std::string::npos) {
            mapping[line.substr(0, tab)] = line.substr(tab + 1);
        }
    }

    if (verbose) {
        std::cerr << "Loaded " << mapping.size() << " mappings from " << path << "\n";
    }
    return mapping;
}

// Estimate d_max from alignment mismatches (first pass)
// Counts damage-consistent substitutions at terminal vs interior positions
struct DamageEstimate {
    size_t terminal_damage = 0;   // Damage-consistent subs within 9nt of terminus
    size_t terminal_total = 0;    // All subs within 9nt of terminus
    size_t interior_damage = 0;   // Damage-consistent subs in interior
    size_t interior_total = 0;    // All subs in interior

    float estimate_d_max() const {
        if (terminal_total == 0) return 0.0f;
        float term_rate = static_cast<float>(terminal_damage) / static_cast<float>(terminal_total);
        float int_rate = (interior_total > 0)
            ? static_cast<float>(interior_damage) / static_cast<float>(interior_total)
            : 0.0f;
        // d_max ≈ terminal excess rate / expected AA change rate from C→T
        // ~50% of C→T cause AA change, so multiply by 2
        float excess = std::max(0.0f, term_rate - int_rate);
        return std::min(0.5f, excess * 2.0f);  // Cap at 50%
    }
};

static void count_damage_for_estimate(
    const std::string& query_id, std::string_view qaln, std::string_view taln,
    size_t qstart, size_t qlen, DamageEstimate& est)
{
    char strand = parse_strand(query_id);
    size_t q_pos = qstart;

    for (size_t i = 0; i < qaln.size() && i < taln.size(); ++i) {
        char q_aa = qaln[i];
        char t_aa = taln[i];

        if (q_aa == '-') continue;
        if (t_aa == '-') { q_pos++; continue; }

        if (q_aa != t_aa) {
            size_t dist_5, dist_3;
            if (strand == '+') {
                dist_5 = q_pos * 3;
                dist_3 = (qlen > 0 && q_pos < qlen) ? (qlen - 1 - q_pos) * 3 : 0;
            } else {
                dist_5 = (qlen > 0 && q_pos < qlen) ? (qlen - 1 - q_pos) * 3 : 0;
                dist_3 = q_pos * 3;
            }
            size_t min_dist = std::min(dist_5, dist_3);

            bool is_terminal = (min_dist <= 9);  // First 3 codons
            const DamageSubstitution* sub = find_damage_sub(t_aa, q_aa);

            if (is_terminal) {
                est.terminal_total++;
                if (sub) est.terminal_damage++;
            } else {
                est.interior_total++;
                if (sub) est.interior_damage++;
            }
        }
        q_pos++;
    }
}

int cmd_damage_annotate(int argc, char* argv[]) {
    const auto run_start = std::chrono::steady_clock::now();
    std::string emi_file;
    std::string output_file;
    std::string sites_file;
    std::string dcx_file;              // .dcx columnar interchange (grain A/B for convergence)
    std::string library_label;         // per-library id stamped into the .dcx calib section
    std::string protein_summary_file;  // Per-protein aggregation
    std::string damage_index_file;     // Binary damage index from predict
    float d_max = -1.0f;  // -1 means "estimate from data"
    float lambda = 0.3f;
    bool lambda_set_by_user = false;
    bool refine_damage = false;
    std::string lib_type = "auto";    // ss, ds, or auto-detect
    const float threshold = 0.7f;           // Classification threshold for is_damaged
    bool verbose = false;

    // Gene summary options (work with or without --em)
    std::string gene_summary_file;
    float min_breadth = 0.10f;
    float min_depth = 0.5f;
    uint32_t min_reads = 3;
    uint32_t min_distinct_starts = 2;
    // Q-CLUST ship-gate: the pooled (K,N) binomial tail is anti-conservative when the eligible
    // sites concentrate in one molecule (measured: 62% of p>0 proteins carry >=90% of N on a
    // single read post-derep). Fail single-molecule-dominated / too-sparse proteins.
    uint32_t min_informative_reads = 3;
    float max_single_molecule_frac = 0.5f;
    float min_positional_score = 0.0f;
    float min_terminal_ratio = 0.0f;
    // Binomial provenance-LR authentication (PRIMARY per-protein ancient call).
    float modern_damage_rate = 0.02f;        // --modern-damage-rate/--eps0: H0 modern damaged-read rate
    float auth_alpha = 0.05f;                // --auth-alpha: level for the LR power gate + ancient call
    float per_read_damage_threshold = 0.5f;  // --per-read-damage-threshold: p_read cutoff for the per-read damaged call (k)
    bool user_set_terminal_ratio = false;  // True when --min-terminal-ratio is explicitly provided
    bool auto_calibrate_spurious = false;  // Auto-derive thresholds from data

    // EM reassignment options (enabled by default for abundance estimation)
    bool use_em = true;
    // Node-local scratch dir for the O_DIRECT EM tape (production streaming path).
    // Must be fast local NVMe supporting O_DIRECT (NOT /tmp, NOT NFS). Env fallback
    // DART_EM_TAPE_DIR; when unset the tape path errors and asks for it explicitly.
    std::string em_tape_dir;
    // EM fixed-point checkpoint (OPT-G). Iterations only; the final tape pass always re-runs.
    std::string em_checkpoint_path;  // --em-checkpoint; default <em_tape_dir>/dart_emtape_<pid>.ckpt
    bool em_resume = false;          // --resume: load the checkpoint when its fingerprint matches
    const uint32_t em_max_iters = 100;
    // Raw-bit temperature for the legacy (non-calibrated) real component. mmseqs
    // bit scores are log2-odds; the principled bits→nats factor is ln2≈0.693.
    // The old 3.0 over-tempered ~4.3x and forced homolog near-ties to gamma≈1.
    const double em_lambda_b = M_LN2;
    const double em_tol = 1e-4;
    double em_min_prob = 1e-6;       // absolute posterior threshold for keeping alignment
    const double em_min_margin = 0.5;      // auth evidence requires gamma_top1 - gamma_top2
    const double em_max_neff = 2.0;        // auth evidence requires exp(H(gamma)) <= cap
    // Fail-loud: a flat targeted DB can force-map absent-protein reads onto their
    // nearest wrong homolog. EM without a null model is refused unless the caller
    // opts in with --no-null (accept the risk explicitly).
    bool no_null = false;
    // Damage-free mapping (removes circularity: damage likelihoods must not both
    // pick the homolog and then be reused in gamma_ancient). OFF => mapping gamma
    // is damage-free; gamma_ancient is emitted post-hoc against frozen weights when the
    // damage index is pi DETECTED (see agd_read_prior).
    bool em_use_damage_mapping = false;
    // Soft multi-mapping emit (DCX 2.2): emit the top-K EM candidates per read (not just the
    // argmax) so downstream can distribute a read's SAAV evidence across homologs by gamma.
    // OFF => byte-identical one-argmax-row-per-read. See convergence/DCX_2.2_MULTIMAP.md.
    bool   emit_multimap = false;
    int    multimap_k = 5;           // keep top-K candidates per read (ranked by gamma)
    double multimap_floor = 0.05;    // AND gamma > floor

    // Null/noise component (spurious short-peptide rejection). Off by default.
    std::string em_null_model_file;  // .null sidecar (or deprecated TSV); default <emi>.null when present
    double em_null_tau = 0.5;        // calibrated path: reject read when gamma_null > tau
    double null_p_max = 0.01;        // uncalibrated path: reject read when null_p > this
    bool null_p_max_set = false;

    // Coherent-LR calibration (f_true + f_0). When loaded, the real EM component
    // becomes the normalized true-hit density log f_true(bit|Lbin) and the null
    // f_0 lives on the raw-bit scale — a coherent likelihood ratio (lambda≡1).
    std::string em_calib_model_file;      // load a frozen calibration model
    const double em_calib_temperature = 1.0;    // frozen misspecification correction
    // Posterior-identifiability abstention (sol): the 3-way call is on the per-read
    // density log-ratio LLR = log f_true - log f_0 of the best EM candidate, NOT on
    // gamma_null (which saturates to 1 via the EM w_null dynamics even for reads whose
    // evidence cannot separate true from null — that mislabels information-poor short
    // reads as contamination). LLR maps to a posterior/FDR floor:
    //   LLR >= keep_llr   → keep (assign);  LLR <= reject_llr → reject (spurious);
    //   between           → abstain (insufficient_information): excluded from
    //   assignment/auth but NOT counted as a contamination rejection.
    const double em_keep_llr = 0.5;      // f_true >= e^0.5 (~1.65x) as likely as f_0 → keep
    const double em_reject_llr = -5.0;   // f_0 >= e^5 (~150x) more likely than f_true → reject

    // Coverage-aware outer EM options

    // Bayesian scoring parameters
    float prior_ancient = 0.10f;     // Prior P(ancient) for Bayesian scorer
    bool  auto_prior_ancient = false; // Auto-calibrate prior from mean p_read distribution
    int   min_damage_sites = -1;     // Min susceptible positions to trust site evidence (-1 = auto)
    float ancient_threshold  = 0.60f; // --ancient-threshold
    float modern_threshold   = 0.25f; // --modern-threshold
    float terminal_threshold = 0.50f; // --terminal-threshold
    float site_cap           = 3.0f;  // --site-cap
    float w0                 = 0.3f;  // --w0
    // Coverage EM tuning

    // Alignment-level pre-filters (applied before any processing)
    const float aln_min_identity = 0.0f;   // 0 = no filter
    const float aln_min_bits = 0.0f;       // 0 = no filter
    const float aln_max_evalue = 1e10f;    // very large = no filter
    int threads = 0;                 // 0 = OpenMP runtime default

    // Functional profiling options
    std::string map_file;            // gene_id -> group mapping
    std::string functional_summary_file;  // Output: per-function stats
    std::string anvio_ko_file;       // Output: Anvi'o-compatible grouped abundance
    std::string annotation_source = "DART";

    for (int i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "--emi") == 0 || strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
            emi_file = argv[++i];
        } else if (strcmp(argv[i], "--d-max") == 0 && i + 1 < argc) {
            d_max = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda") == 0 && i + 1 < argc) {
            lambda = std::stof(argv[++i]);
            lambda_set_by_user = true;
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--sites") == 0 && i + 1 < argc) {
            sites_file = argv[++i];
        } else if (strcmp(argv[i], "--dcx") == 0 && i + 1 < argc) {
            dcx_file = argv[++i];
        } else if (strcmp(argv[i], "--emit-multimap") == 0) {
            emit_multimap = true;
        } else if (strcmp(argv[i], "--multimap-k") == 0 && i + 1 < argc) {
            multimap_k = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--multimap-floor") == 0 && i + 1 < argc) {
            multimap_floor = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--library-label") == 0 && i + 1 < argc) {
            library_label = argv[++i];
        } else if (strcmp(argv[i], "--protein-summary") == 0 && i + 1 < argc) {
            protein_summary_file = argv[++i];
        } else if ((strcmp(argv[i], "--library-type") == 0 || strcmp(argv[i], "--lib-type") == 0) && i + 1 < argc) {
            lib_type = argv[++i];
            if (lib_type != "ss" && lib_type != "ds" && lib_type != "auto") {
                std::cerr << "Error: --library-type must be ss, ds or auto, got '" << lib_type << "'\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--damage-index") == 0 && i + 1 < argc) {
            damage_index_file = argv[++i];
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--refine-damage") == 0) {
            refine_damage = true;
        } else if (strcmp(argv[i], "--em") == 0) {
            use_em = true;
        } else if (strcmp(argv[i], "--no-em") == 0) {
            use_em = false;
        } else if (strcmp(argv[i], "--em-min-prob") == 0 && i + 1 < argc) {
            em_min_prob = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--no-null") == 0) {
            no_null = true;
        } else if (strcmp(argv[i], "--em-use-damage-mapping") == 0) {
            em_use_damage_mapping = true;
        } else if (strcmp(argv[i], "--em-null-model") == 0 && i + 1 < argc) {
            em_null_model_file = argv[++i];
        } else if (strcmp(argv[i], "--em-calib-model") == 0 && i + 1 < argc) {
            em_calib_model_file = argv[++i];
        } else if (strcmp(argv[i], "--em-null-tau") == 0 && i + 1 < argc) {
            em_null_tau = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--null-p-max") == 0 && i + 1 < argc) {
            null_p_max = std::stod(argv[++i]);
            null_p_max_set = true;
        } else if (strcmp(argv[i], "--no-coverage-em") == 0) {
            // Kept for scripts: asks for the only behaviour there is.
        } else if (strcmp(argv[i], "--prior-ancient") == 0 && i + 1 < argc) {
            prior_ancient = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--min-damage-sites") == 0 && i + 1 < argc) {
            min_damage_sites = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--ancient-threshold") == 0 && i + 1 < argc) {
            ancient_threshold = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--modern-threshold") == 0 && i + 1 < argc) {
            modern_threshold = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--terminal-threshold") == 0 && i + 1 < argc) {
            terminal_threshold = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--site-cap") == 0 && i + 1 < argc) {
            site_cap = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--w0") == 0 && i + 1 < argc) {
            w0 = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--auto-prior-ancient") == 0) {
            auto_prior_ancient = true;
        } else if (strcmp(argv[i], "--em-streaming") == 0) {
            std::cerr << "Error: --em-streaming was removed; compact_em_map is the sole "
                         "EM solver (in-memory CSR). Drop --em-streaming.\n";
            return 1;
        } else if (strcmp(argv[i], "--em-max-memory") == 0 || strcmp(argv[i], "-m") == 0 ||
                   strcmp(argv[i], "--em-null-fdr") == 0 ||
                   strcmp(argv[i], "--coverage-em") == 0 ||
                   strncmp(argv[i], "--coverage-em-", 14) == 0 ||
                   strcmp(argv[i], "--depth-gate-lo") == 0 ||
                   strcmp(argv[i], "--depth-gate-hi") == 0 ||
                   strcmp(argv[i], "--w-min") == 0 ||
                   strcmp(argv[i], "--aln-min-identity") == 0 ||
                   strcmp(argv[i], "--aln-min-bits") == 0 ||
                   strcmp(argv[i], "--aln-max-evalue") == 0 ||
                   strcmp(argv[i], "--analysis-prefix") == 0 ||
                   strcmp(argv[i], "--analysis-proteins") == 0 ||
                   strcmp(argv[i], "--analysis-blast8") == 0 ||
                   strcmp(argv[i], "--analysis-categories") == 0 ||
                   strcmp(argv[i], "--blast8-unique") == 0 ||
                   strcmp(argv[i], "--combined-output") == 0 ||
                   strcmp(argv[i], "--corrected") == 0 ||
                   strcmp(argv[i], "--protein-filtered") == 0 ||
                   strcmp(argv[i], "--max-dist") == 0 ||
                   strcmp(argv[i], "--threshold") == 0 ||
                   strcmp(argv[i], "--use-identity") == 0 ||
                   strcmp(argv[i], "--no-identity") == 0 ||
                   strcmp(argv[i], "--em-iters") == 0 ||
                   strcmp(argv[i], "--em-lambda") == 0 ||
                   strcmp(argv[i], "--em-tol") == 0 ||
                   strcmp(argv[i], "--em-min-margin") == 0 ||
                   strcmp(argv[i], "--em-max-neff") == 0 ||
                   strcmp(argv[i], "--em-audit") == 0 ||
                   strcmp(argv[i], "--em-keep-llr") == 0 ||
                   strcmp(argv[i], "--em-reject-llr") == 0 ||
                   strcmp(argv[i], "--em-calib-temperature") == 0 ||
                   strcmp(argv[i], "--preset") == 0) {
            std::cerr << "Error: " << argv[i] << " was removed in DART 2.0. Drop it; runs without it "
                         "behave as before (see 'Removed in 2.0' in the Command Reference).\n";
            return 1;
        } else if (strcmp(argv[i], "--fit-null-from") == 0 || strcmp(argv[i], "--null-out") == 0) {
            std::cerr << "Error: " << argv[i] << " was removed in DART 2.0. Fit the null model with "
                         "'dart hits2emi --decoy-hits decoy_hits.tsv' or 'dart null'; damage-annotate "
                         "loads <emi>.null by itself.\n";
            return 1;
        } else if (strcmp(argv[i], "--em-tape-dir") == 0 && i + 1 < argc) {
            em_tape_dir = argv[++i];
        } else if (strcmp(argv[i], "--em-checkpoint") == 0 && i + 1 < argc) {
            em_checkpoint_path = argv[++i];
        } else if (strcmp(argv[i], "--resume") == 0) {
            em_resume = true;
        } else if (strcmp(argv[i], "--gene-summary") == 0 && i + 1 < argc) {
            gene_summary_file = argv[++i];
        } else if (strcmp(argv[i], "--min-breadth") == 0 && i + 1 < argc) {
            min_breadth = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--min-depth") == 0 && i + 1 < argc) {
            min_depth = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--min-reads") == 0 && i + 1 < argc) {
            min_reads = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--min-distinct-starts") == 0 && i + 1 < argc) {
            min_distinct_starts = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--min-informative-reads") == 0 && i + 1 < argc) {
            min_informative_reads = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--max-single-molecule-frac") == 0 && i + 1 < argc) {
            max_single_molecule_frac = std::stof(argv[++i]);
        } else if ((strcmp(argv[i], "--modern-damage-rate") == 0 || strcmp(argv[i], "--eps0") == 0) && i + 1 < argc) {
            modern_damage_rate = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--auth-alpha") == 0 && i + 1 < argc) {
            auth_alpha = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--per-read-damage-threshold") == 0 && i + 1 < argc) {
            per_read_damage_threshold = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--min-positional-score") == 0 && i + 1 < argc) {
            min_positional_score = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--min-terminal-ratio") == 0 && i + 1 < argc) {
            min_terminal_ratio = std::stof(argv[++i]);
            user_set_terminal_ratio = true;
        } else if (strcmp(argv[i], "--auto-calibrate-spurious") == 0) {
            auto_calibrate_spurious = true;
        } else if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0) && i + 1 < argc) {
            threads = std::stoi(argv[++i]);
            if (threads < 1) {
                std::cerr << "Error: --threads must be >= 1\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--map") == 0 && i + 1 < argc) {
            map_file = argv[++i];
        } else if (strcmp(argv[i], "--functional-summary") == 0 && i + 1 < argc) {
            functional_summary_file = argv[++i];
        } else if (strcmp(argv[i], "--anvio-ko") == 0 && i + 1 < argc) {
            anvio_ko_file = argv[++i];
        } else if (strcmp(argv[i], "--annotation-source") == 0 && i + 1 < argc) {
            annotation_source = argv[++i];
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0 ||
                   strcmp(argv[i], "--help-all") == 0) {
            std::cout <<
                "Usage: dart damage-annotate --emi <hits.emi> [options]\n"
                "\n"
                "Input:\n"
                "  --emi, -i FILE          EMI v7 alignment index\n"
                "  --damage-index FILE     AGD damage index from predict\n"
                "  --library-type TYPE     ss, ds or auto (default: auto)\n"
                "Null model (pick one):\n"
                "  --no-null               Run without a null model\n"
                "  --em-null-model FILE    Empirical decoy null model\n"
                "  --em-calib-model FILE   Calibrated true-hit and null model\n"
                "Output:\n"
                "  -o FILE                 Per-read TSV (default: stdout)\n"
                "  --protein-summary FILE  Per-protein damage TSV\n"
                "  --gene-summary FILE     Per-gene statistics TSV\n"
                "  --dcx FILE              Columnar alignment and site output (requires EM)\n"
                "  --library-label STR     Library ID in DCX (default: DCX basename)\n"
                "Running:\n"
                "  --em-tape-dir DIR       Node-local disk for EM tape; or DART_EM_TAPE_DIR\n"
                "  --threads, -t INT       OpenMP threads (default: runtime default)\n"
                "  -v                      Verbose output\n";
            if (strcmp(argv[i], "--help-all") == 0) {
                std::cout <<
                "\n"
                "Damage:\n"
                "  --sites FILE                  Per-site damage TSV\n"
                "  --refine-damage               Update d_max from EM assignments\n"
                "  --auto-prior-ancient          Estimate the ancient prior from AGD reads\n"
                "  --prior-ancient FLOAT         Ancient prior (default: 0.10)\n"
                "  --d-max RATE                  Damage rate (default: estimated from data)\n"
                "  --lambda FLOAT                Decay constant (default: 0.3)\n"
                "Multimapping:\n"
                "  --emit-multimap               Emit multiple EM candidates per read\n"
                "  --multimap-k INT              Maximum candidates per read (default: 5)\n"
                "  --multimap-floor FLOAT        Minimum candidate gamma (default: 0.05)\n"
                "Protein filters:\n"
                "  --auto-calibrate-spurious     Estimate mapping thresholds from supported proteins\n"
                "  --min-breadth FLOAT           Minimum breadth (default: 0.10)\n"
                "  --min-depth FLOAT             Minimum depth (default: 0.5)\n"
                "  --min-reads INT               Minimum reads (default: 3)\n"
                "Functional profiling:\n"
                "  --map FILE                    Gene-to-group TSV\n"
                "  --functional-summary FILE     Per-function TSV (requires --map)\n"
                "  --anvio-ko FILE               Gene-level group abundance (requires --map)\n"
                "  --annotation-source STR       Anvi'o source label (default: DART)\n"
                "Null model:\n"
                "  --em-null-model FILE          .null sidecar from hits2emi --decoy-hits or dart null\n"
                "                                (default: <emi>.null when it exists)\n"
                "  --null-p-max P                Reject a read when the chance that a decoy scores as\n"
                "                                high in its aligned-length bin exceeds P (default: 0.01)\n"
                "  --em-null-tau FLOAT           Null responsibility cutoff, calibrated null column only\n"
                "                                (--em-calib-model; default: 0.5)\n"
                "Resume:\n"
                "  --em-checkpoint FILE          EM checkpoint (default: tape dir/dart_emtape_<pid>.ckpt)\n"
                "  --resume                      Resume a matching EM checkpoint\n"
                "Model parameters (change only with a reason):\n"
                "  --ancient-threshold FLOAT     Ancient posterior cutoff (default: 0.60)\n"
                "  --modern-threshold FLOAT      Modern posterior cutoff (default: 0.25)\n"
                "  --terminal-threshold FLOAT    Terminal read posterior cutoff (default: 0.50)\n"
                "  --auth-alpha FLOAT            Protein authentication level (default: 0.05)\n"
                "  --modern-damage-rate, --eps0 FLOAT  Modern damaged-read rate (default: 0.02)\n"
                "  --per-read-damage-threshold FLOAT   Damaged-read cutoff (default: 0.50)\n"
                "  --min-damage-sites INT        Minimum susceptible sites (default: auto)\n"
                "  --site-cap FLOAT              Site log-Bayes-factor cap (default: 3.0)\n"
                "  --w0 FLOAT                    Absence-evidence weight (default: 0.3)\n"
                "  --em-min-prob FLOAT           Minimum EM posterior to keep (default: 1e-6)\n"
                "  --min-positional-score FLOAT  Minimum start diversity score (default: 0)\n"
                "  --min-terminal-ratio FLOAT    Minimum terminal/middle coverage ratio (default: 0)\n"
                "  --min-distinct-starts INT     Minimum distinct starts for 3+ reads (default: 2)\n"
                "  --min-informative-reads INT   Minimum informative reads (default: 3)\n"
                "  --max-single-molecule-frac FLOAT  Maximum single-molecule fraction (default: 0.5)\n"
                "EM mode:\n"
                "  --no-em                       Use best hits without EM\n"
                "  --em-use-damage-mapping       Let damage evidence affect mapping (default: off)\n";
            } else {
                std::cout << "\nRun 'dart damage-annotate --help-all' for every option.\n";
            }
            return 0;
        } else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            std::cerr << "Run 'dart damage-annotate --help' for usage.\n";
            return 1;
        }
    }

    if (emi_file.empty()) {
        std::cerr << "Error: No EMI index specified (--emi).\n";
        std::cerr << "Run 'dart damage-annotate --help' for usage.\n";
        return 1;
    }
    if (em_min_prob < 0.0 || em_min_prob > 1.0) {
        std::cerr << "Error: --em-min-prob must be in [0,1]\n";
        return 1;
    }
    if (em_null_tau < 0.0 || em_null_tau > 1.0) {
        std::cerr << "Error: --em-null-tau must be in [0,1]\n";
        return 1;
    }
    if (!(null_p_max >= 0.0 && null_p_max <= 1.0)) {
        std::cerr << "Error: --null-p-max must be in [0,1]\n";
        return 1;
    }

    // Load the empirical null model (Gumbel per length bin) if requested.
    // use_em_null gates the whole null path; disabled → zero behavior change.
    dart::NullModel em_null_model;
    bool use_em_null = false;
    std::string em_null_source = "--em-null-model";
    if (em_null_model_file.empty() && !no_null && em_calib_model_file.empty()) {
        const std::string sidecar = dart::null_sidecar_path(emi_file);
        if (std::ifstream(sidecar).good()) { em_null_model_file = sidecar; em_null_source = "sidecar next to the EMI"; }
    }
    // null_p (the per-read chance rule) reads the model on the scale it was loaded at:
    // em_null_model.mu = raw bits * null_p_scale. N is its decoy denominator.
    dart::NullProvenance em_null_prov;
    double null_p_scale = 1.0;
    uint64_t null_p_n = 0;
    bool null_p_per_submitted = false;
    if (!em_null_model_file.empty()) {
        std::string null_err;
        if (!dart::read_null_model(em_null_model_file, em_null_model, &em_null_prov, &null_err)) {
            std::cerr << "Error: cannot read null model " << em_null_model_file << ": " << null_err << "\n";
            return 1;
        }
        if (em_null_prov.legacy_tsv) {
            std::cerr << "Warning: " << em_null_model_file << " is the deprecated TSV null sidecar; "
                         "refit with dart null or hits2emi --decoy-hits to get the binary .null format.\n";
        } else {
            em_null_model.scale(em_lambda_b);   // sidecar is on the raw bit scale
        }
        null_p_scale = em_lambda_b;   // the legacy TSV is already on the lambda scale
        null_p_n = dart::null_decoy_denominator(em_null_model, em_null_prov, &null_p_per_submitted);
        use_em_null = em_null_model.valid();
        const int fitted = dart::fitted_bins(em_null_model);
        std::cerr << "Null model: " << em_null_model_file << " (" << em_null_source << ", "
                  << dart::null_fit_summary(em_null_model) << "; reject null_p > " << null_p_max
                  << (null_p_per_submitted ? " per submitted decoy" : " per decoy with a hit") << ")\n";
        if (fitted == 0)
            std::cerr << "Warning: no null bin is fitted (the decoy hits gave no finite Gumbel: fewer than two\n"
                         "  queries or zero variance): chance hits are too rare at this database size and search\n"
                         "  stringency to model. The null coordinate stays off and the output equals --no-null.\n";
    }

    // Load the coherent-LR calibration model (f_true + f_0) if requested. This
    // supersedes --em-null-model: it carries its own f_0 and switches the real
    // component to the normalized true-hit density.
    dart::CalibModel em_calib_model;
    bool use_em_calib = false;
    if (!em_calib_model_file.empty()) {
        if (!dart::read_calib_model(em_calib_model_file, em_calib_model) ||
            !em_calib_model.valid()) {
            std::cerr << "Error: cannot read --em-calib-model " << em_calib_model_file << "\n";
            return 1;
        }
        use_em_calib = true;
        // The calibration file supplies f_0 too, so the null path is active.
        em_null_model.mu = em_calib_model.null.mu;
        em_null_model.beta = em_calib_model.null.beta;
        em_null_model.n = em_calib_model.null.n;
        use_em_null = true;
        null_p_scale = 1.0;   // the calibration f_0 is on the raw bit scale
        null_p_n = dart::null_decoy_denominator(em_null_model, dart::NullProvenance{}, &null_p_per_submitted);
        if (verbose) {
            std::cerr << "Loaded calibration model (f_true " << em_calib_model.true_bins.size()
                      << " bins, f_0 " << em_calib_model.null.mu.size() << " bins, censor="
                      << em_calib_model.censor << ") from " << em_calib_model_file << "\n";
        }
    }
    if (null_p_max_set && !use_em_null)
        std::cerr << "Warning: --null-p-max has no effect without a null model (--em-null-model or <emi>.null)\n";
    if (w0 < 0.0f || w0 > 1.0f) {
        std::cerr << "Error: --w0 must be in [0, 1]\n";
        return 1;
    }
    if (site_cap <= 0.0f) {
        std::cerr << "Error: --site-cap must be > 0\n";
        return 1;
    }
    if (terminal_threshold < 0.0f || terminal_threshold > 1.0f) {
        std::cerr << "Error: --terminal-threshold must be in [0, 1]\n";
        return 1;
    }
    if (modern_threshold < 0.0f || modern_threshold >= 1.0f) {
        std::cerr << "Error: --modern-threshold must be in [0, 1)\n";
        return 1;
    }
    if (ancient_threshold <= 0.0f || ancient_threshold > 1.0f) {
        std::cerr << "Error: --ancient-threshold must be in (0, 1]\n";
        return 1;
    }
    if (ancient_threshold <= modern_threshold) {
        std::cerr << "Error: --ancient-threshold must be > --modern-threshold\n";
        return 1;
    }
    if (modern_damage_rate <= 0.0f || modern_damage_rate >= 0.5f) {
        std::cerr << "Error: --modern-damage-rate/--eps0 must be in (0, 0.5)\n";
        return 1;
    }
    if (auth_alpha <= 0.0f || auth_alpha >= 0.5f) {
        std::cerr << "Error: --auth-alpha must be in (0, 0.5)\n";
        return 1;
    }
    if (per_read_damage_threshold <= 0.0f || per_read_damage_threshold >= 1.0f) {
        std::cerr << "Error: --per-read-damage-threshold must be in (0, 1)\n";
        return 1;
    }


    int runtime_threads = 1;
#ifdef _OPENMP
    runtime_threads = (threads > 0) ? threads : std::max(1, omp_get_num_procs());
    omp_set_dynamic(0);
    omp_set_num_threads(runtime_threads);
#else
    if (threads > 0) {
        std::cerr << "Warning: --threads ignored (OpenMP not enabled in this build)\n";
    }
#endif

    if (verbose) {
        std::cerr << "Damage annotation v" << DART_VERSION << "\n";
        std::cerr << "EMI: " << emi_file << "\n";
        if (d_max >= 0.0f) {
            std::cerr << "d_max: " << d_max << "\n";
        } else {
            std::cerr << "d_max: auto\n";
        }
        std::cerr << "lambda: " << lambda << "\n";
        std::cerr << "threads: " << runtime_threads
#ifdef _OPENMP
                  << "\n";
#else
                  << " (OpenMP disabled)\n";
#endif
        std::cerr << "\n";
    }

    // Load damage index if provided (for synonymous damage detection)
    std::unique_ptr<DamageIndexReader> damage_index;
    if (!damage_index_file.empty()) {
        try {
            damage_index = std::make_unique<DamageIndexReader>(damage_index_file);
            if (verbose) {
                std::cerr << "Loaded damage index: " << damage_index->record_count() << " records\n";
                std::cerr << "Index d_max: " << (damage_index->d_max() * 100.0f) << "%\n";
                std::cerr << "Index lambda: " << damage_index->lambda() << "\n";
                std::cerr << "Damage informative: "
                          << (damage_index->damage_informative() ? "yes" : "no") << "\n";
                if (!damage_index->damage_informative()) {
                    std::cerr << "  (damage_validated=" << damage_index->damage_validated()
                              << " damage_artifact=" << damage_index->damage_artifact()
                              << " terminal_shift=" << std::fixed << std::setprecision(4)
                              << damage_index->terminal_shift() << ")\n";
                }
                std::cerr << "\n";
            }

            // Warm up page cache - critical for NFS performance
            // Forces entire .agd file into RAM before random lookups begin
            if (verbose) {
                std::cerr << "Warming up damage index cache ("
                          << (damage_index->record_count() * sizeof(AgdRecord) / 1024 / 1024)
                          << " MB)...\n";
            }
            auto t_warmup = std::chrono::steady_clock::now();
            damage_index->warmup_cache();
            if (verbose) {
                auto warmup_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - t_warmup).count();
                std::cerr << "Cache warmup complete (" << warmup_ms << " ms)\n\n";
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Cannot load damage index: " << e.what() << "\n";
            std::cerr << "Synonymous damage detection disabled.\n\n";
        }
    }

    // Open EMI index (v7 native reader — thread-safe per-block decode + get_qaln)
    dart::emiv7::Reader reader(emi_file);
    if (!reader.is_valid()) {
        std::cerr << "Error: Cannot open EMI index: " << emi_file << "\n";
        return 1;
    }
    if (!reader.has_alignment_strings()) {
        std::cerr << "Error: EMI index is missing alignment strings (qaln/taln).\n";
        std::cerr << "Rebuild with hits2emi (alignment strings are required).\n";
        return 1;
    }

    const uint32_t n_reads = reader.num_reads();
    const uint32_t n_refs = reader.num_refs();

    // Per-read ancient prior for post-hoc gamma_ancient: the AGD p_damaged, used only when the
    // library is pi DETECTED, where libtaph wrote its prior-aware read posterior. Otherwise the
    // AGD byte is sigmoid(prior-free LLR) under an implicit 0.5 prior, so gamma_ancient is left
    // empty instead of being reported from it.
    const bool agd_read_prior = damage_index && damage_index->pi_detected();
    std::vector<uint8_t> agd_dmg_q;
    if (agd_read_prior) {
        agd_dmg_q.assign(n_reads, 0);
        uint64_t agd_missing = 0;
        #pragma omp parallel for schedule(static) reduction(+:agd_missing)
        for (uint32_t r = 0; r < n_reads; ++r) {
            if (const AgdRecord* rec = damage_index->find(reader.read_name(r))) agd_dmg_q[r] = rec->p_damaged_q;
            else ++agd_missing;
        }
        if (agd_missing > 0)
            std::cerr << "Warning: " << agd_missing << " of " << n_reads
                      << " EMI reads are missing from the damage index; their p_damaged is 0\n";
    }
    if (verbose) {
        std::cerr << "Per-read ancient prior: "
                  << (agd_read_prior ? "AGD p_damaged (library pi DETECTED)"
                                     : "none; gamma_ancient left empty (library not pi DETECTED)")
                  << "\n";
    }
    auto read_dmg = [&](const float* col, size_t i, uint32_t ridx) {
        return agd_read_prior ? agd_dmg_q[ridx] * (1.0f / 255.0f) : std::clamp(col[i], 0.0f, 1.0f);
    };

    if (verbose) {
        std::cerr << "EMI metadata:\n";
        std::cerr << "  Version: " << reader.version() << "\n";
        std::cerr << "  Alignments: " << reader.num_alignments() << "\n";
        std::cerr << "  Reads: " << n_reads << "\n";
        std::cerr << "  References: " << n_refs << "\n";
        std::cerr << "  Sorted (read,score): "
                  << (reader.is_sorted_by_read_then_score() ? "yes" : "no") << "\n";
        std::cerr << "\n";
    }

    std::vector<ProteinDamageSummary> summaries;
    std::vector<int32_t> summary_index_by_read(n_reads, -1);

    struct BestHitData {
        bool has = false;
        float bits = 0.0f;
        float bits_second = 0.0f;  // second-best bit score for this query
        float damage_score = 0.0f;
        uint32_t ref_idx = 0;  // For EM: ref_idx corresponding to target_id
        float fident = 0.0f;
        float evalue = 0.0f;
        uint32_t qstart_0 = 0;
        uint32_t tstart_0 = 0;
        uint32_t qlen = 0;
        uint32_t tlen = 0;
        uint32_t aln_len_on_target = 0;
        uint16_t aln_len = 0;  // raw EMI ALN_LEN of the best hit: the null bin (the fit binned on it)
        uint32_t best_rg = UINT32_MAX;
        uint32_t best_row = UINT32_MAX;
        // Soft multi-mapping (--emit-multimap): the EM mapping responsibility for THIS
        // (read,ref), joined from the finalize; candidate_rank assigned 0..k-1 by gamma.
        float    gamma = 0.0f, gamma_ancient = 0.0f;
        uint8_t  candidate_rank = 0;
    };
    std::vector<BestHitData> best_hits(n_reads);
    // Per-read top-K-by-bits candidate POOL for --emit-multimap. Slim element: only the
    // trim key (bits), the row locators (ref_idx/best_rg/best_row) and the post-EM gamma
    // join are retained per candidate; every column-derived field (fident/evalue/qstart/
    // tstart/qlen/tlen/damage_score) is re-fetched byte-identically in the annotation scan
    // (which already re-visits these row groups). This drops the pooled element ~68->~28 B,
    // turning the O(n_reads * soft_cap) resident pool from the peak-RSS driver (measured
    // 150-260 GB at 98M reads) into a fraction of that. gamma is joined post-EM, then the
    // pool is filtered rank<K & gamma>floor.
    struct MmCand {
        float    bits = 0.0f;                   // scan-time top-K-by-bits trim key
        uint32_t ref_idx = 0;                   // finalize one-HSP-per-ref dedup + target_id
        uint32_t best_rg = UINT32_MAX;
        uint32_t best_row = UINT32_MAX;
        float    gamma = 0.0f, gamma_ancient = 0.0f;  // EM-derived, joined post-EM
        uint8_t  candidate_rank = 0;            // assigned 0..k-1 by gamma at finalize
    };
    std::vector<std::vector<MmCand>> mm_pool;
    if (emit_multimap) mm_pool.assign(n_reads, {});
    // Production streaming EM tape path (PI mandate 2026-07-23): when the config needs
    // only ref/gamma/gamma_ancient per read (no calib, damage-free mapping) the EM runs
    // off an O_DIRECT tape on node-local scratch, feeding tape_em_core directly from
    // streaming ingestion. Nothing O(alignments) is materialized in RAM (kills the ~288 GB
    // in-RAM transient measured on LV626), giving aln-independent peak RSS for the 6-10B-aln
    // production target. The decoy null rides along as one double per read (the best null
    // log-density over the read's alignments, computed while the tape is written), so a
    // null model keeps the tape. All predicates here are CLI-known, so the O(A)
    // em_alignments buffers below are skipped entirely on this path.
    const bool tape_path = use_em && !use_em_calib && !em_use_damage_mapping;
    // read-major per-read degree (uint64 prefix: aln offsets exceed 2^32 at 6B+ alns) and
    // per-ref length (== build_alignment_data's ref_lengths). Tape path only.
    std::vector<uint64_t> em_read_deg;    // [n_reads+1], prefix-summed into read_offsets
    std::vector<uint32_t> em_ref_lengths; // max(tlen | aln_end) per ref
    if (tape_path) {
        em_read_deg.assign(static_cast<size_t>(n_reads) + 1, 0);
        em_ref_lengths.assign(n_refs, 0);
    }
    std::vector<dart::CompactAlignment> em_alignments;
    if (use_em && !tape_path) {
        em_alignments.reserve(static_cast<size_t>(reader.num_alignments()));
    }
    std::mutex em_damage_stats_mutex;
    uint64_t em_damage_n = 0;
    uint64_t em_damage_nonzero = 0;
    double em_damage_sum = 0.0;
    double em_damage_sumsq = 0.0;
    uint64_t em_damage_ll_n = 0;
    uint64_t em_damage_ll_nonzero = 0;
    double em_damage_llbf_sum = 0.0;
    double em_damage_llbf_sumsq = 0.0;
    float em_damage_min = 1.0f;
    float em_damage_max = 0.0f;
    struct RgDamageBoundary {
        bool has_rows = false;
        bool intra_read_variation = false;
        uint32_t first_read = 0;
        uint32_t last_read = 0;
        float first_ds = 0.0f;
        float last_ds = 0.0f;
    };
    std::vector<RgDamageBoundary> em_damage_boundaries;
    if (use_em) {
        em_damage_boundaries.resize(reader.num_row_groups());
    }

    std::atomic<size_t> total_alignments{0};
    std::atomic<size_t> filtered_identity{0};
    std::atomic<size_t> filtered_bits{0};
    std::atomic<size_t> filtered_evalue{0};
    constexpr size_t kBestHitShards = 4096;
    std::array<std::mutex, kBestHitShards> best_hit_mutexes;
    const float aln_max_evalue_log10 = (aln_max_evalue > 0.0f)
        ? std::log10(aln_max_evalue)
        : -std::numeric_limits<float>::infinity();

    size_t em_shard_count = 1;
#ifdef _OPENMP
    if (use_em && !tape_path) {
        em_shard_count = static_cast<size_t>(std::max(1, omp_get_max_threads()));
    }
#endif
    std::vector<std::vector<dart::CompactAlignment>> em_alignments_by_thread;
    if (use_em && !tape_path) {
        em_alignments_by_thread.resize(em_shard_count);
        const size_t reserve_per_shard =
            static_cast<size_t>(reader.num_alignments() / std::max<size_t>(1, em_shard_count) + 1);
        for (auto& shard : em_alignments_by_thread) {
            shard.reserve(reserve_per_shard);
        }
    }

    // Pass 1: numeric-only scan for filtering, best-hit selection, and EM payload.
    const auto pass1_start = std::chrono::steady_clock::now();
    reader.set_columns({
        dart::ColumnID::READ_IDX,
        dart::ColumnID::REF_IDX,
        dart::ColumnID::BIT_SCORE,
        dart::ColumnID::DAMAGE_SCORE,
        dart::ColumnID::EVALUE_LOG10,
        dart::ColumnID::DMG_K,
        dart::ColumnID::DMG_M,
        dart::ColumnID::DMG_LL_A,
        dart::ColumnID::DMG_LL_M,
        dart::ColumnID::IDENTITY_Q,
        dart::ColumnID::ALN_LEN,
        dart::ColumnID::QSTART,
        dart::ColumnID::TSTART,
        dart::ColumnID::TLEN,
        dart::ColumnID::QLEN
    });
    reader.parallel_scan([&](
        uint32_t rg_idx, uint32_t num_rows,
        const uint32_t* read_idx, const uint32_t* ref_idx,
        const float* bit_score, const float* damage_score,
        const float* evalue_log10,
        const uint16_t* /*dmg_k*/, const uint16_t* /*dmg_m*/,
        const float* dmg_ll_a, const float* dmg_ll_m,
        const uint16_t* identity_q,
        const uint16_t* aln_len, const uint16_t* qstart, const uint16_t* /*qend*/,
        const uint16_t* tstart, const uint16_t* /*tend*/,
        const uint16_t* tlen, const uint16_t* qlen,
        const uint16_t* /*mismatch*/, const uint16_t* /*gapopen*/)
    {
        if (!read_idx || !ref_idx || !bit_score || !evalue_log10 ||
            !identity_q || !aln_len || !qstart || !tstart || !tlen || !qlen || !damage_score) {
            return;
        }

        struct BestHitCandidate {
            float bits = 0.0f;
            float bits_second = 0.0f;
            float damage_score = 0.0f;
            uint32_t ref_idx = 0;
            float fident = 0.0f;
            float evalue = 0.0f;
            uint32_t qstart_0 = 0;
            uint32_t tstart_0 = 0;
            uint32_t qlen = 0;
            uint32_t tlen = 0;
            uint32_t aln_len_on_target = 0;
            uint32_t best_rg = UINT32_MAX;
            uint32_t best_row = UINT32_MAX;
        };

        size_t local_total = 0;
        size_t local_filtered_identity = 0;
        size_t local_filtered_bits = 0;
        size_t local_filtered_evalue = 0;

        std::vector<dart::CompactAlignment> local_em;
        if (use_em && !tape_path) {
            local_em.reserve(num_rows);
        }
        uint64_t local_em_n = 0;
        uint64_t local_em_nonzero = 0;
        double local_em_sum = 0.0;
        double local_em_sumsq = 0.0;
        uint64_t local_em_ll_n = 0;
        uint64_t local_em_ll_nonzero = 0;
        double local_em_llbf_sum = 0.0;
        double local_em_llbf_sumsq = 0.0;
        float local_em_min = 1.0f;
        float local_em_max = 0.0f;
        bool local_em_has_rows = false;
        bool local_em_intra_read_variation = false;
        uint32_t local_em_first_read = 0;
        uint32_t local_em_last_read = 0;
        float local_em_first_ds = 0.0f;
        float local_em_last_ds = 0.0f;
        std::unordered_map<uint32_t, float> local_em_read_damage;
        if (use_em) {
            local_em_read_damage.reserve(std::max<uint32_t>(1, num_rows / 8));
        }

        bool local_sorted_hits = true;
        for (uint32_t i = 1; i < num_rows; ++i) {
            if (read_idx[i] < read_idx[i - 1] ||
                (read_idx[i] == read_idx[i - 1] && bit_score[i] > bit_score[i - 1])) {
                local_sorted_hits = false;
                break;
            }
        }

        std::vector<std::pair<uint32_t, BestHitCandidate>> local_best_sorted;
        std::unordered_map<uint32_t, BestHitCandidate> local_best_unsorted;
        // --emit-multimap: collect EVERY passing (read,ref) candidate (not collapsed to
        // argmax) so the pool can keep the true top-K; merged into mm_pool after the loop.
        std::vector<std::pair<uint32_t, MmCand>> local_mm;
        if (local_sorted_hits) {
            local_best_sorted.reserve(std::max<uint32_t>(1, num_rows / 8));
        } else {
            local_best_unsorted.reserve(std::max<uint32_t>(1, num_rows / 4));
        }


        auto update_candidate = [&](BestHitCandidate& cand,
                                    float bits,
                                    float damage_score_v,
                                    uint32_t tidx,
                                    float fident,
                                    float evalue,
                                    uint32_t qstart_0,
                                    uint32_t tstart_0,
                                    uint32_t qlen_v,
                                    uint32_t tlen_v,
                                    uint32_t aln_len_on_target_v,
                                    uint32_t row_in_group) {
            if (bits > cand.bits) {
                cand.bits_second = cand.bits;
                cand.bits = bits;
                cand.damage_score = damage_score_v;
                cand.ref_idx = tidx;
                cand.fident = fident;
                cand.evalue = evalue;
                cand.qstart_0 = qstart_0;
                cand.tstart_0 = tstart_0;
                cand.qlen = qlen_v;
                cand.tlen = tlen_v;
                cand.aln_len_on_target = aln_len_on_target_v;
                cand.best_rg = rg_idx;
                cand.best_row = row_in_group;
            } else if (bits > cand.bits_second) {
                cand.bits_second = bits;
            }
        };

        for (uint32_t i = 0; i < num_rows; ++i) {
            local_total++;

            const float fident = static_cast<float>(identity_q[i]) / 65535.0f;
            const float bits = bit_score[i];
            const float evalue_log10_row = evalue_log10[i];
            if (fident < aln_min_identity) {
                local_filtered_identity++;
                continue;
            }
            if (bits < aln_min_bits) {
                local_filtered_bits++;
                continue;
            }
            if (evalue_log10_row > aln_max_evalue_log10) {
                local_filtered_evalue++;
                continue;
            }

            const uint32_t ridx = read_idx[i];
            const uint32_t tidx = ref_idx[i];
            const uint32_t qstart_0 = (qstart[i] > 0) ? static_cast<uint32_t>(qstart[i] - 1) : 0u;
            const uint32_t tstart_0 = (tstart[i] > 0) ? static_cast<uint32_t>(tstart[i] - 1) : 0u;
            const uint32_t tlen_v = static_cast<uint32_t>(tlen[i]);
            const uint32_t qlen_v = static_cast<uint32_t>(qlen[i]);
            const uint16_t aln_len_v = aln_len[i];

            if (emit_multimap) {
                // Slim pooled candidate: only trim key + locators. The column-derived
                // metadata (fident/evalue/qstart/tstart/qlen/tlen/damage_score) is re-read
                // byte-identically from (best_rg,best_row) in the annotation scan below.
                MmCand d;
                d.bits = bits; d.ref_idx = tidx;
                d.best_rg = rg_idx; d.best_row = i;
                local_mm.emplace_back(ridx, d);
            }


            if (use_em) {
                dart::CompactAlignment ca{};
                ca.read_idx = ridx;
                ca.ref_idx = tidx;
                ca.bit_score = bits;
                const float ds = read_dmg(damage_score, i, ridx);
                ca.damage_score = ds;
                const float ll_a = dmg_ll_a ? dmg_ll_a[i] : 0.0f;
                const float ll_m = dmg_ll_m ? dmg_ll_m[i] : 0.0f;
                ca.damage_ll_a = ll_a;
                ca.damage_ll_m = ll_m;
                if (!local_em_has_rows) {
                    local_em_has_rows = true;
                    local_em_first_read = ridx;
                    local_em_first_ds = ds;
                }
                local_em_last_read = ridx;
                local_em_last_ds = ds;
                auto [it_ds, inserted_ds] = local_em_read_damage.try_emplace(ridx, ds);
                if (!inserted_ds && std::fabs(it_ds->second - ds) > 1e-6f) {
                    local_em_intra_read_variation = true;
                }
                ca.aln_start = static_cast<uint16_t>(std::min(tstart_0, 65535u));
                ca.aln_end = static_cast<uint16_t>(
                    std::min(tstart_0 + static_cast<uint32_t>(aln_len_v), 65535u));
                ca.identity_q = identity_q[i];
                ca.flags = tlen[i];  // pass true reference length hint to EM builder
                if (tape_path) {
                    // read-major: each read lives in one row group (one thread), so these
                    // shared vectors are written at disjoint read indices — no atomics needed
                    // (guarded by is_sorted_by_read_then_score() before the solve). The
                    // em_ref_lengths max is racy but benign: tlen is constant per ref, so all
                    // writers store the same value. ceiling: a tlen==0 ref with varying aln_end
                    // could under-max under the race; upgrade: serial reduce (real refs tlen>0).
                    em_read_deg[ridx + 1]++;
                    const uint32_t ref_len =
                        (tlen[i] > 0) ? static_cast<uint32_t>(tlen[i]) : static_cast<uint32_t>(ca.aln_end);
                    if (ref_len > em_ref_lengths[tidx]) em_ref_lengths[tidx] = ref_len;
                } else {
                    local_em.push_back(ca);
                }
                local_em_n++;
                local_em_sum += ds;
                local_em_sumsq += static_cast<double>(ds) * static_cast<double>(ds);
                const double llbf = static_cast<double>(ll_a) - static_cast<double>(ll_m);
                local_em_ll_n++;
                if (std::fabs(llbf) > 1e-9) local_em_ll_nonzero++;
                local_em_llbf_sum += llbf;
                local_em_llbf_sumsq += llbf * llbf;
                if (ds > 0.0f) local_em_nonzero++;
                local_em_min = std::min(local_em_min, ds);
                local_em_max = std::max(local_em_max, ds);
            }

            const float evalue = evalue_from_log10(evalue_log10_row);
            if (local_sorted_hits) {
                if (local_best_sorted.empty() || local_best_sorted.back().first != ridx) {
                    BestHitCandidate cand;
                    cand.bits = bits;
                    cand.bits_second = 0.0f;
                    cand.damage_score = read_dmg(damage_score, i, ridx);
                    cand.ref_idx = tidx;
                    cand.fident = fident;
                    cand.evalue = evalue;
                    cand.qstart_0 = qstart_0;
                    cand.tstart_0 = tstart_0;
                    cand.qlen = qlen_v;
                    cand.tlen = tlen_v;
                    cand.aln_len_on_target = static_cast<uint32_t>(aln_len_v);
                    cand.best_rg = rg_idx;
                    cand.best_row = i;
                    local_best_sorted.emplace_back(ridx, cand);
                } else {
                    update_candidate(local_best_sorted.back().second,
                                     bits, read_dmg(damage_score, i, ridx),
                                     tidx, fident, evalue,
                                     qstart_0, tstart_0, qlen_v, tlen_v,
                                     static_cast<uint32_t>(aln_len_v), i);
                }
            } else {
                auto [it, inserted] = local_best_unsorted.try_emplace(ridx);
                if (inserted) {
                    BestHitCandidate cand;
                    cand.bits = bits;
                    cand.bits_second = 0.0f;
                    cand.damage_score = read_dmg(damage_score, i, ridx);
                    cand.ref_idx = tidx;
                    cand.fident = fident;
                    cand.evalue = evalue;
                    cand.qstart_0 = qstart_0;
                    cand.tstart_0 = tstart_0;
                    cand.qlen = qlen_v;
                    cand.tlen = tlen_v;
                    cand.aln_len_on_target = static_cast<uint32_t>(aln_len_v);
                    cand.best_rg = rg_idx;
                    cand.best_row = i;
                    it->second = cand;
                } else {
                    update_candidate(it->second,
                                     bits, read_dmg(damage_score, i, ridx),
                                     tidx, fident, evalue,
                                     qstart_0, tstart_0, qlen_v, tlen_v,
                                     static_cast<uint32_t>(aln_len_v), i);
                }
            }
        }

        total_alignments.fetch_add(local_total, std::memory_order_relaxed);
        filtered_identity.fetch_add(local_filtered_identity, std::memory_order_relaxed);
        filtered_bits.fetch_add(local_filtered_bits, std::memory_order_relaxed);
        filtered_evalue.fetch_add(local_filtered_evalue, std::memory_order_relaxed);

        auto merge_candidate = [&](uint32_t ridx, const BestHitCandidate& cand) {
            const size_t shard = static_cast<size_t>(ridx) & (kBestHitShards - 1);
            std::lock_guard<std::mutex> lock(best_hit_mutexes[shard]);
            BestHitData& hit = best_hits[ridx];

            if (!hit.has || cand.bits > hit.bits) {
                const float prev_best = hit.has ? hit.bits : 0.0f;
                hit.has = true;
                hit.bits = cand.bits;
                hit.bits_second = std::max(cand.bits_second, std::max(hit.bits_second, prev_best));
                hit.damage_score = cand.damage_score;
                hit.ref_idx = cand.ref_idx;
                hit.fident = cand.fident;
                hit.evalue = cand.evalue;
                hit.qstart_0 = cand.qstart_0;
                hit.tstart_0 = cand.tstart_0;
                hit.qlen = cand.qlen;
                hit.tlen = cand.tlen;
                hit.aln_len_on_target = cand.aln_len_on_target;
                // aln_len_on_target is still the raw aln_len_v here (the annotation scan clips it later).
                hit.aln_len = static_cast<uint16_t>(std::min(cand.aln_len_on_target, 65535u));
                hit.best_rg = cand.best_rg;
                hit.best_row = cand.best_row;
            } else {
                hit.bits_second = std::max(hit.bits_second, std::max(cand.bits, cand.bits_second));
            }
        };

        if (local_sorted_hits) {
            for (const auto& kv : local_best_sorted) {
                merge_candidate(kv.first, kv.second);
            }
        } else {
            for (const auto& kv : local_best_unsorted) {
                merge_candidate(kv.first, kv.second);
            }
        }

        if (emit_multimap && !local_mm.empty()) {
            const size_t soft_cap = static_cast<size_t>(8 * std::max(1, multimap_k));
            for (auto& kv : local_mm) {
                const uint32_t ridx = kv.first;
                const size_t shard = static_cast<size_t>(ridx) & (kBestHitShards - 1);
                std::lock_guard<std::mutex> lock(best_hit_mutexes[shard]);
                auto& pool = mm_pool[ridx];
                pool.push_back(kv.second);
                if (pool.size() > soft_cap) {  // keep memory bounded; final top-K taken post-EM
                    std::sort(pool.begin(), pool.end(),
                              [](const MmCand& a, const MmCand& b) { return a.bits > b.bits; });
                    // Preserve the intended 8*K pre-EM candidate reservoir. Shrinking to K
                    // here discards candidates before gamma exists and makes the later
                    // gamma-ranked top-K depend on row-group ingestion order.
                    pool.resize(soft_cap);
                }
            }
        }


        if (use_em && local_em_has_rows) {
            if (!local_em.empty()) {
                size_t em_shard = 0;
#ifdef _OPENMP
                em_shard = static_cast<size_t>(omp_get_thread_num());
#endif
                if (em_shard >= em_alignments_by_thread.size()) em_shard = 0;
                auto& dst = em_alignments_by_thread[em_shard];
                dst.insert(dst.end(), local_em.begin(), local_em.end());
            }

            std::lock_guard<std::mutex> lock(em_damage_stats_mutex);
            em_damage_n += local_em_n;
            em_damage_nonzero += local_em_nonzero;
            em_damage_sum += local_em_sum;
            em_damage_sumsq += local_em_sumsq;
            em_damage_ll_n += local_em_ll_n;
            em_damage_ll_nonzero += local_em_ll_nonzero;
            em_damage_llbf_sum += local_em_llbf_sum;
            em_damage_llbf_sumsq += local_em_llbf_sumsq;
            em_damage_min = std::min(em_damage_min, local_em_min);
            em_damage_max = std::max(em_damage_max, local_em_max);
            if (rg_idx < em_damage_boundaries.size()) {
                em_damage_boundaries[rg_idx] = {
                    local_em_has_rows,
                    local_em_intra_read_variation,
                    local_em_first_read,
                    local_em_last_read,
                    local_em_first_ds,
                    local_em_last_ds
                };
            }
        }
    });

    bool em_damage_varies_within_read = false;
    if (use_em) {
        RgDamageBoundary prev{};
        bool have_prev = false;
        for (const auto& b : em_damage_boundaries) {
            if (!b.has_rows) continue;
            if (b.intra_read_variation) {
                em_damage_varies_within_read = true;
                break;
            }
            if (have_prev && prev.last_read == b.first_read &&
                std::fabs(prev.last_ds - b.first_ds) > 1e-6f) {
                em_damage_varies_within_read = true;
                break;
            }
            prev = b;
            have_prev = true;
        }
    }

    if (use_em) {
        size_t em_total = 0;
        for (const auto& shard : em_alignments_by_thread) {
            em_total += shard.size();
        }
        em_alignments.clear();
        em_alignments.reserve(em_total);
        for (auto& shard : em_alignments_by_thread) {
            em_alignments.insert(em_alignments.end(), shard.begin(), shard.end());
        }
        // Fix 2: free shard memory immediately after merge
        em_alignments_by_thread.clear();
        std::vector<std::vector<dart::CompactAlignment>>{}.swap(em_alignments_by_thread);
    }
    const auto pass1_end = std::chrono::steady_clock::now();
    if (verbose) {
        std::cerr << "Pass 1 runtime: "
                  << dart::log_utils::format_elapsed(pass1_start, pass1_end) << "\n";
    }

    // Release NFS pages accumulated during Pass 1 before the annotation pass.
    // MADV_SEQUENTIAL readahead holds the entire 114 GB EMI resident even with
    // per-row-group DONTNEED; flush_pages() calls posix_fadvise+madvise DONTNEED
    // on the full range and switches to MADV_RANDOM so the annotation
    // parallel_scan_selected (which uses explicit WILLNEED per selected row group)
    // won't trigger speculative readahead on the full file.
    reader.flush_pages();

    // Combined Pass 2+3: load strings per row group, annotate inline, release pages.
    // Strings are never copied into best_hits — they are string_views from the mmap,
    // valid only during their row group's callback. Peak string memory ≈ 8 threads ×
    // one row group's qaln/taln data (~8 MB each) instead of ~20 GB across all reads.
    const auto pass2_start = std::chrono::steady_clock::now();

    // Build the row-group index: maps rg_idx → sorted list of (row_in_group, ridx).
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> best_rows_by_rg(reader.num_row_groups());
    std::vector<uint32_t> selected_row_groups;
    selected_row_groups.reserve(reader.num_row_groups());
    size_t selected_best_hits = 0;
    for (uint32_t ridx = 0; ridx < n_reads; ++ridx) {
        const auto& hit = best_hits[ridx];
        if (!hit.has || hit.best_rg == UINT32_MAX || hit.best_rg >= best_rows_by_rg.size()) continue;
        if (best_rows_by_rg[hit.best_rg].empty()) {
            selected_row_groups.push_back(hit.best_rg);
        }
        best_rows_by_rg[hit.best_rg].emplace_back(hit.best_row, ridx);
        selected_best_hits++;
    }
    for (auto& rows : best_rows_by_rg) {
        std::sort(rows.begin(), rows.end());
    }
    std::sort(selected_row_groups.begin(), selected_row_groups.end());

    // Determine d_max before annotation (annotate_alignment needs it as a parameter).
    if (d_max < 0.0f) {
        if (damage_index && damage_index->d_max() > 0.0f) {
            d_max = damage_index->d_max();
            if (!lambda_set_by_user && damage_index->lambda() > 0.0f)
                lambda = damage_index->lambda();
            if (verbose) {
                std::cerr << "Using d_max from damage index header: " << std::fixed
                          << std::setprecision(3) << d_max << "\n";
                if (!lambda_set_by_user)
                    std::cerr << "Using lambda from damage index header: " << std::fixed
                              << std::setprecision(3) << lambda << "\n";
            }
        } else if (reader.d_max() > 0.0f) {
            d_max = reader.d_max();
            if (!lambda_set_by_user && reader.lambda() > 0.0f)
                lambda = reader.lambda();
            if (verbose) {
                std::cerr << "Using d_max from EMI header: " << std::fixed
                          << std::setprecision(3) << d_max << "\n";
                if (!lambda_set_by_user)
                    std::cerr << "Using lambda from EMI header: " << std::fixed
                              << std::setprecision(3) << lambda << "\n";
            }
        } else {
            // Rare fallback: estimate d_max from alignment strings.
            // Requires a preliminary sequential pass through strings before the main annotation.
            if (verbose)
                std::cerr << "Estimating d_max from EMI best hits...\n";
            reader.set_columns({dart::ColumnID::QALN, dart::ColumnID::TALN});
            DamageEstimate est;
            for (const uint32_t rg_idx : selected_row_groups) {
                for (const auto& [row_in_group, ridx] : best_rows_by_rg[rg_idx]) {
                    std::string_view qaln_sv = reader.get_qaln(rg_idx, row_in_group);
                    std::string_view taln_sv = reader.get_taln(rg_idx, row_in_group);
                    std::string query_id(reader.read_name(ridx));
                    count_damage_for_estimate(query_id, qaln_sv, taln_sv,
                        best_hits[ridx].qstart_0, best_hits[ridx].qlen, est);
                }
            }
            d_max = est.estimate_d_max();
            if (verbose) {
                std::cerr << "  Terminal subs: " << est.terminal_damage << "/" << est.terminal_total
                          << " (" << (est.terminal_total > 0 ? 100.0f * est.terminal_damage / est.terminal_total : 0) << "%)\n";
                std::cerr << "  Interior subs: " << est.interior_damage << "/" << est.interior_total
                          << " (" << (est.interior_total > 0 ? 100.0f * est.interior_damage / est.interior_total : 0) << "%)\n";
                std::cerr << "  Estimated d_max: " << std::fixed << std::setprecision(3) << d_max << "\n";
            }
        }
        if (verbose) std::cerr << "\n";
    }


    if (verbose) {
        size_t total = total_alignments.load(std::memory_order_relaxed);
        size_t filt_id = filtered_identity.load(std::memory_order_relaxed);
        size_t filt_bits = filtered_bits.load(std::memory_order_relaxed);
        size_t filt_eval = filtered_evalue.load(std::memory_order_relaxed);
        size_t passed = total - filt_id - filt_bits - filt_eval;
        std::cerr << "Alignment pre-filtering:\n";
        std::cerr << "  Total alignments: " << total << "\n";
        if (filt_id > 0 || aln_min_identity > 0)
            std::cerr << "  Filtered (identity < " << aln_min_identity << "): " << filt_id << "\n";
        if (filt_bits > 0 || aln_min_bits > 0)
            std::cerr << "  Filtered (bits < " << aln_min_bits << "): " << filt_bits << "\n";
        if (filt_eval > 0 || aln_max_evalue < 1e9f)
            std::cerr << "  Filtered (evalue > " << aln_max_evalue << "): " << filt_eval << "\n";
        std::cerr << "  Passed: " << passed << " (" << (total > 0 ? 100.0 * passed / total : 0.0) << "%)\n\n";
    }

    // per_read_aa_decay declared here for --refine-damage; populated during
    // annotation pass (which runs after EM to prevent summaries coexisting with
    // EM state in memory, reducing annotation-phase peak RSS by ~25-42 GB).
    std::vector<float> per_read_aa_decay(n_reads, 0.0f);
    // Per-read C→T/G→A opportunity count for EVERY best-hit read, incl. zero-mismatch
    // (clean) reads that are dropped from `summaries`. Feeds the per-opportunity damage
    // fit and q_modern their 0/n_eligible negative evidence (BUG1 fix).
    std::vector<uint32_t> per_read_m_opp(n_reads, 0);

    // CSR read offsets (uint64: aln offsets exceed 2^32 at 6B+ alns) + total kept alns.
    // Tape path: prefix-sum the streaming per-read degrees. In-RAM path: copied from
    // aln_data.read_offsets after build_alignment_data (below). Kept in scope for both the
    // solve (tape_em_core) and the verbose per-read stats.
    std::vector<uint64_t> read_offsets;
    uint64_t em_tape_A = 0;
    if (tape_path) {
        read_offsets.assign(static_cast<size_t>(n_reads) + 1, 0);
        for (uint32_t r = 0; r < n_reads; ++r)
            read_offsets[r + 1] = read_offsets[r] + em_read_deg[r + 1];
        em_tape_A = read_offsets[n_reads];
        std::vector<uint64_t>().swap(em_read_deg);  // degrees folded into read_offsets
    }
    const size_t em_aln_count = tape_path ? static_cast<size_t>(em_tape_A)
                                          : em_alignments.size();

    // EM reassignment
    dart::EMState em_state;
    struct EmReadResult {
        float gamma = 0.0f;             // Responsibility for the best-hit target in summaries
        float gamma_ancient = 0.0f;     // Ancient responsibility for the same target
        bool best_hit = false;          // True if best-hit target == EM best target
        bool em_keep = true;            // Passes bam-filter-like posterior filter
        float em_threshold = 0.0f;      // Applied per-read keep threshold
        float em_margin = 0.0f;         // EM confidence margin: top1 - top2 responsibility
        float em_ref_neff = 1.0f;       // Effective number of references: 1 / sum_j gamma_j^2
        float gamma_null = 0.0f;        // null-component responsibility (spurious-hit mass)
        uint8_t em_reject_reason = 0;   // 0=keep, bit0=null_dominated, bit1=low_margin, bit2=high_neff
        float null_p = std::numeric_limits<float>::quiet_NaN();  // chance probability of the best hit
    };
    std::vector<EmReadResult> em_read_results(n_reads);
    // Null chance rule, per null bin: reads with a best hit and reads flagged by null_p.
    std::array<uint64_t, dart::kNumNullBins> null_bin_reads{}, null_bin_flagged{};
    bool null_p_ran = false;
    // True when damage-aware EM was active and gamma_ancient was actually populated.
    // Used by --refine-damage to avoid falling back to total gamma for modern reads.
    // Atomic because it may be set inside the parallel streaming_em_finalize callback.
    std::atomic<bool> damage_em_active{false};
    std::unordered_map<uint32_t, dart::ProteinCoverageStats> coverage_stats_map;

    // Cached ref_names (NFS-safe pread). Populated by streaming EM right after flush_pages,
    // or lazily at output phase for non-streaming paths.
    std::vector<std::string> cached_ref_names;

    // Fail-loud null-model gate — fires for BOTH the streaming and non-streaming
    // EM branches below. A flat targeted database (barley SAAV / KEGG / viral) has
    // no "other" column, so an absent-protein read is force-mapped onto its nearest
    // wrong homolog and emitted as a false annotation. Refuse EM without a null
    // model unless the caller accepts the risk with --no-null.
    if (use_em && !use_em_null && !no_null) {
        std::cerr << "Error: refusing EM annotation without a null model; "
                     "a flat targeted database can force-map absent-protein reads to "
                     "false annotations. Search the decoy peptides predict wrote next to "
                     "--fasta-aa and pass them to hits2emi --decoy-hits (annotate then finds "
                     "<emi>.null), pass --em-null-model FILE, or --no-null to accept this explicitly.\n";
        return 1;
    }
    if (use_em && em_aln_count > 0) {
        if (verbose) {
            std::cerr << "\nEM reassignment:\n";
            std::cerr << "  Alignments: " << em_aln_count << "\n";
            std::cerr << "  Reads: " << n_reads << "\n";
            std::cerr << "  References: " << n_refs << "\n";
            std::cerr << "  Lambda_b: " << em_lambda_b << "\n";
            std::cerr << "  Max iters: " << em_max_iters << "\n";
            std::cerr << "  Tolerance: " << std::scientific << em_tol << std::defaultfloat << "\n";
            std::cerr << "  Param tolerance: " << std::scientific << std::setprecision(3)
                      << std::sqrt(std::max(em_tol, 1e-15)) << std::defaultfloat << "\n";
            std::cerr << "  Min prob keep: " << em_min_prob << "\n";
            std::cerr << "  Length normalization: on\n";
        }

        // Build CSR-format alignment data (in-RAM path only). On the tape path aln_data
        // stays empty — the EM runs off the O_DIRECT tape and nothing O(A) is materialized;
        // read_offsets was already prefix-summed from the streaming degrees above.
        dart::AlignmentData aln_data;
        if (!tape_path) {
            aln_data = dart::build_alignment_data(
                em_alignments.data(), em_aln_count, n_reads, n_refs);
            // Fix 3: release em_alignments – aln_data has its own copy
            std::vector<dart::CompactAlignment>{}.swap(em_alignments);
            read_offsets.assign(aln_data.read_offsets.begin(), aln_data.read_offsets.end());
        }

        // Set up EM parameters
        dart::EMParams em_params;
        em_params.lambda_b = em_lambda_b;
        em_params.max_iters = em_max_iters;
        em_params.tol = em_tol;
        em_params.use_squarem = true;
        em_params.normalize_by_length = true;
        // Null/noise component: adds one mixture column scored by the decoy Gumbel.
        // Absolute bit-score information is required, so per-query score normalization
        // is forced off when the null is active (it would erase the null's signal).
        if (use_em_null) {
            em_params.use_null = true;
            em_params.null_mu = em_null_model.mu;
            em_params.null_beta = em_null_model.beta;
            em_params.log_w_null = std::log(0.096);   // measured hit-FDR seed
            em_params.alpha_null = 2.0;               // informative floor; can't collapse to 0
            em_params.lambda_null = 1.0;
            em_params.normalize_scores_per_query = false;
        }
        // Coherent-LR: real component becomes log f_true, f_0 on the raw-bit scale.
        if (use_em_calib) {
            em_params.use_calib = true;
            em_params.true_bins = em_calib_model.true_bins;
            em_params.true_censor = em_calib_model.censor;
            em_params.calib_temperature = em_calib_temperature;
            // v2 length-continuous section (sol short-fragment fix): supersedes the
            // discrete f_true/f_0 so short reads are scored against a length-appropriate
            // mean, not the long-read pooled fallback.
            if (em_calib_model.cont) {
                em_params.use_calib_cont = true;
                em_params.ft_mu_a = em_calib_model.ft_mu_a;
                em_params.ft_mu_b = em_calib_model.ft_mu_b;
                em_params.ft_ls_a = em_calib_model.ft_ls_a;
                em_params.ft_ls_b = em_calib_model.ft_ls_b;
                em_params.ft_nu   = em_calib_model.ft_nu;
                em_params.f0_loc_a = em_calib_model.f0_loc_a;
                em_params.f0_loc_b = em_calib_model.f0_loc_b;
                em_params.f0_beta  = em_calib_model.f0_beta;
            }
        }
        if (verbose) {
            std::cerr << "  Min iters: " << em_params.min_iters << "\n";
            std::cerr << "  SQUAREM warmup iters: " << em_params.squarem_warmup_iters << "\n";
            if (use_em_null)
                std::cerr << "  Null component: on (w_null init=0.096, tau=" << em_null_tau << ")\n";
        }
        const double em_damage_mean =
            (em_damage_n > 0) ? (em_damage_sum / static_cast<double>(em_damage_n)) : 0.0;
        const double em_damage_var = (em_damage_n > 0)
            ? std::max(0.0, (em_damage_sumsq / static_cast<double>(em_damage_n)) -
                               (em_damage_mean * em_damage_mean))
            : 0.0;
        const double em_damage_std = std::sqrt(em_damage_var);
        const double em_llbf_mean =
            (em_damage_ll_n > 0) ? (em_damage_llbf_sum / static_cast<double>(em_damage_ll_n)) : 0.0;
        const double em_llbf_var = (em_damage_ll_n > 0)
            ? std::max(0.0, (em_damage_llbf_sumsq / static_cast<double>(em_damage_ll_n)) -
                               (em_llbf_mean * em_llbf_mean))
            : 0.0;
        const double em_llbf_std = std::sqrt(em_llbf_var);
        const float em_damage_min_print = (em_damage_n > 0) ? em_damage_min : 0.0f;
        const float em_damage_max_print = (em_damage_n > 0) ? em_damage_max : 0.0f;
        const bool emi_has_alignment_damage =
            (em_damage_ll_n > 0) && (em_damage_ll_nonzero > 0);
        const bool emi_has_damage_scores =
            (em_damage_n > 0) && (em_damage_nonzero > 0) &&
            ((em_damage_max - em_damage_min) > 1e-6f);
        // Damage-aware reassignment needs within-read alignment-level variation.
        // If DAMAGE_SCORE is constant for all alignments of a read, it only shifts
        // read-level ancientness and cancels from ref responsibilities.
        // Damage-free mapping by default: damage likelihoods must not both choose
        // the homolog and then be reused in gamma_ancient. --em-use-damage-mapping
        // restores the damage-aware mapping. When off, gamma_ancient is produced
        // post-hoc (one damage E-step against the frozen weights) below.
        em_params.use_alignment_damage_likelihood =
            em_use_damage_mapping && emi_has_alignment_damage;
        if (em_use_damage_mapping) {
            if (em_params.use_alignment_damage_likelihood) {
                em_params.use_damage = true;
            } else {
                em_params.use_damage = emi_has_damage_scores && em_damage_varies_within_read;
            }
        } else {
            em_params.use_damage = false;
        }
        if (verbose) {
            std::cerr << "  Damage term in EM: "
                      << (em_params.use_damage ? "enabled" : "disabled")
                      << "\n";
            std::cerr << "  EMI damage_score stats: mean=" << std::fixed << std::setprecision(4)
                      << em_damage_mean
                      << " std=" << std::fixed << std::setprecision(4) << em_damage_std
                      << " min=" << std::fixed << std::setprecision(4) << em_damage_min_print
                      << " max=" << std::fixed << std::setprecision(4) << em_damage_max_print
                      << " nonzero=" << em_damage_nonzero << "/" << em_damage_n << "\n";
            std::cerr << "  Within-read DAMAGE_SCORE variation: "
                      << (em_damage_varies_within_read ? "yes" : "no") << "\n";
            std::cerr << "  Alignment damage logBF stats: mean="
                      << std::fixed << std::setprecision(4) << em_llbf_mean
                      << " std=" << std::fixed << std::setprecision(4) << em_llbf_std
                      << " nonzero=" << em_damage_ll_nonzero << "/" << em_damage_ll_n << "\n";
            if (!em_params.use_damage) {
                if (!emi_has_damage_scores) {
                    std::cerr << "  Reason: EMI damage scores are all zero or constant.\n";
                } else if (!em_damage_varies_within_read) {
                    std::cerr << "  Reason: DAMAGE_SCORE is per-read constant across candidate "
                                 "references, so it cannot alter reassignment.\n";
                } else {
                    std::cerr << "  Reason: damage evidence not informative for reassignment.\n";
                }
            } else {
                if (em_params.use_alignment_damage_likelihood) {
                    std::cerr << "  Source: using alignment-level EMI damage likelihoods.\n";
                } else {
                    std::cerr << "  Source: using alignment-level EMI damage variation.\n";
                }
            }
            std::cerr << "  EM iteration diagnostics:\n";
        }

        // Post-hoc ancient responsibility is produced by the EM final pass itself when the
        // mapping is damage-free (gamma_ancient = gamma * P(ancient|aln) at frozen weights).
        em_params.posthoc_gamma_ancient = !em_use_damage_mapping && agd_read_prior;

        // Streaming (transient-gamma) EM: iterations hold no O(alignments) gamma, the final
        // pass streams per-read gamma straight into the finalization below, and
        // compact_em_map releases aln_data.alignments after decode -> steady RSS ~20 B/aln
        // (vs 20+32). Gated to the production config where per-read finalization needs only
        // ref/gamma/gamma_ancient (+ the per-read null responsibility on the tape path): no
        // calib densities, damage-free mapping. Post-hoc p_anc = sigmoid(base_a - base_m)
        // matches the legacy loop with or without DMG_LL columns (without them it is the AGD
        // read posterior itself).
        const bool can_stream =
            !use_em_calib && !em_params.use_damage && (!use_em_null || tape_path);
        const float em_keep_thr = static_cast<float>(em_min_prob);

        // ---- Tape path: build the O_DIRECT hot tape by a streaming decode (scan-2) ----
        // read_offsets + em_aln_count are already computed. The production finalize callback
        // needs only ref/gamma/gamma_ancient per read, and gamma_ancient is derived in-core
        // (posthoc: gamma*sigmoid(base_a-base_m)), so NO raw {bit,aln_len,ds} tape is needed —
        // only the 3 SoA hot streams tape_em_core consumes. Nothing O(alignments) in RAM.
        // With a null model the tape writer also keeps, per read, the best null log-density
        // over the read's kept alignments (what read_null_log_density returns on the in-RAM
        // CSR); tape_em_core scores the null column from it and returns gamma_null per read.
        std::string tape_dir_path, p_tape_ref, p_tape_ba, p_tape_bm, em_input_fp;
        std::vector<double> tape_null_lp, tape_gnull;
        if (tape_path && use_em_null) {
            tape_null_lp.assign(n_reads, -std::numeric_limits<double>::infinity());
            tape_gnull.assign(n_reads, 0.0);
        }
        struct TapeDirGuard {
            std::string dir, a, b, c;
            ~TapeDirGuard() {
                if (dir.empty()) return;
                ::unlink(a.c_str()); ::unlink(b.c_str()); ::unlink(c.c_str());
                ::rmdir(dir.c_str());
            }
        } tape_guard;
        if (tape_path) {
            if (em_tape_dir.empty()) { const char* e = std::getenv("DART_EM_TAPE_DIR"); if (e) em_tape_dir = e; }
            if (em_tape_dir.empty()) {
                std::cerr << "Error: EM tape path needs a node-local scratch dir. Pass --em-tape-dir "
                             "DIR (fast local NVMe with O_DIRECT; not /tmp, not NFS) or set "
                             "DART_EM_TAPE_DIR.\n";
                return 1;
            }
            // Read-major invariant: the tape is written in scan/CSR order with no reorder, so
            // reads must not span row groups and row groups must be globally read-ordered.
            if (!reader.is_sorted_by_read_then_score()) {
                std::cerr << "Error: EMI is not read-major (is_sorted_by_read_then_score=false); "
                             "the streaming EM tape requires read-major order.\n";
                return 1;
            }
            // Ensure the (node-local, possibly not-yet-present) base dir exists, then a
            // uniquely-named subdir. /scratch is per-node AND may be reaped during the long
            // pre-EM setup (Pass 1 + sort can run ~1h on a giant, leaving emtape empty and
            // eligible for an empty-dir reaper). A single ::mkdir then fails ENOENT because an
            // ANCESTOR is gone — so recreate the WHOLE chain (mkdir -p), best-effort.
            auto mkdir_p = [](const std::string& p) {
                for (size_t i = 1; i <= p.size(); ++i)
                    if (i == p.size() || p[i] == '/') {
                        std::string sub = p.substr(0, i);
                        ::mkdir(sub.c_str(), 0755);  // best-effort; EEXIST is fine
                    }
            };
            mkdir_p(em_tape_dir);  // recursive best-effort; EEXIST is fine
            tape_dir_path = em_tape_dir + "/dart_emtape_" + std::to_string(::getpid());
            // OPT-G: default checkpoint is a SIBLING of tape_dir_path (not a child), so the
            // TapeDirGuard that unlinks the tape does not delete it. Resume across runs needs
            // an explicit shared --em-checkpoint (this default embeds the pid).
            if (em_checkpoint_path.empty())
                em_checkpoint_path = em_tape_dir + "/dart_emtape_" + std::to_string(::getpid()) + ".ckpt";
            struct ::stat emi_st {};
            em_input_fp = emi_file + ":" +
                (::stat(emi_file.c_str(), &emi_st) == 0
                     ? std::to_string(static_cast<long long>(emi_st.st_size)) : "0");
            if (::mkdir(tape_dir_path.c_str(), 0755) != 0 && errno != EEXIST) {
                std::cerr << "Error: cannot create EM tape dir " << tape_dir_path << ": "
                          << std::strerror(errno) << "\n";
                return 1;
            }
            p_tape_ref = tape_dir_path + "/tape.ref";
            p_tape_ba  = tape_dir_path + "/tape.basea";
            p_tape_bm  = tape_dir_path + "/tape.basem";
            // Arm the RAII cleanup by assigning FIELDS directly. Do NOT construct a
            // temporary (tape_guard = TapeDirGuard{...}) — the temporary's destructor would
            // run at end of statement and rmdir the dir we just created, before we open it.
            tape_guard.dir = tape_dir_path;
            tape_guard.a = p_tape_ref;
            tape_guard.b = p_tape_ba;
            tape_guard.c = p_tape_bm;

            // Shell ref metadata — bit-identical to build_alignment_data: ref_log_len_half =
            // (float)(0.5*log(max(len,1))), consumed by length_log_penalty in the decode below.
            dart::AlignmentData shell;
            shell.num_reads = n_reads; shell.num_refs = n_refs;
            shell.ref_lengths = em_ref_lengths;
            shell.ref_log_len_half.assign(n_refs, 0.0f);
            for (uint32_t t = 0; t < n_refs; ++t) {
                const double len = static_cast<double>(std::max(shell.ref_lengths[t], 1u));
                shell.ref_log_len_half[t] = static_cast<float>(0.5 * std::log(len));
            }

            const uint64_t tape_A = static_cast<uint64_t>(em_aln_count);
            // Hot tape: 3 SoA streams (ref u32, base_a f64, base_m f64) in CSR order. Written
            // buffered at CSR offsets (per row group, contiguous under read-major), then fsync'd
            // and tail-padded to 4096 so tape_em_core's O_DIRECT aligned reads never straddle
            // EOF. Buffered writes touch only reclaimable page cache (not RssAnon), then
            // posix_fadvise(DONTNEED) drops it — peak RSS stays O(num_refs).
            int fd_ref = ::open(p_tape_ref.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            int fd_ba  = ::open(p_tape_ba.c_str(),  O_WRONLY | O_CREAT | O_TRUNC, 0644);
            int fd_bm  = ::open(p_tape_bm.c_str(),  O_WRONLY | O_CREAT | O_TRUNC, 0644);
            if (fd_ref < 0 || fd_ba < 0 || fd_bm < 0) {
                std::cerr << "Error: cannot create EM tape files in " << tape_dir_path
                          << ": " << std::strerror(errno) << "\n";
                return 1;
            }
            std::atomic<bool> tape_ok{true};
            const double eps_dec = 1e-12;
            reader.set_columns({
                dart::ColumnID::READ_IDX, dart::ColumnID::REF_IDX, dart::ColumnID::BIT_SCORE,
                dart::ColumnID::DAMAGE_SCORE, dart::ColumnID::EVALUE_LOG10,
                dart::ColumnID::DMG_LL_A, dart::ColumnID::DMG_LL_M,
                dart::ColumnID::IDENTITY_Q, dart::ColumnID::ALN_LEN,
                dart::ColumnID::TSTART, dart::ColumnID::TLEN
            });
            reader.parallel_scan([&](
                uint32_t /*rg_idx*/, uint32_t num_rows,
                const uint32_t* read_idx, const uint32_t* ref_idx,
                const float* bit_score, const float* damage_score,
                const float* evalue_log10,
                const uint16_t* /*dmg_k*/, const uint16_t* /*dmg_m*/,
                const float* dmg_ll_a, const float* dmg_ll_m,
                const uint16_t* identity_q,
                const uint16_t* aln_len, const uint16_t* /*qstart*/, const uint16_t* /*qend*/,
                const uint16_t* tstart, const uint16_t* /*tend*/,
                const uint16_t* tlen, const uint16_t* /*qlen*/,
                const uint16_t* /*mismatch*/, const uint16_t* /*gapopen*/)
            {
                if (!read_idx || !ref_idx || !bit_score || !evalue_log10 || !identity_q ||
                    !aln_len || !tstart || !tlen || !damage_score) return;
                std::vector<uint32_t> lref; std::vector<double> lba, lbm;
                lref.reserve(num_rows); lba.reserve(num_rows); lbm.reserve(num_rows);
                uint64_t base = UINT64_MAX;
                for (uint32_t i = 0; i < num_rows; ++i) {
                    const float fident = static_cast<float>(identity_q[i]) / 65535.0f;
                    const float bits = bit_score[i];
                    if (fident < aln_min_identity) continue;
                    if (bits < aln_min_bits) continue;
                    if (evalue_log10[i] > aln_max_evalue_log10) continue;
                    const uint32_t ridx = read_idx[i];
                    if (base == UINT64_MAX) base = read_offsets[ridx];  // read-major: rg is contiguous
                    const uint32_t tstart_0 = (tstart[i] > 0) ? static_cast<uint32_t>(tstart[i] - 1) : 0u;
                    const uint32_t aln_start = std::min(tstart_0, 65535u);
                    const uint32_t aln_end = std::min(tstart_0 + static_cast<uint32_t>(aln_len[i]), 65535u);
                    const uint16_t aln_len_em = (aln_end > aln_start)
                        ? static_cast<uint16_t>(std::min(aln_end - aln_start, 65535u)) : 0;
                    if (!tape_null_lp.empty()) {
                        // Same function as the in-RAM path, one alignment at a time; the
                        // per-read max is order-exact, so scan order cannot change it.
                        dart::Alignment na{};
                        na.bit_score = bits;
                        na.aln_len = aln_len_em;
                        const double lp = dart::read_null_log_density(em_params, &na, 1);
                        double& slot = tape_null_lp[ridx];   // read-major: this rg owns ridx
                        if (lp > slot) slot = lp;
                    }
                    const uint32_t ref = ref_idx[i];
                    // Decode == tape_decode_read with score_scale=1 (normalize_scores_per_query is
                    // always false here): bc = length_log_penalty(shell) + coverage(0) + lambda_b*bit,
                    // split by em_base_pair (damage-free here: the tape path excludes
                    // --em-use-damage-mapping).
                    const double bc = dart::length_log_penalty(shell, em_params, ref)
                                    + dart::coverage_log_penalty(em_params, ref)
                                    + dart::real_log_density(em_params, static_cast<double>(bits), aln_len_em, 1.0);
                    const double pr = std::clamp(
                        static_cast<double>(read_dmg(damage_score, i, ridx)), eps_dec, 1.0 - eps_dec);
                    const double la = dmg_ll_a ? static_cast<double>(dmg_ll_a[i]) : 0.0;
                    const double lm = dmg_ll_m ? static_cast<double>(dmg_ll_m[i]) : 0.0;
                    double ba, bm;
                    dart::em_base_pair(em_params, bc, pr, la, lm, ba, bm);
                    lref.push_back(ref);
                    lba.push_back(ba);
                    lbm.push_back(bm);
                }
                if (lref.empty()) return;
                const off_t o_ref = static_cast<off_t>(base) * static_cast<off_t>(sizeof(uint32_t));
                const off_t o_d   = static_cast<off_t>(base) * static_cast<off_t>(sizeof(double));
                if (::pwrite(fd_ref, lref.data(), lref.size() * sizeof(uint32_t), o_ref) < 0) tape_ok = false;
                if (::pwrite(fd_ba,  lba.data(),  lba.size()  * sizeof(double),   o_d)   < 0) tape_ok = false;
                if (::pwrite(fd_bm,  lbm.data(),  lbm.size()  * sizeof(double),   o_d)   < 0) tape_ok = false;
            });
            auto pad4096 = [](int fd, uint64_t nbytes) {
                uint64_t padded = (nbytes + 4095) & ~uint64_t(4095);
                return ::ftruncate(fd, static_cast<off_t>(padded)) == 0;
            };
            bool ok = tape_ok.load();
            ok = pad4096(fd_ref, tape_A * sizeof(uint32_t)) && ok;
            ok = pad4096(fd_ba,  tape_A * sizeof(double))   && ok;
            ok = pad4096(fd_bm,  tape_A * sizeof(double))   && ok;
            ::fsync(fd_ref); ::fsync(fd_ba); ::fsync(fd_bm);
#ifdef __linux__
            ::posix_fadvise(fd_ref, 0, 0, POSIX_FADV_DONTNEED);
            ::posix_fadvise(fd_ba,  0, 0, POSIX_FADV_DONTNEED);
            ::posix_fadvise(fd_bm,  0, 0, POSIX_FADV_DONTNEED);
#endif
            ::close(fd_ref); ::close(fd_ba); ::close(fd_bm);
            if (!ok) { std::cerr << "Error: EM tape write failed\n"; return 1; }
            reader.flush_pages();
            if (verbose) {
                const double bh_gb = static_cast<double>(n_reads) *
                    static_cast<double>(sizeof(best_hits[0])) / (1024.0 * 1024.0 * 1024.0);
                std::cerr << "  EM tape: A=" << tape_A << " alns, R=" << n_reads
                          << " reads -> " << tape_dir_path
                          << " | O(R) best_hits ~" << std::fixed << std::setprecision(2) << bh_gb
                          << " GB (" << sizeof(best_hits[0]) << " B/read) + read_offsets "
                          << std::setprecision(2)
                          << (static_cast<double>(n_reads) * 8.0 / (1024.0*1024.0*1024.0)) << " GB\n";
            }
        }

        const auto em_start = std::chrono::steady_clock::now();

        // compact_em_map (fusion dual-b8247d43): the SOLE EM reassignment solver. Decode-once
        // CSR + fused per-ref sufficient stats + deterministic-lane M-step -> bitwise-identical
        // across thread counts. Streaming form additionally makes gamma EM-transient.
        if (can_stream) {
            // Per-read finalization, streamed from the EM final pass. gamma is transient;
            // only O(reads) results survive. Mirrors the legacy per-read loop below for the
            // production config (null off, calib off, audit off). On the tape path (OPT-F)
            // this runs CONCURRENTLY across 32 read lanes: every write below is indexed by r
            // (em_read_results[r], mm_pool[r]) and damage_em_active is atomic, so no lock.
            auto finalize_cb = [&](uint32_t r, uint32_t deg, const uint32_t* ref,
                                   const double* gamma, const double* gamma_anc) {
                uint32_t em_best_ref = 0;
                double best_g = -1.0;
                float top1 = 0.0f, top2 = 0.0f;
                double sum_sq = 0.0;
                // Tape path with a null: tape_em_core stores gamma_null[r] before calling us.
                // Margin/neff use the real responsibilities renormalized by 1/(1-gamma_null),
                // exactly as the legacy loop below; best-hit gamma and the pool stay raw.
                const double gnull = tape_gnull.empty() ? 0.0 : tape_gnull[r];
                const double real_norm = (gnull < 1.0 - 1e-9) ? (1.0 / (1.0 - gnull)) : 0.0;
                for (uint32_t j = 0; j < deg; ++j) {
                    const double g = gamma[j] * real_norm;
                    if (gamma[j] > best_g) { best_g = gamma[j]; em_best_ref = ref[j]; }
                    sum_sq += g * g;
                    if (g > static_cast<double>(top1)) { top2 = top1; top1 = static_cast<float>(g); }
                    else if (g > static_cast<double>(top2)) { top2 = static_cast<float>(g); }
                }
                float best_hit_gamma = 0.0f, best_hit_gamma_ancient = 0.0f;
                bool bitscore_best_matches_em = false;
                bool em_keep = false;
                if (deg > 0 && best_hits[r].has) {
                    const uint32_t best_hit_ref = best_hits[r].ref_idx;
                    bitscore_best_matches_em = (best_hit_ref == em_best_ref);
                    for (uint32_t j = 0; j < deg; ++j) {
                        if (ref[j] == best_hit_ref) {
                            best_hit_gamma = static_cast<float>(gamma[j]);
                            if (gamma_anc) {
                                best_hit_gamma_ancient = static_cast<float>(gamma_anc[j]);
                                damage_em_active = true;
                            }
                            break;
                        }
                    }
                    em_keep = (best_hit_gamma >= em_keep_thr);
                }
                if (deg > 0 && emit_multimap && r < mm_pool.size() && !mm_pool[r].empty()) {
                    for (uint32_t j = 0; j < deg; ++j) {
                        const uint32_t jref = ref[j];
                        for (auto& c : mm_pool[r]) {
                            if (c.ref_idx == jref) {
                                c.gamma = static_cast<float>(gamma[j]);
                                c.gamma_ancient = gamma_anc ? static_cast<float>(gamma_anc[j]) : 0.0f;
                                break;
                            }
                        }
                    }
                }
                const float em_margin = std::max(0.0f, top1 - top2);
                const float em_ref_neff = (sum_sq > 0.0) ? static_cast<float>(1.0 / sum_sq) : 0.0f;
                uint8_t reason = 0;
                if (em_margin < static_cast<float>(em_min_margin)) reason |= 0x2;
                if (em_ref_neff > static_cast<float>(em_max_neff)) reason |= 0x4;
                // Bit 0x1 on this (uncalibrated) path comes from null_p, set after the EM.
                em_read_results[r] = {
                    best_hit_gamma, best_hit_gamma_ancient, bitscore_best_matches_em,
                    em_keep, em_keep_thr, em_margin, em_ref_neff,
                    static_cast<float>(gnull), reason
                };
            };
            if (tape_path) {
                // Stream-solve directly off the O_DIRECT hot tape: base_a/base_m live on
                // scratch NVMe, not RAM, so peak RSS is O(num_refs) + fixed I/O buffers,
                // INDEPENDENT of alignment count. With a null model the per-read null
                // log-density gathered by the tape writer feeds the null column.
                dart::TapeEMResult res = dart::tape_em_core(
                    p_tape_ref, p_tape_ba, p_tape_bm, read_offsets.data(),
                    n_reads, n_refs, static_cast<uint64_t>(em_aln_count), em_params,
                    tape_null_lp.empty() ? nullptr : tape_null_lp.data(),
                    tape_gnull.empty() ? nullptr : tape_gnull.data(),
                    finalize_cb, em_checkpoint_path, em_resume, em_input_fp);
                em_state.weights = std::move(res.weights);
                em_state.pi = res.pi;
                em_state.w_null = res.w_null;
                em_state.log_likelihood = res.log_likelihood;
                em_state.iterations = res.iterations;
                if (!tape_gnull.empty()) em_state.gamma_null = std::move(tape_gnull);
            } else {
                em_state = dart::compact_em_map(aln_data, em_params, finalize_cb);
            }
        } else {
            em_state = dart::compact_em_map(aln_data, em_params);
        }

        const auto em_end = std::chrono::steady_clock::now();

        if (verbose) {
            std::cerr << "  Converged in " << em_state.iterations << " iterations\n";
            std::cerr << "  Log-likelihood: " << std::fixed << std::setprecision(2)
                      << em_state.log_likelihood << "\n";
            if (em_params.use_damage && !em_params.use_alignment_damage_likelihood) {
                std::cerr << "  Estimated pi (ancient fraction): "
                          << std::fixed << std::setprecision(4) << em_state.pi << "\n";
            } else if (em_params.use_alignment_damage_likelihood) {
                std::cerr << "  Estimated pi (ancient fraction): n/a (per-read priors used)\n";
            }
            std::cerr << "  Runtime: "
                      << dart::log_utils::format_elapsed(em_start, em_end) << "\n";
        }

        // The remaining per-read finalization runs only on the NON-streaming path. The
        // streaming path already produced em_read_results / mm_pool joins inside the EM
        // final-pass callback above (and released aln_data.alignments).
        if (!can_stream) {
        // Post-hoc gamma_ancient when the mapping was damage-free (default). The
        // mapping marginal em_state.gamma is left BIT-IDENTICAL; gamma_ancient is the
        // factorized product gamma_map_j * P(ancient | alignment_j), the ancient
        // posterior evaluated against the FROZEN weights only. This mirrors the
        // streaming path exactly. Crucially it is NOT a jointly-normalized damage
        // E-step: that yields a damage-tilted marginal when per-alignment damage
        // varies, which would re-couple mapping and damage in the ancient channel.
        // Damage never re-enters the mapping/M-step; pi stays at its 0.10 prior
        // because damage never drove the mapping loop.
        if (em_params.posthoc_gamma_ancient) {
            const double eps = 1e-12;
            // p_anc = the AGD read posterior; the alignment damage likelihoods would count
            // its terminal C->T sites a second time (see em_base_pair).
            const size_t A = aln_data.alignments.size();
            em_state.gamma_ancient.assign(A, 0.0);
            for (size_t j = 0; j < A; ++j) {
                const auto& aln = aln_data.alignments[j];
                const double p_anc = std::clamp(static_cast<double>(aln.damage_score), eps, 1.0 - eps);
                em_state.gamma_ancient[j] = em_state.gamma[j] * p_anc;
            }
        }

        // Extract best assignments and build per-read lookup
        auto best = dart::reassign_reads(aln_data, em_state, 0.0);


        // Build per-read lookup: extract gamma for the best-hit target used in summaries,
        // while tracking whether EM would have reassigned this read elsewhere.
        for (uint32_t r = 0; r < n_reads; ++r) {
            uint32_t start = aln_data.read_offsets[r];
            uint32_t end = aln_data.read_offsets[r + 1];

            float best_hit_gamma = 0.0f;
            float best_hit_gamma_ancient = 0.0f;
            uint32_t em_best_ref = best[r].first;
            bool bitscore_best_matches_em = false;
            float top1 = 0.0f, top2 = 0.0f;
            double sum_sq = 0.0;
            // EM-assigned candidate's bit / aligned length (for the audit densities).
            double em_best_bit = 0.0; uint16_t em_best_aln = 0;

            // Null responsibility for this read (0 when the null column is off).
            const double gnull = (em_params.use_null && r < em_state.gamma_null.size())
                ? em_state.gamma_null[r] : 0.0;
            // Renormalize real responsibilities by 1/(1-gamma_null) so margin/neff are
            // computed on the REAL simplex — avoids double-penalizing a read that already
            // pays for its null mass. Guard the degenerate gamma_null ≈ 1 case.
            const double real_norm = (gnull < 1.0 - 1e-9) ? (1.0 / (1.0 - gnull)) : 0.0;

            for (uint32_t j = start; j < end; ++j) {
                const double g = em_state.gamma[j] * real_norm;
                sum_sq += g * g;
                if (g > static_cast<double>(top1)) {
                    top2 = top1;
                    top1 = static_cast<float>(g);
                } else if (g > static_cast<double>(top2)) {
                    top2 = static_cast<float>(g);
                }
                if (aln_data.alignments[j].ref_idx == em_best_ref) {
                    em_best_bit = static_cast<double>(aln_data.alignments[j].bit_score);
                    em_best_aln = aln_data.alignments[j].aln_len;
                }
            }
            if (emit_multimap && r < mm_pool.size() && !mm_pool[r].empty()) {  // join gamma onto the pool
                for (uint32_t j = start; j < end; ++j) {
                    const uint32_t jref = aln_data.alignments[j].ref_idx;
                    for (auto& c : mm_pool[r]) {
                        if (c.ref_idx == jref) {
                            c.gamma = static_cast<float>(em_state.gamma[j]);
                            c.gamma_ancient = em_state.gamma_ancient.empty() ? 0.0f
                                : static_cast<float>(em_state.gamma_ancient[j]);
                            break;
                        }
                    }
                }
            }

            // Check if the best-by-bitscore hit matches the EM-assigned best
            const auto& bh = best_hits[r];
            float keep_thr = static_cast<float>(em_min_prob);
            bool em_keep = true;
            if (bh.has) {
                const uint32_t best_hit_ref = bh.ref_idx;
                bitscore_best_matches_em = (best_hit_ref == em_best_ref);
                for (uint32_t j = start; j < end; ++j) {
                    if (aln_data.alignments[j].ref_idx == best_hit_ref) {
                        best_hit_gamma = static_cast<float>(em_state.gamma[j]);
                        // gamma_ancient is filled either by the damage-aware EM
                        // (--em-use-damage-mapping) or the post-hoc frozen-weight pass
                        // above; a non-empty vector is the single source of truth.
                        if (!em_state.gamma_ancient.empty()) {
                            best_hit_gamma_ancient = static_cast<float>(em_state.gamma_ancient[j]);
                            damage_em_active = true;
                        }
                        break;
                    }
                }
                em_keep = (best_hit_gamma >= keep_thr);
            } else {
                em_keep = false;
            }

            const float em_margin = std::max(0.0f, top1 - top2);
            const float em_ref_neff = (sum_sq > 0.0) ? static_cast<float>(1.0 / sum_sq) : 0.0f;

            // Per-read identifiability densities for the best EM candidate (frozen calib).
            double log_ft = std::numeric_limits<double>::quiet_NaN();
            double log_f0 = std::numeric_limits<double>::quiet_NaN();
            if (em_params.use_calib_cont) {
                log_ft = dart::log_f_true_cont(em_params, em_best_bit, em_best_aln);
                log_f0 = dart::log_f0_cont(em_params, em_best_bit, em_best_aln);
            } else if (em_params.use_calib && !em_params.true_bins.empty()) {
                log_ft = dart::log_f_true(em_params, em_best_bit, em_best_aln);
                int nb = dart::null_length_bin(em_best_aln);
                if (nb < static_cast<int>(em_params.null_mu.size()))
                    log_f0 = dart::log_p_null_gumbel(em_best_bit,
                                 em_params.null_mu[nb], em_params.null_beta[nb]);
            }
            const bool have_llr = std::isfinite(log_ft) && std::isfinite(log_f0);
            const double llr = have_llr ? (log_ft - log_f0) : 0.0;

            uint8_t reason = 0;
            if (em_margin < static_cast<float>(em_min_margin)) reason |= 0x2;  // low_margin
            if (em_ref_neff > static_cast<float>(em_max_neff)) reason |= 0x4;  // high_neff

            // 3-way decision (sol). With calibrated densities the floor is on the LLR
            // (identifiability), which — unlike gamma_null — does not saturate to 1 via
            // the EM w_null dynamics, so information-poor reads ABSTAIN rather than being
            // mislabeled contamination. Without an LLR the calibrated path keeps
            // gamma_null > tau. The uncalibrated null gets bit 0x1 from null_p after the EM
            // (gamma_null cannot reach tau there: the real term is unnormalised).
            if (em_params.use_null && use_em_calib && end > start) {
                if (have_llr) {
                    // Abstained/rejected reads leave the assignment/auth pools (em_keep=false).
                    if (llr <= em_reject_llr)     { reason |= 0x1; em_keep = false; }
                    else if (llr < em_keep_llr)   { reason |= 0x8; em_keep = false; }
                } else if (gnull > em_null_tau)   { reason |= 0x1; }
            }


            em_read_results[r] = {
                best_hit_gamma,
                best_hit_gamma_ancient,
                bitscore_best_matches_em,
                em_keep,
                keep_thr,
                em_margin,
                em_ref_neff,
                static_cast<float>(gnull),
                reason
            };
        }

        }  // end if(!can_stream) legacy per-read finalization

        // Chance probability of each read's best hit under the decoy null, the same on the
        // tape and in-RAM paths. On the uncalibrated path it alone sets reject bit 0x1.
        if (use_em_null && null_p_n > 0) {
            null_p_ran = true;
            for (uint32_t r = 0; r < n_reads; ++r) {
                const auto& bh = best_hits[r];
                if (!bh.has) continue;
                const double p = dart::null_chance_p(em_null_model, null_p_scale, null_p_n, bh.bits, bh.aln_len);
                auto& er = em_read_results[r];
                er.null_p = static_cast<float>(p);
                const int b = dart::null_length_bin(bh.aln_len);
                ++null_bin_reads[b];
                if (!use_em_calib && p > null_p_max) {
                    er.em_reject_reason |= 0x1;
                    ++null_bin_flagged[b];
                }
            }
        }

        if (verbose) {
            // Filtering stats
            size_t multi = 0, unique = 0;
            size_t reassigned = 0;
            size_t kept_total = 0, kept_multi = 0;
            size_t fail_abs = 0;
            double margin_sum = 0.0;
            double neff_sum = 0.0;
            size_t margin_n = 0;
            for (uint32_t r = 0; r < n_reads; ++r) {
                uint32_t deg = static_cast<uint32_t>(read_offsets[r + 1] - read_offsets[r]);
                if (deg > 1) multi++;
                else unique++;
                const auto& er = em_read_results[r];
                if (!er.best_hit) reassigned++;
                if (er.em_keep) {
                    kept_total++;
                    if (deg > 1) kept_multi++;
                } else {
                    if (er.gamma < static_cast<float>(em_min_prob)) fail_abs++;
                }
                if (deg > 0) {
                    margin_sum += er.em_margin;
                    neff_sum += er.em_ref_neff;
                    margin_n++;
                }
            }
            std::cerr << "  Unique mappers: " << unique << "\n";
            std::cerr << "  Multi-mappers: " << multi << "\n\n";
            std::cerr << "  EM keep filter (best-hit target): kept "
                      << kept_total << "/" << n_reads
                      << " (" << std::fixed << std::setprecision(2)
                      << (n_reads > 0 ? 100.0 * static_cast<double>(kept_total) / static_cast<double>(n_reads) : 0.0)
                      << "%), kept multi " << kept_multi << "/" << multi << "\n";
            std::cerr << "    failed abs(min_prob): " << fail_abs << "\n";
            if (margin_n > 0) {
                const double reassigned_pct =
                    100.0 * static_cast<double>(reassigned) / static_cast<double>(margin_n);
                std::cerr << "  Reassigned vs bitscore-best: " << std::fixed << std::setprecision(2)
                          << reassigned_pct << "%\n";
                std::cerr << "  Mean EM margin (top1-top2): " << std::fixed << std::setprecision(4)
                          << (margin_sum / static_cast<double>(margin_n)) << "\n";
                std::cerr << "  Mean EM effective refs/read: " << std::fixed << std::setprecision(4)
                          << (neff_sum / static_cast<double>(margin_n)) << "\n\n";
            }
        }
    }

    // Stage covered alignments in bounded sort runs. Posteriors are filled after
    // scoring, while auth picks are retained separately in read order.
    const bool emit_dcx = !dcx_file.empty();
    if (emit_dcx && !use_em) {
        std::cerr << "Error: --dcx requires EM (grain-A gamma/em_keep/posterior); "
                     "remove --no-em\n";
        return 1;
    }
    // .dcx records go to an external sort keyed by target rank (target name byte order),
    // so the hot loop compares integers, never names. Each thread stages into its own sink.
    // Auth picks (the em_keep record of each read) are collected per thread and reduced
    // after the loop.
    struct AuthPick { uint32_t read, rank, ref; uint16_t gamma, prior; };
    std::vector<uint32_t> ref_rank;
    std::unique_ptr<dcx::ExternalSort> dcx_sort;
    std::mutex dcx_pick_mutex;
    std::vector<std::vector<AuthPick>> dcx_pick_lists;
    std::atomic<bool> dcx_stage_failed{false};
    std::exception_ptr dcx_stage_error;
    if (emit_dcx) {
        // Decode the reference names in bulk, one arena per dictionary block, and rank them.
        const auto& dict = reader.ref_dict();
        const uint32_t bs = dict.dict_block_size();
        const uint32_t nblk = n_refs ? (n_refs + bs - 1) / bs : 0;
        std::vector<std::vector<char>> arena(nblk);
        std::vector<std::vector<uint64_t>> ends(nblk);
        #pragma omp parallel for schedule(dynamic, 16)
        for (uint32_t blk = 0; blk < nblk; ++blk) dict.decode_block_bytes(blk, arena[blk], ends[blk]);
        std::vector<std::string_view> names(n_refs);
        for (uint32_t blk = 0; blk < nblk; ++blk)
            for (size_t k = 0; k < ends[blk].size(); ++k) {
                const uint64_t from = k ? ends[blk][k - 1] : 0;
                names[size_t(blk) * bs + k] = std::string_view(arena[blk].data() + from, ends[blk][k] - from);
            }
        std::vector<uint32_t> order(n_refs);
        std::iota(order.begin(), order.end(), 0u);
        auto by_name = [&](uint32_t a, uint32_t b) { return names[a] < names[b]; };
#ifdef _OPENMP
        __gnu_parallel::sort(order.begin(), order.end(), by_name);
#else
        std::sort(order.begin(), order.end(), by_name);
#endif
        ref_rank.assign(n_refs, 0);
        uint32_t rank = 0;
        for (size_t k = 0; k < order.size(); ++k) {
            if (k > 0 && names[order[k]] != names[order[k - 1]]) ++rank;
            ref_rank[order[k]] = rank;
        }
        std::string run_dir = em_tape_dir;
        if (run_dir.empty()) {
            if (const char* dir = std::getenv("DART_EM_TAPE_DIR")) run_dir = dir;
        }
        size_t budget = size_t(8) << 30;
        if (const char* value = std::getenv("DART_DCX_SORT_BYTES")) budget = std::stoull(value);
        size_t fail_run = std::numeric_limits<size_t>::max();
        if (const char* value = std::getenv("DART_DCX_TEST_FAIL_RUN")) fail_run = std::stoull(value);
        dcx_sort = std::make_unique<dcx::ExternalSort>(run_dir, budget, runtime_threads, fail_run);
    }
    // Names for a list of reference indices, decoding each dictionary block once.
    auto ref_names_of = [&](const std::vector<uint32_t>& refs) {
        const auto& dict = reader.ref_dict();
        const uint32_t bs = dict.dict_block_size();
        std::vector<uint32_t> blocks;
        for (uint32_t r : refs) blocks.push_back(r / bs);
        std::sort(blocks.begin(), blocks.end());
        blocks.erase(std::unique(blocks.begin(), blocks.end()), blocks.end());
        std::vector<std::vector<char>> arena(blocks.size());
        std::vector<std::vector<uint64_t>> ends(blocks.size());
        #pragma omp parallel for schedule(dynamic, 16)
        for (size_t k = 0; k < blocks.size(); ++k) dict.decode_block_bytes(blocks[k], arena[k], ends[k]);
        std::vector<std::string> out;
        out.reserve(refs.size());
        for (uint32_t r : refs) {
            const size_t k = size_t(std::lower_bound(blocks.begin(), blocks.end(), r / bs) - blocks.begin());
            const size_t i = r % bs;
            const uint64_t from = i ? ends[k][i - 1] : 0;
            out.emplace_back(arena[k].data() + from, ends[k][i] - from);
        }
        return out;
    };
    auto stage_dcx = [&](const DcxReadRec& rec, std::vector<AuthPick>& picks) {
        if (dcx_stage_failed.load(std::memory_order_relaxed)) return;
        try {
            const uint32_t rank = ref_rank[rec.ref_idx];
            if (rec.em_keep)
                picks.push_back({rec.read_idx, rank, rec.ref_idx,
                                 dcx::enc_unit(rec.gamma), dcx::enc_unit(rec.p_read)});
            dcx_sort->local().append(rec, rank);
        } catch (...) {
            std::lock_guard<std::mutex> lock(dcx_pick_mutex);
            if (!dcx_stage_failed.exchange(true)) dcx_stage_error = std::current_exception();
        }
    };

    // ── Annotation pass (runs after EM) ──────────────────────────────────────
    // Running EM first means summaries (~8 GB compact) never coexist with EM
    // state in memory, eliminating ~25-42 GB of peak RSS vs the old order.
    {
        const auto annotate_start = std::chrono::steady_clock::now();
        std::vector<std::pair<uint32_t, ProteinDamageSummary>> all_results;

        // Open per-site and corrected-protein output files before the parallel scan.
        std::ofstream sites_out_f;
        size_t written_sites = 0;
        if (!sites_file.empty()) {
            sites_out_f.open(sites_file);
            if (!sites_out_f.is_open()) {
                std::cerr << "Error: Cannot open sites file: " << sites_file << "\n";
                return 1;
            }
            sites_out_f << "query_id\ttarget_id\tquery_pos\ttarget_pos\t"
                           "target_aa\tquery_aa\tdamage_class\tconfidence\t"
                           "dist_5prime\tdist_3prime\tp_damage\n";
        }

        // --emit-multimap: finalize the soft top-K pool by gamma and REBUILD best_rows_by_rg
        // from the kept candidates. best_rows was built pre-EM from the argmax; the gamma
        // top-K is only known now. mm_row_cand maps (rg,row)->candidate for the annotation
        // loop. best_rows / selected_row_groups are consumed ONLY in the scan below.
        std::unordered_map<uint64_t, const MmCand*> mm_row_cand;
        if (emit_multimap) {
            for (auto& rg : best_rows_by_rg) rg.clear();
            selected_row_groups.clear();
            selected_best_hits = 0;
            std::vector<uint8_t> rg_seen(best_rows_by_rg.size(), 0);
            const float floorv = static_cast<float>(multimap_floor);
            const int K = std::max(1, multimap_k);
            for (uint32_t r = 0; r < n_reads; ++r) {
                auto& pool = mm_pool[r];
                if (pool.empty()) continue;
                std::sort(pool.begin(), pool.end(),
                          [](const MmCand& a, const MmCand& b) { return a.gamma > b.gamma; });
                std::vector<MmCand> kept;
                std::unordered_set<uint32_t> seen_ref;
                for (auto& c : pool) {
                    if (c.gamma <= floorv) break;                       // gamma-sorted desc
                    if (static_cast<int>(kept.size()) >= K) break;
                    if (!seen_ref.insert(c.ref_idx).second) continue;   // one HSP per ref
                    c.candidate_rank = static_cast<uint8_t>(kept.size());
                    kept.push_back(c);
                }
                pool.swap(kept);
                for (auto& c : pool) {
                    if (c.best_rg >= best_rows_by_rg.size()) continue;
                    if (!rg_seen[c.best_rg]) { rg_seen[c.best_rg] = 1; selected_row_groups.push_back(c.best_rg); }
                    best_rows_by_rg[c.best_rg].emplace_back(c.best_row, r);
                    mm_row_cand[(static_cast<uint64_t>(c.best_rg) << 32) | c.best_row] = &c;
                    selected_best_hits++;
                }
            }
            for (auto& rows : best_rows_by_rg) std::sort(rows.begin(), rows.end());
            std::sort(selected_row_groups.begin(), selected_row_groups.end());
        }

        if (selected_best_hits > 0) {
            // --emit-multimap re-fetches the slim pool's column-derived metadata here (it was
            // dropped from MmCand to bound peak RSS), so the selected scan must load those
            // numeric columns alongside QALN/TALN. The non-multimap path keeps the lean
            // QALN/TALN set (best_hits still carries full metadata).
            if (emit_multimap)
                reader.set_columns({dart::ColumnID::QALN, dart::ColumnID::TALN,
                                    dart::ColumnID::DAMAGE_SCORE, dart::ColumnID::EVALUE_LOG10,
                                    dart::ColumnID::IDENTITY_Q, dart::ColumnID::QSTART,
                                    dart::ColumnID::TSTART, dart::ColumnID::TLEN,
                                    dart::ColumnID::QLEN});
            else
                reader.set_columns({dart::ColumnID::QALN, dart::ColumnID::TALN});

            reader.parallel_scan_selected(selected_row_groups, [&](
                uint32_t rg_idx, uint32_t num_rows,
                const uint32_t*, const uint32_t*,
                const float*, const float* mm_damage_score, const float* mm_evalue_log10,
                const uint16_t*, const uint16_t*,
                const float*, const float*,
                const uint16_t* mm_identity_q, const uint16_t*,
                const uint16_t* mm_qstart, const uint16_t*,
                const uint16_t* mm_tstart, const uint16_t*,
                const uint16_t* mm_tlen, const uint16_t* mm_qlen,
                const uint16_t*, const uint16_t*)
            {
                const auto& rows = best_rows_by_rg[rg_idx];
                if (rows.empty()) return;
                if (emit_multimap && (!mm_damage_score || !mm_evalue_log10 || !mm_identity_q ||
                                      !mm_qstart || !mm_tstart || !mm_tlen || !mm_qlen)) return;

                // Parallel annotation within this row group.
                // string_view references into the mmap are valid until MADV_DONTNEED fires
                // after this callback returns (handled by parallel_scan_selected).
                #pragma omp parallel
                {
                    std::vector<std::pair<uint32_t, ProteinDamageSummary>> local_results;
                    std::ostringstream local_sites;
                    std::vector<AuthPick> local_picks;

                    #pragma omp for schedule(dynamic, 64) nowait
                    for (size_t i = 0; i < rows.size(); ++i) {
                        const auto& [row_in_group, ridx] = rows[i];
                        if (row_in_group >= num_rows) continue;
                        BestHitData mm_cand;  // --emit-multimap: THIS candidate's own metadata
                        BestHitData* hitp = &best_hits[ridx];
                        if (emit_multimap) {
                            auto itc = mm_row_cand.find(
                                (static_cast<uint64_t>(rg_idx) << 32) | row_in_group);
                            if (itc == mm_row_cand.end()) continue;
                            const MmCand& mc = *itc->second;
                            // Reconstruct the column-derived fields dropped from the slim pool,
                            // using the SAME expressions as the pre-EM fill so emitted bytes are
                            // identical. row_in_group indexes this row group's column arrays.
                            mm_cand.has = true;
                            mm_cand.bits = mc.bits;
                            mm_cand.ref_idx = mc.ref_idx;
                            mm_cand.best_rg = mc.best_rg;
                            mm_cand.best_row = mc.best_row;
                            mm_cand.gamma = mc.gamma;
                            mm_cand.gamma_ancient = mc.gamma_ancient;
                            mm_cand.candidate_rank = mc.candidate_rank;
                            mm_cand.damage_score = read_dmg(mm_damage_score, row_in_group, ridx);
                            mm_cand.fident = static_cast<float>(mm_identity_q[row_in_group]) / 65535.0f;
                            mm_cand.evalue = std::pow(10.0f, mm_evalue_log10[row_in_group]);
                            mm_cand.qstart_0 = (mm_qstart[row_in_group] > 0)
                                ? static_cast<uint32_t>(mm_qstart[row_in_group] - 1) : 0u;
                            mm_cand.tstart_0 = (mm_tstart[row_in_group] > 0)
                                ? static_cast<uint32_t>(mm_tstart[row_in_group] - 1) : 0u;
                            mm_cand.qlen = static_cast<uint32_t>(mm_qlen[row_in_group]);
                            mm_cand.tlen = static_cast<uint32_t>(mm_tlen[row_in_group]);
                            hitp = &mm_cand;
                        }
                        BestHitData& hit = *hitp;
                        if (!hit.has) continue;

                        std::string_view qaln_sv = reader.get_qaln(rg_idx, row_in_group);
                        std::string_view taln_sv = reader.get_taln(rg_idx, row_in_group);

                        // Compute aln_len_on_target (previously done in the old Pass 2).
                        uint32_t aln_len_target = 0;
                        for (char c : taln_sv) if (c != '-') ++aln_len_target;
                        hit.aln_len_on_target = aln_len_target;

                        std::string query_id(reader.read_name(ridx));
                        std::string target_id(reader.ref_name(hit.ref_idx));


                        // Annotate (target_id not passed — not stored in struct).
                        auto summary = annotate_alignment(
                            query_id, qaln_sv, taln_sv,
                            hit.qstart_0, hit.tstart_0, hit.qlen,
                            hit.evalue, hit.bits, hit.fident, d_max, lambda,
                            em_state.pi);

                        per_read_aa_decay[ridx] = summary.aa_sum_exp_decay;
                        // Capture opportunity count for ALL reads (clean incl.) before the
                        // total_mismatches>0 gate discards clean-read summaries below.
                        per_read_m_opp[ridx] = summary.aa_m_opportunities;

                        if (summary.total_mismatches > 0) {
                            summary.read_idx = ridx;
                            summary.qlen = hit.qlen;
                            summary.qstart_0 = hit.qstart_0;
                            summary.tstart_0 = hit.tstart_0;
                            summary.tlen = hit.tlen;
                            summary.delta_bits = hit.bits - hit.bits_second;
                            summary.bpa = (summary.alnlen > 0)
                                ? hit.bits / static_cast<float>(summary.alnlen) : 0.0f;
                            summary.length_bin = static_cast<uint8_t>(
                                dart::LengthBinStats::get_bin(hit.qlen));
                            summary.p_read = std::clamp(hit.damage_score, 0.0f, 1.0f);
                            summary.has_p_read = true;
                            if (damage_index) {
                                if (const AgdRecord* rec = damage_index->find(query_id)) {
                                    auto syn_result = detect_synonymous_damage(
                                        *rec, damage_index->d_max(), damage_index->lambda());
                                    summary.syn_5prime = syn_result.synonymous_5prime;
                                    summary.syn_3prime = syn_result.synonymous_3prime;
                                    summary.has_syn_data = true;
                                    summary.p_read = static_cast<float>(rec->p_damaged_q) / 255.0f;
                                    summary.has_p_read = true;
                                }
                            }


                            // Write per-site output inline (before clearing sites).
                            if (sites_out_f.is_open()) {
                                for (const auto& site : summary.sites) {
                                    local_sites << query_id << "\t" << target_id << "\t"
                                        << site.query_pos << "\t" << site.target_pos << "\t"
                                        << site.target_aa << "\t" << site.query_aa << "\t"
                                        << (site.damage_class == 'C' ? "CT" : "GA") << "\t"
                                        << (site.high_confidence ? "high" : "medium") << "\t"
                                        << site.dist_5prime << "\t" << site.dist_3prime << "\t"
                                        << std::fixed << std::setprecision(4)
                                        << site.p_damage << "\n";
                                }
                            }

                            // Aggregate sites into compact scalars (replaces sites loop in
                            // ProteinAggInput and protein aggregate building).
                            constexpr size_t TERM_COD = 3;  // codons from terminus
                            for (const auto& site : summary.sites) {
                                summary.info_sum += site.p_damage;
                                if (site.damage_class == 'C' && site.dist_5prime / 3 < TERM_COD)
                                    ++summary.terminal_5;
                                else if (site.damage_class == 'G' && site.dist_3prime / 3 < TERM_COD)
                                    ++summary.terminal_3;
                            }
                        }

                        // Stage the .dcx grain-A/B record for EVERY best-hit alignment —
                        // INCLUDING zero-substitution (clean) reads AND em_keep=0 reads.
                        // em_keep is a SOFT weight (fusion decision D / M5), not a filter:
                        // the emit must not hard-gate on it or the producer contradicts the
                        // consumer (saav_observation carries em_keep as a no-WHERE column).
                        // Two bugs were fixed here: (1) this block was nested inside the
                        // total_mismatches>0 guard, dropping every clean read; (2) it was
                        // gated on em_keep, dropping every em_keep=0 read (99.8% of best
                        // hits) and starving the per-opportunity denominator. The real
                        // em_keep bit is carried in rec.em_keep. Walk qaln/taln here — the
                        // only place they are live — emitting every covered target position;
                        // p_damage / dist / damage_class come from summary.sites (distance-
                        // filtered for mismatched reads, empty for clean) so grain-B matches
                        // the --sites TSV. posterior is patched in after scoring below.
                        if (emit_dcx && ridx < em_read_results.size()) {
                            std::unordered_map<size_t, const DamageSite*> site_by_pos;
                            site_by_pos.reserve(summary.sites.size() * 2);
                            for (const auto& site : summary.sites)
                                site_by_pos.emplace(site.target_pos, &site);

                            DcxReadRec rec;
                            rec.read_idx = ridx;
                            rec.ref_idx = hit.ref_idx;
                            rec.src_idx = (static_cast<uint64_t>(rg_idx) << 32) | uint64_t(row_in_group);
                            rec.query_id = query_id;
                            rec.tstart = hit.tstart_0;
                            rec.strand = (query_id.find("_-_") != std::string::npos) ? 1 : 0;
                            const auto& em = em_read_results[ridx];
                            rec.gamma = em.gamma;
                            rec.fident = hit.fident;   // recruitment identity (COL_FIDENT)
                            rec.gamma_ancient = em.gamma_ancient;
                            rec.em_margin = em.em_margin;
                            rec.em_ref_neff = em.em_ref_neff;
                            rec.gamma_null = em.gamma_null;
                            rec.p_read = summary.p_read;   // reference-free per-read damage posterior
                            rec.em_keep = em.em_keep ? 1 : 0;
                            if (emit_multimap) {
                                // per-candidate: this (read,target)'s own gamma + rank; em_keep=1
                                // on exactly the argmax (rank 0) so v2.1 one-row-per-read consumers bridge.
                                rec.gamma = hit.gamma;
                                rec.gamma_ancient = hit.gamma_ancient;
                                rec.candidate_rank = hit.candidate_rank;
                                rec.em_keep = (hit.candidate_rank == 0) ? 1 : 0;
                            }
                            // AGD record carries this read's terminal codons (first/last 5),
                            // indexed by terminal AA distance -> the actual read codon at a
                            // terminal substitution (where damage occurs).
                            const AgdRecord* agd = damage_index ? damage_index->find(query_id) : nullptr;

                            std::vector<uint8_t> cover;
                            uint32_t tpos = hit.tstart_0;
                            const size_t naln = std::min(qaln_sv.size(), taln_sv.size());
                            for (size_t k = 0; k < naln; ++k) {
                                const char tc = taln_sv[k], qc = qaln_sv[k];
                                if (tc == '-') continue;   // query insertion: no target pos
                                const bool covered = (qc != '-');
                                cover.push_back(covered ? 1 : 0);
                                if (covered) {
                                    rec.ref_aa.push_back(uint8_t(tc));
                                    rec.obs_aa.push_back(uint8_t(qc));
                                    float pd = 0.0f; uint16_t d5 = 0, d3 = 0;
                                    uint8_t dc = uint8_t(dcx::DamageClass::None);
                                    uint8_t codon = 255;   // 255 = central / not in AGD terminal zone
                                    auto sit = site_by_pos.find(tpos);
                                    if (sit != site_by_pos.end()) {
                                        const DamageSite* s = sit->second;
                                        pd = s->p_damage;
                                        d5 = uint16_t(std::min<size_t>(s->dist_5prime, 65535));
                                        d3 = uint16_t(std::min<size_t>(s->dist_3prime, 65535));
                                        dc = (s->damage_class == 'C')
                                            ? uint8_t(dcx::DamageClass::CT)
                                            : uint8_t(dcx::DamageClass::GA);
                                        // actual read codon at this substitution: terminal codons
                                        // are indexed by terminal AA distance (one codon per AA).
                                        if (agd) {
                                            if (d5 < agd->n_5prime)      codon = agd->codons_5prime[d5];
                                            else if (d3 < agd->n_3prime) codon = agd->codons_3prime[d3];
                                        }
                                    }
                                    rec.p_damage.push_back(pd);
                                    rec.dist5.push_back(d5);
                                    rec.dist3.push_back(d3);
                                    rec.dcls.push_back(dc);
                                    rec.obs_codon.push_back(codon);
                                }
                                ++tpos;
                            }
                            rec.span = uint32_t(cover.size());
                            if (rec.span > 0) {
                                rec.presence.assign((rec.span + 7) / 8, 0);
                                for (uint32_t k = 0; k < rec.span; ++k)
                                    if (cover[k])
                                        rec.presence[k >> 3] |= uint8_t(1u << (k & 7));  // LSB-first
                                stage_dcx(rec, local_picks);
                            }
                        }

                        if (summary.total_mismatches > 0) {
                            // Return sites heap to allocator immediately (~10 GB for 40M reads).
                            summary.sites.clear();
                            summary.sites.shrink_to_fit();
                            // Return corrected-output strings (already written above).

                            local_results.emplace_back(ridx, std::move(summary));
                        }
                    }

                    #pragma omp critical
                    {
                        all_results.insert(all_results.end(),
                            std::make_move_iterator(local_results.begin()),
                            std::make_move_iterator(local_results.end()));
                        if (sites_out_f.is_open()) {
                            std::string s = local_sites.str();
                            if (!s.empty()) {
                                sites_out_f << s;
                                written_sites += static_cast<size_t>(
                                    std::count(s.begin(), s.end(), '\n'));
                            }
                        }
                        if (emit_dcx && !local_picks.empty()) {
                            dcx_pick_lists.push_back(std::move(local_picks));
                        }
                    }
                }
            });
        }

        if (verbose && sites_out_f.is_open()) {
            std::cerr << "Site calls written to: " << sites_file
                      << " (" << written_sites << " sites)\n";
        }

        if (dcx_stage_error) {
            try { std::rethrow_exception(dcx_stage_error); }
            catch (const std::exception& e) {
                std::cerr << "Error: --dcx run generation failed: " << e.what() << "\n";
                return 1;
            }
        }

        // Sort by read index for deterministic output order.
        std::sort(all_results.begin(), all_results.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
        summaries.reserve(all_results.size());
        for (auto& [ridx, summary] : all_results) {
            summary_index_by_read[ridx] = static_cast<int32_t>(summaries.size());
            summaries.push_back(std::move(summary));
        }
        all_results.clear();
        all_results.shrink_to_fit();
        // Return freed heap to OS before Bayesian scoring.
#ifdef __linux__
        malloc_trim(0);
#endif

        const auto annotate_end = std::chrono::steady_clock::now();
        if (verbose) {
            std::cerr << "Annotation runtime: "
                      << dart::log_utils::format_elapsed(annotate_start, annotate_end)
                      << " (" << summaries.size() << " reads with mismatches)\n";
        }
    }

    // Statistics
    size_t total_sites = 0, total_ct = 0, total_ga = 0, total_high = 0;
    size_t proteins_with_damage = 0;
    // Channel D: Oxidative damage statistics (G→T at DNA level)
    uint64_t total_gly_to_cys = 0, total_cys_to_gly = 0;
    uint64_t total_gly_to_trp = 0, total_trp_to_gly = 0;
    uint64_t total_gly = 0;
    for (const auto& s : summaries) {
        total_sites += s.damage_consistent;
        total_ct += s.ct_sites;
        total_ga += s.ga_sites;
        total_high += s.high_conf_sites;
        if (s.damage_consistent > 0) proteins_with_damage++;
        // Oxidative damage accumulation
        total_gly_to_cys += s.gly_to_cys;
        total_cys_to_gly += s.cys_to_gly;
        total_gly_to_trp += s.gly_to_trp;
        total_trp_to_gly += s.trp_to_gly;
        total_gly += s.gly_total;
    }

    // Auto-derive min_damage_sites from median read length if not set
    if (min_damage_sites < 0 && !summaries.empty()) {
        std::vector<size_t> qlens;
        qlens.reserve(summaries.size());
        for (const auto& s : summaries) qlens.push_back(s.qlen);
        std::nth_element(qlens.begin(), qlens.begin() + qlens.size() / 2, qlens.end());
        const size_t median_qlen = qlens[qlens.size() / 2];
        min_damage_sites = std::max(1, static_cast<int>(median_qlen / 10));
        if (verbose) std::cerr << "Auto min_damage_sites: " << min_damage_sites
                               << " (median qlen=" << median_qlen << " aa)\n";
    } else if (min_damage_sites < 0) {
        min_damage_sites = 3;  // fallback if no summaries
    }

    if (null_p_ran && !use_em_calib) {
        // Flagged reads whose best hit has no mismatch get no table row; they stay in the
        // modern-baseline pool but leave the authentication evidence. Count them apart.
        uint64_t flagged = 0, flagged_no_row = 0, with_hit = 0;
        for (uint32_t r = 0; r < n_reads; ++r) {
            if (!(em_read_results[r].null_p > null_p_max)) continue;
            ++flagged;
            if (summary_index_by_read[r] < 0) ++flagged_no_row;
        }
        for (int b = 0; b < dart::kNumNullBins; ++b) {
            with_hit += null_bin_reads[b];
            const double s = dart::null_chance_bits_at(em_null_model, null_p_scale, null_p_n, b, null_p_max);
            std::ostringstream line;
            line << "null reject: bin " << dart::kNullBinEdges[b] << "-";
            if (b + 1 < dart::kNumNullBins) line << (dart::kNullBinEdges[b + 1] - 1); else line << "inf";
            line << " reads=" << null_bin_reads[b] << " flagged=" << null_bin_flagged[b] << " ("
                 << std::fixed << std::setprecision(2)
                 << (null_bin_reads[b] ? 100.0 * null_bin_flagged[b] / null_bin_reads[b] : 0.0)
                 << "%) bits at null_p=" << std::defaultfloat << null_p_max << ": ";
            if (std::isfinite(s)) line << std::fixed << std::setprecision(1) << s; else line << "NA";
            std::cerr << line.str() << "\n";
        }
        std::ostringstream total;
        total << "null reject: " << flagged << " of " << with_hit << " reads with a best hit flagged ("
              << std::fixed << std::setprecision(2) << (with_hit ? 100.0 * flagged / with_hit : 0.0)
              << "%); " << flagged_no_row << " of them have no mismatch (no table row)";
        std::cerr << total.str() << "\n";
    }

    if (verbose) {
        std::cerr << "Reads with mismatches: " << summaries.size() << "\n";
        std::cerr << "Reads with damage-consistent sites: " << proteins_with_damage << "\n";
        std::cerr << "Total damage-consistent sites: " << total_sites << "\n";
        std::cerr << "  C->T class: " << total_ct << "\n";
        std::cerr << "  G->A class: " << total_ga << "\n";
        std::cerr << "  High confidence: " << total_high << "\n";
        // Channel D: Oxidative damage statistics
        const float gly_cys_ratio = total_cys_to_gly > 0
            ? static_cast<float>(total_gly_to_cys) / static_cast<float>(total_cys_to_gly)
            : 0.0f;
        const float gly_trp_ratio = total_trp_to_gly > 0
            ? static_cast<float>(total_gly_to_trp) / static_cast<float>(total_trp_to_gly)
            : 0.0f;
        std::cerr << "Channel D (oxidative damage, G→T at DNA):\n";
        std::cerr << "  Gly positions: " << total_gly << "\n";
        std::cerr << "  Gly→Cys: " << total_gly_to_cys << " | Cys→Gly: " << total_cys_to_gly
                  << " | Ratio: " << std::fixed << std::setprecision(2) << gly_cys_ratio << "x\n";
        std::cerr << "  Gly→Trp: " << total_gly_to_trp << " | Trp→Gly: " << total_trp_to_gly
                  << " | Ratio: " << std::fixed << std::setprecision(2) << gly_trp_ratio << "x\n";
        if (gly_cys_ratio > 1.15f || gly_trp_ratio > 1.15f) {
            std::cerr << "  [!] Elevated oxidative damage detected (ratio > 1.15x)\n";
        }
    }

    // Bayesian scoring setup
    // Estimate q_modern from likely-undamaged reads using AA opportunity accounting.
    std::vector<std::pair<uint32_t, uint32_t>> modern_proxy;
    double p_read_sum_all = 0.0;
    size_t p_read_n_all = 0;
    double p_read_sum_low = 0.0;
    size_t p_read_n_low = 0;
    float p_read_tail_cut = 0.20f;
    std::vector<float> p_read_values;
    p_read_values.reserve(summaries.size());

    for (const auto& s : summaries) {
        if (!s.has_p_read) continue;
        const float p = std::clamp(s.p_read, 0.0f, 1.0f);
        p_read_sum_all += p;
        p_read_n_all++;
        p_read_values.push_back(p);
    }

    float q20 = 0.20f;
    float q35 = 0.35f;
    if (!p_read_values.empty()) {
        std::vector<float> qbuf = p_read_values;
        q20 = quantile_inplace(qbuf, 0.20);
        qbuf = p_read_values;
        q35 = quantile_inplace(qbuf, 0.35);
    }
    p_read_tail_cut = std::clamp(q20, 0.02f, 0.35f);

    auto collect_modern_proxy = [&](float cutoff) {
        p_read_sum_low = 0.0;
        p_read_n_low = 0;
        modern_proxy.clear();
        for (const auto& s : summaries) {
            if (!s.has_p_read) continue;
            const float p = std::clamp(s.p_read, 0.0f, 1.0f);
            if (p > cutoff) continue;
            p_read_sum_low += p;
            p_read_n_low++;
            if (s.aa_m_opportunities > 0) {
                modern_proxy.emplace_back(s.aa_m_opportunities, s.aa_k_hits);
            }
        }
    };

    collect_modern_proxy(p_read_tail_cut);
    if (modern_proxy.empty() && p_read_n_all > 0) {
        const float widened_cut = std::clamp(q35, p_read_tail_cut, 0.50f);
        if (widened_cut > p_read_tail_cut + 1e-6f) {
            p_read_tail_cut = widened_cut;
            collect_modern_proxy(p_read_tail_cut);
        }
    }

    // BUG1 fix: the loop above only sees reads WITH a mismatch (they are the only ones
    // in `summaries`), so q_modern was estimated from a mismatch-only vector and biased
    // HIGH — the perfectly-matched reads (k=0, n>0), which are the cleanest technical
    // baseline, were absent. Pool every zero-mismatch best-hit read's opportunities as a
    // single (Σm_clean, 0) entry; estimate_modern_baseline sums (m,k), so aggregation is
    // exact and O(1) in memory. em-ambiguous reads are excluded to match the fit, except
    // for the uncalibrated null_p bit (0x1): a chance-level best hit says nothing about the
    // technical substitution rate, and dropping those reads would move q_modern and so every
    // other read's posterior.
    {
        const uint8_t pool_reject_mask = use_em_calib ? uint8_t(0xFF) : uint8_t(~0x1u);
        uint64_t clean_m_sum = 0;
        for (uint32_t r = 0; r < n_reads; ++r) {
            if (!best_hits[r].has) continue;
            if (r < summary_index_by_read.size() && summary_index_by_read[r] >= 0) continue;  // mismatched: already counted
            if (per_read_m_opp[r] == 0) continue;
            if (use_em) {
                if (r >= em_read_results.size() || !em_read_results[r].em_keep) continue;
                const auto& er = em_read_results[r];
                const bool ambiguous = ((er.em_reject_reason & pool_reject_mask) != 0) ||
                                       (er.em_margin < static_cast<float>(em_min_margin)) ||
                                       (er.em_ref_neff > static_cast<float>(em_max_neff));
                if (ambiguous) continue;
            }
            clean_m_sum += per_read_m_opp[r];
        }
        // Push the pooled clean denominator as one or more (m, 0) entries, each <= UINT32_MAX,
        // so estimate_modern_baseline's uint64_t M0 accumulates Σm_clean EXACTLY even when a
        // library's total clean opportunities exceed 2^32 (large libraries / long references).
        while (clean_m_sum > 0) {
            const uint32_t chunk = static_cast<uint32_t>(std::min<uint64_t>(clean_m_sum, UINT32_MAX));
            modern_proxy.emplace_back(chunk, static_cast<uint32_t>(0));
            clean_m_sum -= chunk;
        }
    }

    dart::ModernBaseline modern_est = dart::estimate_modern_baseline(modern_proxy);
    // Use fixed pi0 = 0.10 as baseline for terminal evidence transform
    // Empirical estimation fails for heavily damaged samples where all reads have high p_read
    // (no low-tail to establish baseline), causing terminal evidence to collapse
    constexpr float PI0_FIXED = 0.10f;
    const float empirical_pi0 = PI0_FIXED;

    // Self-calibrating protein-level NULL floor (π₀) for the Beta-Binomial damage
    // test (fit_protein_damage): the mean per-read p_damaged in the sample's LOW tail
    // -- the library's own non-ancient / technical baseline. This replaces the old
    // θ≈0 null that pinned p_protein_damaged at 1.0. Clamped to a technical floor
    // (~pi0_technical_floor.tsv) so a clean modern library still has a positive null,
    // capped at 0.5; falls back to the floor when the low tail is empty.
    // Self-calibrated protein null pi0 now lives in the libtaph engine
    // (taph::self_calibrated_pi0): DART supplies the low-tail sum/count (data).
    float sample_pi0_null = taph::self_calibrated_pi0(p_read_sum_low, p_read_n_low);
    if (verbose) std::cerr << "  Protein-null pi0 (self-calibrated): " << std::fixed
                           << std::setprecision(4) << sample_pi0_null << "\n";

    // Auto-calibrate prior_ancient from mean p_read if requested
    // E[p_read] = weighted average of P(damaged|read) ≈ fraction of ancient reads in the sample.
    // Requires AGD index (has_p_read); falls back to default if insufficient data.
    if (auto_prior_ancient) {
        if (p_read_n_all >= 50) {
            const float mean_p_read = static_cast<float>(p_read_sum_all / p_read_n_all);
            prior_ancient = std::clamp(mean_p_read, 0.05f, 0.80f);
            if (verbose) std::cerr << "Auto prior_ancient: " << std::fixed << std::setprecision(3)
                                   << prior_ancient << " (mean p_read=" << mean_p_read
                                   << " from " << p_read_n_all << " reads)\n";
        } else if (verbose) {
            std::cerr << "Auto prior_ancient: skipped (need >=50 reads with p_read, got "
                      << p_read_n_all << "); using default " << prior_ancient << "\n";
        }
    }

    // Set up scoring parameters with empirical Bayes approach
    dart::BayesianScoreParams score_params;
    score_params.pi = prior_ancient;
    score_params.pi0 = empirical_pi0;    // Empirical baseline for terminal evidence transform
    score_params.q_modern = modern_est.q_modern;
    score_params.w0 = w0;
    score_params.terminal_threshold = terminal_threshold;
    score_params.site_cap = site_cap;
    score_params.min_opportunities = min_damage_sites;
    score_params.channel_b_valid = (d_max > 0.0f);  // Trust site evidence if damage detected
    score_params.ancient_threshold = ancient_threshold;
    score_params.modern_threshold = modern_threshold;
    score_params.use_identity = true;
    score_params.k_identity = 20.0f;
    score_params.identity_baseline = 0.90f;
    score_params.d_max = d_max;  // For fixed-average site evidence formula

    // Sample-level env-damage flags (v4+) — captured before damage_index is released.
    bool  sample_ox_detected      = false;
    float sample_ox_rate_terminal = 0.0f;
    bool  sample_depurination     = false;
    float sample_purine_enrich_5p = 0.0f;
    float sample_purine_enrich_3p = 0.0f;

    // Damage informativeness gating from AGD v3 header
    // When uninformative: terminal evidence is zeroed (posterior still computed from sites+identity)
    float damage_detectability = 1.0f;
    // v5 library pre-gate: libtaph read_ancient_posterior only emits per-read damage when the
    // library's pi validated as DETECTED. A v5 index that is not DETECTED fails every protein
    // (surfaced, never silently AA-scored). Captured before the index is released.
    bool library_v5           = false;
    bool library_pi_detected  = false;
    bool library_pi_usable    = false;
    uint8_t library_pi_state  = 255;
    double  library_pi_point  = -1.0;
    // v6: library-level empirical-Bayes damage prior (pi_damaged, δ/d_max). Survives
    // LOW_ABUNDANCE, so it is the prior for per-protein p_protein_damaged — NOT the now-dead
    // EM pi. Captured before the index is released.
    double  library_damage_prior = -1.0;
    uint8_t library_damage_state = 255;
    uint8_t agd_library_type = 0;  // 0=unknown, 1=ss, 2=ds
    if (damage_index) {
        agd_library_type    = damage_index->library_type();
        library_v5          = damage_index->pi_state() != 255;
        library_pi_detected = damage_index->pi_detected();
        library_pi_usable   = damage_index->pi_usable();
        library_pi_state    = damage_index->pi_state();
        library_pi_point    = damage_index->pi_point();
        library_damage_prior = damage_index->damage_prior_point();
        library_damage_state = damage_index->damage_prior_state();
        score_params.damage_informative = damage_index->damage_informative();
        score_params.channel_b_valid = damage_index->channel_b_valid();
        damage_detectability = dart::compute_damage_detectability(
            damage_index->d_max(), damage_index->damage_validated(),
            damage_index->damage_artifact(), damage_index->channel_b_valid(),
            damage_index->stop_decay_llr(), damage_index->terminal_shift());
        score_params.damage_detectability = damage_detectability;
        // Env-damage channels (C/E) not yet implemented in AGD format; default to zero.
        sample_ox_detected      = false;
        sample_ox_rate_terminal = 0.0f;
        sample_depurination     = false;
        sample_purine_enrich_5p = 0.0f;
        sample_purine_enrich_3p = 0.0f;
        // All per-read lookups are done; release the 31 GB index before Bayesian scoring.
        damage_index.reset();
#ifdef __linux__
        malloc_trim(0);  // Return freed AGD pages to OS immediately
#endif
    }
    // State-aware gate. Only the abstain-equivalents with no usable pi point
    // (NOT_DETECTED/UNDETERMINED/BELOW_FLOOR, or a non-positive point) fail every
    // protein. The usable-advisory states (TRACE/ANCIENT_CPG/LOW_ABUNDANCE) carry a
    // real pi.point and are annotated: the per-read path already degrades to the
    // prior-free damage_evidence_llr ranking when read_ancient_posterior refuses, so a
    // library whose deamination is authenticated but whose independent-molecule count
    // lands LOW_ABUNDANCE is annotated (flagged advisory), never erased.
    if (library_v5 && !library_pi_usable) {
        std::cerr << "[libtaph gate] WARNING: library pi.state=" << static_cast<int>(library_pi_state)
                  << " has no usable pi — every protein fails the damage gate (library_not_authenticatable). "
                  << "libtaph emits no per-read damage claim for an abstain-state library.\n";
    } else if (library_v5 && library_pi_usable && !library_pi_detected) {
        std::cerr << "[libtaph gate] ADVISORY: library pi.state=" << static_cast<int>(library_pi_state)
                  << " is usable-advisory (not DETECTED(0)); proteins are annotated with a prior-free "
                  << "per-read ranking and flagged for stratification.\n";
    }

    if (verbose) {
        std::cerr << "\nBayesian scoring (empirical Bayes):\n";
        std::cerr << "  Modern proxy reads: " << modern_est.n_reads << "\n";
        std::cerr << "  Pooled opportunities: " << modern_est.M0 << "\n";
        std::cerr << "  Pooled hits: " << modern_est.K0 << "\n";
        std::cerr << "  p_read baseline (pi0): " << std::fixed << std::setprecision(4)
                  << score_params.pi0 << " (low-tail n=" << p_read_n_low
                  << ", all n=" << p_read_n_all
                  << ", tail_cut=" << std::fixed << std::setprecision(4) << p_read_tail_cut
                  << ")\n";
        std::cerr << "  Estimated q_modern: " << std::fixed << std::setprecision(5)
                  << modern_est.q_modern << "\n";
        std::cerr << "  Tempering w0: " << score_params.w0 << "\n";
        std::cerr << "  Channel B valid: " << (score_params.channel_b_valid ? "yes" : "no") << "\n";
        std::cerr << "  Damage informative: " << (score_params.damage_informative ? "yes" : "no") << "\n";
        std::cerr << "  Damage detectability: " << std::fixed << std::setprecision(3)
                  << damage_detectability << "\n";
        std::cerr << "  Identity evidence: enabled (k=20, baseline=0.90)\n";
    }

    // Phase 2: gamma-weighted MLE update of d_max (--refine-damage)
    if (refine_damage && score_params.d_max > 0.0f && use_em) {
        dart::ThetaDSuffStats theta_ss;

        // Numerator: weighted damage hits from reads that had any mismatches.
        // When damage-aware EM was active, use gamma_ancient directly (even when 0 —
        // that means a modern read and should not contribute to the ancient refit).
        // When damage-aware EM was not active, fall back to total gamma.
        for (const auto& s : summaries) {
            if (!s.has_p_read || s.aa_m_opportunities == 0) continue;
            const auto& er = em_read_results[s.read_idx];
            const double g = damage_em_active
                ? static_cast<double>(er.gamma_ancient)
                : static_cast<double>(er.gamma);
            if (g <= 0.0) continue;
            theta_ss.sum_gamma_k += g * static_cast<double>(s.aa_k_hits);
            ++theta_ss.n_reads;
        }

        // Denominator: weighted decay over ALL reads, including zero-mismatch
        // reads that were excluded from summaries.  Skipping these would bias
        // d_max_hat upward because their decay contributes to the expected
        // denominator without any matching hits in the numerator.
        for (uint32_t r = 0; r < n_reads; ++r) {
            if (per_read_aa_decay[r] <= 0.0f) continue;
            const auto& er = em_read_results[r];
            const double g = damage_em_active
                ? static_cast<double>(er.gamma_ancient)
                : static_cast<double>(er.gamma);
            if (g <= 0.0) continue;
            theta_ss.sum_gamma_decay += g * static_cast<double>(per_read_aa_decay[r]);
        }

        const float d_max_old = score_params.d_max;
        const float d_max_new = theta_ss.refit(score_params.d_max);
        score_params.d_max = d_max_new;
        if (verbose) {
            std::cerr << "  [refine-damage] d_max: " << d_max_old
                      << " -> " << d_max_new
                      << " (n=" << theta_ss.n_reads << " hit-reads, "
                      << n_reads << " total)\n";
        }
    }

    // Compute length-binned BPA z-scores (two-pass: accumulate then normalize)
    {
        dart::LengthBinStats bpa_stats;
        for (const auto& s : summaries) {
            if (s.has_p_read && s.p_read < 0.20f) {
                bpa_stats.add(s.qlen, s.bpa);
            }
        }
        bpa_stats.finalize();
        for (auto& s : summaries) {
            s.z_bpa = bpa_stats.z_score(s.qlen, s.bpa);
        }
    }

    // Write per-protein summary
    std::ostream* out = &std::cout;
    std::ofstream file_out;
    if (!output_file.empty()) {
        file_out.open(output_file);
        if (!file_out.is_open()) {
            std::cerr << "Error: Cannot open output file: " << output_file << "\n";
            return 1;
        }
        out = &file_out;
    }

    // Cache all ref_names via pread (NFS-safe) if not already done.
    // Streaming EM caches early (after flush_pages) when memory is lowest.
    // Non-streaming paths cache here, at output phase.
    if (cached_ref_names.empty()) {
        if (verbose) {
            std::cerr << "Caching ref_names via pread...\n";
        }
        cached_ref_names = reader.ref_names_pread();
        if (cached_ref_names.empty() && reader.num_refs() > 0) {
            std::cerr << "Error: ref_names_pread() returned empty vector\n";
            return 1;
        }
        if (verbose) {
            std::cerr << "Cached " << cached_ref_names.size() << " ref names\n";
        }
    }

    // Output header with Bayesian decomposition columns
    *out << "query_id\ttarget_id\tevalue\tbits\tfident\talnlen\t"
            "total_mismatches\tdamage_consistent\tnon_damage\t"
            "ct_sites\tga_sites\thigh_conf\tdamage_fraction\t"
            "max_p_damage\tp_damaged\tp_read\tposterior\tis_damaged\t"
            "classification\tlogBF_terminal\tlogBF_sites\tlogBF_identity\t"
            "m_opportunities\tk_hits\tq_eff\ttier\tdamage_class\t"
            "syn_5prime\tsyn_3prime\t"
            "delta_bits\tbpa\tlength_bin\tz_bpa\t"
            "damage_informative\t"
            "ox_hits\tox_total\tox_fraction";
    if (use_em) {
        *out << "\tgamma\tgamma_ancient\tbest_hit\tem_keep\tem_keep_threshold\tem_margin\tem_ref_neff\tem_ambiguous\tgamma_null\tnull_p\tem_reject_reason";
    }
    *out << "\n";

    // Compute and write per-read posteriors; store in summaries[si] for ReadRefEntry building.
    // query_id and target_id are reconstructed from read_idx via reader (not stored in struct).
    for (size_t si = 0; si < summaries.size(); ++si) {
        auto& s = summaries[si];
        // Compute Bayesian score with decomposition
        float posterior = s.p_damaged;  // fallback when AGD p_read is unavailable
        dart::BayesianScoreOutput bayes_out{};
        bayes_out.informative = score_params.damage_informative;

        if (s.has_p_read) {
            // AA-level evidence from alignment opportunities and damage-consistent hits.
            dart::SiteEvidence ev{};
            ev.k = s.aa_k_hits;
            ev.m = s.aa_m_opportunities;
            ev.sum_qA = s.aa_sum_qA;
            ev.sum_log_qA_hits = s.aa_sum_log_qA_hits;
            ev.q_eff = (ev.m > 0) ? (ev.sum_qA / static_cast<float>(ev.m)) : 0.0f;
            ev.sum_exp_decay = s.aa_sum_exp_decay;

            // Compute Bayesian score (with optional identity evidence)
            bayes_out = dart::compute_bayesian_score(s.p_read, ev, score_params, s.fident);
            posterior = bayes_out.posterior;
        }

        // Reconstruct names at output time. Query name via mmap (safe for read names).
        // Ref name from cached vector (NFS-safe pread).
        const std::string qname = reader.read_name(s.read_idx);
        const uint32_t ref_idx = best_hits[s.read_idx].ref_idx;
        const std::string& tname = cached_ref_names[ref_idx];

        *out << qname << "\t" << tname << "\t"
             << std::scientific << std::setprecision(2) << s.evalue << "\t"
             << std::fixed << std::setprecision(1) << s.bits << "\t"
             << std::fixed << std::setprecision(3) << s.fident << "\t"
             << s.alnlen << "\t"
             << s.total_mismatches << "\t" << s.damage_consistent << "\t"
             << s.non_damage_mismatches << "\t"
             << s.ct_sites << "\t" << s.ga_sites << "\t" << s.high_conf_sites << "\t"
             << std::fixed << std::setprecision(3) << s.damage_fraction << "\t"
             << std::fixed << std::setprecision(4) << s.max_p_damage << "\t"
             << std::fixed << std::setprecision(4) << s.p_damaged << "\t"
             << std::fixed << std::setprecision(4) << (s.has_p_read ? s.p_read : 0.0f) << "\t";
        // Classify using damage informativeness + posterior + quality metrics
        auto classification = dart::classify_protein(
            bayes_out.informative, posterior, s.fident, s.z_bpa, s.delta_bits,
            score_params.ancient_threshold, score_params.modern_threshold);

        s.posterior = posterior;
        s.classification = static_cast<uint8_t>(classification);

        // Always output posterior (computed from site + identity evidence even when
        // terminal is uninformative). is_damaged classification only applies when informative.
        *out << std::fixed << std::setprecision(4) << posterior << "\t";
        std::string is_damaged = ".";
        if (bayes_out.informative) {
            is_damaged = (posterior >= threshold) ? "1" : "0";
        }
        *out << is_damaged << "\t"
             << dart::classification_name(classification) << "\t";
        *out
             << std::fixed << std::setprecision(3) << bayes_out.logBF_terminal << "\t"
             << std::fixed << std::setprecision(3) << bayes_out.logBF_sites << "\t"
             << std::fixed << std::setprecision(3) << bayes_out.logBF_identity << "\t"
             << bayes_out.m << "\t" << bayes_out.k << "\t"
             << std::fixed << std::setprecision(5) << bayes_out.q_eff << "\t"
             << dart::tier_name(bayes_out.tier) << "\t"
             << dart::damage_class_name(bayes_out.damage_class) << "\t"
             << (s.has_syn_data ? std::to_string(s.syn_5prime) : ".") << "\t"
             << (s.has_syn_data ? std::to_string(s.syn_3prime) : ".") << "\t"
             << std::fixed << std::setprecision(1) << s.delta_bits << "\t"
             << std::fixed << std::setprecision(3) << s.bpa << "\t"
             << static_cast<int>(s.length_bin) << "\t"
             << std::fixed << std::setprecision(3) << s.z_bpa << "\t"
             << (bayes_out.informative ? 1 : 0)
             << "\t" << s.gly_to_cys
             << "\t" << s.gly_total
             << "\t" << std::fixed << std::setprecision(4)
             << (s.gly_total > 0 ? s.gly_to_cys / static_cast<float>(s.gly_total) : 0.0f);

        if (use_em) {
            if (s.read_idx < em_read_results.size()) {
                const auto& em = em_read_results[s.read_idx];
                *out << "\t" << std::scientific << std::setprecision(6) << em.gamma << "\t";
                if (damage_em_active) *out << std::scientific << std::setprecision(6) << em.gamma_ancient;
                *out << "\t" << (em.best_hit ? 1 : 0)
                     << "\t" << (em.em_keep ? 1 : 0)
                     << "\t" << std::fixed << std::setprecision(6) << em.em_threshold
                     << "\t" << std::fixed << std::setprecision(4) << em.em_margin
                     << "\t" << std::fixed << std::setprecision(4) << em.em_ref_neff
                     << "\t" << (((em.em_reject_reason != 0) ||
                                  (em.em_margin < static_cast<float>(em_min_margin)) ||
                                  (em.em_ref_neff > static_cast<float>(em_max_neff))) ? 1 : 0)
                     << "\t" << std::scientific << std::setprecision(6) << em.gamma_null;
                char np[32] = "NA";
                if (std::isfinite(em.null_p)) std::snprintf(np, sizeof np, "%.3g", em.null_p);
                *out << "\t" << np
                     << "\t" << static_cast<int>(em.em_reject_reason);
            } else {
                *out << "\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.";
            }
        }
        *out << "\n";
    }

    // ── .dcx interchange emit (--dcx) ────────────────────────────────────────
    // Merge target-sorted runs and patch the scored per-read posterior. Full and
    // grain-A-only writers cut their blocks independently, as before.
    if (emit_dcx) {
        std::vector<float> post_by_read(n_reads, 0.0f);
        // Deterministic per-read collapse. A multimapping read has one summary per candidate,
        // each with its own posterior; plain emplace() kept whichever appeared first in the
        // parallel-built summaries vector -> the kept value flipped between runs. Take the max
        // (order-independent) so the audit posterior is reproducible. posterior is audit-only
        // (convergence detect reads gamma/em_keep/p_damage, never posterior).
        for (const auto& s : summaries)
            post_by_read[s.read_idx] = std::max(post_by_read[s.read_idx], s.posterior);
        std::vector<std::string> a_targets;
        std::vector<uint32_t> a_target_refs;
        std::vector<uint64_t> a_N, a_k;
        std::vector<double> a_logLR;
        static const double AUTH_TAU[5] = {0.3, 0.4, 0.5, 0.6, 0.7};
        // Per read keep the em_keep record with the highest quantized gamma, ties to the
        // lower target rank; reduce in canonical read order so the floating-point sums
        // match the old in-memory reduction.
        std::vector<AuthPick> auth_picks;
        for (auto& list : dcx_pick_lists) auth_picks.insert(auth_picks.end(), list.begin(), list.end());
        std::vector<std::vector<AuthPick>>().swap(dcx_pick_lists);
        std::sort(auth_picks.begin(), auth_picks.end(), [](const AuthPick& a, const AuthPick& b) {
            if (a.read != b.read) return a.read < b.read;
            if (a.gamma != b.gamma) return a.gamma > b.gamma;
            return a.rank < b.rank;
        });
        std::unordered_map<uint32_t, uint32_t> auth_target_ord;
        for (size_t k = 0; k < auth_picks.size(); ++k) {
            const AuthPick& pick = auth_picks[k];
            if (k > 0 && auth_picks[k - 1].read == pick.read) continue;
            auto [it, inserted] = auth_target_ord.try_emplace(pick.ref, static_cast<uint32_t>(a_target_refs.size()));
            const uint32_t ix = it->second;
            if (inserted) {
                a_target_refs.push_back(pick.ref);
                a_N.push_back(0); a_logLR.push_back(0.0); a_k.insert(a_k.end(), 5, 0);
            }
            const double g = dcx::dec_unit(pick.gamma), pr = dcx::dec_unit(pick.prior);
            ++a_N[ix];
            const double pc = std::min(std::max(pr, 1e-6), 1.0 - 1e-6);
            a_logLR[ix] += g * std::log(pc / (1.0 - pc));
            for (int j = 0; j < 5; ++j) if (pr > AUTH_TAU[j]) ++a_k[ix * 5 + j];
        }
        std::vector<AuthPick>().swap(auth_picks);
        a_targets = ref_names_of(a_target_refs);
        uint32_t dcx_flags;
        dcx_flags = dcx::COL_EM_MARGIN |
                      dcx::COL_EM_REF_NEFF | dcx::COL_GAMMA_NULL |
                      dcx::COL_P_READ | dcx::COL_OBS_CODON |
                      dcx::COL_STRAND | dcx::COL_QID_DICT | dcx::COL_FIDENT;
        if (damage_em_active) dcx_flags |= dcx::COL_GAMMA_ANCIENT;
        if (emit_multimap)
            dcx_flags |= dcx::COL_SOFT_MULTIMAP | dcx::COL_CANDIDATE_RANK;
        dcx::DcxWriteOptions opt;
        opt.col_flags = dcx_flags;
        opt.sample_flags = 0;   // grain-B target_pos is 0-based (DART convention)
        // Test-only: force a small block cap so a modest input spans >1 .dcx block,
        // exercising the writer's block loop and the reader's cross-block walk.
        if (const char* cap = std::getenv("DART_DCX_BLOCK_CAP")) {
            unsigned long v = std::strtoul(cap, nullptr, 10);
            if (v > 0) opt.block_raw_cap = uint32_t(v);
        }
        opt.calib_model = "dart-positional/1";
        // library_label: explicit --library-label, else the .dcx output basename.
        if (!library_label.empty()) {
            opt.library_label = library_label;
        } else {
            std::string base = dcx_file;
            size_t slash = base.find_last_of('/');
            if (slash != std::string::npos) base = base.substr(slash + 1);
            size_t dot = base.rfind(".dcx");
            if (dot != std::string::npos) base = base.substr(0, dot);
            opt.library_label = base;
        }
        // Canonical contract preimage: the damage model identity + its final
        // params (positional hazard d_max*exp(-lambda*dist)). calib_contract_id =
        // two FNV-1a halves over this blob (deterministic, no new dependency).
        // pi_gate carries the library-level authentication verdict so convergence can
        // stratify DETECTED-sourced from advisory-sourced SAAVs from the .dcx alone
        // (pi_state is a library constant; no per-protein column needed). Part of the
        // contract preimage: a LOW_ABUNDANCE model scores reads differently than a
        // DETECTED one, so the calib id must reflect it.
        const char* pi_gate = !library_v5 ? "none"
                            : library_pi_detected ? "detected"
                            : library_pi_usable   ? "advisory"
                                                  : "fatal";
        std::ostringstream blob;
        blob << "model=dart-positional/1;d_max=" << std::fixed << std::setprecision(6)
             << d_max << ";lambda=" << std::setprecision(6) << lambda
             << ";pi_state=" << static_cast<int>(library_pi_state)
             << ";pi_point=" << std::setprecision(6) << library_pi_point
             << ";pi_gate=" << pi_gate;
        opt.calib_contract_blob = blob.str();
        const uint64_t h0 = dart::fnv1a_hash(opt.calib_contract_blob);
        const uint64_t h1 = dart::fnv1a_hash(opt.calib_contract_blob + "\x01");
        std::memcpy(opt.calib_contract_id, &h0, 8);
        std::memcpy(opt.calib_contract_id + 8, &h1, 8);

        const std::string auth_file = dcx_file + ".auth";
        std::vector<std::string> target_dict, query_dict;
        std::vector<uint32_t> target_refs;
        std::vector<uint64_t> query_ord(n_reads, UINT64_MAX);
        uint64_t dcx_records = 0, dcx_observations = 0;
        try {
            dcx::StreamWriter full(dcx_file.c_str(), opt);
            dcx::StreamWriter auth(auth_file.c_str(), opt, true);
            uint32_t last_rank = UINT32_MAX;
            dcx_sort->finish([&](const DcxReadRec& rec, uint32_t rank) {
                if (rank != last_rank) { target_refs.push_back(rec.ref_idx); last_rank = rank; }
                uint64_t& qord = query_ord[rec.read_idx];
                if (qord == UINT64_MAX) { qord = query_dict.size(); query_dict.push_back(rec.query_id); }
                const uint32_t tord = static_cast<uint32_t>(target_refs.size() - 1);
                const float posterior = post_by_read[rec.read_idx];
                full.append(rec, tord, qord, posterior);
                auth.append(rec, tord, qord, posterior);
                ++dcx_records; dcx_observations += rec.ref_aa.size();
            });
            target_dict = ref_names_of(target_refs);
            std::vector<uint32_t>().swap(target_refs);
            full.close(target_dict, query_dict);
            auth.close(target_dict, query_dict);
            dcx_sort.reset();
        } catch (const std::exception& e) {
            std::cerr << "Error: --dcx write failed: " << e.what() << "\n";
            return 1;
        }
        // EMIT the write-time auth-stats table (<dcx>.authstats) computed above.
        const std::string authstats_file = dcx_file + ".authstats";
        try {
            if (!dcx::write_auth_stats(authstats_file.c_str(), opt.library_label,
                                       AUTH_TAU, 5, a_targets, a_N, a_k, a_logLR, 19)) {
                std::cerr << "Error: authstats write returned false for " << authstats_file << "\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: authstats write failed: " << e.what() << "\n";
            return 1;
        }
        if (verbose) {
            std::cerr << ".dcx written to: " << dcx_file << " (" << dcx_records
                      << " alignments, " << dcx_observations << " covered positions, "
                      << target_dict.size() << " targets)\n";
            std::cerr << ".dcx auth sidecar: " << auth_file << " (grain-A only, "
                      << dcx_records << " alignments)\n";
            std::cerr << ".dcx authstats: " << authstats_file << " (write-time reduce, "
                      << a_targets.size() << " targets)\n";
        }
    }

    // Extract compact protein_agg inputs from summaries, then free the full summaries
    // vector before building gene_agg. Uses pre-aggregated site scalars (info_sum,
    // terminal_5, terminal_3) populated during annotation — no per-read sites heap.
    const bool need_protein_agg =
        !protein_summary_file.empty();
    std::vector<ProteinAggInput> protein_agg_inputs;
    if (need_protein_agg) {
        protein_agg_inputs.reserve(summaries.size());
        for (const auto& s : summaries) {
            ProteinAggInput c;
            c.read_idx = s.read_idx;
            c.ref_idx = best_hits[s.read_idx].ref_idx;
            c.tlen = static_cast<uint16_t>(std::min<size_t>(s.tlen, UINT16_MAX));
            c.p_damaged = s.p_damaged;
            c.p_read = s.p_read;
            c.has_p_read = s.has_p_read;
            c.damage_consistent = static_cast<uint32_t>(s.damage_consistent);
            c.info_sum = s.info_sum;
            c.terminal_5 = s.terminal_5;
            c.terminal_3 = s.terminal_3;
            c.posterior = s.posterior;
            c.gly_to_cys = s.gly_to_cys;
            c.gly_total = s.gly_total;
            c.aa_m_opportunities = s.aa_m_opportunities;
            c.aa_k_hits = s.aa_k_hits;
            // EM confidence gate: a read is ambiguous when it is null-dominated
            // (gamma_null > tau), has a weak top1-top2 margin, or too many effective
            // references. Margin/neff are computed on real-only responsibilities (see
            // em_read_results fill), so the null does not double-penalize. Ambiguous
            // reads are excluded from the protein's damage/authentication evidence
            // (below) but still contribute to abundance/coverage via read_ref_entries.
            if (use_em && s.read_idx < em_read_results.size()) {
                const auto& em = em_read_results[s.read_idx];
                c.gamma = em.gamma;   // real mapping responsibility for auth weighting (2b)
                c.em_ambiguous = (em.em_reject_reason != 0) ||
                                 (em.em_margin < static_cast<float>(em_min_margin)) ||
                                 (em.em_ref_neff > static_cast<float>(em_max_neff));
            }
            protein_agg_inputs.push_back(c);
        }
    }
    const size_t n_summaries = summaries.size();

    // =====================================================================
    // Streaming output: sort reads by ref_idx, process one ref at a time.
    // Avoids materializing 9.9M-entry gene_agg / protein_agg maps (~35 GB).
    // =====================================================================

    const bool need_any_aggregation =
        !gene_summary_file.empty() || !functional_summary_file.empty() ||
        !anvio_ko_file.empty() || !protein_summary_file.empty();

    // 1. Build ReadRefEntry list (only em_keep reads), sort by ref_idx.
    // summaries must remain alive here: e.posterior/e.classification are read from
    // summaries[si] (populated during the TSV output loop above). In the old code
    // these were in separate read_posteriors/read_classifications arrays; now they
    // live directly on the compact summaries structs.
    std::vector<ReadRefEntry> read_ref_entries;
    size_t assigned_kept_reads = 0;

    if (need_any_aggregation) {
        read_ref_entries.reserve(std::min<size_t>(n_reads, size_t{1} << 26));
        for (uint32_t ridx = 0; ridx < n_reads; ++ridx) {
            const auto& hit = best_hits[ridx];
            if (!hit.has) continue;
            if (use_em && (ridx >= em_read_results.size() || !em_read_results[ridx].em_keep))
                continue;
            ReadRefEntry e;
            e.ref_idx              = hit.ref_idx;
            e.qstart_0             = static_cast<uint32_t>(hit.qstart_0);
            e.tstart_0             = static_cast<uint32_t>(hit.tstart_0);
            e.gamma                = use_em ? em_read_results[ridx].gamma : 1.0f;
            e.fident               = hit.fident;
            e.bits                 = hit.bits;
            e.bits_second          = hit.bits_second;
            e.aln_len_on_target    = static_cast<uint32_t>(hit.aln_len_on_target);
            e.tlen                 = static_cast<uint16_t>(std::min<uint32_t>(static_cast<uint32_t>(hit.tlen), UINT16_MAX));
            e.length_bin           = static_cast<uint8_t>(dart::LengthBinStats::get_bin(hit.qlen));
            const int32_t si       = (ridx < summary_index_by_read.size())
                                     ? summary_index_by_read[ridx] : -1;
            e.has_damage_obs       = (si >= 0);
            e.auth_keep            = true;
            if (si >= 0 && static_cast<size_t>(si) < summaries.size()) {
                e.posterior        = summaries[static_cast<size_t>(si)].posterior;
                e.classification   = summaries[static_cast<size_t>(si)].classification;
            } else {
                e.posterior        = 0.0f;
                e.classification   = 0;
            }
            // BUG1: opportunity count + EM-ambiguity for the per-opportunity fit. For clean
            // reads (si<0, absent from protein_agg_inputs) this is the ONLY carrier of their
            // 0/n_eligible negative evidence into the protein-level Binomial below.
            e.aa_m_opportunities = (ridx < per_read_m_opp.size()) ? per_read_m_opp[ridx] : 0;
            {
                const auto& er = em_read_results[ridx];
                e.em_ambiguous = (use_em && ((er.em_reject_reason != 0) ||
                                             (er.em_margin < static_cast<float>(em_min_margin)) ||
                                             (er.em_ref_neff > static_cast<float>(em_max_neff)))) ? 1 : 0;
            }
            read_ref_entries.push_back(e);
        }
        std::sort(read_ref_entries.begin(), read_ref_entries.end(),
            [](const ReadRefEntry& a, const ReadRefEntry& b) { return a.ref_idx < b.ref_idx; });
        assigned_kept_reads = read_ref_entries.size();
    } else {
        for (uint32_t ridx = 0; ridx < n_reads; ++ridx) {
            if (!best_hits[ridx].has) continue;
            if (use_em && (ridx >= em_read_results.size() || !em_read_results[ridx].em_keep)) continue;
            ++assigned_kept_reads;
        }
    }

    // Free summaries now that read_ref_entries has extracted posterior/classification.
    { decltype(summaries) tmp; tmp.swap(summaries); }
#ifdef __linux__
    malloc_trim(0);  // Return freed summaries heap to OS before gene_agg/protein_agg
#endif

    // 2. Filter protein_agg_inputs by em_keep, then sort by ref_idx.
    if (!protein_agg_inputs.empty() && use_em) {
        auto it = std::remove_if(protein_agg_inputs.begin(), protein_agg_inputs.end(),
            [&](const ProteinAggInput& c) {
                return c.read_idx >= em_read_results.size() ||
                       !em_read_results[c.read_idx].em_keep;
            });
        protein_agg_inputs.erase(it, protein_agg_inputs.end());
    }
    if (!protein_agg_inputs.empty()) {
        std::sort(protein_agg_inputs.begin(), protein_agg_inputs.end(),
            [](const ProteinAggInput& a, const ProteinAggInput& b) {
                return a.ref_idx < b.ref_idx;
            });
    }

    // 3. Free large per-read vectors — no longer needed.
    { decltype(best_hits) tmp; tmp.swap(best_hits); }
    { decltype(em_read_results) tmp; tmp.swap(em_read_results); }
    { decltype(summary_index_by_read) tmp; tmp.swap(summary_index_by_read); }
#ifdef __linux__
    malloc_trim(0);
#endif

    if (!need_any_aggregation) {
        if (verbose) {
            std::cerr << "\nSummary: assigned_reads=" << assigned_kept_reads
                      << ", mismatch_reads=" << n_summaries
                      << ", reads_with_damage_sites=" << proteins_with_damage << "\n";
            const auto run_end = std::chrono::steady_clock::now();
            std::cerr << "Runtime: "
                      << dart::log_utils::format_elapsed(run_start, run_end) << "\n";
        }
        return 0;
    }

    // 4. Auto-calibration pre-pass (if --auto-calibrate-spurious).
    if (auto_calibrate_spurious) {
        std::vector<float> calib_pos_scores;
        std::vector<float> calib_term_ratios;

        auto ri = read_ref_entries.cbegin();
        while (ri != read_ref_entries.cend()) {
            const uint32_t cur_ref = ri->ref_idx;
            auto ri_end = ri;
            while (ri_end != read_ref_entries.cend() && ri_end->ref_idx == cur_ref) ++ri_end;

            GeneSummaryAccumulator calib_acc;
            calib_acc.tlen = ri->tlen;
            if (calib_acc.tlen > 0) calib_acc.coverage.assign(calib_acc.tlen, 0.0f);
            for (auto it = ri; it != ri_end; ++it)
                calib_acc.add_assignment(it->fident, it->tstart_0, it->aln_len_on_target, it->gamma);

            if (calib_acc.n_reads >= 10 && calib_acc.n_unique_starts() >= 3) {
                calib_pos_scores.push_back(calib_acc.positional_score());
                calib_term_ratios.push_back(calib_acc.terminal_middle_ratio());
            }
            ri = ri_end;
        }

        if (calib_pos_scores.size() >= 100) {
            std::sort(calib_pos_scores.begin(), calib_pos_scores.end());
            std::sort(calib_term_ratios.begin(), calib_term_ratios.end());
            const size_t idx_5pct = calib_pos_scores.size() / 20;
            float calib_min_pos  = calib_pos_scores[idx_5pct];
            float calib_min_term = std::max(0.20f, calib_term_ratios[idx_5pct]);
            if (calib_min_pos > min_positional_score)
                min_positional_score = calib_min_pos;
            if (!user_set_terminal_ratio && calib_min_term > min_terminal_ratio)
                min_terminal_ratio = calib_min_term;
            if (verbose) {
                std::cerr << "Auto-calibration (n=" << calib_pos_scores.size() << " genes):\n"
                          << "  min_positional_score: " << std::fixed << std::setprecision(4)
                          << min_positional_score << "\n"
                          << "  min_terminal_ratio: " << std::fixed << std::setprecision(4)
                          << min_terminal_ratio << "\n";
            }
        } else if (verbose) {
            std::cerr << "Auto-calibration: insufficient calibration set ("
                      << calib_pos_scores.size() << " genes, need >= 100)\n";
        }
    }

    // 5. Library type detection (needs terminal site counts).
    bool use_both_ends = true;
    if (lib_type == "ss") {
        use_both_ends = false;
    } else if (lib_type == "ds") {
        use_both_ends = true;
    } else if (agd_library_type == 1 || agd_library_type == 2) {
        use_both_ends = (agd_library_type == 2);
        if (verbose) {
            std::cerr << "Library type from damage index: " << (use_both_ends ? "dsDNA" : "ssDNA") << "\n";
        }
    } else {
        size_t total_5 = 0, total_3 = 0;
        for (const auto& c : protein_agg_inputs) {
            total_5 += c.terminal_5;
            total_3 += c.terminal_3;
        }
        float ratio = (total_5 > 0) ? static_cast<float>(total_3) / static_cast<float>(total_5) : 0.0f;
        use_both_ends = (ratio > 0.3f);
        if (verbose) {
            std::cerr << "Library type auto-detected: " << (use_both_ends ? "dsDNA" : "ssDNA")
                      << " (3'/5' ratio: " << ratio << ")\n";
        }
    }
    const float n_ends_lib = use_both_ends ? 2.0f : 1.0f;

    // 6. Load mapping file.
    std::unordered_map<std::string, std::string> group_map;
    if (!map_file.empty()) {
        group_map = load_mapping_file(map_file, verbose);
        if (verbose && !group_map.empty()) {
            const auto& first = *group_map.cbegin();
            std::cerr << "  Mapping sample — first key: \"" << first.first
                      << "\" -> \"" << first.second << "\"\n";
        }
    }

    // 7. Open output files.
    std::ofstream gf, pf_os, af_f;
    if (!gene_summary_file.empty()) {
        gf.open(gene_summary_file);
        if (!gf.is_open()) {
            std::cerr << "Error: Cannot open gene summary file: " << gene_summary_file << "\n";
            return 1;
        }
        gf << "target_id\tn_reads\teff_reads\tn_damage_reads\teff_damage_reads\t"
              "n_ancient\tn_modern\tn_undetermined\t"
              "ancient_frac\tci_low\tci_high\tmean_posterior\tmean_identity\t"
              "breadth\tabundance\tdepth_std\tdepth_evenness\t"
              "n_unique_starts\tstart_diversity\tstart_span\tpositional_score\t"
              "terminal_middle_ratio\tavg_aln_len\tstd_aln_len\t"
              "coverage_deviance\tcompleteness_fhat\tcoverage_weight\tabundance_adjusted\n";
    }
    if (!anvio_ko_file.empty()) {
        af_f.open(anvio_ko_file);
        if (!af_f.is_open()) {
            std::cerr << "Error: Cannot open Anvi'o KO file: " << anvio_ko_file << "\n";
            return 1;
        }
        af_f << "gene_callers_id\tenzyme_accession\tsource\tcoverage\tdetection\n";
    }
    constexpr float P_SUSCEPTIBLE_K = 0.25f;
    constexpr float P_NONSYNONYMOUS_K = 0.70f;
    const std::string protein_header =
        "protein_id\tn_reads\teff_reads\tn_damage_reads\tn_damaged\tdamage_sites\t"
        "terminal_5\tterminal_3\t"
        "max_p_damaged\tmean_p_damaged\tfrac_damaged\tfrac_damaged_all\t"
        "protein_damage_confidence\tdamage_confidence_ungated\tdelta_mle\tci_lower\tci_upper\t"
        "lrt_pvalue\tlog_bf\td_aa\t"
        "avg_delta_bits\tavg_bpa\tlength_bin\t"
        "n_unique_starts\tstart_diversity\tstart_span\tpositional_score\t"
        "terminal_middle_ratio\tavg_aln_len\tstd_aln_len\t"
        "assign_n_reads\tassign_eff_reads\tassign_n_damage_reads\tassign_eff_damage_reads\t"
        "assign_n_ancient\tassign_n_modern\tassign_n_undetermined\t"
        "assign_ancient_frac\tassign_ci_low\tassign_ci_high\t"
        "assign_mean_posterior\tassign_mean_identity\t"
        "assign_breadth\tassign_abundance\tassign_depth_std\tassign_depth_evenness\t"
        "coverage_deviance\tcompleteness_fhat\tcoverage_weight\tabundance_adjusted\t"
        "ox_hits\tox_total\tp_ox_protein\tany_damage_pass\t"
        "ox_detected\tox_rate_terminal\t"
        "depurination_detected\tpurine_enrich_5p\tpurine_enrich_3p\t"
        "pass_mapping_filter\tfilter_fail\tmax_ni_frac\tn_informative_reads\t"
        "p_protein_damaged\tp_protein_damaged_band_low\tp_protein_damaged_band_high\t"
        "protein_damage_prior\tlibrary_damage_state\tprotein_auth_tier\t"
        "auth_binom_k\tauth_binom_n\tauth_binom_logLR\tauth_binom_p_ancient\t"
        "auth_binom_pvalue\tauth_binom_call\n";
    if (!protein_summary_file.empty()) {
        pf_os.open(protein_summary_file);
        if (!pf_os.is_open()) {
            std::cerr << "Error: Cannot open protein summary file: " << protein_summary_file << "\n";
            return 1;
        }
        pf_os << protein_header;
    }

    // 8. Functional accumulator (written after streaming pass).
    std::unordered_map<std::string, FunctionalAccumulator> func_acc;

    // 9. Main streaming pass: one ref at a time.
    constexpr size_t TERMINAL_CODONS_S = 3;
    size_t gene_written = 0, protein_written = 0;

    auto ri_cur = read_ref_entries.cbegin();
    // Binomial provenance-LR library constants (uniform across the library → computed once).
    // pi = library damaged-read fraction (LSD-joint pi_damaged = damage_prior_point); eps0 = modern
    // baseline; prov_min_n = min independent-read count for ~80% power at auth_alpha to discriminate
    // pi from eps0. Low-pi libraries yield a huge prov_min_n → most proteins land undetermined.
    const double prov_eps0 = std::clamp(static_cast<double>(modern_damage_rate), 1e-4, 0.499);
    const bool   prov_pi_available = (library_damage_prior > 0.0) && library_pi_usable;
    const double prov_pi = prov_pi_available
        ? std::clamp(library_damage_prior, prov_eps0 + 1e-6, 0.98) : -1.0;
    constexpr int PROV_NMAX = 100000;
    const int prov_min_n = prov_pi_available
        ? min_n_for_power(prov_pi, prov_eps0, static_cast<double>(auth_alpha), 0.80, PROV_NMAX)
        : PROV_NMAX;
    if (verbose && prov_pi_available) {
        std::cerr << "  Provenance-LR: pi=" << prov_pi << " eps0=" << prov_eps0
                  << " alpha=" << auth_alpha << " -> min_N(80% power)=" << prov_min_n << "\n";
    }

    auto pi_cur = protein_agg_inputs.cbegin();

    while (ri_cur != read_ref_entries.cend() || pi_cur != protein_agg_inputs.cend()) {
        uint32_t cur_ref;
        if (ri_cur != read_ref_entries.cend() && pi_cur != protein_agg_inputs.cend())
            cur_ref = std::min(ri_cur->ref_idx, pi_cur->ref_idx);
        else if (ri_cur != read_ref_entries.cend())
            cur_ref = ri_cur->ref_idx;
        else
            cur_ref = pi_cur->ref_idx;

        auto ri_end = ri_cur;
        while (ri_end != read_ref_entries.cend() && ri_end->ref_idx == cur_ref) ++ri_end;
        auto pi_end = pi_cur;
        while (pi_end != protein_agg_inputs.cend() && pi_end->ref_idx == cur_ref) ++pi_end;

        // Use cached ref_names vector (NFS-safe pread).
        const std::string& target_name = cached_ref_names[cur_ref];

        // --- Gene accumulator ---
        GeneSummaryAccumulator gene_acc;
        gene_acc.target_id = target_name;
        if (ri_cur != ri_end) {
            gene_acc.tlen = ri_cur->tlen;
            if (gene_acc.tlen > 0) gene_acc.coverage.assign(gene_acc.tlen, 0.0f);
            for (auto it = ri_cur; it != ri_end; ++it) {
                gene_acc.add_assignment(it->fident, it->tstart_0, it->aln_len_on_target, it->gamma);
                if (it->has_damage_obs) {
                    gene_acc.add_damage_observation(
                        static_cast<dart::AncientClassification>(it->classification),
                        it->posterior, it->gamma);
                }
            }
        }

        // --- Protein aggregate ---
        ProteinAggregate pa;
        pa.protein_id = target_name;
        for (auto it = ri_cur; it != ri_end; ++it) {
            if (pa.tlen == 0 && it->tlen > 0) pa.tlen = it->tlen;
            pa.n_reads++;
            pa.n_reads_eff += it->gamma;
            pa.sum_delta_bits += it->gamma * (it->bits - it->bits_second);
            const float bpa_v = (it->aln_len_on_target > 0)
                ? it->bits / static_cast<float>(it->aln_len_on_target) : 0.0f;
            pa.sum_bpa += it->gamma * bpa_v;
            pa.length_bin_counts[it->length_bin]++;
            pa.unique_starts.insert(it->tstart_0);
            pa.unique_start_pairs.insert((static_cast<uint64_t>(it->qstart_0) << 32) |
                                         static_cast<uint64_t>(it->tstart_0));
            if (it->tstart_0 < pa.min_start) pa.min_start = it->tstart_0;
            if (it->tstart_0 > pa.max_start) pa.max_start = it->tstart_0;
            // Per-position gamma-weighted coverage (mirrors GeneSummaryAccumulator).
            if (pa.tlen > 0 && pa.coverage.empty()) pa.coverage.assign(pa.tlen, 0.0f);
            if (!pa.coverage.empty()) {
                const size_t cs = std::min(static_cast<size_t>(it->tstart_0), pa.coverage.size());
                const size_t ce = std::min(cs + static_cast<size_t>(it->aln_len_on_target),
                                           pa.coverage.size());
                for (size_t p = cs; p < ce; ++p) pa.coverage[p] += it->gamma;
            }
            const size_t aln = it->aln_len_on_target;
            pa.sum_aln_len += it->gamma * static_cast<double>(aln);
            pa.sum_aln_len_sq += it->gamma * static_cast<double>(aln) * static_cast<double>(aln);
            const size_t tlen_sz = std::max<size_t>(1, static_cast<size_t>(it->tlen));
            const size_t window = std::min(std::max<size_t>(1, tlen_sz / 10), tlen_sz / 2);
            const size_t s_pos = std::min(static_cast<size_t>(it->tstart_0), tlen_sz);
            const size_t e_pos = std::min(s_pos + aln, tlen_sz);
            if (e_pos > s_pos) {
                const size_t left_end    = window;
                const size_t right_start = (tlen_sz > window) ? (tlen_sz - window) : 0;
                const size_t lo = (s_pos < left_end) ? (std::min(e_pos, left_end) - s_pos) : 0;
                const size_t ro = (e_pos > right_start && s_pos < tlen_sz)
                    ? (e_pos - std::max(s_pos, right_start)) : 0;
                size_t mo = 0;
                if (right_start > left_end) {
                    const size_t ms = std::max(s_pos, left_end);
                    const size_t me = std::min(e_pos, right_start);
                    if (me > ms) mo = me - ms;
                }
                pa.left_cov   += it->gamma * static_cast<double>(lo);
                pa.right_cov  += it->gamma * static_cast<double>(ro);
                pa.middle_cov += it->gamma * static_cast<double>(mo);
            }
        }
        uint32_t total_ox_hits = 0;
        uint32_t total_ox_total = 0;
        for (auto it = pi_cur; it != pi_end; ++it) {
            if (pa.tlen == 0 && it->tlen > 0) pa.tlen = it->tlen;
            // EM-ambiguous reads are excluded from damage/authentication evidence.
            if (it->em_ambiguous) continue;
            pa.n_damage_reads++;
            // Authenticate on the AGD REFERENCE-FREE per-read posterior it->p_read (the .agd value
            // from predict's prior-free terminal-site LLR), NOT the reference-BASED it->p_damaged
            // (annotate's per-site LLR sum, which folds in the reference-identity term #5 and is
            // capped near ~0.11 for ancient reads). p_read keeps p_protein_damaged free of the
            // reference-identity term (#5) and short-read capable (#1) — the reason libtaph was
            // wired in. A read WITHOUT an AGD p_read (has_p_read==false) carries no reference-free
            // damage evidence: it is EXCLUDED from authentication entirely (no fallback injected),
            // though it still counts for abundance/coverage via n_damage_reads and for the
            // per-opportunity fit via read_obs.
            pa.total_damage_sites += it->damage_consistent;
            const float p_eff = it->p_read;
            if (it->has_p_read) {
                pa.n_auth_reads++;
                // Per-read "damaged" call for the binomial-LR k: the reference-FREE MAP call on
                // p_read (identity EXCLUDED). For a DETECTED library p_read=sigmoid(llr+logit(pi)),
                // so p_read>0.5 <=> llr>-logit(pi): the threshold adapts to the library prior, which
                // is what keeps library-wide k/N tied to pi (per-read-damage-threshold default 0.5).
                if (p_eff > per_read_damage_threshold) pa.n_damaged++;
                if (p_eff > pa.max_p_damaged) pa.max_p_damaged = p_eff;
                pa.sum_p_damaged += p_eff;
                // Stage 2b: bank the per-read PRIOR-FREE damage LLR weighted by the REAL EM
                // responsibility. p_read is the .agd reference-free terminal-site posterior:
                //   read_ancient_posterior fires IFF pi.state==DETECTED, storing
                //     p_read = sigmoid(llr + logit(pi.point)) ⇒ logit(p_read) = llr + logit(pi.point);
                //   otherwise (LOW_ABUNDANCE etc.) the fallback stores
                //     p_read = sigmoid(llr) ⇒ logit(p_read) = llr (already prior-free).
                // Recover the raw prior-free llr by subtracting the baked-in logit(pi.point) ONLY
                // for DETECTED libraries; else subtract 0. Without this a DETECTED library
                // re-injects the logit(pi.point)·Σγ coverage bias the Stage-2b fix removed. (Only
                // 8-bit p_read quantization is lost — invariant under the rank/sign validation.)
                const double g = static_cast<double>(it->gamma);
                const double p = std::clamp(static_cast<double>(p_eff), 1e-4, 1.0 - 1e-4);
                double llr = std::log(p / (1.0 - p));
                if (library_pi_detected && library_pi_point > 0.0) {
                    const double pp = std::clamp(library_pi_point, 1e-4, 1.0 - 1e-4);
                    llr -= std::log(pp / (1.0 - pp));   // strip the DETECTED posterior's baked prior
                }
                pa.auth_sum_g_llr += g * llr;
                pa.auth_sum_g2    += g * g;
            }
            ReadDamageObs obs;
            obs.p_damaged = p_eff;
            obs.log_lr    = 0.0f;
            obs.info      = it->info_sum;
            obs.n_sites   = static_cast<int>(it->damage_consistent);
            obs.is_damaged = p_eff > 0.5f;
            // #2: per-opportunity counts for fit_protein_damage_binomial (k damaged / n eligible).
            // Ceiling: aa_m_opportunities/aa_k_hits are PER-READ counts bounded by the read's AA
            // length (qlen/3, tens of codons) — orders of magnitude below INT_MAX, so this cast to
            // int cannot overflow. The pooled per-protein N=Σn_eligible, K=Σk_damaged accumulate in
            // double inside fit_protein_damage_binomial (protein_damage.hpp) and in uint64 in
            // estimate_modern_baseline — no pooled overflow either. The whole path is width-safe.
            obs.n_eligible = static_cast<int>(it->aa_m_opportunities);
            obs.k_damaged  = static_cast<int>(it->aa_k_hits);
            pa.read_obs.push_back(obs);
            pa.terminal_5 += it->terminal_5;
            pa.terminal_3 += it->terminal_3;
            total_ox_hits  += it->gly_to_cys;
            total_ox_total += it->gly_total;
        }

        // BUG1: eligible-clean reads (zero mismatch → no summary → absent from
        // protein_agg_inputs above) still carry damage OPPORTUNITIES. Feed their
        // 0/n_eligible negative evidence into the per-opportunity fit from
        // read_ref_entries, which covers ALL best-hit reads. Deliberately NOT added to
        // n_damage_reads / sum_p_damaged (legacy soft-posterior columns) — only to the
        // opportunity likelihood, so p_protein_damaged sees the true denominator without
        // perturbing the legacy per-read damage stats.
        for (auto it = ri_cur; it != ri_end; ++it) {
            if (it->has_damage_obs) continue;      // mismatched: already in pa.read_obs
            if (it->em_ambiguous) continue;        // EM-ambiguous excluded from evidence
            if (it->aa_m_opportunities == 0) continue;  // info-free
            ReadDamageObs obs;
            obs.p_damaged  = 0.0f;
            obs.log_lr     = 0.0f;
            obs.info       = 0.0f;
            obs.n_sites    = 0;
            obs.is_damaged = false;
            obs.n_eligible = static_cast<int>(it->aa_m_opportunities);
            obs.k_damaged  = 0;
            pa.read_obs.push_back(obs);
        }

        // --- Coverage EM stats for this ref ---
        float cov_deviance = 0.0f, cov_fhat = 1.0f, cov_weight = 1.0f;
        {
            auto cit = coverage_stats_map.find(cur_ref);
            if (cit != coverage_stats_map.end()) {
                cov_deviance = cit->second.coverage_deviance;
                cov_fhat     = cit->second.completeness_fhat;
                cov_weight   = cit->second.coverage_weight;
            }
        }

        // --- Gene filter ---
        const float gene_breadth   = gene_acc.breadth();
        const float gene_abundance = gene_acc.depth_mean();
        const float gene_pos_score = gene_acc.positional_score();
        const float gene_term_ratio= gene_acc.terminal_middle_ratio();
        const bool gene_passes =
            gene_acc.n_reads >= min_reads &&
            gene_breadth  >= min_breadth &&
            gene_abundance >= min_depth &&
            gene_pos_score >= min_positional_score &&
            gene_term_ratio >= min_terminal_ratio;

        // --- Gene output ---
        if (gf.is_open() && gene_passes) {
            const float n_ancient = static_cast<float>(gene_acc.n_ancient_conf + gene_acc.n_ancient_likely);
            const float n_modern  = static_cast<float>(gene_acc.n_modern_conf);
            const float n_undet   = static_cast<float>(gene_acc.n_undetermined);
            const double n_eff_class = static_cast<double>(n_ancient + n_modern);
            auto [ci_lo, ci_hi] = wilson_ci_effective(n_ancient, n_eff_class);
            const float assign_eff = static_cast<float>(gene_acc.eff_reads());
            gf << target_name
               << '\t' << gene_acc.n_reads
               << '\t' << std::fixed << std::setprecision(2) << assign_eff
               << '\t' << gene_acc.n_damage_reads
               << '\t' << std::setprecision(2) << gene_acc.n_damage_reads_eff
               << '\t' << std::setprecision(2) << n_ancient
               << '\t' << std::setprecision(2) << n_modern
               << '\t' << std::setprecision(2) << n_undet
               << '\t' << std::setprecision(4) << gene_acc.ancient_frac()
               << '\t' << std::setprecision(4) << ci_lo
               << '\t' << std::setprecision(4) << ci_hi
               << '\t' << std::setprecision(4) << gene_acc.avg_posterior()
               << '\t' << std::setprecision(4) << gene_acc.avg_fident()
               << '\t' << std::setprecision(4) << gene_breadth
               << '\t' << std::setprecision(4) << gene_abundance
               << '\t' << std::setprecision(4) << gene_acc.depth_std()
               << '\t' << std::setprecision(4) << gene_acc.depth_evenness()
               << '\t' << gene_acc.n_unique_starts()
               << '\t' << std::setprecision(4) << gene_acc.start_diversity()
               << '\t' << std::setprecision(4) << gene_acc.start_span()
               << '\t' << std::setprecision(4) << gene_pos_score
               << '\t' << std::setprecision(4) << gene_term_ratio
               << '\t' << std::setprecision(2) << gene_acc.avg_aln_len()
               << '\t' << std::setprecision(2) << gene_acc.std_aln_len()
               << '\t' << std::setprecision(4) << cov_deviance
               << '\t' << std::setprecision(4) << cov_fhat
               << '\t' << std::setprecision(4) << cov_weight
               << '\t' << std::setprecision(2) << (assign_eff * cov_weight)
               << '\n';
            ++gene_written;
        }

        // Functional accumulation
        if (!group_map.empty() && gene_passes) {
            auto map_it = group_map.find(target_name);
            if (map_it != group_map.end()) {
                auto& fa = func_acc[map_it->second];
                if (fa.function_id.empty()) { fa.function_id = map_it->second; fa.db_type = "MAP"; }
                const float eff = static_cast<float>(gene_acc.eff_reads());
                const float n_anc = static_cast<float>(gene_acc.n_ancient_conf + gene_acc.n_ancient_likely);
                const float n_mod = static_cast<float>(gene_acc.n_modern_conf);
                const float n_und = static_cast<float>(gene_acc.n_undetermined);
                const float mean_post = gene_acc.avg_posterior();
                const bool gene_damaged = (mean_post >= threshold);
                fa.add_gene(static_cast<uint32_t>(std::round(eff)), n_anc, n_mod, n_und, mean_post, gene_damaged,
                            gene_acc.depth_mean(),
                            gene_acc.sum_fident, gene_acc.sum_fident_sq,
                            gene_acc.sum_aln_len, gene_acc.sum_aln_len_sq,
                            gene_acc.n_reads_eff);
            } else if (verbose && func_acc.empty()) {
                std::cerr << "  Warning: first gene_passes target not found in map: \""
                          << target_name << "\"\n";
            }
        }

        // Anvio KO
        if (af_f.is_open() && gene_passes && !group_map.empty()) {
            auto it = group_map.find(target_name);
            if (it != group_map.end()) {
                af_f << target_name << '\t' << it->second << '\t' << annotation_source
                     << '\t' << std::fixed << std::setprecision(4) << gene_abundance
                     << '\t' << std::setprecision(4) << gene_breadth << '\n';
            }
        }

        // --- Protein output ---
        const bool have_protein = (pa.n_reads > 0 || pa.n_damage_reads > 0);
        const float assign_eff_reads = static_cast<float>(pa.n_reads_eff);
        const float abundance_adj    = assign_eff_reads * cov_weight;
        bool protein_passes = false;
        std::string filter_fail_p;
        if (pf_os.is_open() && have_protein) {
            const float mean_p  = pa.n_auth_reads > 0
                ? (pa.sum_p_damaged / static_cast<float>(pa.n_auth_reads)) : 0.0f;
            const float frac    = pa.n_auth_reads > 0
                ? static_cast<float>(pa.n_damaged) / static_cast<float>(pa.n_auth_reads) : 0.0f;
            const float frac_all= pa.n_reads > 0
                ? static_cast<float>(pa.n_damaged) / static_cast<float>(pa.n_reads) : 0.0f;

            // #2: per-opportunity Beta-Binomial (k damaged / n eligible); null = per-opportunity
            // technical floor q_modern (baseline substitution rate at eligible sites).
            ProteinDamageResult result = fit_protein_damage_binomial(pa.read_obs, modern_est.q_modern, 0.5, 2.0);
            const float p_prot  = static_cast<float>(result.p_damaged);
            const float d_mle   = static_cast<float>(result.delta_max);
            const float ci_lo2  = static_cast<float>(result.ci_lower);
            const float ci_hi2  = static_cast<float>(result.ci_upper);
            const float lrt_p   = static_cast<float>(result.p_value);
            const float log_bf  = static_cast<float>(result.log_bayes_factor);

            // Oxidation posterior: Beta(1+hits, 9+misses) with ~10% prior oxidation rate
            float p_ox_protein = 0.0f;
            if (total_ox_total > 0) {
                const float ox_alpha = 1.0f + static_cast<float>(total_ox_hits);
                const float ox_beta  = 9.0f + (static_cast<float>(total_ox_total) - static_cast<float>(total_ox_hits));
                p_ox_protein = ox_alpha / (ox_alpha + ox_beta);
            }
            // p_prot (= 1 - one-sided LRT p) is a FREQUENTIST detection confidence, p_ox_protein is a
            // Bayesian posterior mean — different quantities, NOT complement-combinable. The old
            // noisy-OR 1-(1-p_prot)(1-p_ox) was illegitimate (treated a 1-p confidence as
            // P(damaged|data)). Report each channel separately and combine only their PASS states:
            // deamination authenticated iff protein_passes (below); oxidation detected iff the
            // sample-level oxidation channel fired AND this protein's posterior clears the prior.
            // any_damage_pass is a boolean any-channel authentication flag, not a probability.
            // (Provisional oxidation gate: sample_ox_detected && p_ox_protein > prior 0.1 — flagged.)

            const float observed= use_both_ends
                ? static_cast<float>(pa.terminal_5 + pa.terminal_3)
                : static_cast<float>(pa.terminal_5);
            const float expected= assign_eff_reads * n_ends_lib *
                static_cast<float>(TERMINAL_CODONS_S) * P_SUSCEPTIBLE_K * P_NONSYNONYMOUS_K;
            const float d_aa    = (expected > 0.0f) ? std::min(1.0f, observed / expected) : 0.0f;

            const float pos_score  = pa.positional_score();
            const float term_ratio = pa.terminal_middle_ratio();
            auto add_fp = [&filter_fail_p](const char* t) {
                if (!filter_fail_p.empty()) filter_fail_p.push_back(',');
                filter_fail_p += t;
            };
            const float prot_breadth = pa.breadth();
            const float prot_depth   = (pa.tlen > 0) ? pa.depth_mean() : 0.0f;
            if (pa.n_reads < min_reads)          add_fp("min_reads");
            if (prot_breadth < min_breadth)      add_fp("min_breadth");
            if (prot_depth < min_depth)          add_fp("min_depth");
            if (min_distinct_starts > 0 && pa.n_reads >= 3 &&
                pa.unique_start_pairs.size() < min_distinct_starts)
                add_fp("min_distinct_starts");
            if (pos_score < min_positional_score) add_fp("min_positional_score");
            if (term_ratio < min_terminal_ratio)  add_fp("min_terminal_ratio");
            // v5 library gate (state-aware): a protein fails only when the library has NO usable
            // pi -- the abstain states NOT_DETECTED/UNDETERMINED/BELOW_FLOOR (or a non-positive
            // point). Usable-advisory states (TRACE/ANCIENT_CPG/LOW_ABUNDANCE) carry a real
            // pi.point and DO pass: the per-read path degrades to the prior-free
            // damage_evidence_llr ranking when read_ancient_posterior refuses, so an authenticated
            // library whose independent-molecule count lands LOW_ABUNDANCE is annotated, not erased.
            // Stratification (DETECTED vs advisory) is a library-level property (pi_state, uniform
            // across the library) and travels to convergence via the .dcx calib metadata.
            if (library_v5 && !library_pi_usable) add_fp("library_not_authenticatable");
            // Q-CLUST ship-gate: fail proteins whose per-opportunity binomial rests on too few
            // informative reads or is dominated by a single molecule (pooled-tail
            // anti-conservative regime). Override via --min-informative-reads /
            // --max-single-molecule-frac.
            if (result.n_informative_reads < min_informative_reads) add_fp("min_informative_reads");
            if (result.max_ni_frac > max_single_molecule_frac)      add_fp("single_molecule_dominated");
            protein_passes = filter_fail_p.empty();
            if (filter_fail_p.empty()) filter_fail_p = ".";

            // Stage 2b: per-protein authentication probability. Prior = library pi_damaged
            // (empirical Bayes), NEVER zeroed by filter_fail — that is a TIER, not a gate.
            // L = logit(prior) + Σ gamma_i·llr_i. The prior enters ONCE (library level); the
            // per-read llr is already prior-free, so there is NO per-read −logit(prior) term.
            // A protein whose reads all sit at the baseline nets Σ=0 → lands on the prior (spike);
            // net-positive terminal deamination pushes above, anti-damage below. Fallback prior of
            // 0.5 MUST surface loudly as auth_tier=no_library (never a silent authentication); we
            // do NOT fall back to pi_point (the downgraded, gated value).
            const bool auth_prior_available = (library_damage_prior > 0.0);
            double auth_prior = auth_prior_available ? library_damage_prior : 0.5;
            auth_prior = std::clamp(auth_prior, 0.02, 0.98);
            const double logit_prior_auth = std::log(auth_prior / (1.0 - auth_prior));
            double auth_L = logit_prior_auth + pa.auth_sum_g_llr;
            auth_L = std::clamp(auth_L, -10.0, 10.0);
            const float p_protein_damaged = 1.0f / (1.0f + std::exp(-static_cast<float>(auth_L)));
            // Heuristic width from the effective sample size (Σ gamma²), NOT a calibrated
            // likelihood CI — reported as a band so it is not read as a confidence interval.
            const double auth_band = std::sqrt(pa.auth_sum_g2);
            const float p_auth_lo = 1.0f / (1.0f + std::exp(
                -static_cast<float>(std::clamp(auth_L - 1.96 * auth_band, -10.0, 10.0))));
            const float p_auth_hi = 1.0f / (1.0f + std::exp(
                -static_cast<float>(std::clamp(auth_L + 1.96 * auth_band, -10.0, 10.0))));
            // Authentication tier: filter status is a TIER modifier, never a zero. The 0.5
            // prior fallback surfaces LOUDLY as no_library (an un-anchored score, never a silent
            // authentication). Strongest to weakest library evidence otherwise; the protein
            // probability above is emitted regardless of the filter.
            std::string auth_tier;
            if (!auth_prior_available)   auth_tier = "no_library";
            else if (!library_pi_usable) auth_tier = "unauthenticatable";
            else                         auth_tier = library_pi_detected ? "detected" : "advisory";
            if (auth_prior_available && library_pi_usable && !protein_passes) auth_tier += "_filtered";

            // --- Binomial provenance LR: PRIMARY per-protein ancient authentication ---
            // logLR = k*log(pi/eps0) + (N-k)*log((1-pi)/(1-eps0)); H1:rate=pi vs H0:rate=eps0.
            // N = independent (em_keep, non-ambiguous) reads carrying reference-free evidence
            // (pa.n_auth_reads); k = those the reference-free MAP call flags damaged (pa.n_damaged).
            // gamma weights only the SECONDARY p_protein_damaged above; the LR uses INTEGER
            // molecule counts (em_keep independence), never gamma. p_ancient is the LR posterior
            // at a flat prior; pvalue is the exact one-sided binomial tail under eps0.
            const int prov_k = static_cast<int>(pa.n_damaged);
            const int prov_n = static_cast<int>(pa.n_auth_reads);
            double prov_logLR = 0.0;
            if (prov_pi_available && prov_n > 0) {
                prov_logLR = prov_k * std::log(prov_pi / prov_eps0)
                           + (prov_n - prov_k) * std::log((1.0 - prov_pi) / (1.0 - prov_eps0));
            }
            const float prov_p_ancient = prov_pi_available
                ? 1.0f / (1.0f + std::exp(-static_cast<float>(std::clamp(prov_logLR, -40.0, 40.0))))
                : 0.0f;
            const double prov_pvalue = prov_pi_available
                ? binom_upper_tail(prov_k, prov_n, prov_eps0) : 1.0;
            // Coverage-aware call: underpowered (N below the min for ~80% power, incl. all low-pi
            // libraries where prov_min_n is huge) → undetermined; else reject H0 (modern) at
            // auth_alpha → ancient, otherwise modern.
            std::string prov_call;
            if (!prov_pi_available || prov_n < prov_min_n) prov_call = "undetermined_low_coverage";
            else if (prov_pvalue <= static_cast<double>(auth_alpha)) prov_call = "ancient";
            else prov_call = "modern";

            // assign_* columns come from the gene accumulator
            const float assign_n_ancient = static_cast<float>(gene_acc.n_ancient_conf + gene_acc.n_ancient_likely);
            const float assign_n_modern  = static_cast<float>(gene_acc.n_modern_conf);
            const float assign_n_undet   = static_cast<float>(gene_acc.n_undetermined);
            const float assign_anc_frac  = gene_acc.ancient_frac();
            const double n_eff_cls = static_cast<double>(assign_n_ancient + assign_n_modern);
            auto [a_ci_lo, a_ci_hi] = wilson_ci_effective(assign_n_ancient, n_eff_cls);

            std::ostringstream row;
            row << pa.protein_id
                << '\t' << pa.n_reads
                << '\t' << std::fixed << std::setprecision(2) << assign_eff_reads
                << '\t' << pa.n_damage_reads
                << '\t' << pa.n_damaged
                << '\t' << pa.total_damage_sites
                << '\t' << pa.terminal_5
                << '\t' << pa.terminal_3
                << '\t' << std::setprecision(4) << pa.max_p_damaged
                << '\t' << std::setprecision(4) << mean_p
                << '\t' << std::setprecision(4) << frac
                << '\t' << std::setprecision(4) << frac_all
                // col13 protein_damage_confidence: GATED — 0 when the protein fails authentication
                // (incl. the Q-CLUST single-molecule gate), so thresholding the primary score cannot
                // resurrect a suppressed verdict. col14 damage_confidence_ungated: raw diagnostic.
                << '\t' << std::setprecision(4) << (protein_passes ? p_prot : 0.0f)
                << '\t' << std::setprecision(4) << p_prot
                << '\t' << std::setprecision(4) << d_mle
                << '\t' << std::setprecision(4) << ci_lo2
                << '\t' << std::setprecision(4) << ci_hi2
                << '\t' << std::scientific << std::setprecision(3) << lrt_p
                << '\t' << std::fixed << std::setprecision(2) << log_bf
                << '\t' << std::setprecision(4) << d_aa
                << '\t' << std::setprecision(1)
                << (assign_eff_reads > 0.0f ? pa.sum_delta_bits / assign_eff_reads : 0.0f)
                << '\t' << std::setprecision(3)
                << (assign_eff_reads > 0.0f ? pa.sum_bpa / assign_eff_reads : 0.0f)
                << '\t' << static_cast<int>(std::max_element(
                    pa.length_bin_counts.begin(), pa.length_bin_counts.end())
                    - pa.length_bin_counts.begin())
                << '\t' << pa.unique_starts.size()
                << '\t' << std::setprecision(4) << pa.start_diversity()
                << '\t' << std::setprecision(4) << pa.start_span()
                << '\t' << std::setprecision(4) << pos_score
                << '\t' << std::setprecision(4) << term_ratio
                << '\t' << std::setprecision(2) << pa.avg_aln_len()
                << '\t' << std::setprecision(2) << pa.std_aln_len()
                << '\t' << gene_acc.n_reads
                << '\t' << std::setprecision(2) << static_cast<float>(gene_acc.eff_reads())
                << '\t' << gene_acc.n_damage_reads
                << '\t' << std::setprecision(2) << static_cast<float>(gene_acc.n_damage_reads_eff)
                << '\t' << std::setprecision(2) << assign_n_ancient
                << '\t' << std::setprecision(2) << assign_n_modern
                << '\t' << std::setprecision(2) << assign_n_undet
                << '\t' << std::setprecision(4) << assign_anc_frac
                << '\t' << std::setprecision(4) << a_ci_lo
                << '\t' << std::setprecision(4) << a_ci_hi
                << '\t' << std::setprecision(4) << gene_acc.avg_posterior()
                << '\t' << std::setprecision(4) << gene_acc.avg_fident()
                << '\t' << std::setprecision(4) << gene_breadth
                << '\t' << std::setprecision(4) << gene_abundance
                << '\t' << std::setprecision(4) << gene_acc.depth_std()
                << '\t' << std::setprecision(4) << gene_acc.depth_evenness()
                << '\t' << std::setprecision(4) << cov_deviance
                << '\t' << std::setprecision(4) << cov_fhat
                << '\t' << std::setprecision(4) << cov_weight
                << '\t' << std::setprecision(2) << abundance_adj
                << '\t' << total_ox_hits
                << '\t' << total_ox_total
                << '\t' << std::setprecision(4) << p_ox_protein
                << '\t' << ((protein_passes || (sample_ox_detected && p_ox_protein > 0.1f)) ? 1 : 0)
                << '\t' << (sample_ox_detected ? 1 : 0)
                << '\t' << std::setprecision(4) << sample_ox_rate_terminal
                << '\t' << (sample_depurination ? 1 : 0)
                << '\t' << std::setprecision(4) << sample_purine_enrich_5p
                << '\t' << std::setprecision(4) << sample_purine_enrich_3p
                << '\t' << (protein_passes ? 1 : 0)
                << '\t' << filter_fail_p
                << '\t' << std::setprecision(4) << result.max_ni_frac
                << '\t' << result.n_informative_reads
                // Stage 2b: authentication probability (never gated to 0), an effective-N band
                // (NOT a calibrated CI), the library prior used, the library DamageConfidence
                // state, and the auth tier.
                << '\t' << std::setprecision(4) << p_protein_damaged
                << '\t' << std::setprecision(4) << p_auth_lo
                << '\t' << std::setprecision(4) << p_auth_hi
                << '\t' << std::setprecision(4) << static_cast<float>(auth_prior)
                << '\t' << static_cast<int>(library_damage_state)
                << '\t' << auth_tier
                // Binomial provenance-LR authentication (PRIMARY): k/N independent-read counts,
                // logLR, LR posterior P(ancient), exact binomial tail p-value, coverage-aware call.
                << '\t' << prov_k
                << '\t' << prov_n
                << '\t' << std::fixed << std::setprecision(3) << prov_logLR
                << '\t' << std::setprecision(4) << prov_p_ancient
                << '\t' << std::scientific << std::setprecision(3) << prov_pvalue
                << '\t' << prov_call
                << '\n';
            if (pf_os.is_open()) { pf_os << row.str(); ++protein_written; }
        }
        ri_cur = ri_end;
        pi_cur = pi_end;
    }

    // 10. Write functional summary.
    if (!func_acc.empty() && !functional_summary_file.empty()) {
        std::ofstream ff_fs(functional_summary_file);
        if (!ff_fs.is_open()) {
            std::cerr << "Error: Cannot open functional summary file: " << functional_summary_file << "\n";
        } else {
            ff_fs << "db\tfunction_id\tlevel\tn_genes\tn_reads\tn_ancient\tn_modern\tn_undetermined\t"
                     "ancient_frac\tci_low\tci_high\tmean_posterior\t"
                     "n_damaged_genes\tdamaged_gene_frac\tdamage_enrichment\t"
                     "coverage_mean\tcoverage_stdev\t"
                     "avg_identity\tstdev_identity\t"
                     "avg_read_length\tstdev_read_length\n";
            uint32_t global_genes = 0, global_damaged = 0;
            for (const auto& [fid, fa] : func_acc) {
                global_genes   += fa.n_genes;
                global_damaged += fa.n_damaged_genes;
            }
            const float global_dam_frac = (global_genes > 0)
                ? static_cast<float>(global_damaged) / static_cast<float>(global_genes) : 0.0f;
            for (const auto& [fid, fa] : func_acc) {
                const float frac_f = fa.ancient_frac();
                auto [ci_l, ci_h] = wilson_ci_effective(
                    static_cast<double>(fa.n_ancient),
                    static_cast<double>(fa.n_ancient + fa.n_modern));
                const float dam_frac = fa.damaged_gene_frac();
                const float enrich   = (global_dam_frac > 0.0f) ? (dam_frac / global_dam_frac) : 0.0f;
                ff_fs << fa.db_type << '\t' << fa.function_id << '\t' << "group" << '\t'
                      << fa.n_genes << '\t' << fa.n_reads << '\t'
                      << std::fixed << std::setprecision(1) << fa.n_ancient << '\t'
                      << fa.n_modern << '\t' << fa.n_undetermined << '\t'
                      << std::setprecision(4) << frac_f << '\t' << ci_l << '\t' << ci_h << '\t'
                      << fa.mean_posterior() << '\t'
                      << fa.n_damaged_genes << '\t'
                      << dam_frac << '\t'
                      << enrich << '\t'
                      << std::setprecision(4) << fa.coverage_mean() << '\t' << fa.coverage_stdev() << '\t'
                      << fa.avg_identity() << '\t' << fa.stdev_identity() << '\t'
                      << std::setprecision(2) << fa.avg_read_length() << '\t' << fa.stdev_read_length() << '\n';
            }
            if (verbose) std::cerr << "Functional summary written to: " << functional_summary_file
                                   << " (" << func_acc.size() << " functions)\n";
        }
    }

    if (verbose) {
        if (gf.is_open())
            std::cerr << "Gene summary written to: " << gene_summary_file << " (" << gene_written << " genes)\n";
        if (pf_os.is_open())
            std::cerr << "Protein summary written to: " << protein_summary_file << " (" << protein_written << " proteins)\n";
        if (af_f.is_open())
            std::cerr << "Anvi'o gene abundance written to: " << anvio_ko_file << "\n";
        std::cerr << "\nSummary: assigned_reads=" << assigned_kept_reads
                  << ", mismatch_reads=" << n_summaries
                  << ", reads_with_damage_sites=" << proteins_with_damage << "\n";
        const auto run_end = std::chrono::steady_clock::now();
        std::cerr << "Runtime: "
                  << dart::log_utils::format_elapsed(run_start, run_end) << "\n";
    }

    return 0;
}

// Register subcommand
namespace {
    struct DamageAnnotateRegistrar {
        DamageAnnotateRegistrar() {
            SubcommandRegistry::instance().register_command(
                "damage-annotate",
                "Post-mapping damage annotation from EMI alignments",
                cmd_damage_annotate, 40);
        }
    } damage_annotate_registrar;
}

}  // namespace cli
}  // namespace dart
