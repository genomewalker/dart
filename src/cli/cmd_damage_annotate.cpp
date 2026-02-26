// dart damage-annotate: Post-mapping damage annotation
//
// Given EMI alignments (with qaln/taln), identifies amino acid positions
// where the query (observed protein) differs from the target (reference)
// in a damage-consistent way. C->T deamination (5' end) and G->A (3' end)
// cause specific amino acid substitution patterns that can be identified
// with ~90% precision when a reference protein is available.

#include "subcommand.hpp"
#include "dart/version.h"
#include "dart/damage_stats.hpp"
#include "dart/damage_index_reader.hpp"
#include "dart/bayesian_damage_score.hpp"
#include "dart/em_reassign.hpp"
#include "dart/theta_d_coordinator.hpp"
#include "dart/coverage_em.hpp"
#include "dart/columnar_index.hpp"
#include "dart/log_utils.hpp"
#include "dart/mmap_array.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
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

// Damage-consistent amino acid substitutions
// target = reference (undamaged), query = observed (potentially damaged)
struct DamageSubstitution {
    char target_aa;
    char query_aa;
    char damage_class;  // 'C' for C->T, 'G' for G->A
    bool high_confidence;
};

static constexpr DamageSubstitution DAMAGE_SUBS[] = {
    // C->T damage (5' end)
    // Codon position 1: C->T
    {'R', 'W', 'C', true},    // CGG->TGG
    {'R', 'C', 'C', true},    // CGC->TGC, CGT->TGT
    {'R', '*', 'C', true},    // CGA->TGA
    {'Q', '*', 'C', true},    // CAA->TAA, CAG->TAG
    {'H', 'Y', 'C', true},    // CAC->TAC, CAT->TAT
    {'P', 'S', 'C', false},   // CCN->TCN
    {'P', 'L', 'C', false},   // CCG->CTG
    // Codon position 2: C->T
    {'T', 'I', 'C', false},   // ACN->ATN (except ACG->ATG)
    {'T', 'M', 'C', false},   // ACG->ATG
    {'A', 'V', 'C', false},   // GCN->GTN
    {'S', 'F', 'C', false},   // TCC->TTC, TCT->TTT
    {'S', 'L', 'C', false},   // TCA->TTA, TCG->TTG

    // G->A damage (3' end)
    // Codon position 1: G->A
    {'E', 'K', 'G', true},    // GAA/GAG->AAA/AAG
    {'D', 'N', 'G', true},    // GAC/GAT->AAC/AAT
    {'A', 'T', 'G', false},   // GCN->ACN
    {'G', 'R', 'G', false},   // GGA/GGG->AGA/AGG
    {'G', 'S', 'G', false},   // GGC/GGT->AGC/AGT
    {'V', 'I', 'G', false},   // GTA/GTC/GTT->ATA/ATC/ATT
    {'V', 'M', 'G', false},   // GTG->ATG
    // Codon position 2: G->A
    {'R', 'K', 'G', true},    // AGA/AGG->AAA/AAG
    {'R', 'Q', 'G', true},    // CGA/CGG->CAA/CAG
    {'R', 'H', 'G', true},    // CGC/CGT->CAC/CAT
    {'G', 'E', 'G', true},    // GGA/GGG->GAA/GAG
    {'G', 'D', 'G', true},    // GGC/GGT->GAC/GAT
    {'S', 'N', 'G', false},   // AGC/AGT->AAC/AAT
    {'C', 'Y', 'G', false},   // TGC/TGT->TAC/TAT
    {'W', '*', 'G', false},   // TGG->TAG/TGA

    // X-masked positions (from search_protein): X was originally damaged AA
    // DART X-masks: (1) damage-convertible stops, (2) W/C likely from R at 5' terminus
    // High-confidence because DART only X-masks when damage probability is high
    {'Q', 'X', 'C', true},    // CAA/CAG->TAA/TAG (stop X-masked)
    {'R', 'X', 'C', true},    // CGA->TGA (stop) OR CGG->TGG (W) OR CGC/CGT->TGC/TGT (C)
    {'W', 'X', 'G', true},    // TGG->TGA/TAG (stop X-masked)
};

static constexpr size_t NUM_DAMAGE_SUBS = sizeof(DAMAGE_SUBS) / sizeof(DAMAGE_SUBS[0]);

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
    // Stored for corrected protein output (empty unless --corrected-proteins requested)
    std::string qaln;
    std::string taln;
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
    float p_damaged;
    uint32_t damage_consistent;
    float info_sum;       // sum(site.p_damage) → ReadDamageObs.info
    uint32_t terminal_5;  // C-terminal 5' site count → pa.terminal_5
    uint32_t terminal_3;  // G-terminal 3' site count → pa.terminal_3
};

// Compact per-read assignment data sorted by ref_idx for streaming gene/protein output.
// Replaces the gene_agg and protein_agg unordered_maps with a sort-then-stream pass.
struct ReadRefEntry {
    uint32_t ref_idx;
    uint32_t tstart_0;
    float gamma;
    float fident;
    float bits;
    float bits_second;
    uint32_t aln_len_on_target;
    uint16_t tlen;
    uint8_t length_bin;
    bool has_damage_obs;
    float posterior;
    uint8_t classification;  // dart::AncientClassification cast to uint8_t
    uint8_t pad[2];
};  // 40 bytes

// Parse strand and frame from DART header suffix (e.g., readname_+_1 -> '+', 1)
// Returns {strand, frame} where frame is 0, 1, or 2
static std::pair<char, int> parse_strand_and_frame(const std::string& query_id) {
    char strand = '+';
    int frame = 0;

    if (query_id.size() >= 4) {
        size_t n = query_id.size();
        for (size_t i = n - 1; i >= 2; --i) {
            if (query_id[i] >= '0' && query_id[i] <= '2') {
                frame = query_id[i] - '0';
                continue;
            }
            if (query_id[i] == '_' && i >= 2 &&
                (query_id[i-1] == '+' || query_id[i-1] == '-') &&
                query_id[i-2] == '_') {
                strand = query_id[i-1];
                break;
            }
            if (query_id[i] != '_') break;
        }
    }
    return {strand, frame};
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
    float d_max, float lambda)
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
            // NOTE: AA_SCALE is NOT applied here - it's only for the fixed-average fallback.
            // The position-weighted formula uses raw d_max because we're computing
            // P(damage at this position | ancient) directly from the damage model.
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
    constexpr float PRIOR_DAMAGED = 0.5f;          // Prior P(read is damaged)

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
    return out;
}

static void recompute_summary_site_metrics(ProteinDamageSummary& summary) {
    constexpr float Q_BASELINE = 0.005f;
    // Bayesian LLR parameters (must match compute_damage_summary)
    constexpr float SEQ_ERROR_RATE = 0.01f;
    constexpr float MIN_DAMAGE_RATE = 0.01f;       // Interior sites: LLR=0 (no evidence)
    constexpr float PRIOR_DAMAGED = 0.5f;

    float llr_sum = 0.0f;
    summary.ct_sites = 0;
    summary.ga_sites = 0;
    summary.high_conf_sites = 0;
    summary.max_p_damage = 0.0f;
    summary.aa_k_hits = 0;
    summary.aa_sum_log_qA_hits = 0.0f;

    for (const auto& s : summary.sites) {
        if (s.damage_class == 'C') summary.ct_sites++;
        else summary.ga_sites++;
        if (s.high_confidence) summary.high_conf_sites++;
        summary.max_p_damage = std::max(summary.max_p_damage, s.p_damage);
        summary.aa_k_hits++;
        summary.aa_sum_log_qA_hits += std::log((s.p_damage + Q_BASELINE) / Q_BASELINE);

        // Bayesian LLR for each site
        float p_site_given_damaged = std::max(s.p_damage, MIN_DAMAGE_RATE);
        float llr = std::log(p_site_given_damaged / SEQ_ERROR_RATE);
        llr = std::clamp(llr, -3.0f, 5.0f);
        llr_sum += llr;
    }

    if (!summary.sites.empty()) {
        float logit_prior = std::log(PRIOR_DAMAGED / (1.0f - PRIOR_DAMAGED));
        float logit_posterior = logit_prior + llr_sum;
        logit_posterior = std::clamp(logit_posterior, -10.0f, 10.0f);
        summary.p_damaged = 1.0f / (1.0f + std::exp(-logit_posterior));
    } else {
        summary.p_damaged = 0.0f;
    }
}

// Filter sites by distance from relevant terminus
static void filter_sites_by_distance(ProteinDamageSummary& summary, int max_dist) {
    if (max_dist < 0) return;

    std::vector<DamageSite> filtered;
    for (auto& site : summary.sites) {
        size_t dist = (site.damage_class == 'C') ? site.dist_5prime : site.dist_3prime;
        if (dist <= static_cast<size_t>(max_dist)) {
            filtered.push_back(std::move(site));
        }
    }

    summary.sites = std::move(filtered);
    summary.damage_consistent = summary.sites.size();
    recompute_summary_site_metrics(summary);
    summary.non_damage_mismatches = summary.total_mismatches - summary.damage_consistent;
    summary.damage_fraction = summary.total_mismatches > 0
        ? static_cast<float>(summary.damage_consistent) / static_cast<float>(summary.total_mismatches)
        : 0.0f;
}

// Generate corrected protein by replacing damage sites with reference AAs
// Uses sorted index walk instead of hash set (sites are ordered by query_pos)
static std::string generate_corrected_protein(
    const std::string& qaln, const std::string& taln,
    const std::vector<DamageSite>& sites,
    size_t qstart)
{
    size_t site_idx = 0;
    std::string corrected;
    corrected.reserve(qaln.size());
    size_t q_pos = qstart;

    for (size_t i = 0; i < qaln.size() && i < taln.size(); ++i) {
        char q_aa = qaln[i];
        char t_aa = taln[i];

        if (q_aa == '-') continue;
        if (t_aa == '-') {
            corrected += q_aa;
            q_pos++;
            continue;
        }

        // Advance site_idx to match current position
        while (site_idx < sites.size() && sites[site_idx].query_pos < q_pos) {
            site_idx++;
        }

        if (site_idx < sites.size() && sites[site_idx].query_pos == q_pos) {
            corrected += t_aa;  // Replace with reference AA
            site_idx++;
        } else {
            corrected += q_aa;  // Keep observed
        }
        q_pos++;
    }

    return corrected;
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
    size_t n_damage_reads = 0;    // Assigned reads with mismatch evidence
    size_t n_damaged = 0;         // Damage-evidence reads with p_damaged > 0.5
    size_t total_damage_sites = 0; // Sum of damage_consistent across reads
    size_t terminal_5 = 0;        // C→T sites at 5' end (first 3 codons)
    size_t terminal_3 = 0;        // G→A sites at 3' end (first 3 codons)
    float max_p_damaged = 0.0f;   // Maximum p_damaged across reads
    float sum_p_damaged = 0.0f;   // Sum of p_damaged (for averaging)

    // Per-read observations for Beta-Binomial model
    std::vector<ReadDamageObs> read_obs;  // Collected for proper statistical inference

    // Short protein filtering aggregates
    float sum_delta_bits = 0.0f;
    float sum_bpa = 0.0f;
    std::array<uint32_t, 4> length_bin_counts{};  // per-bin read counts

    // Mapping-pattern aggregates for spurious-hit filtering
    std::unordered_set<uint32_t> unique_starts;
    uint32_t min_start = UINT32_MAX;
    uint32_t max_start = 0;
    double left_cov = 0.0;   // Alignment overlap in left terminal window
    double right_cov = 0.0;  // Alignment overlap in right terminal window
    double middle_cov = 0.0; // Alignment overlap in middle region
    double sum_aln_len = 0.0;
    double sum_aln_len_sq = 0.0;

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
        sum_fident += gamma * fident;
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

    void add_gene(uint32_t gene_reads, float gene_ancient, float gene_modern,
                  float gene_undetermined, float gene_mean_posterior,
                  bool gene_is_damaged) {
        ++n_genes;
        if (gene_is_damaged) ++n_damaged_genes;
        n_reads += gene_reads;
        n_ancient += gene_ancient;
        n_modern += gene_modern;
        n_undetermined += gene_undetermined;
        sum_posterior += gene_mean_posterior * gene_reads;
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
    std::string corrected_file;
    std::string protein_summary_file;  // Per-protein aggregation
    std::string protein_filtered_file; // Per-protein filtered subset
    std::string combined_output_file;  // Combined per-target protein+gene summary
    std::string blast8_unique_file;    // BLAST8-like export with unique mappers only
    std::string analysis_prefix;       // Output bundle prefix
    std::string analysis_proteins_file;   // Explicit name for analysis proteins output
    std::string analysis_blast8_file;     // Explicit name for analysis blast8 unique output
    std::string analysis_categories_file; // Explicit name for analysis categories output
    std::string damage_index_file;     // Binary damage index from predict
    float d_max = -1.0f;  // -1 means "estimate from data"
    float lambda = 0.3f;
    bool lambda_set_by_user = false;
    bool refine_damage = false;
    int max_dist = -1;
    std::string lib_type = "auto";    // ss, ds, or auto-detect
    float threshold = 0.7f;           // Classification threshold for is_damaged
    bool use_identity = true;         // Identity evidence enabled by default (+7.6% AUC)
    bool verbose = false;

    // Gene summary options (work with or without --em)
    std::string gene_summary_file;
    float min_breadth = 0.10f;
    float min_depth = 0.5f;
    uint32_t min_reads = 3;
    float min_positional_score = 0.0f;
    float min_terminal_ratio = 0.0f;
    bool user_set_terminal_ratio = false;  // True when --min-terminal-ratio is explicitly provided
    bool auto_calibrate_spurious = false;  // Auto-derive thresholds from data

    // EM reassignment options (enabled by default for abundance estimation)
    bool use_em = true;
    bool em_streaming = false;       // Use streaming EM (O(num_refs) memory instead of O(num_alns))
    size_t em_max_memory_mb = 4096;  // Max memory for EM (MB), 0 = no limit, auto-switch to streaming if exceeded
    uint32_t em_max_iters = 100;
    double em_lambda_b = 3.0;
    double em_tol = 1e-4;
    double em_min_prob = 1e-6;       // absolute posterior threshold for keeping alignment

    // Coverage-aware outer EM options
    bool use_coverage_em = false;
    uint32_t coverage_em_iters = 5;
    float coverage_em_tau = 0.13f;
    uint32_t coverage_em_bins = 6;

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
    float depth_gate_lo = 8.0f;       // --depth-gate-lo
    float depth_gate_hi = 20.0f;      // --depth-gate-hi
    float w_min         = 0.20f;      // --w-min

    // Alignment-level pre-filters (applied before any processing)
    float aln_min_identity = 0.0f;   // 0 = no filter
    float aln_min_bits = 0.0f;       // 0 = no filter
    float aln_max_evalue = 1e10f;    // very large = no filter
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
        } else if (strcmp(argv[i], "--max-dist") == 0 && i + 1 < argc) {
            max_dist = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--sites") == 0 && i + 1 < argc) {
            sites_file = argv[++i];
        } else if (strcmp(argv[i], "--corrected") == 0 && i + 1 < argc) {
            corrected_file = argv[++i];
        } else if (strcmp(argv[i], "--protein-summary") == 0 && i + 1 < argc) {
            protein_summary_file = argv[++i];
        } else if (strcmp(argv[i], "--protein-filtered") == 0 && i + 1 < argc) {
            protein_filtered_file = argv[++i];
        } else if (strcmp(argv[i], "--combined-output") == 0 && i + 1 < argc) {
            combined_output_file = argv[++i];
        } else if (strcmp(argv[i], "--blast8-unique") == 0 && i + 1 < argc) {
            blast8_unique_file = argv[++i];
        } else if (strcmp(argv[i], "--analysis-prefix") == 0 && i + 1 < argc) {
            analysis_prefix = argv[++i];
        } else if (strcmp(argv[i], "--analysis-proteins") == 0 && i + 1 < argc) {
            analysis_proteins_file = argv[++i];
        } else if (strcmp(argv[i], "--analysis-blast8") == 0 && i + 1 < argc) {
            analysis_blast8_file = argv[++i];
        } else if (strcmp(argv[i], "--analysis-categories") == 0 && i + 1 < argc) {
            analysis_categories_file = argv[++i];
        } else if ((strcmp(argv[i], "--library-type") == 0 || strcmp(argv[i], "--lib-type") == 0) && i + 1 < argc) {
            lib_type = argv[++i];
        } else if (strcmp(argv[i], "--damage-index") == 0 && i + 1 < argc) {
            damage_index_file = argv[++i];
        } else if (strcmp(argv[i], "--threshold") == 0 && i + 1 < argc) {
            threshold = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--use-identity") == 0) {
            use_identity = true;
        } else if (strcmp(argv[i], "--no-identity") == 0) {
            use_identity = false;
        } else if (strcmp(argv[i], "--refine-damage") == 0) {
            refine_damage = true;
        } else if (strcmp(argv[i], "--em") == 0) {
            use_em = true;
        } else if (strcmp(argv[i], "--no-em") == 0) {
            use_em = false;
        } else if (strcmp(argv[i], "--em-iters") == 0 && i + 1 < argc) {
            em_max_iters = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--em-lambda") == 0 && i + 1 < argc) {
            em_lambda_b = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--em-tol") == 0 && i + 1 < argc) {
            em_tol = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--em-min-prob") == 0 && i + 1 < argc) {
            em_min_prob = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--coverage-em") == 0) {
            use_coverage_em = true;
        } else if (strcmp(argv[i], "--no-coverage-em") == 0) {
            // Backward-compatible alias; coverage EM is now opt-in.
            use_coverage_em = false;
        } else if (strcmp(argv[i], "--coverage-em-iters") == 0 && i + 1 < argc) {
            coverage_em_iters = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--coverage-em-tau") == 0 && i + 1 < argc) {
            coverage_em_tau = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--coverage-em-bins") == 0 && i + 1 < argc) {
            coverage_em_bins = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--depth-gate-lo") == 0 && i + 1 < argc) {
            depth_gate_lo = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--depth-gate-hi") == 0 && i + 1 < argc) {
            depth_gate_hi = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--w-min") == 0 && i + 1 < argc) {
            w_min = std::stof(argv[++i]);
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
        } else if (strcmp(argv[i], "--preset") == 0 && i + 1 < argc) {
            const char* preset = argv[++i];
            if (strcmp(preset, "loose") == 0) {
                min_reads = 1; min_breadth = 0.0f; min_depth = 0.0f; min_damage_sites = 1;
            } else if (strcmp(preset, "strict") == 0) {
                min_reads = 5; min_breadth = 0.20f; min_depth = 1.0f;
                auto_calibrate_spurious = true; min_damage_sites = 5;
            } else {
                std::cerr << "Error: --preset must be 'loose' or 'strict'\n"; return 1;
            }
        } else if (strcmp(argv[i], "--auto-prior-ancient") == 0) {
            auto_prior_ancient = true;
        } else if (strcmp(argv[i], "--em-streaming") == 0) {
            em_streaming = true;
        } else if ((strcmp(argv[i], "--em-max-memory") == 0 || strcmp(argv[i], "-m") == 0) && i + 1 < argc) {
            em_max_memory_mb = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--gene-summary") == 0 && i + 1 < argc) {
            gene_summary_file = argv[++i];
        } else if (strcmp(argv[i], "--min-breadth") == 0 && i + 1 < argc) {
            min_breadth = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--min-depth") == 0 && i + 1 < argc) {
            min_depth = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--min-reads") == 0 && i + 1 < argc) {
            min_reads = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--min-positional-score") == 0 && i + 1 < argc) {
            min_positional_score = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--min-terminal-ratio") == 0 && i + 1 < argc) {
            min_terminal_ratio = std::stof(argv[++i]);
            user_set_terminal_ratio = true;
        } else if (strcmp(argv[i], "--auto-calibrate-spurious") == 0) {
            auto_calibrate_spurious = true;
        } else if (strcmp(argv[i], "--aln-min-identity") == 0 && i + 1 < argc) {
            aln_min_identity = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--aln-min-bits") == 0 && i + 1 < argc) {
            aln_min_bits = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--aln-max-evalue") == 0 && i + 1 < argc) {
            aln_max_evalue = std::stof(argv[++i]);
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
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            std::cout << "Usage: dart damage-annotate --emi <hits.emi> [options]\n\n";
            std::cout << "Post-mapping damage annotation using EMI alignments.\n";
            std::cout << "Compares observed proteins against reference proteins to identify\n";
            std::cout << "damage-consistent amino acid substitutions with ~90% precision.\n\n";
            std::cout << "Required:\n";
            std::cout << "  --emi, -i FILE    EMI v5 index with qaln/taln + damage evidence columns\n\n";
            std::cout << "Options:\n";
            std::cout << "  --d-max RATE      Override damage rate (default: estimated from data)\n";
            std::cout << "  --lambda FLOAT    Decay constant (default: 0.3)\n";
            std::cout << "  --max-dist INT    Positional filter: max nt from relevant terminus (-1=off)\n";
            std::cout << "  -o FILE           Per-read damage summary TSV (default: stdout)\n";
            std::cout << "  --protein-summary FILE  Per-protein aggregated damage TSV\n";
            std::cout << "  --protein-filtered FILE Per-protein filtered subset TSV\n";
            std::cout << "  --combined-output FILE  Per-target combined protein+gene summary TSV\n";
            std::cout << "  --blast8-unique FILE  BLAST8-like hits for unique mappers only\n";
            std::cout << "  --analysis-prefix STR  Output prefix: <STR>.reads.tsv, <STR>.protein.tsv,\n";
            std::cout << "                         <STR>.blast8.unique.tsv, and <STR>.categories.tsv when --map is provided\n";
            std::cout << "  --analysis-proteins FILE   Proteins output file\n";
            std::cout << "  --analysis-blast8 FILE     Unique BLAST8 output file\n";
            std::cout << "  --analysis-categories FILE Categories output file\n";
            std::cout << "  --library-type TYPE     Library type: ss, ds, or auto (default: auto)\n";
            std::cout << "  --sites FILE      Per-site damage calls TSV\n";
            std::cout << "  --corrected FILE  Reference-guided corrected proteins FASTA\n";
            std::cout << "  --damage-index FILE  Binary damage index (.agd) from predict\n";
            std::cout << "                       Enables synonymous damage detection\n";
            std::cout << "  --threshold FLOAT Damage classification threshold (default: 0.7)\n";
            std::cout << "                       Adds is_damaged column (1 if posterior >= threshold)\n";
            std::cout << "  --no-identity     Disable identity evidence (enabled by default)\n";
            std::cout << "  --refine-damage       Jointly update d_max from EM assignments (Phase 2)\n";
            std::cout << "  -v                Verbose output\n\n";
            std::cout << "Gene summary (per-gene aggregation, works with or without --em):\n";
            std::cout << "  --gene-summary FILE  Per-gene statistics TSV\n";
            std::cout << "  --min-breadth FLOAT  Minimum breadth filter (default: 0.1)\n";
            std::cout << "  --min-depth FLOAT    Minimum depth filter (default: 0.5)\n";
            std::cout << "  --min-reads INT      Minimum read count (default: 3)\n\n";
            std::cout << "Functional profiling (requires --map):\n";
            std::cout << "  --map FILE           Gene-to-group mapping TSV (gene_id<TAB>group)\n";
            std::cout << "  --functional-summary FILE  Per-function stats TSV\n";
            std::cout << "  --anvio-ko FILE      Gene-level group abundance for anvi-estimate-metabolism\n";
            std::cout << "  --annotation-source STR  Source label for anvi output (default: DART)\n\n";
            std::cout << "EM reassignment (multi-mapping resolution, enabled by default):\n";
            std::cout << "  --no-em           Disable EM (use best-hit only)\n";
            std::cout << "  --em-iters INT    Max EM iterations (default: 100)\n";
            std::cout << "  --em-lambda FLOAT Temperature parameter (default: 3.0)\n";
            std::cout << "  --em-tol FLOAT    Convergence tolerance (default: 1e-4)\n";
            std::cout << "  --em-min-prob FLOAT     Min EM posterior to keep (default: 1e-6)\n";
            std::cout << "\nCoverage-aware EM (outer loop, opt-in):\n";
            std::cout << "  --coverage-em            Enable coverage-aware outer EM loop (default: off)\n";
            std::cout << "  --no-coverage-em         Explicitly disable coverage-aware outer EM (kept for scripts)\n";
            std::cout << "  --coverage-em-iters INT  Max outer EM iterations (default: 5)\n";
            std::cout << "  --coverage-em-tau FLOAT  KL penalty temperature (default: 0.13)\n";
            std::cout << "  --coverage-em-bins INT   Coverage bins per protein (default: 6)\n";
            std::cout << "  --depth-gate-lo FLOAT    Coverage depth gate lower bound (default: 8.0)\n";
            std::cout << "  --depth-gate-hi FLOAT    Coverage depth gate upper bound (default: 20.0)\n";
            std::cout << "  --w-min FLOAT            Coverage weight floor (default: 0.20)\n";
            std::cout << "\nBayesian scoring:\n";
            std::cout << "  --prior-ancient FLOAT    Prior P(ancient) for Bayesian scorer (default: 0.10)\n";
            std::cout << "  --auto-prior-ancient     Auto-calibrate prior from mean p_read distribution\n";
            std::cout << "                           Requires AGD index (--damage-index). Uses E[p_read]\n";
            std::cout << "                           as an empirical Bayes estimate of fraction ancient.\n";
            std::cout << "  --min-damage-sites INT   Min susceptible AA positions to use site evidence\n";
            std::cout << "                           (default: auto = max(1, median_qlen/10))\n";
            std::cout << "  --ancient-threshold FLOAT  Ancient posterior cutoff (default: 0.60)\n";
            std::cout << "  --modern-threshold FLOAT   Modern posterior cutoff (default: 0.25)\n";
            std::cout << "  --terminal-threshold FLOAT Min p_read for terminal evidence (default: 0.50)\n";
            std::cout << "  --site-cap FLOAT           Max |logBF| from AA site evidence (default: 3.0)\n";
            std::cout << "  --w0 FLOAT                 Absence-evidence tempering weight 0-1 (default: 0.30)\n";
            std::cout << "\nPresets (applied immediately; individual flags after --preset override):\n";
            std::cout << "  --preset loose   Maximize recall: min-reads=1, breadth=0, depth=0, min-damage-sites=1\n";
            std::cout << "  --preset strict  Maximize precision: min-reads=5, breadth=0.2, depth=1.0,\n";
            std::cout << "                   auto-calibrate-spurious, min-damage-sites=5\n";
            std::cout << "  --em-streaming    Use streaming EM (low memory, ~2.5x slower)\n";
            std::cout << "  --em-max-memory, -m MB  Max EM memory in MB (default: 4096)\n";
            std::cout << "                       Auto-switches to streaming if estimate exceeds limit\n";
            std::cout << "  --min-positional-score FLOAT  Min read start diversity (default: 0)\n";
            std::cout << "                       Filters spurious matches where all reads hit same position\n";
            std::cout << "                       Score = sqrt(diversity * span), range 0-1\n";
            std::cout << "  --min-terminal-ratio FLOAT  Min terminal/middle coverage ratio (default: 0)\n";
            std::cout << "                       Filters matches with reads only in middle (spurious motif)\n";
            std::cout << "  --auto-calibrate-spurious   Auto-derive thresholds from data\n";
            std::cout << "                       Uses 5th percentile of well-supported proteins (n_reads>=10, n_unique_starts>=3)\n\n";
            std::cout << "Alignment-level pre-filters (applied before any processing):\n";
            std::cout << "  --aln-min-identity FLOAT  Min identity fraction (default: 0, no filter)\n";
            std::cout << "  --aln-min-bits FLOAT      Min bit score (default: 0, no filter)\n";
            std::cout << "  --aln-max-evalue FLOAT    Max e-value (default: 1e10, no filter)\n\n";
            std::cout << "Parallelism:\n";
            std::cout << "  --threads, -t INT    OpenMP threads for EMI scans (default: runtime default)\n\n";
            std::cout << "Pipeline:\n";
            std::cout << "  hits2emi -> damage-annotate\n";
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
    if (coverage_em_bins < 1) {
        std::cerr << "Error: --coverage-em-bins must be >= 1\n";
        return 1;
    }
    if (coverage_em_iters < 1) {
        std::cerr << "Error: --coverage-em-iters must be >= 1\n";
        return 1;
    }
    if (coverage_em_tau < 0.0f) {
        std::cerr << "Error: --coverage-em-tau must be >= 0\n";
        return 1;
    }
    if (depth_gate_lo < 0.0f || depth_gate_hi <= depth_gate_lo) {
        std::cerr << "Error: --depth-gate-lo must be >= 0 and < --depth-gate-hi\n";
        return 1;
    }
    if (w_min < 0.0f || w_min > 1.0f) {
        std::cerr << "Error: --w-min must be in [0, 1]\n";
        return 1;
    }
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

    if (!analysis_proteins_file.empty()) {
        protein_summary_file = analysis_proteins_file;
    }
    if (!analysis_blast8_file.empty()) {
        blast8_unique_file = analysis_blast8_file;
    }
    if (!analysis_categories_file.empty()) {
        functional_summary_file = analysis_categories_file;
    }

    if (!analysis_prefix.empty()) {
        if (output_file.empty()) {
            output_file = analysis_prefix + ".reads.tsv";
        }
        if (protein_summary_file.empty()) {
            protein_summary_file = analysis_prefix + ".protein.tsv";
        }
        if (blast8_unique_file.empty()) {
            blast8_unique_file = analysis_prefix + ".blast8.unique.tsv";
        }
        if (!map_file.empty() && functional_summary_file.empty()) {
            functional_summary_file = analysis_prefix + ".categories.tsv";
        }
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
        if (max_dist >= 0) {
            std::cerr << "Positional filter: <=" << max_dist << " nt from terminus\n";
        }
        if (!analysis_prefix.empty()) {
            std::cerr << "analysis-prefix: " << analysis_prefix << "\n";
        }
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

    // Open EMI index
    dart::ColumnarIndexReader reader(emi_file);
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

    if (verbose) {
        std::cerr << "EMI metadata:\n";
        std::cerr << "  Version: " << reader.version() << "\n";
        std::cerr << "  Alignments: " << reader.num_alignments() << "\n";
        std::cerr << "  Reads: " << n_reads << "\n";
        std::cerr << "  References: " << n_refs << "\n";
        std::cerr << "  Sorted (read,score): "
                  << (reader.is_sorted_by_read_then_score() ? "yes" : "no") << "\n";
        if (max_dist >= 0) {
            std::cerr << "  Positional filter: <=" << max_dist << " nt from terminus\n";
        }
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
        uint32_t best_rg = UINT32_MAX;
        uint32_t best_row = UINT32_MAX;
    };
    std::vector<BestHitData> best_hits(n_reads);
    const bool need_unique_degree = !blast8_unique_file.empty();
    std::vector<uint32_t> read_degree;
    if (need_unique_degree) {
        read_degree.assign(n_reads, 0);
    }
    std::vector<dart::CompactAlignment> em_alignments;
    if (use_em && !em_streaming) {
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
    if (use_em && !em_streaming) {
        em_shard_count = static_cast<size_t>(std::max(1, omp_get_max_threads()));
    }
#endif
    std::vector<std::vector<dart::CompactAlignment>> em_alignments_by_thread;
    if (use_em && !em_streaming) {
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
        if (use_em && !em_streaming) {
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
        }  // Note: local_em_read_damage always needed for damage-variation detection

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
        if (local_sorted_hits) {
            local_best_sorted.reserve(std::max<uint32_t>(1, num_rows / 8));
        } else {
            local_best_unsorted.reserve(std::max<uint32_t>(1, num_rows / 4));
        }

        std::unordered_map<uint32_t, uint32_t> local_degree_counts;
        if (need_unique_degree) {
            local_degree_counts.reserve(std::max<uint32_t>(1, num_rows / 8));
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

            if (need_unique_degree) {
                auto it_deg = local_degree_counts.find(ridx);
                if (it_deg == local_degree_counts.end()) {
                    local_degree_counts.emplace(ridx, 1u);
                } else {
                    it_deg->second += 1u;
                }
            }

            if (use_em) {
                dart::CompactAlignment ca{};
                ca.read_idx = ridx;
                ca.ref_idx = tidx;
                ca.bit_score = bits;
                const float ds = std::clamp(damage_score[i], 0.0f, 1.0f);
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
                if (!em_streaming) local_em.push_back(ca);
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
                    cand.damage_score = std::clamp(damage_score[i], 0.0f, 1.0f);
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
                                     bits, std::clamp(damage_score[i], 0.0f, 1.0f),
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
                    cand.damage_score = std::clamp(damage_score[i], 0.0f, 1.0f);
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
                                     bits, std::clamp(damage_score[i], 0.0f, 1.0f),
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

        if (need_unique_degree) {
            for (const auto& kv : local_degree_counts) {
                const uint32_t ridx = kv.first;
                const uint32_t add = kv.second;
                const size_t shard = static_cast<size_t>(ridx) & (kBestHitShards - 1);
                std::lock_guard<std::mutex> lock(best_hit_mutexes[shard]);
                read_degree[ridx] += add;
            }
        }

        if (use_em && local_em_has_rows) {
            if (!em_streaming && !local_em.empty()) {
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

    if (use_em && !em_streaming) {
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

    // Cache all ref names and read names to heap BEFORE flush_pages() evicts NFS pages.
    // flush_pages() calls posix_fadvise+madvise DONTNEED on the entire file; on NFS +
    // Linux 4.18 MAP_PRIVATE, re-faulting evicted pages returns garbage instead of file
    // content under heavy memory pressure. Materialise once while pages are still clean.
    std::vector<std::string> cached_ref_names(reader.num_refs());
    for (uint32_t i = 0; i < reader.num_refs(); ++i)
        cached_ref_names[i] = std::string(reader.ref_name(i));

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

    // Open BLAST8 output file if requested.
    std::ofstream blast8_out;
    size_t written_unique = 0;
    if (!blast8_unique_file.empty()) {
        blast8_out.open(blast8_unique_file);
        if (!blast8_out.is_open()) {
            std::cerr << "Error: Cannot open BLAST8 unique file: " << blast8_unique_file << "\n";
            return 1;
        }
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

    const size_t em_aln_count = em_alignments.size();

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
    };
    std::vector<EmReadResult> em_read_results(n_reads);
    // True when damage-aware EM was active and gamma_ancient was actually populated.
    // Used by --refine-damage to avoid falling back to total gamma for modern reads.
    // Atomic because it may be set inside the parallel streaming_em_finalize callback.
    std::atomic<bool> damage_em_active{false};
    std::unordered_map<uint32_t, dart::ProteinCoverageStats> coverage_stats_map;

    // Auto-switch to streaming EM if estimated memory exceeds limit
    // Memory estimate: ~72 bytes per alignment at SQUAREM peak:
    //   Alignment struct (40B) + gamma (8B) + gamma_ancient (8B)
    //   + gamma_scratch (8B) + gamma_ancient_scratch (8B) = 72B
    // Plus ~4 bytes per read for offsets
    if (use_em && !em_streaming && em_max_memory_mb > 0) {
        const size_t num_alignments = reader.num_alignments();
        const size_t estimated_mb = (num_alignments * 72 + n_reads * 4) / (1024 * 1024);
        if (estimated_mb > em_max_memory_mb) {
            em_streaming = true;
            if (verbose) {
                std::cerr << "\nAuto-switching to streaming EM:\n";
                std::cerr << "  Estimated memory: " << estimated_mb << " MB\n";
                std::cerr << "  Memory limit: " << em_max_memory_mb << " MB\n";
                std::cerr << "  Use --em-max-memory to adjust limit, or --em-streaming to force\n";
            }
        }
    }

    if (use_em && em_streaming) {
        // Streaming EM path - O(num_refs) memory instead of O(num_alignments)
        if (verbose) {
            std::cerr << "\nEM reassignment (streaming mode):\n";
            std::cerr << "  Alignments: " << reader.num_alignments() << " (from EMI)\n";
            std::cerr << "  Reads: " << n_reads << "\n";
            std::cerr << "  References: " << n_refs << "\n";
            std::cerr << "  Lambda_b: " << em_lambda_b << "\n";
            std::cerr << "  Max iters: " << em_max_iters << "\n";
            std::cerr << "  Note: Streaming mode uses all EMI alignments (ignores pre-filters)\n";
        }

        dart::EMParams em_params;
        em_params.lambda_b = em_lambda_b;
        em_params.max_iters = em_max_iters;
        em_params.tol = em_tol;
        em_params.use_squarem = false;  // Streaming doesn't support SQUAREM yet
        em_params.use_damage = true;
        em_params.normalize_by_length = true;
        // Enable alignment damage likelihoods if available (checked via NULL in E-step)
        em_params.use_alignment_damage_likelihood = true;

        const auto em_start = std::chrono::steady_clock::now();
        // Flush all NFS-backed EMI pages accumulated during Pass 1 and scoring.
        // Without this, ~80-150 GB of page-cache pages from the sequential scan
        // stay resident, consuming headroom before streaming_em even starts.
        reader.flush_pages();

        // Restrict to the 7 columns streaming_em actually reads.
        // set_all_columns() would include QALN/TALN (~13 MB/row-group) in
        // MADV_WILLNEED prefetch, inflating NFS RSS by ~50-100 GB over 100+
        // EM iterations. Scoring passes that need QALN/TALN run before this
        // block and restore their own column set; streaming_em_finalize uses
        // the same 7 columns. No restore needed — reader not used after finalize.
        reader.set_columns({
            dart::ColumnID::READ_IDX,
            dart::ColumnID::REF_IDX,
            dart::ColumnID::BIT_SCORE,
            dart::ColumnID::DAMAGE_SCORE,
            dart::ColumnID::DMG_LL_A,
            dart::ColumnID::DMG_LL_M,
            dart::ColumnID::TLEN,
        });

        if (verbose) {
            std::cerr << "  Starting streaming EM"
                      << (use_coverage_em ? " (with coverage outer loop)" : "")
                      << "...\n";
            std::cerr << std::flush;
        }

        // Coverage EM accumulation now happens after streaming_em_finalize (below).
        // Always run standalone streaming_em; build coverage stats from sorted CovEntry
        // after per-read gammas are available. Eliminates per-thread unordered_maps
        // (source of ~400 GB heap fragmentation in coverage_em_outer_loop).
        dart::CoverageEMParams cov_params_streaming;
        if (use_coverage_em) {
            cov_params_streaming.max_outer_iters = coverage_em_iters;
            cov_params_streaming.tau = coverage_em_tau;
            cov_params_streaming.default_bins = coverage_em_bins;
            cov_params_streaming.depth_gate_lo = depth_gate_lo;
            cov_params_streaming.depth_gate_hi = depth_gate_hi;
            cov_params_streaming.w_min = w_min;
        }

        dart::StreamingEMResult streaming_result;
        streaming_result = dart::streaming_em(reader, em_params,
            [verbose](const dart::EMIterationDiagnostics& d) {
                if (verbose && (d.iteration % 10 == 0 || d.iteration < 5)) {
                    std::cerr << "    iter " << d.iteration
                              << " ll=" << std::fixed << std::setprecision(2) << d.log_likelihood
                              << "\n";
                    std::cerr << std::flush;
                }
                // Return per-row-group Parquet heap back to OS after each EM iteration.
                // Without this, 100 sequential EMI scans accumulate ~150 GB of fragmented
                // glibc arena pages that MADV_DONTNEED alone cannot reclaim on NFS mounts.
#ifdef __linux__
                malloc_trim(0);
#endif
            });

        const auto em_end = std::chrono::steady_clock::now();

        if (verbose) {
            std::cerr << "  Converged in " << streaming_result.iterations << " iterations\n";
            std::cerr << "  Log-likelihood: " << std::fixed << std::setprecision(2)
                      << streaming_result.log_likelihood << "\n";
            std::cerr << "  Pi (ancient): " << std::fixed << std::setprecision(4)
                      << streaming_result.pi << "\n";
            std::cerr << "  Memory: O(" << n_refs << " refs) vs O("
                      << reader.num_alignments() << " alignments)\n";
            std::cerr << "  Runtime: "
                      << dart::log_utils::format_elapsed(em_start, em_end) << "\n";
        }

        // Copy weights to em_state for downstream use
        em_state.weights = streaming_result.weights;
        em_state.pi = streaming_result.pi;
        em_state.log_likelihood = streaming_result.log_likelihood;
        em_state.iterations = streaming_result.iterations;

        // Compute per-read stats via streaming finalize
        // Process gamma per row group without storing all values
        // Note: parallel_scan runs in parallel, need mutex for shared updates
        std::vector<float> read_top1(n_reads, 0.0f);
        std::vector<float> read_top2(n_reads, 0.0f);
        std::vector<double> read_sum_sq(n_reads, 0.0);
        constexpr size_t kFinalizeShards = 4096;
        std::array<std::mutex, kFinalizeShards> finalize_mutexes;

        dart::streaming_em_finalize(reader, streaming_result, em_params,
            [&](uint32_t /*rg_idx*/, uint32_t num_rows,
                const uint32_t* read_idx, const uint32_t* ref_idx,
                const float* gamma, const float* gamma_ancient)
            {
                if (gamma_ancient) damage_em_active = true;
                for (uint32_t i = 0; i < num_rows; ++i) {
                    const uint32_t r = read_idx[i];
                    if (r >= n_reads) continue;
                    const float g = gamma[i];
                    const float g_a = gamma_ancient ? gamma_ancient[i] : 0.0f;
                    const uint32_t ref = ref_idx[i];

                    const size_t shard = static_cast<size_t>(r) & (kFinalizeShards - 1);
                    std::lock_guard<std::mutex> lock(finalize_mutexes[shard]);

                    read_sum_sq[r] += static_cast<double>(g) * g;
                    if (g > read_top1[r]) {
                        read_top2[r] = read_top1[r];
                        read_top1[r] = g;
                    } else if (g > read_top2[r]) {
                        read_top2[r] = g;
                    }

                    const auto& bh = best_hits[r];
                    if (bh.has && ref == bh.ref_idx) {
                        em_read_results[r].gamma = g;
                        em_read_results[r].gamma_ancient = g_a;
                    }
                }
            });

        // Finalize per-read stats
        for (uint32_t r = 0; r < n_reads; ++r) {
            float keep_thr = static_cast<float>(em_min_prob);

            em_read_results[r].em_threshold = keep_thr;
            em_read_results[r].em_margin = std::max(0.0f, read_top1[r] - read_top2[r]);
            em_read_results[r].em_ref_neff = (read_sum_sq[r] > 0.0)
                ? static_cast<float>(1.0 / read_sum_sq[r]) : 0.0f;
            em_read_results[r].em_keep = (em_read_results[r].gamma >= keep_thr);
            // best_hit: we don't track EM-best ref in streaming mode, leave as false
        }

        if (verbose) {
            std::cerr << "  Streaming finalize complete\n\n";
        }

        // Post-finalize coverage accumulation: build sorted CovEntry from
        // best_hits + em_read_results.gamma, then accumulate per-ref coverage
        // stats in a single linear scan (O(1) memory per ref, no fragmentation).
        if (use_coverage_em) {
            const auto cov_start = std::chrono::steady_clock::now();
            if (verbose) {
                std::cerr << "Coverage EM (sorted in-memory pass)...\n";
                std::cerr << std::flush;
            }

            // Build compact CovEntry vector from best_hits + em_read_results
            std::vector<dart::CovEntry> cov_entries;
            cov_entries.reserve(n_reads);
            for (uint32_t r = 0; r < n_reads; ++r) {
                const auto& bh = best_hits[r];
                if (!bh.has) continue;
                const float g = em_read_results[r].gamma;
                if (g < 1e-6f) continue;
                dart::CovEntry e;
                e.ref_idx           = bh.ref_idx;
                e.gamma             = g;
                e.tstart_0          = static_cast<uint16_t>(std::min(bh.tstart_0,          65535u));
                e.aln_len_on_target = static_cast<uint16_t>(std::min(bh.aln_len_on_target, 65535u));
                e.tlen              = static_cast<uint16_t>(std::min(bh.tlen,              65535u));
                cov_entries.push_back(e);
            }
            std::sort(cov_entries.begin(), cov_entries.end(),
                [](const dart::CovEntry& a, const dart::CovEntry& b) {
                    return a.ref_idx < b.ref_idx;
                });

            dart::coverage_accumulate_from_sorted(
                cov_entries, cov_params_streaming, coverage_stats_map);

            // Apply coverage weights to em_state.weights (one reweighting step)
            double w_sum = 0.0;
            for (uint32_t j = 0; j < static_cast<uint32_t>(em_state.weights.size()); ++j) {
                auto it = coverage_stats_map.find(j);
                if (it != coverage_stats_map.end()) {
                    const double cw = static_cast<double>(it->second.coverage_weight);
                    em_state.weights[j] *= cw;
                }
                w_sum += em_state.weights[j];
            }
            if (w_sum > 0.0) {
                const double inv = 1.0 / w_sum;
                for (auto& w : em_state.weights) w *= inv;
            }

            const auto cov_end = std::chrono::steady_clock::now();
            if (verbose) {
                std::cerr << "  Coverage stats: " << coverage_stats_map.size() << " proteins\n";
                std::cerr << "  Runtime: "
                          << dart::log_utils::format_elapsed(cov_start, cov_end) << "\n\n";
            }
        }

    } else if (use_em && em_aln_count > 0) {
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

        // Build CSR-format alignment data
        auto aln_data = dart::build_alignment_data(
            em_alignments.data(), em_aln_count, n_reads, n_refs);
        // Fix 3: release em_alignments – aln_data has its own copy
        std::vector<dart::CompactAlignment>{}.swap(em_alignments);

        // Set up EM parameters
        dart::EMParams em_params;
        em_params.lambda_b = em_lambda_b;
        em_params.max_iters = em_max_iters;
        em_params.tol = em_tol;
        em_params.use_squarem = true;
        em_params.normalize_by_length = true;
        if (verbose) {
            std::cerr << "  Min iters: " << em_params.min_iters << "\n";
            std::cerr << "  SQUAREM warmup iters: " << em_params.squarem_warmup_iters << "\n";
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
        em_params.use_alignment_damage_likelihood = emi_has_alignment_damage;
        if (em_params.use_alignment_damage_likelihood) {
            em_params.use_damage = true;
        } else {
            em_params.use_damage = emi_has_damage_scores && em_damage_varies_within_read;
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

        dart::EMProgressCallback em_progress;
        if (verbose) {
            em_progress = [&](const dart::EMIterationDiagnostics& d) {
                const bool rel_valid = std::isfinite(d.rel_change) && d.iteration > 1;
                std::cerr << "    iter " << d.iteration
                          << " ll=" << std::fixed << std::setprecision(2) << d.log_likelihood
                          << " obj=" << std::fixed << std::setprecision(2) << d.objective
                          << " rel=";
                if (rel_valid) {
                    std::cerr << std::scientific << std::setprecision(3) << d.rel_change;
                } else {
                    std::cerr << "NA";
                }
                std::cerr << " dphi=" << std::scientific << std::setprecision(3)
                          << d.param_change;
                std::cerr
                          << " mode=" << (d.squarem_active ? "SQUAREM" : "EM");
                if (d.squarem_active) {
                    std::cerr << " alpha=" << std::fixed << std::setprecision(3) << d.alpha
                              << " fallback=" << (d.fallback_to_em ? "1" : "0");
                }
                std::cerr << std::defaultfloat << "\n";
            };
        }

        const auto em_start = std::chrono::steady_clock::now();

        if (use_coverage_em) {
            dart::CoverageEMParams cov_params;
            cov_params.max_outer_iters = coverage_em_iters;
            cov_params.tau = coverage_em_tau;
            cov_params.default_bins = coverage_em_bins;
            cov_params.depth_gate_lo = depth_gate_lo;
            cov_params.depth_gate_hi = depth_gate_hi;
            cov_params.w_min = w_min;
            if (verbose) {
                std::cerr << "  Coverage-EM outer iterations: " << cov_params.max_outer_iters << "\n";
            }
            em_state = dart::coverage_em_squarem(
                aln_data, em_params, cov_params, coverage_stats_map, em_progress, verbose);
        } else {
            em_state = dart::squarem_em(aln_data, em_params, nullptr, em_progress);
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
            if (use_coverage_em) {
                std::cerr << "  Coverage stats: " << coverage_stats_map.size() << " proteins\n";
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

            for (uint32_t j = start; j < end; ++j) {
                const double g = em_state.gamma[j];
                sum_sq += g * g;
                if (g > static_cast<double>(top1)) {
                    top2 = top1;
                    top1 = static_cast<float>(g);
                } else if (g > static_cast<double>(top2)) {
                    top2 = static_cast<float>(g);
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
                        if (em_params.use_damage && !em_state.gamma_ancient.empty()) {
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

            em_read_results[r] = {
                best_hit_gamma,
                best_hit_gamma_ancient,
                bitscore_best_matches_em,
                em_keep,
                keep_thr,
                std::max(0.0f, top1 - top2),
                (sum_sq > 0.0) ? static_cast<float>(1.0 / sum_sq) : 0.0f
            };
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
                uint32_t deg = aln_data.read_offsets[r + 1] - aln_data.read_offsets[r];
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

    // ── Annotation pass (runs after EM) ──────────────────────────────────────
    // Running EM first means summaries (~8 GB compact) never coexist with EM
    // state in memory, eliminating ~25-42 GB of peak RSS vs the old order.
    {
        const auto annotate_start = std::chrono::steady_clock::now();
        std::vector<std::pair<uint32_t, ProteinDamageSummary>> all_results;

        // Open per-site and corrected-protein output files before the parallel scan.
        std::ofstream sites_out_f, corrected_out_f;
        size_t written_sites = 0, written_corrected = 0;
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
        if (!corrected_file.empty()) {
            corrected_out_f.open(corrected_file);
            if (!corrected_out_f.is_open()) {
                std::cerr << "Error: Cannot open corrected file: " << corrected_file << "\n";
                return 1;
            }
        }

        if (selected_best_hits > 0) {
            reader.set_columns({dart::ColumnID::QALN, dart::ColumnID::TALN});

            reader.parallel_scan_selected(selected_row_groups, [&](
                uint32_t rg_idx, uint32_t num_rows,
                const uint32_t*, const uint32_t*,
                const float*, const float*, const float*,
                const uint16_t*, const uint16_t*,
                const float*, const float*,
                const uint16_t*, const uint16_t*,
                const uint16_t*, const uint16_t*,
                const uint16_t*, const uint16_t*,
                const uint16_t*, const uint16_t*,
                const uint16_t*, const uint16_t*)
            {
                const auto& rows = best_rows_by_rg[rg_idx];
                if (rows.empty()) return;

                // Parallel annotation within this row group.
                // string_view references into the mmap are valid until MADV_DONTNEED fires
                // after this callback returns (handled by parallel_scan_selected).
                #pragma omp parallel
                {
                    std::vector<std::pair<uint32_t, ProteinDamageSummary>> local_results;
                    std::ostringstream local_blast8;
                    std::ostringstream local_sites;
                    std::ostringstream local_corrected;
                    size_t local_corrected_count = 0;

                    #pragma omp for schedule(dynamic, 64) nowait
                    for (size_t i = 0; i < rows.size(); ++i) {
                        const auto& [row_in_group, ridx] = rows[i];
                        if (row_in_group >= num_rows) continue;
                        BestHitData& hit = best_hits[ridx];
                        if (!hit.has) continue;

                        std::string_view qaln_sv = reader.get_qaln(rg_idx, row_in_group);
                        std::string_view taln_sv = reader.get_taln(rg_idx, row_in_group);

                        // Compute aln_len_on_target (previously done in the old Pass 2).
                        uint32_t aln_len_target = 0;
                        for (char c : taln_sv) if (c != '-') ++aln_len_target;
                        hit.aln_len_on_target = aln_len_target;

                        std::string query_id(reader.read_name(ridx));
                        std::string target_id(reader.ref_name(hit.ref_idx));

                        // BLAST8 unique-mapper output.
                        if (!blast8_unique_file.empty() && need_unique_degree &&
                            ridx < read_degree.size() && read_degree[ridx] == 1) {
                            int aln_len = 0, mismatches = 0, gapopen = 0;
                            int q_aln_non_gap = 0, t_aln_non_gap = 0;
                            bool in_gap_q = false, in_gap_t = false;
                            const size_t n = std::min(qaln_sv.size(), taln_sv.size());
                            for (size_t j = 0; j < n; ++j) {
                                const char q = qaln_sv[j];
                                const char t = taln_sv[j];
                                aln_len++;
                                if (q == '-') { if (!in_gap_q) { gapopen++; in_gap_q = true; } }
                                else { in_gap_q = false; q_aln_non_gap++; }
                                if (t == '-') { if (!in_gap_t) { gapopen++; in_gap_t = true; } }
                                else { in_gap_t = false; t_aln_non_gap++; }
                                if (q != '-' && t != '-' && q != t) mismatches++;
                            }
                            const int qstart1 = static_cast<int>(hit.qstart_0 + 1);
                            const int sstart1 = static_cast<int>(hit.tstart_0 + 1);
                            const int qend1 = (q_aln_non_gap > 0) ? (qstart1 + q_aln_non_gap - 1) : qstart1;
                            const int send1 = (t_aln_non_gap > 0) ? (sstart1 + t_aln_non_gap - 1) : sstart1;
                            local_blast8 << query_id << '\t' << target_id << '\t'
                                << std::fixed << std::setprecision(3) << (hit.fident * 100.0f) << '\t'
                                << aln_len << '\t' << mismatches << '\t' << gapopen << '\t'
                                << qstart1 << '\t' << qend1 << '\t'
                                << sstart1 << '\t' << send1 << '\t'
                                << std::scientific << std::setprecision(2) << hit.evalue << '\t'
                                << std::fixed << std::setprecision(1) << hit.bits << '\n';
                        }

                        // Annotate (target_id not passed — not stored in struct).
                        auto summary = annotate_alignment(
                            query_id, qaln_sv, taln_sv,
                            hit.qstart_0, hit.tstart_0, hit.qlen,
                            hit.evalue, hit.bits, hit.fident, d_max, lambda);

                        per_read_aa_decay[ridx] = summary.aa_sum_exp_decay;

                        if (summary.total_mismatches > 0) {
                            summary.read_idx = ridx;
                            if (!corrected_file.empty()) {
                                summary.qaln = std::string(qaln_sv);
                                summary.taln = std::string(taln_sv);
                            }
                            filter_sites_by_distance(summary, max_dist);
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

                            // Write corrected proteins inline (before clearing sites).
                            if (!summary.sites.empty() && corrected_out_f.is_open() &&
                                !summary.qaln.empty()) {
                                std::string corrected = generate_corrected_protein(
                                    summary.qaln, summary.taln, summary.sites, summary.qstart_0);
                                local_corrected << ">" << query_id
                                    << " corrected_sites=" << summary.sites.size()
                                    << " damage_frac=" << std::fixed << std::setprecision(3)
                                    << summary.damage_fraction << "\n";
                                for (size_t j = 0; j < corrected.size(); j += 60)
                                    local_corrected << corrected.substr(j, 60) << "\n";
                                ++local_corrected_count;
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

                            // Return sites heap to allocator immediately (~10 GB for 40M reads).
                            summary.sites.clear();
                            summary.sites.shrink_to_fit();
                            // Return corrected-output strings (already written above).
                            std::string{}.swap(summary.qaln);
                            std::string{}.swap(summary.taln);

                            local_results.emplace_back(ridx, std::move(summary));
                        }
                    }

                    #pragma omp critical
                    {
                        all_results.insert(all_results.end(),
                            std::make_move_iterator(local_results.begin()),
                            std::make_move_iterator(local_results.end()));
                        if (!blast8_unique_file.empty()) {
                            std::string s = local_blast8.str();
                            if (!s.empty()) {
                                blast8_out << s;
                                written_unique += static_cast<size_t>(
                                    std::count(s.begin(), s.end(), '\n'));
                            }
                        }
                        if (sites_out_f.is_open()) {
                            std::string s = local_sites.str();
                            if (!s.empty()) {
                                sites_out_f << s;
                                written_sites += static_cast<size_t>(
                                    std::count(s.begin(), s.end(), '\n'));
                            }
                        }
                        if (corrected_out_f.is_open() && local_corrected_count > 0) {
                            corrected_out_f << local_corrected.str();
                            written_corrected += local_corrected_count;
                        }
                    }
                }
            });
        }

        if (!blast8_unique_file.empty() && verbose) {
            std::cerr << "Unique-mapper BLAST8 written to: " << blast8_unique_file
                      << " (" << written_unique << " rows)\n";
        }
        if (verbose && sites_out_f.is_open()) {
            std::cerr << "Site calls written to: " << sites_file
                      << " (" << written_sites << " sites)\n";
        }
        if (verbose && corrected_out_f.is_open()) {
            std::cerr << "Corrected proteins written: " << written_corrected
                      << " to " << corrected_file << "\n";
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
    for (const auto& s : summaries) {
        total_sites += s.damage_consistent;
        total_ct += s.ct_sites;
        total_ga += s.ga_sites;
        total_high += s.high_conf_sites;
        if (s.damage_consistent > 0) proteins_with_damage++;
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

    if (verbose) {
        std::cerr << "Reads with mismatches: " << summaries.size() << "\n";
        std::cerr << "Reads with damage-consistent sites: " << proteins_with_damage << "\n";
        std::cerr << "Total damage-consistent sites: " << total_sites << "\n";
        std::cerr << "  C->T class: " << total_ct << "\n";
        std::cerr << "  G->A class: " << total_ga << "\n";
        std::cerr << "  High confidence: " << total_high << "\n";
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

    dart::ModernBaseline modern_est = dart::estimate_modern_baseline(modern_proxy);
    // Use fixed pi0 = 0.10 as baseline for terminal evidence transform
    // Empirical estimation fails for heavily damaged samples where all reads have high p_read
    // (no low-tail to establish baseline), causing terminal evidence to collapse
    constexpr float PI0_FIXED = 0.10f;
    const float empirical_pi0 = PI0_FIXED;

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
    score_params.use_identity = use_identity;
    score_params.k_identity = 20.0f;
    score_params.identity_baseline = 0.90f;
    score_params.d_max = d_max;  // For fixed-average site evidence formula

    // Damage informativeness gating from AGD v3 header
    // When uninformative: terminal evidence is zeroed (posterior still computed from sites+identity)
    float damage_detectability = 1.0f;
    if (damage_index) {
        score_params.damage_informative = damage_index->damage_informative();
        score_params.channel_b_valid = damage_index->channel_b_valid();
        damage_detectability = dart::compute_damage_detectability(
            damage_index->d_max(), damage_index->damage_validated(),
            damage_index->damage_artifact(), damage_index->channel_b_valid(),
            damage_index->stop_decay_llr(), damage_index->terminal_shift());
        score_params.damage_detectability = damage_detectability;
        // All per-read lookups are done; release the 31 GB index before Bayesian scoring.
        damage_index.reset();
#ifdef __linux__
        malloc_trim(0);  // Return freed AGD pages to OS immediately
#endif
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
        std::cerr << "  Identity evidence: " << (use_identity ? "enabled (k=20, baseline=0.90)" : "disabled") << "\n";
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

    // Note: per-site and corrected-protein output are written during the annotation
    // pass (above), before sites are cleared. Nothing to do here.

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

    // Output header with Bayesian decomposition columns
    *out << "query_id\ttarget_id\tevalue\tbits\tfident\talnlen\t"
            "total_mismatches\tdamage_consistent\tnon_damage\t"
            "ct_sites\tga_sites\thigh_conf\tdamage_fraction\t"
            "max_p_damage\tp_damaged\tp_read\tposterior\tis_damaged\t"
            "classification\tlogBF_terminal\tlogBF_sites\tlogBF_identity\t"
            "m_opportunities\tk_hits\tq_eff\ttier\tdamage_class\t"
            "syn_5prime\tsyn_3prime\t"
            "delta_bits\tbpa\tlength_bin\tz_bpa\t"
            "damage_informative";
    if (use_em) {
        *out << "\tgamma\tgamma_ancient\tbest_hit\tem_keep\tem_keep_threshold\tem_margin\tem_ref_neff";
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

        // Reconstruct names at output time (O(1), no heap allocation).
        const std::string_view qname = reader.read_name(s.read_idx);
        const std::string_view tname = cached_ref_names[best_hits[s.read_idx].ref_idx];

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
             << (bayes_out.informative ? 1 : 0);

        if (use_em) {
            if (s.read_idx < em_read_results.size()) {
                const auto& em = em_read_results[s.read_idx];
                *out << "\t" << std::scientific << std::setprecision(6) << em.gamma
                     << "\t" << std::scientific << std::setprecision(6) << em.gamma_ancient
                     << "\t" << (em.best_hit ? 1 : 0)
                     << "\t" << (em.em_keep ? 1 : 0)
                     << "\t" << std::fixed << std::setprecision(6) << em.em_threshold
                     << "\t" << std::fixed << std::setprecision(4) << em.em_margin
                     << "\t" << std::fixed << std::setprecision(4) << em.em_ref_neff;
            } else {
                *out << "\t.\t.\t.\t.\t.\t.\t.";
            }
        }
        *out << "\n";
    }

    // Extract compact protein_agg inputs from summaries, then free the full summaries
    // vector before building gene_agg. Uses pre-aggregated site scalars (info_sum,
    // terminal_5, terminal_3) populated during annotation — no per-read sites heap.
    const bool need_protein_agg =
        !protein_summary_file.empty() || !protein_filtered_file.empty() || !combined_output_file.empty();
    std::vector<ProteinAggInput> protein_agg_inputs;
    if (need_protein_agg) {
        protein_agg_inputs.reserve(summaries.size());
        for (const auto& s : summaries) {
            ProteinAggInput c;
            c.read_idx = s.read_idx;
            c.ref_idx = best_hits[s.read_idx].ref_idx;
            c.tlen = static_cast<uint16_t>(std::min<size_t>(s.tlen, UINT16_MAX));
            c.p_damaged = s.p_damaged;
            c.damage_consistent = static_cast<uint32_t>(s.damage_consistent);
            c.info_sum = s.info_sum;
            c.terminal_5 = s.terminal_5;
            c.terminal_3 = s.terminal_3;
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
        !anvio_ko_file.empty() || !combined_output_file.empty() ||
        !protein_summary_file.empty() || !protein_filtered_file.empty();

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
            if (si >= 0 && static_cast<size_t>(si) < summaries.size()) {
                e.posterior        = summaries[static_cast<size_t>(si)].posterior;
                e.classification   = summaries[static_cast<size_t>(si)].classification;
            } else {
                e.posterior        = 0.0f;
                e.classification   = 0;
            }
            e.pad[0] = 0; e.pad[1] = 0;
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
    std::ofstream gf, pf_os, pff_os, cof_os, af_f;
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
        "p_protein_damaged\tcombined_score\tdelta_mle\tci_lower\tci_upper\t"
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
        "pass_mapping_filter\tfilter_fail\n";
    if (!protein_summary_file.empty()) {
        pf_os.open(protein_summary_file);
        if (!pf_os.is_open()) {
            std::cerr << "Error: Cannot open protein summary file: " << protein_summary_file << "\n";
            return 1;
        }
        pf_os << protein_header;
    }
    if (!protein_filtered_file.empty()) {
        pff_os.open(protein_filtered_file);
        if (!pff_os.is_open()) {
            std::cerr << "Error: Cannot open protein filtered file: " << protein_filtered_file << "\n";
            return 1;
        }
        pff_os << protein_header;
    }
    if (!combined_output_file.empty()) {
        cof_os.open(combined_output_file);
        if (!cof_os.is_open()) {
            std::cerr << "Error: Cannot open combined output file: " << combined_output_file << "\n";
            return 1;
        }
        cof_os << "target_id\tgroup_id\tannotation_source\tgene_present\tprotein_present\t"
                  "n_reads\teff_reads\tabundance\tbreadth\tmean_identity\tmean_posterior\t"
                  "protein_n_damage_reads\tprotein_mean_p_damaged\tprotein_p_protein_damaged\tprotein_d_aa\t"
                  "positional_score\tterminal_middle_ratio\tpass_mapping_filter\tfilter_fail\n";
    }

    // 8. Functional accumulator (written after streaming pass).
    std::unordered_map<std::string, FunctionalAccumulator> func_acc;

    // 9. Main streaming pass: one ref at a time.
    constexpr size_t TERMINAL_CODONS_S = 3;
    size_t gene_written = 0, protein_written = 0, filtered_written = 0, combined_written = 0;

    auto ri_cur = read_ref_entries.cbegin();
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
            if (it->tstart_0 < pa.min_start) pa.min_start = it->tstart_0;
            if (it->tstart_0 > pa.max_start) pa.max_start = it->tstart_0;
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
        for (auto it = pi_cur; it != pi_end; ++it) {
            if (pa.tlen == 0 && it->tlen > 0) pa.tlen = it->tlen;
            pa.n_damage_reads++;
            if (it->p_damaged > 0.5f) pa.n_damaged++;
            pa.total_damage_sites += it->damage_consistent;
            if (it->p_damaged > pa.max_p_damaged) pa.max_p_damaged = it->p_damaged;
            pa.sum_p_damaged += it->p_damaged;
            ReadDamageObs obs;
            obs.p_damaged = it->p_damaged;
            obs.log_lr    = 0.0f;
            obs.info      = it->info_sum;
            obs.n_sites   = static_cast<int>(it->damage_consistent);
            obs.is_damaged = it->p_damaged > 0.5f;
            pa.read_obs.push_back(obs);
            pa.terminal_5 += it->terminal_5;
            pa.terminal_3 += it->terminal_3;
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
                fa.add_gene(static_cast<uint32_t>(std::round(eff)), n_anc, n_mod, n_und, mean_post, gene_damaged);
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

        if ((pf_os.is_open() || pff_os.is_open() || cof_os.is_open()) && have_protein) {
            const float mean_p  = pa.n_damage_reads > 0
                ? (pa.sum_p_damaged / static_cast<float>(pa.n_damage_reads)) : 0.0f;
            const float frac    = pa.n_damage_reads > 0
                ? static_cast<float>(pa.n_damaged) / static_cast<float>(pa.n_damage_reads) : 0.0f;
            const float frac_all= pa.n_reads > 0
                ? static_cast<float>(pa.n_damaged) / static_cast<float>(pa.n_reads) : 0.0f;

            ProteinDamageResult result = fit_protein_damage(pa.read_obs);
            const float p_prot  = static_cast<float>(result.p_damaged);
            const float d_mle   = static_cast<float>(result.delta_max);
            const float ci_lo2  = static_cast<float>(result.ci_lower);
            const float ci_hi2  = static_cast<float>(result.ci_upper);
            const float lrt_p   = static_cast<float>(result.p_value);
            const float log_bf  = static_cast<float>(result.log_bayes_factor);

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
            if (pa.n_reads < min_reads)          add_fp("min_reads");
            if (pos_score < min_positional_score) add_fp("min_positional_score");
            if (term_ratio < min_terminal_ratio)  add_fp("min_terminal_ratio");
            protein_passes = filter_fail_p.empty();
            if (filter_fail_p.empty()) filter_fail_p = ".";

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
                << '\t' << std::setprecision(4) << p_prot
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
                << '\t' << (protein_passes ? 1 : 0)
                << '\t' << filter_fail_p
                << '\n';
            if (pf_os.is_open()) { pf_os << row.str(); ++protein_written; }
            if (pff_os.is_open() && protein_passes) { pff_os << row.str(); ++filtered_written; }
        }

        // --- Combined output ---
        if (cof_os.is_open()) {
            const bool gene_present    = (gene_acc.n_reads > 0);
            const bool protein_present = have_protein;

            std::string gene_fail_str;
            if (gene_present && !gene_passes) {
                auto add_gf = [&gene_fail_str](const char* t) {
                    if (!gene_fail_str.empty()) gene_fail_str.push_back(',');
                    gene_fail_str += t;
                };
                if (gene_acc.n_reads < min_reads)          add_gf("min_reads");
                if (gene_breadth < min_breadth)            add_gf("min_breadth");
                if (gene_abundance < min_depth)            add_gf("min_depth");
                if (gene_pos_score < min_positional_score) add_gf("min_positional_score");
                if (gene_term_ratio < min_terminal_ratio)  add_gf("min_terminal_ratio");
            }
            const std::string gene_filter_fail = gene_present
                ? (gene_passes ? "." : gene_fail_str)
                : "missing";
            const std::string prot_filter_fail = protein_present
                ? (filter_fail_p.empty() ? "." : filter_fail_p)
                : "missing";

            std::string comb_fail;
            auto append_cf = [&comb_fail](const std::string& t) {
                if (t.empty() || t == ".") return;
                if (!comb_fail.empty()) comb_fail.push_back(',');
                comb_fail += t;
            };
            if (!gene_present)    append_cf("no_gene");
            if (!protein_present) append_cf("no_protein");
            if (gene_filter_fail != "." && gene_filter_fail != "missing") append_cf(gene_filter_fail);
            if (prot_filter_fail != "." && prot_filter_fail != "missing") append_cf(prot_filter_fail);
            if (comb_fail.empty()) comb_fail = ".";

            const bool pass_final = gene_present && protein_present && gene_passes && protein_passes;
            const uint32_t comb_n_reads = protein_present ? static_cast<uint32_t>(pa.n_reads) : gene_acc.n_reads;
            const float comb_pos_score  = protein_present ? pa.positional_score() : gene_pos_score;
            const float comb_term_ratio = protein_present ? pa.terminal_middle_ratio() : gene_term_ratio;
            const float protein_mean_p  = pa.n_damage_reads > 0
                ? (pa.sum_p_damaged / static_cast<float>(pa.n_damage_reads)) : 0.0f;
            ProteinDamageResult comb_result = protein_present
                ? fit_protein_damage(pa.read_obs)
                : ProteinDamageResult{};
            const float protein_d_aa_c = protein_present ? [&]() {
                const float obs = use_both_ends
                    ? static_cast<float>(pa.terminal_5 + pa.terminal_3)
                    : static_cast<float>(pa.terminal_5);
                const float exp_val = static_cast<float>(pa.n_reads_eff) * n_ends_lib *
                    3.0f * 0.25f * 0.70f;
                return (exp_val > 0.0f) ? std::min(1.0f, obs / exp_val) : 0.0f;
            }() : 0.0f;

            std::string grp_id = ".";
            auto gm_it = group_map.find(target_name);
            if (gm_it != group_map.end()) grp_id = gm_it->second;

            cof_os << target_name << '\t' << grp_id << '\t' << annotation_source
                   << '\t' << (gene_present ? 1 : 0)
                   << '\t' << (protein_present ? 1 : 0)
                   << '\t' << comb_n_reads
                   << '\t' << std::fixed << std::setprecision(2) << static_cast<float>(gene_acc.eff_reads())
                   << '\t' << std::setprecision(4) << gene_abundance
                   << '\t' << std::setprecision(4) << gene_breadth
                   << '\t' << std::setprecision(4) << gene_acc.avg_fident()
                   << '\t' << std::setprecision(4) << gene_acc.avg_posterior()
                   << '\t' << static_cast<uint32_t>(pa.n_damage_reads)
                   << '\t' << std::setprecision(4) << protein_mean_p
                   << '\t' << std::setprecision(4) << static_cast<float>(comb_result.p_damaged)
                   << '\t' << std::setprecision(4) << protein_d_aa_c
                   << '\t' << std::setprecision(4) << comb_pos_score
                   << '\t' << std::setprecision(4) << comb_term_ratio
                   << '\t' << (pass_final ? 1 : 0)
                   << '\t' << comb_fail
                   << '\n';
            ++combined_written;
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
                     "n_damaged_genes\tdamaged_gene_frac\tdamage_enrichment\n";
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
                      << enrich << '\n';
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
        if (pff_os.is_open())
            std::cerr << "Filtered proteins written to: " << protein_filtered_file << " (" << filtered_written << " proteins passed mapping filters)\n";
        if (cof_os.is_open())
            std::cerr << "Combined output written to: " << combined_output_file << " (" << combined_written << " targets)\n";
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
