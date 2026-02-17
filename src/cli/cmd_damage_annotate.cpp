// agp damage-annotate: Post-mapping damage annotation
//
// Given EMI alignments (with qaln/taln), identifies amino acid positions
// where the query (observed protein) differs from the target (reference)
// in a damage-consistent way. C->T deamination (5' end) and G->A (3' end)
// cause specific amino acid substitution patterns that can be identified
// with ~90% precision when a reference protein is available.

#include "subcommand.hpp"
#include "agp/version.h"
#include "agp/damage_stats.hpp"
#include "agp/damage_index_reader.hpp"
#include "agp/bayesian_damage_score.hpp"
#include "agp/em_reassign.hpp"
#include "agp/columnar_index.hpp"
#include "agp/log_utils.hpp"
#include "agp/mmap_array.hpp"
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
#include <atomic>
#include <mutex>
#include <limits>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace agp {
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
    // AGP X-masks: (1) damage-convertible stops, (2) W/C likely from R at 5' terminus
    // High-confidence because AGP only X-masks when damage probability is high
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
    std::string query_id;
    std::string target_id;
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
    std::string query_id;
    std::string target_id;
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
    std::vector<DamageSite> sites;
    // Synonymous damage (from damage index, invisible at protein level)
    size_t syn_5prime = 0;      // Synonymous C→T at 5' terminus
    size_t syn_3prime = 0;      // Synonymous G→A at 3' terminus
    bool has_syn_data = false;  // True if damage index was used
    // Per-read p_damaged from AGP predict (stored in AGD index)
    float p_read = 0.0f;        // Per-read damage probability from predict
    bool has_p_read = false;    // True if p_read was retrieved from index
    // AA-level Bayesian evidence terms (derived from alignment + site calls)
    uint32_t aa_m_opportunities = 0;
    uint32_t aa_k_hits = 0;
    float aa_sum_qA = 0.0f;
    float aa_sum_log_qA_hits = 0.0f;
    // Full Bernoulli LLR (hits + non-hits, position-weighted)
    float aa_llr_bernoulli = 0.0f;
    float aa_sum_exp_decay = 0.0f;
    bool aa_has_bernoulli = false;
    // Stored for single-pass corrected output
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

// Parse strand and frame from AGP header suffix (e.g., readname_+_1 -> '+', 1)
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
    const std::string& query_id, const std::string& target_id,
    const std::string& qaln, const std::string& taln,
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

    // Fix 3 & 4: Full Bernoulli LLR including non-hits
    // LLR = sum_j [ y_j*log(q1_j/q0_j) + (1-y_j)*log((1-q1_j)/(1-q0_j)) ]
    float llr_bernoulli = 0.0f;
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

            // Position-weighted LLR contribution (hits only)
            // Absence evidence removed: most ancient reads have only 1-2 damage events,
            // so absence at a position doesn't mean "not ancient" - just means that
            // particular C/G wasn't converted. Including absence evidence is too aggressive.
            if (is_hit) {
                // Positive evidence: log(q1/q0)
                llr_bernoulli += std::log(q1_safe / q0_safe);
            }
            // Non-hits: no penalty (one-sided test)

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
                    query_id, target_id,
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
        p_damaged = 0.0f;  // No damage-consistent sites = no AA-level evidence
    }

    // damage_consistent now includes ALL damage-consistent sites (not just terminal)
    // terminal_damage is the subset with high p_damage
    size_t all_damage_consistent = sites.size();

    ProteinDamageSummary out{};
    out.read_idx = 0;
    out.query_id = query_id;
    out.target_id = target_id;
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
    // Full Bernoulli LLR (Fixes 3 & 4)
    out.aa_llr_bernoulli = llr_bernoulli;
    out.aa_sum_exp_decay = sum_exp_decay;
    out.aa_has_bernoulli = (aa_m_opportunities > 0);
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
    double terminal_cov = 0.0;  // Alignment overlap in terminal windows
    double middle_cov = 0.0;    // Alignment overlap in middle region
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
        const double term_mean = terminal_cov / static_cast<double>(term_len);
        if (mid_len == 0) {
            return term_mean > 0.0 ? std::numeric_limits<float>::infinity() : 0.0f;
        }
        const double mid_mean = middle_cov / static_cast<double>(mid_len);
        if (mid_mean <= 0.0) {
            return term_mean > 0.0 ? std::numeric_limits<float>::infinity() : 0.0f;
        }
        return static_cast<float>(term_mean / mid_mean);
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

    // Terminal/middle coverage ratio.
    // Low values indicate reads concentrated in the middle (typical spurious motif hits).
    float terminal_middle_ratio() const {
        if (coverage.empty()) return 0.0f;
        const size_t n = coverage.size();
        size_t edge = std::max<size_t>(1, n / 10);
        edge = std::min(edge, n / 2);
        if (edge == 0) return 0.0f;

        double term_sum = 0.0;
        double mid_sum = 0.0;
        size_t term_n = 0;
        size_t mid_n = 0;
        for (size_t i = 0; i < n; ++i) {
            if (i < edge || i >= n - edge) {
                term_sum += coverage[i];
                ++term_n;
            } else {
                mid_sum += coverage[i];
                ++mid_n;
            }
        }

        if (term_n == 0) return 0.0f;
        const double term_mean = term_sum / static_cast<double>(term_n);
        if (mid_n == 0) {
            return term_mean > 0.0 ? std::numeric_limits<float>::infinity() : 0.0f;
        }
        const double mid_mean = mid_sum / static_cast<double>(mid_n);
        if (mid_mean <= 0.0) {
            return term_mean > 0.0 ? std::numeric_limits<float>::infinity() : 0.0f;
        }
        return static_cast<float>(term_mean / mid_mean);
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
    const std::string& query_id, const std::string& qaln, const std::string& taln,
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
    bool auto_calibrate_spurious = false;  // Auto-derive thresholds from data

    // EM reassignment options (enabled by default for abundance estimation)
    bool use_em = true;
    bool em_streaming = false;       // Use streaming EM (O(num_refs) memory instead of O(num_alns))
    size_t em_max_memory_mb = 4096;  // Max memory for EM (MB), 0 = no limit, auto-switch to streaming if exceeded
    uint32_t em_max_iters = 100;
    double em_lambda_b = 3.0;
    double em_tol = 1e-4;
    double em_min_prob = 1e-6;       // absolute posterior threshold for keeping alignment
    double em_prob_fraction = 0.3;   // keep if gamma >= fraction * read_max_gamma

    // Alignment-level pre-filters (applied before any processing)
    float aln_min_identity = 0.0f;   // 0 = no filter
    float aln_min_bits = 0.0f;       // 0 = no filter
    float aln_max_evalue = 1e10f;    // very large = no filter
    int threads = 0;                 // 0 = OpenMP runtime default

    // Functional profiling options
    std::string map_file;            // gene_id -> group mapping
    std::string functional_summary_file;  // Output: per-function stats
    std::string anvio_ko_file;       // Output: Anvi'o-compatible grouped abundance
    std::string annotation_source = "AGP";

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
        } else if (strcmp(argv[i], "--em-prob-fraction") == 0 && i + 1 < argc) {
            em_prob_fraction = std::stod(argv[++i]);
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
            std::cout << "Usage: agp damage-annotate --emi <hits.emi> [options]\n\n";
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
            std::cout << "  --annotation-source STR  Source label for anvi output (default: AGP)\n\n";
            std::cout << "EM reassignment (multi-mapping resolution, enabled by default):\n";
            std::cout << "  --no-em           Disable EM (use best-hit only)\n";
            std::cout << "  --em-iters INT    Max EM iterations (default: 100)\n";
            std::cout << "  --em-lambda FLOAT Temperature parameter (default: 3.0)\n";
            std::cout << "  --em-tol FLOAT    Convergence tolerance (default: 1e-4)\n";
            std::cout << "  --em-min-prob FLOAT     Min EM posterior to keep (default: 1e-6)\n";
            std::cout << "  --em-prob-fraction FLOAT Keep if gamma >= fraction*max_gamma(read) (default: 0.3)\n";
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
            std::cerr << "Run 'agp damage-annotate --help' for usage.\n";
            return 1;
        }
    }

    if (emi_file.empty()) {
        std::cerr << "Error: No EMI index specified (--emi).\n";
        std::cerr << "Run 'agp damage-annotate --help' for usage.\n";
        return 1;
    }
    if (em_min_prob < 0.0 || em_min_prob > 1.0) {
        std::cerr << "Error: --em-min-prob must be in [0,1]\n";
        return 1;
    }
    if (em_prob_fraction < 0.0 || em_prob_fraction > 1.0) {
        std::cerr << "Error: --em-prob-fraction must be in [0,1]\n";
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
        std::cerr << "Damage annotation v" << AGP_VERSION << "\n";
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
    agp::ColumnarIndexReader reader(emi_file);
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
        size_t qstart_0 = 0;
        size_t tstart_0 = 0;
        size_t qlen = 0;
        size_t tlen = 0;
        size_t aln_len_on_target = 0;
        uint32_t best_rg = UINT32_MAX;
        uint32_t best_row = UINT32_MAX;
        std::string qaln;
        std::string taln;
    };
    std::vector<BestHitData> best_hits(n_reads);
    const bool need_unique_degree = !blast8_unique_file.empty();
    std::vector<uint32_t> read_degree;
    if (need_unique_degree) {
        read_degree.assign(n_reads, 0);
    }
    std::vector<agp::CompactAlignment> em_alignments;
    if (use_em) {
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
    if (use_em) {
        em_shard_count = static_cast<size_t>(std::max(1, omp_get_max_threads()));
    }
#endif
    std::vector<std::vector<agp::CompactAlignment>> em_alignments_by_thread;
    if (use_em) {
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
        agp::ColumnID::READ_IDX,
        agp::ColumnID::REF_IDX,
        agp::ColumnID::BIT_SCORE,
        agp::ColumnID::DAMAGE_SCORE,
        agp::ColumnID::EVALUE_LOG10,
        agp::ColumnID::DMG_K,
        agp::ColumnID::DMG_M,
        agp::ColumnID::DMG_LL_A,
        agp::ColumnID::DMG_LL_M,
        agp::ColumnID::IDENTITY_Q,
        agp::ColumnID::ALN_LEN,
        agp::ColumnID::QSTART,
        agp::ColumnID::TSTART,
        agp::ColumnID::TLEN,
        agp::ColumnID::QLEN
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
            size_t qstart_0 = 0;
            size_t tstart_0 = 0;
            size_t qlen = 0;
            size_t tlen = 0;
            size_t aln_len_on_target = 0;
            uint32_t best_rg = UINT32_MAX;
            uint32_t best_row = UINT32_MAX;
        };

        size_t local_total = 0;
        size_t local_filtered_identity = 0;
        size_t local_filtered_bits = 0;
        size_t local_filtered_evalue = 0;

        std::vector<agp::CompactAlignment> local_em;
        if (use_em) {
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
                                    size_t qstart_0,
                                    size_t tstart_0,
                                    size_t qlen_v,
                                    size_t tlen_v,
                                    size_t aln_len_on_target_v,
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
            const size_t qstart_0 = (qstart[i] > 0) ? static_cast<size_t>(qstart[i] - 1) : 0;
            const size_t tstart_0 = (tstart[i] > 0) ? static_cast<size_t>(tstart[i] - 1) : 0;
            const size_t tlen_v = static_cast<size_t>(tlen[i]);
            const size_t qlen_v = static_cast<size_t>(qlen[i]);
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
                agp::CompactAlignment ca{};
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
                ca.aln_start = static_cast<uint16_t>(std::min(tstart_0, size_t(65535)));
                ca.aln_end = static_cast<uint16_t>(
                    std::min(tstart_0 + static_cast<size_t>(aln_len_v), size_t(65535)));
                ca.identity_q = identity_q[i];
                ca.flags = tlen[i];  // pass true reference length hint to EM builder
                local_em.push_back(ca);
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
                    cand.aln_len_on_target = static_cast<size_t>(aln_len_v);
                    cand.best_rg = rg_idx;
                    cand.best_row = i;
                    local_best_sorted.emplace_back(ridx, cand);
                } else {
                    update_candidate(local_best_sorted.back().second,
                                     bits, std::clamp(damage_score[i], 0.0f, 1.0f),
                                     tidx, fident, evalue,
                                     qstart_0, tstart_0, qlen_v, tlen_v,
                                     static_cast<size_t>(aln_len_v), i);
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
                    cand.aln_len_on_target = static_cast<size_t>(aln_len_v);
                    cand.best_rg = rg_idx;
                    cand.best_row = i;
                    it->second = cand;
                } else {
                    update_candidate(it->second,
                                     bits, std::clamp(damage_score[i], 0.0f, 1.0f),
                                     tidx, fident, evalue,
                                     qstart_0, tstart_0, qlen_v, tlen_v,
                                     static_cast<size_t>(aln_len_v), i);
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

        if (use_em && !local_em.empty()) {
            size_t em_shard = 0;
#ifdef _OPENMP
            em_shard = static_cast<size_t>(omp_get_thread_num());
#endif
            if (em_shard >= em_alignments_by_thread.size()) em_shard = 0;
            auto& dst = em_alignments_by_thread[em_shard];
            dst.insert(dst.end(), local_em.begin(), local_em.end());

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
    }
    const auto pass1_end = std::chrono::steady_clock::now();
    if (verbose) {
        std::cerr << "Pass 1 runtime: "
                  << agp::log_utils::format_elapsed(pass1_start, pass1_end) << "\n";
    }

    // Pass 2: fetch alignment strings only for selected best hits.
    const auto pass2_start = std::chrono::steady_clock::now();
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

    if (selected_best_hits > 0) {
        reader.set_columns({agp::ColumnID::QALN, agp::ColumnID::TALN});
        std::atomic<size_t> loaded_best_hit_strings{0};
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
            size_t local_loaded = 0;

            for (const auto& row_pair : rows) {
                const uint32_t row_in_group = row_pair.first;
                const uint32_t ridx = row_pair.second;
                if (row_in_group >= num_rows) continue;

                std::string_view qaln_sv = reader.get_qaln(rg_idx, row_in_group);
                std::string_view taln_sv = reader.get_taln(rg_idx, row_in_group);

                BestHitData& hit = best_hits[ridx];
                hit.qaln.assign(qaln_sv.begin(), qaln_sv.end());
                hit.taln.assign(taln_sv.begin(), taln_sv.end());
                size_t aln_len_target = 0;
                for (char c : hit.taln) {
                    if (c != '-') ++aln_len_target;
                }
                hit.aln_len_on_target = aln_len_target;
                local_loaded++;
            }
            loaded_best_hit_strings.fetch_add(local_loaded, std::memory_order_relaxed);
        });

        if (verbose && loaded_best_hit_strings.load(std::memory_order_relaxed) != selected_best_hits) {
            std::cerr << "Warning: Loaded " << loaded_best_hit_strings.load(std::memory_order_relaxed)
                      << "/" << selected_best_hits << " best-hit alignment strings from EMI.\n";
        }
    }
    const auto pass2_end = std::chrono::steady_clock::now();
    if (verbose) {
        std::cerr << "Pass 2 runtime: "
                  << agp::log_utils::format_elapsed(pass2_start, pass2_end)
                  << " (" << selected_best_hits << " best-hit alignments)\n";
    }

    if (!blast8_unique_file.empty()) {
        std::ofstream bf(blast8_unique_file);
        if (!bf.is_open()) {
            std::cerr << "Error: Cannot open BLAST8 unique file: " << blast8_unique_file << "\n";
            return 1;
        }

        size_t written_unique = 0;
        for (uint32_t ridx = 0; ridx < n_reads; ++ridx) {
            if (need_unique_degree) {
                if (ridx >= read_degree.size() || read_degree[ridx] != 1) continue;
            }
            const auto& hit = best_hits[ridx];
            if (!hit.has) continue;

            const std::string qid(reader.read_name(ridx));
            const std::string sid(reader.ref_name(hit.ref_idx));

            int aln_len = 0;
            int mismatches = 0;
            int gapopen = 0;
            int q_aln_non_gap = 0;
            int t_aln_non_gap = 0;
            bool in_gap_q = false;
            bool in_gap_t = false;

            const size_t n = std::min(hit.qaln.size(), hit.taln.size());
            for (size_t i = 0; i < n; ++i) {
                const char q = hit.qaln[i];
                const char t = hit.taln[i];
                aln_len++;

                if (q == '-') {
                    if (!in_gap_q) {
                        gapopen++;
                        in_gap_q = true;
                    }
                } else {
                    in_gap_q = false;
                    q_aln_non_gap++;
                }

                if (t == '-') {
                    if (!in_gap_t) {
                        gapopen++;
                        in_gap_t = true;
                    }
                } else {
                    in_gap_t = false;
                    t_aln_non_gap++;
                }

                if (q != '-' && t != '-' && q != t) mismatches++;
            }

            const int qstart1 = static_cast<int>(hit.qstart_0 + 1);
            const int sstart1 = static_cast<int>(hit.tstart_0 + 1);
            const int qend1 = (q_aln_non_gap > 0) ? (qstart1 + q_aln_non_gap - 1) : qstart1;
            const int send1 = (t_aln_non_gap > 0) ? (sstart1 + t_aln_non_gap - 1) : sstart1;
            const float pident = hit.fident * 100.0f;

            // BLAST outfmt 6 style (qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore)
            bf << qid << '\t' << sid << '\t'
               << std::fixed << std::setprecision(3) << pident << '\t'
               << aln_len << '\t' << mismatches << '\t' << gapopen << '\t'
               << qstart1 << '\t' << qend1 << '\t'
               << sstart1 << '\t' << send1 << '\t'
               << std::scientific << std::setprecision(2) << hit.evalue << '\t'
               << std::fixed << std::setprecision(1) << hit.bits << '\n';
            written_unique++;
        }

        if (verbose) {
            std::cerr << "Unique-mapper BLAST8 written to: " << blast8_unique_file
                      << " (" << written_unique << " rows)\n";
        }
    }

    if (d_max < 0.0f) {
        if (damage_index && damage_index->d_max() > 0.0f) {
            d_max = damage_index->d_max();
            if (!lambda_set_by_user && damage_index->lambda() > 0.0f) {
                lambda = damage_index->lambda();
            }
            if (verbose) {
                std::cerr << "Using d_max from damage index header: " << std::fixed
                          << std::setprecision(3) << d_max << "\n";
                if (!lambda_set_by_user) {
                    std::cerr << "Using lambda from damage index header: " << std::fixed
                              << std::setprecision(3) << lambda << "\n";
                }
            }
        } else if (reader.d_max() > 0.0f) {
            d_max = reader.d_max();
            if (!lambda_set_by_user && reader.lambda() > 0.0f) {
                lambda = reader.lambda();
            }
            if (verbose) {
                std::cerr << "Using d_max from EMI header: " << std::fixed
                          << std::setprecision(3) << d_max << "\n";
                if (!lambda_set_by_user) {
                    std::cerr << "Using lambda from EMI header: " << std::fixed
                              << std::setprecision(3) << lambda << "\n";
                }
            }
        } else {
            if (verbose) {
                std::cerr << "Estimating d_max from EMI best hits...\n";
            }
            DamageEstimate est;
            for (uint32_t ridx = 0; ridx < n_reads; ++ridx) {
                const BestHitData& hit = best_hits[ridx];
                if (!hit.has) continue;
                std::string query_id(reader.read_name(ridx));
                count_damage_for_estimate(query_id, hit.qaln, hit.taln, hit.qstart_0, hit.qlen, est);
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

    // Process best hits into summaries
    const auto annotate_start = std::chrono::steady_clock::now();
    for (uint32_t ridx = 0; ridx < n_reads; ++ridx) {
        const auto& hit = best_hits[ridx];
        if (!hit.has) continue;

        std::string query_id(reader.read_name(ridx));
        std::string target_id(reader.ref_name(hit.ref_idx));
        auto summary = annotate_alignment(
            query_id, target_id, hit.qaln, hit.taln,
            hit.qstart_0, hit.tstart_0, hit.qlen,
            hit.evalue, hit.bits, hit.fident, d_max, lambda);

        if (summary.total_mismatches > 0) {
            summary.read_idx = ridx;
            filter_sites_by_distance(summary, max_dist);
            // Store alignment data for Bayesian scoring and corrected output
            summary.qlen = hit.qlen;
            summary.qaln = hit.qaln;
            summary.taln = hit.taln;
            summary.qstart_0 = hit.qstart_0;
            summary.tstart_0 = hit.tstart_0;
            summary.tlen = hit.tlen;
            // Short protein filtering metrics
            summary.delta_bits = hit.bits - hit.bits_second;
            summary.bpa = (summary.alnlen > 0)
                ? hit.bits / static_cast<float>(summary.alnlen) : 0.0f;
            summary.length_bin = static_cast<uint8_t>(
                agp::LengthBinStats::get_bin(hit.qlen));
            // EMI always carries DAMAGE_SCORE. Use it as per-read p_read fallback so
            // Bayesian scoring works in EMI-only workflows.
            summary.p_read = std::clamp(hit.damage_score, 0.0f, 1.0f);
            summary.has_p_read = true;
            // Look up synonymous damage and per-read p_damaged from index if available
            if (damage_index) {
                if (const AgdRecord* rec = damage_index->find(query_id)) {
                    auto syn_result = detect_synonymous_damage(
                        *rec, damage_index->d_max(), damage_index->lambda());
                    summary.syn_5prime = syn_result.synonymous_5prime;
                    summary.syn_3prime = syn_result.synonymous_3prime;
                    summary.has_syn_data = true;
                    // Retrieve per-read p_damaged (dequantize from uint8)
                    summary.p_read = static_cast<float>(rec->p_damaged_q) / 255.0f;
                    summary.has_p_read = true;
                }
            }
            summary_index_by_read[ridx] = static_cast<int32_t>(summaries.size());
            summaries.push_back(std::move(summary));
        }
    }
    const auto annotate_end = std::chrono::steady_clock::now();
    if (verbose) {
        std::cerr << "Annotation runtime: "
                  << agp::log_utils::format_elapsed(annotate_start, annotate_end)
                  << " (" << summaries.size() << " reads with mismatches)\n";
    }

    const size_t em_aln_count = em_alignments.size();

    // EM reassignment
    agp::EMState em_state;
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

    // Auto-switch to streaming EM if estimated memory exceeds limit
    // Memory estimate: ~48 bytes per alignment (gamma double + Alignment struct overhead)
    // Plus ~4 bytes per read for offsets
    if (use_em && !em_streaming && em_max_memory_mb > 0) {
        const size_t num_alignments = reader.num_alignments();
        const size_t estimated_mb = (num_alignments * 48 + n_reads * 4) / (1024 * 1024);
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

        agp::EMParams em_params;
        em_params.lambda_b = em_lambda_b;
        em_params.max_iters = em_max_iters;
        em_params.tol = em_tol;
        em_params.use_squarem = false;  // Streaming doesn't support SQUAREM yet
        em_params.use_damage = true;
        // Enable alignment damage likelihoods if available (checked via NULL in E-step)
        em_params.use_alignment_damage_likelihood = true;

        const auto em_start = std::chrono::steady_clock::now();
        // Reset column selection - previous code may have limited columns
        reader.set_all_columns();

        if (verbose) {
            std::cerr << "  Starting streaming EM...\n";
            std::cerr << std::flush;
        }

        auto streaming_result = agp::streaming_em(reader, em_params,
            [verbose](const agp::EMIterationDiagnostics& d) {
                if (verbose && (d.iteration % 10 == 0 || d.iteration < 5)) {
                    std::cerr << "    iter " << d.iteration
                              << " ll=" << std::fixed << std::setprecision(2) << d.log_likelihood
                              << "\n";
                    std::cerr << std::flush;
                }
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
                      << agp::log_utils::format_elapsed(em_start, em_end) << "\n";
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
        std::mutex finalize_mtx;

        agp::streaming_em_finalize(reader, streaming_result, em_params,
            [&](uint32_t /*rg_idx*/, uint32_t num_rows,
                const uint32_t* read_idx, const uint32_t* ref_idx,
                const float* gamma, const float* gamma_ancient)
            {
                // Collect updates locally first, then apply with lock
                struct ReadUpdate {
                    uint32_t r;
                    float g;
                    float g_a;
                    uint32_t ref;
                };
                std::vector<ReadUpdate> updates;
                updates.reserve(num_rows);

                for (uint32_t i = 0; i < num_rows; ++i) {
                    updates.push_back({
                        read_idx[i],
                        gamma[i],
                        gamma_ancient ? gamma_ancient[i] : 0.0f,
                        ref_idx[i]
                    });
                }

                // Apply updates under lock
                std::lock_guard<std::mutex> lock(finalize_mtx);
                for (const auto& u : updates) {
                    if (u.r >= n_reads) continue;  // bounds check

                    read_sum_sq[u.r] += static_cast<double>(u.g) * u.g;
                    if (u.g > read_top1[u.r]) {
                        read_top2[u.r] = read_top1[u.r];
                        read_top1[u.r] = u.g;
                    } else if (u.g > read_top2[u.r]) {
                        read_top2[u.r] = u.g;
                    }

                    // Check if this is the best-hit ref for this read
                    const auto& bh = best_hits[u.r];
                    if (bh.has && u.ref == bh.ref_idx) {
                        em_read_results[u.r].gamma = u.g;
                        em_read_results[u.r].gamma_ancient = u.g_a;
                    }
                }
            });

        // Finalize per-read stats
        for (uint32_t r = 0; r < n_reads; ++r) {
            float keep_thr = std::max(
                static_cast<float>(em_min_prob),
                static_cast<float>(em_prob_fraction) * read_top1[r]);

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
            std::cerr << "  Prob fraction keep: " << em_prob_fraction << "\n";
            std::cerr << "  Length normalization: on\n";
        }

        // Build CSR-format alignment data
        auto aln_data = agp::build_alignment_data(
            em_alignments.data(), em_aln_count, n_reads, n_refs);

        // Set up EM parameters
        agp::EMParams em_params;
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

        agp::EMProgressCallback em_progress;
        if (verbose) {
            em_progress = [&](const agp::EMIterationDiagnostics& d) {
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

        // Run EM with SQUAREM acceleration
        em_state = agp::squarem_em(aln_data, em_params, nullptr, em_progress);
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
                      << agp::log_utils::format_elapsed(em_start, em_end) << "\n";
        }

        // Extract best assignments and build per-read lookup
        auto best = agp::reassign_reads(aln_data, em_state, 0.0);

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
            float keep_thr = std::max(
                static_cast<float>(em_min_prob),
                static_cast<float>(em_prob_fraction) * top1);
            bool em_keep = true;
            if (bh.has) {
                const uint32_t best_hit_ref = bh.ref_idx;
                bitscore_best_matches_em = (best_hit_ref == em_best_ref);
                for (uint32_t j = start; j < end; ++j) {
                    if (aln_data.alignments[j].ref_idx == best_hit_ref) {
                        best_hit_gamma = static_cast<float>(em_state.gamma[j]);
                        if (em_params.use_damage && !em_state.gamma_ancient.empty()) {
                            best_hit_gamma_ancient = static_cast<float>(em_state.gamma_ancient[j]);
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
            size_t fail_abs = 0, fail_frac = 0;
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
                    else fail_frac++;
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
            std::cerr << "    failed abs(min_prob): " << fail_abs
                      << ", failed frac(max*prob_fraction): " << fail_frac << "\n";
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

    agp::ModernBaseline modern_est = agp::estimate_modern_baseline(modern_proxy);
    // Use fixed pi0 = 0.10 as baseline for terminal evidence transform
    // Empirical estimation fails for heavily damaged samples where all reads have high p_read
    // (no low-tail to establish baseline), causing terminal evidence to collapse
    constexpr float PI0_FIXED = 0.10f;
    const float empirical_pi0 = PI0_FIXED;

    // Set up scoring parameters with empirical Bayes approach
    agp::BayesianScoreParams score_params;
    score_params.pi = 0.10f;             // Prior P(ancient) - conservative
    score_params.pi0 = empirical_pi0;    // Empirical baseline for terminal evidence transform
    score_params.q_modern = modern_est.q_modern;
    score_params.w0 = 0.3f;              // Temper absence evidence (robustness)
    score_params.terminal_threshold = 0.50f;
    score_params.site_cap = 3.0f;        // Cap site contribution (matches original)
    score_params.min_opportunities = 3;  // Minimum m to trust site evidence
    score_params.channel_b_valid = (d_max > 0.0f);  // Trust site evidence if damage detected
    score_params.ancient_threshold = 0.60f;
    score_params.modern_threshold = 0.25f;
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
        damage_detectability = agp::compute_damage_detectability(
            damage_index->d_max(), damage_index->damage_validated(),
            damage_index->damage_artifact(), damage_index->channel_b_valid(),
            damage_index->stop_decay_llr(), damage_index->terminal_shift());
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

    // Compute length-binned BPA z-scores (two-pass: accumulate then normalize)
    {
        agp::LengthBinStats bpa_stats;
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

    // Write per-site output
    if (!sites_file.empty()) {
        std::ofstream sf(sites_file);
        if (!sf.is_open()) {
            std::cerr << "Error: Cannot open sites file: " << sites_file << "\n";
            return 1;
        }
        sf << "query_id\ttarget_id\tquery_pos\ttarget_pos\t"
              "target_aa\tquery_aa\tdamage_class\tconfidence\t"
              "dist_5prime\tdist_3prime\tp_damage\n";
        for (const auto& s : summaries) {
            for (const auto& site : s.sites) {
                sf << site.query_id << "\t" << site.target_id << "\t"
                   << site.query_pos << "\t" << site.target_pos << "\t"
                   << site.target_aa << "\t" << site.query_aa << "\t"
                   << (site.damage_class == 'C' ? "CT" : "GA") << "\t"
                   << (site.high_confidence ? "high" : "medium") << "\t"
                   << site.dist_5prime << "\t" << site.dist_3prime << "\t"
                   << std::fixed << std::setprecision(4) << site.p_damage << "\n";
            }
        }
        if (verbose) {
            std::cerr << "Site calls written to: " << sites_file << "\n";
        }
    }

    // Write corrected proteins (single pass — alignment data stored during parsing)
    if (!corrected_file.empty()) {
        std::ofstream cf(corrected_file);
        if (!cf.is_open()) {
            std::cerr << "Error: Cannot open corrected file: " << corrected_file << "\n";
            return 1;
        }

        size_t corrected_count = 0;
        for (const auto& s : summaries) {
            if (s.sites.empty() || s.qaln.empty()) continue;

            std::string corrected = generate_corrected_protein(
                s.qaln, s.taln, s.sites, s.qstart_0);

            cf << ">" << s.query_id << " corrected_sites=" << s.sites.size()
               << " damage_frac=" << std::fixed << std::setprecision(3)
               << s.damage_fraction << "\n";

            for (size_t j = 0; j < corrected.size(); j += 60) {
                cf << corrected.substr(j, 60) << "\n";
            }
            corrected_count++;
        }

        if (verbose) {
            std::cerr << "Corrected proteins written: " << corrected_count
                      << " to " << corrected_file << "\n";
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

    // Store per-read posteriors and classifications for gene summary aggregation
    std::vector<float> read_posteriors(summaries.size());
    std::vector<agp::AncientClassification> read_classifications(summaries.size());

    for (size_t si = 0; si < summaries.size(); ++si) {
        const auto& s = summaries[si];
        // Compute Bayesian score with decomposition
        float posterior = s.p_damaged;  // fallback when AGD p_read is unavailable
        agp::BayesianScoreOutput bayes_out{};
        bayes_out.informative = score_params.damage_informative;

        if (s.has_p_read) {
            // AA-level evidence from alignment opportunities and damage-consistent hits.
            agp::SiteEvidence ev{};
            ev.k = s.aa_k_hits;
            ev.m = s.aa_m_opportunities;
            ev.sum_qA = s.aa_sum_qA;
            ev.sum_log_qA_hits = s.aa_sum_log_qA_hits;
            ev.q_eff = (ev.m > 0) ? (ev.sum_qA / static_cast<float>(ev.m)) : 0.0f;
            // Full Bernoulli LLR (Fixes 3 & 4)
            ev.llr_bernoulli = s.aa_llr_bernoulli;
            ev.sum_exp_decay = s.aa_sum_exp_decay;
            ev.has_bernoulli = s.aa_has_bernoulli;

            // Compute Bayesian score (with optional identity evidence)
            bayes_out = agp::compute_bayesian_score(s.p_read, ev, score_params, s.fident);
            posterior = bayes_out.posterior;
        }

        *out << s.query_id << "\t" << s.target_id << "\t"
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
        auto classification = agp::classify_protein(
            bayes_out.informative, posterior, s.fident, s.z_bpa, s.delta_bits,
            score_params.ancient_threshold, score_params.modern_threshold);

        read_posteriors[si] = posterior;
        read_classifications[si] = classification;

        // Always output posterior (computed from site + identity evidence even when
        // terminal is uninformative). is_damaged classification only applies when informative.
        *out << std::fixed << std::setprecision(4) << posterior << "\t";
        std::string is_damaged = ".";
        if (bayes_out.informative) {
            is_damaged = (posterior >= threshold) ? "1" : "0";
        }
        *out << is_damaged << "\t"
             << agp::classification_name(classification) << "\t";
        *out
             << std::fixed << std::setprecision(3) << bayes_out.logBF_terminal << "\t"
             << std::fixed << std::setprecision(3) << bayes_out.logBF_sites << "\t"
             << std::fixed << std::setprecision(3) << bayes_out.logBF_identity << "\t"
             << bayes_out.m << "\t" << bayes_out.k << "\t"
             << std::fixed << std::setprecision(5) << bayes_out.q_eff << "\t"
             << agp::tier_name(bayes_out.tier) << "\t"
             << agp::damage_class_name(bayes_out.damage_class) << "\t"
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

    std::unordered_map<std::string, GeneSummaryAccumulator> gene_agg;
    bool have_gene_agg = false;

    // Build/write gene summary (supports both best-hit and EM-weighted modes)
    const bool need_gene_aggregation =
        !gene_summary_file.empty() || !functional_summary_file.empty() ||
        !anvio_ko_file.empty() || !combined_output_file.empty() ||
        !protein_summary_file.empty() || !protein_filtered_file.empty();
    if (need_gene_aggregation) {
        gene_agg.clear();

        for (uint32_t ridx = 0; ridx < n_reads; ++ridx) {
            const auto& hit = best_hits[ridx];
            if (!hit.has) continue;
            if (use_em) {
                if (ridx >= em_read_results.size() || !em_read_results[ridx].em_keep) {
                    continue;
                }
            }
            const std::string target_id(reader.ref_name(hit.ref_idx));
            auto& acc = gene_agg[target_id];
            if (acc.tlen == 0) {
                acc.target_id = target_id;
                acc.tlen = static_cast<uint32_t>(hit.tlen);
                if (acc.tlen > 0) {
                    acc.coverage.assign(acc.tlen, 0.0f);
                }
            }

            const float gamma = use_em ? em_read_results[ridx].gamma : 1.0f;
            acc.add_assignment(hit.fident, hit.tstart_0, hit.aln_len_on_target, gamma);

            const int32_t summary_idx = (ridx < summary_index_by_read.size())
                ? summary_index_by_read[ridx]
                : -1;
            if (summary_idx >= 0) {
                const size_t si = static_cast<size_t>(summary_idx);
                if (si < read_posteriors.size() && si < read_classifications.size()) {
                    acc.add_damage_observation(read_classifications[si], read_posteriors[si], gamma);
                }
            }
        }

        // Auto-calibrate spurious hit thresholds if requested
        if (auto_calibrate_spurious) {
            std::vector<float> calib_pos_scores;
            std::vector<float> calib_term_ratios;
            calib_pos_scores.reserve(gene_agg.size() / 10);
            calib_term_ratios.reserve(gene_agg.size() / 10);

            // Collect scores from well-supported genes (calibration set)
            for (const auto& [tid, acc] : gene_agg) {
                if (acc.n_reads >= 10 && acc.n_unique_starts() >= 3) {
                    calib_pos_scores.push_back(acc.positional_score());
                    calib_term_ratios.push_back(acc.terminal_middle_ratio());
                }
            }

            if (calib_pos_scores.size() >= 100) {
                // Compute 5th percentile for each metric
                std::sort(calib_pos_scores.begin(), calib_pos_scores.end());
                std::sort(calib_term_ratios.begin(), calib_term_ratios.end());

                const size_t idx_5pct = calib_pos_scores.size() / 20;  // 5th percentile
                float calib_min_pos = calib_pos_scores[idx_5pct];
                float calib_min_term = calib_term_ratios[idx_5pct];

                // Apply floor to terminal ratio (min 0.05)
                calib_min_term = std::max(0.05f, calib_min_term);

                // Override if auto-calibration yields stricter thresholds
                if (calib_min_pos > min_positional_score) {
                    min_positional_score = calib_min_pos;
                }
                if (calib_min_term > min_terminal_ratio) {
                    min_terminal_ratio = calib_min_term;
                }

                if (verbose) {
                    std::cerr << "Auto-calibration (n=" << calib_pos_scores.size() << " genes):\n"
                              << "  min_positional_score: " << std::fixed << std::setprecision(4) << min_positional_score << "\n"
                              << "  min_terminal_ratio: " << std::fixed << std::setprecision(4) << min_terminal_ratio << "\n";
                }
            } else if (verbose) {
                std::cerr << "Auto-calibration: insufficient calibration set ("
                          << calib_pos_scores.size() << " genes, need >= 100)\n";
            }
        }

        if (!gene_summary_file.empty()) {
            std::ofstream gf(gene_summary_file);
            if (!gf.is_open()) {
                std::cerr << "Error: Cannot open gene summary file: " << gene_summary_file << "\n";
                return 1;
            }

            // Full gene-level stats including coverage, damage, alignment, and diversity metrics
            gf << "target_id\tn_reads\teff_reads\tn_damage_reads\teff_damage_reads\t"
                  "n_ancient\tn_modern\tn_undetermined\t"
                  "ancient_frac\tci_low\tci_high\tmean_posterior\tmean_identity\t"
                  "breadth\tabundance\tdepth_std\tdepth_evenness\t"
                  "n_unique_starts\tstart_diversity\tstart_span\tpositional_score\t"
                  "terminal_middle_ratio\tavg_aln_len\tstd_aln_len\n";

            for (const auto& [tid, acc] : gene_agg) {
                if (acc.n_reads < min_reads) continue;
                float b = acc.breadth();
                if (b < min_breadth) continue;
                float abundance = acc.depth_mean();  // Gamma-weighted depth = abundance
                if (abundance < min_depth) continue;
                const float pos_score = acc.positional_score();
                if (pos_score < min_positional_score) continue;
                const float term_ratio = acc.terminal_middle_ratio();
                if (term_ratio < min_terminal_ratio) continue;

                // Gamma-weighted ancient/modern counts
                float n_ancient = static_cast<float>(acc.n_ancient_conf + acc.n_ancient_likely);
                float n_modern = static_cast<float>(acc.n_modern_conf);
                float n_undet = static_cast<float>(acc.n_undetermined);

                const double n_eff_class = static_cast<double>(n_ancient + n_modern);
                auto [ci_lo, ci_hi] = wilson_ci_effective(n_ancient, n_eff_class);

                gf << acc.target_id
                   << '\t' << acc.n_reads
                   << '\t' << std::fixed << std::setprecision(2) << acc.eff_reads()
                   << '\t' << acc.n_damage_reads
                   << '\t' << std::setprecision(2) << acc.n_damage_reads_eff
                   << '\t' << std::setprecision(2) << n_ancient
                   << '\t' << std::setprecision(2) << n_modern
                   << '\t' << std::setprecision(2) << n_undet
                   << '\t' << std::setprecision(4) << acc.ancient_frac()
                   << '\t' << std::setprecision(4) << ci_lo
                   << '\t' << std::setprecision(4) << ci_hi
                   << '\t' << std::setprecision(4) << acc.avg_posterior()
                   << '\t' << std::setprecision(4) << acc.avg_fident()
                   << '\t' << std::setprecision(4) << b
                   << '\t' << std::setprecision(4) << abundance
                   << '\t' << std::setprecision(4) << acc.depth_std()
                   << '\t' << std::setprecision(4) << acc.depth_evenness()
                   << '\t' << acc.n_unique_starts()
                   << '\t' << std::setprecision(4) << acc.start_diversity()
                   << '\t' << std::setprecision(4) << acc.start_span()
                   << '\t' << std::setprecision(4) << pos_score
                   << '\t' << std::setprecision(4) << term_ratio
                   << '\t' << std::setprecision(2) << acc.avg_aln_len()
                   << '\t' << std::setprecision(2) << acc.std_aln_len()
                   << '\n';
            }

            if (verbose) {
                size_t total = gene_agg.size();
                size_t passed = 0;
                for (const auto& [tid, acc] : gene_agg) {
                    if (acc.n_reads >= min_reads && acc.breadth() >= min_breadth
                        && acc.depth_mean() >= min_depth
                        && acc.positional_score() >= min_positional_score
                        && acc.terminal_middle_ratio() >= min_terminal_ratio)
                        ++passed;
                }
                std::cerr << "Gene summary written to: " << gene_summary_file
                          << " (" << passed << "/" << total << " genes after filtering)\n";
            }
        }

        have_gene_agg = true;

        // Functional profiling (aggregate gene stats by mapped group)
        bool do_functional = !map_file.empty();
        if (do_functional) {
            // Load mapping file: gene_id -> group
            auto group_map = load_mapping_file(map_file, verbose);

            // Accumulator by group
            std::unordered_map<std::string, FunctionalAccumulator> group_acc;

            // Aggregate gene stats by group
            for (const auto& [tid, acc] : gene_agg) {
                // Skip filtered genes
                float abundance = acc.depth_mean();  // EM-weighted abundance
                if (acc.n_reads < min_reads || acc.breadth() < min_breadth || abundance < min_depth
                    || acc.positional_score() < min_positional_score
                    || acc.terminal_middle_ratio() < min_terminal_ratio)
                    continue;

                // Use effective reads (gamma-weighted) for proper abundance estimation
                uint32_t eff_reads = static_cast<uint32_t>(std::round(acc.eff_reads()));
                float gene_ancient = static_cast<float>(acc.n_ancient_conf + acc.n_ancient_likely);
                float gene_modern = static_cast<float>(acc.n_modern_conf);
                float gene_undetermined = static_cast<float>(acc.n_undetermined);

                auto map_it = group_map.find(tid);
                if (map_it != group_map.end()) {
                    auto& fa = group_acc[map_it->second];
                    if (fa.function_id.empty()) {
                        fa.function_id = map_it->second;
                        fa.db_type = "MAP";
                    }
                    const float gene_mean_post = acc.avg_posterior();
                    const bool gene_is_damaged = (gene_mean_post >= threshold);
                    fa.add_gene(eff_reads, gene_ancient, gene_modern, gene_undetermined,
                                gene_mean_post, gene_is_damaged);
                }
            }

            // Write functional summary
            if (!functional_summary_file.empty()) {
                std::ofstream ff(functional_summary_file);
                if (!ff.is_open()) {
                    std::cerr << "Error: Cannot open functional summary file: " << functional_summary_file << "\n";
                } else {
                    ff << "db\tfunction_id\tlevel\tn_genes\tn_reads\tn_ancient\tn_modern\tn_undetermined\t"
                       << "ancient_frac\tci_low\tci_high\tmean_posterior\t"
                       << "n_damaged_genes\tdamaged_gene_frac\tdamage_enrichment\n";

                    uint32_t global_genes = 0;
                    uint32_t global_damaged_genes = 0;
                    for (const auto& [fid, fa] : group_acc) {
                        global_genes += fa.n_genes;
                        global_damaged_genes += fa.n_damaged_genes;
                    }
                    const float global_damaged_frac = (global_genes > 0)
                        ? static_cast<float>(global_damaged_genes) / static_cast<float>(global_genes)
                        : 0.0f;

                    auto write_acc = [&ff, global_damaged_frac](
                        const std::unordered_map<std::string, FunctionalAccumulator>& accs,
                        const std::string& level) {
                        for (const auto& [fid, fa] : accs) {
                            float frac = fa.ancient_frac();
                            auto [ci_low, ci_high] = wilson_ci_effective(
                                static_cast<double>(fa.n_ancient),
                                static_cast<double>(fa.n_ancient + fa.n_modern));
                            const float damaged_frac = fa.damaged_gene_frac();
                            const float enrichment = (global_damaged_frac > 0.0f)
                                ? (damaged_frac / global_damaged_frac)
                                : 0.0f;
                            ff << fa.db_type << '\t' << fa.function_id << '\t' << level << '\t'
                               << fa.n_genes << '\t' << fa.n_reads << '\t'
                               << std::fixed << std::setprecision(1) << fa.n_ancient << '\t'
                               << fa.n_modern << '\t' << fa.n_undetermined << '\t'
                               << std::setprecision(4) << frac << '\t' << ci_low << '\t' << ci_high << '\t'
                               << fa.mean_posterior() << '\t'
                               << fa.n_damaged_genes << '\t'
                               << damaged_frac << '\t'
                               << enrichment << '\n';
                        }
                    };

                    write_acc(group_acc, "group");

                    if (verbose) {
                        size_t total_funcs = group_acc.size();
                        std::cerr << "Functional summary written to: " << functional_summary_file
                                  << " (" << total_funcs << " functions)\n";
                    }
                }
            }

            // Write Anvi'o-compatible gene abundance (for anvi-estimate-metabolism --enzymes-txt)
            // Format: gene_id  enzyme_accession  source  coverage  detection
            // Where coverage = EM-weighted abundance (gamma-weighted depth), detection = breadth
            if (!anvio_ko_file.empty() && !group_map.empty()) {
                std::ofstream af(anvio_ko_file);
                if (!af.is_open()) {
                    std::cerr << "Error: Cannot open Anvi'o KO file: " << anvio_ko_file << "\n";
                } else {
                    af << "gene_id\tenzyme_accession\tsource\tcoverage\tdetection\n";
                    size_t n_written = 0;
                    for (const auto& [tid, acc] : gene_agg) {
                        // Skip filtered genes
                        float abundance = acc.depth_mean();  // EM-weighted abundance
                        if (acc.n_reads < min_reads || acc.breadth() < min_breadth || abundance < min_depth
                            || acc.positional_score() < min_positional_score
                            || acc.terminal_middle_ratio() < min_terminal_ratio)
                            continue;
                        // Look up mapped group for this gene
                        auto it = group_map.find(tid);
                        if (it == group_map.end()) continue;
                        // coverage = EM-weighted abundance (sum of gamma-weighted depth)
                        af << tid << '\t' << it->second << '\t' << annotation_source << '\t'
                           << std::fixed << std::setprecision(4) << abundance << '\t'
                           << std::setprecision(4) << acc.breadth() << '\n';
                        ++n_written;
                    }
                    if (verbose) {
                        std::cerr << "Anvi'o gene abundance written to: " << anvio_ko_file
                                  << " (" << n_written << " genes with mapping annotations)\n";
                    }
                }
            }
        }
    }

    std::unordered_map<std::string, ProteinAggregate> protein_agg;
    bool have_protein_agg = false;

    // Build/write per-protein aggregated summary (aggregates across all reads per target)
    if (!protein_summary_file.empty() || !protein_filtered_file.empty() || !combined_output_file.empty()) {
        protein_agg.clear();

        // Assignment channel: aggregate all EM-kept best-hit reads for abundance/patterns.
        // Damage channel: aggregate mismatch-evidence summaries for damage statistics.
        constexpr size_t TERMINAL_CODONS = 3;  // First 3 codons = 9 nt

        for (uint32_t ridx = 0; ridx < n_reads; ++ridx) {
            const auto& hit = best_hits[ridx];
            if (!hit.has) continue;
            if (use_em) {
                if (ridx >= em_read_results.size() || !em_read_results[ridx].em_keep) {
                    continue;
                }
            }

            const std::string target_id(reader.ref_name(hit.ref_idx));
            auto& agg = protein_agg[target_id];
            if (agg.protein_id.empty()) agg.protein_id = target_id;
            if (agg.tlen == 0) agg.tlen = static_cast<uint32_t>(hit.tlen);
            const float gamma = use_em ? em_read_results[ridx].gamma : 1.0f;
            agg.n_reads++;
            agg.n_reads_eff += gamma;
            agg.sum_delta_bits += gamma * (hit.bits - hit.bits_second);
            const float bpa = (hit.aln_len_on_target > 0)
                ? (hit.bits / static_cast<float>(hit.aln_len_on_target))
                : 0.0f;
            agg.sum_bpa += gamma * bpa;
            const uint8_t length_bin = static_cast<uint8_t>(agp::LengthBinStats::get_bin(hit.qlen));
            agg.length_bin_counts[length_bin]++;
            agg.unique_starts.insert(static_cast<uint32_t>(hit.tstart_0));
            if (hit.tstart_0 < agg.min_start) agg.min_start = static_cast<uint32_t>(hit.tstart_0);
            if (hit.tstart_0 > agg.max_start) agg.max_start = static_cast<uint32_t>(hit.tstart_0);

            // Mapping pattern metrics (lightweight, no per-position vector needed)
            const size_t aln_on_target = hit.aln_len_on_target;
            agg.sum_aln_len += static_cast<double>(gamma) * static_cast<double>(aln_on_target);
            agg.sum_aln_len_sq += static_cast<double>(gamma) * static_cast<double>(aln_on_target) *
                                  static_cast<double>(aln_on_target);
            const size_t tlen = std::max<size_t>(1, static_cast<size_t>(hit.tlen));
            const size_t window = std::min(std::max<size_t>(1, tlen / 10), tlen / 2);
            const size_t start = std::min(static_cast<size_t>(hit.tstart_0), tlen);
            const size_t end = std::min(start + aln_on_target, tlen);
            if (end > start) {
                const size_t left_end = window;
                const size_t right_start = (tlen > window) ? (tlen - window) : 0;
                const size_t left_overlap = (start < left_end) ? (std::min(end, left_end) - start) : 0;
                const size_t right_overlap = (end > right_start && start < tlen)
                    ? (end - std::max(start, right_start)) : 0;
                size_t middle_overlap = 0;
                if (right_start > left_end) {
                    const size_t mid_start = std::max(start, left_end);
                    const size_t mid_end = std::min(end, right_start);
                    if (mid_end > mid_start) middle_overlap = mid_end - mid_start;
                }
                agg.terminal_cov += static_cast<double>(gamma) *
                                    static_cast<double>(left_overlap + right_overlap);
                agg.middle_cov += static_cast<double>(gamma) * static_cast<double>(middle_overlap);
            }
        }

        for (const auto& s : summaries) {
            if (use_em) {
                if (s.read_idx >= em_read_results.size() || !em_read_results[s.read_idx].em_keep) {
                    continue;
                }
            }
            auto& agg = protein_agg[s.target_id];
            if (agg.protein_id.empty()) agg.protein_id = s.target_id;
            if (agg.tlen == 0) agg.tlen = static_cast<uint32_t>(s.tlen);
            agg.n_damage_reads++;
            if (s.p_damaged > 0.5f) agg.n_damaged++;
            agg.total_damage_sites += s.damage_consistent;
            if (s.p_damaged > agg.max_p_damaged) agg.max_p_damaged = s.p_damaged;
            agg.sum_p_damaged += s.p_damaged;

            // Collect per-read observation for Beta-Binomial model
            // info = sum of p_damage at sites (informativeness for weighting)
            float info = 0.0f;
            for (const auto& site : s.sites) {
                info += site.p_damage;
            }
            ReadDamageObs read_obs;
            read_obs.p_damaged = s.p_damaged;
            read_obs.log_lr = 0.0f;  // Could compute from damage_fraction vs expected
            read_obs.info = info;
            read_obs.n_sites = static_cast<int>(s.damage_consistent);
            read_obs.is_damaged = s.p_damaged > 0.5f;
            agg.read_obs.push_back(read_obs);

            // Count terminal sites from per-read sites
            for (const auto& site : s.sites) {
                size_t codon_pos_5 = site.dist_5prime / 3;
                size_t codon_pos_3 = site.dist_3prime / 3;
                if (site.damage_class == 'C' && codon_pos_5 < TERMINAL_CODONS) {
                    agg.terminal_5++;
                } else if (site.damage_class == 'G' && codon_pos_3 < TERMINAL_CODONS) {
                    agg.terminal_3++;
                }
            }
        }

        // Auto-detect library type from terminal site ratio
        size_t total_5 = 0, total_3 = 0;
        for (const auto& [id, agg] : protein_agg) {
            total_5 += agg.terminal_5;
            total_3 += agg.terminal_3;
        }
        bool use_both_ends = true;  // default: dsDNA
        if (lib_type == "ss") {
            use_both_ends = false;
        } else if (lib_type == "ds") {
            use_both_ends = true;
        } else {
            // auto-detect: if 3'/5' ratio < 0.3, assume ssDNA
            float ratio = (total_5 > 0) ? static_cast<float>(total_3) / total_5 : 0.0f;
            use_both_ends = (ratio > 0.3f);
            if (verbose) {
                std::cerr << "Library type auto-detected: " << (use_both_ends ? "dsDNA" : "ssDNA")
                          << " (3'/5' ratio: " << ratio << ")\n";
            }
        }

        std::ofstream pf;
        if (!protein_summary_file.empty()) {
            pf.open(protein_summary_file);
            if (!pf.is_open()) {
                std::cerr << "Error: Cannot open protein summary file: " << protein_summary_file << "\n";
                return 1;
            }
        }
        std::ofstream pff;
        if (!protein_filtered_file.empty()) {
            pff.open(protein_filtered_file);
            if (!pff.is_open()) {
                std::cerr << "Error: Cannot open protein filtered file: " << protein_filtered_file << "\n";
                return 1;
            }
        }

        // D_aa calculation constants (derived from observed AA susceptibility):
        // P_SUSCEPTIBLE: ~25% of codons can produce visible damage substitutions
        // P_NONSYNONYMOUS: ~70% of those produce amino acid changes
        // These are not magic numbers - derived from codon table analysis
        constexpr float P_SUSCEPTIBLE = 0.25f;
        constexpr float P_NONSYNONYMOUS = 0.70f;

        // Extended header with proper statistical outputs
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
            "pass_mapping_filter\tfilter_fail\n";
        if (pf.is_open()) pf << protein_header;
        if (pff.is_open()) pff << protein_header;

        size_t filtered_written = 0;
        for (const auto& [id, agg] : protein_agg) {
            float mean_p = agg.n_damage_reads > 0
                ? (agg.sum_p_damaged / static_cast<float>(agg.n_damage_reads))
                : 0.0f;
            float frac = agg.n_damage_reads > 0
                ? static_cast<float>(agg.n_damaged) / static_cast<float>(agg.n_damage_reads)
                : 0.0f;
            float frac_all = agg.n_reads > 0
                ? static_cast<float>(agg.n_damaged) / static_cast<float>(agg.n_reads)
                : 0.0f;
            uint32_t assign_n_reads = 0;
            float assign_eff_reads = 0.0f;
            uint32_t assign_n_damage_reads = 0;
            float assign_eff_damage_reads = 0.0f;
            float assign_n_ancient = 0.0f;
            float assign_n_modern = 0.0f;
            float assign_n_undetermined = 0.0f;
            float assign_ancient_frac = 0.0f;
            float assign_ci_low = 0.0f;
            float assign_ci_high = 0.0f;
            float assign_mean_posterior = 0.0f;
            float assign_mean_identity = 0.0f;
            float assign_breadth = 0.0f;
            float assign_abundance = 0.0f;
            float assign_depth_std = 0.0f;
            float assign_depth_evenness = 0.0f;

            auto git = gene_agg.find(id);
            if (git != gene_agg.end()) {
                const auto& g = git->second;
                assign_n_reads = g.n_reads;
                assign_eff_reads = g.eff_reads();
                assign_n_damage_reads = g.n_damage_reads;
                assign_eff_damage_reads = static_cast<float>(g.n_damage_reads_eff);
                assign_n_ancient = static_cast<float>(g.n_ancient_conf + g.n_ancient_likely);
                assign_n_modern = static_cast<float>(g.n_modern_conf);
                assign_n_undetermined = static_cast<float>(g.n_undetermined);
                assign_ancient_frac = g.ancient_frac();
                const double n_eff_class = static_cast<double>(assign_n_ancient + assign_n_modern);
                auto [ci_lo, ci_hi] = wilson_ci_effective(assign_n_ancient, n_eff_class);
                assign_ci_low = ci_lo;
                assign_ci_high = ci_hi;
                assign_mean_posterior = g.avg_posterior();
                assign_mean_identity = g.avg_fident();
                assign_breadth = g.breadth();
                assign_abundance = g.depth_mean();
                assign_depth_std = g.depth_std();
                assign_depth_evenness = g.depth_evenness();
            }

            // ================================================================
            // PROPER STATISTICAL INFERENCE: Beta-Binomial model
            // Replaces naive "any read damaged → p=1.0" with principled approach
            // ================================================================
            ProteinDamageResult result = fit_protein_damage(agg.read_obs);

            float p_protein_damaged = static_cast<float>(result.p_damaged);
            float delta_mle = static_cast<float>(result.delta_max);
            float ci_lower = static_cast<float>(result.ci_lower);
            float ci_upper = static_cast<float>(result.ci_upper);
            float lrt_pvalue = static_cast<float>(result.p_value);
            float log_bf = static_cast<float>(result.log_bayes_factor);
            float combined_score = p_protein_damaged;

            // D_aa: AA-level damage amplitude (0-1)
            // D_aa = observed_terminal_sites / expected_if_fully_damaged
            // This is an empirical estimate complementing the Bayesian model
            float observed = use_both_ends
                ? static_cast<float>(agg.terminal_5 + agg.terminal_3)
                : static_cast<float>(agg.terminal_5);
            float n_ends = use_both_ends ? 2.0f : 1.0f;
            float expected = static_cast<float>(agg.n_reads_eff) * n_ends *
                            TERMINAL_CODONS * P_SUSCEPTIBLE * P_NONSYNONYMOUS;
            float d_aa = (expected > 0.0f) ? std::min(1.0f, observed / expected) : 0.0f;

            const float pos_score = agg.positional_score();
            const float start_div = agg.start_diversity();
            const float term_ratio = agg.terminal_middle_ratio();
            const bool pass_mapping_filter =
                (agg.n_reads >= min_reads) &&
                (pos_score >= min_positional_score) &&
                (term_ratio >= min_terminal_ratio);
            std::string filter_fail = ".";
            if (!pass_mapping_filter) {
                filter_fail.clear();
                auto add_fail = [&filter_fail](const char* token) {
                    if (!filter_fail.empty()) filter_fail.push_back(',');
                    filter_fail += token;
                };
                if (agg.n_reads < min_reads) add_fail("min_reads");
                if (pos_score < min_positional_score) add_fail("min_positional_score");
                if (term_ratio < min_terminal_ratio) add_fail("min_terminal_ratio");
                if (filter_fail.empty()) filter_fail = "unknown";
            }

            std::ostringstream row;
            row << agg.protein_id << "\t"
                << agg.n_reads << "\t"
                << std::fixed << std::setprecision(2) << agg.n_reads_eff << "\t"
                << agg.n_damage_reads << "\t"
                << agg.n_damaged << "\t"
                << agg.total_damage_sites << "\t"
                << agg.terminal_5 << "\t"
                << agg.terminal_3 << "\t"
                << std::fixed << std::setprecision(4) << agg.max_p_damaged << "\t"
                << std::fixed << std::setprecision(4) << mean_p << "\t"
                << std::fixed << std::setprecision(4) << frac << "\t"
                << std::fixed << std::setprecision(4) << frac_all << "\t"
                << std::fixed << std::setprecision(4) << p_protein_damaged << "\t"
                << std::fixed << std::setprecision(4) << combined_score << "\t"
                << std::fixed << std::setprecision(4) << delta_mle << "\t"
                << std::fixed << std::setprecision(4) << ci_lower << "\t"
                << std::fixed << std::setprecision(4) << ci_upper << "\t"
                << std::scientific << std::setprecision(3) << lrt_pvalue << "\t"
                << std::fixed << std::setprecision(2) << log_bf << "\t"
                << std::fixed << std::setprecision(4) << d_aa << "\t"
                << std::fixed << std::setprecision(1)
                << (agg.n_reads_eff > 0.0 ? agg.sum_delta_bits / static_cast<float>(agg.n_reads_eff) : 0.0f) << "\t"
                << std::fixed << std::setprecision(3)
                << (agg.n_reads_eff > 0.0 ? agg.sum_bpa / static_cast<float>(agg.n_reads_eff) : 0.0f) << "\t"
                << static_cast<int>(std::max_element(
                       agg.length_bin_counts.begin(),
                       agg.length_bin_counts.end()) - agg.length_bin_counts.begin()) << "\t"
                << agg.unique_starts.size() << "\t"
                << std::fixed << std::setprecision(4) << start_div << "\t"
                << std::fixed << std::setprecision(4) << agg.start_span() << "\t"
                << std::fixed << std::setprecision(4) << pos_score << "\t"
                << std::fixed << std::setprecision(4) << term_ratio << "\t"
                << std::fixed << std::setprecision(2) << agg.avg_aln_len() << "\t"
                << std::fixed << std::setprecision(2) << agg.std_aln_len() << "\t"
                << assign_n_reads << "\t"
                << std::fixed << std::setprecision(2) << assign_eff_reads << "\t"
                << assign_n_damage_reads << "\t"
                << std::fixed << std::setprecision(2) << assign_eff_damage_reads << "\t"
                << std::fixed << std::setprecision(2) << assign_n_ancient << "\t"
                << std::fixed << std::setprecision(2) << assign_n_modern << "\t"
                << std::fixed << std::setprecision(2) << assign_n_undetermined << "\t"
                << std::fixed << std::setprecision(4) << assign_ancient_frac << "\t"
                << std::fixed << std::setprecision(4) << assign_ci_low << "\t"
                << std::fixed << std::setprecision(4) << assign_ci_high << "\t"
                << std::fixed << std::setprecision(4) << assign_mean_posterior << "\t"
                << std::fixed << std::setprecision(4) << assign_mean_identity << "\t"
                << std::fixed << std::setprecision(4) << assign_breadth << "\t"
                << std::fixed << std::setprecision(4) << assign_abundance << "\t"
                << std::fixed << std::setprecision(4) << assign_depth_std << "\t"
                << std::fixed << std::setprecision(4) << assign_depth_evenness << "\t"
                << (pass_mapping_filter ? 1 : 0) << "\t"
                << filter_fail << "\n";

            if (pf.is_open()) pf << row.str();
            if (pff.is_open() && pass_mapping_filter) {
                pff << row.str();
                ++filtered_written;
            }
        }

        if (verbose) {
            if (pf.is_open()) {
                std::cerr << "Protein summary written to: " << protein_summary_file
                          << " (" << protein_agg.size() << " proteins)\n";
            }
            if (pff.is_open()) {
                std::cerr << "Filtered proteins written to: " << protein_filtered_file
                          << " (" << filtered_written << " proteins passed mapping filters)\n";
            }
        }

        have_protein_agg = true;
    }

    if (!combined_output_file.empty()) {
        std::ofstream cof(combined_output_file);
        if (!cof.is_open()) {
            std::cerr << "Error: Cannot open combined output file: " << combined_output_file << "\n";
            return 1;
        }

        std::unordered_map<std::string, std::string> combined_group_map;
        if (!map_file.empty()) {
            combined_group_map = load_mapping_file(map_file, false);
        }

        bool use_both_ends_combined = true;
        if (have_protein_agg) {
            size_t total_5 = 0;
            size_t total_3 = 0;
            for (const auto& kv : protein_agg) {
                total_5 += kv.second.terminal_5;
                total_3 += kv.second.terminal_3;
            }
            if (lib_type == "ss") {
                use_both_ends_combined = false;
            } else if (lib_type == "ds") {
                use_both_ends_combined = true;
            } else {
                float ratio = (total_5 > 0) ? static_cast<float>(total_3) / total_5 : 0.0f;
                use_both_ends_combined = (ratio > 0.3f);
            }
        }

        cof << "target_id\tgroup_id\tannotation_source\tgene_present\tprotein_present\t"
               "n_reads\teff_reads\tabundance\tbreadth\tmean_identity\tmean_posterior\t"
               "protein_n_damage_reads\tprotein_mean_p_damaged\tprotein_p_protein_damaged\tprotein_d_aa\t"
               "positional_score\tterminal_middle_ratio\tpass_mapping_filter\tfilter_fail\n";

        std::vector<std::string> target_ids;
        target_ids.reserve(gene_agg.size() + protein_agg.size());
        std::unordered_set<std::string> seen_ids;
        seen_ids.reserve(gene_agg.size() + protein_agg.size());

        if (have_protein_agg) {
            for (const auto& kv : protein_agg) {
                if (seen_ids.insert(kv.first).second) target_ids.push_back(kv.first);
            }
        }
        if (have_gene_agg) {
            for (const auto& kv : gene_agg) {
                if (seen_ids.insert(kv.first).second) target_ids.push_back(kv.first);
            }
        }

        size_t combined_written = 0;
        for (const auto& tid : target_ids) {
            const auto git = gene_agg.find(tid);
            const auto pit = protein_agg.find(tid);
            const bool gene_present = (git != gene_agg.end());
            const bool protein_present = (pit != protein_agg.end());

            uint32_t gene_n_reads = 0;
            float gene_eff_reads = 0.0f;
            float gene_n_ancient = 0.0f;
            float gene_n_modern = 0.0f;
            float gene_n_undetermined = 0.0f;
            float gene_ancient_frac = 0.0f;
            float gene_ci_low = 0.0f;
            float gene_ci_high = 0.0f;
            float gene_mean_posterior = 0.0f;
            float gene_mean_identity = 0.0f;
            float gene_breadth = 0.0f;
            float gene_abundance = 0.0f;
            float gene_depth_std = 0.0f;
            float gene_depth_evenness = 0.0f;
            uint32_t gene_n_unique_starts = 0;
            float gene_start_diversity = 0.0f;
            float gene_start_span = 0.0f;
            float gene_pos_score = 0.0f;
            float gene_term_ratio = 0.0f;
            float gene_avg_aln_len = 0.0f;
            float gene_std_aln_len = 0.0f;
            bool pass_gene_filter = false;
            std::string gene_filter_fail = "missing";
            if (gene_present) {
                const auto& g = git->second;
                gene_n_reads = g.n_reads;
                gene_eff_reads = g.eff_reads();
                gene_n_ancient = static_cast<float>(g.n_ancient_conf + g.n_ancient_likely);
                gene_n_modern = static_cast<float>(g.n_modern_conf);
                gene_n_undetermined = static_cast<float>(g.n_undetermined);
                gene_ancient_frac = g.ancient_frac();
                const double n_eff_class = static_cast<double>(gene_n_ancient + gene_n_modern);
                auto ci = wilson_ci_effective(gene_n_ancient, n_eff_class);
                gene_ci_low = ci.first;
                gene_ci_high = ci.second;
                gene_mean_posterior = g.avg_posterior();
                gene_mean_identity = g.avg_fident();
                gene_breadth = g.breadth();
                gene_abundance = g.depth_mean();
                gene_depth_std = g.depth_std();
                gene_depth_evenness = g.depth_evenness();
                gene_n_unique_starts = g.n_unique_starts();
                gene_start_diversity = g.start_diversity();
                gene_start_span = g.start_span();
                gene_pos_score = g.positional_score();
                gene_term_ratio = g.terminal_middle_ratio();
                gene_avg_aln_len = g.avg_aln_len();
                gene_std_aln_len = g.std_aln_len();

                std::string fail;
                auto add_fail = [&fail](const char* token) {
                    if (!fail.empty()) fail.push_back(',');
                    fail += token;
                };
                if (g.n_reads < min_reads) add_fail("min_reads");
                if (gene_breadth < min_breadth) add_fail("min_breadth");
                if (gene_abundance < min_depth) add_fail("min_depth");
                if (gene_pos_score < min_positional_score) add_fail("min_positional_score");
                if (gene_term_ratio < min_terminal_ratio) add_fail("min_terminal_ratio");

                pass_gene_filter =
                    fail.empty();
                gene_filter_fail = fail.empty() ? "." : fail;
            }

            uint32_t protein_n_reads = 0;
            uint32_t protein_n_damage_reads = 0;
            uint32_t protein_n_damaged = 0;
            uint32_t protein_damage_sites = 0;
            uint32_t protein_terminal_5 = 0;
            uint32_t protein_terminal_3 = 0;
            float protein_mean_p = 0.0f;
            float protein_max_p = 0.0f;
            float protein_frac_damaged = 0.0f;
            float protein_p_protein_damaged = 0.0f;
            float protein_delta_mle = 0.0f;
            float protein_ci_lower = 0.0f;
            float protein_ci_upper = 0.0f;
            float protein_lrt_pvalue = 1.0f;
            float protein_log_bf = 0.0f;
            float protein_d_aa = 0.0f;
            float protein_avg_delta_bits = 0.0f;
            float protein_avg_bpa = 0.0f;
            int protein_length_bin = 0;
            uint32_t protein_n_unique_starts = 0;
            float protein_start_diversity = 0.0f;
            float protein_start_span = 0.0f;
            float protein_pos_score = 0.0f;
            float protein_term_ratio = 0.0f;
            float protein_avg_aln_len = 0.0f;
            float protein_std_aln_len = 0.0f;
            bool pass_mapping_filter = false;
            std::string protein_filter_fail = "missing";
            if (protein_present) {
                const auto& p = pit->second;
                protein_n_reads = static_cast<uint32_t>(p.n_reads);
                protein_n_damage_reads = static_cast<uint32_t>(p.n_damage_reads);
                protein_n_damaged = static_cast<uint32_t>(p.n_damaged);
                protein_damage_sites = static_cast<uint32_t>(p.total_damage_sites);
                protein_terminal_5 = static_cast<uint32_t>(p.terminal_5);
                protein_terminal_3 = static_cast<uint32_t>(p.terminal_3);
                protein_mean_p = (p.n_damage_reads > 0)
                    ? (p.sum_p_damaged / static_cast<float>(p.n_damage_reads))
                    : 0.0f;
                protein_max_p = p.max_p_damaged;
                protein_frac_damaged = (p.n_damage_reads > 0)
                    ? static_cast<float>(p.n_damaged) / static_cast<float>(p.n_damage_reads)
                    : 0.0f;
                ProteinDamageResult result = fit_protein_damage(p.read_obs);
                protein_p_protein_damaged = static_cast<float>(result.p_damaged);
                protein_delta_mle = static_cast<float>(result.delta_max);
                protein_ci_lower = static_cast<float>(result.ci_lower);
                protein_ci_upper = static_cast<float>(result.ci_upper);
                protein_lrt_pvalue = static_cast<float>(result.p_value);
                protein_log_bf = static_cast<float>(result.log_bayes_factor);
                constexpr float TERMINAL_CODONS_COMBINED = 3.0f;
                constexpr float P_SUSCEPTIBLE_COMBINED = 0.25f;
                constexpr float P_NONSYNONYMOUS_COMBINED = 0.70f;
                const float observed = use_both_ends_combined
                    ? static_cast<float>(p.terminal_5 + p.terminal_3)
                    : static_cast<float>(p.terminal_5);
                const float n_ends = use_both_ends_combined ? 2.0f : 1.0f;
                const float expected = static_cast<float>(p.n_reads_eff) * n_ends *
                    TERMINAL_CODONS_COMBINED * P_SUSCEPTIBLE_COMBINED * P_NONSYNONYMOUS_COMBINED;
                protein_d_aa = (expected > 0.0f) ? std::min(1.0f, observed / expected) : 0.0f;
                protein_avg_delta_bits = (p.n_reads_eff > 0.0)
                    ? (p.sum_delta_bits / static_cast<float>(p.n_reads_eff)) : 0.0f;
                protein_avg_bpa = (p.n_reads_eff > 0.0)
                    ? (p.sum_bpa / static_cast<float>(p.n_reads_eff)) : 0.0f;
                protein_length_bin = static_cast<int>(std::max_element(
                    p.length_bin_counts.begin(), p.length_bin_counts.end()) - p.length_bin_counts.begin());
                protein_n_unique_starts = static_cast<uint32_t>(p.unique_starts.size());
                protein_start_diversity = p.start_diversity();
                protein_start_span = p.start_span();
                protein_pos_score = p.positional_score();
                protein_term_ratio = p.terminal_middle_ratio();
                protein_avg_aln_len = p.avg_aln_len();
                protein_std_aln_len = p.std_aln_len();

                std::string fail;
                auto add_fail = [&fail](const char* token) {
                    if (!fail.empty()) fail.push_back(',');
                    fail += token;
                };
                if (p.n_reads < min_reads) add_fail("min_reads");
                if (protein_pos_score < min_positional_score) add_fail("min_positional_score");
                if (protein_term_ratio < min_terminal_ratio) add_fail("min_terminal_ratio");
                pass_mapping_filter = fail.empty();
                protein_filter_fail = fail.empty() ? "." : fail;
            }

            std::string group_id = ".";
            auto gm_it = combined_group_map.find(tid);
            if (gm_it != combined_group_map.end()) {
                group_id = gm_it->second;
            }
            std::string filter_fail;
            auto append_fail = [&filter_fail](const std::string& token) {
                if (token.empty()) return;
                if (!filter_fail.empty()) filter_fail.push_back(',');
                filter_fail += token;
            };
            if (!gene_present) append_fail("no_gene");
            if (!protein_present) append_fail("no_protein");
            if (gene_filter_fail != "." && gene_filter_fail != "missing") append_fail(gene_filter_fail);
            if (protein_filter_fail != "." && protein_filter_fail != "missing") append_fail(protein_filter_fail);
            if (filter_fail.empty()) filter_fail = ".";
            const bool pass_mapping_filter_final =
                gene_present && protein_present && pass_gene_filter && pass_mapping_filter;
            const uint32_t combined_n_reads = protein_present ? protein_n_reads : gene_n_reads;
            const float combined_pos_score = protein_present ? protein_pos_score : gene_pos_score;
            const float combined_term_ratio = protein_present ? protein_term_ratio : gene_term_ratio;

            cof << tid << '\t' << group_id << '\t' << annotation_source << '\t'
                << (gene_present ? 1 : 0) << '\t'
                << (protein_present ? 1 : 0) << '\t'
                << combined_n_reads << '\t'
                << std::fixed << std::setprecision(2) << gene_eff_reads << '\t'
                << std::setprecision(4) << gene_abundance << '\t'
                << std::setprecision(4) << gene_breadth << '\t'
                << std::setprecision(4) << gene_mean_identity << '\t'
                << std::setprecision(4) << gene_mean_posterior << '\t'
                << protein_n_damage_reads << '\t'
                << std::setprecision(4) << protein_mean_p << '\t'
                << std::setprecision(4) << protein_p_protein_damaged << '\t'
                << std::setprecision(4) << protein_d_aa << '\t'
                << std::setprecision(4) << combined_pos_score << '\t'
                << std::setprecision(4) << combined_term_ratio << '\t'
                << (pass_mapping_filter_final ? 1 : 0) << '\t'
                << filter_fail << '\n';
            combined_written++;
        }

        if (verbose) {
            std::cerr << "Combined output written to: " << combined_output_file
                      << " (" << combined_written << " targets)\n";
        }
    }

    if (verbose) {
        size_t assigned_kept_reads = 0;
        for (uint32_t ridx = 0; ridx < n_reads; ++ridx) {
            if (!best_hits[ridx].has) continue;
            if (use_em) {
                if (ridx >= em_read_results.size() || !em_read_results[ridx].em_keep) continue;
            }
            assigned_kept_reads++;
        }
        std::cerr << "\nSummary: assigned_reads=" << assigned_kept_reads
                  << ", mismatch_reads=" << summaries.size()
                  << ", reads_with_damage_sites=" << proteins_with_damage << "\n";
        const auto run_end = std::chrono::steady_clock::now();
        std::cerr << "Runtime: "
                  << agp::log_utils::format_elapsed(run_start, run_end) << "\n";
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
}  // namespace agp
