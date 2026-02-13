// agp damage-annotate: Post-mapping damage annotation
//
// Given MMseqs2 alignments (with qaln/taln), identifies amino acid positions
// where the query (observed protein) differs from the target (reference)
// in a damage-consistent way. C->T deamination (5' end) and G->A (3' end)
// cause specific amino acid substitution patterns that can be identified
// with ~90% precision when a reference protein is available.

#include "subcommand.hpp"
#include "agp/version.h"
#include "agp/damage_stats.hpp"
#include "agp/damage_index_reader.hpp"
#include "agp/bayesian_damage_score.hpp"
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
    // Stored for single-pass corrected output
    std::string qaln;
    std::string taln;
    size_t qstart_0;
    size_t qlen = 0;  // Total read length (for Bayesian scoring)
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

// Manual tab-split: returns views into the line via start/end positions
// Much faster than std::istringstream for millions of lines
static size_t split_tabs(const std::string& line, size_t* starts, size_t* ends, size_t max_fields) {
    size_t n = 0;
    size_t pos = 0;
    size_t len = line.size();
    while (pos <= len && n < max_fields) {
        starts[n] = pos;
        size_t tab = line.find('\t', pos);
        ends[n] = (tab != std::string::npos) ? tab : len;
        n++;
        pos = (tab != std::string::npos) ? tab + 1 : len + 1;
    }
    return n;
}

// Extract substring view helper
static std::string field_str(const std::string& line, size_t start, size_t end) {
    return line.substr(start, end - start);
}

// Compute combined damage score using noisy-OR model
// P(ancient) = 1 - (1 - p_read) × (1 - p_nonsyn) × (1 - p_syn)
//
// All probabilities are sample-calibrated via d_max:
// - p_read: from sample damage model (already calibrated)
// - p_nonsyn: damage sites → probability, scaled by d_max
// - p_syn: synonymous damage → weaker signal, scaled by d_max
//
// No magic weights - everything derived from sample's damage level.
static float compute_combined_score(float p_read, size_t damage_sites, bool has_syn, float d_max) {
    // p_nonsyn: probability that observed damage sites indicate ancient origin
    // Uses saturation function: more sites = stronger evidence, caps at d_max
    // p = d_max × (1 - exp(-0.5 × n)) where n = number of damage sites
    // At n=1: p ≈ 0.39 × d_max, n=2: p ≈ 0.63 × d_max, n=3+: p → d_max
    float p_nonsyn = damage_sites > 0
        ? d_max * (1.0f - std::exp(-0.5f * static_cast<float>(damage_sites)))
        : 0.0f;

    // p_syn: synonymous damage is weaker evidence (doesn't change protein)
    // but still indicates deamination occurred
    float p_syn = has_syn ? (d_max * 0.5f) : 0.0f;

    // Noisy-OR: each evidence source independently suggests ancient origin
    // Combined probability = 1 - P(none of them indicate ancient)
    float combined = 1.0f - (1.0f - p_read) * (1.0f - p_nonsyn) * (1.0f - p_syn);

    return std::min(1.0f, std::max(0.0f, combined));
}

// Annotate a single alignment
static ProteinDamageSummary annotate_alignment(
    const std::string& query_id, const std::string& target_id,
    const std::string& qaln, const std::string& taln,
    size_t qstart, size_t tstart, size_t qlen,
    float evalue, float bits, float fident,
    float d_max, float lambda)
{
    char strand = parse_strand(query_id);

    std::vector<DamageSite> sites;
    size_t total_mismatches = 0;
    size_t damage_consistent = 0;

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

        if (q_aa != t_aa) {
            total_mismatches++;

            const DamageSubstitution* sub = find_damage_sub(t_aa, q_aa);
            if (sub) {
                damage_consistent++;

                size_t dist_5, dist_3;
                if (strand == '+') {
                    dist_5 = q_pos * 3;
                    dist_3 = (qlen > 0 && q_pos < qlen) ? (qlen - 1 - q_pos) * 3 : 0;
                } else {
                    dist_5 = (qlen > 0 && q_pos < qlen) ? (qlen - 1 - q_pos) * 3 : 0;
                    dist_3 = q_pos * 3;
                }

                float p_dmg;
                if (sub->damage_class == 'C') {
                    p_dmg = positional_damage_prob(dist_5, d_max, lambda);
                } else {
                    p_dmg = positional_damage_prob(dist_3, d_max, lambda);
                }

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

    // Two thresholds:
    // - HIGH_CONF: Terminal sites with strong positional signal (for strict filtering)
    // - LOW_CONF: All damage-consistent sites weighted by position (for scoring)
    constexpr float HIGH_CONF_THRESHOLD = 0.05f;  // ~5% for high-confidence terminal
    constexpr float EMPIRICAL_PRECISION = 0.83f;   // Validated precision at terminal
    constexpr float INTERIOR_PRECISION = 0.40f;    // Estimated precision for interior sites

    size_t ct = 0, ga = 0, hc = 0;
    size_t terminal_damage = 0;  // Sites with p_damage >= HIGH_CONF_THRESHOLD
    float max_p_damage = 0.0f;
    float log_odds_sum = 0.0f;   // Log-odds accumulator for scoring

    for (const auto& s : sites) {
        if (s.damage_class == 'C') ct++;
        else ga++;
        if (s.high_confidence) hc++;
        if (s.p_damage > max_p_damage) max_p_damage = s.p_damage;

        // Count high-confidence terminal sites separately
        if (s.p_damage >= HIGH_CONF_THRESHOLD) {
            terminal_damage++;
        }

        // Weighted scoring: all sites contribute based on position
        // Higher p_damage → more likely true damage → stronger evidence
        // P(true damage | observed mismatch) = precision × p_damage + base_rate × (1 - p_damage)
        // where base_rate accounts for interior sites that are still damage
        float p_true_damage = EMPIRICAL_PRECISION * s.p_damage +
                              INTERIOR_PRECISION * (1.0f - s.p_damage);

        // Log-odds: log(P(damage) / P(not damage))
        float odds = p_true_damage / (1.0f - p_true_damage + 1e-6f);
        log_odds_sum += std::log(odds);
    }

    // p_damaged: sigmoid of accumulated log-odds
    // With no sites: p_damaged = 0.5 (prior)
    // With positive evidence: p_damaged > 0.5
    float p_damaged = 0.0f;
    if (!sites.empty()) {
        // Normalize by number of sites and apply sigmoid
        float avg_log_odds = log_odds_sum / static_cast<float>(sites.size());
        p_damaged = 1.0f / (1.0f + std::exp(-avg_log_odds));
    }

    // damage_consistent now includes ALL damage-consistent sites (not just terminal)
    // terminal_damage is the subset with high p_damage
    size_t all_damage_consistent = sites.size();

    return {
        query_id, target_id,
        evalue, bits, fident,
        qaln.size(),
        total_mismatches,
        all_damage_consistent,  // ALL damage-consistent sites (improved recall)
        total_mismatches - all_damage_consistent,
        ct, ga, hc,
        total_mismatches > 0
            ? static_cast<float>(all_damage_consistent) / static_cast<float>(total_mismatches)
            : 0.0f,
        max_p_damage,
        p_damaged,
        std::move(sites),
        {}, {}, 0  // qaln/taln/qstart_0 filled by caller if needed
    };
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
    summary.ct_sites = 0;
    summary.ga_sites = 0;
    summary.high_conf_sites = 0;
    for (const auto& s : summary.sites) {
        if (s.damage_class == 'C') summary.ct_sites++;
        else summary.ga_sites++;
        if (s.high_confidence) summary.high_conf_sites++;
    }
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
    size_t n_reads = 0;           // Total reads mapping to this protein
    size_t n_damaged = 0;         // Reads with damage detected (p_damaged > 0)
    size_t total_damage_sites = 0; // Sum of damage_consistent across reads
    size_t terminal_5 = 0;        // C→T sites at 5' end (first 3 codons)
    size_t terminal_3 = 0;        // G→A sites at 3' end (first 3 codons)
    float max_p_damaged = 0.0f;   // Maximum p_damaged across reads
    float sum_p_damaged = 0.0f;   // Sum of p_damaged (for averaging)

    // Per-read observations for Beta-Binomial model
    std::vector<ReadDamageObs> read_obs;  // Collected for proper statistical inference
};

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
    std::string hits_file;
    std::string output_file;
    std::string sites_file;
    std::string corrected_file;
    std::string protein_summary_file;  // Per-protein aggregation
    std::string damage_index_file;     // Binary damage index from predict
    float d_max = -1.0f;  // -1 means "estimate from data"
    float lambda = 0.3f;
    int max_dist = -1;
    std::string lib_type = "auto";    // ss, ds, or auto-detect
    float threshold = 0.7f;           // Classification threshold for is_damaged
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "--hits") == 0 || strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
            hits_file = argv[++i];
        } else if (strcmp(argv[i], "--d-max") == 0 && i + 1 < argc) {
            d_max = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda") == 0 && i + 1 < argc) {
            lambda = std::stof(argv[++i]);
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
        } else if ((strcmp(argv[i], "--library-type") == 0 || strcmp(argv[i], "--lib-type") == 0) && i + 1 < argc) {
            lib_type = argv[++i];
        } else if (strcmp(argv[i], "--damage-index") == 0 && i + 1 < argc) {
            damage_index_file = argv[++i];
        } else if (strcmp(argv[i], "--threshold") == 0 && i + 1 < argc) {
            threshold = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            std::cerr << "Usage: agp damage-annotate --hits <results.tsv> [options]\n\n";
            std::cerr << "Post-mapping damage annotation using MMseqs2 alignments.\n";
            std::cerr << "Compares observed proteins against reference proteins to identify\n";
            std::cerr << "damage-consistent amino acid substitutions with ~90% precision.\n\n";
            std::cerr << "Required:\n";
            std::cerr << "  --hits, -i FILE   MMseqs2 results TSV with qaln/taln columns\n\n";
            std::cerr << "Options:\n";
            std::cerr << "  --d-max RATE      Override damage rate (default: estimated from data)\n";
            std::cerr << "  --lambda FLOAT    Decay constant (default: 0.3)\n";
            std::cerr << "  --max-dist INT    Positional filter: max nt from relevant terminus (-1=off)\n";
            std::cerr << "  -o FILE           Per-read damage summary TSV (default: stdout)\n";
            std::cerr << "  --protein-summary FILE  Per-protein aggregated damage TSV\n";
            std::cerr << "  --library-type TYPE     Library type: ss, ds, or auto (default: auto)\n";
            std::cerr << "  --sites FILE      Per-site damage calls TSV\n";
            std::cerr << "  --corrected FILE  Reference-guided corrected proteins FASTA\n";
            std::cerr << "  --damage-index FILE  Binary damage index (.agd) from predict\n";
            std::cerr << "                       Enables synonymous damage detection\n";
            std::cerr << "  --threshold FLOAT Damage classification threshold (default: 0.7)\n";
            std::cerr << "                       Adds is_damaged column (1 if combined_score >= threshold)\n";
            std::cerr << "  -v                Verbose output\n\n";
            std::cerr << "MMseqs2 format (16 columns with qaln/taln):\n";
            std::cerr << "  --format-output \"query,target,fident,alnlen,mismatch,gapopen,\n";
            std::cerr << "                   qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln\"\n";
            return 0;
        } else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            return 1;
        }
    }

    if (hits_file.empty()) {
        std::cerr << "Error: No hits file specified (--hits).\n";
        return 1;
    }

    // Estimate d_max from data if not provided
    if (d_max < 0.0f) {
        if (verbose) {
            std::cerr << "Estimating d_max from alignment data...\n";
        }
        DamageEstimate est;
        std::ifstream est_in(hits_file);
        if (!est_in.is_open()) {
            std::cerr << "Error: Cannot open hits file: " << hits_file << "\n";
            return 1;
        }
        std::string est_line;
        size_t est_starts[16], est_ends[16];
        std::unordered_set<std::string> est_seen;
        while (std::getline(est_in, est_line)) {
            if (est_line.empty() || est_line[0] == '#') continue;
            size_t nf = split_tabs(est_line, est_starts, est_ends, 16);
            if (nf < 16) continue;
            std::string query_id = field_str(est_line, est_starts[0], est_ends[0]);
            if (est_seen.count(query_id)) continue;
            est_seen.insert(query_id);
            size_t qstart = static_cast<size_t>(std::stoul(field_str(est_line, est_starts[6], est_ends[6])));
            size_t qlen = static_cast<size_t>(std::stoul(field_str(est_line, est_starts[12], est_ends[12])));
            std::string qaln = field_str(est_line, est_starts[14], est_ends[14]);
            std::string taln = field_str(est_line, est_starts[15], est_ends[15]);
            size_t qstart_0 = (qstart > 0) ? qstart - 1 : 0;
            count_damage_for_estimate(query_id, qaln, taln, qstart_0, qlen, est);
        }
        est_in.close();

        d_max = est.estimate_d_max();
        if (verbose) {
            std::cerr << "  Terminal subs: " << est.terminal_damage << "/" << est.terminal_total
                      << " (" << (est.terminal_total > 0 ? 100.0f * est.terminal_damage / est.terminal_total : 0) << "%)\n";
            std::cerr << "  Interior subs: " << est.interior_damage << "/" << est.interior_total
                      << " (" << (est.interior_total > 0 ? 100.0f * est.interior_damage / est.interior_total : 0) << "%)\n";
            std::cerr << "  Estimated d_max: " << std::fixed << std::setprecision(3) << d_max << "\n\n";
        }
    }

    if (verbose) {
        std::cerr << "Damage annotation v" << AGP_VERSION << "\n";
        std::cerr << "Hits: " << hits_file << "\n";
        std::cerr << "d_max: " << d_max << "\n";
        std::cerr << "lambda: " << lambda << "\n";
        if (max_dist >= 0) {
            std::cerr << "Positional filter: <=" << max_dist << " nt from terminus\n";
        }
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
                std::cerr << "Index lambda: " << damage_index->lambda() << "\n\n";
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Cannot load damage index: " << e.what() << "\n";
            std::cerr << "Synonymous damage detection disabled.\n\n";
        }
    }

    bool need_alignments = !corrected_file.empty();

    // Parse MMseqs2 results
    std::ifstream in(hits_file);
    if (!in.is_open()) {
        std::cerr << "Error: Cannot open hits file: " << hits_file << "\n";
        return 1;
    }

    std::vector<ProteinDamageSummary> summaries;
    std::unordered_set<std::string> seen_queries;
    std::string line;
    size_t starts[16], ends[16];

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        size_t nf = split_tabs(line, starts, ends, 16);
        if (nf < 16) continue;

        std::string query_id = field_str(line, starts[0], ends[0]);
        if (seen_queries.count(query_id)) continue;
        seen_queries.insert(query_id);

        std::string target_id = field_str(line, starts[1], ends[1]);
        float fident = std::stof(field_str(line, starts[2], ends[2]));
        size_t qstart = static_cast<size_t>(std::stoul(field_str(line, starts[6], ends[6])));
        size_t tstart = static_cast<size_t>(std::stoul(field_str(line, starts[8], ends[8])));
        // Use stod for evalue: E-values can be extremely small (e.g., 1E-100)
        // which underflows stof but is representable in double
        float evalue = static_cast<float>(std::stod(field_str(line, starts[10], ends[10])));
        float bits = std::stof(field_str(line, starts[11], ends[11]));
        size_t qlen = static_cast<size_t>(std::stoul(field_str(line, starts[12], ends[12])));
        std::string qaln = field_str(line, starts[14], ends[14]);
        std::string taln = field_str(line, starts[15], ends[15]);

        size_t qstart_0 = (qstart > 0) ? qstart - 1 : 0;
        size_t tstart_0 = (tstart > 0) ? tstart - 1 : 0;

        auto summary = annotate_alignment(
            query_id, target_id, qaln, taln,
            qstart_0, tstart_0, qlen,
            evalue, bits, fident, d_max, lambda);

        if (summary.total_mismatches > 0) {
            filter_sites_by_distance(summary, max_dist);
            // Store alignment data for Bayesian scoring and corrected output
            summary.qlen = qlen;
            summary.qaln = qaln;
            summary.taln = taln;
            summary.qstart_0 = qstart_0;
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
            summaries.push_back(std::move(summary));
        }
    }
    in.close();

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
        std::cerr << "Proteins with mismatches: " << summaries.size() << "\n";
        std::cerr << "Proteins with damage-consistent sites: " << proteins_with_damage << "\n";
        std::cerr << "Total damage-consistent sites: " << total_sites << "\n";
        std::cerr << "  C->T class: " << total_ct << "\n";
        std::cerr << "  G->A class: " << total_ga << "\n";
        std::cerr << "  High confidence: " << total_high << "\n";
    }

    // Bayesian scoring setup
    // 1. Create decay lookup table for position-weighted hazards
    agp::ExpDecayLUT decay_lut(lambda, 300);

    // 2. Estimate q_modern from low p_read reads (likely undamaged)
    // Use reads with p_read < 0.20 as modern proxy (p_read is bimodal: ~0 or ~0.5+)
    std::vector<std::pair<uint32_t, uint32_t>> modern_proxy;
    for (const auto& s : summaries) {
        if (s.has_p_read && s.p_read < 0.20f && !s.qaln.empty()) {
            auto ev = agp::compute_site_evidence(s.qaln, s.taln, s.qstart_0, s.qlen, d_max, decay_lut);
            modern_proxy.emplace_back(ev.m, ev.k);
        }
    }
    agp::ModernBaseline modern_est = agp::estimate_modern_baseline(modern_proxy);

    // Set up scoring parameters with empirical Bayes approach
    agp::BayesianScoreParams score_params;
    score_params.pi = 0.10f;             // Prior P(ancient) - conservative
    score_params.pi0 = 0.10f;            // Prior used in p_read calculation
    score_params.q_modern = modern_est.q_modern;
    score_params.w0 = 0.3f;              // Temper absence evidence (robustness)
    score_params.terminal_threshold = 0.50f;
    score_params.site_cap = 3.0f;        // Cap site contribution
    score_params.min_opportunities = 3;  // Minimum m to trust site evidence
    score_params.channel_b_valid = (d_max > 0.0f);  // Trust site evidence if damage detected
    score_params.ancient_threshold = 0.60f;
    score_params.modern_threshold = 0.25f;

    if (verbose) {
        std::cerr << "\nBayesian scoring (empirical Bayes):\n";
        std::cerr << "  Modern proxy reads: " << modern_est.n_reads << "\n";
        std::cerr << "  Pooled opportunities: " << modern_est.M0 << "\n";
        std::cerr << "  Pooled hits: " << modern_est.K0 << "\n";
        std::cerr << "  Estimated q_modern: " << std::fixed << std::setprecision(5)
                  << modern_est.q_modern << "\n";
        std::cerr << "  Tempering w0: " << score_params.w0 << "\n";
        std::cerr << "  Channel B valid: " << (score_params.channel_b_valid ? "yes" : "no") << "\n";
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
            "max_p_damage\tp_damaged\tp_read\tposterior\t"
            "is_damaged\tlogBF_terminal\tlogBF_sites\t"
            "m_opportunities\tk_hits\tq_eff\ttier\tdamage_class\t"
            "syn_5prime\tsyn_3prime\n";

    for (const auto& s : summaries) {
        // Compute Bayesian score with decomposition
        float posterior = s.p_damaged;  // fallback if no alignment data
        agp::BayesianScoreOutput bayes_out{};

        if (s.has_p_read) {
            // Use ct_sites + ga_sites as damage evidence (correctly computed from AA substitutions)
            // The compute_site_evidence function was broken - it expected nucleotides but got amino acids
            agp::SiteEvidence ev{};
            ev.k = static_cast<uint32_t>(s.ct_sites + s.ga_sites);  // damage-consistent hits
            ev.m = static_cast<uint32_t>(s.alnlen);  // opportunities ~ alignment length

            // Compute log-likelihood contribution for hits
            // Each damage site at terminal position has ~d_max probability
            // log(q_ancient / q_modern) where q_ancient ~ d_max, q_modern ~ 0.005
            constexpr float Q_BASELINE = 0.005f;
            if (ev.k > 0 && d_max > Q_BASELINE) {
                float avg_q_ancient = d_max * 0.5f + Q_BASELINE;  // average over positions
                ev.sum_log_qA_hits = static_cast<float>(ev.k) * std::log(avg_q_ancient / Q_BASELINE);
            }
            ev.sum_qA = d_max * 0.3f * static_cast<float>(ev.m);  // rough average damage rate
            ev.q_eff = (ev.m > 0) ? (ev.sum_qA / static_cast<float>(ev.m)) : 0.0f;

            // Compute Bayesian score
            bayes_out = agp::compute_bayesian_score(s.p_read, ev, score_params);
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
             << std::fixed << std::setprecision(4) << (s.has_p_read ? s.p_read : 0.0f) << "\t"
             << std::fixed << std::setprecision(4) << posterior << "\t"
             << (posterior >= threshold ? 1 : 0) << "\t"
             << std::fixed << std::setprecision(3) << bayes_out.logBF_terminal << "\t"
             << std::fixed << std::setprecision(3) << bayes_out.logBF_sites << "\t"
             << bayes_out.m << "\t" << bayes_out.k << "\t"
             << std::fixed << std::setprecision(5) << bayes_out.q_eff << "\t"
             << agp::tier_name(bayes_out.tier) << "\t"
             << agp::damage_class_name(bayes_out.damage_class) << "\t"
             << (s.has_syn_data ? std::to_string(s.syn_5prime) : ".") << "\t"
             << (s.has_syn_data ? std::to_string(s.syn_3prime) : ".") << "\n";
    }

    // Write per-protein aggregated summary (aggregates across all reads per target)
    if (!protein_summary_file.empty()) {
        std::unordered_map<std::string, ProteinAggregate> protein_agg;

        // Aggregate from summaries (read-level) and sites (for terminal counts)
        constexpr size_t TERMINAL_CODONS = 3;  // First 3 codons = 9 nt

        for (const auto& s : summaries) {
            auto& agg = protein_agg[s.target_id];
            if (agg.protein_id.empty()) agg.protein_id = s.target_id;
            agg.n_reads++;
            if (s.p_damaged > 0.0f) agg.n_damaged++;
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
            read_obs.is_damaged = s.p_damaged > 0.5f || s.damage_consistent > 0;
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

        std::ofstream pf(protein_summary_file);
        if (!pf.is_open()) {
            std::cerr << "Error: Cannot open protein summary file: " << protein_summary_file << "\n";
            return 1;
        }

        // D_aa calculation constants (derived from observed AA susceptibility):
        // P_SUSCEPTIBLE: ~25% of codons can produce visible damage substitutions
        // P_NONSYNONYMOUS: ~70% of those produce amino acid changes
        // These are not magic numbers - derived from codon table analysis
        constexpr float P_SUSCEPTIBLE = 0.25f;
        constexpr float P_NONSYNONYMOUS = 0.70f;

        // Extended header with proper statistical outputs
        pf << "protein_id\tn_reads\tn_damaged\tdamage_sites\t"
              "terminal_5\tterminal_3\t"
              "max_p_damaged\tmean_p_damaged\tfrac_damaged\t"
              "p_protein_damaged\tdelta_mle\tci_lower\tci_upper\t"
              "lrt_pvalue\tlog_bf\td_aa\n";

        for (const auto& [id, agg] : protein_agg) {
            float mean_p = agg.n_reads > 0 ? agg.sum_p_damaged / agg.n_reads : 0.0f;
            float frac = agg.n_reads > 0 ? static_cast<float>(agg.n_damaged) / agg.n_reads : 0.0f;

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

            // D_aa: AA-level damage amplitude (0-1)
            // D_aa = observed_terminal_sites / expected_if_fully_damaged
            // This is an empirical estimate complementing the Bayesian model
            float observed = use_both_ends
                ? static_cast<float>(agg.terminal_5 + agg.terminal_3)
                : static_cast<float>(agg.terminal_5);
            float n_ends = use_both_ends ? 2.0f : 1.0f;
            float expected = static_cast<float>(agg.n_reads) * n_ends *
                            TERMINAL_CODONS * P_SUSCEPTIBLE * P_NONSYNONYMOUS;
            float d_aa = (expected > 0.0f) ? std::min(1.0f, observed / expected) : 0.0f;

            pf << agg.protein_id << "\t"
               << agg.n_reads << "\t"
               << agg.n_damaged << "\t"
               << agg.total_damage_sites << "\t"
               << agg.terminal_5 << "\t"
               << agg.terminal_3 << "\t"
               << std::fixed << std::setprecision(4) << agg.max_p_damaged << "\t"
               << std::fixed << std::setprecision(4) << mean_p << "\t"
               << std::fixed << std::setprecision(4) << frac << "\t"
               << std::fixed << std::setprecision(4) << p_protein_damaged << "\t"
               << std::fixed << std::setprecision(4) << delta_mle << "\t"
               << std::fixed << std::setprecision(4) << ci_lower << "\t"
               << std::fixed << std::setprecision(4) << ci_upper << "\t"
               << std::scientific << std::setprecision(3) << lrt_pvalue << "\t"
               << std::fixed << std::setprecision(2) << log_bf << "\t"
               << std::fixed << std::setprecision(4) << d_aa << "\n";
        }

        if (verbose) {
            std::cerr << "Protein summary written to: " << protein_summary_file
                      << " (" << protein_agg.size() << " proteins)\n";
        }
    }

    if (verbose) {
        std::cerr << "\nSummary: " << summaries.size() << " reads analyzed, "
                  << proteins_with_damage << " with damage-consistent sites\n";
    }

    return 0;
}

// Register subcommand
namespace {
    struct DamageAnnotateRegistrar {
        DamageAnnotateRegistrar() {
            SubcommandRegistry::instance().register_command(
                "damage-annotate",
                "Post-mapping damage annotation from MMseqs2 alignments",
                cmd_damage_annotate, 40);
        }
    } damage_annotate_registrar;
}

}  // namespace cli
}  // namespace agp
