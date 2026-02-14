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
#include "agp/em_reassign.hpp"
#include "agp/mmap_array.hpp"
#include "agp/reference_stats.hpp"
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

    // Short protein filtering aggregates
    float sum_delta_bits = 0.0f;
    float sum_bpa = 0.0f;
    std::array<uint32_t, 4> length_bin_counts{};  // per-bin read counts
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

// Per-gene accumulator for functional profiling
// Aggregates read-level damage annotations into gene-level summaries
struct GeneSummaryAccumulator {
    std::string target_id;
    uint32_t tlen = 0;
    uint32_t n_reads = 0;
    uint32_t n_ancient_conf = 0;
    uint32_t n_ancient_likely = 0;
    uint32_t n_undetermined = 0;
    uint32_t n_modern_conf = 0;
    double sum_posterior = 0.0;
    double sum_fident = 0.0;
    std::vector<float> coverage;  // Per-position coverage depth

    void add_read(AncientClassification cls, float posterior, float fident,
                  size_t tstart, size_t aln_len_on_target) {
        ++n_reads;
        sum_posterior += posterior;
        sum_fident += fident;

        switch (cls) {
            case AncientClassification::AncientConfident: ++n_ancient_conf; break;
            case AncientClassification::AncientLikely:    ++n_ancient_likely; break;
            case AncientClassification::Undetermined:     ++n_undetermined; break;
            case AncientClassification::ModernConfident:  ++n_modern_conf; break;
        }

        // Update per-position coverage
        if (tlen > 0 && coverage.empty()) {
            coverage.resize(tlen, 0.0f);
        }
        size_t end = std::min(tstart + aln_len_on_target, static_cast<size_t>(tlen));
        for (size_t p = tstart; p < end; ++p) {
            coverage[p] += 1.0f;
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
        return static_cast<float>(total / coverage.size());
    }

    float avg_posterior() const {
        return n_reads > 0 ? static_cast<float>(sum_posterior / n_reads) : 0.0f;
    }

    float avg_fident() const {
        return n_reads > 0 ? static_cast<float>(sum_fident / n_reads) : 0.0f;
    }
};

// Functional profiling accumulator
struct FunctionalAccumulator {
    std::string function_id;
    std::string db_type;  // "KEGG", "CAZyme", "Viral"
    uint32_t n_genes = 0;
    uint32_t n_reads = 0;
    float n_ancient = 0.0f;  // Can be fractional from EM
    float n_modern = 0.0f;
    float n_undetermined = 0.0f;
    double sum_posterior = 0.0;

    void add_gene(uint32_t gene_reads, float gene_ancient, float gene_modern,
                  float gene_undetermined, float gene_mean_posterior) {
        ++n_genes;
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

// Extract CAZyme class from family (e.g., "GH13" -> "GH", "GT2" -> "GT")
static std::string cazyme_class_from_family(const std::string& family) {
    std::string cls;
    for (char c : family) {
        if (std::isalpha(c)) {
            cls += c;
        } else {
            break;  // Stop at first non-alpha (digit)
        }
    }
    return cls;
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
    bool use_identity = true;         // Identity evidence enabled by default (+7.6% AUC)
    bool verbose = false;

    // Gene summary options (work with or without --em)
    std::string gene_summary_file;
    float min_breadth = 0.10f;
    float min_depth = 0.5f;
    uint32_t min_reads = 3;
    float min_positional_score = 0.0f;
    float min_terminal_ratio = 0.0f;

    // EM reassignment options
    bool use_em = false;
    uint32_t em_max_iters = 100;
    double em_lambda_b = 3.0;
    double em_tol = 1e-4;

    // Alignment-level pre-filters (applied before any processing)
    float aln_min_identity = 0.0f;   // 0 = no filter
    float aln_min_bits = 0.0f;       // 0 = no filter
    float aln_max_evalue = 1e10f;    // very large = no filter

    // Functional profiling options
    std::string kegg_map_file;       // gene_id -> KO mapping
    std::string cazyme_map_file;     // gene_id -> CAZyme family mapping
    std::string viral_map_file;      // gene_id -> VOG mapping
    std::string functional_summary_file;  // Output: per-function stats
    std::string anvio_ko_file;       // Output: Anvi'o-compatible KO abundance

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
        } else if (strcmp(argv[i], "--use-identity") == 0) {
            use_identity = true;
        } else if (strcmp(argv[i], "--no-identity") == 0) {
            use_identity = false;
        } else if (strcmp(argv[i], "--em") == 0) {
            use_em = true;
        } else if (strcmp(argv[i], "--em-iters") == 0 && i + 1 < argc) {
            em_max_iters = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (strcmp(argv[i], "--em-lambda") == 0 && i + 1 < argc) {
            em_lambda_b = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--em-tol") == 0 && i + 1 < argc) {
            em_tol = std::stod(argv[++i]);
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
        } else if (strcmp(argv[i], "--aln-min-identity") == 0 && i + 1 < argc) {
            aln_min_identity = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--aln-min-bits") == 0 && i + 1 < argc) {
            aln_min_bits = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--aln-max-evalue") == 0 && i + 1 < argc) {
            aln_max_evalue = std::stof(argv[++i]);
        } else if (strcmp(argv[i], "--kegg-map") == 0 && i + 1 < argc) {
            kegg_map_file = argv[++i];
        } else if (strcmp(argv[i], "--cazyme-map") == 0 && i + 1 < argc) {
            cazyme_map_file = argv[++i];
        } else if (strcmp(argv[i], "--viral-map") == 0 && i + 1 < argc) {
            viral_map_file = argv[++i];
        } else if (strcmp(argv[i], "--functional-summary") == 0 && i + 1 < argc) {
            functional_summary_file = argv[++i];
        } else if (strcmp(argv[i], "--anvio-ko") == 0 && i + 1 < argc) {
            anvio_ko_file = argv[++i];
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
            std::cerr << "  --no-identity     Disable identity evidence (enabled by default)\n";
            std::cerr << "  -v                Verbose output\n\n";
            std::cerr << "Gene summary (per-gene aggregation, works with or without --em):\n";
            std::cerr << "  --gene-summary FILE  Per-gene statistics TSV\n";
            std::cerr << "  --min-breadth FLOAT  Minimum breadth filter (default: 0.1)\n";
            std::cerr << "  --min-depth FLOAT    Minimum depth filter (default: 0.5)\n";
            std::cerr << "  --min-reads INT      Minimum read count (default: 3)\n\n";
            std::cerr << "Functional profiling (requires --gene-summary):\n";
            std::cerr << "  --kegg-map FILE      Gene-to-KO mapping TSV (gene_id<TAB>KO)\n";
            std::cerr << "  --cazyme-map FILE    Gene-to-CAZyme family mapping TSV\n";
            std::cerr << "  --viral-map FILE     Gene-to-VOG mapping TSV\n";
            std::cerr << "  --functional-summary FILE  Per-function stats TSV\n";
            std::cerr << "  --anvio-ko FILE      Gene-level KO abundance for anvi-estimate-metabolism\n\n";
            std::cerr << "EM reassignment (multi-mapping resolution):\n";
            std::cerr << "  --em              Enable EM multi-mapping resolution\n";
            std::cerr << "  --em-iters INT    Max EM iterations (default: 100)\n";
            std::cerr << "  --em-lambda FLOAT Temperature parameter (default: 3.0)\n";
            std::cerr << "  --em-tol FLOAT    Convergence tolerance (default: 1e-4)\n";
            std::cerr << "  --min-positional-score FLOAT  Min read start diversity (default: 0)\n";
            std::cerr << "                       Filters spurious matches where all reads hit same position\n";
            std::cerr << "                       Score = sqrt(diversity * span), range 0-1\n";
            std::cerr << "  --min-terminal-ratio FLOAT  Min terminal/middle coverage ratio (default: 0)\n";
            std::cerr << "                       Filters matches with reads only in middle (spurious motif)\n\n";
            std::cerr << "Alignment-level pre-filters (applied before any processing):\n";
            std::cerr << "  --aln-min-identity FLOAT  Min identity fraction (default: 0, no filter)\n";
            std::cerr << "  --aln-min-bits FLOAT      Min bit score (default: 0, no filter)\n";
            std::cerr << "  --aln-max-evalue FLOAT    Max e-value (default: 1e10, no filter)\n\n";
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
    std::string line;
    size_t starts[16], ends[16];

    // Alignment filter counters
    size_t total_alignments = 0;
    size_t filtered_identity = 0;
    size_t filtered_bits = 0;
    size_t filtered_evalue = 0;

    // EM mode: also collect CompactAlignment records for all hits
    agp::StringTable read_names, ref_names;
    agp::MmapArray<agp::CompactAlignment> em_alignments;
    size_t em_aln_count = 0;
    if (use_em) {
        em_alignments.resize(1024 * 1024);  // pre-allocate 1M, grows as needed
        em_aln_count = 0;
    }

    // Best hit tracking: query_id -> (bit_score, line_data)
    // MMseqs2 doesn't guarantee best-first ordering, so we keep track of the
    // best hit (highest bit_score) for each query and process them after.
    struct BestHitData {
        float bits;
        float bits_second = 0.0f;  // second-best bit score for this query
        std::string target_id;
        uint32_t ref_idx;  // For EM: ref_idx corresponding to target_id
        float fident;
        float evalue;
        size_t qstart_0;
        size_t tstart_0;
        size_t qlen;
        size_t tlen;
        std::string qaln;
        std::string taln;
    };
    std::unordered_map<std::string, BestHitData> best_hits;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        size_t nf = split_tabs(line, starts, ends, 16);
        if (nf < 16) continue;

        std::string query_id = field_str(line, starts[0], ends[0]);
        std::string target_id = field_str(line, starts[1], ends[1]);
        float fident = std::stof(field_str(line, starts[2], ends[2]));
        size_t qstart = static_cast<size_t>(std::stoul(field_str(line, starts[6], ends[6])));
        size_t tstart = static_cast<size_t>(std::stoul(field_str(line, starts[8], ends[8])));
        // Use stod for evalue: E-values can be extremely small (e.g., 1E-100)
        // which underflows stof but is representable in double
        float evalue = static_cast<float>(std::stod(field_str(line, starts[10], ends[10])));
        float bits = std::stof(field_str(line, starts[11], ends[11]));
        size_t qlen = static_cast<size_t>(std::stoul(field_str(line, starts[12], ends[12])));
        size_t tlen = static_cast<size_t>(std::stoul(field_str(line, starts[13], ends[13])));
        std::string qaln = field_str(line, starts[14], ends[14]);
        std::string taln = field_str(line, starts[15], ends[15]);

        size_t qstart_0 = (qstart > 0) ? qstart - 1 : 0;
        size_t tstart_0 = (tstart > 0) ? tstart - 1 : 0;

        // Alignment-level pre-filtering
        ++total_alignments;
        if (fident < aln_min_identity) { ++filtered_identity; continue; }
        if (bits < aln_min_bits) { ++filtered_bits; continue; }
        if (evalue > aln_max_evalue) { ++filtered_evalue; continue; }

        // EM mode: collect CompactAlignment for every hit
        uint32_t tidx = 0;  // Track ref_idx for best hit comparison
        if (use_em) {
            uint32_t ridx = read_names.insert(query_id);
            tidx = ref_names.insert(target_id);

            // Look up damage score from AGD index
            float ds = 0.0f;
            if (damage_index) {
                if (const AgdRecord* rec = damage_index->find(query_id)) {
                    ds = static_cast<float>(rec->p_damaged_q) / 255.0f;
                }
            }

            // Grow if needed
            if (em_aln_count >= em_alignments.size()) {
                em_alignments.resize(em_alignments.size() * 2);
            }

            agp::CompactAlignment& ca = em_alignments[em_aln_count++];
            ca.read_idx = ridx;
            ca.ref_idx = tidx;
            ca.bit_score = bits;
            ca.damage_score = ds;
            ca.aln_start = static_cast<uint16_t>(std::min(tstart_0, size_t(65535)));
            ca.aln_end = static_cast<uint16_t>(std::min(tstart_0 + taln.size(), size_t(65535)));
            ca.identity_q = static_cast<uint16_t>(fident * 65535.0f);
            ca.flags = 0;
        }

        // Track best hit per query (highest bit_score) and second-best
        auto it = best_hits.find(query_id);
        if (it == best_hits.end()) {
            best_hits[query_id] = BestHitData{
                bits, 0.0f, target_id, tidx, fident, evalue,
                qstart_0, tstart_0, qlen, tlen,
                qaln, taln
            };
        } else if (bits > it->second.bits) {
            it->second.bits_second = it->second.bits;  // demote current best
            it->second.bits = bits;
            it->second.target_id = target_id;
            it->second.ref_idx = tidx;
            it->second.fident = fident;
            it->second.evalue = evalue;
            it->second.qstart_0 = qstart_0;
            it->second.tstart_0 = tstart_0;
            it->second.qlen = qlen;
            it->second.tlen = tlen;
            it->second.qaln = qaln;
            it->second.taln = taln;
        } else if (bits > it->second.bits_second) {
            it->second.bits_second = bits;  // new second-best
        }
    }
    in.close();

    if (verbose) {
        size_t passed = total_alignments - filtered_identity - filtered_bits - filtered_evalue;
        std::cerr << "Alignment pre-filtering:\n";
        std::cerr << "  Total alignments: " << total_alignments << "\n";
        if (filtered_identity > 0 || aln_min_identity > 0)
            std::cerr << "  Filtered (identity < " << aln_min_identity << "): " << filtered_identity << "\n";
        if (filtered_bits > 0 || aln_min_bits > 0)
            std::cerr << "  Filtered (bits < " << aln_min_bits << "): " << filtered_bits << "\n";
        if (filtered_evalue > 0 || aln_max_evalue < 1e9f)
            std::cerr << "  Filtered (evalue > " << aln_max_evalue << "): " << filtered_evalue << "\n";
        std::cerr << "  Passed: " << passed << " (" << (100.0 * passed / total_alignments) << "%)\n\n";
    }

    // Process best hits into summaries
    for (auto& [query_id, hit] : best_hits) {
        auto summary = annotate_alignment(
            query_id, hit.target_id, hit.qaln, hit.taln,
            hit.qstart_0, hit.tstart_0, hit.qlen,
            hit.evalue, hit.bits, hit.fident, d_max, lambda);

        if (summary.total_mismatches > 0) {
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

    // EM reassignment
    agp::EMState em_state;
    agp::ReferenceStatsCollector stats_collector;
    // Per-read EM results: query_id -> (gamma, gamma_ancient, best_hit)
    std::unordered_map<std::string, std::tuple<float, float, bool>> em_read_results;

    if (use_em && em_aln_count > 0) {
        uint32_t n_reads = static_cast<uint32_t>(read_names.size());
        uint32_t n_refs = static_cast<uint32_t>(ref_names.size());

        if (verbose) {
            std::cerr << "\nEM reassignment:\n";
            std::cerr << "  Alignments: " << em_aln_count << "\n";
            std::cerr << "  Reads: " << n_reads << "\n";
            std::cerr << "  References: " << n_refs << "\n";
            std::cerr << "  Lambda_b: " << em_lambda_b << "\n";
            std::cerr << "  Max iters: " << em_max_iters << "\n";
            std::cerr << "  Tolerance: " << em_tol << "\n";
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
        em_params.use_damage = (damage_index != nullptr);

        // Run EM with SQUAREM acceleration
        em_state = agp::squarem_em(aln_data, em_params, &stats_collector);

        if (verbose) {
            std::cerr << "  Converged in " << em_state.iterations << " iterations\n";
            std::cerr << "  Log-likelihood: " << std::fixed << std::setprecision(2)
                      << em_state.log_likelihood << "\n";
            if (em_params.use_damage) {
                std::cerr << "  Estimated pi (ancient fraction): "
                          << std::fixed << std::setprecision(4) << em_state.pi << "\n";
            }
        }

        // Extract best assignments and build per-read lookup
        auto best = agp::reassign_reads(aln_data, em_state, 0.01);

        // Build per-read lookup: for each read, find its gamma for the best hit
        // Summaries use best-hit-per-read (by bit_score), so we map read_idx -> gamma of best assignment
        for (uint32_t r = 0; r < n_reads; ++r) {
            std::string_view rname = read_names.get(r);
            std::string rname_str(rname);

            uint32_t start = aln_data.read_offsets[r];
            uint32_t end = aln_data.read_offsets[r + 1];

            // Find gamma for the best alignment per read (by EM assignment)
            float best_gamma = 0.0f;
            float best_gamma_ancient = 0.0f;
            uint32_t best_ref = best[r].first;
            bool bitscore_best_matches_em = false;

            for (uint32_t j = start; j < end; ++j) {
                if (aln_data.alignments[j].ref_idx == best_ref) {
                    best_gamma = static_cast<float>(em_state.gamma[j]);
                    if (em_params.use_damage && !em_state.gamma_ancient.empty()) {
                        best_gamma_ancient = static_cast<float>(em_state.gamma_ancient[j]);
                    }
                    break;
                }
            }

            // Check if the best-by-bitscore hit matches the EM-assigned best
            auto bh_it = best_hits.find(rname_str);
            if (bh_it != best_hits.end() && bh_it->second.ref_idx == best_ref) {
                bitscore_best_matches_em = true;
            }

            em_read_results[rname_str] = {best_gamma, best_gamma_ancient, bitscore_best_matches_em};
        }

        // Aggregate per-target short protein metrics for gene summary
        struct GeneMetrics {
            float sum_delta_bits = 0.0f;
            float sum_bpa = 0.0f;
            uint32_t count = 0;
            std::array<uint32_t, 4> bin_counts{};
        };
        std::unordered_map<std::string, GeneMetrics> gene_metrics;
        for (const auto& s : summaries) {
            auto& gm = gene_metrics[s.target_id];
            gm.sum_delta_bits += s.delta_bits;
            gm.sum_bpa += s.bpa;
            gm.count++;
            gm.bin_counts[s.length_bin]++;
        }

        // Write gene summary
        if (!gene_summary_file.empty()) {
            agp::FilterThresholds thresh;
            thresh.min_reads = min_reads;
            thresh.min_breadth = min_breadth;
            thresh.min_depth = min_depth;
            thresh.min_positional_score = min_positional_score;
            thresh.min_terminal_ratio = min_terminal_ratio;

            auto ref_stats = stats_collector.finalize_filtered(
                static_cast<float>(em_state.pi), thresh);

            // Replace numeric ref_idx IDs with actual names
            for (auto& s : ref_stats) {
                uint32_t idx = static_cast<uint32_t>(std::stoul(s.ref_id));
                if (idx < ref_names.size()) {
                    s.ref_id = std::string(ref_names.get(idx));
                }
            }

            std::ofstream gf(gene_summary_file);
            if (!gf.is_open()) {
                std::cerr << "Error: Cannot open gene summary file: " << gene_summary_file << "\n";
                return 1;
            }

            gf << "gene_id\tn_reads\tn_effective\tn_ancient\tn_modern\t"
                  "breadth\tdepth_mean\tavg_identity\tenrichment\t"
                  "start_diversity\tstart_span\tpositional_score\t"
                  "terminal_ratio\tn_unique_starts\tis_ancient\t"
                  "avg_delta_bits\tavg_bpa\tlength_bin\n";

            for (const auto& s : ref_stats) {
                bool is_ancient = s.damage_enrichment > 0.0f;
                gf << s.ref_id
                   << '\t' << s.n_reads
                   << '\t' << std::fixed << std::setprecision(2) << s.n_effective
                   << '\t' << std::setprecision(2) << s.n_ancient
                   << '\t' << std::setprecision(2) << s.n_modern
                   << '\t' << std::setprecision(4) << s.breadth
                   << '\t' << std::setprecision(2) << s.depth_mean
                   << '\t' << std::setprecision(4) << s.avg_identity
                   << '\t' << std::setprecision(4) << s.damage_enrichment
                   << '\t' << std::setprecision(4) << s.start_diversity
                   << '\t' << std::setprecision(4) << s.start_span
                   << '\t' << std::setprecision(4) << s.positional_score
                   << '\t' << std::setprecision(4) << s.terminal_ratio
                   << '\t' << s.n_unique_starts
                   << '\t' << (is_ancient ? 1 : 0);
                auto gm_it = gene_metrics.find(s.ref_id);
                if (gm_it != gene_metrics.end() && gm_it->second.count > 0) {
                    const auto& gm = gm_it->second;
                    gf << '\t' << std::setprecision(1) << (gm.sum_delta_bits / gm.count)
                       << '\t' << std::setprecision(3) << (gm.sum_bpa / gm.count)
                       << '\t' << static_cast<int>(std::max_element(
                              gm.bin_counts.begin(),
                              gm.bin_counts.end()) - gm.bin_counts.begin());
                } else {
                    gf << "\t.\t.\t.";
                }
                gf << '\n';
            }

            if (verbose) {
                std::cerr << "Gene summary written to: " << gene_summary_file
                          << " (" << ref_stats.size() << " genes after filtering)\n";
            }
        }

        if (verbose) {
            // Filtering stats
            size_t multi = 0, unique = 0;
            for (uint32_t r = 0; r < n_reads; ++r) {
                uint32_t deg = aln_data.read_offsets[r + 1] - aln_data.read_offsets[r];
                if (deg > 1) multi++;
                else unique++;
            }
            std::cerr << "  Unique mappers: " << unique << "\n";
            std::cerr << "  Multi-mappers: " << multi << "\n\n";
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
    score_params.use_identity = use_identity;
    score_params.k_identity = 20.0f;
    score_params.identity_baseline = 0.90f;

    // Damage informativeness gating from AGD v3 header
    // When uninformative: terminal evidence is zeroed, posterior output is NA
    float damage_detectability = 1.0f;
    if (damage_index) {
        score_params.damage_informative = damage_index->damage_informative();
        score_params.channel_b_valid = damage_index->channel_b_valid();
        damage_detectability = agp::compute_damage_detectability(
            damage_index->d_max(), damage_index->stop_decay_llr(),
            damage_index->terminal_shift(), damage_index->damage_validated(),
            damage_index->damage_artifact());
    }

    if (verbose) {
        std::cerr << "\nBayesian scoring (empirical Bayes):\n";
        std::cerr << "  Modern proxy reads: " << modern_est.n_reads << "\n";
        std::cerr << "  Pooled opportunities: " << modern_est.M0 << "\n";
        std::cerr << "  Pooled hits: " << modern_est.K0 << "\n";
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
            "max_p_damage\tp_damaged\tp_read\tposterior\t"
            "classification\tlogBF_terminal\tlogBF_sites\tlogBF_identity\t"
            "m_opportunities\tk_hits\tq_eff\ttier\tdamage_class\t"
            "syn_5prime\tsyn_3prime\t"
            "delta_bits\tbpa\tlength_bin\tz_bpa\t"
            "damage_informative";
    if (use_em) {
        *out << "\tgamma\tgamma_ancient\tbest_hit";
    }
    *out << "\n";

    // Store per-read posteriors and classifications for gene summary aggregation
    std::vector<float> read_posteriors(summaries.size());
    std::vector<agp::AncientClassification> read_classifications(summaries.size());

    for (size_t si = 0; si < summaries.size(); ++si) {
        const auto& s = summaries[si];
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

        if (bayes_out.informative) {
            *out << std::fixed << std::setprecision(4) << posterior << "\t";
        } else {
            *out << "NA\t";
        }
        *out << agp::classification_name(classification) << "\t";
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
            auto it = em_read_results.find(s.query_id);
            if (it != em_read_results.end()) {
                auto& [gamma, gamma_anc, best] = it->second;
                *out << "\t" << std::fixed << std::setprecision(4) << gamma
                     << "\t" << std::fixed << std::setprecision(4) << gamma_anc
                     << "\t" << (best ? 1 : 0);
            } else {
                *out << "\t.\t.\t.";
            }
        }
        *out << "\n";
    }

    // Write gene summary (works without EM — uses best-hit per read)
    if (!gene_summary_file.empty() && !use_em) {
        std::unordered_map<std::string, GeneSummaryAccumulator> gene_agg;

        for (size_t i = 0; i < summaries.size(); ++i) {
            const auto& s = summaries[i];
            auto& acc = gene_agg[s.target_id];
            if (acc.tlen == 0) {
                acc.target_id = s.target_id;
                acc.tlen = static_cast<uint32_t>(s.tlen);
                if (acc.tlen > 0) {
                    acc.coverage.assign(acc.tlen, 0.0f);
                }
            }
            size_t aln_on_target = 0;
            for (char c : s.taln) {
                if (c != '-') ++aln_on_target;
            }
            acc.add_read(read_classifications[i], read_posteriors[i],
                         s.fident, s.tstart_0, aln_on_target);
        }

        std::ofstream gf(gene_summary_file);
        if (!gf.is_open()) {
            std::cerr << "Error: Cannot open gene summary file: " << gene_summary_file << "\n";
            return 1;
        }

        gf << "target_id\tn_reads\tn_ancient_conf\tn_ancient_likely\t"
              "n_undetermined\tn_modern_conf\tancient_frac\tci_low\tci_high\t"
              "mean_posterior\tmean_identity\tbreadth\tdepth\n";

        for (const auto& [tid, acc] : gene_agg) {
            if (acc.n_reads < min_reads) continue;
            float b = acc.breadth();
            if (b < min_breadth) continue;
            float d = acc.depth_mean();
            if (d < min_depth) continue;

            uint32_t n_ancient = acc.n_ancient_conf + acc.n_ancient_likely;
            auto [ci_lo, ci_hi] = wilson_ci(n_ancient, acc.n_reads);

            gf << acc.target_id
               << '\t' << acc.n_reads
               << '\t' << acc.n_ancient_conf
               << '\t' << acc.n_ancient_likely
               << '\t' << acc.n_undetermined
               << '\t' << acc.n_modern_conf
               << '\t' << std::fixed << std::setprecision(4)
               << (acc.n_reads > 0 ? static_cast<float>(n_ancient) / acc.n_reads : 0.0f)
               << '\t' << std::setprecision(4) << ci_lo
               << '\t' << std::setprecision(4) << ci_hi
               << '\t' << std::setprecision(4) << acc.avg_posterior()
               << '\t' << std::setprecision(4) << acc.avg_fident()
               << '\t' << std::setprecision(4) << b
               << '\t' << std::setprecision(2) << d
               << '\n';
        }

        if (verbose) {
            size_t total = gene_agg.size();
            size_t passed = 0;
            for (const auto& [tid, acc] : gene_agg) {
                if (acc.n_reads >= min_reads && acc.breadth() >= min_breadth
                    && acc.depth_mean() >= min_depth)
                    ++passed;
            }
            std::cerr << "Gene summary written to: " << gene_summary_file
                      << " (" << passed << "/" << total << " genes after filtering)\n";
        }

        // Functional profiling (aggregate gene stats by function)
        bool do_functional = !kegg_map_file.empty() || !cazyme_map_file.empty() || !viral_map_file.empty();
        if (do_functional) {
            // Load mapping files
            auto kegg_map = kegg_map_file.empty() ? std::unordered_map<std::string, std::string>{}
                          : load_mapping_file(kegg_map_file, verbose);
            auto cazyme_map = cazyme_map_file.empty() ? std::unordered_map<std::string, std::string>{}
                            : load_mapping_file(cazyme_map_file, verbose);
            auto viral_map = viral_map_file.empty() ? std::unordered_map<std::string, std::string>{}
                           : load_mapping_file(viral_map_file, verbose);

            // Accumulators by function
            std::unordered_map<std::string, FunctionalAccumulator> kegg_acc;    // KO -> stats
            std::unordered_map<std::string, FunctionalAccumulator> cazyme_fam_acc;  // Family -> stats
            std::unordered_map<std::string, FunctionalAccumulator> cazyme_cls_acc;  // Class -> stats
            std::unordered_map<std::string, FunctionalAccumulator> viral_acc;   // VOG -> stats

            // Aggregate gene stats by function
            for (const auto& [tid, acc] : gene_agg) {
                // Skip filtered genes
                if (acc.n_reads < min_reads || acc.breadth() < min_breadth || acc.depth_mean() < min_depth)
                    continue;

                float gene_ancient = static_cast<float>(acc.n_ancient_conf + acc.n_ancient_likely);
                float gene_modern = static_cast<float>(acc.n_modern_conf);
                float gene_undetermined = static_cast<float>(acc.n_undetermined);

                // KEGG
                auto kegg_it = kegg_map.find(tid);
                if (kegg_it != kegg_map.end()) {
                    auto& fa = kegg_acc[kegg_it->second];
                    if (fa.function_id.empty()) {
                        fa.function_id = kegg_it->second;
                        fa.db_type = "KEGG";
                    }
                    fa.add_gene(acc.n_reads, gene_ancient, gene_modern, gene_undetermined, acc.avg_posterior());
                }

                // CAZyme
                auto cazyme_it = cazyme_map.find(tid);
                if (cazyme_it != cazyme_map.end()) {
                    const std::string& family = cazyme_it->second;
                    // Family level
                    auto& fam_acc = cazyme_fam_acc[family];
                    if (fam_acc.function_id.empty()) {
                        fam_acc.function_id = family;
                        fam_acc.db_type = "CAZyme";
                    }
                    fam_acc.add_gene(acc.n_reads, gene_ancient, gene_modern, gene_undetermined, acc.avg_posterior());

                    // Class level (GH, GT, PL, CE, AA, CBM)
                    std::string cls = cazyme_class_from_family(family);
                    if (!cls.empty()) {
                        auto& cls_acc = cazyme_cls_acc[cls];
                        if (cls_acc.function_id.empty()) {
                            cls_acc.function_id = cls;
                            cls_acc.db_type = "CAZyme_class";
                        }
                        cls_acc.add_gene(acc.n_reads, gene_ancient, gene_modern, gene_undetermined, acc.avg_posterior());
                    }
                }

                // Viral
                auto viral_it = viral_map.find(tid);
                if (viral_it != viral_map.end()) {
                    auto& va = viral_acc[viral_it->second];
                    if (va.function_id.empty()) {
                        va.function_id = viral_it->second;
                        va.db_type = "Viral";
                    }
                    va.add_gene(acc.n_reads, gene_ancient, gene_modern, gene_undetermined, acc.avg_posterior());
                }
            }

            // Write functional summary
            if (!functional_summary_file.empty()) {
                std::ofstream ff(functional_summary_file);
                if (!ff.is_open()) {
                    std::cerr << "Error: Cannot open functional summary file: " << functional_summary_file << "\n";
                } else {
                    ff << "db\tfunction_id\tlevel\tn_genes\tn_reads\tn_ancient\tn_modern\tn_undetermined\t"
                       << "ancient_frac\tci_low\tci_high\tmean_posterior\n";

                    auto write_acc = [&ff](const std::unordered_map<std::string, FunctionalAccumulator>& accs,
                                           const std::string& level) {
                        for (const auto& [fid, fa] : accs) {
                            float frac = fa.ancient_frac();
                            auto [ci_low, ci_high] = wilson_ci(
                                static_cast<uint32_t>(fa.n_ancient),
                                static_cast<uint32_t>(fa.n_ancient + fa.n_modern));
                            ff << fa.db_type << '\t' << fa.function_id << '\t' << level << '\t'
                               << fa.n_genes << '\t' << fa.n_reads << '\t'
                               << std::fixed << std::setprecision(1) << fa.n_ancient << '\t'
                               << fa.n_modern << '\t' << fa.n_undetermined << '\t'
                               << std::setprecision(4) << frac << '\t' << ci_low << '\t' << ci_high << '\t'
                               << fa.mean_posterior() << '\n';
                        }
                    };

                    write_acc(kegg_acc, "KO");
                    write_acc(cazyme_fam_acc, "family");
                    write_acc(cazyme_cls_acc, "class");
                    write_acc(viral_acc, "VOG");

                    if (verbose) {
                        size_t total_funcs = kegg_acc.size() + cazyme_fam_acc.size() +
                                             cazyme_cls_acc.size() + viral_acc.size();
                        std::cerr << "Functional summary written to: " << functional_summary_file
                                  << " (" << total_funcs << " functions)\n";
                    }
                }
            }

            // Write Anvi'o-compatible gene abundance (for anvi-estimate-metabolism --enzymes-txt)
            // Format: gene_id  enzyme_accession  source  coverage  detection
            // Where coverage = depth (abundance estimate), detection = breadth
            if (!anvio_ko_file.empty() && !kegg_map.empty()) {
                std::ofstream af(anvio_ko_file);
                if (!af.is_open()) {
                    std::cerr << "Error: Cannot open Anvi'o KO file: " << anvio_ko_file << "\n";
                } else {
                    af << "gene_id\tenzyme_accession\tsource\tcoverage\tdetection\n";
                    size_t n_written = 0;
                    for (const auto& [tid, acc] : gene_agg) {
                        // Skip filtered genes
                        if (acc.n_reads < min_reads || acc.breadth() < min_breadth || acc.depth_mean() < min_depth)
                            continue;
                        // Look up KO for this gene
                        auto it = kegg_map.find(tid);
                        if (it == kegg_map.end()) continue;
                        std::string ko = it->second;
                        // Strip "ko:" prefix if present
                        if (ko.size() > 3 && ko.substr(0, 3) == "ko:") {
                            ko = ko.substr(3);
                        }
                        af << tid << '\t' << ko << '\t' << "AGP" << '\t'
                           << std::fixed << std::setprecision(4) << acc.depth_mean() << '\t'
                           << std::setprecision(4) << acc.breadth() << '\n';
                        ++n_written;
                    }
                    if (verbose) {
                        std::cerr << "Anvi'o gene abundance written to: " << anvio_ko_file
                                  << " (" << n_written << " genes with KO annotations)\n";
                    }
                }
            }
        }
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
            agg.sum_delta_bits += s.delta_bits;
            agg.sum_bpa += s.bpa;
            agg.length_bin_counts[s.length_bin]++;

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
              "lrt_pvalue\tlog_bf\td_aa\t"
              "avg_delta_bits\tavg_bpa\tlength_bin\n";

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
               << std::fixed << std::setprecision(4) << d_aa << "\t"
               << std::fixed << std::setprecision(1)
               << (agg.n_reads > 0 ? agg.sum_delta_bits / agg.n_reads : 0.0f) << "\t"
               << std::fixed << std::setprecision(3)
               << (agg.n_reads > 0 ? agg.sum_bpa / agg.n_reads : 0.0f) << "\t"
               << static_cast<int>(std::max_element(
                      agg.length_bin_counts.begin(),
                      agg.length_bin_counts.end()) - agg.length_bin_counts.begin())
               << "\n";
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
