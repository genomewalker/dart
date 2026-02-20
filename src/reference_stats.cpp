#include "dart/reference_stats.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace dart {

// Write per-reference stats as TSV
void write_reference_stats_tsv(std::ostream& out,
                               const std::vector<ReferenceStats>& stats) {
    out << "ref_id\tref_length\tn_reads\tbreadth\tdepth_mean\tdepth_std"
        << "\tdepth_evenness\tstart_diversity\tstart_span\tpositional_score"
        << "\tterminal_ratio\tn_unique_starts"
        << "\tavg_identity\tavg_aln_len\tavg_read_len"
        << "\tn_effective\tn_ancient\tn_modern\tdamage_enrichment\n";

    out << std::fixed;
    for (const auto& s : stats) {
        out << s.ref_id
            << '\t' << s.ref_length
            << '\t' << s.n_reads
            << '\t' << std::setprecision(4) << s.breadth
            << '\t' << std::setprecision(2) << s.depth_mean
            << '\t' << std::setprecision(2) << s.depth_std
            << '\t' << std::setprecision(4) << s.depth_evenness
            << '\t' << std::setprecision(4) << s.start_diversity
            << '\t' << std::setprecision(4) << s.start_span
            << '\t' << std::setprecision(4) << s.positional_score
            << '\t' << std::setprecision(4) << s.terminal_ratio
            << '\t' << s.n_unique_starts
            << '\t' << std::setprecision(4) << s.avg_identity
            << '\t' << std::setprecision(1) << s.avg_alignment_length
            << '\t' << std::setprecision(1) << s.avg_read_length
            << '\t' << std::setprecision(2) << s.n_effective
            << '\t' << std::setprecision(2) << s.n_ancient
            << '\t' << std::setprecision(2) << s.n_modern
            << '\t' << std::setprecision(4) << s.damage_enrichment
            << '\n';
    }
}

// Write JSON summary of sample-level aggregates
void write_reference_stats_json(std::ostream& out,
                                const std::vector<ReferenceStats>& stats,
                                float pi) {
    // Aggregate across all references
    uint64_t total_reads = 0;
    double total_effective = 0.0;
    double total_ancient = 0.0;
    double total_modern = 0.0;
    uint32_t n_refs = static_cast<uint32_t>(stats.size());
    uint32_t n_ancient_enriched = 0;

    double sum_breadth = 0.0;
    double sum_depth = 0.0;
    double sum_identity = 0.0;
    uint32_t n_with_reads = 0;

    for (const auto& s : stats) {
        total_reads += s.n_reads;
        total_effective += s.n_effective;
        total_ancient += s.n_ancient;
        total_modern += s.n_modern;
        if (s.damage_enrichment > 0.0f) ++n_ancient_enriched;
        if (s.n_reads > 0) {
            sum_breadth += s.breadth;
            sum_depth += s.depth_mean;
            sum_identity += s.avg_identity;
            ++n_with_reads;
        }
    }

    float avg_breadth = (n_with_reads > 0)
        ? static_cast<float>(sum_breadth / n_with_reads) : 0.0f;
    float avg_depth = (n_with_reads > 0)
        ? static_cast<float>(sum_depth / n_with_reads) : 0.0f;
    float avg_identity = (n_with_reads > 0)
        ? static_cast<float>(sum_identity / n_with_reads) : 0.0f;
    float ancient_frac = (total_effective > 0.0)
        ? static_cast<float>(total_ancient / total_effective) : 0.0f;

    out << std::fixed;
    out << "{\n";
    out << "  \"references\": {\n";
    out << "    \"total\": " << n_refs << ",\n";
    out << "    \"with_reads\": " << n_with_reads << ",\n";
    out << "    \"ancient_enriched\": " << n_ancient_enriched << "\n";
    out << "  },\n";
    out << "  \"reads\": {\n";
    out << "    \"total\": " << total_reads << ",\n";
    out << "    \"n_effective\": " << std::setprecision(1) << total_effective << ",\n";
    out << "    \"n_ancient\": " << std::setprecision(1) << total_ancient << ",\n";
    out << "    \"n_modern\": " << std::setprecision(1) << total_modern << ",\n";
    out << "    \"ancient_fraction\": " << std::setprecision(4) << ancient_frac << "\n";
    out << "  },\n";
    out << "  \"quality\": {\n";
    out << "    \"avg_breadth\": " << std::setprecision(4) << avg_breadth << ",\n";
    out << "    \"avg_depth\": " << std::setprecision(2) << avg_depth << ",\n";
    out << "    \"avg_identity\": " << std::setprecision(4) << avg_identity << "\n";
    out << "  },\n";
    out << "  \"prior_pi\": " << std::setprecision(4) << pi << "\n";
    out << "}\n";
}

} // namespace dart
