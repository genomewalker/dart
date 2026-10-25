#pragma once

#include <cstdint>

namespace dart {

inline bool null_abundance_rejected(bool enabled, uint8_t reason) {
    return enabled && (reason & 0x1) != 0;
}

// Shared by best-hit aggregation and the summary-only fallback.
template<class Hits, class Decisions, class Keep, class Drop>
void visit_best_abundance_reads(const Hits& hits, const Decisions& decisions,
                               bool use_em, bool drop_null_flagged, Keep keep, Drop drop) {
    for (uint32_t r = 0; r < hits.size(); ++r) {
        if (!hits[r].has) continue;
        if (use_em && (r >= decisions.size() || !decisions[r].em_keep)) continue;
        if (r < decisions.size() &&
            null_abundance_rejected(drop_null_flagged, decisions[r].em_reject_reason)) {
            drop(r);
            continue;
        }
        keep(r);
    }
}

} // namespace dart
