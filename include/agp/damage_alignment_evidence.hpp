#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string_view>

namespace agp {

struct DamageEvidenceSub {
    char target_aa;
    char query_aa;
    char damage_class;  // 'C' for C->T, 'G' for G->A
    bool high_confidence;
};

// Keep this table aligned with damage-annotate substitution logic.
static constexpr DamageEvidenceSub DAMAGE_EVIDENCE_SUBS[] = {
    // C->T damage (5' end)
    {'R', 'W', 'C', true}, {'R', 'C', 'C', true}, {'R', '*', 'C', true},
    {'Q', '*', 'C', true}, {'H', 'Y', 'C', true}, {'P', 'S', 'C', false},
    {'P', 'L', 'C', false}, {'T', 'I', 'C', false}, {'T', 'M', 'C', false},
    {'A', 'V', 'C', false}, {'S', 'F', 'C', false}, {'S', 'L', 'C', false},

    // G->A damage (3' end)
    {'E', 'K', 'G', true}, {'D', 'N', 'G', true}, {'A', 'T', 'G', false},
    {'G', 'R', 'G', false}, {'G', 'S', 'G', false}, {'V', 'I', 'G', false},
    {'V', 'M', 'G', false}, {'R', 'K', 'G', true}, {'R', 'Q', 'G', true},
    {'R', 'H', 'G', true}, {'G', 'E', 'G', true}, {'G', 'D', 'G', true},
    {'S', 'N', 'G', false}, {'C', 'Y', 'G', false}, {'W', '*', 'G', false},

    // X-masked positions
    {'Q', 'X', 'C', true}, {'R', 'X', 'C', true}, {'W', 'X', 'G', true},
};

static constexpr size_t NUM_DAMAGE_EVIDENCE_SUBS =
    sizeof(DAMAGE_EVIDENCE_SUBS) / sizeof(DAMAGE_EVIDENCE_SUBS[0]);

struct DamageEvidenceLookup {
    std::array<std::array<uint8_t, 128>, 128> table{};
    constexpr DamageEvidenceLookup() {
        for (size_t i = 0; i < NUM_DAMAGE_EVIDENCE_SUBS; ++i) {
            const auto t = static_cast<unsigned char>(DAMAGE_EVIDENCE_SUBS[i].target_aa);
            const auto q = static_cast<unsigned char>(DAMAGE_EVIDENCE_SUBS[i].query_aa);
            if (t < 128 && q < 128) {
                table[t][q] = static_cast<uint8_t>(i + 1);
            }
        }
    }
};

struct DamageEvidenceTargetClassLookup {
    std::array<uint8_t, 128> class_mask{};  // bit0: C->T, bit1: G->A
    constexpr DamageEvidenceTargetClassLookup() {
        for (size_t i = 0; i < NUM_DAMAGE_EVIDENCE_SUBS; ++i) {
            const auto t = static_cast<unsigned char>(DAMAGE_EVIDENCE_SUBS[i].target_aa);
            if (t >= 128) continue;
            if (DAMAGE_EVIDENCE_SUBS[i].damage_class == 'C') class_mask[t] |= 0x1;
            if (DAMAGE_EVIDENCE_SUBS[i].damage_class == 'G') class_mask[t] |= 0x2;
        }
    }
};

static constexpr DamageEvidenceLookup DAMAGE_EVIDENCE_LOOKUP{};
static constexpr DamageEvidenceTargetClassLookup DAMAGE_EVIDENCE_TARGET_CLASSES{};

inline const DamageEvidenceSub* find_damage_evidence_sub(char target_aa, char query_aa) {
    const auto t = static_cast<unsigned char>(target_aa);
    const auto q = static_cast<unsigned char>(query_aa);
    if (t >= 128 || q >= 128) return nullptr;
    const uint8_t idx = DAMAGE_EVIDENCE_LOOKUP.table[t][q];
    return idx ? &DAMAGE_EVIDENCE_SUBS[idx - 1] : nullptr;
}

inline bool query_id_is_reverse(std::string_view query_id) {
    // Expected AGP suffix: ..._<+|->_<0..9>
    if (query_id.size() < 4) return false;
    for (size_t i = query_id.size() - 1; i >= 2; --i) {
        if (query_id[i] == '_' &&
            query_id[i - 2] == '_' &&
            (query_id[i - 1] == '+' || query_id[i - 1] == '-')) {
            return query_id[i - 1] == '-';
        }
        if (query_id[i] != '_' && query_id[i] != '+' && query_id[i] != '-' &&
            (query_id[i] < '0' || query_id[i] > '9')) {
            break;
        }
    }
    return false;
}

inline float positional_damage_prob(size_t dist_from_terminus, float d_max, float lambda) {
    if (d_max <= 0.0f) return 0.0f;
    return d_max * std::exp(-lambda * static_cast<float>(dist_from_terminus));
}

struct AlignmentDamageEvidence {
    uint16_t k_hits = 0;   // damage-consistent hits among opportunities
    uint16_t m_sites = 0;  // opportunities in aligned columns
    float ll_ancient = 0.0f;
    float ll_modern = 0.0f;
};

inline AlignmentDamageEvidence compute_alignment_damage_evidence(
    std::string_view query_id,
    std::string_view qaln,
    std::string_view taln,
    uint16_t qstart_1based,
    uint16_t qlen_aa,
    float d_max,
    float lambda,
    float q_modern = 0.005f)
{
    AlignmentDamageEvidence out{};
    const bool reverse = query_id_is_reverse(query_id);
    const double eps = 1e-8;
    const double qM = std::clamp(static_cast<double>(q_modern), eps, 1.0 - eps);
    size_t q_pos = (qstart_1based > 0) ? static_cast<size_t>(qstart_1based - 1) : 0;

    const size_t n = std::min(qaln.size(), taln.size());
    for (size_t i = 0; i < n; ++i) {
        const char q_aa = qaln[i];
        const char t_aa = taln[i];
        if (q_aa == '-') continue;
        if (t_aa == '-') {
            q_pos++;
            continue;
        }

        size_t dist_5 = 0;
        size_t dist_3 = 0;
        if (!reverse) {
            dist_5 = q_pos * 3;
            dist_3 = (qlen_aa > 0 && q_pos < qlen_aa) ? (static_cast<size_t>(qlen_aa) - 1 - q_pos) * 3 : 0;
        } else {
            dist_5 = (qlen_aa > 0 && q_pos < qlen_aa) ? (static_cast<size_t>(qlen_aa) - 1 - q_pos) * 3 : 0;
            dist_3 = q_pos * 3;
        }

        const auto t_code = static_cast<unsigned char>(t_aa);
        const uint8_t class_mask =
            (t_code < 128) ? DAMAGE_EVIDENCE_TARGET_CLASSES.class_mask[t_code] : 0;
        if (class_mask != 0) {
            const double q_ct = (class_mask & 0x1)
                ? static_cast<double>(positional_damage_prob(dist_5, d_max, lambda))
                : 0.0;
            const double q_ga = (class_mask & 0x2)
                ? static_cast<double>(positional_damage_prob(dist_3, d_max, lambda))
                : 0.0;

            out.m_sites = static_cast<uint16_t>(std::min<uint32_t>(65535u, static_cast<uint32_t>(out.m_sites) + 1u));

            const DamageEvidenceSub* sub = nullptr;
            if (q_aa != t_aa) {
                sub = find_damage_evidence_sub(t_aa, q_aa);
            }

            if (sub) {
                out.k_hits = static_cast<uint16_t>(std::min<uint32_t>(65535u, static_cast<uint32_t>(out.k_hits) + 1u));
                const double qA = (sub->damage_class == 'C') ? q_ct : q_ga;
                const double pA = std::clamp(qA + qM, eps, 1.0 - eps);
                out.ll_ancient += static_cast<float>(std::log(pA));
                out.ll_modern += static_cast<float>(std::log(qM));
            } else {
                const double qA = std::max(q_ct, q_ga);
                const double pA_no = std::clamp(1.0 - (qA + qM), eps, 1.0 - eps);
                out.ll_ancient += static_cast<float>(std::log(pA_no));
                out.ll_modern += static_cast<float>(std::log(1.0 - qM));
            }
        }

        q_pos++;
    }

    return out;
}

}  // namespace agp

