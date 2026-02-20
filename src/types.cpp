#include "dart/types.hpp"
#include <cmath>
#include <algorithm>

namespace dart {

float DamageProfile::get_ct_damage(Position pos, Position seq_len) const {
    // Calculate distance from 5' end
    float dist_5prime = static_cast<float>(pos);

    // Exponential decay from 5' end
    float damage_5prime = delta_max * std::exp(-lambda_5prime * dist_5prime) + delta_background;

    // Use precomputed profile if available
    if (!ct_prob_5prime.empty() && pos < ct_prob_5prime.size()) {
        return ct_prob_5prime[pos];
    }

    return std::clamp(damage_5prime, 0.0f, 1.0f);
}

float DamageProfile::get_ga_damage(Position pos, Position seq_len) const {
    // Calculate distance from 3' end
    float dist_3prime = static_cast<float>(seq_len - pos - 1);

    // Exponential decay from 3' end
    float damage_3prime = delta_max * std::exp(-lambda_3prime * dist_3prime) + delta_background;

    // Use precomputed profile if available
    if (!ga_prob_3prime.empty() && (seq_len - pos - 1) < ga_prob_3prime.size()) {
        return ga_prob_3prime[seq_len - pos - 1];
    }

    return std::clamp(damage_3prime, 0.0f, 1.0f);
}

} // namespace dart
