// DamageProfile creation and caching

#include "agp/damage_model.hpp"
#include "agp/frame_selector.hpp"
#include <cmath>
#include <algorithm>

namespace agp {

// Forward declaration from analyzer.cpp
DamageModel DamageAnalyzer_analyze(const std::vector<std::string>& sequences);

DamageModel DamageModel::from_parameters(float lambda_5prime,
                                        float lambda_3prime,
                                        float delta_max,
                                        float delta_background) {
    DamageModel model;
    model.lambda_5prime_ = lambda_5prime;
    model.lambda_3prime_ = lambda_3prime;
    model.delta_max_ = delta_max;
    model.delta_background_ = delta_background;
    return model;
}

void DamageModel::update_from_sample_profile(const SampleDamageProfile& profile) {
    // Update decay rates from Pass 1 observations
    lambda_5prime_ = profile.lambda_5prime;
    lambda_3prime_ = profile.lambda_3prime;

    // Update maximum damage rates from terminal positions
    delta_max_ = std::max(profile.max_damage_5prime, profile.max_damage_3prime);

    // Update background rate from internal positions (use minimum observed)
    // Internal positions should have minimal damage
    float sum_internal_5 = 0.0f;
    float sum_internal_3 = 0.0f;
    int count_internal_5 = 0;
    int count_internal_3 = 0;
    for (size_t i = 10; i < 15; ++i) {
        if (profile.t_freq_5prime[i] + profile.c_freq_5prime[i] > 0) {
            sum_internal_5 += profile.t_freq_5prime[i] / (profile.t_freq_5prime[i] + profile.c_freq_5prime[i]);
            count_internal_5++;
        }
        if (profile.a_freq_3prime[i] + profile.g_freq_3prime[i] > 0) {
            sum_internal_3 += profile.a_freq_3prime[i] / (profile.a_freq_3prime[i] + profile.g_freq_3prime[i]);
            count_internal_3++;
        }
    }
    float avg_internal_5 = (count_internal_5 > 0) ? (sum_internal_5 / count_internal_5) : 0.0f;
    float avg_internal_3 = (count_internal_3 > 0) ? (sum_internal_3 / count_internal_3) : 0.0f;

    // Use a conservative floor to avoid overconfidence when no internal data
    const float background_floor = 0.002f;
    delta_background_ = std::max(std::min(avg_internal_5, avg_internal_3), background_floor);

    // Clear profile cache so new profiles reflect updated parameters
    std::lock_guard<std::mutex> lock(cache_mutex_);
    profile_cache_ = {};
}

DamageProfile DamageModel::create_profile(size_t seq_len) const {
    DamageProfile profile;
    profile.lambda_5prime = lambda_5prime_;
    profile.lambda_3prime = lambda_3prime_;
    profile.delta_max = delta_max_;
    profile.delta_background = delta_background_;

    // Precompute damage probabilities for each position
    profile.ct_prob_5prime.resize(seq_len);
    profile.ct_prob_3prime.resize(seq_len);
    profile.ga_prob_5prime.resize(seq_len);
    profile.ga_prob_3prime.resize(seq_len);

    for (size_t pos = 0; pos < seq_len; ++pos) {
        // C->T from 5' end
        float dist_5prime = static_cast<float>(pos);
        profile.ct_prob_5prime[pos] = delta_max_ * std::exp(-lambda_5prime_ * dist_5prime)
                                    + delta_background_;

        // G->A from 3' end (reverse complement of C->T)
        float dist_3prime = static_cast<float>(seq_len - pos - 1);
        profile.ga_prob_3prime[pos] = delta_max_ * std::exp(-lambda_3prime_ * dist_3prime)
                                    + delta_background_;
    }

    // Use empirical data if available
    if (empirical_ct_5prime_) {
        for (size_t i = 0; i < std::min(empirical_ct_5prime_->size(), seq_len); ++i) {
            profile.ct_prob_5prime[i] = (*empirical_ct_5prime_)[i];
        }
    }

    if (empirical_ga_3prime_) {
        for (size_t i = 0; i < std::min(empirical_ga_3prime_->size(), seq_len); ++i) {
            profile.ga_prob_3prime[seq_len - i - 1] = (*empirical_ga_3prime_)[i];
        }
    }

    return profile;
}

const DamageProfile& DamageModel::create_profile_cached(size_t seq_len) const {
    // Clamp to cache size
    if (seq_len >= profile_cache_.size()) {
        seq_len = profile_cache_.size() - 1;
    }

    // Check cache first (double-checked locking for performance)
    if (profile_cache_[seq_len].has_value()) {
        return *profile_cache_[seq_len];
    }

    // Lock and create profile if not cached
    std::lock_guard<std::mutex> lock(cache_mutex_);

    // Double-check after acquiring lock
    if (!profile_cache_[seq_len].has_value()) {
        profile_cache_[seq_len] = create_profile(seq_len);
    }

    return *profile_cache_[seq_len];
}

bool DamageModel::is_ancient_compatible() const {
    return delta_max_ > 0.1f && delta_max_ < 0.6f &&
           lambda_5prime_ > 0.01f && lambda_5prime_ < 1.0f;
}

}  // namespace agp
