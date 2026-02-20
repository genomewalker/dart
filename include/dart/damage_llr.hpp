#pragma once

/**
 * @file damage_llr.hpp
 * @brief Unified damage log-likelihood ratio computation
 *
 * Consolidates all LLR computations for ancient DNA damage:
 * - Base-level: P(obs_base | damaged) / P(obs_base | undamaged)
 * - Codon-level: P(obs_codon | damaged, is_stop) / P(obs_codon | undamaged)
 * - Hexamer-level: P(obs_hexamer | damaged) / P(obs_hexamer | undamaged)
 *
 * Uses domain-specific hexamer tables for proper priors.
 *
 * From OpenCode collaboration:
 * - STOP_LOG_ODDS = 1.0 threshold for stop codon correction
 * - NONSTOP_LOG_ODDS = 2.5 threshold for non-stop correction
 * - D_MIN_STOP = 0.005 minimum damage rate to consider
 */

#include "hexamer_tables.hpp"
#include "gtdb_dicodon_phase.hpp"
#include "fungi_dicodon_phase.hpp"
#include "plant_dicodon_phase.hpp"
#include "viral_dicodon_phase.hpp"
#include "protozoa_dicodon_phase.hpp"
#include "invertebrate_dicodon_phase.hpp"
#include "vertebrate_mammalian_dicodon_phase.hpp"
#include "vertebrate_other_dicodon_phase.hpp"
#include <cmath>
#include <array>
#include <string>
#include <algorithm>

namespace dart {

/**
 * Unified damage log-likelihood ratio computer
 *
 * Handles both C→T (5' end) and G→A (3' end) damage patterns with
 * proper domain-specific priors from hexamer tables.
 */
class DamageLLR {
public:
    // Damage parameters
    float d_max_5p = 0.0f;    // Maximum C→T rate at 5' end
    float d_max_3p = 0.0f;    // Maximum G→A rate at 3' end
    float lambda_5p = 0.3f;   // Decay rate from 5' end
    float lambda_3p = 0.3f;   // Decay rate from 3' end
    Domain domain = Domain::GTDB;

    // Correction thresholds (from OpenCode review)
    static constexpr float STOP_LOG_ODDS = 1.0f;    // ~73% posterior for stop correction
    static constexpr float NONSTOP_LOG_ODDS = 2.5f; // ~92% posterior for non-stop correction
    static constexpr float D_MIN_STOP = 0.005f;     // Minimum damage rate to consider

    DamageLLR() = default;

    DamageLLR(float d5, float d3, float l5, float l3, Domain dom = Domain::GTDB)
        : d_max_5p(d5), d_max_3p(d3), lambda_5p(l5), lambda_3p(l3), domain(dom) {}

    /**
     * Check if damage model is active (has meaningful damage rates)
     */
    bool is_active() const {
        return d_max_5p >= D_MIN_STOP || d_max_3p >= D_MIN_STOP;
    }

    // =========================================================================
    // Position-dependent damage rates
    // =========================================================================

    /**
     * Get C→T damage rate at position from 5' end
     */
    float ct_rate_at_5prime(size_t pos) const {
        return d_max_5p * std::exp(-lambda_5p * static_cast<float>(pos));
    }

    /**
     * Get G→A damage rate at position from 3' end
     */
    float ga_rate_at_3prime(size_t pos) const {
        return d_max_3p * std::exp(-lambda_3p * static_cast<float>(pos));
    }

    /**
     * Get C→T damage rate at absolute position for given sequence length
     */
    float ct_rate(size_t abs_pos, size_t seq_len, bool is_reverse) const {
        if (is_reverse) {
            // Reverse strand: C→T is from original 3' (now 5' after RC)
            return ga_rate_at_3prime(abs_pos);
        } else {
            // Forward strand: C→T is from 5'
            return ct_rate_at_5prime(abs_pos);
        }
    }

    /**
     * Get G→A damage rate at absolute position for given sequence length
     */
    float ga_rate(size_t abs_pos, size_t seq_len, bool is_reverse) const {
        size_t dist_from_3prime = seq_len - 1 - abs_pos;
        if (is_reverse) {
            // Reverse strand: G→A is from original 5' (now 3' after RC)
            return ct_rate_at_5prime(dist_from_3prime);
        } else {
            // Forward strand: G→A is from 3'
            return ga_rate_at_3prime(dist_from_3prime);
        }
    }

    // =========================================================================
    // Base-level LLR
    // =========================================================================

    /**
     * Compute log-likelihood ratio for observed base at position
     *
     * LLR = log[ P(obs | damaged) / P(obs | undamaged) ]
     *
     * For 5' C→T damage:
     *   - Observe T: LLR = log[ (pT + pC*D) / pT ]
     *   - Observe C: LLR = log[ (1-D) ]  (C survived, evidence against damage)
     *   - Other bases: LLR = 0
     *
     * @param obs Observed base (a/c/g/t)
     * @param abs_pos Position in sequence
     * @param seq_len Total sequence length
     * @param is_reverse Whether analyzing reverse strand
     * @param phase Reading frame phase (for hexamer prior)
     * @return Log-likelihood ratio (positive = evidence for damage)
     */
    float base_llr(char obs, size_t abs_pos, size_t seq_len, bool is_reverse, int phase = 0) const {
        if (!is_active()) return 0.0f;

        char b = obs | 0x20;  // lowercase

        // Get damage rates
        float d_ct = ct_rate(abs_pos, seq_len, is_reverse);
        float d_ga = ga_rate(abs_pos, seq_len, is_reverse);

        // Get base priors from domain-specific tables
        float pT = 0.25f, pC = 0.25f, pA = 0.25f, pG = 0.25f;
        // TODO: Could use hexamer context for more accurate priors

        float llr = 0.0f;

        // 5' end: C→T damage
        if (abs_pos < 15) {
            if (b == 't') {
                // T observed: could be original T or damaged C
                float p_obs_damaged = pT + pC * d_ct;
                float p_obs_undamaged = pT;
                llr += std::log(p_obs_damaged / std::max(p_obs_undamaged, 1e-6f));
            } else if (b == 'c') {
                // C observed: survived damage = evidence against damage
                llr += std::log(std::max(1.0f - d_ct, 1e-6f));
            }
        }

        // 3' end: G→A damage
        size_t dist_from_3p = seq_len - 1 - abs_pos;
        if (dist_from_3p < 15) {
            if (b == 'a') {
                // A observed: could be original A or damaged G
                float p_obs_damaged = pA + pG * d_ga;
                float p_obs_undamaged = pA;
                llr += std::log(p_obs_damaged / std::max(p_obs_undamaged, 1e-6f));
            } else if (b == 'g') {
                // G observed: survived damage = evidence against damage
                llr += std::log(std::max(1.0f - d_ga, 1e-6f));
            }
        }

        return llr;
    }

    // =========================================================================
    // Codon-level LLR
    // =========================================================================

    /**
     * Compute log-odds that a stop codon is damage-induced
     *
     * Uses Bayesian reasoning:
     * - Prior: P(real_stop) ~ 0.01 for internal codons (stop=0.01 vs coding=0.99)
     * - Likelihood: P(obs_stop | damage, original_was_coding)
     *
     * Damage-induced stop codons:
     * - TAA from CAA (C→T at pos 1): P(stop|damage) = d_ct
     * - TAG from CAG (C→T at pos 1): P(stop|damage) = d_ct
     * - TGA from CGA (C→T at pos 1): P(stop|damage) = d_ct
     * - TGA from TGG (G→A at pos 3): P(stop|damage) = d_ga
     *
     * @param codon 3-character codon string
     * @param codon_pos Position of first nucleotide in sequence
     * @param seq_len Total sequence length
     * @param is_reverse Whether analyzing reverse strand
     * @return Log-odds that stop is damage-induced (positive = likely damage)
     */
    float stop_codon_log_odds(const char* codon, size_t codon_pos, size_t seq_len, bool is_reverse) const {
        if (!is_active()) return -10.0f;  // No damage, stop is real

        char c0 = codon[0] | 0x20;
        char c1 = codon[1] | 0x20;
        char c2 = codon[2] | 0x20;

        // Check if this is a stop codon
        bool is_stop = (c0 == 't' && c1 == 'a' && c2 == 'a') ||  // TAA
                       (c0 == 't' && c1 == 'a' && c2 == 'g') ||  // TAG
                       (c0 == 't' && c1 == 'g' && c2 == 'a');    // TGA

        if (!is_stop) return -10.0f;

        // Get damage rates at codon positions
        float d_ct_pos1 = ct_rate(codon_pos, seq_len, is_reverse);
        float d_ga_pos3 = ga_rate(codon_pos + 2, seq_len, is_reverse);

        // Prior: P(internal stop) is very low in coding sequences
        constexpr float prior_stop = 0.01f;
        float log_prior_odds = std::log(prior_stop / (1.0f - prior_stop));

        // Likelihood ratio for damage
        float damage_prob = 0.0f;

        // TAA from CAA, TAG from CAG: C→T at position 1
        if ((c0 == 't' && c1 == 'a' && c2 == 'a') ||
            (c0 == 't' && c1 == 'a' && c2 == 'g')) {
            damage_prob = d_ct_pos1;
        }

        // TGA can come from CGA (C→T at pos 1) or TGG (G→A at pos 3)
        if (c0 == 't' && c1 == 'g' && c2 == 'a') {
            damage_prob = std::max(d_ct_pos1, d_ga_pos3);
        }

        if (damage_prob < D_MIN_STOP) {
            return -10.0f;  // Too far from terminal, stop is real
        }

        // Log-odds that this is damage-induced
        // P(damage-induced | obs_stop) / P(real_stop | obs_stop)
        // = P(obs_stop | damage) * P(damage) / P(obs_stop | real) * P(real)
        // Approximation: log_odds = log(damage_prob / prior_stop)
        float log_likelihood_ratio = std::log(damage_prob / std::max(prior_stop, 1e-6f));

        return log_likelihood_ratio;
    }

    /**
     * Check if a codon correction should be made
     *
     * @param codon 3-character codon string
     * @param codon_pos Position of first nucleotide
     * @param seq_len Total sequence length
     * @param is_reverse Whether analyzing reverse strand
     * @return True if codon should be corrected
     */
    bool should_correct_codon(const char* codon, size_t codon_pos, size_t seq_len, bool is_reverse) const {
        float log_odds = stop_codon_log_odds(codon, codon_pos, seq_len, is_reverse);
        return log_odds >= STOP_LOG_ODDS;
    }

    // =========================================================================
    // Hexamer-level LLR
    // =========================================================================

    /**
     * Compute damage LLR for a hexamer at a position
     *
     * Uses domain-specific hexamer tables to compute:
     * LLR = log[ P(obs_hex | damaged) / P(obs_hex | undamaged) ]
     *
     * For 5' hexamers: considers C→T damage at each T position
     * For 3' hexamers: considers G→A damage at each A position
     *
     * @param hex_code Encoded hexamer (0-4095)
     * @param start_pos Position of hexamer start in sequence
     * @param seq_len Total sequence length
     * @param phase Reading frame phase (0, 1, 2)
     * @return Log-likelihood ratio for damage
     */
    float hexamer_llr(uint32_t hex_code, size_t start_pos, size_t seq_len, int phase = 0) const {
        if (!is_active() || hex_code >= 4096) return 0.0f;

        float llr = 0.0f;

        // Check if in 5' or 3' terminal region
        bool is_5prime = (start_pos < 15);
        size_t end_pos = start_pos + 5;
        bool is_3prime = (seq_len - end_pos < 15);

        if (!is_5prime && !is_3prime) return 0.0f;

        // Get hexamer frequencies
        float p_undamaged = get_phase_freq(hex_code, phase);

        // Compute P(observed | damaged) by summing over possible source hexamers
        float p_damaged = compute_damaged_hexamer_prob(hex_code, start_pos, seq_len, is_5prime);

        if (p_undamaged > 1e-10f && p_damaged > 1e-10f) {
            llr = std::log(p_damaged / p_undamaged);
        }

        return llr;
    }

    /**
     * Compute total damage LLR for a sequence by scanning hexamers
     *
     * @param seq DNA sequence
     * @param terminal_window Number of positions to scan at each end
     * @return Total log-likelihood ratio
     */
    float sequence_llr(const std::string& seq, int terminal_window = 15) const {
        if (!is_active() || seq.length() < 12) return 0.0f;

        float total_llr = 0.0f;
        const size_t len = seq.length();

        // 5' terminal hexamers
        for (size_t i = 0; i + 5 < len && i < static_cast<size_t>(terminal_window); i += 3) {
            uint32_t code = encode_hexamer(seq.c_str() + i);
            if (code < 4096) {
                total_llr += hexamer_llr(code, i, len, i % 3);
            }
        }

        // 3' terminal hexamers
        for (size_t dist = 0; dist + 5 < len && dist < static_cast<size_t>(terminal_window); dist += 3) {
            size_t pos = len - 6 - dist;
            if (pos < len) {
                uint32_t code = encode_hexamer(seq.c_str() + pos);
                if (code < 4096) {
                    total_llr += hexamer_llr(code, pos, len, pos % 3);
                }
            }
        }

        return total_llr;
    }

private:
    /**
     * Get domain-specific phase hexamer frequency
     */
    float get_phase_freq(uint32_t code, int phase) const {
        if (code >= 4096) return 1.0f / 4096.0f;

        switch (domain) {
            case Domain::FUNGI:
                switch (phase) {
                    case 0: return FUNGI_PHASE0_HEXAMER_FREQ[code];
                    case 1: return FUNGI_PHASE1_HEXAMER_FREQ[code];
                    case 2: return FUNGI_PHASE2_HEXAMER_FREQ[code];
                }
                break;
            case Domain::PLANT:
                switch (phase) {
                    case 0: return PLANT_PHASE0_HEXAMER_FREQ[code];
                    case 1: return PLANT_PHASE1_HEXAMER_FREQ[code];
                    case 2: return PLANT_PHASE2_HEXAMER_FREQ[code];
                }
                break;
            case Domain::VIRAL:
                switch (phase) {
                    case 0: return VIRAL_PHASE0_HEXAMER_FREQ[code];
                    case 1: return VIRAL_PHASE1_HEXAMER_FREQ[code];
                    case 2: return VIRAL_PHASE2_HEXAMER_FREQ[code];
                }
                break;
            case Domain::PROTOZOA:
                switch (phase) {
                    case 0: return PROTOZOA_PHASE0_HEXAMER_FREQ[code];
                    case 1: return PROTOZOA_PHASE1_HEXAMER_FREQ[code];
                    case 2: return PROTOZOA_PHASE2_HEXAMER_FREQ[code];
                }
                break;
            case Domain::INVERTEBRATE:
                switch (phase) {
                    case 0: return INVERTEBRATE_PHASE0_HEXAMER_FREQ[code];
                    case 1: return INVERTEBRATE_PHASE1_HEXAMER_FREQ[code];
                    case 2: return INVERTEBRATE_PHASE2_HEXAMER_FREQ[code];
                }
                break;
            case Domain::VERTEBRATE_MAMMALIAN:
                switch (phase) {
                    case 0: return VERTEBRATE_MAMMALIAN_PHASE0_HEXAMER_FREQ[code];
                    case 1: return VERTEBRATE_MAMMALIAN_PHASE1_HEXAMER_FREQ[code];
                    case 2: return VERTEBRATE_MAMMALIAN_PHASE2_HEXAMER_FREQ[code];
                }
                break;
            case Domain::VERTEBRATE_OTHER:
                switch (phase) {
                    case 0: return VERTEBRATE_OTHER_PHASE0_HEXAMER_FREQ[code];
                    case 1: return VERTEBRATE_OTHER_PHASE1_HEXAMER_FREQ[code];
                    case 2: return VERTEBRATE_OTHER_PHASE2_HEXAMER_FREQ[code];
                }
                break;
            case Domain::GTDB:
            case Domain::META:
            default:
                switch (phase) {
                    case 0: return GTDB_PHASE0_HEXAMER_FREQ[code];
                    case 1: return GTDB_PHASE1_HEXAMER_FREQ[code];
                    case 2: return GTDB_PHASE2_HEXAMER_FREQ[code];
                }
                break;
        }
        return 1.0f / 4096.0f;
    }

    /**
     * Compute P(observed hexamer | damaged) by marginalizing over source hexamers
     *
     * For 5' damage (C→T): each T could have come from C
     * For 3' damage (G→A): each A could have come from G
     */
    float compute_damaged_hexamer_prob(uint32_t obs_code, size_t start_pos, size_t seq_len, bool is_5prime) const {
        if (obs_code >= 4096) return 0.0f;

        // Decode observed hexamer
        char obs_hex[7];
        for (int i = 5; i >= 0; --i) {
            static const char bases[] = "ACGT";
            obs_hex[i] = bases[obs_code & 3];
            obs_code >>= 2;
        }
        obs_hex[6] = '\0';

        // Restore obs_code (we destroyed it above)
        obs_code = 0;
        for (int i = 0; i < 6; ++i) {
            obs_code <<= 2;
            switch (obs_hex[i]) {
                case 'A': break;
                case 'C': obs_code |= 1; break;
                case 'G': obs_code |= 2; break;
                case 'T': obs_code |= 3; break;
            }
        }

        float total_prob = 0.0f;

        if (is_5prime) {
            // 5' C→T damage: find T positions that could be from C
            int t_positions[6];
            int n_t = 0;
            for (int i = 0; i < 6; ++i) {
                if (obs_hex[i] == 'T') t_positions[n_t++] = i;
            }

            // Enumerate all 2^n_t combinations
            for (int mask = 0; mask < (1 << n_t); ++mask) {
                char src_hex[7];
                std::copy(obs_hex, obs_hex + 7, src_hex);

                float trans_prob = 1.0f;
                for (int j = 0; j < n_t; ++j) {
                    int pos = t_positions[j];
                    float d = ct_rate_at_5prime(start_pos + pos);

                    if (mask & (1 << j)) {
                        // This T was originally C
                        src_hex[pos] = 'C';
                        trans_prob *= d;
                    }
                    // Else: T stays T (prob = 1, no change)
                }

                // C positions must survive
                for (int i = 0; i < 6; ++i) {
                    if (obs_hex[i] == 'C') {
                        float d = ct_rate_at_5prime(start_pos + i);
                        trans_prob *= (1.0f - d);
                    }
                }

                // Get source hexamer frequency
                uint32_t src_code = encode_hexamer(src_hex);
                if (src_code < 4096) {
                    float src_freq = get_phase_freq(src_code, 0);
                    total_prob += src_freq * trans_prob;
                }
            }
        } else {
            // 3' G→A damage: find A positions that could be from G
            size_t dist_from_end = seq_len - start_pos - 6;

            int a_positions[6];
            int n_a = 0;
            for (int i = 0; i < 6; ++i) {
                if (obs_hex[i] == 'A') a_positions[n_a++] = i;
            }

            for (int mask = 0; mask < (1 << n_a); ++mask) {
                char src_hex[7];
                std::copy(obs_hex, obs_hex + 7, src_hex);

                float trans_prob = 1.0f;
                for (int j = 0; j < n_a; ++j) {
                    int pos = a_positions[j];
                    float d = ga_rate_at_3prime(dist_from_end + (5 - pos));

                    if (mask & (1 << j)) {
                        // This A was originally G
                        src_hex[pos] = 'G';
                        trans_prob *= d;
                    }
                }

                // G positions must survive
                for (int i = 0; i < 6; ++i) {
                    if (obs_hex[i] == 'G') {
                        float d = ga_rate_at_3prime(dist_from_end + (5 - i));
                        trans_prob *= (1.0f - d);
                    }
                }

                uint32_t src_code = encode_hexamer(src_hex);
                if (src_code < 4096) {
                    float src_freq = get_phase_freq(src_code, 0);
                    total_prob += src_freq * trans_prob;
                }
            }
        }

        return total_prob;
    }
};

}  // namespace dart
