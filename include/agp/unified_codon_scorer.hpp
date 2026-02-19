/**
 * Unified Codon-Marginalization Frame Scorer for Ancient DNA
 *
 * Joint probabilistic model for frame, strand, and damage correction:
 * P(frame, strand | seq, damage) via codon-level marginalization
 *
 * Scoring components:
 * 1. Codon marginal likelihood (damage-aware)
 * 2. X circular code enrichment (frame-specific)
 * 3. AA bigram directional bias (strand-specific)
 * 4. Stop codon positional explainability (damage-aware)
 */

#pragma once

#include <string>
#include <array>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <string_view>
#include "types.hpp"
#include "codon_tables.hpp"
#include "dicodon_phase.hpp"
#include "strand_scoring.hpp"

namespace agp {
namespace codon {

// ============================================================================
// Thread-local buffers for zero-allocation hot path
// ============================================================================

// Maximum sequence length we support without allocation
constexpr size_t MAX_SEQ_LEN = 2048;
constexpr size_t MAX_CODONS = MAX_SEQ_LEN / 3;

// Thread-local buffers to avoid heap allocations in hot path
struct alignas(64) ThreadLocalBuffers {
    char rc_buffer[MAX_SEQ_LEN];           // Reverse complement
    char protein_buffer[6][MAX_CODONS];    // Proteins for 6 hypotheses
    char corrected_buffer[6][MAX_CODONS];  // Corrected proteins
    float codon_scores[6][MAX_CODONS];     // Per-codon log marginals
    float stop_probs[6][MAX_CODONS];       // Per-codon stop probabilities
    size_t protein_lens[6];                // Actual lengths
    bool initialized = false;
};

inline ThreadLocalBuffers& get_thread_buffers() {
    static thread_local ThreadLocalBuffers buffers;
    return buffers;
}

// Fast reverse complement into pre-allocated buffer
// Returns string_view into thread-local buffer
inline std::string_view revcomp_fast(const char* seq, size_t len) {
    auto& buf = get_thread_buffers();
    if (len > MAX_SEQ_LEN) {
        // Fallback for very long sequences (rare)
        static thread_local std::string fallback;
        fallback.resize(len);
        for (size_t i = 0; i < len; ++i) {
            char c = seq[len - 1 - i];
            switch (c) {
                case 'A': fallback[i] = 'T'; break;
                case 'T': fallback[i] = 'A'; break;
                case 'C': fallback[i] = 'G'; break;
                case 'G': fallback[i] = 'C'; break;
                case 'a': fallback[i] = 't'; break;
                case 't': fallback[i] = 'a'; break;
                case 'c': fallback[i] = 'g'; break;
                case 'g': fallback[i] = 'c'; break;
                default:  fallback[i] = 'N'; break;
            }
        }
        return std::string_view(fallback);
    }

    // Fast path: use thread-local buffer
    for (size_t i = 0; i < len; ++i) {
        char c = seq[len - 1 - i];
        switch (c) {
            case 'A': buf.rc_buffer[i] = 'T'; break;
            case 'T': buf.rc_buffer[i] = 'A'; break;
            case 'C': buf.rc_buffer[i] = 'G'; break;
            case 'G': buf.rc_buffer[i] = 'C'; break;
            case 'a': buf.rc_buffer[i] = 't'; break;
            case 't': buf.rc_buffer[i] = 'a'; break;
            case 'c': buf.rc_buffer[i] = 'g'; break;
            case 'g': buf.rc_buffer[i] = 'c'; break;
            default:  buf.rc_buffer[i] = 'N'; break;
        }
    }
    return std::string_view(buf.rc_buffer, len);
}

// Damage parameters (physical read coordinates)
struct DamageParams {
    float d_max_5p = 0.0f;
    float d_max_3p = 0.0f;
    float lambda_5p = 0.3f;
    float lambda_3p = 0.3f;

    // CpG context boost: methylated CpG → C→T is elevated
    // Typical boost is 2-5x in ancient DNA
    float cpg_boost = 5.0f;

    bool has_damage() const { return d_max_5p > 0.01f || d_max_3p > 0.01f; }

    float delta_CT(int phys_pos) const {
        return d_max_5p * std::exp(-lambda_5p * static_cast<float>(phys_pos));
    }

    // CpG-aware C→T rate: boosted when next base is G
    float delta_CT_context(int phys_pos, bool next_is_G) const {
        float base_rate = delta_CT(phys_pos);
        if (next_is_G) {
            // CpG context: boost the damage rate
            return std::min(1.0f, base_rate * cpg_boost);
        }
        return base_rate;
    }

    float delta_GA(int phys_pos, int L) const {
        int dist_from_3p = L - 1 - phys_pos;
        return d_max_3p * std::exp(-lambda_3p * static_cast<float>(dist_from_3p));
    }

    // CpG-aware G→A rate: boosted when previous base is C (CpG on opposite strand)
    float delta_GA_context(int phys_pos, int L, bool prev_is_C) const {
        float base_rate = delta_GA(phys_pos, L);
        if (prev_is_C) {
            // CpG context on opposite strand: boost the damage rate
            return std::min(1.0f, base_rate * cpg_boost);
        }
        return base_rate;
    }
};

// ============================================================================
// SELF-TRAINED PERIODIC MARKOV MODEL (EM-learned from sample)
// ============================================================================

/**
 * 3-periodic Markov model with order-2 (dinucleotide) context.
 * Learned via EM from sample reads to adapt to metagenome composition.
 *
 * Parameters: θ[phase][context][base] = P(base | context, phase)
 *   - phase ∈ {0,1,2} = codon position
 *   - context = 16 dinucleotides (4^2)
 *   - base = 4 nucleotides
 *
 * Total: 3 × 16 × 4 = 192 parameters
 */
struct PeriodicMarkovModel {
    // θ[phase][context][base] - probabilities
    std::array<std::array<std::array<float, 4>, 16>, 3> theta;

    // Counts for accumulation during M-step
    std::array<std::array<std::array<float, 4>, 16>, 3> counts;
    std::array<std::array<float, 16>, 3> totals;

    bool trained = false;
    float alpha = 1.0f;  // Dirichlet smoothing parameter

    PeriodicMarkovModel() {
        reset();
    }

    // Local base index function (A=0, C=1, G=2, T=3)
    static int base_idx_local(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    }

    // Local complement function
    static char complement_local(char c) {
        switch (c) {
            case 'A': return 'T'; case 'T': return 'A';
            case 'C': return 'G'; case 'G': return 'C';
            case 'a': return 't'; case 't': return 'a';
            case 'c': return 'g'; case 'g': return 'c';
            default: return 'N';
        }
    }

    // Local revcomp function
    static std::string revcomp_local(const std::string& seq) {
        std::string rc;
        rc.resize(seq.size());
        for (size_t i = 0; i < seq.size(); ++i) {
            rc[i] = complement_local(seq[seq.size() - 1 - i]);
        }
        return rc;
    }

    void reset() {
        // Initialize with uniform distribution
        for (int p = 0; p < 3; ++p) {
            for (int c = 0; c < 16; ++c) {
                for (int b = 0; b < 4; ++b) {
                    theta[p][c][b] = 0.25f;
                    counts[p][c][b] = 0.0f;
                }
                totals[p][c] = 0.0f;
            }
        }
        trained = false;
    }

    // Convert dinucleotide to context index (0-15)
    static int ctx_index(char c1, char c2) {
        int i1 = base_idx_local(c1);
        int i2 = base_idx_local(c2);
        if (i1 < 0 || i2 < 0) return -1;
        return i1 * 4 + i2;
    }

    // Accumulate weighted counts from a read (M-step contribution)
    void accumulate(const std::string& seq, int frame, bool forward, float weight) {
        std::string oriented = forward ? seq : revcomp_local(seq);
        int L = static_cast<int>(oriented.size());

        for (int j = 2; j < L; ++j) {
            int phase = (j - frame + 300) % 3;
            int ctx = ctx_index(oriented[j-2], oriented[j-1]);
            int base = base_idx_local(oriented[j]);

            if (ctx >= 0 && base >= 0) {
                counts[phase][ctx][base] += weight;
                totals[phase][ctx] += weight;
            }
        }
    }

    // Normalize counts to probabilities (end of M-step)
    void normalize() {
        for (int p = 0; p < 3; ++p) {
            for (int c = 0; c < 16; ++c) {
                float total = totals[p][c] + 4.0f * alpha;
                for (int b = 0; b < 4; ++b) {
                    theta[p][c][b] = (counts[p][c][b] + alpha) / total;
                }
            }
        }
        trained = true;
    }

    // Score a sequence under this frame hypothesis
    float score(const std::string& seq, int frame, bool forward) const {
        if (!trained) return 0.0f;

        std::string oriented = forward ? seq : revcomp_local(seq);
        int L = static_cast<int>(oriented.size());

        float ll = 0.0f;
        int count = 0;

        for (int j = 2; j < L; ++j) {
            int phase = (j - frame + 300) % 3;
            int ctx = ctx_index(oriented[j-2], oriented[j-1]);
            int base = base_idx_local(oriented[j]);

            if (ctx >= 0 && base >= 0) {
                ll += std::log(theta[phase][ctx][base] + 1e-10f);
                count++;
            }
        }

        // Normalize by length for comparability
        return count > 0 ? ll / count : 0.0f;
    }
};

// Global model instance (set during Pass 1, used during Pass 2)
inline PeriodicMarkovModel& get_learned_periodic_model() {
    static PeriodicMarkovModel model;
    return model;
}

// ============================================================================
// GC-conditional periodic Markov model
// ============================================================================

/**
 * GC-conditional periodic Markov model with multiple GC bins.
 *
 * Learns separate transition parameters for different GC content ranges,
 * with shrinkage to a global model for sparse bins.
 */
struct GCConditionalPeriodicModel {
    static constexpr int N_BINS = 4;  // <40%, 40-50%, 50-60%, >60%
    static constexpr float BIN_CENTERS[4] = {0.30f, 0.45f, 0.55f, 0.70f};
    static constexpr float BIN_WIDTH = 0.15f;

    std::array<PeriodicMarkovModel, N_BINS> bin_models;
    PeriodicMarkovModel global_model;
    std::array<int, N_BINS> bin_counts;
    bool trained = false;

    // Compute GC from interior (avoid terminal damage bias)
    static float compute_interior_gc(const std::string& seq, int skip = 10) {
        int gc = 0, total = 0;
        int L = static_cast<int>(seq.size());
        int start = std::min(skip, L / 3);
        int end = std::max(L - skip, 2 * L / 3);

        for (int i = start; i < end; ++i) {
            char c = seq[i];
            if (c == 'G' || c == 'C' || c == 'g' || c == 'c') gc++;
            if (c == 'A' || c == 'T' || c == 'G' || c == 'C' ||
                c == 'a' || c == 't' || c == 'g' || c == 'c') total++;
        }
        return total > 0 ? static_cast<float>(gc) / total : 0.5f;
    }

    // Get hard bin index (0-3)
    static int get_gc_bin(float gc) {
        if (gc < 0.40f) return 0;
        if (gc < 0.50f) return 1;
        if (gc < 0.60f) return 2;
        return 3;
    }

    // Get soft bin weights (triangular kernel)
    static void get_soft_weights(float gc, float* weights) {
        float sum = 0.0f;
        for (int b = 0; b < N_BINS; ++b) {
            float dist = std::abs(gc - BIN_CENTERS[b]);
            weights[b] = std::max(0.0f, 1.0f - dist / BIN_WIDTH);
            sum += weights[b];
        }
        // Normalize
        if (sum > 0.0f) {
            for (int b = 0; b < N_BINS; ++b) weights[b] /= sum;
        } else {
            // Fallback to hard assignment
            int bin = get_gc_bin(gc);
            for (int b = 0; b < N_BINS; ++b) weights[b] = (b == bin) ? 1.0f : 0.0f;
        }
    }

    void reset() {
        for (int b = 0; b < N_BINS; ++b) {
            bin_models[b].reset();
            bin_counts[b] = 0;
        }
        global_model.reset();
        trained = false;
    }

    // Accumulate with soft binning
    void accumulate(const std::string& seq, int frame, bool forward, float weight) {
        float gc = compute_interior_gc(seq);
        float bin_weights[N_BINS];
        get_soft_weights(gc, bin_weights);

        for (int b = 0; b < N_BINS; ++b) {
            if (bin_weights[b] > 0.01f) {
                bin_models[b].accumulate(seq, frame, forward, weight * bin_weights[b]);
            }
        }
        global_model.accumulate(seq, frame, forward, weight);

        // Track primary bin for count
        int primary_bin = get_gc_bin(gc);
        bin_counts[primary_bin]++;
    }

    // Normalize with shrinkage to global model
    void normalize() {
        global_model.normalize();

        for (int b = 0; b < N_BINS; ++b) {
            bin_models[b].normalize();

            // Shrinkage: blend with global based on sample size
            // Small bins use more global prior
            float lambda = 1.0f / (1.0f + bin_counts[b] / 1000.0f);

            for (int p = 0; p < 3; ++p) {
                for (int c = 0; c < 16; ++c) {
                    for (int base = 0; base < 4; ++base) {
                        bin_models[b].theta[p][c][base] =
                            (1.0f - lambda) * bin_models[b].theta[p][c][base] +
                            lambda * global_model.theta[p][c][base];
                    }
                }
            }
            bin_models[b].trained = true;
        }
        trained = true;
    }

    // Score with marginalization over GC bins
    float score(const std::string& seq, int frame, bool forward) const {
        if (!trained) return 0.0f;

        float gc = compute_interior_gc(seq);
        float bin_weights[N_BINS];
        get_soft_weights(gc, bin_weights);

        // Log-sum-exp for marginalization
        float log_scores[N_BINS];
        float max_log = -1e30f;

        for (int b = 0; b < N_BINS; ++b) {
            if (bin_weights[b] > 0.01f) {
                log_scores[b] = std::log(bin_weights[b]) +
                                bin_models[b].score(seq, frame, forward) *
                                static_cast<float>(seq.size());  // Undo normalization
            } else {
                log_scores[b] = -1e30f;
            }
            if (log_scores[b] > max_log) max_log = log_scores[b];
        }

        // Compute log-sum-exp
        float sum = 0.0f;
        for (int b = 0; b < N_BINS; ++b) {
            if (log_scores[b] > -1e20f) {
                sum += std::exp(log_scores[b] - max_log);
            }
        }

        return (max_log + std::log(sum + 1e-30f)) / static_cast<float>(seq.size());
    }
};

// Global GC-conditional model
inline GCConditionalPeriodicModel& get_gc_conditional_model() {
    static GCConditionalPeriodicModel model;
    return model;
}

// Scoring weights
// After empirical testing:
// - Base codon marginalization: main signal (2x random baseline)
// - Strand-only LLR: NEW - uses all positions for strand discrimination
// - Self-trained periodic: EM-learned from sample
// - X code: DISABLED - hurts on metagenome (designed for single organisms)
// - AA bigram: DISABLED - partial table was not robust
// - Stop positional: marginal benefit, keep at low weight
// IMPORTANT: Strand discrimination is FUNDAMENTALLY LIMITED for double-stranded
// ancient DNA libraries without reference alignment:
// - Hexamer patterns are symmetric (both strands look equally coding)
// - Damage patterns are symmetric (both strands show C→T at 5', G→A at 3')
// For metagenomics, use --both-strands and let MMseqs2 determine strand via alignment.
// Frame selection (86% accuracy) works well; strand selection (~50%) is random.
struct ScoringWeights {
    float w_codon = 1.0f;      // Base codon marginal (stop codon penalty)
    float w_strand_only = 0.0f;// Strand-only LLR - DISABLED (doesn't help for DS libraries)
    float w_strand_hex = 0.0f; // Strand hexamer LLR - DISABLED (cannot distinguish gene strand)
    float w_periodic = 1.0f;   // 3-periodic hexamer (dicodon phase) - PRIMARY for frame
    float w_self_trained = 0.0f; // Self-trained periodic
    float w_gc_conditional = 0.0f; // GC-conditional
    float w_fourier = 1.0f;    // RNY bias score - KEY IMPROVEMENT for frame
    float w_xcode = 0.0f;      // X circular code
    float w_bigram = 0.0f;     // AA bigram
    float w_stoppos = 1.0f;    // Stop positional - KEY IMPROVEMENT for frame
};

// Frame+strand result
struct FrameStrandResult {
    int frame;
    bool forward;
    float log_likelihood;
    float posterior;
    float score_codon;
    float score_strand_only;
    float score_strand_hex;   // Strand hexamer LLR
    float score_periodic;
    float score_self_trained;
    float score_fourier;
    float score_xcode;
    float score_bigram;
    float score_stoppos;
    std::string protein;
    std::string corrected_protein;
};

// ============================================================================
// Sample-estimated base priors (from interior positions to avoid damage bias)
// ============================================================================

struct BasePriors {
    std::array<float, 4> pi;  // A=0, C=1, G=2, T=3 (ACGT order)

    BasePriors() : pi{{0.25f, 0.25f, 0.25f, 0.25f}} {}

    BasePriors(float pA, float pC, float pG, float pT) {
        float sum = pA + pC + pG + pT;
        pi[0] = pA / sum;
        pi[1] = pC / sum;
        pi[2] = pG / sum;
        pi[3] = pT / sum;
    }

    float operator[](char base) const {
        switch (base) {
            case 'A': case 'a': return pi[0];
            case 'C': case 'c': return pi[1];
            case 'G': case 'g': return pi[2];
            case 'T': case 't': return pi[3];
            default: return 0.25f;
        }
    }

    static int base_idx(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    }
};

// ============================================================================
// Base encoding utilities
// ============================================================================

// Fast base-to-index lookup table (TCAG order: T=0, C=1, A=2, G=3)
alignas(64) constexpr int8_t BASE_TO_TCAG[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 0-15
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 16-31
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 32-47
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 48-63
    -1, 2,-1, 1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1, // 64-79  (A=65->2, C=67->1, G=71->3)
    -1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 80-95  (T=84->0)
    -1, 2,-1, 1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1, // 96-111 (a=97->2, c=99->1, g=103->3)
    -1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 112-127 (t=116->0)
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
};

inline int base_to_tcag(char c) {
    return BASE_TO_TCAG[static_cast<unsigned char>(c)];
}

inline char tcag_to_base(int idx) {
    constexpr char bases[] = "TCAG";
    return (idx >= 0 && idx < 4) ? bases[idx] : 'N';
}

// Fast codon to index - inline the lookup
inline int codon_to_idx(char c0, char c1, char c2) {
    int i0 = BASE_TO_TCAG[static_cast<unsigned char>(c0)];
    int i1 = BASE_TO_TCAG[static_cast<unsigned char>(c1)];
    int i2 = BASE_TO_TCAG[static_cast<unsigned char>(c2)];
    // Use bitwise OR to check all at once (if any is -1, result has high bit set)
    if ((i0 | i1 | i2) < 0) return -1;
    return (i0 << 4) | (i1 << 2) | i2;
}

inline void idx_to_codon(int idx, char& c0, char& c1, char& c2) {
    c0 = tcag_to_base(idx / 16);
    c1 = tcag_to_base((idx / 4) % 4);
    c2 = tcag_to_base(idx % 4);
}

// Fast complement using lookup table
alignas(64) constexpr char COMPLEMENT_TABLE[256] = {
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', // 0-15
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', // 16-31
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', // 32-47
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', // 48-63
    'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N', // 64-79  (A=65->T, C=67->G, G=71->C)
    'N','N','N','N','A','N','N','N','N','N','N','N','N','N','N','N', // 80-95  (T=84->A)
    'N','t','N','g','N','N','N','c','N','N','N','N','N','N','N','N', // 96-111 (a=97->t, c=99->g, g=103->c)
    'N','N','N','N','a','N','N','N','N','N','N','N','N','N','N','N', // 112-127 (t=116->a)
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N', // 128+
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
};

inline char complement(char c) {
    return COMPLEMENT_TABLE[static_cast<unsigned char>(c)];
}

// ============================================================================
// Genetic code
// ============================================================================

inline char translate_codon_idx(int idx) {
    static constexpr char AA_TABLE[64] = {
        'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S',
        'Y', 'Y', '*', '*', 'C', 'C', '*', 'W',
        'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P',
        'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R',
        'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T',
        'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R',
        'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A',
        'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G',
    };
    return (idx >= 0 && idx < 64) ? AA_TABLE[idx] : 'X';
}

inline int aa_to_idx(char aa) {
    switch (aa) {
        case 'A': return 0;  case 'C': return 1;  case 'D': return 2;
        case 'E': return 3;  case 'F': return 4;  case 'G': return 5;
        case 'H': return 6;  case 'I': return 7;  case 'K': return 8;
        case 'L': return 9;  case 'M': return 10; case 'N': return 11;
        case 'P': return 12; case 'Q': return 13; case 'R': return 14;
        case 'S': return 15; case 'T': return 16; case 'V': return 17;
        case 'W': return 18; case 'Y': return 19;
        default: return -1;
    }
}

// ============================================================================
// X Circular Code (Arquès & Michel)
// 20 codons with frame-retrieval property, enriched in coding sequences
// ============================================================================

inline bool is_x_code_codon(int codon_idx) {
    // X = {AAC, AAT, ACC, ATC, ATT, CAG, CCT, CGA, CTA, CTC,
    //      CTG, GAA, GAC, GAG, GAT, GCC, GGC, GGT, GTA, GTC}
    // Convert to TCAG encoding indices
    static constexpr bool X_SET[64] = {
        // TTx TCx TAx TGx CTx CCx CAx CGx ATx ACx AAx AGx GTx GCx GAx GGx
        0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  // Txx (none)
        0,1,0,1, 0,0,0,0, 0,1,0,0, 0,0,1,0,  // Cxx: CCT(5), CTG(19), CTA(18), CTC(17), CAG(27), CGA(30)
        0,0,1,0, 1,0,0,0, 0,1,1,0, 0,0,0,0,  // Axx: AAC(41), AAT(40), ACC(37), ATC(33), ATT(32)
        0,0,0,0, 1,1,0,0, 0,0,1,1, 0,1,1,0,  // Gxx: GAA(58), GAC(57), GAG(59), GAT(56), GCC(53), GGC(61), GGT(60), GTA(50), GTC(49)
    };
    // Recompute based on TCAG encoding: idx = T*16 + C*4 + A*G
    // AAC = A(2)*16 + A(2)*4 + C(1) = 32+8+1 = 41
    // Let me build this properly
    static const bool X_CODON[64] = {
        // idx: base0*16 + base1*4 + base2, T=0,C=1,A=2,G=3
        false, false, false, false,  // 0-3:   TTT TTC TTA TTG
        false, false, false, false,  // 4-7:   TCT TCC TCA TCG
        false, false, false, false,  // 8-11:  TAT TAC TAA TAG
        false, false, false, false,  // 12-15: TGT TGC TGA TGG
        false, false, true,  true,   // 16-19: CTT CTC CTA CTG  (CTA=18, CTG=19)
        true,  false, false, false,  // 20-23: CCT CCC CCA CCG  (CCT=20)
        false, false, false, true,   // 24-27: CAT CAC CAA CAG  (CAG=27)
        false, false, true,  false,  // 28-31: CGT CGC CGA CGG  (CGA=30)
        true,  true,  false, false,  // 32-35: ATT ATC ATA ATG  (ATT=32, ATC=33)
        false, true,  false, false,  // 36-39: ACT ACC ACA ACG  (ACC=37)
        true,  true,  false, false,  // 40-43: AAT AAC AAA AAG  (AAT=40, AAC=41)
        false, false, false, false,  // 44-47: AGT AGC AGA AGG
        false, true,  true,  false,  // 48-51: GTT GTC GTA GTG  (GTC=49, GTA=50)
        false, true,  false, false,  // 52-55: GCT GCC GCA GCG  (GCC=53)
        true,  true,  true,  true,   // 56-59: GAT GAC GAA GAG  (all 4)
        true,  true,  false, false,  // 60-63: GGT GGC GGA GGG  (GGT=60, GGC=61)
    };
    return (codon_idx >= 0 && codon_idx < 64) && X_CODON[codon_idx];
}

// ============================================================================
// Default AA bigram frequencies (bacterial proteome)
// ============================================================================

inline float default_aa_freq(int aa_idx) {
    static constexpr float AA_FREQ[20] = {
        0.082f, 0.013f, 0.054f, 0.067f, 0.039f,  // A C D E F
        0.071f, 0.023f, 0.057f, 0.058f, 0.096f,  // G H I K L
        0.024f, 0.045f, 0.047f, 0.039f, 0.055f,  // M N P Q R
        0.066f, 0.054f, 0.068f, 0.011f, 0.029f,  // S T V W Y
    };
    return (aa_idx >= 0 && aa_idx < 20) ? AA_FREQ[aa_idx] : 0.05f;
}

// ============================================================================
// Codon prior with stop penalty
// ============================================================================

inline float codon_prior(int codon_idx) {
    if (codon_idx == 10 || codon_idx == 11 || codon_idx == 14) {
        return 1e-6f;  // Stop codons (TAA, TAG, TGA)
    }
    return 1.0f / 61.0f;
}

/**
 * Joint codon prior using periodic Markov model.
 *
 * For a codon c = (c0, c1, c2) at positions (j, j+1, j+2) in frame f:
 *   P(c | context, frame) = P(c0 | ctx, phase0) × P(c1 | ctx', phase1) × P(c2 | ctx'', phase2)
 *
 * Where:
 *   - phase_k = (j+k - frame + 300) mod 3
 *   - ctx is the dinucleotide context before the codon
 *   - ctx' = (ctx[1], c0) for second base
 *   - ctx'' = (c0, c1) for third base
 *
 * @param codon_idx Index of true codon (0-63)
 * @param prev_ctx Dinucleotide context before codon (chars)
 * @param codon_start Position j of first codon base
 * @param frame Frame hypothesis (0, 1, 2)
 * @param model Self-trained periodic Markov model
 */
inline float codon_prior_joint(
    int codon_idx,
    char prev0, char prev1,
    int codon_start,
    int frame,
    const PeriodicMarkovModel& model) {

    // If model not trained, fall back to uniform prior
    if (!model.trained) {
        return codon_prior(codon_idx);
    }

    // Stop codon penalty
    if (codon_idx == 10 || codon_idx == 11 || codon_idx == 14) {
        return 1e-6f;
    }

    // Get codon bases
    char c0, c1, c2;
    idx_to_codon(codon_idx, c0, c1, c2);

    // Phases for each position
    int phase0 = (codon_start - frame + 300) % 3;
    int phase1 = (codon_start + 1 - frame + 300) % 3;
    int phase2 = (codon_start + 2 - frame + 300) % 3;

    // Context indices
    int ctx0 = PeriodicMarkovModel::ctx_index(prev0, prev1);
    int ctx1 = PeriodicMarkovModel::ctx_index(prev1, c0);
    int ctx2 = PeriodicMarkovModel::ctx_index(c0, c1);

    // Base indices
    int b0 = BasePriors::base_idx(c0);
    int b1 = BasePriors::base_idx(c1);
    int b2 = BasePriors::base_idx(c2);

    // Skip if any index is invalid
    if (ctx0 < 0 || ctx1 < 0 || ctx2 < 0 || b0 < 0 || b1 < 0 || b2 < 0) {
        return codon_prior(codon_idx);
    }

    // Joint probability
    float p0 = model.theta[phase0][ctx0][b0];
    float p1 = model.theta[phase1][ctx1][b1];
    float p2 = model.theta[phase2][ctx2][b2];

    return p0 * p1 * p2;
}

// ============================================================================
// Damage emission model
// ============================================================================

// Sparse candidate enumeration for observed base
// Returns viable true bases and their emission probabilities
struct ViableBases {
    char bases[2];       // Up to 2 viable true bases
    float emissions[2];  // Their emission probabilities
    int count;           // 1 or 2
};

inline ViableBases get_viable_bases_physical(char obs_phys, float delta_CT, float delta_GA) {
    ViableBases v;
    delta_CT = std::clamp(delta_CT, 0.0f, 1.0f);
    delta_GA = std::clamp(delta_GA, 0.0f, 1.0f);

    switch (obs_phys) {
        case 'A':
            // Could be true A (1.0) or true G that became A (delta_GA)
            v.bases[0] = 'A'; v.emissions[0] = 1.0f;
            if (delta_GA > 1e-6f) {
                v.bases[1] = 'G'; v.emissions[1] = delta_GA;
                v.count = 2;
            } else {
                v.count = 1;
            }
            break;
        case 'T':
            // Could be true T (1.0) or true C that became T (delta_CT)
            v.bases[0] = 'T'; v.emissions[0] = 1.0f;
            if (delta_CT > 1e-6f) {
                v.bases[1] = 'C'; v.emissions[1] = delta_CT;
                v.count = 2;
            } else {
                v.count = 1;
            }
            break;
        case 'C':
            // Only true C can produce observed C (with prob 1-delta_CT)
            v.bases[0] = 'C'; v.emissions[0] = 1.0f - delta_CT;
            v.count = 1;
            break;
        case 'G':
            // Only true G can produce observed G (with prob 1-delta_GA)
            v.bases[0] = 'G'; v.emissions[0] = 1.0f - delta_GA;
            v.count = 1;
            break;
        default:
            // Unknown base - uniform
            v.bases[0] = 'N'; v.emissions[0] = 0.25f;
            v.count = 1;
            break;
    }
    return v;
}

inline float damage_emission_physical(char true_phys, char obs_phys, float delta_CT, float delta_GA) {
    delta_CT = std::clamp(delta_CT, 0.0f, 1.0f);
    delta_GA = std::clamp(delta_GA, 0.0f, 1.0f);

    switch (true_phys) {
        case 'A': return (obs_phys == 'A') ? 1.0f : 0.0f;
        case 'T': return (obs_phys == 'T') ? 1.0f : 0.0f;
        case 'C':
            if (obs_phys == 'C') return 1.0f - delta_CT;
            if (obs_phys == 'T') return delta_CT;
            return 0.0f;
        case 'G':
            if (obs_phys == 'G') return 1.0f - delta_GA;
            if (obs_phys == 'A') return delta_GA;
            return 0.0f;
        default: return 0.25f;
    }
}

inline float damage_emission_oriented(
    char true_oriented, char obs_oriented,
    int oriented_pos, bool forward, int L,
    const DamageParams& dmg) {

    int phys_pos = forward ? oriented_pos : (L - 1 - oriented_pos);
    float delta_CT = dmg.delta_CT(phys_pos);
    float delta_GA = dmg.delta_GA(phys_pos, L);

    char true_phys = forward ? true_oriented : complement(true_oriented);
    char obs_phys = forward ? obs_oriented : complement(obs_oriented);

    return damage_emission_physical(true_phys, obs_phys, delta_CT, delta_GA);
}

/**
 * CpG-aware damage emission with context.
 *
 * For C→T damage: boosted when next TRUE base is G (CpG dinucleotide)
 * For G→A damage: boosted when prev TRUE base is C (CpG on opposite strand)
 *
 * @param true_oriented True base in oriented coordinates
 * @param obs_oriented Observed base in oriented coordinates
 * @param oriented_pos Position in oriented sequence
 * @param forward True if forward strand hypothesis
 * @param L Sequence length
 * @param dmg Damage parameters with CpG boost
 * @param next_true_oriented Next TRUE base (for CpG context), or 'N' if unknown
 * @param prev_true_oriented Previous TRUE base (for CpG context), or 'N' if unknown
 */
inline float damage_emission_oriented_cpg(
    char true_oriented, char obs_oriented,
    int oriented_pos, bool forward, int L,
    const DamageParams& dmg,
    char next_true_oriented = 'N',
    char prev_true_oriented = 'N') {

    int phys_pos = forward ? oriented_pos : (L - 1 - oriented_pos);

    char true_phys = forward ? true_oriented : complement(true_oriented);
    char obs_phys = forward ? obs_oriented : complement(obs_oriented);

    // Check CpG context in physical coordinates
    char next_true_phys = forward ? next_true_oriented : complement(next_true_oriented);
    char prev_true_phys = forward ? prev_true_oriented : complement(prev_true_oriented);

    // For forward strand: CpG = true_phys='C' and next_true_phys='G'
    // For reverse strand: coordinates are flipped, so check accordingly
    bool next_is_G = (next_true_phys == 'G');
    bool prev_is_C = (prev_true_phys == 'C');

    float delta_CT = dmg.delta_CT_context(phys_pos, next_is_G);
    float delta_GA = dmg.delta_GA_context(phys_pos, L, prev_is_C);

    return damage_emission_physical(true_phys, obs_phys, delta_CT, delta_GA);
}

// ============================================================================
// STRAND-ONLY LLR: Uses ALL positions (not just codon-aligned)
// ============================================================================

// Per-position emission probability marginalizing over true base
inline float position_marginal_ll(
    char obs, int phys_pos, int L,
    const DamageParams& dmg,
    const BasePriors& pi) {

    float delta_CT = dmg.delta_CT(phys_pos);
    float delta_GA = dmg.delta_GA(phys_pos, L);

    // Sum over all possible true physical bases
    float prob = 0.0f;

    // True = A: P(obs|A) = 1 if obs=A, else 0
    if (obs == 'A') prob += pi['A'];

    // True = T: P(obs|T) = 1 if obs=T, else 0
    if (obs == 'T') prob += pi['T'];

    // True = C: P(obs=C|C) = 1-delta_CT, P(obs=T|C) = delta_CT
    if (obs == 'C') prob += pi['C'] * (1.0f - delta_CT);
    if (obs == 'T') prob += pi['C'] * delta_CT;

    // True = G: P(obs=G|G) = 1-delta_GA, P(obs=A|G) = delta_GA
    if (obs == 'G') prob += pi['G'] * (1.0f - delta_GA);
    if (obs == 'A') prob += pi['G'] * delta_GA;

    return std::log(prob + 1e-30f);
}

/**
 * Strand-only log-likelihood ratio using ALL positions.
 *
 * DAMAGE-ONLY formulation (difference-of-differences):
 *   LLR_damage_only = LLR_strand(d) - LLR_strand(0)
 *
 * This subtracts the composition baseline, isolating damage signal.
 * At d=0, returns exactly 0 (no composition leakage).
 *
 * Forward hypothesis: observed seq is in physical orientation
 *   - Position i maps to physical position i
 *   - Observed base at i is the physical base
 *
 * Reverse hypothesis: observed seq is reverse complement of physical
 *   - Position i maps to physical position L-1-i
 *   - Observed base at i is the complement of physical base
 */
inline float strand_only_llr(
    const std::string& seq,
    const DamageParams& dmg,
    const BasePriors& pi = BasePriors()) {

    int L = static_cast<int>(seq.size());
    if (L < 10 || !dmg.has_damage()) return 0.0f;

    // Compute with damage
    float ll_F_d = 0.0f, ll_R_d = 0.0f;
    // Compute without damage (baseline)
    float ll_F_0 = 0.0f, ll_R_0 = 0.0f;

    DamageParams no_dmg;  // d_max = 0

    for (int i = 0; i < L; ++i) {
        char obs = seq[i];
        int idx = BasePriors::base_idx(obs);
        if (idx < 0) continue;  // Skip non-ACGT

        // Forward: physical position = i, observed = obs
        // For forward strand, we interpret the observed base directly
        ll_F_d += position_marginal_ll(obs, i, L, dmg, pi);
        ll_F_0 += position_marginal_ll(obs, i, L, no_dmg, pi);

        // Reverse: SAME physical position i (damage is physical, not biological)
        // But the TRUE base prior should be for the complement
        // This captures the asymmetry: if genome has GC skew, strand matters
        char obs_for_rev = complement(obs);
        ll_R_d += position_marginal_ll(obs_for_rev, i, L, dmg, pi);
        ll_R_0 += position_marginal_ll(obs_for_rev, i, L, no_dmg, pi);
    }

    // Difference-of-differences: isolate damage signal
    float llr_d = ll_F_d - ll_R_d;
    float llr_0 = ll_F_0 - ll_R_0;

    return llr_d - llr_0;
}

/**
 * Strand-only LLR focusing on terminal positions where damage is strongest.
 * Includes CpG context boost for stronger strand signal.
 */
inline float strand_only_llr_terminal(
    const std::string& seq,
    const DamageParams& dmg,
    const BasePriors& pi,
    int terminal_window = 15) {

    int L = static_cast<int>(seq.size());
    if (L < 2 * terminal_window || !dmg.has_damage()) return 0.0f;

    DamageParams no_dmg;
    float score = 0.0f;

    // Only process terminal positions
    for (int i = 0; i < L; ++i) {
        char obs = seq[i];
        int idx = BasePriors::base_idx(obs);
        if (idx < 0) continue;

        // Check CpG context
        bool is_cpg_C = (obs == 'C' || obs == 'c') && (i + 1 < L) &&
                        (seq[i+1] == 'G' || seq[i+1] == 'g');
        bool is_cpg_G = (obs == 'G' || obs == 'g') && (i > 0) &&
                        (seq[i-1] == 'C' || seq[i-1] == 'c');

        // Weight by damage probability (higher at terminals)
        // Boost weight for CpG sites
        float weight_5p = is_cpg_C ? dmg.delta_CT_context(i, true) : dmg.delta_CT(i);
        float weight_3p = is_cpg_G ? dmg.delta_GA_context(i, L, true) : dmg.delta_GA(i, L);
        float weight = std::max(weight_5p, weight_3p);
        if (weight < 0.01f) continue;  // Skip interior

        // FIXED: Damage position is always i (physical), not L-1-i
        char obs_for_rev = complement(obs);

        float ll_F_d = position_marginal_ll(obs, i, L, dmg, pi);
        float ll_F_0 = position_marginal_ll(obs, i, L, no_dmg, pi);
        float ll_R_d = position_marginal_ll(obs_for_rev, i, L, dmg, pi);
        float ll_R_0 = position_marginal_ll(obs_for_rev, i, L, no_dmg, pi);

        // Weighted difference-of-differences
        // Extra weight for CpG positions (they're more informative)
        float cpg_multiplier = (is_cpg_C || is_cpg_G) ? dmg.cpg_boost : 1.0f;
        float delta_F = ll_F_d - ll_F_0;
        float delta_R = ll_R_d - ll_R_0;
        score += weight * cpg_multiplier * (delta_F - delta_R);
    }

    return score;
}

/**
 * Estimate base priors from interior of reads (avoid damage bias).
 * Skip first/last skip_ends bases.
 */
inline BasePriors estimate_base_priors(
    const std::string& seq,
    int skip_ends = 10) {

    int L = static_cast<int>(seq.size());
    if (L <= 2 * skip_ends) {
        return BasePriors();  // Default uniform
    }

    std::array<int, 4> counts = {0, 0, 0, 0};
    for (int i = skip_ends; i < L - skip_ends; ++i) {
        int idx = BasePriors::base_idx(seq[i]);
        if (idx >= 0) counts[idx]++;
    }

    int total = counts[0] + counts[1] + counts[2] + counts[3];
    if (total < 10) return BasePriors();

    return BasePriors(
        static_cast<float>(counts[0]) / total,
        static_cast<float>(counts[1]) / total,
        static_cast<float>(counts[2]) / total,
        static_cast<float>(counts[3]) / total
    );
}

// ============================================================================
// FOURIER PERIODOGRAM: Detect 3-nt periodicity using spectral analysis
// ============================================================================

/**
 * Compute RNY bias score for a specific frame hypothesis.
 *
 * In coding sequences, there's a well-known RNY bias:
 * - Position 0 (codon pos 1): enriched for purines (R = A or G)
 * - Position 1 (codon pos 2): most variable (N)
 * - Position 2 (codon pos 3/wobble): varies by GC content
 *
 * We measure how well the base composition at each codon position matches
 * the expected RNY pattern for this frame hypothesis.
 *
 * Score = R_freq at pos0 + Y_freq at pos2 (both should be > 0.5 for correct frame)
 *
 * @param seq DNA sequence
 * @param frame Frame to test (0, 1, or 2)
 * @param skip_ends Skip terminal bases to avoid damage artifacts
 * @return RNY bias score (higher = better match to coding pattern)
 */
inline float compute_rny_bias_score(
    const std::string& seq,
    int frame,
    int skip_ends = 5) {

    int L = static_cast<int>(seq.size());
    if (L < 2 * skip_ends + 9) return 0.0f;

    int start = skip_ends;
    int end = L - skip_ends;

    // Count R (A,G) vs Y (C,T) at each codon position
    std::array<int, 3> R_count = {0, 0, 0};  // Purine count at each codon pos
    std::array<int, 3> Y_count = {0, 0, 0};  // Pyrimidine count at each codon pos

    for (int j = start; j < end; ++j) {
        char base = seq[j];
        int codon_pos = (j - frame + 300) % 3;

        if (base == 'A' || base == 'G' || base == 'a' || base == 'g') {
            R_count[codon_pos]++;
        } else if (base == 'C' || base == 'T' || base == 'c' || base == 't') {
            Y_count[codon_pos]++;
        }
    }

    // Compute RNY score
    // Position 0: expect R enrichment
    int total0 = R_count[0] + Y_count[0];
    float R_freq0 = total0 > 0 ? R_count[0] / static_cast<float>(total0) : 0.5f;

    // Position 2 dependence on GC is intentionally omitted in this score.

    // Score is just R enrichment at position 0 vs baseline 0.5
    // LLR form: log(R_freq0 / 0.5) - log((1-R_freq0) / 0.5)
    // = log(R_freq0) - log(1-R_freq0)

    float score = std::log(R_freq0 + 0.01f) - std::log(1.0f - R_freq0 + 0.01f);

    return score;
}

// Alias for backward compatibility
inline float fourier_periodicity_power(
    const std::string& seq,
    int frame,
    int skip_ends = 5) {
    return compute_rny_bias_score(seq, frame, skip_ends);
}

/**
 * RNY-based frame score.
 *
 * For each frame, compute the R enrichment at codon position 0.
 * Higher score = better match to coding RNY pattern.
 */
inline float score_fourier_periodicity(
    const std::string& seq,
    int frame) {

    // The RNY bias score is already in log-odds form
    return compute_rny_bias_score(seq, frame, 5);
}

// OPTIMIZED: Fast version with raw pointer
inline float score_fourier_periodicity_fast(
    const char* seq,
    int L,
    int frame,
    int skip_ends = 5) {

    if (L < 2 * skip_ends + 9) return 0.0f;

    int start = skip_ends;
    int end = L - skip_ends;

    // Use lookup table for R/Y classification
    alignas(64) static constexpr int8_t IS_PURINE[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // A=65, G=71
        0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // a=97, g=103
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    };
    alignas(64) static constexpr int8_t IS_PYRIMIDINE[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0, // C=67, T=84
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0, // c=99, t=116
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    };

    int R_count_0 = 0, Y_count_0 = 0;

    // Only count position 0 (the most discriminative)
    for (int j = start; j < end; ++j) {
        int codon_pos = (j - frame + 300) % 3;
        if (codon_pos == 0) {
            unsigned char c = static_cast<unsigned char>(seq[j]);
            R_count_0 += IS_PURINE[c];
            Y_count_0 += IS_PYRIMIDINE[c];
        }
    }

    int total0 = R_count_0 + Y_count_0;
    if (total0 == 0) return 0.0f;

    float R_freq0 = R_count_0 / static_cast<float>(total0);
    return std::log(R_freq0 + 0.01f) - std::log(1.0f - R_freq0 + 0.01f);
}

// ============================================================================
// Utility functions
// ============================================================================

inline float logsumexp(const float* logv, int n) {
    float maxv = logv[0];
    for (int i = 1; i < n; ++i) {
        if (logv[i] > maxv) maxv = logv[i];
    }
    if (!std::isfinite(maxv)) return -1e30f;

    float sum = 0.0f;
    for (int i = 0; i < n; ++i) {
        if (std::isfinite(logv[i])) {
            sum += std::exp(logv[i] - maxv);
        }
    }
    return maxv + std::log(sum);
}

// ============================================================================
// Codon posterior computation (needed for all scoring components)
// ============================================================================

struct CodonPosteriors {
    std::array<float, 64> post;  // P(true_codon = c | observed, damage)
    float x_membership;          // P(codon in X | observed, damage)
    std::array<float, 20> aa_post;  // P(AA = a | observed, damage)
    float stop_prob;             // P(stop | observed, damage)
    float log_marginal;          // log P(observed | frame, damage) - the marginal likelihood
    int best_codon;              // MAP codon (avoids 64-iteration search)
    float best_post;             // Posterior of best codon
};

/**
 * Compute codon posteriors with optional joint periodic Markov prior.
 *
 * If the periodic model is trained, uses:
 *   P(true_codon | observed, damage, context, frame) ∝
 *     P(true_codon | context, frame) × P(observed | true_codon, damage)
 *
 * Otherwise falls back to uniform prior.
 */
inline CodonPosteriors compute_codon_posteriors_joint(
    const std::string& oriented,
    int codon_start,
    bool forward,
    int L,
    int frame,
    const DamageParams& dmg,
    const PeriodicMarkovModel& model) {

    CodonPosteriors result;
    result.post.fill(0.0f);
    result.aa_post.fill(0.0f);
    result.x_membership = 0.0f;
    result.stop_prob = 0.0f;
    result.log_marginal = 0.0f;

    if (codon_start + 2 >= static_cast<int>(oriented.size())) {
        for (int c = 0; c < 64; ++c) result.post[c] = 1.0f / 64.0f;
        result.log_marginal = std::log(1.0f / 64.0f);  // Uniform marginal
        return result;
    }

    char o0 = oriented[codon_start];
    char o1 = oriented[codon_start + 1];
    char o2 = oriented[codon_start + 2];

    if (base_to_tcag(o0) < 0 || base_to_tcag(o1) < 0 || base_to_tcag(o2) < 0) {
        for (int c = 0; c < 64; ++c) result.post[c] = 1.0f / 64.0f;
        return result;
    }

    // Get context (previous bases for periodic model, next base for CpG)
    char prev0 = codon_start >= 2 ? oriented[codon_start - 2] : 'N';
    char prev1 = codon_start >= 1 ? oriented[codon_start - 1] : 'N';
    char next_base = codon_start + 3 < L ? oriented[codon_start + 3] : 'N';

    float log_terms[64];
    for (int c = 0; c < 64; ++c) {
        // Use joint prior if model is trained
        float prior = model.trained
            ? codon_prior_joint(c, prev0, prev1, codon_start, frame, model)
            : codon_prior(c);

        if (prior <= 0.0f) {
            log_terms[c] = -1e30f;
            continue;
        }

        char t0, t1, t2;
        idx_to_codon(c, t0, t1, t2);

        // CpG-aware emissions: use context from within the codon
        // Position 0: next true is t1, prev true is prev1 (from context)
        // Position 1: next true is t2, prev true is t0
        // Position 2: next true is next_base (outside codon), prev true is t1
        float e0 = damage_emission_oriented_cpg(t0, o0, codon_start, forward, L, dmg, t1, prev1);
        float e1 = damage_emission_oriented_cpg(t1, o1, codon_start + 1, forward, L, dmg, t2, t0);
        float e2 = damage_emission_oriented_cpg(t2, o2, codon_start + 2, forward, L, dmg, next_base, t1);

        float likelihood = e0 * e1 * e2;
        log_terms[c] = (likelihood > 0.0f) ? (std::log(prior) + std::log(likelihood)) : -1e30f;
    }

    float log_z = logsumexp(log_terms, 64);
    result.log_marginal = log_z;  // Store marginal likelihood for scoring

    for (int c = 0; c < 64; ++c) {
        result.post[c] = std::exp(log_terms[c] - log_z);

        if (is_x_code_codon(c)) {
            result.x_membership += result.post[c];
        }

        char aa = translate_codon_idx(c);
        if (aa == '*') {
            result.stop_prob += result.post[c];
        } else {
            int aa_idx = aa_to_idx(aa);
            if (aa_idx >= 0) {
                result.aa_post[aa_idx] += result.post[c];
            }
        }
    }

    return result;
}

// OPTIMIZED: Fast version taking raw pointer (no string overhead)
inline CodonPosteriors compute_codon_posteriors_joint_fast(
    const char* oriented,
    int codon_start,
    bool forward,
    int L,
    int frame,
    const DamageParams& dmg,
    const PeriodicMarkovModel& model) {

    CodonPosteriors result;
    result.x_membership = 0.0f;
    result.stop_prob = 0.0f;
    result.log_marginal = 0.0f;

    if (codon_start + 2 >= L) {
        result.post.fill(1.0f / 64.0f);
        result.aa_post.fill(0.0f);
        result.log_marginal = std::log(1.0f / 64.0f);
        return result;
    }

    const char o0 = oriented[codon_start];
    const char o1 = oriented[codon_start + 1];
    const char o2 = oriented[codon_start + 2];

    // Fast validation using lookup table
    if ((BASE_TO_TCAG[static_cast<unsigned char>(o0)] |
         BASE_TO_TCAG[static_cast<unsigned char>(o1)] |
         BASE_TO_TCAG[static_cast<unsigned char>(o2)]) < 0) {
        result.post.fill(1.0f / 64.0f);
        result.aa_post.fill(0.0f);
        return result;
    }

    // Context for periodic model
    const char prev0 = codon_start >= 2 ? oriented[codon_start - 2] : 'N';
    const char prev1 = codon_start >= 1 ? oriented[codon_start - 1] : 'N';
    const char next_base = codon_start + 3 < L ? oriented[codon_start + 3] : 'N';

    // Pre-compute damage rates at these positions (avoid recomputing per codon)
    const int phys0 = forward ? codon_start : (L - 1 - codon_start);
    const int phys1 = forward ? codon_start + 1 : (L - 2 - codon_start);
    const int phys2 = forward ? codon_start + 2 : (L - 3 - codon_start);

    const float delta_CT_0 = dmg.delta_CT(phys0);
    const float delta_CT_1 = dmg.delta_CT(phys1);
    const float delta_CT_2 = dmg.delta_CT(phys2);
    const float delta_GA_0 = dmg.delta_GA(phys0, L);
    const float delta_GA_1 = dmg.delta_GA(phys1, L);
    const float delta_GA_2 = dmg.delta_GA(phys2, L);

    alignas(64) float log_terms[64];
    float max_log = -1e30f;

    for (int c = 0; c < 64; ++c) {
        // Prior (uniform or from model)
        float prior = model.trained
            ? codon_prior_joint(c, prev0, prev1, codon_start, frame, model)
            : (1.0f / 64.0f);

        if (prior <= 0.0f) {
            log_terms[c] = -1e30f;
            continue;
        }

        // Decode codon
        char t0 = tcag_to_base(c >> 4);
        char t1 = tcag_to_base((c >> 2) & 3);
        char t2 = tcag_to_base(c & 3);

        // Fast damage emission (without CpG context)
        float e0 = damage_emission_physical(
            forward ? t0 : COMPLEMENT_TABLE[static_cast<unsigned char>(t0)],
            forward ? o0 : COMPLEMENT_TABLE[static_cast<unsigned char>(o0)],
            delta_CT_0, delta_GA_0);
        float e1 = damage_emission_physical(
            forward ? t1 : COMPLEMENT_TABLE[static_cast<unsigned char>(t1)],
            forward ? o1 : COMPLEMENT_TABLE[static_cast<unsigned char>(o1)],
            delta_CT_1, delta_GA_1);
        float e2 = damage_emission_physical(
            forward ? t2 : COMPLEMENT_TABLE[static_cast<unsigned char>(t2)],
            forward ? o2 : COMPLEMENT_TABLE[static_cast<unsigned char>(o2)],
            delta_CT_2, delta_GA_2);

        float likelihood = e0 * e1 * e2;
        log_terms[c] = (likelihood > 0.0f) ? (std::log(prior) + std::log(likelihood)) : -1e30f;
        if (log_terms[c] > max_log) max_log = log_terms[c];
    }

    // Log-sum-exp for normalization
    float sum = 0.0f;
    for (int c = 0; c < 64; ++c) {
        if (log_terms[c] > -1e20f) {
            sum += std::exp(log_terms[c] - max_log);
        }
    }
    float log_z = max_log + std::log(sum + 1e-30f);
    result.log_marginal = log_z;

    // Compute posteriors and derived quantities
    result.post.fill(0.0f);
    result.aa_post.fill(0.0f);

    for (int c = 0; c < 64; ++c) {
        float post_c = std::exp(log_terms[c] - log_z);
        result.post[c] = post_c;

        if (is_x_code_codon(c)) {
            result.x_membership += post_c;
        }

        char aa = translate_codon_idx(c);
        if (aa == '*') {
            result.stop_prob += post_c;
        } else {
            int aa_idx = aa_to_idx(aa);
            if (aa_idx >= 0) {
                result.aa_post[aa_idx] += post_c;
            }
        }
    }

    return result;
}

// SPARSE OPTIMIZATION: Only iterate over viable codons (1-8 instead of 64)
// This exploits the fact that damage emission returns 0 for most base combinations
inline CodonPosteriors compute_codon_posteriors_sparse(
    const char* oriented,
    int codon_start,
    bool forward,
    int L,
    int frame,
    const DamageParams& dmg,
    const PeriodicMarkovModel& model) {

    CodonPosteriors result;
    result.x_membership = 0.0f;
    result.stop_prob = 0.0f;
    result.log_marginal = 0.0f;
    result.post.fill(0.0f);
    result.aa_post.fill(0.0f);

    if (codon_start + 2 >= L) {
        result.post.fill(1.0f / 64.0f);
        result.log_marginal = std::log(1.0f / 64.0f);
        return result;
    }

    const char o0 = oriented[codon_start];
    const char o1 = oriented[codon_start + 1];
    const char o2 = oriented[codon_start + 2];

    // Fast validation
    if ((BASE_TO_TCAG[static_cast<unsigned char>(o0)] |
         BASE_TO_TCAG[static_cast<unsigned char>(o1)] |
         BASE_TO_TCAG[static_cast<unsigned char>(o2)]) < 0) {
        result.post.fill(1.0f / 64.0f);
        return result;
    }

    // Context for periodic model
    const char prev0 = codon_start >= 2 ? oriented[codon_start - 2] : 'N';
    const char prev1 = codon_start >= 1 ? oriented[codon_start - 1] : 'N';

    // Pre-compute damage rates
    const int phys0 = forward ? codon_start : (L - 1 - codon_start);
    const int phys1 = forward ? codon_start + 1 : (L - 2 - codon_start);
    const int phys2 = forward ? codon_start + 2 : (L - 3 - codon_start);

    const float delta_CT_0 = dmg.delta_CT(phys0);
    const float delta_CT_1 = dmg.delta_CT(phys1);
    const float delta_CT_2 = dmg.delta_CT(phys2);
    const float delta_GA_0 = dmg.delta_GA(phys0, L);
    const float delta_GA_1 = dmg.delta_GA(phys1, L);
    const float delta_GA_2 = dmg.delta_GA(phys2, L);

    // Get viable bases for each position (1-2 candidates each)
    // For forward strand, use observed bases directly
    // For reverse strand, complement the observed bases for physical interpretation
    char eff_o0 = forward ? o0 : COMPLEMENT_TABLE[static_cast<unsigned char>(o0)];
    char eff_o1 = forward ? o1 : COMPLEMENT_TABLE[static_cast<unsigned char>(o1)];
    char eff_o2 = forward ? o2 : COMPLEMENT_TABLE[static_cast<unsigned char>(o2)];

    ViableBases v0 = get_viable_bases_physical(eff_o0, delta_CT_0, delta_GA_0);
    ViableBases v1 = get_viable_bases_physical(eff_o1, delta_CT_1, delta_GA_1);
    ViableBases v2 = get_viable_bases_physical(eff_o2, delta_CT_2, delta_GA_2);

    // Maximum 2×2×2 = 8 viable codons
    float log_terms[8];
    int codon_indices[8];
    int n_viable = 0;
    float max_log = -1e30f;

    for (int i0 = 0; i0 < v0.count; ++i0) {
        char t0_phys = v0.bases[i0];
        float e0 = v0.emissions[i0];
        // Convert physical base back to oriented for codon index
        char t0 = forward ? t0_phys : COMPLEMENT_TABLE[static_cast<unsigned char>(t0_phys)];
        int b0 = base_to_tcag(t0);
        if (b0 < 0) continue;

        for (int i1 = 0; i1 < v1.count; ++i1) {
            char t1_phys = v1.bases[i1];
            float e1 = v1.emissions[i1];
            char t1 = forward ? t1_phys : COMPLEMENT_TABLE[static_cast<unsigned char>(t1_phys)];
            int b1 = base_to_tcag(t1);
            if (b1 < 0) continue;

            for (int i2 = 0; i2 < v2.count; ++i2) {
                char t2_phys = v2.bases[i2];
                float e2 = v2.emissions[i2];
                char t2 = forward ? t2_phys : COMPLEMENT_TABLE[static_cast<unsigned char>(t2_phys)];
                int b2 = base_to_tcag(t2);
                if (b2 < 0) continue;

                int c = (b0 << 4) | (b1 << 2) | b2;

                float prior = model.trained
                    ? codon_prior_joint(c, prev0, prev1, codon_start, frame, model)
                    : (1.0f / 64.0f);

                if (prior <= 0.0f) continue;

                float likelihood = e0 * e1 * e2;
                if (likelihood <= 0.0f) continue;

                // OPTIMIZATION: Single log instead of two
                float log_term = std::log(prior * likelihood);
                log_terms[n_viable] = log_term;
                codon_indices[n_viable] = c;
                if (log_term > max_log) max_log = log_term;
                ++n_viable;
            }
        }
    }

    if (n_viable == 0) {
        result.post.fill(1.0f / 64.0f);
        result.log_marginal = std::log(1.0f / 64.0f);
        result.best_codon = 0;
        result.best_post = 1.0f / 64.0f;
        return result;
    }

    // OPTIMIZATION: Fused exp loop - compute exp once, then normalize
    // This saves n_viable exp calls by computing unnormalized posteriors first
    float unnorm[8];
    float sum = 0.0f;
    for (int i = 0; i < n_viable; ++i) {
        float u = std::exp(log_terms[i] - max_log);
        unnorm[i] = u;
        sum += u;
    }

    float log_z = max_log + std::log(sum + 1e-30f);
    result.log_marginal = log_z;
    float inv_sum = 1.0f / (sum + 1e-30f);

    // Normalize and compute derived quantities in single pass
    result.best_codon = codon_indices[0];
    result.best_post = 0.0f;

    for (int i = 0; i < n_viable; ++i) {
        int c = codon_indices[i];
        float post_c = unnorm[i] * inv_sum;  // Division instead of exp
        result.post[c] = post_c;

        // Track MAP codon
        if (post_c > result.best_post) {
            result.best_post = post_c;
            result.best_codon = c;
        }

        if (is_x_code_codon(c)) {
            result.x_membership += post_c;
        }

        char aa = translate_codon_idx(c);
        if (aa == '*') {
            result.stop_prob += post_c;
        } else {
            int aa_idx = aa_to_idx(aa);
            if (aa_idx >= 0) {
                result.aa_post[aa_idx] += post_c;
            }
        }
    }

    return result;
}

// Legacy version without joint model
inline CodonPosteriors compute_codon_posteriors(
    const std::string& oriented,
    int codon_start,
    bool forward,
    int L,
    const DamageParams& dmg) {

    CodonPosteriors result;
    result.post.fill(0.0f);
    result.aa_post.fill(0.0f);
    result.x_membership = 0.0f;
    result.stop_prob = 0.0f;
    result.log_marginal = 0.0f;

    if (codon_start + 2 >= static_cast<int>(oriented.size())) {
        // Return uniform if out of bounds
        for (int c = 0; c < 64; ++c) result.post[c] = 1.0f / 64.0f;
        result.log_marginal = std::log(1.0f / 64.0f);
        return result;
    }

    char o0 = oriented[codon_start];
    char o1 = oriented[codon_start + 1];
    char o2 = oriented[codon_start + 2];

    if (base_to_tcag(o0) < 0 || base_to_tcag(o1) < 0 || base_to_tcag(o2) < 0) {
        for (int c = 0; c < 64; ++c) result.post[c] = 1.0f / 64.0f;
        result.log_marginal = std::log(1.0f / 64.0f);
        return result;
    }

    float log_terms[64];
    for (int c = 0; c < 64; ++c) {
        float prior = codon_prior(c);
        if (prior <= 0.0f) {
            log_terms[c] = -1e30f;
            continue;
        }

        char t0, t1, t2;
        idx_to_codon(c, t0, t1, t2);

        float e0 = damage_emission_oriented(t0, o0, codon_start, forward, L, dmg);
        float e1 = damage_emission_oriented(t1, o1, codon_start + 1, forward, L, dmg);
        float e2 = damage_emission_oriented(t2, o2, codon_start + 2, forward, L, dmg);

        float likelihood = e0 * e1 * e2;
        log_terms[c] = (likelihood > 0.0f) ? (std::log(prior) + std::log(likelihood)) : -1e30f;
    }

    float log_z = logsumexp(log_terms, 64);
    result.log_marginal = log_z;  // Store marginal likelihood for scoring

    for (int c = 0; c < 64; ++c) {
        result.post[c] = std::exp(log_terms[c] - log_z);

        if (is_x_code_codon(c)) {
            result.x_membership += result.post[c];
        }

        char aa = translate_codon_idx(c);
        if (aa == '*') {
            result.stop_prob += result.post[c];
        } else {
            int aa_idx = aa_to_idx(aa);
            if (aa_idx >= 0) {
                result.aa_post[aa_idx] += result.post[c];
            }
        }
    }

    return result;
}

// ============================================================================
// SCORING COMPONENT 1: Base codon marginal LLR
// ============================================================================

inline float score_codon_marginal(
    const std::vector<CodonPosteriors>& posteriors) {

    float ll = 0.0f;
    for (const auto& p : posteriors) {
        // This is already computed as part of posteriors, but we sum log-marginals
        float marginal = 0.0f;
        for (int c = 0; c < 64; ++c) {
            marginal += p.post[c];
        }
        if (marginal > 0.0f) {
            ll += std::log(marginal);
        }
    }
    return ll;
}

// ============================================================================
// SCORING COMPONENT 2: X Circular Code Enrichment
// Frame-specific: true frame enriched in X code codons
// ============================================================================

inline float score_x_code_llr(
    const std::vector<CodonPosteriors>& posteriors,
    float p1 = 0.35f,  // P(X | correct frame) - empirical ~35%
    float p0 = 0.20f)  // P(X | wrong frame) - empirical ~20%
{
    if (posteriors.empty()) return 0.0f;

    float log_ratio_in = std::log(p1 / p0);
    float log_ratio_out = std::log((1.0f - p1) / (1.0f - p0));

    float score = 0.0f;
    for (const auto& p : posteriors) {
        float wX = p.x_membership;
        score += wX * log_ratio_in + (1.0f - wX) * log_ratio_out;
    }

    return score;
}

// ============================================================================
// SCORING COMPONENT 3: AA Bigram Directional Bias
// Strand-specific: real coding has directional AA adjacency patterns
// ============================================================================

inline float score_aa_bigram_llr(
    const std::vector<CodonPosteriors>& posteriors) {

    if (posteriors.size() < 2) return 0.0f;

    // Simple model: prefer common AA pairs, penalize rare ones
    // This is directional because (A,B) != (B,A) in general

    // Log background frequency for denominator
    float score = 0.0f;

    for (size_t i = 0; i + 1 < posteriors.size(); ++i) {
        const auto& p1 = posteriors[i];
        const auto& p2 = posteriors[i + 1];

        // Skip if high stop probability
        if (p1.stop_prob > 0.5f || p2.stop_prob > 0.5f) continue;

        // Compute expected bigram log-likelihood
        for (int a = 0; a < 20; ++a) {
            float pa = p1.aa_post[a];
            if (pa < 1e-6f) continue;

            for (int b = 0; b < 20; ++b) {
                float pb = p2.aa_post[b];
                if (pb < 1e-6f) continue;

                // Simple directional model: some pairs are preferred
                // Use product of marginals as baseline
                float freq_a = default_aa_freq(a);
                float freq_b = default_aa_freq(b);

                // Directional bonus: certain pairs are more common in proteins
                float pair_bonus = 0.0f;

                // Hydrophobic pairs tend to cluster
                bool hydro_a = (a == 0 || a == 4 || a == 7 || a == 9 || a == 10 || a == 17 || a == 18);
                bool hydro_b = (b == 0 || b == 4 || b == 7 || b == 9 || b == 10 || b == 17 || b == 18);
                if (hydro_a && hydro_b) pair_bonus += 0.1f;

                // Charged pairs have specific patterns
                bool pos_a = (a == 8 || a == 14);  // K, R
                bool neg_b = (b == 2 || b == 3);   // D, E
                if (pos_a && neg_b) pair_bonus += 0.05f;

                float log_joint = std::log(freq_a * freq_b + 0.001f) + pair_bonus;
                float log_indep = std::log(freq_a) + std::log(freq_b);

                score += pa * pb * (log_joint - log_indep);
            }
        }
    }

    return score;
}

// ============================================================================
// SCORING COMPONENT 4: Stop Codon Positional Explainability
// Damage-aware: stops near physical ends are explainable by damage
// ============================================================================

inline float score_stop_position_llr(
    const std::string& oriented,
    int frame,
    bool forward,
    int L,
    const DamageParams& dmg,
    const std::vector<CodonPosteriors>& posteriors,
    float r0 = 0.02f)  // Background stop rate
{
    if (posteriors.empty()) return 0.0f;

    float score = 0.0f;

    for (size_t i = 0; i < posteriors.size(); ++i) {
        int codon_start = frame + static_cast<int>(i) * 3;
        float stop_mass = posteriors[i].stop_prob;

        // Compute explainability based on position
        // Stop is explainable if:
        // 1. Near physical 5' end AND could be CAA/CAG/CGA -> TAA/TAG/TGA
        // 2. Pattern matches C->T damage

        int phys_pos = forward ? codon_start : (L - 1 - codon_start);
        float delta = dmg.delta_CT(phys_pos);

        // Convertible precursors: CAA(26), CAG(27), CGA(30)
        // -> TAA(10), TAG(11), TGA(14)
        // Probability a sense codon could become stop via damage
        float p_convert = 0.0f;

        // CAA -> TAA: prior(CAA) * delta
        p_convert += codon_prior(26) * delta;  // CAA -> TAA
        // CAG -> TAG: prior(CAG) * delta
        p_convert += codon_prior(27) * delta;  // CAG -> TAG
        // CGA -> TGA: prior(CGA) * delta
        p_convert += codon_prior(30) * delta;  // CGA -> TGA

        // Total explainable stop probability
        float p_stop_explainable = r0 + (1.0f - r0) * p_convert;
        p_stop_explainable = std::clamp(p_stop_explainable, 1e-6f, 1.0f - 1e-6f);

        float log_p = std::log(p_stop_explainable);
        float log_1mp = std::log(1.0f - p_stop_explainable);
        float log_r0 = std::log(r0);
        float log_1mr0 = std::log(1.0f - r0);

        // Bernoulli LLR with soft stop indicator
        score += stop_mass * (log_p - log_r0);
        score += (1.0f - stop_mass) * (log_1mp - log_1mr0);
    }

    // Normalize by codon count to avoid length bias
    // (longer reads would otherwise accumulate higher scores)
    return score / static_cast<float>(posteriors.size());
}

// ============================================================================
// Main scoring function: combine all components
// OPTIMIZED: Zero allocations in hot path, compute revcomp once
// ============================================================================

inline std::array<FrameStrandResult, 6> score_all_hypotheses(
    const std::string& seq,
    const DamageParams& dmg,
    const ScoringWeights& weights = ScoringWeights(),
    bool populate_all_proteins = false) {

    std::array<FrameStrandResult, 6> results;
    float total_scores[6];

    const int L = static_cast<int>(seq.size());
    const char* fwd_seq = seq.c_str();

    // OPTIMIZATION: Compute reverse complement ONCE (not 3x in loop)
    std::string_view rc_view = revcomp_fast(fwd_seq, seq.size());
    const char* rev_seq = rc_view.data();

    // Compute strand-only LLR once (same for all frames of same strand)
    BasePriors pi = estimate_base_priors(seq, 10);
    float strand_llr = strand_only_llr_terminal(seq, dmg, pi, 15);

    // Compute strand hexamer LLR
    float strand_hex_fwd = 0.0f, strand_hex_rev = 0.0f;
    bool use_interior_only = (dmg.d_max_5p > 0.05f && seq.size() >= 45);
    strand::compute_both_llrs(fwd_seq, seq.size(), strand_hex_fwd, strand_hex_rev, use_interior_only);
    size_t effective_len = use_interior_only ? (seq.size() > 30 ? seq.size() - 30 : 1) : seq.size();
    float norm = static_cast<float>(effective_len > 6 ? effective_len - 5 : 1);
    strand_hex_fwd /= norm;
    strand_hex_rev /= norm;

    // Get learned model
    const auto& periodic_model = get_learned_periodic_model();

    // OPTIMIZATION: Pre-compute number of codons per frame
    const int n_codons[3] = {
        (L - 0) / 3,
        (L - 1) / 3,
        (L - 2) / 3
    };

    // OPTIMIZATION: Use thread-local buffers for proteins
    auto& buf = get_thread_buffers();

    // OPTIMIZATION: Stack-allocated codon posteriors (max ~680 codons for 2KB sequence)
    // Use fixed-size array instead of vector
    alignas(64) CodonPosteriors posteriors_buf[MAX_CODONS];

    for (int h = 0; h < 6; ++h) {
        const bool forward = (h < 3);
        const int frame = h % 3;
        const char* oriented = forward ? fwd_seq : rev_seq;
        // Clamp to MAX_CODONS to prevent buffer overflow for long sequences
        const int num_codons = std::min(n_codons[frame], static_cast<int>(MAX_CODONS));

        // OPTIMIZATION: Direct pointer arithmetic for codon iteration
        float s_codon = 0.0f;
        float s_xcode = 0.0f;
        float s_stoppos = 0.0f;
        int stop_count = 0;
        int interior_stops = 0;

        // Single pass: compute codon posteriors and all derived scores
        for (int c = 0; c < num_codons; ++c) {
            const int i = frame + c * 3;
            if (i + 2 >= L) break;

            // Compute codon posterior
            auto& cp = posteriors_buf[c];
            cp = compute_codon_posteriors_sparse(oriented, i, forward, L, frame, dmg, periodic_model);

            // Accumulate log marginal and store per-codon score for frameshift detection
            s_codon += cp.log_marginal;
            buf.codon_scores[h][c] = cp.log_marginal;

            // Observed codon -> protein
            const int obs_idx = codon_to_idx(oriented[i], oriented[i+1], oriented[i+2]);
            buf.protein_buffer[h][c] = translate_codon_idx(obs_idx);

            // MAP correction (but force 'X' for ambiguous/invalid codons)
            if (obs_idx < 0) {
                // Codon contains N or invalid base - cannot infer original
                buf.corrected_buffer[h][c] = 'X';
            } else {
                // Use pre-computed best codon from sparse computation (avoids 64-iteration search)
                buf.corrected_buffer[h][c] = translate_codon_idx(cp.best_codon);
            }

            // X-code membership (fused into loop)
            s_xcode += cp.x_membership;

            // Stop tracking for positional score
            if (cp.stop_prob > 0.5f) {
                stop_count++;
                if (c > 0 && c < num_codons - 1) interior_stops++;
            }
        }
        buf.protein_lens[h] = static_cast<size_t>(num_codons);

        // Normalize scores
        if (num_codons > 0) {
            s_codon /= static_cast<float>(num_codons);
            s_xcode = (s_xcode / static_cast<float>(num_codons)) - 0.31f; // X-code baseline
        }

        // Stop position score: penalize interior stops
        s_stoppos = -0.5f * static_cast<float>(interior_stops);

        // Strand scores
        float s_strand_only = forward ? strand_llr : -strand_llr;
        float s_strand_hex = forward ? strand_hex_fwd : strand_hex_rev;

        // Dicodon phase score - pass raw pointer for speed
        float s_periodic = dicodon::calculate_frame_score_fast(oriented, L, frame, Domain::GTDB);

        // Self-trained and GC-conditional (these use the original seq)
        float s_self_trained = periodic_model.trained ? periodic_model.score(seq, frame, forward) : 0.0f;
        float s_gc_conditional = get_gc_conditional_model().trained ?
            get_gc_conditional_model().score(seq, frame, forward) : 0.0f;

        // Fourier periodogram
        float s_fourier = score_fourier_periodicity_fast(oriented, L, frame);

        // AA bigram term (currently disabled)
        float s_bigram = 0.0f;

        // Combined score
        float total = weights.w_codon * s_codon
                    + weights.w_strand_only * s_strand_only
                    + weights.w_strand_hex * s_strand_hex
                    + weights.w_periodic * s_periodic
                    + weights.w_self_trained * s_self_trained
                    + weights.w_gc_conditional * s_gc_conditional
                    + weights.w_fourier * s_fourier
                    + weights.w_xcode * s_xcode
                    + weights.w_bigram * s_bigram
                    + weights.w_stoppos * s_stoppos;

        total_scores[h] = total;

        results[h].frame = frame;
        results[h].forward = forward;
        results[h].log_likelihood = total;
        results[h].score_codon = s_codon;
        results[h].score_strand_only = s_strand_only;
        results[h].score_strand_hex = s_strand_hex;
        results[h].score_periodic = s_periodic;
        results[h].score_self_trained = s_self_trained;
        results[h].score_fourier = s_fourier;
        results[h].score_xcode = s_xcode;
        results[h].score_bigram = s_bigram;
        results[h].score_stoppos = s_stoppos;
        // Defer string assignment until we know which hypothesis wins
        results[h].protein.clear();
        results[h].corrected_protein.clear();
    }

    // Compute posteriors via softmax
    float log_z = logsumexp(total_scores, 6);
    for (int h = 0; h < 6; ++h) {
        results[h].posterior = std::exp(total_scores[h] - log_z);
    }

    // Find best indices without full sort (faster for 6 elements)
    int best_idx = 0;
    float best_post = results[0].posterior;
    for (int h = 1; h < 6; ++h) {
        if (results[h].posterior > best_post) {
            best_post = results[h].posterior;
            best_idx = h;
        }
    }

    // Copy proteins - either just the winner or all 6 if requested
    if (populate_all_proteins) {
        // Copy all 6 proteins (needed for adaptive/all frames mode)
        for (int h = 0; h < 6; ++h) {
            results[h].protein.assign(buf.protein_buffer[h], buf.protein_lens[h]);
            results[h].corrected_protein.assign(buf.corrected_buffer[h], buf.protein_lens[h]);
        }
    } else {
        // Only copy the winning hypothesis (faster for normal mode)
        results[best_idx].protein.assign(buf.protein_buffer[best_idx],
                                          buf.protein_lens[best_idx]);
        results[best_idx].corrected_protein.assign(buf.corrected_buffer[best_idx],
                                                    buf.protein_lens[best_idx]);
    }

    // Sort by posterior (needed for API compatibility)
    std::sort(results.begin(), results.end(),
              [](const FrameStrandResult& a, const FrameStrandResult& b) {
                  return a.posterior > b.posterior;
              });

    return results;
}

// Legacy wrapper for compatibility
inline std::string revcomp(const std::string& seq) {
    auto view = revcomp_fast(seq.c_str(), seq.size());
    return std::string(view);
}

inline std::pair<FrameStrandResult, FrameStrandResult> select_best_per_strand(
    const std::string& seq,
    const DamageParams& dmg) {

    auto all = score_all_hypotheses(seq, dmg, ScoringWeights());

    FrameStrandResult best_fwd, best_rev;
    best_fwd.posterior = -1.0f;
    best_rev.posterior = -1.0f;

    for (const auto& r : all) {
        if (r.forward) {
            if (r.posterior > best_fwd.posterior) best_fwd = r;
        } else {
            if (r.posterior > best_rev.posterior) best_rev = r;
        }
    }

    return {best_fwd, best_rev};
}

inline float compute_read_damage_probability(
    const std::string& seq,
    const DamageParams& dmg,
    float prior_damaged = 0.5f) {

    if (!dmg.has_damage()) return 0.0f;

    auto results_dmg = score_all_hypotheses(seq, dmg, ScoringWeights());
    float ll_damaged = results_dmg[0].log_likelihood;

    DamageParams no_dmg;
    auto results_nodmg = score_all_hypotheses(seq, no_dmg, ScoringWeights());
    float ll_nodamaged = results_nodmg[0].log_likelihood;

    float log_bf = ll_damaged - ll_nodamaged;
    float log_prior_ratio = std::log(prior_damaged / (1.0f - prior_damaged));
    float log_posterior_odds = log_bf + log_prior_ratio;

    return 1.0f / (1.0f + std::exp(-std::clamp(log_posterior_odds, -20.0f, 20.0f)));
}

// ============================================================================
// EM training for periodic Markov model
// ============================================================================

/**
 * Run one EM iteration on a batch of sequences.
 *
 * E-step: Score each sequence with current model + other features to get posteriors
 * M-step: Accumulate weighted counts based on posteriors, then normalize
 *
 * @param sequences Vector of DNA sequences
 * @param dmg Damage parameters for scoring
 * @param model Model to update
 */
inline void em_iteration(
    const std::vector<std::string>& sequences,
    const DamageParams& dmg,
    PeriodicMarkovModel& model) {

    // Reset counts for M-step
    for (int p = 0; p < 3; ++p) {
        for (int c = 0; c < 16; ++c) {
            for (int b = 0; b < 4; ++b) {
                model.counts[p][c][b] = 0.0f;
            }
            model.totals[p][c] = 0.0f;
        }
    }

    // E-step: compute posteriors and accumulate
    ScoringWeights weights;
    for (const auto& seq : sequences) {
        if (seq.length() < 20) continue;

        auto results = score_all_hypotheses(seq, dmg, weights);

        // Accumulate weighted counts for each hypothesis
        for (const auto& r : results) {
            model.accumulate(seq, r.frame, r.forward, r.posterior);
        }
    }

    // M-step: normalize
    model.normalize();
}

/**
 * Train periodic Markov model via EM on a sample of sequences.
 *
 * @param sequences Vector of DNA sequences
 * @param dmg Damage parameters
 * @param n_iterations Number of EM iterations (typically 2-5)
 */
inline void train_periodic_model(
    const std::vector<std::string>& sequences,
    const DamageParams& dmg,
    int n_iterations = 3) {

    auto& model = get_learned_periodic_model();
    model.reset();

    for (int iter = 0; iter < n_iterations; ++iter) {
        em_iteration(sequences, dmg, model);
    }

    // Also train GC-conditional model
    auto& gc_model = get_gc_conditional_model();
    gc_model.reset();

    // Single pass to accumulate GC-binned counts
    ScoringWeights weights;
    for (const auto& seq : sequences) {
        if (seq.length() < 20) continue;

        auto results = score_all_hypotheses(seq, dmg, weights);

        for (const auto& r : results) {
            gc_model.accumulate(seq, r.frame, r.forward, r.posterior);
        }
    }
    gc_model.normalize();
}

}  // namespace codon
}  // namespace agp
