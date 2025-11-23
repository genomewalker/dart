/**
 * Frame selection for ancient DNA reads
 */

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include "agp/scoring.hpp"
#include "agp/simd_utils.hpp"
#include <cmath>
#include <algorithm>

namespace agp {

// Thread-local buffer for reverse complement to avoid allocations
static thread_local std::string tl_rc_buffer;

// Compute reverse complement into thread-local buffer using SIMD
static const std::string& reverse_complement_fast(const std::string& seq) {
    tl_rc_buffer.resize(seq.length());

#ifdef USE_AVX512
    simd::reverse_complement_avx512(seq.data(), tl_rc_buffer.data(), seq.length());
#else
    // Scalar fallback
    for (size_t i = 0; i < seq.length(); ++i) {
        tl_rc_buffer[i] = fast_complement(seq[seq.length() - 1 - i]);
    }
#endif
    return tl_rc_buffer;
}

std::string FrameSelector::reverse_complement(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.length());

    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        rc += fast_complement(*it);
    }
    return rc;
}

// Optimized frame scoring - no vector allocations
// Now includes damage-frame consistency for ancient DNA
static float score_frame_fast(
    const std::string& seq,
    int frame,
    std::string& protein_out) {

    if (seq.length() < static_cast<size_t>(3 + frame)) {
        return -1000.0f;
    }

    // Translate directly
    protein_out = translate_sequence(seq, frame);

    if (protein_out.empty()) {
        return -1000.0f;
    }

    // Calculate scores without allocations
    float codon = calculate_codon_usage_score_direct(seq, frame);
    float stop = FrameSelector::calculate_stop_penalty(protein_out);
    float aa_comp = FrameSelector::calculate_aa_composition_score(protein_out);
    float dicodon = calculate_dicodon_score(seq, frame);
    float dipeptide = calculate_dipeptide_score(protein_out);

    // GC3
    int gc3_count = 0, total3 = 0;
    for (size_t i = frame + 2; i < seq.length(); i += 3) {
        char c = fast_upper(seq[i]);
        if (c != 'N') {
            if (c == 'G' || c == 'C') gc3_count++;
            total3++;
        }
    }
    float gc3 = total3 > 0 ? static_cast<float>(gc3_count) / total3 : 0.5f;
    float gc3_bonus = (gc3 > 0.5f) ? 0.1f : 0.0f;

    float len_bonus = std::min(0.2f, static_cast<float>(protein_out.length()) / 100.0f);

    // Damage-frame consistency score
    // If sequence shows damage pattern, check if it's consistent with this frame
    float damage_consistency = FrameSelector::score_damage_frame_consistency(seq, frame);
    // Only apply if there's evidence of damage (T at 5' or A at 3')
    float damage_bonus = 0.0f;
    if (seq.length() > 0) {
        bool has_damage_signal = (fast_upper(seq[0]) == 'T') ||
                                  (fast_upper(seq[seq.length()-1]) == 'A');
        if (has_damage_signal) {
            damage_bonus = (damage_consistency - 0.5f) * 0.15f;  // +/- 0.075 max
        }
    }

    return 0.15f * codon + 0.28f * stop + 0.10f * aa_comp +
           0.13f * dicodon + 0.13f * dipeptide +
           0.05f * gc3_bonus + 0.05f * len_bonus +
           damage_bonus;  // Add damage consistency
}

FrameScore FrameSelector::score_frame(
    const std::string& seq,
    int frame,
    bool forward,
    const DamageProfile* damage) {

    FrameScore score;
    score.frame = frame;
    score.forward = forward;

    if (seq.length() < static_cast<size_t>(3 + frame)) {
        score.total_score = -1000.0f;
        return score;
    }

    // Get working sequence
    std::string working_seq;
    if (forward) {
        working_seq = seq;
    } else {
        working_seq = reverse_complement(seq);
    }

    // Use optimized scoring
    score.total_score = score_frame_fast(working_seq, frame, score.protein);

    if (score.total_score < -100.0f) {
        return score;
    }

    // Extract individual scores for reporting (optional)
    score.codon_score = calculate_codon_usage_score_direct(working_seq, frame);
    score.stop_codon_penalty = calculate_stop_penalty(score.protein);
    score.aa_composition_score = calculate_aa_composition_score(score.protein);

    if (damage) {
        score.damage_consistency = calculate_damage_consistency(
            seq, working_seq, frame, forward, *damage);
    } else {
        score.damage_consistency = 0.5f;
    }

    return score;
}

FrameScore FrameSelector::select_best_frame(
    const std::string& seq,
    const DamageProfile* damage) {

    if (seq.length() < 3) {
        FrameScore empty;
        empty.total_score = -1000.0f;
        return empty;
    }

    // Compute reverse complement ONCE using thread-local buffer (no allocation)
    const std::string& rc_seq = reverse_complement_fast(seq);

    FrameScore best;
    best.total_score = -1000.0f;

    std::string protein;

    // Score all 6 frames using optimized function
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        float total = score_frame_fast(seq, frame, protein);
        if (total > best.total_score) {
            best.total_score = total;
            best.frame = frame;
            best.forward = true;
            best.codon_score = calculate_codon_usage_score_direct(seq, frame);
            best.stop_codon_penalty = calculate_stop_penalty(protein);
            best.aa_composition_score = calculate_aa_composition_score(protein);
            best.damage_consistency = 0.5f;
            best.protein = std::move(protein);
            protein.clear();
        }

        // Reverse frame
        total = score_frame_fast(rc_seq, frame, protein);
        if (total > best.total_score) {
            best.total_score = total;
            best.frame = frame;
            best.forward = false;
            best.codon_score = calculate_codon_usage_score_direct(rc_seq, frame);
            best.stop_codon_penalty = calculate_stop_penalty(protein);
            best.aa_composition_score = calculate_aa_composition_score(protein);
            best.damage_consistency = 0.5f;
            best.protein = std::move(protein);
            protein.clear();
        }
    }

    return best;
}

std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand(
    const std::string& seq,
    const DamageProfile* damage) {

    if (seq.length() < 3) {
        FrameScore empty;
        empty.total_score = -1000.0f;
        return {empty, empty};
    }

    // Compute reverse complement ONCE using thread-local buffer (no allocation)
    const std::string& rc_seq = reverse_complement_fast(seq);

    FrameScore best_fwd, best_rev;
    best_fwd.total_score = -1000.0f;
    best_fwd.forward = true;
    best_rev.total_score = -1000.0f;
    best_rev.forward = false;

    std::string protein;

    // Score all 6 frames
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        float total = score_frame_fast(seq, frame, protein);
        if (total > best_fwd.total_score) {
            best_fwd.total_score = total;
            best_fwd.frame = frame;
            best_fwd.forward = true;
            best_fwd.codon_score = calculate_codon_usage_score_direct(seq, frame);
            best_fwd.stop_codon_penalty = calculate_stop_penalty(protein);
            best_fwd.aa_composition_score = calculate_aa_composition_score(protein);
            best_fwd.damage_consistency = 0.5f;
            best_fwd.protein = std::move(protein);
            protein.clear();
        }

        // Reverse frame
        total = score_frame_fast(rc_seq, frame, protein);
        if (total > best_rev.total_score) {
            best_rev.total_score = total;
            best_rev.frame = frame;
            best_rev.forward = false;
            best_rev.codon_score = calculate_codon_usage_score_direct(rc_seq, frame);
            best_rev.stop_codon_penalty = calculate_stop_penalty(protein);
            best_rev.aa_composition_score = calculate_aa_composition_score(protein);
            best_rev.damage_consistency = 0.5f;
            best_rev.protein = std::move(protein);
            protein.clear();
        }
    }

    return {best_fwd, best_rev};
}

std::vector<FrameScore> FrameSelector::score_all_frames(
    const std::string& seq,
    const DamageProfile* damage) {

    std::vector<FrameScore> scores;
    scores.reserve(6);

    for (int frame = 0; frame < 3; ++frame) {
        scores.push_back(score_frame(seq, frame, true, damage));
        scores.push_back(score_frame(seq, frame, false, damage));
    }

    std::sort(scores.begin(), scores.end(),
              [](const FrameScore& a, const FrameScore& b) {
                  return a.total_score > b.total_score;
              });

    return scores;
}

} // namespace agp
