/**
 * Frame selection for ancient DNA reads
 */

#include "agp/frame_selector.hpp"
#include "agp/codon_tables.hpp"
#include "agp/scoring.hpp"
#include "agp/simd_utils.hpp"
#include "agp/damage_model.hpp"
#include "agp/gtdb_hexamer_table.hpp"
#include "agp/hexamer_damage_lookup.hpp"
#include "agp/multi_domain_positional_hexamer.hpp"
#include <cmath>
#include <algorithm>
#include <climits>

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
// Now uses unified positional hexamer scoring with damage awareness
static float score_frame_fast(
    const std::string& seq,
    int frame,
    std::string& protein_out,
    float delta_5 = 0.0f,  // Optional damage parameters for damage-aware scoring
    float delta_3 = 0.0f) {

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
    float dipeptide = calculate_dipeptide_score(protein_out);

    // Unified positional hexamer score (replaces calculate_dicodon_score + positional_bonus)
    // Uses position-specific tables (START/INTERNAL/END) and damage-weighted scoring
    float hexamer = positional::calculate_hexamer_score_damage_aware(seq, frame, delta_5, delta_3);

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
    // Smooth sigmoid bonus: centered at 0.5, max bonus 0.1 at gc3=0.7+
    float gc3_bonus = 0.1f / (1.0f + std::exp(-10.0f * (gc3 - 0.5f)));

    float len_bonus = std::min(0.2f, static_cast<float>(protein_out.length()) / 100.0f);

    // Damage-frame consistency score (simple T at 5' / A at 3' check)
    float damage_bonus = 0.0f;
    if (seq.length() > 0 && (delta_5 > 0.01f || delta_3 > 0.01f)) {
        float damage_consistency = FrameSelector::score_damage_frame_consistency(seq, frame);
        bool has_damage_signal = (fast_upper(seq[0]) == 'T') ||
                                  (fast_upper(seq[seq.length()-1]) == 'A');
        if (has_damage_signal) {
            damage_bonus = (damage_consistency - 0.5f) * 0.15f;  // +/- 0.075 max
        }
    }

    // Dynamic weights: when no stop codons (stop=1.0), redistribute weight to hexamer
    // Base weights: codon=0.15, stop=0.28, aa=0.10, hexamer=0.18, dipeptide=0.13, gc3=0.05, len=0.05
    // Note: hexamer weight increased since it now includes positional information
    float stop_weight = 0.28f;
    float hexamer_weight = 0.18f;  // Increased from 0.13 to account for positional info
    float codon_weight = 0.15f;

    // If stop penalty doesn't discriminate (no internal stops), boost hexamer weight
    if (stop > 0.99f) {
        // No stop codons detected - rely more heavily on sequence patterns
        stop_weight = 0.0f;
        hexamer_weight = 0.32f;  // +0.14: hexamer is most discriminative for frame
        codon_weight = 0.29f;    // +0.14: codon usage also frame-specific
    }

    return codon_weight * codon + stop_weight * stop + 0.10f * aa_comp +
           hexamer_weight * hexamer + 0.13f * dipeptide +
           0.05f * gc3_bonus + 0.05f * len_bonus +
           damage_bonus;
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

// Helper: Calculate hexamer log-likelihood for strand comparison
// Returns sum of log(P_coding / P_uniform) for all hexamers in sequence
static float calculate_hexamer_strand_score(const std::string& seq, int frame) {
    if (seq.length() < static_cast<size_t>(frame + 6)) return 0.0f;

    float log_prob_sum = 0.0f;
    int hexamer_count = 0;

    // Background frequency (uniform = 1/4096)
    static constexpr float LOG_BACKGROUND = -8.317766f;  // log(1/4096)

    char hexbuf[7];
    hexbuf[6] = '\0';

    for (size_t i = frame; i + 5 < seq.length(); i += 3) {
        hexbuf[0] = fast_upper(seq[i]);
        hexbuf[1] = fast_upper(seq[i + 1]);
        hexbuf[2] = fast_upper(seq[i + 2]);
        hexbuf[3] = fast_upper(seq[i + 3]);
        hexbuf[4] = fast_upper(seq[i + 4]);
        hexbuf[5] = fast_upper(seq[i + 5]);

        // Skip if any N
        if (hexbuf[0] == 'N' || hexbuf[1] == 'N' || hexbuf[2] == 'N' ||
            hexbuf[3] == 'N' || hexbuf[4] == 'N' || hexbuf[5] == 'N') continue;

        uint32_t code = encode_hexamer(hexbuf);
        if (code == UINT32_MAX) continue;

        float freq = GTDB_HEXAMER_FREQ[code];
        if (freq < 1e-8f) freq = 1e-8f;

        log_prob_sum += std::log(freq) - LOG_BACKGROUND;
        hexamer_count++;
    }

    if (hexamer_count == 0) return 0.0f;
    return log_prob_sum;  // Return total, not average (to weight by length)
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

    // Pre-compute best hexamer scores for strand comparison
    // This captures the strand-specific codon usage signal
    float best_fwd_hexamer = -1000.0f;
    float best_rev_hexamer = -1000.0f;
    for (int frame = 0; frame < 3; ++frame) {
        float fwd_hex = calculate_hexamer_strand_score(seq, frame);
        float rev_hex = calculate_hexamer_strand_score(rc_seq, frame);
        if (fwd_hex > best_fwd_hexamer) best_fwd_hexamer = fwd_hex;
        if (rev_hex > best_rev_hexamer) best_rev_hexamer = rev_hex;
    }

    // Calculate strand comparison bonus
    // If one strand has much better hexamer evidence, boost its scores
    float hexamer_diff = best_fwd_hexamer - best_rev_hexamer;
    float strand_bonus_fwd = 0.0f;
    float strand_bonus_rev = 0.0f;

    // Apply sigmoid-like scaling: max bonus +/- 0.15
    // Threshold at diff > 2.0 (about 7x likelihood difference)
    if (hexamer_diff > 0.5f) {
        strand_bonus_fwd = std::min(0.15f, hexamer_diff * 0.05f);
    } else if (hexamer_diff < -0.5f) {
        strand_bonus_rev = std::min(0.15f, -hexamer_diff * 0.05f);
    }

    std::string protein;

    // Score all 6 frames with strand comparison bonus
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        float total = score_frame_fast(seq, frame, protein);
        total += strand_bonus_fwd;  // Add strand comparison bonus
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
        total += strand_bonus_rev;  // Add strand comparison bonus
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

float FrameSelector::calculate_wobble_enrichment(
    const std::string& seq,
    int frame,
    size_t terminal_region) {

    if (seq.length() < static_cast<size_t>(frame + 3)) {
        return 0.5f;  // Neutral if too short
    }

    // Limit analysis to terminal region
    size_t analyze_len = std::min(terminal_region, seq.length());

    // Count T's at each codon position within the terminal region
    int t_at_wobble = 0;    // T's at position 3 (wobble)
    int t_at_other = 0;     // T's at positions 1 and 2
    int total_pos3 = 0;     // Total bases at wobble position
    int total_other = 0;    // Total bases at other positions

    for (size_t i = frame; i < analyze_len; ++i) {
        char c = std::toupper(seq[i]);
        if (c == 'N') continue;

        int codon_pos = (i - frame) % 3;  // 0, 1, 2

        if (codon_pos == 2) {  // Wobble position (3rd in codon)
            total_pos3++;
            if (c == 'T') t_at_wobble++;
        } else {  // Position 1 or 2
            total_other++;
            if (c == 'T') t_at_other++;
        }
    }

    // Calculate enrichment: how much more T's at wobble vs other positions
    // Expected: if damage is evenly distributed, T freq should be similar
    // In correct frame: T's from C→T damage more likely at wobble (synonymous)
    if (total_pos3 == 0 || total_other == 0) {
        return 0.5f;
    }

    float t_freq_wobble = static_cast<float>(t_at_wobble) / total_pos3;
    float t_freq_other = static_cast<float>(t_at_other) / total_other;

    // Enrichment ratio normalized to 0-1
    // If t_freq_wobble > t_freq_other, this frame is more likely correct
    // Use logistic function to bound between 0 and 1
    float ratio = (t_freq_wobble + 0.01f) / (t_freq_other + 0.01f);
    float enrichment = 1.0f / (1.0f + std::exp(-(ratio - 1.0f) * 3.0f));

    return enrichment;
}

float FrameSelector::detect_damage_induced_stops(
    const std::string& seq,
    int frame) {

    if (seq.length() < static_cast<size_t>(frame + 3)) {
        return 1.0f;  // No stops = good
    }

    // Check first few codons for potential damage-induced stops
    // TAA could be from CAA (Gln), TAG from CAG (Gln), TGA from CGA (Arg)
    int suspicious_stops = 0;
    int codons_checked = 0;

    // Check first 5 codons (15 bp) at 5' end
    for (size_t i = frame; i + 2 < seq.length() && codons_checked < 5; i += 3) {
        char c1 = std::toupper(seq[i]);
        char c2 = std::toupper(seq[i + 1]);
        char c3 = std::toupper(seq[i + 2]);
        codons_checked++;

        // Check for damage-induced stop codons
        // TAA: could be from CAA (C→T damage at position 1)
        // TAG: could be from CAG (C→T damage at position 1)
        // TGA: could be from CGA (C→T damage at position 1)
        if (c1 == 'T') {
            if ((c2 == 'A' && c3 == 'A') ||  // TAA from CAA
                (c2 == 'A' && c3 == 'G') ||  // TAG from CAG
                (c2 == 'G' && c3 == 'A')) {  // TGA from CGA
                // This stop could be damage-induced
                // Weight more heavily if at position 0 (highest damage rate)
                float pos_weight = 1.0f - static_cast<float>(i - frame) / 15.0f;
                suspicious_stops += static_cast<int>(pos_weight * 2 + 0.5f);
            }
        }
    }

    // Return score: 1.0 = no suspicious stops, 0.0 = many suspicious stops
    // Each suspicious stop reduces confidence
    float score = 1.0f / (1.0f + suspicious_stops * 0.5f);
    return score;
}

float FrameSelector::calculate_strand_confidence(
    const std::string& seq,
    float fwd_score,
    float rev_score,
    int fwd_frame,
    int rev_frame,
    const SampleDamageProfile& sample_profile) {

    // Base confidence from score difference
    float score_diff = std::abs(fwd_score - rev_score);
    float base_conf = 0.5f + score_diff * 0.5f;  // 0.5 to 1.0

    // Check if sample has significant damage (only then can we use damage signals)
    bool has_damage = sample_profile.max_damage_5prime > 0.05f ||
                      sample_profile.max_damage_3prime > 0.05f;

    if (!has_damage) {
        // Without damage signal, strand prediction is based only on coding score
        return std::min(1.0f, base_conf);
    }

    // Compute reverse complement for reverse strand analysis
    std::string rc = reverse_complement(seq);

    // Calculate wobble enrichment for both strands
    float wobble_fwd = calculate_wobble_enrichment(seq, fwd_frame);
    float wobble_rev = calculate_wobble_enrichment(rc, rev_frame);

    // Detect damage-induced stops
    float stop_fwd = detect_damage_induced_stops(seq, fwd_frame);
    float stop_rev = detect_damage_induced_stops(rc, rev_frame);

    // Combine signals
    // If forward has better wobble enrichment AND fewer suspicious stops
    float fwd_signal = wobble_fwd * stop_fwd;
    float rev_signal = wobble_rev * stop_rev;

    // Calculate strand signal confidence
    float signal_ratio = (fwd_signal + 0.01f) / (rev_signal + 0.01f);
    float signal_conf;
    if (fwd_score > rev_score) {
        // We're predicting forward - check if signals agree
        signal_conf = (signal_ratio > 1.0f) ? 0.5f + 0.5f * (1.0f - 1.0f/signal_ratio) : 0.5f * signal_ratio;
    } else {
        // We're predicting reverse - check if signals agree
        signal_conf = (signal_ratio < 1.0f) ? 0.5f + 0.5f * (1.0f - signal_ratio) : 0.5f / signal_ratio;
    }

    // Weight sample damage level: stronger damage = more reliable signal
    float damage_weight = std::min(1.0f, (sample_profile.max_damage_5prime + sample_profile.max_damage_3prime) * 2.0f);

    // Combined confidence
    float combined = base_conf * (1.0f - damage_weight * 0.3f) + signal_conf * damage_weight * 0.3f;

    return std::max(0.5f, std::min(1.0f, combined));
}

// Damage-aware frame scoring with Bayesian correction
static float score_frame_damage_aware(
    const std::string& seq,
    int frame,
    const DamageModel& damage_model,
    std::string& protein_out,
    std::string& corrected_dna_out,
    std::string& corrected_protein_out,
    size_t& stops_corrected_out,
    size_t& aa_changes_out) {

    if (seq.length() < static_cast<size_t>(3 + frame)) {
        return -1000.0f;
    }

    // Step 1: Apply Bayesian stop codon correction to DNA
    auto [corrected_dna, stops_fixed] = damage_model.correct_stop_codons_bayesian(seq, frame);
    corrected_dna_out = corrected_dna;
    stops_corrected_out = stops_fixed;

    // Step 2: Translate BOTH original and corrected DNA
    protein_out = translate_sequence(seq, frame);

    // Step 3: Use Bayesian AA probability translation for corrected protein
    auto [bayesian_protein, confidences] = damage_model.translate_with_confidence(corrected_dna, frame);
    corrected_protein_out = bayesian_protein;

    // Count AA changes
    aa_changes_out = 0;
    for (size_t i = 0; i < protein_out.length() && i < corrected_protein_out.length(); ++i) {
        if (protein_out[i] != corrected_protein_out[i]) {
            aa_changes_out++;
        }
    }

    if (corrected_protein_out.empty()) {
        return -1000.0f;
    }

    // Step 4: Score the CORRECTED sequences (this is the key insight!)
    // The correct frame will score better after correction
    float codon = calculate_codon_usage_score_direct(corrected_dna, frame);
    float stop = FrameSelector::calculate_stop_penalty(corrected_protein_out);  // Use corrected!
    float aa_comp = FrameSelector::calculate_aa_composition_score(corrected_protein_out);
    float dicodon = calculate_dicodon_score(corrected_dna, frame);
    float dipeptide = calculate_dipeptide_score(corrected_protein_out);

    // GC3
    int gc3_count = 0, total3 = 0;
    for (size_t i = frame + 2; i < corrected_dna.length(); i += 3) {
        char c = fast_upper(corrected_dna[i]);
        if (c != 'N') {
            if (c == 'G' || c == 'C') gc3_count++;
            total3++;
        }
    }
    float gc3 = total3 > 0 ? static_cast<float>(gc3_count) / total3 : 0.5f;
    float gc3_bonus = (gc3 > 0.5f) ? 0.1f : 0.0f;

    float len_bonus = std::min(0.2f, static_cast<float>(corrected_protein_out.length()) / 100.0f);

    // Bonus for stop codon corrections (indicates this frame might be correct)
    // Damage-induced stops only get corrected if they look like artifacts
    float correction_bonus = 0.0f;
    if (stops_fixed > 0) {
        // Each stop corrected adds confidence this is the right frame
        // But cap it to avoid over-weighting
        correction_bonus = std::min(0.15f, stops_fixed * 0.05f);
    }

    return 0.15f * codon + 0.28f * stop + 0.10f * aa_comp +
           0.13f * dicodon + 0.13f * dipeptide +
           0.05f * gc3_bonus + 0.05f * len_bonus +
           correction_bonus;
}

/**
 * Calculate damage-adjusted stop codon penalty multiplier
 * Returns a multiplier (0.0-1.0) that's more lenient for damage-induced stops
 *
 * This replaces calculate_stop_penalty() for damage-aware frame selection.
 * Instead of treating all internal stops equally, it weights them by damage probability.
 */
static float calculate_damage_adjusted_stop_multiplier(
    const std::string& seq,
    const std::string& protein,
    int frame,
    const DamageModel& damage_model) {

    if (protein.empty()) return 1.0f;

    size_t seq_len = seq.length();

    // Helper to get hexamer and calculate log probability
    auto get_hexamer_log_prob = [](const std::string& s, size_t pos) -> float {
        if (pos + 6 > s.length()) return -10.0f;
        uint32_t code = agp::encode_hexamer(s.c_str() + pos);
        if (code == UINT32_MAX) return -10.0f;
        float freq = agp::GTDB_HEXAMER_FREQ[code];
        if (freq < 1e-8f) freq = 1e-8f;
        return std::log(freq);
    };

    // Count internal stops and their damage probabilities
    int total_internal_stops = 0;
    float damage_weighted_stops = 0.0f;

    for (size_t i = 0; i + 1 < protein.length(); ++i) {  // i + 1 to exclude terminal stop
        if (protein[i] == '*') {
            total_internal_stops++;

            // Get codon position in DNA
            size_t codon_start = frame + (i * 3);
            if (codon_start + 3 > seq_len) {
                damage_weighted_stops += 1.0f;
                continue;
            }

            // Get stop codon
            char c1 = seq[codon_start];
            char c2 = seq[codon_start + 1];
            char c3 = seq[codon_start + 2];

            // Determine if stop is near terminus (damage-prone region)
            // Be very conservative: only first/last 2 codons
            bool near_5prime = (codon_start < 6);  // First 2 codons only
            bool near_3prime = (codon_start + 3 > seq_len - 6);  // Last 2 codons only

            // Check if this is a potential C→T damage (TAA from CAA, TAG from CAG, TGA from CGA)
            bool potential_ct_damage = ((c1 == 'T' || c1 == 't') &&
                                       ((c2 == 'A' || c2 == 'a') || (c2 == 'G' || c2 == 'g')));

            // Check if this is a potential G→A damage (TAA from TAG, TGA from TGG on RC)
            bool potential_ga_damage = ((c2 == 'A' || c2 == 'a') && (c3 == 'A' || c3 == 'a'));

            float stop_weight = 1.0f;  // Default: full penalty

            // Only apply hexamer-based evaluation if stop is in damage-prone region
            if ((near_5prime && potential_ct_damage) || (near_3prime && potential_ga_damage)) {
                // Find hexamer context (align to codon boundary)
                size_t hexamer_start = (codon_start / 3) * 3;
                if (hexamer_start + 6 > seq_len && seq_len >= 6) {
                    hexamer_start = seq_len - 6;
                }

                if (hexamer_start + 6 <= seq_len) {
                    // Get observed (stop-containing) hexamer probability
                    float log_p_stop = get_hexamer_log_prob(seq, hexamer_start);

                    // Create corrected hexamer (stop → sense codon)
                    std::string corrected_hex = seq.substr(hexamer_start, 6);
                    size_t stop_pos_in_hex = codon_start - hexamer_start;

                    // Try most likely correction based on context
                    if (potential_ct_damage && stop_pos_in_hex < 6) {
                        corrected_hex[stop_pos_in_hex] = 'C';  // T→C
                    } else if (potential_ga_damage && stop_pos_in_hex + 2 < 6) {
                        corrected_hex[stop_pos_in_hex + 2] = 'G';  // A→G at position 3
                    }

                    float log_p_corrected = get_hexamer_log_prob(corrected_hex, 0);
                    float log_likelihood_ratio = log_p_corrected - log_p_stop;

                    // If corrected hexamer is much more likely, this stop is probably damaged
                    // LLR > 3.0 means corrected is 20x more likely
                    // Be ULTRA conservative: only discount with overwhelming evidence
                    if (log_likelihood_ratio > 3.0f) {
                        // Overwhelming evidence (20x more likely): reduce weight to 0.4
                        stop_weight = 0.4f;
                    } else if (log_likelihood_ratio > 2.3f) {
                        // Very strong evidence (10x more likely): reduce weight to 0.6
                        stop_weight = 0.6f;
                    }
                    // else: not strong enough evidence, keep full penalty (stop_weight = 1.0)
                }
            }

            damage_weighted_stops += stop_weight;
        }
    }

    // Return multiplier based on weighted stop count
    // Use standard penalty scale: each weighted stop has strong impact
    // 0 stops → 1.0, 1 full stop → 0.2, 2+ full stops → 0.05
    if (damage_weighted_stops < 0.01f) return 1.0f;
    if (damage_weighted_stops < 1.0f) return 1.0f - (0.8f * damage_weighted_stops);  // Linear decay
    return 0.2f / damage_weighted_stops;  // Strong penalty for multiple stops
}

std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand_damage_aware(
    const std::string& seq,
    const DamageModel& damage_model) {

    if (seq.length() < 3) {
        FrameScore empty;
        empty.total_score = -1000.0f;
        return {empty, empty};
    }

    // Get damage parameters from model for hexamer weighting
    float delta_5 = damage_model.get_delta_5();
    float delta_3 = damage_model.get_delta_3();

    // Compute reverse complement ONCE
    const std::string& rc_seq = reverse_complement_fast(seq);

    FrameScore best_fwd, best_rev;
    best_fwd.total_score = -1000.0f;
    best_fwd.forward = true;
    best_rev.total_score = -1000.0f;
    best_rev.forward = false;

    std::string protein;

    // Score all 6 frames with damage-adjusted scoring
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        protein = translate_sequence(seq, frame);
        if (!protein.empty()) {
            // Calculate all score components
            float codon = calculate_codon_usage_score_direct(seq, frame);
            float stop = calculate_damage_adjusted_stop_multiplier(seq, protein, frame, damage_model);
            float aa_comp = calculate_aa_composition_score(protein);
            float dipeptide = calculate_dipeptide_score(protein);

            // Unified positional hexamer score with damage weighting
            float hexamer = positional::calculate_hexamer_score_damage_aware(seq, frame, delta_5, delta_3);

            // GC3 bonus
            int gc3_count = 0, total3 = 0;
            for (size_t i = frame + 2; i < seq.length(); i += 3) {
                char c = fast_upper(seq[i]);
                if (c != 'N') {
                    if (c == 'G' || c == 'C') gc3_count++;
                    total3++;
                }
            }
            float gc3 = total3 > 0 ? static_cast<float>(gc3_count) / total3 : 0.5f;
            float gc3_bonus = 0.1f / (1.0f + std::exp(-10.0f * (gc3 - 0.5f)));

            float len_bonus = std::min(0.2f, static_cast<float>(protein.length()) / 100.0f);

            // Damage consistency
            float damage_consistency = score_damage_frame_consistency(seq, frame);
            float damage_bonus = 0.0f;
            if (seq.length() > 0) {
                bool has_damage_signal = (fast_upper(seq[0]) == 'T') ||
                                          (fast_upper(seq[seq.length()-1]) == 'A');
                if (has_damage_signal) {
                    damage_bonus = (damage_consistency - 0.5f) * 0.15f;
                }
            }

            // Wobble position T-enrichment (helps identify correct frame in damaged DNA)
            float wobble_enrich = calculate_wobble_enrichment(seq, frame, 15);
            float wobble_bonus = (wobble_enrich - 0.5f) * 0.10f;  // +/- 0.05 max

            // Dynamic weights
            float stop_weight = 0.28f;
            float hexamer_weight = 0.18f;
            float codon_weight = 0.15f;

            if (stop > 0.99f) {
                stop_weight = 0.0f;
                hexamer_weight = 0.32f;
                codon_weight = 0.29f;
            }

            float total = codon_weight * codon + stop_weight * stop + 0.10f * aa_comp +
                         hexamer_weight * hexamer + 0.13f * dipeptide +
                         0.05f * gc3_bonus + 0.05f * len_bonus +
                         damage_bonus + wobble_bonus;

            if (total > best_fwd.total_score) {
                best_fwd.total_score = total;
                best_fwd.frame = frame;
                best_fwd.forward = true;
                best_fwd.codon_score = codon;
                best_fwd.stop_codon_penalty = stop;
                best_fwd.aa_composition_score = aa_comp;
                best_fwd.damage_consistency = damage_consistency;
                best_fwd.protein = std::move(protein);
            }
        }

        // Reverse frame
        protein = translate_sequence(rc_seq, frame);
        if (!protein.empty()) {
            // Calculate all score components
            float codon = calculate_codon_usage_score_direct(rc_seq, frame);
            float stop = calculate_damage_adjusted_stop_multiplier(rc_seq, protein, frame, damage_model);
            float aa_comp = calculate_aa_composition_score(protein);
            float dipeptide = calculate_dipeptide_score(protein);

            // Unified positional hexamer score with damage weighting
            float hexamer = positional::calculate_hexamer_score_damage_aware(rc_seq, frame, delta_5, delta_3);

            // GC3 bonus
            int gc3_count = 0, total3 = 0;
            for (size_t i = frame + 2; i < rc_seq.length(); i += 3) {
                char c = fast_upper(rc_seq[i]);
                if (c != 'N') {
                    if (c == 'G' || c == 'C') gc3_count++;
                    total3++;
                }
            }
            float gc3 = total3 > 0 ? static_cast<float>(gc3_count) / total3 : 0.5f;
            float gc3_bonus = 0.1f / (1.0f + std::exp(-10.0f * (gc3 - 0.5f)));

            float len_bonus = std::min(0.2f, static_cast<float>(protein.length()) / 100.0f);

            // Damage consistency
            float damage_consistency = score_damage_frame_consistency(rc_seq, frame);
            float damage_bonus = 0.0f;
            if (rc_seq.length() > 0) {
                bool has_damage_signal = (fast_upper(rc_seq[0]) == 'T') ||
                                          (fast_upper(rc_seq[rc_seq.length()-1]) == 'A');
                if (has_damage_signal) {
                    damage_bonus = (damage_consistency - 0.5f) * 0.15f;
                }
            }

            // Wobble position T-enrichment (helps identify correct frame in damaged DNA)
            float wobble_enrich = calculate_wobble_enrichment(rc_seq, frame, 15);
            float wobble_bonus = (wobble_enrich - 0.5f) * 0.10f;  // +/- 0.05 max

            // Dynamic weights
            float stop_weight = 0.28f;
            float hexamer_weight = 0.18f;
            float codon_weight = 0.15f;

            if (stop > 0.99f) {
                stop_weight = 0.0f;
                hexamer_weight = 0.32f;
                codon_weight = 0.29f;
            }

            float total = codon_weight * codon + stop_weight * stop + 0.10f * aa_comp +
                         hexamer_weight * hexamer + 0.13f * dipeptide +
                         0.05f * gc3_bonus + 0.05f * len_bonus +
                         damage_bonus + wobble_bonus;

            if (total > best_rev.total_score) {
                best_rev.total_score = total;
                best_rev.frame = frame;
                best_rev.forward = false;
                best_rev.codon_score = codon;
                best_rev.stop_codon_penalty = stop;
                best_rev.aa_composition_score = aa_comp;
                best_rev.damage_consistency = damage_consistency;
                best_rev.protein = std::move(protein);
            }
        }
    }

    return {best_fwd, best_rev};
}

std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand_stop_only(
    const std::string& seq) {

    if (seq.length() < 3) {
        FrameScore empty;
        empty.total_score = -1000.0f;
        return {empty, empty};
    }

    // Compute reverse complement ONCE
    const std::string& rc_seq = reverse_complement_fast(seq);

    FrameScore best_fwd, best_rev;
    best_fwd.total_score = -1000.0f;
    best_fwd.forward = true;
    best_rev.total_score = -1000.0f;
    best_rev.forward = false;

    int best_fwd_stops = INT_MAX;
    int best_rev_stops = INT_MAX;

    std::string protein;

    // Score all 6 frames using ONLY stop codon count
    // Score is based on stop penalty multiplier (0-1 range) for compatibility
    // with min_coding_prob filter: 0 stops → 1.0, 1 stop → 0.2, 2+ stops → 0.05
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        protein = translate_sequence(seq, frame);
        if (!protein.empty()) {
            // Count internal stop codons (exclude terminal)
            int stops = 0;
            for (size_t i = 0; i + 1 < protein.length(); ++i) {
                if (protein[i] == '*') stops++;
            }

            // Use stop penalty as score (0-1 range)
            // Add small frame bonus for tie-breaking: prefer frame 0 > 1 > 2
            float stop_penalty = calculate_stop_penalty(protein);
            float score = stop_penalty - 0.0001f * frame;

            // Update best if fewer stops (or same stops but lower frame)
            if (stops < best_fwd_stops || (stops == best_fwd_stops && frame < best_fwd.frame)) {
                best_fwd_stops = stops;
                best_fwd.total_score = score;
                best_fwd.frame = frame;
                best_fwd.forward = true;
                best_fwd.codon_score = 0.0f;
                best_fwd.stop_codon_penalty = stop_penalty;
                best_fwd.aa_composition_score = 0.0f;
                best_fwd.damage_consistency = 0.5f;
                best_fwd.protein = std::move(protein);
            }
        }

        // Reverse frame
        protein = translate_sequence(rc_seq, frame);
        if (!protein.empty()) {
            // Count internal stop codons (exclude terminal)
            int stops = 0;
            for (size_t i = 0; i + 1 < protein.length(); ++i) {
                if (protein[i] == '*') stops++;
            }

            // Use stop penalty as score (0-1 range)
            float stop_penalty = calculate_stop_penalty(protein);
            float score = stop_penalty - 0.0001f * frame;

            // Update best if fewer stops (or same stops but lower frame)
            if (stops < best_rev_stops || (stops == best_rev_stops && frame < best_rev.frame)) {
                best_rev_stops = stops;
                best_rev.total_score = score;
                best_rev.frame = frame;
                best_rev.forward = false;
                best_rev.codon_score = 0.0f;
                best_rev.stop_codon_penalty = stop_penalty;
                best_rev.aa_composition_score = 0.0f;
                best_rev.damage_consistency = 0.5f;
                best_rev.protein = std::move(protein);
            }
        }
    }

    return {best_fwd, best_rev};
}

std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand_stop_priority(
    const std::string& seq) {

    if (seq.length() < 3) {
        FrameScore empty;
        empty.total_score = -1000.0f;
        return {empty, empty};
    }

    // Compute reverse complement ONCE
    const std::string& rc_seq = reverse_complement_fast(seq);

    // Store all 6 frame results with stop counts and hexamer scores
    struct FrameInfo {
        int frame;
        bool forward;
        int stop_count;
        float hexamer_llr;
        std::string protein;
        float stop_penalty;
    };

    std::vector<FrameInfo> fwd_frames(3);
    std::vector<FrameInfo> rev_frames(3);

    // Score all 6 frames
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        std::string protein = translate_sequence(seq, frame);
        int stops = 0;
        for (size_t i = 0; i + 1 < protein.length(); ++i) {
            if (protein[i] == '*') stops++;
        }
        float hexamer_llr = calculate_hexamer_llr_score(seq, frame);
        float stop_pen = calculate_stop_penalty(protein);
        fwd_frames[frame] = {frame, true, stops, hexamer_llr, std::move(protein), stop_pen};

        // Reverse frame
        protein = translate_sequence(rc_seq, frame);
        stops = 0;
        for (size_t i = 0; i + 1 < protein.length(); ++i) {
            if (protein[i] == '*') stops++;
        }
        hexamer_llr = calculate_hexamer_llr_score(rc_seq, frame);
        stop_pen = calculate_stop_penalty(protein);
        rev_frames[frame] = {frame, false, stops, hexamer_llr, std::move(protein), stop_pen};
    }

    // Sort by: (1) stop_count ascending, (2) hexamer_llr descending
    auto compare = [](const FrameInfo& a, const FrameInfo& b) {
        if (a.stop_count != b.stop_count) return a.stop_count < b.stop_count;
        return a.hexamer_llr > b.hexamer_llr;  // Higher LLR is better
    };

    std::sort(fwd_frames.begin(), fwd_frames.end(), compare);
    std::sort(rev_frames.begin(), rev_frames.end(), compare);

    // Build result FrameScores
    FrameScore best_fwd, best_rev;

    // Best forward
    const auto& bf = fwd_frames[0];
    best_fwd.frame = bf.frame;
    best_fwd.forward = true;
    best_fwd.protein = bf.protein;
    best_fwd.stop_codon_penalty = bf.stop_penalty;
    best_fwd.codon_score = calculate_codon_usage_score_direct(seq, bf.frame);
    best_fwd.aa_composition_score = calculate_aa_composition_score(bf.protein);
    best_fwd.damage_consistency = 0.5f;
    // Composite score: primarily stop penalty, with hexamer bonus
    // Normalize hexamer LLR to 0-0.3 range for small bonus
    float hexamer_bonus_fwd = std::max(0.0f, std::min(0.3f, bf.hexamer_llr * 0.05f + 0.15f));
    best_fwd.total_score = bf.stop_penalty * 0.7f + hexamer_bonus_fwd;

    // Best reverse
    const auto& br = rev_frames[0];
    best_rev.frame = br.frame;
    best_rev.forward = false;
    best_rev.protein = br.protein;
    best_rev.stop_codon_penalty = br.stop_penalty;
    best_rev.codon_score = calculate_codon_usage_score_direct(rc_seq, br.frame);
    best_rev.aa_composition_score = calculate_aa_composition_score(br.protein);
    best_rev.damage_consistency = 0.5f;
    float hexamer_bonus_rev = std::max(0.0f, std::min(0.3f, br.hexamer_llr * 0.05f + 0.15f));
    best_rev.total_score = br.stop_penalty * 0.7f + hexamer_bonus_rev;

    return {best_fwd, best_rev};
}

std::vector<FrameScore> FrameSelector::score_all_frames_full(
    const std::string& seq,
    const DamageProfile* damage) {

    std::vector<FrameScore> scores;
    scores.reserve(6);

    if (seq.length() < 3) {
        return scores;
    }

    // Compute reverse complement ONCE
    const std::string& rc_seq = reverse_complement_fast(seq);

    std::string protein;

    // Score all 6 frames
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        FrameScore fwd;
        fwd.frame = frame;
        fwd.forward = true;
        fwd.total_score = score_frame_fast(seq, frame, protein);
        fwd.protein = protein;
        if (fwd.total_score > -100.0f) {
            fwd.codon_score = calculate_codon_usage_score_direct(seq, frame);
            fwd.stop_codon_penalty = calculate_stop_penalty(protein);
            fwd.aa_composition_score = calculate_aa_composition_score(protein);
            fwd.damage_consistency = score_damage_frame_consistency(seq, frame);
        }
        scores.push_back(std::move(fwd));

        // Reverse frame
        FrameScore rev;
        rev.frame = frame;
        rev.forward = false;
        rev.total_score = score_frame_fast(rc_seq, frame, protein);
        rev.protein = protein;
        if (rev.total_score > -100.0f) {
            rev.codon_score = calculate_codon_usage_score_direct(rc_seq, frame);
            rev.stop_codon_penalty = calculate_stop_penalty(protein);
            rev.aa_composition_score = calculate_aa_composition_score(protein);
            rev.damage_consistency = score_damage_frame_consistency(rc_seq, frame);
        }
        scores.push_back(std::move(rev));
    }

    // Sort by score descending
    std::sort(scores.begin(), scores.end(),
              [](const FrameScore& a, const FrameScore& b) {
                  return a.total_score > b.total_score;
              });

    return scores;
}

std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand_hexamer_aware(
    const std::string& seq,
    size_t terminal_length) {

    if (seq.length() < 3) {
        FrameScore empty;
        empty.total_score = -1000.0f;
        return {empty, empty};
    }

    // Compute reverse complement ONCE
    const std::string& rc_seq = reverse_complement_fast(seq);

    FrameScore best_fwd, best_rev;
    best_fwd.total_score = -1000.0f;
    best_fwd.forward = true;
    best_rev.total_score = -1000.0f;
    best_rev.forward = false;

    std::string protein;

    // Score all 6 frames with hexamer-aware scoring
    for (int frame = 0; frame < 3; ++frame) {
        // Forward frame
        protein = translate_sequence(seq, frame);
        if (!protein.empty()) {
            // Core scores
            float codon = calculate_codon_usage_score_direct(seq, frame);
            float aa_comp = calculate_aa_composition_score(protein);
            float dicodon = calculate_dicodon_score(seq, frame);
            float dipeptide = calculate_dipeptide_score(protein);

            // Hexamer-aware stop penalty
            float stop = calculate_stop_penalty_hexamer_aware(seq, protein, frame, terminal_length);

            // Hexamer log-likelihood ratio score (damage-corrected)
            float hexamer_llr = calculate_hexamer_llr_score(seq, frame);
            // Normalize to 0-1 range (typical LLR range: -2 to +4)
            float hexamer_score = 0.5f + 0.1f * hexamer_llr;
            hexamer_score = std::max(0.0f, std::min(1.0f, hexamer_score));

            // GC3 bonus
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

            float len_bonus = std::min(0.2f, static_cast<float>(protein.length()) / 100.0f);

            // Damage consistency from base pattern
            float damage_consistency = score_damage_frame_consistency(seq, frame);
            float damage_bonus = 0.0f;
            if (seq.length() > 0) {
                bool has_damage_signal = (fast_upper(seq[0]) == 'T') ||
                                          (fast_upper(seq[seq.length()-1]) == 'A');
                if (has_damage_signal) {
                    damage_bonus = (damage_consistency - 0.5f) * 0.10f;
                }
            }

            // Combined score with hexamer evidence
            float total = 0.12f * codon + 0.25f * stop + 0.08f * aa_comp +
                         0.15f * hexamer_score +  // New: hexamer LLR contributes to score
                         0.10f * dicodon + 0.10f * dipeptide +
                         0.05f * gc3_bonus + 0.05f * len_bonus +
                         damage_bonus;

            if (total > best_fwd.total_score) {
                best_fwd.total_score = total;
                best_fwd.frame = frame;
                best_fwd.forward = true;
                best_fwd.codon_score = codon;
                best_fwd.stop_codon_penalty = stop;
                best_fwd.aa_composition_score = aa_comp;
                best_fwd.damage_consistency = damage_consistency;
                best_fwd.protein = std::move(protein);
            }
        }

        // Reverse frame
        protein = translate_sequence(rc_seq, frame);
        if (!protein.empty()) {
            // Core scores
            float codon = calculate_codon_usage_score_direct(rc_seq, frame);
            float aa_comp = calculate_aa_composition_score(protein);
            float dicodon = calculate_dicodon_score(rc_seq, frame);
            float dipeptide = calculate_dipeptide_score(protein);

            // Hexamer-aware stop penalty
            float stop = calculate_stop_penalty_hexamer_aware(rc_seq, protein, frame, terminal_length);

            // Hexamer log-likelihood ratio score (damage-corrected)
            float hexamer_llr = calculate_hexamer_llr_score(rc_seq, frame);
            float hexamer_score = 0.5f + 0.1f * hexamer_llr;
            hexamer_score = std::max(0.0f, std::min(1.0f, hexamer_score));

            // GC3 bonus
            int gc3_count = 0, total3 = 0;
            for (size_t i = frame + 2; i < rc_seq.length(); i += 3) {
                char c = fast_upper(rc_seq[i]);
                if (c != 'N') {
                    if (c == 'G' || c == 'C') gc3_count++;
                    total3++;
                }
            }
            float gc3 = total3 > 0 ? static_cast<float>(gc3_count) / total3 : 0.5f;
            float gc3_bonus = (gc3 > 0.5f) ? 0.1f : 0.0f;

            float len_bonus = std::min(0.2f, static_cast<float>(protein.length()) / 100.0f);

            // Damage consistency from base pattern
            float damage_consistency = score_damage_frame_consistency(rc_seq, frame);
            float damage_bonus = 0.0f;
            if (rc_seq.length() > 0) {
                bool has_damage_signal = (fast_upper(rc_seq[0]) == 'T') ||
                                          (fast_upper(rc_seq[rc_seq.length()-1]) == 'A');
                if (has_damage_signal) {
                    damage_bonus = (damage_consistency - 0.5f) * 0.10f;
                }
            }

            // Combined score with hexamer evidence
            float total = 0.12f * codon + 0.25f * stop + 0.08f * aa_comp +
                         0.15f * hexamer_score +  // New: hexamer LLR contributes to score
                         0.10f * dicodon + 0.10f * dipeptide +
                         0.05f * gc3_bonus + 0.05f * len_bonus +
                         damage_bonus;

            if (total > best_rev.total_score) {
                best_rev.total_score = total;
                best_rev.frame = frame;
                best_rev.forward = false;
                best_rev.codon_score = codon;
                best_rev.stop_codon_penalty = stop;
                best_rev.aa_composition_score = aa_comp;
                best_rev.damage_consistency = damage_consistency;
                best_rev.protein = std::move(protein);
            }
        }
    }

    return {best_fwd, best_rev};
}

} // namespace agp
