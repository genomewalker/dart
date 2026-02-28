/**
 * Frame selection for ancient DNA reads - UNIFIED CODON MARGINALIZATION
 *
 * Joint probabilistic model for frame, strand, and damage:
 * P(frame, strand | seq, damage) via codon-level marginalization
 *
 * For each observed codon, marginalizes over all 64 true codons:
 *   Σ_c P(c) * P(observed | c, damage_at_position)
 *
 * This naturally handles:
 * - Damage-induced stops (high P(CAA→TAA) near 5' explains observed TAA)
 * - Strand discrimination (damage zones map to different codon positions)
 * - Damage correction (MAP codon differs from observed when damage likely)
 */

#include "dart/frame_selector.hpp"
#include "dart/unified_codon_scorer.hpp"
#include "dart/codon_tables.hpp"
#include "dart/hexamer_tables.hpp"
#include "dart/damage_aa_table.hpp"
#include "dart/damage_model.hpp"
#include "dart/damage_context_tables.hpp"
#include <cmath>
#include <algorithm>

namespace dart {

// Convert codon scorer result to FrameScore
static FrameScore codon_to_frame_score(const codon::FrameStrandResult& cr) {
    FrameScore fs;
    fs.frame = cr.frame;
    fs.forward = cr.forward;

    // Use posterior directly as total_score
    fs.total_score = cr.posterior;

    fs.protein = cr.protein;
    fs.corrected_protein = cr.corrected_protein;
    fs.log_likelihood = cr.log_likelihood;
    fs.stop_codon_penalty = 0.0f;  // Handled via marginalization
    fs.aa_composition_score = 0.0f;
    fs.damage_consistency = cr.posterior;
    return fs;
}

// Create codon damage params from sample profile
static codon::DamageParams make_codon_damage_params(const SampleDamageProfile* sample) {
    codon::DamageParams dmg;
    if (sample) {
        dmg.d_max_5p = sample->d_max_5prime;
        dmg.d_max_3p = sample->d_max_3prime;
        dmg.lambda_5p = sample->lambda_5prime;
        dmg.lambda_3p = sample->lambda_3prime;
    }
    return dmg;
}

static codon::DamageParams make_codon_damage_params(const DamageProfile* profile) {
    codon::DamageParams dmg;
    if (profile) {
        dmg.d_max_5p = profile->delta_max;
        dmg.d_max_3p = profile->delta_max;
        dmg.lambda_5p = profile->lambda_5prime;
        dmg.lambda_3p = profile->lambda_3prime;
    }
    return dmg;
}

FrameScore FrameSelector::score_frame(
    const std::string& seq,
    int frame,
    bool forward,
    const DamageProfile* damage) {

    auto dmg = make_codon_damage_params(damage);
    auto results = codon::score_all_hypotheses(seq, dmg);

    // Find the matching frame
    for (const auto& r : results) {
        if (r.frame == frame && r.forward == forward) {
            return codon_to_frame_score(r);
        }
    }

    // Fallback
    FrameScore fs;
    fs.frame = frame;
    fs.forward = forward;
    fs.total_score = 0.0f;
    return fs;
}

FrameScore FrameSelector::select_best_frame(
    const std::string& seq,
    const DamageProfile* damage) {

    auto dmg = make_codon_damage_params(damage);
    auto results = codon::score_all_hypotheses(seq, dmg);

    // Results are sorted by posterior, best first
    return codon_to_frame_score(results[0]);
}

std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand(
    const std::string& seq,
    const DamageProfile* damage) {

    auto dmg = make_codon_damage_params(damage);
    auto [best_fwd, best_rev] = codon::select_best_per_strand(seq, dmg);

    return {codon_to_frame_score(best_fwd), codon_to_frame_score(best_rev)};
}

std::pair<FrameScore, FrameScore> FrameSelector::select_best_per_strand(
    const std::string& seq,
    const DamageProfile* damage,
    const SampleDamageProfile* sample_profile) {

    // Use sample profile if available
    auto dmg = sample_profile ? make_codon_damage_params(sample_profile) : make_codon_damage_params(damage);

    // If damage not validated by Channel B, disable damage model
    if (sample_profile && !sample_profile->damage_validated) {
        dmg.d_max_5p = 0.0f;
        dmg.d_max_3p = 0.0f;
    }

    auto [best_fwd, best_rev] = codon::select_best_per_strand(seq, dmg);

    return {codon_to_frame_score(best_fwd), codon_to_frame_score(best_rev)};
}

std::vector<FrameScore> FrameSelector::score_all_frames(
    const std::string& seq,
    const DamageProfile* damage) {

    auto dmg = make_codon_damage_params(damage);
    // Request all proteins since caller needs all 6 frames
    auto results = codon::score_all_hypotheses(seq, dmg, codon::ScoringWeights(), true);

    std::vector<FrameScore> scores;
    scores.reserve(6);
    for (const auto& r : results) {
        scores.push_back(codon_to_frame_score(r));
    }
    return scores;
}

float FrameSelector::calculate_strand_confidence(
    const std::string& seq,
    float fwd_score,
    float rev_score,
    int fwd_frame,
    int rev_frame,
    const SampleDamageProfile& sample_profile) {

    (void)seq;
    (void)fwd_frame;
    (void)rev_frame;
    (void)sample_profile;

    // Confidence from posterior difference
    float score_diff = std::abs(fwd_score - rev_score);
    return std::clamp(0.5f + score_diff * 0.5f, 0.5f, 1.0f);
}

// Damage strand preference calculation
float calculate_damage_strand_preference(
    const std::string& seq,
    const std::string& rc_seq,
    const SampleDamageProfile& profile) {

    if (profile.max_damage_5prime < 0.05f) return 0.0f;

    float fwd_signal = FrameSelector::estimate_damage_signal(seq);
    float rev_signal = FrameSelector::estimate_damage_signal(rc_seq);

    return fwd_signal - rev_signal;
}

// Damage probability at a nucleotide position using exponential decay model
static float compute_position_damage_prob(
    size_t nt_pos,
    size_t seq_len,
    const SampleDamageProfile& sample_profile,
    bool is_5prime_damage) {  // true for C->T (5'), false for G->A (3')

    // Tri-state damage validation: VALIDATED=full, CONTRADICTED=suppress, UNVALIDATED=soft
    const auto state = get_damage_validation_state(sample_profile);
    if (state == DamageValidationState::CONTRADICTED) {
        return 0.0f;
    }
    // For UNVALIDATED state with low d_max, also suppress
    const float d_max_check = std::max(sample_profile.d_max_5prime, sample_profile.d_max_3prime);
    if (d_max_check < 0.05f && state == DamageValidationState::UNVALIDATED) {
        return 0.0f;
    }

    float d_max, lambda;
    size_t dist;

    if (is_5prime_damage) {
        // C->T damage: strongest at 5' end
        d_max = sample_profile.d_max_5prime;
        lambda = sample_profile.lambda_5prime;
        dist = nt_pos;  // Distance from 5' end
    } else {
        // G->A damage: strongest at 3' end
        d_max = sample_profile.d_max_3prime;
        lambda = sample_profile.lambda_3prime;
        dist = (seq_len > nt_pos) ? seq_len - 1 - nt_pos : 0;  // Distance from 3' end
    }

    // For single-stranded libraries, C->T can happen at both ends
    if (sample_profile.library_type == SampleDamageProfile::LibraryType::SINGLE_STRANDED) {
        // Use minimum distance to either terminal for C->T
        if (is_5prime_damage) {
            size_t dist_3p = (seq_len > nt_pos) ? seq_len - 1 - nt_pos : 0;
            dist = std::min(dist, dist_3p);
            d_max = std::max(sample_profile.d_max_5prime, sample_profile.d_max_3prime);
        }
    }

    if (d_max < 0.01f) {
        return 0.0f;
    }

    float p_damage = d_max * std::exp(-lambda * static_cast<float>(dist));
    return std::clamp(p_damage, 0.0f, 1.0f);
}

float FrameSelector::compute_stop_damage_probability(
    const char* codon,
    size_t codon_nt_start,
    size_t seq_len,
    const SampleDamageProfile& sample_profile) {

    char c0 = fast_upper(codon[0]);
    char c1 = fast_upper(codon[1]);
    char c2 = fast_upper(codon[2]);

    // Check for C->T convertible stops (5' damage)
    // TAA <- CAA (Gln), TAG <- CAG (Gln), TGA <- CGA (Arg)
    bool is_ct_convertible = (c0 == 'T' && c1 == 'A' && c2 == 'A') ||  // TAA from CAA
                             (c0 == 'T' && c1 == 'A' && c2 == 'G') ||  // TAG from CAG
                             (c0 == 'T' && c1 == 'G' && c2 == 'A');    // TGA from CGA

    // Check for G->A convertible stops (3' damage)
    // TGA <- TGG (Trp) via G->A at position 3
    // TAG <- TGG (Trp) via G->A at position 2 (less common pattern)
    bool is_ga_convertible_pos3 = (c0 == 'T' && c1 == 'G' && c2 == 'A');  // TGA from TGG
    bool is_ga_convertible_pos2 = (c0 == 'T' && c1 == 'A' && c2 == 'G');  // TAG from TGG

    float p_damage = 0.0f;

    if (is_ct_convertible) {
        // C->T at position 0 of codon
        p_damage = compute_position_damage_prob(codon_nt_start, seq_len, sample_profile, true);
    }

    if (is_ga_convertible_pos3) {
        // G->A at position 2 of codon (TGG->TGA)
        float p_ga = compute_position_damage_prob(codon_nt_start + 2, seq_len, sample_profile, false);
        p_damage = std::max(p_damage, p_ga);
    }

    if (is_ga_convertible_pos2 && !is_ct_convertible) {
        // G->A at position 1 of codon (TGG->TAG) - only if not already C->T convertible
        float p_ga = compute_position_damage_prob(codon_nt_start + 1, seq_len, sample_profile, false);
        p_damage = std::max(p_damage, p_ga);
    }

    return p_damage;
}

// Check if a stop codon could be from C→T damage
static bool is_damage_convertible_stop(char c0, char c1, char c2) {
    return (c0 == 'T' && c1 == 'A' && c2 == 'A') ||  // TAA ← CAA
           (c0 == 'T' && c1 == 'A' && c2 == 'G') ||  // TAG ← CAG
           (c0 == 'T' && c1 == 'G' && c2 == 'A');    // TGA ← CGA or TGG
}

// Infer the original amino acid from a damage-convertible stop codon
// Returns the most likely pre-damage AA based on codon, position, and hexamer context
// - TAA: Q (CAA→TAA via C→T)
// - TAG: Q (CAG→TAG via C→T)
// - TGA: R or W based on position + hexamer context
//
// For TGA disambiguation, we use hexamer frequencies to compare:
// - [prev_codon + CGA] → R (Arginine)
// - [prev_codon + TGG] → W (Tryptophan)
// The higher scoring hexamer indicates the more likely original codon.
static char infer_damaged_stop_aa(char c0, char c1, char c2,
                                   size_t codon_pos, size_t seq_len,
                                   bool is_forward,
                                   const char* prev_codon = nullptr) {
    // TAA and TAG are unambiguous: both from C→T at position 1
    if (c0 == 'T' && c1 == 'A' && c2 == 'A') return 'Q';  // CAA→TAA
    if (c0 == 'T' && c1 == 'A' && c2 == 'G') return 'Q';  // CAG→TAG

    // TGA is ambiguous: could be CGA→TGA (R) or TGG→TGA (W)
    if (c0 == 'T' && c1 == 'G' && c2 == 'A') {
        // Calculate distance from terminals
        size_t dist_5prime = is_forward ? codon_pos : (seq_len - codon_pos - 3);
        size_t dist_3prime = is_forward ? (seq_len - codon_pos - 3) : codon_pos;

        // Damage zone is ~15-20 nt from terminus
        constexpr size_t DAMAGE_ZONE = 20;

        // Clear cases: only one terminal in damage zone
        if (dist_5prime < DAMAGE_ZONE && dist_3prime >= DAMAGE_ZONE) {
            return 'R';  // 5' terminal: C→T damage → CGA→TGA
        } else if (dist_3prime < DAMAGE_ZONE && dist_5prime >= DAMAGE_ZONE) {
            return 'W';  // 3' terminal: G→A damage → TGG→TGA
        }

        // Ambiguous case: use hexamer context if available
        if (prev_codon != nullptr) {
            // Construct hexamers: [prev_codon + hypothetical_original_codon]
            char hex_cga[7] = {prev_codon[0], prev_codon[1], prev_codon[2], 'C', 'G', 'A', '\0'};
            char hex_tgg[7] = {prev_codon[0], prev_codon[1], prev_codon[2], 'T', 'G', 'G', '\0'};

            // Look up hexamer frequencies
            uint32_t code_cga = encode_hexamer(hex_cga);
            uint32_t code_tgg = encode_hexamer(hex_tgg);

            if (code_cga != UINT32_MAX && code_tgg != UINT32_MAX) {
                float freq_cga = GTDB_HEXAMER_FREQ[code_cga];
                float freq_tgg = GTDB_HEXAMER_FREQ[code_tgg];

                // Return the amino acid corresponding to higher frequency hexamer
                if (freq_tgg > freq_cga * 1.5f) {
                    return 'W';  // TGG significantly more common in this context
                }
                // Otherwise default to R (CGA is generally more common than TGG)
            }
        }

        // Default to R (Arginine) - more common in coding sequences
        return 'R';
    }

    // Fallback (shouldn't reach here if is_damage_convertible_stop was true)
    return 'X';
}

// Generate search-optimized protein by extending past damage-convertible stops.
// This improves MMseqs2 alignment sensitivity for ancient DNA.
static std::string generate_search_protein(
    const std::string& protein,
    const std::string& oriented_seq,
    int frame,
    size_t orf_nt_start,
    size_t orf_aa_end,  // Codon position where ORF ended (stop codon position)
    float d_max,
    const SampleDamageProfile* sample_profile = nullptr,
    bool is_forward = true,
    float orf_coding_score = 0.0f) {  // Average hexamer score for the ORF

    // If no significant damage, return observed protein as-is
    if (d_max < 0.05f) {
        return protein;
    }

    std::string search_protein = protein;

    // If the ORF ended at a stop codon (not at sequence end), check if we should extend
    // Only extend if the ORF is short (early stop) and the stop is damage-convertible
    constexpr size_t EARLY_STOP_THRESHOLD_AA = 5;  // Only extend if stop was in first 5 AAs
    constexpr size_t MAX_EXTENSION = 3;  // Maximum stops to extend through

    if (protein.length() < EARLY_STOP_THRESHOLD_AA) {
        // ORF was terminated early - try to extend through damage-convertible stops
        size_t extensions = 0;
        size_t current_codon = orf_aa_end;  // Start at the stop codon that ended the ORF

        while (extensions < MAX_EXTENSION) {
            size_t codon_nt = frame + current_codon * 3;
            if (codon_nt + 2 >= oriented_seq.length()) break;

            char c0 = fast_upper(oriented_seq[codon_nt]);
            char c1 = fast_upper(oriented_seq[codon_nt + 1]);
            char c2 = fast_upper(oriented_seq[codon_nt + 2]);

            bool is_stop = (c0 == 'T' && c1 == 'A' && (c2 == 'A' || c2 == 'G')) ||
                           (c0 == 'T' && c1 == 'G' && c2 == 'A');

            if (is_stop && is_damage_convertible_stop(c0, c1, c2)) {
                // Extend with inferred original AA instead of truncating
                // Pass previous codon for hexamer-based TGA disambiguation
                const char* prev = (codon_nt >= 3) ? (oriented_seq.data() + codon_nt - 3) : nullptr;
                char inferred = infer_damaged_stop_aa(c0, c1, c2, codon_nt,
                                                       oriented_seq.length(), is_forward, prev);
                search_protein += inferred;
                extensions++;
                current_codon++;
            } else if (is_stop) {
                // Non-convertible stop - terminate
                break;
            } else {
                // Regular codon - translate and add
                char aa = CodonTable::translate_codon(c0, c1, c2);
                search_protein += aa;
                current_codon++;
            }
        }
    }

    // We intentionally do NOT modify non-stop amino acids in search_protein.
    //
    // All damage-induced substitutions (W←R, C←R, Y←H, I←T, V←A, K←E, N←D, etc.)
    // are kept as-observed because:
    //
    // 1. X-masking breaks MMseqs2 prefilter k-mer matching → hurts sensitivity
    // 2. Correction (e.g., W→R) loses damage evidence → can't detect post-mapping
    // 3. Keeping amino acids as-is:
    //    - Prefilter k-mers work normally
    //    - Alignment shows mismatches (e.g., W vs R)
    //    - damage-annotate detects all patterns via DAMAGE_SUBS table
    //
    // The ONLY modification is X-masking of damage-convertible STOPS (above),
    // because stops truncate proteins and X allows alignment through.

    return search_protein;
}


float FrameSelector::compute_per_read_damage_prior(
    const std::string& seq,
    const SampleDamageProfile& sample_profile) {

    constexpr size_t WINDOW = 20;
    constexpr float SMOOTHING = 0.5f;
    constexpr float KAPPA = 20.0f;  // Prior concentration
    constexpr float EPS = 1e-6f;

    const size_t seq_len = seq.length();
    if (seq_len < 60) {
        return 0.0f;
    }

    // Use d_max_combined which accounts for Channel A/B arbitration
    // This respects artifact detection and Channel B fallback when Channel A is unreliable
    const float d_max = sample_profile.d_max_combined > 0.0f
        ? sample_profile.d_max_combined
        : std::max(sample_profile.d_max_5prime, sample_profile.d_max_3prime);
    if (d_max < 0.01f) {
        return 0.0f;
    }

    const float lambda = (sample_profile.lambda_5prime + sample_profile.lambda_3prime) * 0.5f;

    // Compute interior base frequencies
    int t_int = 0, c_int = 0, a_int = 0, g_int = 0;
    for (size_t i = WINDOW; i < seq_len - WINDOW; ++i) {
        char b = fast_upper(seq[i]);
        switch (b) {
            case 'T': ++t_int; break;
            case 'C': ++c_int; break;
            case 'A': ++a_int; break;
            case 'G': ++g_int; break;
        }
    }

    float n = static_cast<float>(seq_len - 2 * WINDOW);
    float p_T = (t_int + SMOOTHING) / (n + 4 * SMOOTHING);
    float p_C = (c_int + SMOOTHING) / (n + 4 * SMOOTHING);
    float p_A = (a_int + SMOOTHING) / (n + 4 * SMOOTHING);
    float p_G = (g_int + SMOOTHING) / (n + 4 * SMOOTHING);

    // Baseline product probabilities (P(T|C/T, no damage), P(A|A/G, no damage))
    float b5 = std::clamp((p_T + 1e-3f) / (p_T + p_C + 2e-3f), 0.01f, 0.99f);
    float b3 = std::clamp((p_A + 1e-3f) / (p_A + p_G + 2e-3f), 0.01f, 0.99f);

    // Collect informative sites: y=1 if product base (T at 5', A at 3')
    struct SiteObs { float b, c; uint8_t y; };
    std::array<SiteObs, 64> sites;
    size_t m = 0;

    // Avoid 5'/3' overlap for short reads
    const size_t w = std::min(WINDOW, seq_len / 2);

    auto push_site = [&](float b, float dist, bool is_product) {
        if (m >= sites.size()) return;
        const float c = (1.0f - b) * std::exp(-lambda * dist);
        if (c < 1e-7f) return;
        sites[m++] = SiteObs{b, c, static_cast<uint8_t>(is_product ? 1 : 0)};
    };

    // Skip position 0 if sample shows adapter artifact (pos0 depleted, pos1 enriched)
    // In artifact samples, position 0 measures adapter ligation bias, not damage
    const size_t skip_5 = sample_profile.position_0_artifact_5prime ? 1 : 0;
    const size_t skip_3 = sample_profile.position_0_artifact_3prime ? 1 : 0;

    // 5' channel: C/T sites
    for (size_t i = skip_5; i < w; ++i) {
        char x = fast_upper(seq[i]);
        // When skipping, adjust effective position for decay calculation
        const float eff_pos = static_cast<float>(i - skip_5);
        if (x == 'C') push_site(b5, eff_pos, false);
        else if (x == 'T') push_site(b5, eff_pos, true);
    }

    // 3' channel: G/A sites
    for (size_t i = skip_3; i < w; ++i) {
        char x = fast_upper(seq[seq_len - 1 - i]);
        const float eff_pos = static_cast<float>(i - skip_3);
        if (x == 'G') push_site(b3, eff_pos, false);
        else if (x == 'A') push_site(b3, eff_pos, true);
    }

    if (m == 0) {
        return d_max * 0.2f;  // Uninformative read, return weak prior
    }

    // ---- MAP-Newton for θ (2 iterations) ----
    // Prior: θ ~ Beta(α, β) with mean = d_max
    float alpha = std::max(1e-3f, d_max * KAPPA);
    float beta_prior = std::max(1e-3f, (1.0f - d_max) * KAPPA);
    if (alpha <= 1.0f || beta_prior <= 1.0f) {
        alpha += 1.0f;
        beta_prior += 1.0f;
    }

    float theta = std::clamp(d_max, 1e-4f, 1.0f - 1e-4f);

    for (int iter = 0; iter < 2; ++iter) {
        float g = (alpha - 1.0f) / theta - (beta_prior - 1.0f) / (1.0f - theta);
        float h = -(alpha - 1.0f) / (theta * theta)
                  - (beta_prior - 1.0f) / ((1.0f - theta) * (1.0f - theta));

        for (size_t k = 0; k < m; ++k) {
            const float b = sites[k].b;
            const float c = sites[k].c;
            const float y = static_cast<float>(sites[k].y);

            const float p1 = std::max(EPS, b + c * theta);        // P(y=1|θ)
            const float p0 = std::max(EPS, 1.0f - b - c * theta); // P(y=0|θ)

            g += y * (c / p1) - (1.0f - y) * (c / p0);
            h -= y * (c * c) / (p1 * p1) + (1.0f - y) * (c * c) / (p0 * p0);
        }

        if (std::fabs(h) < 1e-8f) break;
        float step = std::clamp(g / h, -0.25f, 0.25f);
        theta = std::clamp(theta - step, 1e-5f, 1.0f - 1e-5f);
    }

    // ---- Compute P(read damaged) via opportunity-weighted model ----
    float C_eff = 0.0f;
    float log_p_no_damage = 0.0f;

    for (size_t k = 0; k < m; ++k) {
        const float c = sites[k].c;
        C_eff += c;
        const float qk = std::clamp(theta * c, 0.0f, 0.999999f);
        log_p_no_damage += std::log1p(-qk);
    }

    float p_read_damaged = 1.0f - std::exp(log_p_no_damage);

    // Reliability shrinkage for low-information reads
    float reliability = C_eff / (C_eff + 1.0f);
    p_read_damaged = reliability * p_read_damaged + (1.0f - reliability) * d_max * 0.2f;

    return std::clamp(p_read_damaged, 0.0f, 1.0f);
}


std::vector<FrameSelector::ORFFragment> FrameSelector::enumerate_orf_fragments(
    const std::string& seq,
    const SampleDamageProfile& sample_profile,
    size_t min_aa,
    bool adaptive,
    float per_read_damage) {

    // =========================================================================
    // UNIFIED DAMAGE-AWARE ORF ENUMERATION
    //
    // Architecture:
    // 1. Use pre-computed buffers from score_all_hypotheses():
    //    - protein_buffer[h][c]: observed amino acid at codon c in hypothesis h
    //    - corrected_buffer[h][c]: MAP-inferred amino acid (damage-corrected)
    //    - codon_scores[h][c]: hexamer/coding log-likelihood score
    //
    // 2. Stop classification:
    //    - Damage-induced stop (p_damage >= threshold): continue ORF with inferred AA
    //    - Real stop (p_damage < threshold): split ORF, evidence of wrong frame
    //
    // 3. Frame scoring: S_h = Σ(hexamer_scores) - w_real × real_stops
    //    Real stops penalize the hypothesis as they indicate incorrect frame/strand
    //
    // 4. Per-read damage integration (NEW):
    //    - High per-read damage → more lenient stop penalty, prefer longer X-masked ORFs
    //    - Low per-read damage → standard stop penalty
    // =========================================================================

    constexpr float BASE_STOP_PENALTY = 3.0f;  // Base penalty per real stop

    // Auto-compute per-read damage if not provided
    float p_read_dmg = per_read_damage;
    if (p_read_dmg < 0.0f) {
        p_read_dmg = compute_per_read_damage_prior(seq, sample_profile);
    }

    // Adjust stop penalty based on per-read damage:
    // High damage (p=1.0) → 50% penalty reduction (more lenient)
    // Low damage (p=0.0) → full penalty
    const float REAL_STOP_PENALTY = BASE_STOP_PENALTY * (1.0f - 0.5f * p_read_dmg);

    std::vector<ORFFragment> fwd_fragments;
    std::vector<ORFFragment> rev_fragments;
    fwd_fragments.reserve(18);
    rev_fragments.reserve(18);

    // Compute hexamer/codon scores for all 6 hypotheses
    // This populates the thread-local buffers
    auto dmg = make_codon_damage_params(&sample_profile);
    auto all_scores = codon::score_all_hypotheses(seq, dmg, codon::ScoringWeights(), true);
    auto& buf = codon::get_thread_buffers();

    // STOP RESCUE STRATEGY:
    // AA-level correction (guessing Q, R, W) was disabled due to ~2% precision.
    // Instead: X-mask damage-convertible stops in search_protein (MMseqs2-safe),
    // keep '*' in protein (truthful). This prevents ORF truncation at damage stops
    // without making false AA claims. Post-mapping annotation (dart damage-annotate)
    // achieves ~90% precision using reference proteins.
    // Use d_max_combined which respects Channel A/B arbitration
    const float d_max = sample_profile.d_max_combined > 0.0f
        ? sample_profile.d_max_combined
        : std::max(sample_profile.d_max_5prime, sample_profile.d_max_3prime);
    (void)adaptive;

    // Process both strands
    for (bool is_forward : {true, false}) {
        const std::string& oriented = is_forward ? seq : reverse_complement_cached(seq);
        auto& strand_fragments = is_forward ? fwd_fragments : rev_fragments;

        for (int frame = 0; frame < 3; ++frame) {
            // Hypothesis index: fwd frames 0,1,2 -> h=0,1,2; rev frames 0,1,2 -> h=3,4,5
            int h = is_forward ? frame : (3 + frame);
            size_t num_codons = buf.protein_lens[h];

            if (num_codons == 0) continue;

            // Build ORFs by iterating through codons, splitting at real stops
            std::string current_protein;
            std::string current_search_protein;
            std::string current_observed;
            std::string current_corrected_nt;
            current_protein.reserve(num_codons);
            current_search_protein.reserve(num_codons);
            current_observed.reserve(num_codons);
            current_corrected_nt.reserve(num_codons * 3);

            size_t orf_start = 0;          // Start codon index of current ORF
            size_t frame_real_stops = 0;   // Total real stops in this frame (for scoring)
            size_t orf_passed_stops = 0;   // Damage stops continued through in current ORF (legacy)
            size_t orf_rescued_stops = 0;  // Damage stops X-masked in search_protein
            size_t orf_aa_corrections = 0; // Non-stop AA corrections in current ORF
            float orf_damage_evidence = 0.0f;
            float orf_hexamer_score = 0.0f;

            for (size_t c = 0; c < num_codons; ++c) {
                char observed_aa = buf.protein_buffer[h][c];
                float codon_score = buf.codon_scores[h][c];

                if (observed_aa == '*') {
                    size_t codon_nt_start = frame + c * 3;
                    const char* codon_ptr = oriented.data() + codon_nt_start;

                    char c0 = fast_upper(codon_ptr[0]);
                    char c1 = fast_upper(codon_ptr[1]);
                    char c2 = fast_upper(codon_ptr[2]);

                    // Stop rescue: continue ORF through damage-convertible stops
                    // protein keeps '*' (truthful), search_protein gets inferred original AA
                    // Cap at 2 rescued stops per ORF to avoid noise accumulation
                    //
                    // Per-read damage integration: lower threshold for high-damage reads
                    // d_max threshold: 0.05 (baseline) → 0.02 (high damage, p=1.0)
                    const float x_mask_threshold = 0.05f - 0.03f * p_read_dmg;
                    if (d_max >= x_mask_threshold && orf_rescued_stops < 2 &&
                        codon_nt_start + 2 < oriented.length() &&
                        is_damage_convertible_stop(c0, c1, c2)) {
                        // Infer the original AA from the damage pattern
                        // Pass previous codon for hexamer-based TGA disambiguation
                        const char* prev = (codon_nt_start >= 3) ? (oriented.data() + codon_nt_start - 3) : nullptr;
                        char inferred = infer_damaged_stop_aa(c0, c1, c2, codon_nt_start,
                                                               oriented.length(), is_forward, prev);
                        current_protein += '*';
                        current_search_protein += inferred;
                        current_observed += inferred;
                        current_corrected_nt += c0;
                        current_corrected_nt += c1;
                        current_corrected_nt += c2;
                        orf_rescued_stops++;
                        orf_hexamer_score += codon_score;
                    } else {
                        // Real stop (or rescue cap reached): emit current ORF and start new one
                        frame_real_stops++;

                        size_t orf_nt_start = frame + orf_start * 3;
                        float avg_hexamer_score = current_protein.length() > 0 ?
                            orf_hexamer_score / static_cast<float>(current_protein.length()) : 0.0f;
                        std::string search_protein = generate_search_protein(
                            current_search_protein, oriented, frame, orf_nt_start, c, d_max, &sample_profile, is_forward, avg_hexamer_score);

                        size_t effective_length = std::max(current_protein.length(), search_protein.length());

                        if (effective_length >= min_aa) {
                            ORFFragment orf;
                            orf.search_protein = std::move(search_protein);
                            orf.protein = std::move(current_protein);
                            orf.observed_protein = std::move(current_observed);
                            orf.corrected_nt = std::move(current_corrected_nt);
                            orf.aa_start = orf_start;
                            orf.aa_end = c;
                            orf.nt_start = orf_nt_start;
                            orf.nt_end = frame + c * 3;
                            orf.is_forward = is_forward;
                            orf.frame = frame;
                            orf.length = orf.protein.length();
                            orf.x_count = std::count(orf.search_protein.begin(), orf.search_protein.end(), 'X');
                            orf.aa_corrections = orf_aa_corrections;
                            orf.passed_stops = orf_passed_stops;
                            orf.rescued_stops = orf_rescued_stops;
                            orf.real_stops = frame_real_stops;
                            orf.damage_evidence = orf_damage_evidence;
                            float avg_hexamer = orf.length > 0 ? orf_hexamer_score / static_cast<float>(orf.length) : 0.0f;
                            float length_bonus = std::log(static_cast<float>(orf.length) + 1.0f);
                            // Rescued stop bonus: reward X-masked stops on high-damage reads
                            // High damage (p=1.0) → +0.5 per rescued stop
                            // Low damage (p=0.0) → no bonus
                            float rescue_bonus = p_read_dmg * 0.5f * static_cast<float>(orf_rescued_stops);
                            orf.score = avg_hexamer + length_bonus + rescue_bonus - REAL_STOP_PENALTY * frame_real_stops;
                            strand_fragments.push_back(std::move(orf));
                        }

                        // Reset for next ORF
                        current_protein.clear();
                        current_search_protein.clear();
                        current_observed.clear();
                        current_corrected_nt.clear();
                        orf_start = c + 1;
                        orf_passed_stops = 0;
                        orf_rescued_stops = 0;
                        orf_aa_corrections = 0;
                        orf_damage_evidence = 0.0f;
                        orf_hexamer_score = 0.0f;
                    }
                } else {
                    // Non-stop codon
                    size_t codon_nt_start = frame + c * 3;
                    char c0 = fast_upper(oriented[codon_nt_start]);
                    char c1 = fast_upper(oriented[codon_nt_start + 1]);
                    char c2 = fast_upper(oriented[codon_nt_start + 2]);

                    current_corrected_nt += c0;
                    current_corrected_nt += c1;
                    current_corrected_nt += c2;
                    current_protein += observed_aa;
                    current_search_protein += observed_aa;
                    current_observed += observed_aa;
                    orf_hexamer_score += codon_score;
                }
            }

            // Emit final ORF if long enough
            // For final ORF, there's no stop to extend past (ended at sequence end)
            size_t final_orf_nt_start = frame + orf_start * 3;
            float final_avg_hexamer = current_protein.length() > 0 ?
                orf_hexamer_score / static_cast<float>(current_protein.length()) : 0.0f;
            std::string final_search_protein = generate_search_protein(
                current_search_protein, oriented, frame, final_orf_nt_start, num_codons, d_max, &sample_profile, is_forward, final_avg_hexamer);
            size_t final_effective_length = std::max(current_protein.length(), final_search_protein.length());

            if (final_effective_length >= min_aa) {
                ORFFragment orf;
                orf.search_protein = std::move(final_search_protein);
                orf.protein = std::move(current_protein);
                orf.observed_protein = std::move(current_observed);
                orf.corrected_nt = std::move(current_corrected_nt);
                size_t orf_nt_start = final_orf_nt_start;
                orf.aa_start = orf_start;
                orf.aa_end = num_codons;
                orf.nt_start = orf_nt_start;
                orf.nt_end = frame + num_codons * 3;
                orf.is_forward = is_forward;
                orf.frame = frame;
                orf.length = orf.protein.length();
                orf.x_count = std::count(orf.search_protein.begin(), orf.search_protein.end(), 'X');
                orf.aa_corrections = orf_aa_corrections;
                orf.passed_stops = orf_passed_stops;
                orf.rescued_stops = orf_rescued_stops;
                orf.real_stops = frame_real_stops;
                orf.damage_evidence = orf_damage_evidence;
                float avg_hexamer = orf.length > 0 ? orf_hexamer_score / static_cast<float>(orf.length) : 0.0f;
                float length_bonus = std::log(static_cast<float>(orf.length) + 1.0f);
                float rescue_bonus = p_read_dmg * 0.5f * static_cast<float>(orf_rescued_stops);
                orf.score = avg_hexamer + length_bonus + rescue_bonus - REAL_STOP_PENALTY * frame_real_stops;
                strand_fragments.push_back(std::move(orf));
            }
        }
    }

    // Sort by score (coding signal - stop penalties + damage-aware rescue bonus)
    auto by_score = [](const ORFFragment& a, const ORFFragment& b) {
        return a.score > b.score;
    };
    std::sort(fwd_fragments.begin(), fwd_fragments.end(), by_score);
    std::sort(rev_fragments.begin(), rev_fragments.end(), by_score);

    // =========================================================================
    // SELECTION MODES:
    //
    // 1. DEFAULT (6-FRAME MODE): Output best ORF from each frame
    //    - Guarantees frame coverage (never drops a frame due to hexamer score)
    //    - Equivalent to 6-frame translation but with damage correction
    //    - Up to 6 ORFs total (one per frame)
    //
    // 2. ADAPTIVE MODE: Confidence-based selection using posteriors
    //    - Converts scores to posteriors via softmax
    //    - Emits fewer ORFs when model is confident
    //    - Guarantees at least 1 ORF per strand (strand ambiguity in ancient DNA)
    //    - Falls back to 6-frame on low-complexity sequences
    // =========================================================================

    std::vector<ORFFragment> result;
    result.reserve(6);

    // Check if ORF is already in result (linear scan, max 6 elements)
    auto already_added = [&result](const ORFFragment& orf) {
        for (const auto& r : result) {
            if (r.is_forward == orf.is_forward && r.frame == orf.frame && r.aa_start == orf.aa_start)
                return true;
        }
        return false;
    };
    auto frame_index = [](const ORFFragment& orf) {
        return orf.is_forward ? orf.frame : (3 + orf.frame);
    };

    if (adaptive) {
        // =====================================================================
        // ADAPTIVE MODE FOR LARGE-SCALE PROCESSING (billions of reads)
        //
        // Strategy: Best-per-strand (2 ORFs) + margin-gated extras
        //
        // 1. Always emit best forward + best reverse (2 ORFs baseline)
        // 2. Add 3rd ORF if within-strand ambiguity is high:
        //    - If |fwd_best - fwd_2nd| < MARGIN_THRESHOLD, add fwd_2nd
        //    - If |rev_best - rev_2nd| < MARGIN_THRESHOLD, add rev_2nd
        // 3. Low-complexity burst: if sequence is low-complexity (AT-rich),
        //    fall back to 6-frame to avoid hexamer model failures
        //
        // Target: ~2-3 ORFs/read average, ≥70% coverage
        // =====================================================================
        constexpr float MARGIN_THRESHOLD = 2.0f;    // Score margin for ambiguity
        constexpr float LOW_COMPLEXITY_ENTROPY = 1.5f;  // Shannon entropy threshold

        // Compute sequence complexity (Shannon entropy of dinucleotides)
        auto compute_complexity = [](const std::string& s) -> float {
            if (s.length() < 10) return 2.0f;  // Assume complex if too short
            std::array<int, 16> di_counts = {};
            for (size_t i = 0; i + 1 < s.length(); ++i) {
                int b1 = (s[i] == 'A' || s[i] == 'a') ? 0 :
                         (s[i] == 'C' || s[i] == 'c') ? 1 :
                         (s[i] == 'G' || s[i] == 'g') ? 2 :
                         (s[i] == 'T' || s[i] == 't') ? 3 : -1;
                int b2 = (s[i+1] == 'A' || s[i+1] == 'a') ? 0 :
                         (s[i+1] == 'C' || s[i+1] == 'c') ? 1 :
                         (s[i+1] == 'G' || s[i+1] == 'g') ? 2 :
                         (s[i+1] == 'T' || s[i+1] == 't') ? 3 : -1;
                if (b1 >= 0 && b2 >= 0) di_counts[b1 * 4 + b2]++;
            }
            float total = static_cast<float>(s.length() - 1);
            float entropy = 0.0f;
            for (int c : di_counts) {
                if (c > 0) {
                    float p = static_cast<float>(c) / total;
                    entropy -= p * std::log2(p);
                }
            }
            return entropy;  // Max ~4.0 for uniform, low for repetitive
        };

        float complexity = compute_complexity(seq);
        bool low_complexity = (complexity < LOW_COMPLEXITY_ENTROPY);

        if (low_complexity) {
            // LOW-COMPLEXITY FALLBACK: Use 6-frame mode
            // Hexamer models fail on AT-rich/repetitive sequences
            std::array<bool, 6> frame_represented = {false, false, false, false, false, false};
            for (const auto& orf : fwd_fragments) {
                int fi = frame_index(orf);
                if (!frame_represented[fi]) {
                    result.push_back(orf);
                    frame_represented[fi] = true;
                }
            }
            for (const auto& orf : rev_fragments) {
                int fi = frame_index(orf);
                if (!frame_represented[fi]) {
                    result.push_back(orf);
                    frame_represented[fi] = true;
                }
            }
        } else {
            // NORMAL ADAPTIVE: Best-per-strand + margin-gated extras

            // Step 1: Always emit best from each strand (2 ORFs)
            if (!fwd_fragments.empty()) {
                result.push_back(fwd_fragments[0]);
            }
            if (!rev_fragments.empty()) {
                result.push_back(rev_fragments[0]);
            }

            // Step 2: Add extras if within-strand ambiguity is high
            // Check forward strand ambiguity
            if (fwd_fragments.size() >= 2) {
                float fwd_margin = fwd_fragments[0].score - fwd_fragments[1].score;
                if (fwd_margin < MARGIN_THRESHOLD &&
                    !already_added(fwd_fragments[1])) {
                    result.push_back(fwd_fragments[1]);
                }
            }

            // Check reverse strand ambiguity
            if (rev_fragments.size() >= 2) {
                float rev_margin = rev_fragments[0].score - rev_fragments[1].score;
                if (rev_margin < MARGIN_THRESHOLD &&
                    !already_added(rev_fragments[1])) {
                    result.push_back(rev_fragments[1]);
                }
            }
        }
    } else {
        // 6-FRAME MODE (DEFAULT): Best ORF per frame
        // This ensures we never drop a frame solely due to hexamer score
        std::array<bool, 6> frame_represented = {false, false, false, false, false, false};

        for (const auto& orf : fwd_fragments) {
            int fi = frame_index(orf);
            if (!frame_represented[fi]) {
                result.push_back(orf);
                frame_represented[fi] = true;
            }
        }
        for (const auto& orf : rev_fragments) {
            int fi = frame_index(orf);
            if (!frame_represented[fi]) {
                result.push_back(orf);
                frame_represented[fi] = true;
            }
        }
    }

    // Final sort by score
    std::sort(result.begin(), result.end(), by_score);

    return result;
}

float FrameSelector::infer_per_read_aa_damage(
    const std::string& seq,
    const std::vector<ORFFragment>& orfs,
    const SampleDamageProfile& profile,
    PerReadDamageEvidence* evidence) {

    constexpr float RHO_STOP = 20.0f;
    constexpr int DAMAGE_WINDOW_AA = 15;
    constexpr size_t DAMAGE_WINDOW_NT = 45;

    if (orfs.empty() || seq.length() < 30) {
        if (evidence) *evidence = {};
        return 0.0f;
    }

    if (profile.damage_artifact) {
        if (evidence) {
            *evidence = {};
            evidence->p_read_damaged = 0.05f;
        }
        return 0.05f;
    }

    const float d_max = profile.d_max_combined;
    const float lambda_5 = profile.lambda_5prime;
    const float lambda_3 = profile.lambda_3prime > 0 ? profile.lambda_3prime : lambda_5;

    if (d_max < 0.005f) {
        if (evidence) *evidence = {};
        return 0.0f;
    }

    const auto& orf = orfs[0];
    const size_t seq_len = seq.length();
    const std::string& oriented = orf.is_forward ? seq : reverse_complement_cached(seq);
    const size_t n_codons = orf.length;
    const size_t orf_nt_start = orf.nt_start;

    // ================================================================
    // 1c. Bayesian Posterior Damage Score (ρ=0.157 vs baseline ρ=0.077)
    //
    // Key insight: Ask "P(damaged | T observed)" not "count T at terminal"
    // Uses per-read interior composition as baseline (self-normalized).
    //
    // Formula: P(damaged|T) = P(T|damaged) × P(damage) / P(T)
    // where P(T|damaged) = P_C_interior (C became T via damage)
    //
    // Validated: AUC=0.591, ρ=0.157 on 50K reads (103% improvement over baseline)
    // ================================================================
    float bayesian_damage_score = 0.0f;
    {
        constexpr int WINDOW = 20;
        constexpr float SMOOTHING = 0.5f;  // Optimized (was 4.0)

        // Compute interior base frequencies (per-read self-normalization)
        if (seq_len >= 2 * WINDOW + 10) {
            int t_int = 0, c_int = 0, a_int = 0, g_int = 0;
            for (size_t i = WINDOW; i < seq_len - WINDOW; ++i) {
                char b = fast_upper(oriented[i]);
                switch (b) {
                    case 'T': ++t_int; break;
                    case 'C': ++c_int; break;
                    case 'A': ++a_int; break;
                    case 'G': ++g_int; break;
                }
            }

            float n = static_cast<float>(seq_len - 2 * WINDOW);
            float p_T = (t_int + SMOOTHING) / (n + 4 * SMOOTHING);
            float p_C = (c_int + SMOOTHING) / (n + 4 * SMOOTHING);
            float p_A = (a_int + SMOOTHING) / (n + 4 * SMOOTHING);
            float p_G = (g_int + SMOOTHING) / (n + 4 * SMOOTHING);

            // 5' end: compute P(damaged | T) for each T
            for (size_t i = 0; i < static_cast<size_t>(WINDOW) && i < seq_len; ++i) {
                char b = fast_upper(oriented[i]);
                if (b != 'T') continue;

                float p_damage = d_max * std::exp(-lambda_5 * static_cast<float>(i));
                float p_t_from_damage = p_C * p_damage;
                float p_t_natural = p_T * (1.0f - p_damage);
                float p_t_total = p_t_from_damage + p_t_natural;

                if (p_t_total > 1e-8f) {
                    bayesian_damage_score += p_t_from_damage / p_t_total;
                }
            }

            // 3' end: compute P(damaged | A) for each A
            for (size_t i = 0; i < static_cast<size_t>(WINDOW) && i < seq_len; ++i) {
                size_t pos = seq_len - 1 - i;
                if (pos < static_cast<size_t>(WINDOW)) break;

                char b = fast_upper(oriented[pos]);
                if (b != 'A') continue;

                float p_damage = d_max * std::exp(-lambda_3 * static_cast<float>(i));
                float p_a_from_damage = p_G * p_damage;
                float p_a_natural = p_A * (1.0f - p_damage);
                float p_a_total = p_a_from_damage + p_a_natural;

                if (p_a_total > 1e-8f) {
                    bayesian_damage_score += p_a_from_damage / p_a_total;
                }
            }
        } else {
            // Fallback for short reads: simple position-weighted counting
            for (size_t i = 0; i < seq_len; ++i) {
                char c = fast_upper(oriented[i]);
                if (c == 'T') {
                    float p_damage = d_max * std::exp(-lambda_5 * static_cast<float>(i));
                    bayesian_damage_score += p_damage;
                } else if (c == 'A') {
                    float dist_3 = static_cast<float>(seq_len - 1 - i);
                    float p_damage = d_max * std::exp(-lambda_3 * dist_3);
                    bayesian_damage_score += p_damage;
                }
            }
        }
    }

    // ================================================================
    // Score = bayesian_damage_score directly
    //
    // No hand-tuned combination weights. The score IS the position-
    // weighted T/A sum using only sample-inferred d_max and λ.
    // ================================================================
    float score_logit = bayesian_damage_score;

    float p_read_damaged = 1.0f / (1.0f + std::exp(-score_logit));

    // Apply tri-state damage suppression factor
    const float suppression = get_damage_suppression_factor(profile);
    p_read_damaged *= suppression;

    p_read_damaged = std::clamp(p_read_damaged, 0.0f, 1.0f);

    // ================================================================
    // Diagnostic evidence (only computed when caller needs it)
    // ================================================================
    if (evidence) {
        // Hexamer damage LLR at terminal dicodon positions
        float hex_llr_5p = 0.0f;
        float hex_llr_3p = 0.0f;
        uint16_t n_hex_5p = 0;
        uint16_t n_hex_3p = 0;

        // 5' scan: C→T damage channel
        for (size_t ci = 0; ci < n_codons; ++ci) {
            size_t nt = orf_nt_start + ci * 3;
            if (nt + 5 >= oriented.length()) break;
            if (nt >= DAMAGE_WINDOW_NT) break;

            uint32_t code = encode_hexamer(oriented.c_str() + nt);
            if (code >= 4096) continue;

            float p_damage = orf.is_forward
                ? d_max * std::exp(-lambda_5 * static_cast<float>(nt))
                : d_max * std::exp(-lambda_3 * static_cast<float>(nt));
            if (p_damage < 0.005f) break;

            hex_llr_5p += compute_hexamer_damage_llr(code, p_damage, true, GTDB_HEXAMER_FREQ);
            ++n_hex_5p;
        }

        // 3' scan: G→A damage channel
        for (size_t ci = n_codons; ci > 0; --ci) {
            size_t nt = orf_nt_start + (ci - 1) * 3;
            if (nt + 5 >= oriented.length()) continue;

            size_t dist_3 = seq_len - 1 - (nt + 5);
            if (dist_3 >= DAMAGE_WINDOW_NT) continue;

            uint32_t code = encode_hexamer(oriented.c_str() + nt);
            if (code >= 4096) continue;

            float p_damage = orf.is_forward
                ? d_max * std::exp(-lambda_3 * static_cast<float>(dist_3))
                : d_max * std::exp(-lambda_5 * static_cast<float>(dist_3));
            if (p_damage < 0.005f) break;

            hex_llr_3p += compute_hexamer_damage_llr(code, p_damage, false, GTDB_HEXAMER_FREQ);
            ++n_hex_3p;
        }

        uint16_t n_informative = n_hex_5p + n_hex_3p;

        // Stop evidence: convertible stops in selected ORF
        float logbf_stop = 0.0f;
        uint8_t n_conv_stops = 0;

        for (size_t ci = 0; ci < n_codons; ++ci) {
            size_t nt_pos = orf.nt_start + ci * 3;
            if (nt_pos + 2 >= oriented.length()) break;

            char c0 = fast_upper(oriented[nt_pos]);
            char c1 = fast_upper(oriented[nt_pos + 1]);
            char c2 = fast_upper(oriented[nt_pos + 2]);

            if (!CodonTable::is_stop_codon(c0, c1, c2)) continue;

            bool ct_convertible = (c0 == 'T') &&
                ((c1 == 'A' && c2 == 'A') ||
                 (c1 == 'A' && c2 == 'G') ||
                 (c1 == 'G' && c2 == 'A'));
            bool ga_convertible = (c0 == 'T' && c1 == 'G' && c2 == 'A');

            if (!ct_convertible && !ga_convertible) continue;

            float p_j = 0.0f;
            if (ct_convertible) {
                p_j = compute_position_damage_prob(nt_pos, seq_len, profile, true);
            }
            if (ga_convertible) {
                float p_ga = compute_position_damage_prob(nt_pos + 2, seq_len, profile, false);
                p_j = std::max(p_j, p_ga);
            }

            if (p_j < 0.001f) continue;

            size_t dist_5 = nt_pos;
            size_t dist_3 = (seq_len > nt_pos) ? seq_len - 1 - nt_pos : 0;
            size_t min_dist = std::min(dist_5, dist_3);
            float w_j = std::exp(-0.3f * static_cast<float>(min_dist / 3));

            logbf_stop += w_j * std::log(1.0f + RHO_STOP * p_j);
            ++n_conv_stops;
        }

        // Pre-stop counter-evidence: surviving CAA/CAG/CGA/TGG
        float logbf_prestop = 0.0f;
        uint8_t n_prestops = 0;

        auto scan_prestops = [&](size_t ci_start, size_t ci_end, bool check_5p, bool check_3p) {
            for (size_t ci = ci_start; ci < ci_end; ++ci) {
                size_t nt_pos = orf.nt_start + ci * 3;
                if (nt_pos + 2 >= oriented.length()) break;

                char c0 = fast_upper(oriented[nt_pos]);
                char c1 = fast_upper(oriented[nt_pos + 1]);
                char c2 = fast_upper(oriented[nt_pos + 2]);

                if (check_5p) {
                    bool is_prestop = (c0 == 'C' && c1 == 'A' && c2 == 'A') ||
                                      (c0 == 'C' && c1 == 'A' && c2 == 'G') ||
                                      (c0 == 'C' && c1 == 'G' && c2 == 'A');
                    if (is_prestop) {
                        float p_j = compute_position_damage_prob(nt_pos, seq_len, profile, true);
                        if (p_j > 0.001f) {
                            logbf_prestop += std::log(1.0f - p_j);
                            ++n_prestops;
                        }
                    }
                }
                if (check_3p) {
                    if (c0 == 'T' && c1 == 'G' && c2 == 'G') {
                        float p_j = compute_position_damage_prob(nt_pos + 2, seq_len, profile, false);
                        if (p_j > 0.001f) {
                            logbf_prestop += std::log(1.0f - p_j);
                            ++n_prestops;
                        }
                    }
                }
            }
        };

        scan_prestops(0, std::min(n_codons, static_cast<size_t>(DAMAGE_WINDOW_AA)), true, false);
        if (n_codons > static_cast<size_t>(DAMAGE_WINDOW_AA)) {
            scan_prestops(n_codons - DAMAGE_WINDOW_AA, n_codons, false, true);
        }

        evidence->score_logit = score_logit;
        evidence->p_read_damaged = p_read_damaged;
        evidence->mu_aa_prior = bayesian_damage_score;
        evidence->logbf_terminal = 0.0f;
        evidence->logbf_stop = logbf_stop;
        evidence->logbf_prestop = logbf_prestop;
        evidence->n_informative = n_informative;
        evidence->n_conv_stops = n_conv_stops;
        evidence->n_prestops = n_prestops;
    }

    return p_read_damaged;
}

FrameSelector::FrameshiftResult FrameSelector::detect_frameshifts(
    const std::string& seq,
    const SampleDamageProfile& sample_profile,
    float frameshift_penalty,
    size_t min_segment_codons,
    float min_score_improvement) {

    FrameshiftResult result;
    result.has_frameshift = false;
    result.viterbi_score = -1e9f;
    result.best_single_frame_score = -1e9f;
    result.best_single_frame = 0;
    result.frameshift_position = 0;

    // Need at least 2*min_segment_codons to detect a frameshift
    const size_t min_codons = 2 * min_segment_codons;
    if (seq.length() < min_codons * 3 + 3) {
        return result;  // Too short for frameshift detection
    }

    // Get damage parameters
    auto dmg = make_codon_damage_params(&sample_profile);

    // Score all 6 hypotheses - this populates thread-local buffers
    auto all_scores = codon::score_all_hypotheses(seq, dmg, codon::ScoringWeights(), true);

    // Access thread-local buffers with per-codon data
    auto& buf = codon::get_thread_buffers();

    const size_t num_codons = buf.protein_lens[0];  // All frames have similar codon count
    if (num_codons < min_codons) {
        return result;
    }

    // =========================================================================
    // Step 1: Find best single-frame score (baseline for comparison)
    // =========================================================================
    for (int h = 0; h < 6; ++h) {
        if (all_scores[h].log_likelihood > result.best_single_frame_score) {
            result.best_single_frame_score = all_scores[h].log_likelihood;
            result.best_single_frame = h;
        }
    }

    // =========================================================================
    // Step 2: Viterbi DP for optimal frame path
    // States: 6 frame/strand combinations (fwd 0,1,2 and rev 0,1,2)
    // Emissions: per-codon log marginals from codon_scores
    // Transitions: 0 for same frame, frameshift_penalty for different frame
    // =========================================================================

    // dp[h][i] = best score ending at codon i in hypothesis h
    // backtrack[h][i] = previous hypothesis for backtracking
    std::vector<std::array<float, 6>> dp(num_codons);
    std::vector<std::array<int, 6>> backtrack(num_codons);

    // Initialize first codon
    for (int h = 0; h < 6; ++h) {
        dp[0][h] = buf.codon_scores[h][0];
        backtrack[0][h] = h;  // Start in own state
    }

    // Forward pass
    for (size_t i = 1; i < num_codons; ++i) {
        for (int h = 0; h < 6; ++h) {
            float emission = buf.codon_scores[h][i];
            float best_prev = -1e9f;
            int best_prev_h = h;

            for (int prev_h = 0; prev_h < 6; ++prev_h) {
                float trans = (prev_h == h) ? 0.0f : frameshift_penalty;
                float score = dp[i-1][prev_h] + trans;
                if (score > best_prev) {
                    best_prev = score;
                    best_prev_h = prev_h;
                }
            }

            dp[i][h] = best_prev + emission;
            backtrack[i][h] = best_prev_h;
        }
    }

    // Find best ending state
    float best_final = -1e9f;
    int best_final_h = 0;
    for (int h = 0; h < 6; ++h) {
        if (dp[num_codons - 1][h] > best_final) {
            best_final = dp[num_codons - 1][h];
            best_final_h = h;
        }
    }
    result.viterbi_score = best_final;

    // =========================================================================
    // Step 3: Backtrack to find the frame path
    // =========================================================================
    std::vector<int> path(num_codons);
    path[num_codons - 1] = best_final_h;
    for (int i = static_cast<int>(num_codons) - 2; i >= 0; --i) {
        path[i] = backtrack[i + 1][path[i + 1]];
    }

    // =========================================================================
    // Step 4: Check if frameshift is worthwhile
    // =========================================================================
    float score_improvement = result.viterbi_score - result.best_single_frame_score;

    // Count frame transitions
    size_t n_transitions = 0;
    for (size_t i = 1; i < num_codons; ++i) {
        if (path[i] != path[i-1]) {
            n_transitions++;
        }
    }

    // If no transitions or insufficient improvement, return single-frame result
    if (n_transitions == 0 || score_improvement < min_score_improvement) {
        // Output single best frame as one region
        int h = result.best_single_frame;
        FrameshiftRegion region;
        region.codon_start = 0;
        region.codon_end = num_codons;
        region.frame = h % 3;
        region.forward = (h < 3);
        region.nt_start = 0;
        region.nt_end = seq.length();
        region.score = result.best_single_frame_score;

        // Extract protein from buffer
        region.protein = std::string(buf.protein_buffer[h], buf.protein_lens[h]);

        result.regions.push_back(std::move(region));
        return result;
    }

    // =========================================================================
    // Step 5: Extract regions from path, enforcing minimum segment length
    // =========================================================================
    std::vector<std::pair<size_t, size_t>> raw_regions;  // (start, end) codon indices
    std::vector<int> region_frames;

    size_t region_start = 0;
    for (size_t i = 1; i <= num_codons; ++i) {
        if (i == num_codons || path[i] != path[i-1]) {
            raw_regions.emplace_back(region_start, i);
            region_frames.push_back(path[i-1]);
            region_start = i;
        }
    }

    // Filter out regions shorter than min_segment_codons
    // Merge adjacent short regions with their neighbors
    std::vector<std::pair<size_t, size_t>> valid_regions;
    std::vector<int> valid_frames;

    for (size_t r = 0; r < raw_regions.size(); ++r) {
        size_t len = raw_regions[r].second - raw_regions[r].first;
        if (len >= min_segment_codons) {
            valid_regions.push_back(raw_regions[r]);
            valid_frames.push_back(region_frames[r]);
        } else if (!valid_regions.empty()) {
            // Extend previous region to include this short one
            valid_regions.back().second = raw_regions[r].second;
        }
        // If first region is too short, it will be merged into next valid one
    }

    // If we end up with only one region after filtering, no frameshift
    if (valid_regions.size() <= 1) {
        int h = result.best_single_frame;
        FrameshiftRegion region;
        region.codon_start = 0;
        region.codon_end = num_codons;
        region.frame = h % 3;
        region.forward = (h < 3);
        region.nt_start = 0;
        region.nt_end = seq.length();
        region.score = result.best_single_frame_score;
        region.protein = std::string(buf.protein_buffer[h], buf.protein_lens[h]);
        result.regions.push_back(std::move(region));
        return result;
    }

    // =========================================================================
    // Step 6: Build output regions with proteins
    // =========================================================================
    result.has_frameshift = true;
    result.frameshift_position = valid_regions[0].second;  // First transition point

    for (size_t r = 0; r < valid_regions.size(); ++r) {
        int h = valid_frames[r];
        size_t start = valid_regions[r].first;
        size_t end = valid_regions[r].second;

        FrameshiftRegion region;
        region.codon_start = start;
        region.codon_end = end;
        region.frame = h % 3;
        region.forward = (h < 3);
        region.nt_start = start * 3;
        region.nt_end = end * 3;

        // Compute region score from per-codon scores
        region.score = 0.0f;
        for (size_t i = start; i < end; ++i) {
            region.score += buf.codon_scores[h][i];
        }

        // Extract protein substring for this region
        if (end <= buf.protein_lens[h]) {
            region.protein = std::string(
                buf.protein_buffer[h] + start,
                end - start);
        }

        result.regions.push_back(std::move(region));
    }

    return result;
}

}  // namespace dart
