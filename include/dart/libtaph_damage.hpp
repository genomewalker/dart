#pragma once
/**
 * @file libtaph_damage.hpp
 * @brief Per-read damage authentication via libtaph's reference-free terminal-site LLR.
 *
 * Replaces the frame_selector AA-window path (infer_per_read_aa_damage), which
 * floored short reads and folded reference identity into the score. libtaph's
 * damage_evidence_llr is PRIOR-FREE, reference-free, and has NO length floor:
 * it scores terminal C->T (5') / G->A (3' ds) / C->T (3' ss) excess against the
 * library's fitted amplitude/decay and interior baseline, so genuine composition
 * is nulled and short fragments (which still have termini) get a finite score.
 */

#include <cmath>
#include <string_view>
#include <vector>

#include "taph/read_ancient_llr.hpp"
#include "taph/sample_damage_profile.hpp"

namespace dart {

/// Owning storage for a read's eligible terminal sites plus a non-owning view.
struct ReadTerminalObs {
    std::vector<taph::TerminalMismatch> five;
    std::vector<taph::TerminalMismatch> three;
    int length = 0;

    taph::ReadDamageObs view() const {
        taph::ReadDamageObs o;
        o.length = length;
        o.five = five.empty() ? nullptr : five.data();
        o.n_five = static_cast<std::uint32_t>(five.size());
        o.three = three.empty() ? nullptr : three.data();
        o.n_three = static_cast<std::uint32_t>(three.size());
        return o;
    }
};

/// Reference-free terminal observations from the read itself. A 5' C->T site must
/// currently read C (not deaminated) or T (deaminated); ineligible if A/G. The 3'
/// channel is G->A for ds libraries, C->T for ss. `pos` is 0-indexed from the terminus.
inline ReadTerminalObs terminal_obs_reffree(std::string_view seq, bool ss_library,
                                            int window = taph::SampleDamageProfile::P_PI) {
    ReadTerminalObs obs;
    obs.length = static_cast<int>(seq.size());
    const int L = obs.length;
    const int W = window < L ? window : L;
    for (int p = 0; p < W; ++p) {
        const char b = seq[p];  // 5' C->T
        if (b == 'C' || b == 'c') obs.five.push_back({p, false});
        else if (b == 'T' || b == 't') obs.five.push_back({p, true});
    }
    // 3' starts at p=1: the profile's 3' amplitude/baseline fit excludes terminal
    // position 0, so scoring the highest-leverage terminal base against an amplitude
    // never fit for it biases the per-read LLR. Align the scoring domain with the fit.
    for (int p = 1; p < W; ++p) {
        const char b = seq[L - 1 - p];  // 3'
        if (ss_library) {  // ss 3' C->T
            if (b == 'C' || b == 'c') obs.three.push_back({p, false});
            else if (b == 'T' || b == 't') obs.three.push_back({p, true});
        } else {           // ds 3' G->A
            if (b == 'G' || b == 'g') obs.three.push_back({p, false});
            else if (b == 'A' || b == 'a') obs.three.push_back({p, true});
        }
    }
    return obs;
}

/// Per-read P(damaged) via libtaph. Prior-aware posterior when the library validated
/// as DETECTED; otherwise the prior-free LLR through a sigmoid (a finite ranking score
/// for every read). NEVER the AA-window path — short reads keep a real terminal score.
inline float per_read_damage_libtaph(std::string_view seq,
                                     const taph::SampleDamageProfile& profile) {
    const bool ss = profile.library_type ==
                    taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
    const ReadTerminalObs obs = terminal_obs_reffree(seq, ss);
    const taph::ReadDamageObs v = obs.view();
    if (auto post = taph::read_ancient_posterior(v, profile))
        return static_cast<float>(*post);
    const double llr = taph::damage_evidence_llr(v, profile);
    return static_cast<float>(1.0 / (1.0 + std::exp(-llr)));
}

/// Raw PRIOR-FREE per-read damage LLR: log P(obs|damaged) − log P(obs|non-damaged) via
/// libtaph. Same terminal-site construction as per_read_damage_libtaph, but returns
/// taph::damage_evidence_llr directly — NOT through a sigmoid and NEVER the prior-aware
/// posterior — so per-read evidence can be POOLED at the protein level without folding a
/// per-read prior into each term: 0 for an evidence-free read, >0 for terminal-deamination
/// excess, <0 for a non-deaminated (anti-damage) terminus. Prefer this over logit() of the
/// 8-bit .agd p_damaged_q, which loses resolution near 0/1. (Requires seq + profile, i.e. the
/// predict-time path; damage-annotate reads the .agd and has neither, so it must carry a
/// per-read LLR forward if full resolution is needed there.)
inline float per_read_damage_llr_libtaph(std::string_view seq,
                                         const taph::SampleDamageProfile& profile) {
    const bool ss = profile.library_type ==
                    taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
    const ReadTerminalObs obs = terminal_obs_reffree(seq, ss);
    return static_cast<float>(taph::damage_evidence_llr(obs.view(), profile));
}

}  // namespace dart
