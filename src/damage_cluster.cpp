// src/damage_cluster.cpp
// Damage-aware read clustering for per-read damage estimation
//
// Hard-EM-lite 2-round algorithm:
// 1. Stream reads -> compact store + seed index
// 2. Round 1: Score frames using coding + cluster support + terminal damage
// 3. Round 2: Rescore with neighbors restricted to best slots
// 4. Compute per-read damage + X-mask positions
//
// Design from GPT-5.3-codex discussion (2026-02-06)

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "agp/damage_cluster.hpp"
#include "agp/damage_cluster_seeds.hpp"

namespace agp {
namespace damage_cluster {

// -----------------------------
// Config / result structs
// -----------------------------
struct BuildConfig {
    uint8_t strict_k = 7;
    uint8_t tolerant_k = 5;
    uint8_t terminal_window = 12;

    uint8_t max_strict_seeds_per_frame = 8;
    uint8_t max_tolerant_seeds_per_frame = 4;

    uint32_t min_df = 2;
    uint32_t max_df = 0;  // 0 = derived from IDF prior

    // IDF prior (learned from warmup)
    SeedIndex::IDFPrior idf_prior = SeedIndex::IDFPrior::jeffreys();
};

// Learned parameters container (calibrated from warmup, not hardcoded)
struct ClusterParams {
    FrameWeights frame_weights;          // Weight for code/support/terminal
    SoftmaxTemperature softmax_temp;     // Temperature for frame posteriors
    float kappa_mis_over_fp = 1.0f;      // Cost ratio for X-mask decision
    float d_max = 0.0f;                  // Sample damage rate
    float lambda = 0.0f;                 // Damage decay constant

    static ClusterParams defaults(float d_max_, float lambda_) {
        ClusterParams p;
        p.frame_weights = FrameWeights::defaults_for_damage(d_max_);
        p.softmax_temp = SoftmaxTemperature::defaults_for_damage(d_max_);
        p.kappa_mis_over_fp = 1.0f;
        p.d_max = d_max_;
        p.lambda = lambda_;
        return p;
    }
};

struct NormalizationStats {
    float mu_code = 0.0f, sigma_code = 1.0f;
    float mu_sup = 0.0f, sigma_sup = 1.0f;
    float mu_term = 0.0f, sigma_term = 1.0f;
};

struct BuildArtifacts {
    ReadCompactStore reads;
    SeedIndex index;
};

template <size_t TOPK = 64>
struct FrameSlotEvidence {
    float code_score = 0.0f;
    float cluster_support = 0.0f;
    float terminal_llr = 0.0f;
    float combined_score = 0.0f;
    float posterior = 0.0f;

    uint8_t n_neighbors = 0;
    std::array<uint32_t, TOPK> neighbor_nodes{};
    std::array<float, TOPK> neighbor_scores{};
    std::array<int16_t, TOPK> neighbor_delta{}; // anchored offset: neigh_pos - query_pos
};

template <size_t TOPK = 64>
struct FrameScoreResult {
    uint8_t n_slots = 0;
    std::array<FrameSlotEvidence<TOPK>, 3> slot{};
};

struct PerReadDamageResult {
    float damage_score = 0.0f;               // mean posterior across informative terminal mismatches
    std::vector<uint16_t> x_mask_positions;  // AA positions in search protein to X-mask
};

// -----------------------------
// Internal helpers
// -----------------------------
inline uint32_t make_node_id(uint32_t read_id, uint8_t slot) {
    return (read_id << 2) | static_cast<uint32_t>(slot & 0x3u);
}
inline uint32_t node_read_id(uint32_t node_id) { return node_id >> 2; }
inline uint8_t node_slot(uint32_t node_id) { return static_cast<uint8_t>(node_id & 0x3u); }

inline float safe_norm_at(const std::vector<float>& v, uint32_t i) {
    return (i < v.size()) ? v[i] : 1.0f;
}

inline void keep_lowest_hashes(std::vector<SeedHit>& seeds, uint8_t keep_n) {
    if (seeds.empty() || seeds.size() <= keep_n) return;
    std::nth_element(
        seeds.begin(),
        seeds.begin() + keep_n,
        seeds.end(),
        [](const SeedHit& a, const SeedHit& b) {
            if (a.hash != b.hash) return a.hash < b.hash;
            return a.aa_pos < b.aa_pos;
        }
    );
    seeds.resize(keep_n);
}

inline void dedup_seeds(std::vector<SeedHit>& seeds) {
    std::sort(seeds.begin(), seeds.end(), [](const SeedHit& a, const SeedHit& b) {
        if (a.hash != b.hash) return a.hash < b.hash;
        if (a.aa_pos != b.aa_pos) return a.aa_pos < b.aa_pos;
        if (a.terminal != b.terminal) return a.terminal < b.terminal;
        return a.tolerant < b.tolerant;
    });
    seeds.erase(std::unique(seeds.begin(), seeds.end(), [](const SeedHit& a, const SeedHit& b) {
        return a.hash == b.hash &&
               a.aa_pos == b.aa_pos &&
               a.terminal == b.terminal &&
               a.tolerant == b.tolerant;
    }), seeds.end());
}

inline void collect_query_seeds(
    const std::string& protein,
    const BuildConfig& cfg,
    std::vector<SeedHit>& out
) {
    std::vector<SeedHit> strict_seeds;
    std::vector<SeedHit> tol_seeds;
    emit_strict_seeds(protein, strict_seeds, cfg.strict_k);
    emit_tolerant_seeds(protein, tol_seeds, cfg.tolerant_k, cfg.terminal_window);

    keep_lowest_hashes(strict_seeds, cfg.max_strict_seeds_per_frame);
    keep_lowest_hashes(tol_seeds, cfg.max_tolerant_seeds_per_frame);

    out.clear();
    out.reserve(strict_seeds.size() + tol_seeds.size());
    out.insert(out.end(), strict_seeds.begin(), strict_seeds.end());
    out.insert(out.end(), tol_seeds.begin(), tol_seeds.end());

    dedup_seeds(out);
}

template <size_t TOPK>
inline int find_neighbor_idx(
    const std::array<uint32_t, TOPK>& nodes,
    uint8_t n,
    uint32_t node
) {
    for (uint8_t i = 0; i < n; ++i) {
        if (nodes[i] == node) return static_cast<int>(i);
    }
    return -1;
}

// Re-scan postings for selected neighbors to estimate a single offset per neighbor.
template <size_t TOPK>
inline void estimate_neighbor_offsets(
    const std::vector<SeedHit>& query_seeds,
    const SeedIndex& index,
    FrameSlotEvidence<TOPK>& evidence
) {
    std::array<float, TOPK> delta_sum{};
    std::array<float, TOPK> delta_w{};
    delta_sum.fill(0.0f);
    delta_w.fill(0.0f);

    for (const auto& qs : query_seeds) {
        const auto range = index.lookup(qs.hash);
        if (range.empty()) continue;

        for (const auto* p = range.begin(); p != range.end(); ++p) {
            const uint32_t n = PackedPosting::node_id(*p);
            const int j = find_neighbor_idx(evidence.neighbor_nodes, evidence.n_neighbors, n);
            if (j < 0) continue;

            const int d = static_cast<int>(PackedPosting::aa_pos(*p)) - static_cast<int>(qs.aa_pos);
            delta_sum[static_cast<size_t>(j)] += range.idf * static_cast<float>(d);
            delta_w[static_cast<size_t>(j)] += range.idf;
        }
    }

    for (uint8_t i = 0; i < evidence.n_neighbors; ++i) {
        if (delta_w[i] > 0.0f) {
            const float d = delta_sum[i] / delta_w[i];
            const int di = static_cast<int>(std::lround(d));
            evidence.neighbor_delta[i] = static_cast<int16_t>(
                std::max(-1024, std::min(1024, di))
            );
        } else {
            evidence.neighbor_delta[i] = 0;
        }
    }
}

template <size_t TOPK, typename TranslateFn>
inline void build_terminal_consensus_from_evidence(
    const std::string& query_protein,
    const ReadCompactStore& reads,
    const FrameSlotEvidence<TOPK>& evidence,
    TranslateFn&& translate_frame,
    uint8_t terminal_window,
    TerminalConsensus<>& consensus
) {
    consensus.reset(query_protein.size(), terminal_window);
    const size_t qlen = query_protein.size();
    if (qlen == 0) return;

    const uint8_t w = static_cast<uint8_t>(std::min<size_t>(terminal_window, qlen));

    for (uint8_t ni = 0; ni < evidence.n_neighbors; ++ni) {
        const uint32_t node = evidence.neighbor_nodes[ni];
        const uint32_t nrid = node_read_id(node);
        const uint8_t nslot = node_slot(node);
        if (nrid >= reads.size()) continue;

        const auto& nr = reads.read(nrid);
        if (nslot >= nr.n_candidates) continue;

        const std::string nt = reads.decode_read(nrid);
        const uint8_t nframe = nr.cand[nslot].frame;
        const std::string np = translate_frame(nt, nframe);
        if (np.empty()) continue;

        const int delta = static_cast<int>(evidence.neighbor_delta[ni]);
        const float wgt = evidence.neighbor_scores[ni];

        // 5'
        for (uint8_t i = 0; i < w; ++i) {
            const int qi = static_cast<int>(i);
            const int nj = qi + delta;
            if (nj >= 0 && nj < static_cast<int>(np.size())) {
                consensus.add(true, i, np[static_cast<size_t>(nj)], wgt);
            }
        }

        // 3'
        for (uint8_t i = 0; i < w; ++i) {
            const int qi = static_cast<int>(qlen - 1 - i);
            const int nj = qi + delta;
            if (nj >= 0 && nj < static_cast<int>(np.size())) {
                consensus.add(false, i, np[static_cast<size_t>(nj)], wgt);
            }
        }
    }
}

// -----------------------------
// 1) build_seed_index_from_reads()
// -----------------------------
// select_candidates signature:
//   void(const std::string& nt, std::array<uint8_t,3>& frames,
//        std::array<float,3>& code_scores, uint8_t& n_candidates)
//
// translate_frame signature:
//   std::string(const std::string& nt, uint8_t frame)
template <typename CandidateSelectorFn, typename TranslateFn>
BuildArtifacts build_seed_index_from_reads(
    const std::vector<std::string>& nt_reads,
    CandidateSelectorFn&& select_candidates,
    TranslateFn&& translate_frame,
    const BuildConfig& cfg = {}
) {
    BuildArtifacts out;
    out.reads.reserve_reads(nt_reads.size());

    SeedIndex::Builder builder(static_cast<uint32_t>(nt_reads.size() * 4u));
    // Rough reserve: 12 seeds/frame * ~2.3 frames/read
    builder.reserve(nt_reads.size() * 28u);

    for (const auto& nt : nt_reads) {
        std::array<uint8_t, 3> frames{0, 1, 2};
        std::array<float, 3> scores{0.0f, 0.0f, 0.0f};
        uint8_t n_cand = 0;
        select_candidates(nt, frames, scores, n_cand);
        n_cand = static_cast<uint8_t>(std::max<uint8_t>(1, std::min<uint8_t>(n_cand, 3)));

        const uint32_t rid = out.reads.add_read(
            nt.c_str(),
            static_cast<uint16_t>(nt.size()),
            frames,
            scores,
            n_cand
        );

        for (uint8_t slot = 0; slot < n_cand; ++slot) {
            const std::string protein = translate_frame(nt, frames[slot]);
            if (protein.empty()) continue;

            std::vector<SeedHit> seeds;
            collect_query_seeds(protein, cfg, seeds);

            for (const auto& sh : seeds) {
                builder.append(
                    sh.hash,
                    PackedPosting::pack(
                        rid,
                        slot,
                        sh.aa_pos,
                        sh.terminal,
                        sh.tolerant
                    )
                );
            }
        }
    }

    builder.build(out.index, cfg.min_df, cfg.max_df, cfg.idf_prior);
    return out;
}

// -----------------------------
// 2) score_read_frames()
// -----------------------------
template <size_t ACC_CAP = 8192, size_t TOPK = 64, typename TranslateFn>
FrameScoreResult<TOPK> score_read_frames(
    uint32_t read_id,
    ReadCompactStore& reads,
    const SeedIndex& index,
    const std::vector<float>& node_norm_sq,
    const ClusterParams& params,
    const NormalizationStats& stats,
    TranslateFn&& translate_frame,
    const BuildConfig& cfg = {}
) {
    FrameScoreResult<TOPK> res;
    if (read_id >= reads.size()) return res;

    auto& r = reads.read_mut(read_id);
    res.n_slots = r.n_candidates;

    const std::string nt = reads.decode_read(read_id);
    if (nt.empty() || r.n_candidates == 0) return res;

    NeighborAccumulator<ACC_CAP, TOPK> acc;
    std::array<typename NeighborAccumulator<ACC_CAP, TOPK>::Neighbor, TOPK> topn{};

    std::array<float, 3> combined{0.0f, 0.0f, 0.0f};

    for (uint8_t slot = 0; slot < r.n_candidates; ++slot) {
        auto& ev = res.slot[slot];
        ev.code_score = r.cand[slot].code_score;

        const std::string qp = translate_frame(nt, r.cand[slot].frame);
        if (qp.empty()) {
            ev.combined_score = ev.code_score;
            combined[slot] = ev.combined_score;
            continue;
        }

        std::vector<SeedHit> qseeds;
        collect_query_seeds(qp, cfg, qseeds);
        if (qseeds.empty()) {
            ev.combined_score = ev.code_score;
            combined[slot] = ev.combined_score;
            continue;
        }

        const uint32_t qnode = make_node_id(read_id, slot);
        acc.reset(qnode, safe_norm_at(node_norm_sq, qnode));

        for (const auto& qs : qseeds) {
            const auto range = index.lookup(qs.hash);
            if (range.empty()) continue;

            for (const auto* p = range.begin(); p != range.end(); ++p) {
                const uint32_t nid = PackedPosting::node_id(*p);
                acc.add(nid, range.idf);
            }
        }

        const size_t n_top = acc.topk(
            [&](uint32_t nid) { return safe_norm_at(node_norm_sq, nid); },
            topn,
            0.0f
        );

        ev.n_neighbors = static_cast<uint8_t>(std::min<size_t>(n_top, TOPK));
        for (uint8_t i = 0; i < ev.n_neighbors; ++i) {
            ev.neighbor_nodes[i] = topn[i].node_id;
            ev.neighbor_scores[i] = topn[i].score;
        }

        // Cluster support = mean of top-L neighbor scores
        constexpr size_t L = 16;
        const size_t use_n = std::min<size_t>(ev.n_neighbors, L);
        float sup = 0.0f;
        for (size_t i = 0; i < use_n; ++i) sup += ev.neighbor_scores[i];
        ev.cluster_support = (use_n > 0) ? (sup / static_cast<float>(use_n)) : 0.0f;

        // Estimate neighbor offsets from shared seed anchors
        estimate_neighbor_offsets(qseeds, index, ev);

        // Terminal consensus + terminal LLR
        TerminalConsensus<> cons;
        build_terminal_consensus_from_evidence(
            qp, reads, ev, translate_frame, cfg.terminal_window, cons
        );
        ev.terminal_llr = compute_terminal_damage_llr(qp, cons, params.d_max, params.lambda, 0.10f);

        ev.combined_score = score_frame(
            ev.code_score,
            ev.cluster_support,
            ev.terminal_llr,
            stats.mu_code, stats.sigma_code,
            stats.mu_sup, stats.sigma_sup,
            stats.mu_term, stats.sigma_term,
            params.frame_weights
        );
        combined[slot] = ev.combined_score;
    }

    // Softmax posteriors using learned temperature
    const auto post = softmax_frame_posteriors(combined, r.n_candidates, params.softmax_temp.tau);
    for (uint8_t slot = 0; slot < r.n_candidates; ++slot) {
        res.slot[slot].posterior = post[slot];
        r.cand[slot].combined_score = res.slot[slot].combined_score;
        r.cand[slot].posterior = post[slot];
    }

    return res;
}

// -----------------------------
// 3) compute_per_read_damage()
// -----------------------------
// best_slot is candidate slot [0..n_candidates-1] for this read.
template <size_t TOPK = 64, typename TranslateFn>
PerReadDamageResult compute_per_read_damage(
    uint32_t read_id,
    uint8_t best_slot,
    const ReadCompactStore& reads,
    const SeedIndex& index, // kept for signature symmetry / future use
    const FrameSlotEvidence<TOPK>& neighbors,
    const ClusterParams& params,
    float frame_posterior,
    TranslateFn&& translate_frame,
    uint8_t terminal_window = 12
) {
    (void)index; // not needed in current implementation

    PerReadDamageResult out;
    if (read_id >= reads.size()) return out;

    const auto& r = reads.read(read_id);
    if (best_slot >= r.n_candidates) return out;

    const std::string nt = reads.decode_read(read_id);
    const std::string qp = translate_frame(nt, r.cand[best_slot].frame);
    if (qp.empty()) return out;

    TerminalConsensus<> cons;
    build_terminal_consensus_from_evidence(
        qp, reads, neighbors, translate_frame, terminal_window, cons
    );

    DamagePosterior dp;
    const uint8_t w = static_cast<uint8_t>(std::min<size_t>(terminal_window, qp.size()));

    float sum_post = 0.0f;
    uint32_t n_post = 0;

    auto eval_site = [&](bool five_prime, uint8_t i, size_t qpos) {
        const auto s = cons.site(five_prime, i);
        if (s.depth < 1.0f || s.consensus == 'X') return;

        const char obs = qp[qpos];
        if (obs == s.consensus) return; // not informative for damage mismatch posterior

        DamagePosterior::Input in;
        in.consensus_aa = s.consensus;
        in.observed_aa = obs;
        in.five_prime = five_prime;
        in.three_prime = !five_prime;
        in.dist_to_end = static_cast<float>(i);
        in.d_max = params.d_max;
        in.lambda = params.lambda;
        in.consensus_freq = s.freq;
        in.consensus_depth = s.depth;
        in.frame_posterior = frame_posterior;

        const float p = dp.posterior(in);
        sum_post += p;
        ++n_post;

        // X-mask decision using principled risk minimization
        if (DamagePosterior::should_xmask(p, s.depth, s.freq, frame_posterior,
                                          params.d_max, params.lambda,
                                          static_cast<float>(i), params.kappa_mis_over_fp)) {
            out.x_mask_positions.push_back(static_cast<uint16_t>(qpos));
        }
    };

    // 5'
    for (uint8_t i = 0; i < w; ++i) {
        eval_site(true, i, static_cast<size_t>(i));
    }
    // 3' (avoid overlap double-count)
    for (uint8_t i = 0; i < w; ++i) {
        const size_t qpos = qp.size() - 1 - static_cast<size_t>(i);
        if (qpos < static_cast<size_t>(w)) continue;
        eval_site(false, i, qpos);
    }

    out.damage_score = (n_post > 0) ? (sum_post / static_cast<float>(n_post)) : 0.0f;
    std::sort(out.x_mask_positions.begin(), out.x_mask_positions.end());
    out.x_mask_positions.erase(
        std::unique(out.x_mask_positions.begin(), out.x_mask_positions.end()),
        out.x_mask_positions.end()
    );
    return out;
}

// -----------------------------
// 4) Driver: run_damage_clustering_rounds()
// -----------------------------
struct DamageClusteringResult {
    std::vector<uint8_t> best_frames;                    // biological frame id (0..5) per read
    std::vector<float> frame_posteriors;                 // posterior of best slot
    std::vector<float> damage_scores;                    // per-read damage score
    std::vector<std::vector<uint16_t>> x_mask_positions; // AA positions to X-mask per read
};

namespace driver_detail {

struct RunningMoment {
    double sum = 0.0;
    double sumsq = 0.0;
    uint64_t n = 0;

    inline void add(float x) {
        sum += static_cast<double>(x);
        sumsq += static_cast<double>(x) * static_cast<double>(x);
        ++n;
    }

    inline void mean_sigma(float& mu, float& sigma) const {
        if (n == 0) {
            mu = 0.0f;
            sigma = 1.0f;
            return;
        }
        const double dn = static_cast<double>(n);
        const double m = sum / dn;
        const double var = std::max(0.0, (sumsq / dn) - (m * m));
        mu = static_cast<float>(m);
        sigma = static_cast<float>(std::sqrt(std::max(var, 1e-8)));
        if (sigma < 1e-4f) sigma = 1.0f;
    }
};

inline NormalizationStats finalize_stats(
    const RunningMoment& code_m,
    const RunningMoment& sup_m,
    const RunningMoment& term_m
) {
    NormalizationStats s{};
    code_m.mean_sigma(s.mu_code, s.sigma_code);
    sup_m.mean_sigma(s.mu_sup, s.sigma_sup);
    term_m.mean_sigma(s.mu_term, s.sigma_term);
    return s;
}

template <size_t TOPK>
inline uint8_t argmax_posterior_slot(const FrameScoreResult<TOPK>& fs) {
    if (fs.n_slots == 0) return 0;
    uint8_t best = 0;
    float bestv = fs.slot[0].posterior;
    for (uint8_t i = 1; i < fs.n_slots; ++i) {
        if (fs.slot[i].posterior > bestv) {
            bestv = fs.slot[i].posterior;
            best = i;
        }
    }
    return best;
}

template <size_t TOPK>
inline FrameSlotEvidence<TOPK> filter_neighbors_to_best_slots(
    const FrameSlotEvidence<TOPK>& in,
    const std::vector<uint8_t>& best_slot_by_read
) {
    FrameSlotEvidence<TOPK> out{};
    out.code_score = in.code_score;
    out.cluster_support = 0.0f;
    out.terminal_llr = 0.0f;
    out.combined_score = 0.0f;
    out.posterior = 0.0f;
    out.n_neighbors = 0;

    for (uint8_t i = 0; i < in.n_neighbors; ++i) {
        const uint32_t node = in.neighbor_nodes[i];
        const uint32_t nrid = node_read_id(node);
        const uint8_t nslot = node_slot(node);
        if (nrid >= best_slot_by_read.size()) continue;
        if (nslot != best_slot_by_read[nrid]) continue;

        const uint8_t j = out.n_neighbors++;
        out.neighbor_nodes[j] = node;
        out.neighbor_scores[j] = in.neighbor_scores[i];
        out.neighbor_delta[j] = in.neighbor_delta[i];
        if (out.n_neighbors >= TOPK) break;
    }

    // mean top-L cluster support from retained neighbors
    constexpr uint8_t L = 16;
    const uint8_t use_n = std::min<uint8_t>(out.n_neighbors, L);
    float acc = 0.0f;
    for (uint8_t i = 0; i < use_n; ++i) acc += out.neighbor_scores[i];
    out.cluster_support = (use_n > 0) ? (acc / static_cast<float>(use_n)) : 0.0f;

    return out;
}

} // namespace driver_detail

template <typename CandidateSelectorFn, typename TranslateFn>
DamageClusteringResult run_damage_clustering_rounds(
    const std::vector<std::string>& nt_reads,
    CandidateSelectorFn&& select_candidates,
    TranslateFn&& translate_frame,
    float d_max,
    float lambda,
    const BuildConfig& cfg = {},
    const ClusterParams* params_override = nullptr  // If null, use defaults from d_max/lambda
) {
    using driver_detail::RunningMoment;
    using driver_detail::finalize_stats;

    DamageClusteringResult out;
    const size_t n_reads = nt_reads.size();

    out.best_frames.assign(n_reads, 0);
    out.frame_posteriors.assign(n_reads, 0.0f);
    out.damage_scores.assign(n_reads, 0.0f);
    out.x_mask_positions.assign(n_reads, {});

    if (n_reads == 0) return out;

    // Use provided params or derive defaults from d_max/lambda
    const ClusterParams params = params_override
        ? *params_override
        : ClusterParams::defaults(d_max, lambda);

    // 1) Build compact state + seed index
    BuildArtifacts artifacts = build_seed_index_from_reads(
        nt_reads,
        std::forward<CandidateSelectorFn>(select_candidates),
        std::forward<TranslateFn>(translate_frame),
        cfg
    );

    const std::vector<float> node_norm_sq = artifacts.index.compute_node_idf_norm_sq();

    // 2) Compute normalization stats from probe pass (un-normalized scoring components)
    const NormalizationStats identity_stats{
        0.0f, 1.0f,
        0.0f, 1.0f,
        0.0f, 1.0f
    };

    std::vector<FrameScoreResult<64>> probe_round1(n_reads);
    RunningMoment m_code_r1, m_sup_r1, m_term_r1;

    for (uint32_t rid = 0; rid < n_reads; ++rid) {
        auto fs = score_read_frames<8192, 64>(
            rid,
            artifacts.reads,
            artifacts.index,
            node_norm_sq,
            params,
            identity_stats,
            std::forward<TranslateFn>(translate_frame),
            cfg
        );
        probe_round1[rid] = fs;
        for (uint8_t s = 0; s < fs.n_slots; ++s) {
            m_code_r1.add(fs.slot[s].code_score);
            m_sup_r1.add(fs.slot[s].cluster_support);
            m_term_r1.add(fs.slot[s].terminal_llr);
        }
    }

    const NormalizationStats stats_r1 = finalize_stats(m_code_r1, m_sup_r1, m_term_r1);

    // 3) Round 1 final scoring (all neighbors allowed)
    std::vector<FrameScoreResult<64>> round1(n_reads);
    std::vector<uint8_t> best_slot_r1(n_reads, 0);

    for (uint32_t rid = 0; rid < n_reads; ++rid) {
        round1[rid] = score_read_frames<8192, 64>(
            rid,
            artifacts.reads,
            artifacts.index,
            node_norm_sq,
            params,
            stats_r1,
            std::forward<TranslateFn>(translate_frame),
            cfg
        );
        best_slot_r1[rid] = driver_detail::argmax_posterior_slot(round1[rid]);
    }

    // 4) Round 2: restrict neighbors to each neighbor-read's best slot from round 1
    std::vector<FrameScoreResult<64>> round2(n_reads);
    RunningMoment m_code_r2, m_sup_r2, m_term_r2;

    for (uint32_t rid = 0; rid < n_reads; ++rid) {
        const auto& rstate = artifacts.reads.read(rid);
        FrameScoreResult<64> fs2{};
        fs2.n_slots = rstate.n_candidates;

        const std::string nt = artifacts.reads.decode_read(rid);

        for (uint8_t slot = 0; slot < fs2.n_slots; ++slot) {
            // Start from round-1 evidence and prune neighbors by best_slot_r1
            FrameSlotEvidence<64> ev = driver_detail::filter_neighbors_to_best_slots(
                round1[rid].slot[slot], best_slot_r1
            );

            const std::string qp = translate_frame(nt, rstate.cand[slot].frame);
            if (!qp.empty() && ev.n_neighbors > 0) {
                TerminalConsensus<> cons;
                build_terminal_consensus_from_evidence(
                    qp, artifacts.reads, ev,
                    std::forward<TranslateFn>(translate_frame),
                    cfg.terminal_window, cons
                );
                ev.terminal_llr = compute_terminal_damage_llr(qp, cons, params.d_max, params.lambda, 0.10f);
            } else {
                ev.terminal_llr = 0.0f;
            }

            fs2.slot[slot] = ev;
            m_code_r2.add(ev.code_score);
            m_sup_r2.add(ev.cluster_support);
            m_term_r2.add(ev.terminal_llr);
        }

        round2[rid] = fs2;
    }

    const NormalizationStats stats_r2 = finalize_stats(m_code_r2, m_sup_r2, m_term_r2);

    // Apply normalized frame score + softmax for round 2
    std::vector<uint8_t> best_slot_r2(n_reads, 0);

    for (uint32_t rid = 0; rid < n_reads; ++rid) {
        auto& fs = round2[rid];
        std::array<float, 3> scores{0.0f, 0.0f, 0.0f};

        for (uint8_t slot = 0; slot < fs.n_slots; ++slot) {
            auto& ev = fs.slot[slot];
            ev.combined_score = score_frame(
                ev.code_score,
                ev.cluster_support,
                ev.terminal_llr,
                stats_r2.mu_code, stats_r2.sigma_code,
                stats_r2.mu_sup, stats_r2.sigma_sup,
                stats_r2.mu_term, stats_r2.sigma_term,
                params.frame_weights
            );
            scores[slot] = ev.combined_score;
        }

        const auto post = softmax_frame_posteriors(scores, fs.n_slots, params.softmax_temp.tau);
        for (uint8_t slot = 0; slot < fs.n_slots; ++slot) {
            fs.slot[slot].posterior = post[slot];
        }

        best_slot_r2[rid] = driver_detail::argmax_posterior_slot(fs);
    }

    // 5) Final per-read damage + X-mask using best round-2 slot
    for (uint32_t rid = 0; rid < n_reads; ++rid) {
        const auto& rstate = artifacts.reads.read(rid);
        if (rstate.n_candidates == 0) continue;

        const uint8_t bslot = best_slot_r2[rid];
        const float bpost = round2[rid].slot[bslot].posterior;
        const uint8_t bframe = rstate.cand[bslot].frame;

        out.best_frames[rid] = bframe;
        out.frame_posteriors[rid] = bpost;

        const PerReadDamageResult dres = compute_per_read_damage<64>(
            rid,
            bslot,
            artifacts.reads,
            artifacts.index,
            round2[rid].slot[bslot],  // already restricted neighbors
            params,
            bpost,
            std::forward<TranslateFn>(translate_frame),
            cfg.terminal_window
        );

        out.damage_scores[rid] = dres.damage_score;
        out.x_mask_positions[rid] = dres.x_mask_positions;
    }

    return out;
}

} // namespace damage_cluster
} // namespace agp
