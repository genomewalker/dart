// tests/test_formal_model.cpp
//
// Phase 4 CI tests for the Formal Model and Proof Roadmap theorems:
//
//   T1: ELBO/MAP-objective monotonicity
//       (a) plain EM using inline e_step / m_step / em_objective_from_ll
//       (b) SQUAREM path via squarem_em with per-iteration progress callback
//
//   T2: Stationary-point convergence (parameter recovery as proxy)
//       (a) unique read-to-reference mappings  → EM recovers ground-truth weights
//       (b) symmetric ambiguous mappings       → EM converges to [0.5, 0.5]
//       (c) asymmetric ambiguous mappings      → dominant reference wins
//
//   T3: FDR control contract
//       For every threshold returned by find_threshold_for_fdr, the estimated
//       FDR at that threshold is ≤ the requested target.
//
//   T4 (direction test / output validity): damage-aware EM pi
//       All-ancient data pulls pi > 0.5; all-modern data pulls pi < 0.5.
//       Full SBC calibration is deferred to Phase 4 integration tests.

#include "dart/em_reassign.hpp"
#include "dart/fdr_estimator.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

void expect(bool ok, const std::string& msg, int& failed) {
    if (!ok) {
        std::cerr << "  FAIL: " << msg << "\n";
        ++failed;
    }
}

// Build AlignmentData from a list of reads.
// Each read is a list of (ref_idx, bit_score, damage_score) tuples.
dart::AlignmentData make_data(
    uint32_t num_refs,
    const std::vector<std::vector<std::tuple<uint32_t, float, float>>>& reads)
{
    dart::AlignmentData d;
    d.num_refs  = num_refs;
    d.num_reads = static_cast<uint32_t>(reads.size());
    d.read_offsets.resize(reads.size() + 1);

    uint32_t off = 0;
    for (uint32_t r = 0; r < static_cast<uint32_t>(reads.size()); ++r) {
        d.read_offsets[r] = off;
        for (auto& [ref, score, dam] : reads[r]) {
            dart::Alignment a{};
            a.read_idx     = r;
            a.ref_idx      = ref;
            a.bit_score    = score;
            a.damage_score = dam;
            a.damage_ll_a  = 0.0f;
            a.damage_ll_m  = 0.0f;
            a.fident       = 0.9f;
            a.tstart       = 0;
            a.aln_len      = 60;
            d.alignments.push_back(a);
        }
        off += static_cast<uint32_t>(reads[r].size());
    }
    d.read_offsets[reads.size()] = off;
    d.ref_lengths.assign(num_refs, 200u);
    d.ref_log_len_half.resize(num_refs, 0.5f * std::log(200.0f));
    return d;
}

// ---------------------------------------------------------------------------
// T1a: Plain-EM objective monotonicity (header-only inline path)
// ---------------------------------------------------------------------------
int test_elbo_monotone_plain_em() {
    std::cout << "[T1a] plain-EM objective monotonicity\n";
    int failed = 0;

    // 4 refs, 40 reads, each mapping to 2 refs with different bit scores
    std::vector<std::vector<std::tuple<uint32_t, float, float>>> reads;
    for (int i = 0; i < 20; ++i) reads.push_back({{0, 1.5f, 0.1f}, {1, 1.0f, 0.1f}});
    for (int i = 0; i < 10; ++i) reads.push_back({{1, 1.2f, 0.1f}, {2, 0.8f, 0.1f}});
    for (int i = 0; i < 10; ++i) reads.push_back({{2, 1.0f, 0.1f}, {3, 1.3f, 0.1f}});

    auto data = make_data(4, reads);
    const uint32_t T = data.num_refs;
    const size_t   A = data.alignments.size();

    constexpr double eps = 1e-8;

    for (double alpha : {1.0, 2.0}) {           // MLE and MAP
        dart::EMParams p;
        p.use_squarem = false;
        p.use_damage  = false;
        p.alpha_prior = alpha;

        std::vector<double> w(T, 1.0 / T);
        std::vector<double> gamma(A, 0.0);
        double prev = -std::numeric_limits<double>::infinity();

        for (int iter = 0; iter < 60; ++iter) {
            double ll  = dart::e_step(data, p, w.data(), gamma.data());
            double obj = dart::em_objective_from_ll(p, w.data(), T, ll);
            expect(obj >= prev - eps,
                   "alpha=" + std::to_string(alpha) +
                   " objective decreased at iter " + std::to_string(iter) +
                   " (" + std::to_string(prev) + " -> " + std::to_string(obj) + ")",
                   failed);
            dart::m_step(data, p, gamma.data(), w.data());
            prev = obj;
        }
    }

    if (!failed) std::cout << "  PASS\n";
    return failed;
}

// ---------------------------------------------------------------------------
// T1b: SQUAREM objective monotonicity via progress callback
// ---------------------------------------------------------------------------
int test_elbo_monotone_squarem() {
    std::cout << "[T1b] SQUAREM objective monotonicity\n";
    int failed = 0;

    // 3 refs, 50 reads with multi-mapping — forces the extrapolation/safeguard path
    std::vector<std::vector<std::tuple<uint32_t, float, float>>> reads;
    for (int i = 0; i < 25; ++i) reads.push_back({{0, 1.8f, 0.1f}, {1, 1.0f, 0.1f}});
    for (int i = 0; i < 15; ++i) reads.push_back({{1, 1.4f, 0.1f}, {2, 0.9f, 0.1f}});
    for (int i = 0; i < 10; ++i) reads.push_back({{2, 1.0f, 0.1f}, {0, 1.2f, 0.1f}});

    auto data = make_data(3, reads);

    dart::EMParams p;
    p.use_squarem          = true;
    p.use_damage           = false;
    p.squarem_warmup_iters = 2;
    p.max_iters            = 80;
    p.alpha_prior          = 1.0;

    std::vector<double> objs;
    auto cb = [&](const dart::EMIterationDiagnostics& d) {
        objs.push_back(d.objective);
    };

    dart::squarem_em(data, p, nullptr, cb);

    expect(!objs.empty(), "SQUAREM: no iterations recorded", failed);

    constexpr double eps = 1e-7;
    for (size_t i = 1; i < objs.size(); ++i) {
        expect(objs[i] >= objs[i - 1] - eps,
               "SQUAREM objective decreased at iter " + std::to_string(i) +
               " (" + std::to_string(objs[i - 1]) + " -> " + std::to_string(objs[i]) + ")",
               failed);
    }

    if (!failed)
        std::cout << "  PASS (" << objs.size() << " iterations)\n";
    return failed;
}

// ---------------------------------------------------------------------------
// T2: EM parameter recovery (T2 / T3 in roadmap: convergence)
// ---------------------------------------------------------------------------
int test_em_parameter_recovery() {
    std::cout << "[T2] EM parameter recovery\n";
    int failed = 0;

    // Sub-test A: unique mappings — EM must exactly recover ground-truth weights
    {
        std::vector<std::vector<std::tuple<uint32_t, float, float>>> reads;
        for (int i = 0; i < 60; ++i) reads.push_back({{0, 1.0f, 0.1f}});
        for (int i = 0; i < 30; ++i) reads.push_back({{1, 1.0f, 0.1f}});
        for (int i = 0; i < 10; ++i) reads.push_back({{2, 1.0f, 0.1f}});

        dart::EMParams p;
        p.use_damage  = false;
        p.use_squarem = false;
        auto s = dart::em_solve(make_data(3, reads), p);

        constexpr double tol = 1e-4;
        expect(std::abs(s.weights[0] - 0.60) < tol, "unique: w[0] != 0.60", failed);
        expect(std::abs(s.weights[1] - 0.30) < tol, "unique: w[1] != 0.30", failed);
        expect(std::abs(s.weights[2] - 0.10) < tol, "unique: w[2] != 0.10", failed);
    }

    // Sub-test B: symmetric ambiguous — EM must converge to [0.5, 0.5]
    {
        std::vector<std::vector<std::tuple<uint32_t, float, float>>> reads;
        for (int i = 0; i < 40; ++i)
            reads.push_back({{0, 1.0f, 0.1f}, {1, 1.0f, 0.1f}});

        dart::EMParams p;
        p.use_damage  = false;
        p.use_squarem = false;
        auto s = dart::em_solve(make_data(2, reads), p);

        constexpr double tol = 1e-3;
        expect(std::abs(s.weights[0] - 0.5) < tol, "symmetric: w[0] != 0.5", failed);
        expect(std::abs(s.weights[1] - 0.5) < tol, "symmetric: w[1] != 0.5", failed);
    }

    // Sub-test C: skewed ambiguous — dominant ref wins
    // 60 unique reads → ref 0; 40 ambiguous reads → ref 0 or ref 1
    // With equal bit scores on ambiguous reads, ref 0 dominates by unique support.
    {
        std::vector<std::vector<std::tuple<uint32_t, float, float>>> reads;
        for (int i = 0; i < 60; ++i) reads.push_back({{0, 1.0f, 0.1f}});
        for (int i = 0; i < 40; ++i) reads.push_back({{0, 1.0f, 0.1f}, {1, 1.0f, 0.1f}});

        dart::EMParams p;
        p.use_damage  = false;
        p.use_squarem = false;
        p.max_iters   = 200;
        auto s = dart::em_solve(make_data(2, reads), p);

        expect(s.weights[0] > s.weights[1],
               "skewed: ref 0 should outweigh ref 1 (w0=" + std::to_string(s.weights[0]) +
               " w1=" + std::to_string(s.weights[1]) + ")",
               failed);
    }

    if (!failed) std::cout << "  PASS\n";
    return failed;
}

// ---------------------------------------------------------------------------
// T3: FDR control contract (T5 in roadmap)
// ---------------------------------------------------------------------------
int test_fdr_control_contract() {
    std::cout << "[T3] FDR control contract\n";
    int failed = 0;

    // Contract: the threshold returned by find_threshold_for_fdr always satisfies
    //   estimate_fdr(threshold, targets, decoys) <= target_fdr
    // by construction of the scan. We verify this holds for multiple target FDR
    // levels and score distributions.

    // Sub-test A: targets well-separated from decoys
    {
        std::vector<float> targets = {0.90f, 0.80f, 0.70f, 0.60f, 0.50f, 0.40f, 0.30f, 0.20f, 0.10f};
        std::vector<float> decoys  = {0.35f, 0.15f};

        float t   = dart::FDREstimator::find_threshold_for_fdr(0.05f, targets, decoys);
        float fdr = dart::FDREstimator::estimate_fdr(t, targets, decoys);
        expect(fdr <= 0.05f,
               "separated: FDR=" + std::to_string(fdr) + " > 0.05 at t=" + std::to_string(t),
               failed);
        // At minimum, some targets should pass (the high-scoring ones are clean)
        size_t n_pass = 0;
        for (float s : targets) if (s >= t) ++n_pass;
        expect(n_pass >= 4, "separated: fewer than 4 targets pass", failed);
    }

    // Sub-test B: contract holds across multiple FDR levels on a mixed distribution
    {
        // 50 targets linearly spaced [0.4, 1.0]; 50 decoys linearly spaced [0.0, 0.5]
        std::vector<float> targets(50), decoys(50);
        for (int i = 0; i < 50; ++i) {
            targets[i] = 0.4f + 0.6f * (i / 49.0f);
            decoys[i]  = 0.5f * (i / 49.0f);
        }
        for (float tgt : {0.01f, 0.05f, 0.10f, 0.20f}) {
            float t   = dart::FDREstimator::find_threshold_for_fdr(tgt, targets, decoys);
            float fdr = dart::FDREstimator::estimate_fdr(t, targets, decoys);
            expect(fdr <= tgt + 1e-5f,
                   "mixed: FDR=" + std::to_string(fdr) + " > target=" + std::to_string(tgt),
                   failed);
        }
    }

    // Sub-test C: when decoys fully dominate (all decoys > all targets), no threshold
    // in the scan satisfies FDR ≤ target.  The function returns max_threshold; at that
    // point no targets pass, which is the correct conservative behaviour.
    {
        std::vector<float> targets = {0.3f, 0.2f, 0.1f};
        std::vector<float> decoys  = {0.9f, 0.8f, 0.7f, 0.6f};

        float t = dart::FDREstimator::find_threshold_for_fdr(0.05f, targets, decoys);
        size_t n_pass = 0;
        for (float s : targets) if (s >= t) ++n_pass;
        expect(n_pass == 0,
               "dominated: expected 0 targets to pass at t=" + std::to_string(t) +
               " but got " + std::to_string(n_pass),
               failed);
    }

    if (!failed) std::cout << "  PASS\n";
    return failed;
}

// ---------------------------------------------------------------------------
// T4 (simplified): damage-aware EM pi direction test (T4 in roadmap)
//
// Note on the heuristic model: the E-step uses log(pi) + log(p_read) for
// the ancient branch, so the MLE of pi under the heuristic model is NOT
// the true ancient fraction.  Full posterior calibration (SBC) is deferred
// to Phase 4 integration tests; here we test the weaker but verifiable claim:
//   - All-ancient data (damage_score ≈ 1) should pull pi toward 1.
//   - All-modern  data (damage_score ≈ 0) should pull pi toward 0.
//   - The damage-aware path produces valid, finite EM output.
// ---------------------------------------------------------------------------
int test_em_pi_calibration() {
    std::cout << "[T4] damage-aware EM pi direction + output validity\n";
    int failed = 0;

    auto run = [&](float p_read_val, uint32_t n_reads) {
        std::vector<std::vector<std::tuple<uint32_t, float, float>>> reads;
        for (uint32_t i = 0; i < n_reads; ++i)
            reads.push_back({{i % 3, 1.0f, p_read_val}});
        auto data = make_data(3, reads);
        dart::EMParams p;
        p.use_damage  = true;
        p.use_squarem = true;
        p.max_iters   = 150;
        return dart::squarem_em(data, p, nullptr);
    };

    // All-ancient (damage_score=0.99): pi must converge above 0.5
    {
        auto s = run(0.99f, 120);
        expect(s.pi > 0.5,
               "all-ancient: pi=" + std::to_string(s.pi) + " should be > 0.5",
               failed);
        double wsum = 0.0;
        for (double w : s.weights) wsum += w;
        expect(std::abs(wsum - 1.0) < 1e-6, "all-ancient: weights don't sum to 1", failed);
        expect(std::isfinite(s.log_likelihood), "all-ancient: non-finite LL", failed);
    }

    // All-modern (damage_score=0.01): pi must converge below 0.5
    {
        auto s = run(0.01f, 120);
        expect(s.pi < 0.5,
               "all-modern: pi=" + std::to_string(s.pi) + " should be < 0.5",
               failed);
        double wsum = 0.0;
        for (double w : s.weights) wsum += w;
        expect(std::abs(wsum - 1.0) < 1e-6, "all-modern: weights don't sum to 1", failed);
        expect(std::isfinite(s.log_likelihood), "all-modern: non-finite LL", failed);
    }

    if (!failed) std::cout << "  PASS\n";
    return failed;
}

} // anonymous namespace

int main() {
    int total = 0;
    total += test_elbo_monotone_plain_em();
    total += test_elbo_monotone_squarem();
    total += test_em_parameter_recovery();
    total += test_fdr_control_contract();
    total += test_em_pi_calibration();

    if (total == 0) {
        std::cout << "\nAll formal model tests passed.\n";
        return 0;
    }
    std::cerr << "\n" << total << " test(s) FAILED.\n";
    return 1;
}
