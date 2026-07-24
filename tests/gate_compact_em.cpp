// Bitwise gate for compact_em_map (fusion dual-b8247d43).
// (a) DETERMINISM: compact_em_map gamma/gamma_ancient bitwise-identical at 1 vs 32 threads.
// (b) SAME-MATH vs em_solve (OpenMP reference): max|Δgamma| ~float-order, ZERO winner flips.
// Synthetic AlignmentData (deterministic seed; multi-mapping + frozen damage), the pilots'
// config (use_damage=use_alignment_damage_likelihood=true).
#include "dart/em_reassign.hpp"
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace dart;

static AlignmentData make_synth(uint32_t R, uint32_t T) {
    AlignmentData d;
    d.num_reads = R; d.num_refs = T;
    d.read_offsets.assign(R + 1, 0);
    d.ref_lengths.assign(T, 0);
    d.ref_log_len_half.assign(T, 0.0f);
    std::mt19937_64 rng(12345);
    std::uniform_int_distribution<uint32_t> reflen(120, 900), refpick(0, T - 1), degd(1, 6);
    std::uniform_real_distribution<float> bs(20.f, 200.f), ds(0.02f, 0.98f), ll(-3.f, 0.f), al(30.f, 60.f);
    for (uint32_t t = 0; t < T; ++t) { d.ref_lengths[t] = reflen(rng); d.ref_log_len_half[t] = 0.5f * std::log((float)std::max<uint32_t>(d.ref_lengths[t], 1)); }
    std::vector<Alignment> alns;
    for (uint32_t r = 0; r < R; ++r) {
        uint32_t deg = degd(rng);
        for (uint32_t j = 0; j < deg; ++j) {
            Alignment a{};
            a.read_idx = r; a.ref_idx = refpick(rng);
            a.bit_score = bs(rng); a.aln_len = (uint16_t)al(rng);
            a.damage_score = ds(rng); a.damage_ll_a = ll(rng); a.damage_ll_m = ll(rng);
            a.fident = 0.9f;
            alns.push_back(a);
        }
    }
    // CSR: already grouped by read (built in read order)
    d.alignments = std::move(alns);
    for (const auto& a : d.alignments) d.read_offsets[a.read_idx + 1]++;
    for (uint32_t r = 1; r <= R; ++r) d.read_offsets[r] += d.read_offsets[r - 1];
    return d;
}

static void set_threads(int n) {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
}

int main() {
    AlignmentData data = make_synth(200000, 50000);
    EMParams p;
    p.use_damage = true;
    p.use_alignment_damage_likelihood = true;
    p.use_squarem = true;
    p.max_iters = 60; p.min_iters = 8; p.tol = 1e-6;
    printf("gate: A=%zu R=%u T=%u\n", data.alignments.size(), data.num_reads, data.num_refs);

    // (a) determinism 1 vs 32
    set_threads(1);  EMState c1 = compact_em_map(data, p);
    set_threads(32); EMState c32 = compact_em_map(data, p);
    double det_g = 0, det_ga = 0;
    for (size_t i = 0; i < c1.gamma.size(); ++i) {
        det_g = std::max(det_g, std::fabs(c1.gamma[i] - c32.gamma[i]));
        if (!c1.gamma_ancient.empty()) det_ga = std::max(det_ga, std::fabs(c1.gamma_ancient[i] - c32.gamma_ancient[i]));
    }
    bool det_bit = (det_g == 0.0 && det_ga == 0.0);
    printf("(a) DETERMINISM 1v32: max|dgamma|=%.3e max|dgamma_anc|=%.3e -> %s\n",
           det_g, det_ga, det_bit ? "BITWISE-IDENTICAL" : "NOT DETERMINISTIC");

    // (b) same-math vs em_solve (reference, OpenMP)
    set_threads(32);
    EMState ref = em_solve(data, p);
    double dg = 0, dga = 0; size_t win_flip = 0;
    // per-read argmax-gamma winner flip (discrete auth-decision proxy)
    for (uint32_t r = 0; r < data.num_reads; ++r) {
        uint32_t s = data.read_offsets[r], e = data.read_offsets[r + 1];
        if (e == s) continue;
        size_t cw = s, rw = s; double cb = -1, rb = -1;
        for (uint32_t j = s; j < e; ++j) {
            dg = std::max(dg, std::fabs(c32.gamma[j] - ref.gamma[j]));
            if (!c32.gamma_ancient.empty()) dga = std::max(dga, std::fabs(c32.gamma_ancient[j] - ref.gamma_ancient[j]));
            if (c32.gamma[j] > cb) { cb = c32.gamma[j]; cw = j; }
            if (ref.gamma[j] > rb) { rb = ref.gamma[j]; rw = j; }
        }
        if (data.alignments[cw].ref_idx != data.alignments[rw].ref_idx) win_flip++;
    }
    printf("(b) SAME-MATH vs em_solve: max|dgamma|=%.3e max|dgamma_anc|=%.3e winner_flips=%zu/%u\n",
           dg, dga, win_flip, data.num_reads);
    printf("    compact iters=%u pi=%.6f | ref iters=%u pi=%.6f\n",
           c32.iterations, c32.pi, ref.iterations, ref.pi);

    // (c) TAPE == compact (bit-for-bit) + determinism 1v32 for the O_DIRECT tape solver.
    // scratch dir: $GATE_TAPE_DIR, else /scratch/$USER (node-local NVMe), else ".".
    std::string tdir;
    if (const char* e = std::getenv("GATE_TAPE_DIR")) tdir = e;
    else if (const char* u = std::getenv("USER")) tdir = std::string("/scratch/") + u;
    else tdir = ".";
    { std::string mk = "mkdir -p '" + tdir + "' 2>/dev/null"; if (system(mk.c_str())) {} }
    set_threads(1);  EMState t1 = compact_em_map_tape(data, p, tdir);
    set_threads(32); EMState t32 = compact_em_map_tape(data, p, tdir);
    double tdet_g = 0, tdet_ga = 0, tvc_g = 0, tvc_ga = 0;
    for (size_t i = 0; i < t1.gamma.size(); ++i) {
        tdet_g = std::max(tdet_g, std::fabs(t1.gamma[i] - t32.gamma[i]));
        tvc_g  = std::max(tvc_g,  std::fabs(t32.gamma[i] - c32.gamma[i]));
        if (!t1.gamma_ancient.empty()) {
            tdet_ga = std::max(tdet_ga, std::fabs(t1.gamma_ancient[i] - t32.gamma_ancient[i]));
            tvc_ga  = std::max(tvc_ga,  std::fabs(t32.gamma_ancient[i] - c32.gamma_ancient[i]));
        }
    }
    bool tape_det = (tdet_g == 0.0 && tdet_ga == 0.0);
    bool tape_eq  = (tvc_g == 0.0 && tvc_ga == 0.0);
    printf("(c) TAPE det 1v32: max|dgamma|=%.3e max|dgamma_anc|=%.3e -> %s\n",
           tdet_g, tdet_ga, tape_det ? "BITWISE-IDENTICAL" : "NOT DETERMINISTIC");
    printf("    TAPE vs compact (bit-for-bit): max|dgamma|=%.3e max|dgamma_anc|=%.3e -> %s\n",
           tvc_g, tvc_ga, tape_eq ? "IDENTICAL" : "DIVERGENT");

    bool pass = det_bit && (dg < 1e-9) && (win_flip == 0) && tape_det && tape_eq;
    printf("GATE: %s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
