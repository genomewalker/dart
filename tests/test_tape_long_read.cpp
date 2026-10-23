// The O_DIRECT EM tape refills a fixed window of records; a read with more alignments than
// the window (GPU short-query searches reach 371k per read) must still be read whole.
// Property: compact_em_map_tape is bit-for-bit identical to the in-RAM compact_em_map on
// data where one read has kTapeChunkRecs + 37,856 alignments.
#include "dart/em_reassign.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using namespace dart;

int main() {
    const uint32_t long_deg = (1u << 18) + 37856;  // past the 262,144-record window
    const uint32_t R = 2001, T = long_deg + 500;
    AlignmentData d;
    d.num_reads = R; d.num_refs = T;
    d.read_offsets.assign(R + 1, 0);
    d.ref_lengths.assign(T, 300);
    d.ref_log_len_half.assign(T, 0.5f * std::log(300.0f));
    std::mt19937_64 rng(7);
    std::uniform_int_distribution<uint32_t> refpick(0, T - 1), degd(1, 6);
    std::uniform_real_distribution<float> bs(20.f, 200.f);
    for (uint32_t r = 0; r < R; ++r) {
        const uint32_t deg = (r == R / 2) ? long_deg : degd(rng);  // long read in the middle
        for (uint32_t j = 0; j < deg; ++j) {
            Alignment a{};
            a.read_idx = r;
            a.ref_idx = (r == R / 2) ? j : refpick(rng);
            a.bit_score = bs(rng); a.aln_len = 40; a.fident = 0.9f;
            d.alignments.push_back(a);
        }
        d.read_offsets[r + 1] = d.alignments.size();
    }
    EMParams p;
    p.use_squarem = true;
    p.max_iters = 12; p.min_iters = 4; p.tol = 1e-6;

    const char* tmp = std::getenv("TMPDIR");
    const std::string dir = (tmp && *tmp) ? tmp : ".";
    AlignmentData d_ram = d;
    EMState ram = compact_em_map(d_ram, p);
    EMState tape = compact_em_map_tape(d, p, dir);

    if (ram.gamma.size() != tape.gamma.size() || ram.gamma.size() != d.alignments.size()) {
        std::printf("FAIL: gamma sizes ram=%zu tape=%zu alns=%zu\n",
                    ram.gamma.size(), tape.gamma.size(), d.alignments.size());
        return 1;
    }
    size_t diff = 0;
    for (size_t i = 0; i < ram.gamma.size(); ++i) diff += (ram.gamma[i] != tape.gamma[i]);
    for (size_t i = 0; i < ram.weights.size(); ++i) diff += (ram.weights[i] != tape.weights[i]);
    std::printf("tape long read: deg=%u alns=%zu differing gamma/weights=%zu\n",
                long_deg, d.alignments.size(), diff);
    return diff == 0 ? 0 : 1;
}
