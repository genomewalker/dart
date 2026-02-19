#pragma once
// Damage-aware read clustering for per-read damage estimation
//
// Hard-EM-lite approach with 2 rounds:
// Round 0: Frame from coding score alone
// Round 1: Refine using cluster support + terminal damage consistency
//
#include <array>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace agp {
namespace damage_cluster {

// -----------------------------
// AA helpers (shared)
// -----------------------------
static constexpr uint8_t AA_DIM = 21; // 20 AA + stop(*)

inline uint8_t aa_index(char aa) {
    switch (aa) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'D': case 'd': return 2;
        case 'E': case 'e': return 3;
        case 'F': case 'f': return 4;
        case 'G': case 'g': return 5;
        case 'H': case 'h': return 6;
        case 'I': case 'i': return 7;
        case 'K': case 'k': return 8;
        case 'L': case 'l': return 9;
        case 'M': case 'm': return 10;
        case 'N': case 'n': return 11;
        case 'P': case 'p': return 12;
        case 'Q': case 'q': return 13;
        case 'R': case 'r': return 14;
        case 'S': case 's': return 15;
        case 'T': case 't': return 16;
        case 'V': case 'v': return 17;
        case 'W': case 'w': return 18;
        case 'Y': case 'y': return 19;
        case '*': return 20;
        default: return 255;
    }
}

inline char aa_from_index(uint8_t i) {
    static constexpr char LUT[AA_DIM] = {
        'A','C','D','E','F','G','H','I','K','L','M',
        'N','P','Q','R','S','T','V','W','Y','*'
    };
    return (i < AA_DIM) ? LUT[i] : 'X';
}

inline float clampf(float x, float lo, float hi) {
    return std::max(lo, std::min(hi, x));
}

// -----------------------------
// 1) PackedPosting (64-bit)
// -----------------------------
struct PackedPosting {
    using word_t = uint64_t;

    // Bit layout:
    // [0..23]   read_id       (24 bits, up to 16.7M reads)
    // [24..25]  frame_slot    (2 bits, 0..3)
    // [26..35]  aa_pos        (10 bits, 0..1023 AA)
    // [36]      terminal_flag (1 bit)
    // [37]      tolerant_flag (1 bit)
    // [38..63]  reserved      (26 bits)

    static constexpr uint64_t READ_BITS = 24;
    static constexpr uint64_t FRAME_BITS = 2;
    static constexpr uint64_t POS_BITS = 10;

    static constexpr uint64_t READ_SHIFT = 0;
    static constexpr uint64_t FRAME_SHIFT = 24;
    static constexpr uint64_t POS_SHIFT = 26;
    static constexpr uint64_t TERM_SHIFT = 36;
    static constexpr uint64_t TOL_SHIFT = 37;

    static constexpr uint64_t READ_MASK = ((1ULL << READ_BITS) - 1ULL) << READ_SHIFT;
    static constexpr uint64_t FRAME_MASK = ((1ULL << FRAME_BITS) - 1ULL) << FRAME_SHIFT;
    static constexpr uint64_t POS_MASK = ((1ULL << POS_BITS) - 1ULL) << POS_SHIFT;
    static constexpr uint64_t TERM_MASK = 1ULL << TERM_SHIFT;
    static constexpr uint64_t TOL_MASK = 1ULL << TOL_SHIFT;

    static inline word_t pack(
        uint32_t read_id,
        uint8_t frame_slot,
        uint16_t aa_pos,
        bool terminal_flag,
        bool tolerant_flag
    ) {
        return (static_cast<word_t>(read_id) << READ_SHIFT) |
               (static_cast<word_t>(frame_slot & 0x3u) << FRAME_SHIFT) |
               (static_cast<word_t>(aa_pos & 0x03FFu) << POS_SHIFT) |
               (terminal_flag ? TERM_MASK : 0ULL) |
               (tolerant_flag ? TOL_MASK : 0ULL);
    }

    static inline uint32_t read_id(word_t w) {
        return static_cast<uint32_t>((w & READ_MASK) >> READ_SHIFT);
    }

    static inline uint8_t frame_slot(word_t w) {
        return static_cast<uint8_t>((w & FRAME_MASK) >> FRAME_SHIFT);
    }

    static inline uint16_t aa_pos(word_t w) {
        return static_cast<uint16_t>((w & POS_MASK) >> POS_SHIFT);
    }

    static inline bool terminal_flag(word_t w) {
        return (w & TERM_MASK) != 0ULL;
    }

    static inline bool tolerant_flag(word_t w) {
        return (w & TOL_MASK) != 0ULL;
    }

    // node_id for neighbor graph (read x frame-slot)
    static inline uint32_t node_id(word_t w) {
        return (read_id(w) << 2) | frame_slot(w);
    }
};

// -----------------------------
// 2) SeedIndex (CSR-like)
// -----------------------------
struct SeedOccurrence {
    uint64_t seed_hash;
    PackedPosting::word_t posting;
};

class SeedIndex {
public:
    struct PostingRange {
        const PackedPosting::word_t* data = nullptr;
        uint32_t size = 0;
        float idf = 0.0f;
        inline bool empty() const { return size == 0; }
        inline const PackedPosting::word_t* begin() const { return data; }
        inline const PackedPosting::word_t* end() const { return data + size; }
    };

    // Beta-Binomial IDF prior parameters (learned from warmup, not hardcoded)
    //
    // Derivation: Instead of +1/+1 smoothing, use empirical Bayes:
    //   idf(s) = log((N + α + β) / (df_s + α))
    //
    // Estimate (α, β) from warmup seed-frequency moments:
    //   - If overdispersed (var > mean), use α = β = 1/ratio
    //   - Fallback: Jeffreys prior α = β = 0.5
    //
    struct IDFPrior {
        float alpha = 0.5f;  // Jeffreys prior default
        float beta = 0.5f;

        static IDFPrior jeffreys() { return {0.5f, 0.5f}; }

        static IDFPrior from_moments(double mean_df, double var_df) {
            IDFPrior p;
            if (var_df > mean_df && mean_df > 0) {
                // Overdispersion suggests informative prior
                const double ratio = var_df / (mean_df + 1.0);
                p.alpha = static_cast<float>(std::max(0.1, 1.0 / ratio));
                p.beta = p.alpha;
            }
            return p;
        }
    };

    class Builder {
    public:
        explicit Builder(uint32_t num_nodes = 0) : num_nodes_(num_nodes) {}

        inline void set_num_nodes(uint32_t n) { num_nodes_ = n; }
        inline void reserve(size_t n_occ) { occ_.reserve(n_occ); }

        inline void append(uint64_t seed_hash, PackedPosting::word_t posting) {
            occ_.push_back({seed_hash, posting});
        }

        inline size_t size() const { return occ_.size(); }

        // Frequency thresholds:
        //   min_df >= 2 (df=1 cannot form pairwise edges)
        //   max_df derived from idf_min: max_df = (N+α+β)*exp(-idf_min) - α
        void build(SeedIndex& out, uint32_t min_df = 2, uint32_t max_df = 0,
                   IDFPrior prior = IDFPrior::jeffreys()) {
            if (occ_.empty()) {
                out.clear();
                out.num_nodes_ = num_nodes_;
                return;
            }

            // Derive max_df from informativeness floor if not specified
            if (max_df == 0) {
                const float idf_min = 0.5f;  // Seeds with idf < 0.5 are too common
                const float N = static_cast<float>(num_nodes_);
                max_df = static_cast<uint32_t>(std::max(10.0f,
                    (N + prior.alpha + prior.beta) * std::exp(-idf_min) - prior.alpha));
            }

            radix_sort_by_seed_(occ_);

            out.clear();
            out.num_nodes_ = num_nodes_;
            out.idf_prior_ = prior;
            out.seed_hashes_.reserve(occ_.size() / 4);
            out.begins_.reserve(occ_.size() / 4);
            out.lens_.reserve(occ_.size() / 4);
            out.idfs_.reserve(occ_.size() / 4);
            out.postings_.reserve(occ_.size());

            size_t i = 0;
            while (i < occ_.size()) {
                size_t j = i + 1;
                const uint64_t h = occ_[i].seed_hash;
                while (j < occ_.size() && occ_[j].seed_hash == h) ++j;

                const uint32_t df = static_cast<uint32_t>(j - i);
                if (df >= min_df && df <= max_df) {
                    const uint32_t begin = static_cast<uint32_t>(out.postings_.size());
                    for (size_t k = i; k < j; ++k) {
                        out.postings_.push_back(occ_[k].posting);
                    }
                    out.seed_hashes_.push_back(h);
                    out.begins_.push_back(begin);
                    out.lens_.push_back(df);
                    out.idfs_.push_back(idf_(num_nodes_, df, prior));
                }
                i = j;
            }
        }

    private:
        // Beta-Binomial IDF: log((N + α + β) / (df + α))
        static inline float idf_(uint32_t n_nodes, uint32_t df, const IDFPrior& p) {
            return std::log((static_cast<float>(n_nodes) + p.alpha + p.beta) /
                            (static_cast<float>(df) + p.alpha));
        }

        // LSD radix sort on 64-bit seed_hash
        static void radix_sort_by_seed_(std::vector<SeedOccurrence>& v) {
            std::vector<SeedOccurrence> tmp(v.size());
            constexpr size_t RAD = 256;
            std::array<size_t, RAD> count{};
            std::array<size_t, RAD> offset{};

            for (int pass = 0; pass < 8; ++pass) {
                count.fill(0);
                const int shift = pass * 8;
                for (const auto& x : v) {
                    ++count[(x.seed_hash >> shift) & 0xFFu];
                }

                size_t run = 0;
                for (size_t b = 0; b < RAD; ++b) {
                    offset[b] = run;
                    run += count[b];
                }

                for (const auto& x : v) {
                    const uint8_t bucket = static_cast<uint8_t>((x.seed_hash >> shift) & 0xFFu);
                    tmp[offset[bucket]++] = x;
                }

                v.swap(tmp);
            }
        }

        uint32_t num_nodes_ = 0;
        std::vector<SeedOccurrence> occ_;
    };

    inline void clear() {
        seed_hashes_.clear();
        begins_.clear();
        lens_.clear();
        idfs_.clear();
        postings_.clear();
    }

    inline uint32_t num_nodes() const { return num_nodes_; }
    inline size_t num_buckets() const { return seed_hashes_.size(); }
    inline size_t num_postings() const { return postings_.size(); }

    inline PostingRange lookup(uint64_t seed_hash) const {
        const auto it = std::lower_bound(seed_hashes_.begin(), seed_hashes_.end(), seed_hash);
        if (it == seed_hashes_.end() || *it != seed_hash) return {};

        const size_t idx = static_cast<size_t>(it - seed_hashes_.begin());
        return PostingRange{
            postings_.data() + begins_[idx],
            lens_[idx],
            idfs_[idx]
        };
    }

    // Precompute node IDF norms for cosine-like overlap normalization
    std::vector<float> compute_node_idf_norm_sq() const {
        std::vector<float> norm_sq(num_nodes_, 0.0f);
        for (size_t b = 0; b < seed_hashes_.size(); ++b) {
            const float w2 = idfs_[b] * idfs_[b];
            const uint32_t begin = begins_[b];
            const uint32_t len = lens_[b];
            for (uint32_t i = 0; i < len; ++i) {
                const PackedPosting::word_t p = postings_[begin + i];
                const uint32_t n = PackedPosting::node_id(p);
                if (n < norm_sq.size()) norm_sq[n] += w2;
            }
        }
        return norm_sq;
    }

private:
    uint32_t num_nodes_ = 0; // typically n_reads * 4 (2 bits slot space)
    IDFPrior idf_prior_{};   // Beta-Binomial prior used during build
    std::vector<uint64_t> seed_hashes_;
    std::vector<uint32_t> begins_;
    std::vector<uint32_t> lens_;
    std::vector<float> idfs_;
    std::vector<PackedPosting::word_t> postings_;
};

// -----------------------------
// 3) ReadCompactState + 2-bit NT pool
// -----------------------------
struct FrameCandidate {
    uint8_t frame = 0;          // biological frame id (0-5)
    float code_score = 0.0f;    // existing coding score
    float combined_score = 0.0f;
    float posterior = 0.0f;
};

struct ReadCompactState {
    uint32_t nt_word_offset = 0; // offset in global 2-bit pool
    uint16_t nt_len = 0;
    uint8_t n_candidates = 0;    // 2 or 3
    std::array<FrameCandidate, 3> cand{};
};

class ReadCompactStore {
public:
    inline void reserve_reads(size_t n) { reads_.reserve(n); }
    inline void reserve_nt_words(size_t n) { nt_words_.reserve(n); }

    uint32_t add_read(
        const char* seq,
        uint16_t len,
        const std::array<uint8_t,3>& frames,
        const std::array<float,3>& code_scores,
        uint8_t n_candidates
    ) {
        ReadCompactState st;
        st.nt_word_offset = static_cast<uint32_t>(nt_words_.size());
        st.nt_len = len;
        st.n_candidates = static_cast<uint8_t>(std::min<uint8_t>(n_candidates, 3));

        encode_2bit_append_(seq, len, nt_words_);
        for (uint8_t i = 0; i < st.n_candidates; ++i) {
            st.cand[i].frame = frames[i];
            st.cand[i].code_score = code_scores[i];
            st.cand[i].combined_score = code_scores[i];
            st.cand[i].posterior = 0.0f;
        }

        reads_.push_back(st);
        return static_cast<uint32_t>(reads_.size() - 1);
    }

    inline const ReadCompactState& read(uint32_t read_id) const { return reads_[read_id]; }
    inline ReadCompactState& read_mut(uint32_t read_id) { return reads_[read_id]; }
    inline size_t size() const { return reads_.size(); }

    inline char base_at(uint32_t read_id, uint16_t pos) const {
        const auto& r = reads_[read_id];
        if (pos >= r.nt_len) return 'N';
        const uint32_t bit_index = static_cast<uint32_t>(pos) * 2u;
        const uint32_t word_index = r.nt_word_offset + (bit_index >> 6);
        const uint32_t bit_off = bit_index & 63u;
        const uint64_t w = nt_words_[word_index];
        const uint8_t v = static_cast<uint8_t>((w >> bit_off) & 0x3u);
        return decode_base_(v);
    }

    std::string decode_read(uint32_t read_id) const {
        const auto& r = reads_[read_id];
        std::string out;
        out.resize(r.nt_len);
        for (uint16_t i = 0; i < r.nt_len; ++i) out[i] = base_at(read_id, i);
        return out;
    }

private:
    static inline uint8_t encode_base_(char c) {
        switch (c) {
            case 'A': case 'a': return 0u;
            case 'C': case 'c': return 1u;
            case 'G': case 'g': return 2u;
            case 'T': case 't': return 3u;
            default: return 0u; // map N->A in compact store
        }
    }

    static inline char decode_base_(uint8_t v) {
        static constexpr char LUT[4] = {'A','C','G','T'};
        return LUT[v & 0x3u];
    }

    static void encode_2bit_append_(const char* seq, uint16_t len, std::vector<uint64_t>& out_words) {
        const uint32_t n_bits = static_cast<uint32_t>(len) * 2u;
        const uint32_t n_words = (n_bits + 63u) >> 6;
        out_words.resize(out_words.size() + n_words, 0ULL);
        uint64_t* base = out_words.data() + (out_words.size() - n_words);

        for (uint16_t i = 0; i < len; ++i) {
            const uint32_t bit_index = static_cast<uint32_t>(i) * 2u;
            const uint32_t wi = bit_index >> 6;
            const uint32_t bo = bit_index & 63u;
            const uint64_t val = static_cast<uint64_t>(encode_base_(seq[i]));
            base[wi] |= (val << bo);
        }
    }

    std::vector<uint64_t> nt_words_;
    std::vector<ReadCompactState> reads_;
};

// -----------------------------
// 4) NeighborAccumulator (fixed-size, no allocations in hot path)
// -----------------------------
template <size_t CAPACITY = 8192, size_t TOPK = 64>
class NeighborAccumulator {
    static_assert((CAPACITY & (CAPACITY - 1)) == 0, "CAPACITY must be power of 2");

public:
    struct Neighbor {
        uint32_t node_id = 0;
        float score = 0.0f; // normalized overlap score
    };

    NeighborAccumulator() { reset(0, 1.0f); }

    inline void reset(uint32_t query_node_id, float query_norm_sq) {
        query_node_id_ = query_node_id;
        query_norm_sq_ = std::max(1e-6f, query_norm_sq);
        overflow_ = 0;
        keys_.fill(kEmpty);
        vals_.fill(0.0f);
    }

    // Add one shared-seed contribution (typically +idf(seed))
    inline void add(uint32_t node_id, float w_idf) {
        if (node_id == query_node_id_) return;
        uint32_t slot = mix_(node_id) & (CAPACITY - 1);

        for (size_t probe = 0; probe < CAPACITY; ++probe) {
            if (keys_[slot] == kEmpty) {
                keys_[slot] = node_id;
                vals_[slot] = w_idf;
                return;
            }
            if (keys_[slot] == node_id) {
                vals_[slot] += w_idf;
                return;
            }
            slot = (slot + 1u) & (CAPACITY - 1u);
        }
        ++overflow_; // table saturated
    }

    template <typename NeighborNormFn>
    size_t topk(
        NeighborNormFn&& neighbor_norm_sq_fn,
        std::array<Neighbor, TOPK>& out,
        float min_score = 0.0f
    ) const {
        size_t n_out = 0;
        float min_keep = std::numeric_limits<float>::infinity();
        size_t min_idx = 0;

        for (size_t i = 0; i < CAPACITY; ++i) {
            const uint32_t key = keys_[i];
            if (key == kEmpty) continue;

            const float raw = vals_[i];
            if (raw <= 0.0f) continue;

            const float n2 = std::max(1e-6f, neighbor_norm_sq_fn(key));
            const float score = raw / std::sqrt(query_norm_sq_ * n2); // cosine-like
            if (score < min_score) continue;

            if (n_out < TOPK) {
                out[n_out++] = Neighbor{key, score};
                if (score < min_keep) {
                    min_keep = score;
                    min_idx = n_out - 1;
                }
            } else if (score > min_keep) {
                out[min_idx] = Neighbor{key, score};
                // recompute min in out
                min_keep = out[0].score;
                min_idx = 0;
                for (size_t j = 1; j < TOPK; ++j) {
                    if (out[j].score < min_keep) {
                        min_keep = out[j].score;
                        min_idx = j;
                    }
                }
            }
        }

        std::sort(out.begin(), out.begin() + n_out, [](const Neighbor& a, const Neighbor& b) {
            return a.score > b.score;
        });
        return n_out;
    }

    inline uint32_t overflow_count() const { return overflow_; }

private:
    static constexpr uint32_t kEmpty = 0xFFFFFFFFu;
    std::array<uint32_t, CAPACITY> keys_{};
    std::array<float, CAPACITY> vals_{};
    uint32_t query_node_id_ = 0;
    float query_norm_sq_ = 1.0f;
    uint32_t overflow_ = 0;

    static inline uint32_t mix_(uint32_t x) {
        // 32-bit integer mix
        x ^= x >> 16;
        x *= 0x7feb352dU;
        x ^= x >> 15;
        x *= 0x846ca68bU;
        x ^= x >> 16;
        return x;
    }
};

// -----------------------------
// 5) TerminalConsensus (fixed window)
// -----------------------------
template <size_t MAX_W = 24>
class TerminalConsensus {
public:
    struct Site {
        char consensus = 'X';
        float freq = 0.0f;   // consensus count / depth
        float depth = 0.0f;
    };

    inline void reset(size_t query_len_aa, uint8_t window) {
        query_len_aa_ = query_len_aa;
        w_ = static_cast<uint8_t>(std::min<size_t>(window, MAX_W));
        counts_.fill(0.0f);
        depth_.fill(0.0f);
    }

    // side: true=5', false=3'
    inline void add(bool five_prime, uint8_t pos_from_end, char aa, float weight = 1.0f) {
        if (pos_from_end >= w_) return;
        const uint8_t ai = aa_index(aa);
        if (ai == 255) return;

        const uint32_t side = five_prime ? 0u : 1u;
        const uint32_t site = side * MAX_W + pos_from_end;
        depth_[site] += weight;
        counts_[site * AA_DIM + ai] += weight;
    }

    inline Site site(bool five_prime, uint8_t pos_from_end) const {
        Site s;
        if (pos_from_end >= w_) return s;

        const uint32_t side = five_prime ? 0u : 1u;
        const uint32_t site_idx = side * MAX_W + pos_from_end;
        const float d = depth_[site_idx];
        s.depth = d;
        if (d <= 0.0f) return s;

        const uint32_t base = site_idx * AA_DIM;
        float best = -1.0f;
        uint8_t best_i = 255;
        for (uint8_t i = 0; i < AA_DIM; ++i) {
            const float c = counts_[base + i];
            if (c > best) {
                best = c;
                best_i = i;
            }
        }

        s.consensus = aa_from_index(best_i);
        s.freq = (best > 0.0f) ? (best / d) : 0.0f;
        return s;
    }

    inline uint8_t window() const { return w_; }
    inline size_t query_len_aa() const { return query_len_aa_; }

private:
    size_t query_len_aa_ = 0;
    uint8_t w_ = 12;
    // [2 * MAX_W * AA_DIM]
    std::array<float, 2 * MAX_W * AA_DIM> counts_{};
    // [2 * MAX_W]
    std::array<float, 2 * MAX_W> depth_{};
};

// -----------------------------
// 6) DamagePosterior
// -----------------------------
class DamagePosterior {
public:
    struct Input {
        char consensus_aa = 'X';
        char observed_aa = 'X';
        bool five_prime = false;      // position in 5' terminal window
        bool three_prime = false;     // position in 3' terminal window
        float dist_to_end = 0.0f;     // AA distance from nearest relevant end
        float d_max = 0.0f;
        float lambda = 0.0f;
        float consensus_freq = 0.0f;  // from TerminalConsensus
        float consensus_depth = 0.0f;
        float frame_posterior = 1.0f;
    };

    DamagePosterior() {
        init_bg_();
        init_damage_tables_();
    }

    // P(D_i = 1 | read, consensus, priors)
    float posterior(const Input& in) const {
        if (in.observed_aa == in.consensus_aa) return 0.0f;
        if (!in.five_prime && !in.three_prime) return 0.0f;

        const uint8_t c = aa_index(in.consensus_aa);
        const uint8_t o = aa_index(in.observed_aa);
        if (c == 255 || o == 255) return 0.0f;

        const float q_depth = std::min(1.0f, in.consensus_depth / 6.0f);
        const float q_freq = std::max(0.0f, (in.consensus_freq - 0.5f) / 0.5f);
        const float q = q_depth * q_freq;

        float pi = in.d_max * std::exp(-in.lambda * std::max(0.0f, in.dist_to_end));
        pi *= q * clampf(in.frame_posterior, 0.0f, 1.0f);
        pi = clampf(pi, 1e-4f, 0.90f);

        const float m_ct = table_get_(ct_, c, o);
        const float m_ga = table_get_(ga_, c, o);
        float m = 1e-5f;
        if (in.five_prime) m = std::max(m, m_ct);
        if (in.three_prime) m = std::max(m, m_ga);

        const float b = std::max(1e-6f, table_get_(bg_, c, o));

        const float num = pi * m;
        const float den = num + (1.0f - pi) * b;
        return (den > 0.0f) ? (num / den) : 0.0f;
    }

    // X-mask threshold via risk minimization (replaces hardcoded 0.82 - 0.25*d_max)
    //
    // Derivation: Minimize expected risk E[cost] = P(FP)*C_FP + P(FN)*C_FN
    // At threshold τ:
    //   - Mask if posterior > τ
    //   - P(FP) ~ (1-posterior) when we mask
    //   - P(FN) ~ posterior when we don't mask
    //
    // Optimal threshold: τ* = C_FP / (C_FP + C_FN)
    //
    // where C_FN depends on:
    //   - Position-specific damage rate: d_i = d_max * exp(-λ * dist)
    //   - Consensus quality: q_i = min(1, depth/6) * max(0, (freq-0.5)/0.5)
    //   - Frame confidence: frame_post
    //   - Cost ratio: κ_mis (cost of missed damage) / κ_FP (cost of false X)
    //
    struct XMaskParams {
        float d_max = 0.0f;
        float lambda = 0.0f;
        float dist_from_end = 0.0f;
        float depth = 0.0f;
        float freq = 0.0f;
        float frame_post = 0.0f;
        float kappa_mis_over_fp = 1.0f;  // Learned from calibration
    };

    static float xmask_threshold(const XMaskParams& p) {
        // Position-specific damage probability
        const float d_i = p.d_max * std::exp(-std::max(0.0f, p.lambda) * p.dist_from_end);

        // Consensus quality
        const float q_depth = std::min(1.0f, p.depth / 6.0f);
        const float q_freq = std::max(0.0f, (p.freq - 0.5f) / 0.5f);
        const float q_i = q_depth * q_freq;

        // Cost of false negative
        const float c_fn = p.kappa_mis_over_fp * d_i * q_i * std::max(0.0f, p.frame_post);

        // Cost of false positive = 1 (reference)
        constexpr float c_fp = 1.0f;

        // Optimal threshold: τ = C_FP / (C_FP + C_FN)
        return c_fp / (c_fp + c_fn + 1e-6f);
    }

    static bool should_xmask(float post, const XMaskParams& p) {
        return post > xmask_threshold(p);
    }

    // Backward-compatible overload (uses default kappa ratio)
    static bool should_xmask(float post, float depth, float freq, float frame_post, float d_max,
                             float lambda = 0.3f, float dist = 0.0f, float kappa = 1.0f) {
        XMaskParams p;
        p.d_max = d_max;
        p.lambda = lambda;
        p.dist_from_end = dist;
        p.depth = depth;
        p.freq = freq;
        p.frame_post = frame_post;
        p.kappa_mis_over_fp = kappa;
        return should_xmask(post, p);
    }

private:
    // flat 21x21
    std::array<float, AA_DIM * AA_DIM> ct_{};
    std::array<float, AA_DIM * AA_DIM> ga_{};
    std::array<float, AA_DIM * AA_DIM> bg_{};

    static inline size_t idx_(uint8_t c, uint8_t o) {
        return static_cast<size_t>(c) * AA_DIM + static_cast<size_t>(o);
    }

    static inline float table_get_(const std::array<float, AA_DIM * AA_DIM>& t, uint8_t c, uint8_t o) {
        return t[idx_(c, o)];
    }

    static inline void table_set_(std::array<float, AA_DIM * AA_DIM>& t, char c, char o, float v) {
        const uint8_t ci = aa_index(c);
        const uint8_t oi = aa_index(o);
        if (ci != 255 && oi != 255) t[idx_(ci, oi)] = v;
    }

    void init_bg_() {
        bg_.fill(0.002f); // baseline mismatch likelihood
        for (uint8_t i = 0; i < AA_DIM; ++i) {
            bg_[idx_(i, i)] = 1e-6f; // exact match not used as mismatch evidence
        }
    }

    void init_damage_tables_() {
        ct_.fill(1e-5f);
        ga_.fill(1e-5f);

        // CT-like (5' damage patterns)
        table_set_(ct_, 'R', 'W', 0.95f);
        table_set_(ct_, 'R', 'C', 0.50f);
        table_set_(ct_, 'R', '*', 0.90f);
        table_set_(ct_, 'Q', '*', 0.90f);
        table_set_(ct_, 'H', 'Y', 0.95f);
        table_set_(ct_, 'P', 'S', 0.60f);
        table_set_(ct_, 'P', 'L', 0.40f);
        table_set_(ct_, 'T', 'I', 0.85f);
        table_set_(ct_, 'T', 'M', 0.65f);
        table_set_(ct_, 'A', 'V', 0.60f);
        table_set_(ct_, 'S', 'F', 0.80f);
        table_set_(ct_, 'S', 'L', 0.40f);

        // GA-like (3' damage patterns)
        table_set_(ga_, 'E', 'K', 0.95f);
        table_set_(ga_, 'D', 'N', 0.95f);
        table_set_(ga_, 'A', 'T', 0.70f);
        table_set_(ga_, 'G', 'R', 0.55f);
        table_set_(ga_, 'G', 'S', 0.55f);
        table_set_(ga_, 'G', 'E', 0.35f);
        table_set_(ga_, 'G', 'D', 0.35f);
        table_set_(ga_, 'V', 'I', 0.85f);
        table_set_(ga_, 'V', 'M', 0.60f);
        table_set_(ga_, 'R', 'K', 0.50f);
        table_set_(ga_, 'R', 'Q', 0.45f);
        table_set_(ga_, 'R', 'H', 0.45f);
        table_set_(ga_, 'S', 'N', 0.55f);
        table_set_(ga_, 'C', 'Y', 0.70f);
        table_set_(ga_, 'W', '*', 0.85f);
    }
};

} // namespace damage_cluster
} // namespace agp
