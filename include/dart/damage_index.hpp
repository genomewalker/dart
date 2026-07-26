#pragma once
/**
 * @file damage_index.hpp
 * @brief Binary damage index format (.agd) for efficient post-mapping annotation.
 *
 * Format overview:
 *   [Header: 64 bytes]
 *   [Hash Table: num_buckets × 4 bytes]
 *   [Records: num_records × 32 bytes]
 *   [Chain: num_records × 4 bytes]
 *
 * Provides O(1) lookup of terminal codon information by read ID hash.
 * Enables detection of synonymous damage (C→T in wobble position).
 */

#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>

namespace dart {

// Magic bytes: "AGD\x01" (AGP Damage index, version 1)
constexpr uint32_t AGD_MAGIC = 0x01444741;  // Little-endian "AGD\x01"
constexpr uint32_t AGD_VERSION = 6;

/**
 * @brief File header (64 bytes, fixed size).
 *
 * v3 adds damage informativeness fields (offsets 33-43) from the reserved area.
 * v4 adds environmental damage channels C/D (offsets 44-49).
 * v2/v3 files are still readable via version-gated accessors in DamageIndexReader.
 */
struct AgdHeader {
    uint32_t magic;              // offset 0, 4 bytes
    uint32_t version;            // offset 4, 4 bytes
    uint64_t num_records;        // offset 8, 8 bytes
    uint64_t num_buckets;        // offset 16, 8 bytes
    float d_max;                 // offset 24, 4 bytes
    float lambda;                // offset 28, 4 bytes
    uint8_t library_type;        // offset 32, 1 byte (0=unknown, 1=ss, 2=ds)
    uint8_t damage_validated;    // offset 33, 1 byte (v3+: 0 or 1)
    uint8_t damage_artifact;     // offset 34, 1 byte (v3+: 0 or 1)
    uint8_t channel_b_valid;     // offset 35, 1 byte (v3+: 0 or 1)
    float stop_decay_llr;              // offset 36, 4 bytes (v3+: Channel B LLR)
    float terminal_shift;              // offset 40, 4 bytes (v3+: terminal T/(T+C) shift)
    uint8_t _reserved1[4];             // offset 44, 4 bytes (alignment pad)
    uint64_t records_compressed_size;  // offset 48, 8 bytes (v4+: ZSTD size of records+chain; 0 = uncompressed)
    uint8_t _reserved[8];              // offset 56, 8 bytes (pad to 64)

    // v5: per-read damage-authentication params for libtaph read_ancient_llr /
    // read_ancient_posterior. These are the ONLY SampleDamageProfile scalars
    // damage_evidence_llr() reads (the pi_pos cubes / tri tables it never touches
    // are consumed earlier by finalize_sample_profile). Annotate reconstructs a
    // minimal SampleDamageProfile from this block; no cube/JSON serialization.
    float    d_max_5prime;             // offset 64
    float    d_max_3prime;             // offset 68
    double   d_max_cooc_point;         // offset 72  (co-occurrence d_max point; <=0 ⇒ non-identifiable)
    float    lambda_5prime;            // offset 80
    float    lambda_3prime;            // offset 84
    float    fit_baseline_5prime;      // offset 88
    float    fit_baseline_3prime;      // offset 92
    double   pi_shape_lambda;          // offset 96  (ref-free per-position decay; preferred when fitted)
    double   pi_shape_baseline;        // offset 104
    double   pi_point;                 // offset 112 (validated damaged fraction; for the posterior prior)
    uint8_t  pi_shape_fitted;          // offset 120 (0/1)
    uint8_t  pi_state;                 // offset 121 (DamageConfidence enum value; DETECTED gates the posterior)
    uint8_t  _reserved_v5[6];          // offset 122, pad to 128

    // v6: library-level empirical-Bayes damage prior. damage_prior_point is the damaged-read
    // FRACTION (pi_damaged = δ/d_max); unlike pi_point it is NOT gated by the DamageConfidence
    // verdict, so it survives LOW_ABUNDANCE and is the prior for per-protein p_protein_damaged.
    // damage_prior_state carries the library's DamageConfidence so the consumer can tier the source.
    float    damage_prior_point;       // offset 128 (library pi_damaged; <0 ⇒ unavailable)
    uint8_t  damage_prior_state;       // offset 132 (DamageConfidence enum value)
    uint8_t  _reserved_v6[3];          // offset 133, pad to 136
};
static_assert(sizeof(AgdHeader) == 136, "AgdHeader must be 136 bytes (v6)");

/**
 * @brief Per-gene record (32 bytes, fixed size).
 *
 * Stores terminal codon information for damage detection.
 * Codons are stored as indices (0-63) using standard genetic code ordering.
 */
struct AgdRecord {
    uint64_t id_hash;            // FNV-1a hash of query ID (without frame suffix)
    uint16_t seq_len;            // DNA sequence length
    uint8_t frame_strand;        // bits[0:1] = frame (0-2), bit[7] = strand (0=fwd, 1=rev)
    uint8_t damage_pct_q;        // Quantized damage_pct: round(damage_pct * 2.0), max 200
    uint8_t p_damaged_q;         // Quantized p_read_damaged: round(p * 255)
    uint8_t n_5prime;            // Number of valid 5' codons stored (0-5)
    uint8_t n_3prime;            // Number of valid 3' codons stored (0-5)
    uint8_t _pad;                // Alignment padding
    uint8_t codons_5prime[5];    // Codon indices from 5' end (0-63, 255=invalid)
    uint8_t codons_3prime[5];    // Codon indices from 3' end (0-63, 255=invalid)
    uint8_t nt_5prime[3];        // Raw packed nucleotides (12 nts, 2 bits each)
    uint8_t nt_3prime[3];        // Raw packed nucleotides from 3' end
};
static_assert(sizeof(AgdRecord) == 32, "AgdRecord must be 32 bytes");

/**
 * @brief Hash table bucket entry (4 bytes).
 *
 * Only stores head of chain; collision resolution uses separate chain array.
 */
struct AgdBucket {
    uint32_t record_offset;      // Index into record array (0xFFFFFFFF = empty)
};
static_assert(sizeof(AgdBucket) == 4, "AgdBucket must be 4 bytes");

// -----------------------------------------------------------------------------
// FNV-1a Hash Implementation
// -----------------------------------------------------------------------------

constexpr uint64_t FNV_OFFSET_BASIS = 14695981039346656037ULL;
constexpr uint64_t FNV_PRIME = 1099511628211ULL;

/**
 * @brief Compute FNV-1a hash of a string.
 */
inline uint64_t fnv1a_hash(std::string_view str) {
    uint64_t hash = FNV_OFFSET_BASIS;
    for (char c : str) {
        hash ^= static_cast<uint8_t>(c);
        hash *= FNV_PRIME;
    }
    return hash;
}

/**
 * @brief One parsed DART ORF id.
 *
 * Grammar, emitted by `dart predict` (src/main.cpp):
 *   <read>_<strand>_<frame>              legacy, before ORF ordinals
 *   <read>_<strand>_<frame>_<ordinal>    current
 * <strand> is '+' or '-'; <ordinal> is the per-read gene number, which makes
 * the id unique when one frame of one read yields more than one ORF.
 *
 * Parsing anchors on the `_<strand>_` marker and reads FORWARD. It must not
 * scan backwards from the end: every previous implementation did, and appending
 * a field silently changed their answer instead of failing.
 */
struct DartOrfId {
    std::string_view read;         ///< the id itself when !ok
    char             strand  = '+';
    int              frame   = 0;
    int              ordinal = -1; ///< -1 when absent (legacy id)
    bool             ok      = false;
};

inline DartOrfId parse_dart_orf_id(std::string_view id) {
    DartOrfId r;
    r.read = id;

    constexpr size_t npos = std::string_view::npos;
    const size_t p = id.rfind("_+_");
    const size_t m = id.rfind("_-_");
    const size_t s = (p == npos) ? m : (m == npos ? p : (p > m ? p : m));
    if (s == npos || s == 0) return r;

    size_t i = s + 3, j = i;
    while (j < id.size() && id[j] >= '0' && id[j] <= '9') ++j;
    if (j == i) return r;
    int frame = 0;
    for (size_t k = i; k < j; ++k) frame = frame * 10 + (id[k] - '0');

    int ordinal = -1;
    if (j < id.size()) {
        if (id[j] != '_') return r;
        const size_t o = j + 1;
        size_t e = o;
        while (e < id.size() && id[e] >= '0' && id[e] <= '9') ++e;
        if (e == o || e != id.size()) return r;
        ordinal = 0;
        for (size_t k = o; k < e; ++k) ordinal = ordinal * 10 + (id[k] - '0');
    }

    r.read    = id.substr(0, s);
    r.strand  = id[s + 1];
    r.frame   = frame;
    r.ordinal = ordinal;
    r.ok      = true;
    return r;
}

/**
 * @brief Strip the strand/frame/ordinal suffix from an ORF id.
 *
 * "read_+_1" and "read_+_1_3" both return "read". Anything that is not an ORF
 * id is returned unchanged.
 */
inline std::string_view strip_dart_suffix(std::string_view id) {
    return parse_dart_orf_id(id).read;
}

// -----------------------------------------------------------------------------
// Codon Encoding
// -----------------------------------------------------------------------------

/**
 * @brief Standard genetic code codon-to-index mapping.
 *
 * Order: TTT=0, TTC=1, TTA=2, TTG=3, TCT=4, ..., GGG=63
 * (First base varies slowest, third base varies fastest)
 */
inline int encode_nucleotide(char nt) {
    switch (nt) {
        case 'T': case 't': return 0;
        case 'C': case 'c': return 1;
        case 'A': case 'a': return 2;
        case 'G': case 'g': return 3;
        default: return -1;  // N or other ambiguous
    }
}

/**
 * @brief Encode a 3-letter codon to index (0-63).
 * @return Codon index, or 255 if invalid (contains N or other).
 */
inline uint8_t encode_codon(const char* codon) {
    int b0 = encode_nucleotide(codon[0]);
    int b1 = encode_nucleotide(codon[1]);
    int b2 = encode_nucleotide(codon[2]);
    if (b0 < 0 || b1 < 0 || b2 < 0) return 255;
    return static_cast<uint8_t>((b0 << 4) | (b1 << 2) | b2);
}

/**
 * @brief Decode codon index back to 3-letter string.
 */
inline void decode_codon(uint8_t idx, char* out) {
    static const char BASES[] = "TCAG";
    if (idx > 63) {
        out[0] = out[1] = out[2] = 'N';
        return;
    }
    out[0] = BASES[(idx >> 4) & 3];
    out[1] = BASES[(idx >> 2) & 3];
    out[2] = BASES[idx & 3];
}

/**
 * @brief Pack 4 nucleotides into a single byte (2 bits each).
 */
inline uint8_t pack_nucleotides_4(const char* nts) {
    uint8_t packed = 0;
    for (int i = 0; i < 4; ++i) {
        int val = encode_nucleotide(nts[i]);
        if (val < 0) val = 0;  // Treat N as T for packing
        packed |= (val << (6 - 2 * i));
    }
    return packed;
}

/**
 * @brief Unpack byte to 4 nucleotides.
 */
inline void unpack_nucleotides_4(uint8_t packed, char* out) {
    static const char BASES[] = "TCAG";
    for (int i = 0; i < 4; ++i) {
        out[i] = BASES[(packed >> (6 - 2 * i)) & 3];
    }
}

// -----------------------------------------------------------------------------
// Frame/Strand Encoding
// -----------------------------------------------------------------------------

/**
 * @brief Encode frame (0-2) and strand into single byte.
 */
inline uint8_t encode_frame_strand(int frame, bool is_reverse) {
    return static_cast<uint8_t>((frame & 0x03) | (is_reverse ? 0x80 : 0x00));
}

/**
 * @brief Decode frame from frame_strand byte.
 */
inline int decode_frame(uint8_t fs) {
    return fs & 0x03;
}

/**
 * @brief Decode strand from frame_strand byte.
 */
inline bool decode_is_reverse(uint8_t fs) {
    return (fs & 0x80) != 0;
}

// -----------------------------------------------------------------------------
// Quantization Helpers
// -----------------------------------------------------------------------------

/**
 * @brief Quantize damage_pct (0-100) to byte (0-200, resolution 0.5%).
 */
inline uint8_t quantize_damage_pct(float pct) {
    if (pct <= 0.0f) return 0;
    if (pct >= 100.0f) return 200;
    return static_cast<uint8_t>(pct * 2.0f + 0.5f);
}

/**
 * @brief Dequantize byte back to damage_pct.
 */
inline float dequantize_damage_pct(uint8_t q) {
    return q * 0.5f;
}

/**
 * @brief Quantize probability (0-1) to byte (0-255).
 */
inline uint8_t quantize_probability(float p) {
    if (p <= 0.0f) return 0;
    if (p >= 1.0f) return 255;
    return static_cast<uint8_t>(p * 255.0f + 0.5f);
}

/**
 * @brief Dequantize byte back to probability.
 */
inline float dequantize_probability(uint8_t q) {
    return q / 255.0f;
}

}  // namespace dart
