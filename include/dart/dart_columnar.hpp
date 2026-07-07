#pragma once
// dart_columnar.hpp — .dcx on-disk columnar format for DART damage-annotation
// results, consumed zero-copy by the convergence Rust reader
// (memmap2 + zstd -> Arrow C Data Interface -> DuckDB).
//
// Dependencies: C++17 + zstd only. No arrow/parquet in DART.
// Endianness: little-endian throughout.
//
// DATA MODEL — two joined grains, both sparse/variable shape:
//   Grain A  one row per read x target alignment kept by EM
//   Grain B  one row per COVERED target position inside a record's span
// The A->B join is IMPLICIT: records are sorted by (target_ord, tstart);
// n_obs_i = popcount(presence_i); Grain B columns are concatenated in that
// record order. The prefix-sum of n_obs IS Arrow's List<Struct> offsets
// buffer — recomputed at load from the presence bitmaps, never stored.
//
// Section order:
//   FileHeader | data blocks | BlockIndexEntry[] | dict section | calib section | Footer
// A reader seeks the Footer first (only seek target), validates the index
// adler32, prunes blocks by the trailer zone map, then lazily zstd-decodes.
//
// MULTI-FILE POOLING: one .dcx is one DART run = one library (e.g. one KapK
// library). A formation is many libraries; convergence reads N files and pools
// observations at the same (target_ord, tstart+k) to recover faint per-library
// signal. Each file carries its own calib_contract_id (identity of the damage
// model+params+code that produced p_damage) and a library_label. Convergence
// pools the RAW sufficient statistics (obs_codon, dist5/3, strand) across files
// and may refit one joint damage curve; it MUST NOT average p_damage/p_post
// across differing calib_contract_id without an explicit recalibration. Same
// id => directly poolable. This is why the id lives in the artifact, not a
// filename: renames and concatenation must not silently mix calibrations.

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

namespace dcx {

constexpr uint8_t  FILE_MAGIC[8]  = {'D','A','R','T','D','C','X','\0'};
constexpr uint8_t  BLOCK_MAGIC[4] = {'D','C','X','B'};
constexpr uint8_t  INDEX_MAGIC[4] = {'D','C','X','I'};
constexpr uint8_t  FOOTER_MAGIC[4]= {'D','C','X','E'};

constexpr uint8_t  VERSION_MAJOR = 2;  // bump => bytes reordered/reinterpreted
constexpr uint8_t  VERSION_MINOR = 0;  // bump => only skippable sections / optional columns added
// v2: mandatory Grain-B residue columns RefAa/ObsAa (ASCII); codon columns
// (BestCodon/BestPostQ/ObsCodon) are now capability-gated, not the residue
// carrier. v1 carried the observed residue ONLY via BestCodon, so a codon-free
// write left the consumer unable to recover alt_aa; v2 fixes that at the core.
constexpr uint32_t CODEC_ZSTD    = 1;

// Block = size-bounded pre-zstd, flushed at whichever comes first:
//   BLOCK_TARGET_RAW_BYTES of raw Grain-A+B payload, or a target_ord change
//   that would split a target across blocks past the soft cap. A target is
//   never split unless it alone exceeds the cap. ~8 MiB raw balances zstd
//   ratio against per-target pushdown selectivity at billions of rows.
constexpr uint32_t BLOCK_TARGET_RAW_BYTES = 8u << 20;

// ---- domain quantization (mirror hoki encode_pident/encode_evalue) ----
inline uint16_t enc_unit(float x) {           // [0,1] -> u16
    return static_cast<uint16_t>(std::lround(std::clamp(x, 0.0f, 1.0f) * 65535.0f));
}
inline float dec_unit(uint16_t q) {           // u16 -> [0,1]
    return q * (1.0f / 65535.0f);
}

// ---- codon packing: A=0 C=1 G=2 T=3, b0 | b1<<2 | b2<<4, range 0..63 ----
// Stored as a full u8 (top 2 bits zero) so Arrow adopts it as UInt8 with no
// unpack pass; zstd reclaims the wasted bits.
inline uint8_t pack_codon(uint8_t b0, uint8_t b1, uint8_t b2) {
    return static_cast<uint8_t>(b0 | (b1 << 2) | (b2 << 4));
}

enum class DamageClass : uint8_t { None = 0, CT = 1, GA = 2, Other = 3 };  // 2-bit packed

// Optional-column / layout flags (FileHeader.col_flags and BlockHeader.flags).
enum ColFlag : uint32_t {
    COL_GAMMA_ANCIENT = 1u << 0,   // Grain A: gamma_ancient present (u16)
    COL_EM_MARGIN     = 1u << 1,   // Grain A: em_margin present (u16)
    COL_EM_REF_NEFF   = 1u << 2,   // Grain A: em_ref_neff present (f32)
    COL_QID_DICT      = 1u << 3,   // query_id dictionary-encoded (else delta-varint ordinal)
    COL_DIST3_DERIVED = 1u << 4,   // Grain B: dist_3prime elided = readlen - dist5 - 1
    COL_OBS_CODON     = 1u << 5,   // Grain B: raw observed codon present (u8) alongside best_codon (MAP).
                                   // The calling grain reads obs_codon so DART's damage-corrected MAP
                                   // does not pre-empt the consumer's own real/damage split.
    COL_STRAND        = 1u << 6,   // Grain A: read strand present (1 bit) for the strand-symmetry test
    COL_BEST_CODON    = 1u << 7,   // Grain B: MAP (damage-corrected) codon + its posterior present (panbarley).
                                   // Needs the marginalized-out codon posterior DART discards for KapK, so
                                   // off by default; RefAa/ObsAa (mandatory) carry the residues regardless.
};

// Sample-level semantics (FileHeader.sample_flags). Deamination discrimination
// and coordinate joins are both biased if these are implicit.
enum SampleFlag : uint16_t {
    SAMPLE_LIB_DS       = 1u << 0, // double-stranded library (legit complementary G->A at 3'); else ss
    SAMPLE_COORD_1BASED = 1u << 1, // target_pos is 1-based; else 0-based
};

// Stable column ids for the self-describing block directory. Minor versions
// append new ids; readers skip unknown ids by ColumnDirEntry.byte_size.
enum class ColumnId : uint16_t {
    QueryId = 1, TargetId, TStart, Span,
    GammaQ, PosteriorQ, EmKeep, GammaAncientQ, EmMarginQ, EmRefNeff,
    Presence, BestCodon, BestPostQ, PDamageQ, Dist5, Dist3, DamageClass2,
    Strand, ObsCodon,   // appended (minor bump): read strand (Grain A), raw observed codon (Grain B)
    RefAa, ObsAa,       // v2 mandatory Grain B: reference and observed residue as ASCII byte (Raw u8)
};

enum class ColEncoding : uint16_t {
    Raw = 0, DeltaVarint, Varint, U16Quant, Bitset1, Bits2, CodonU8, F32,
};

#pragma pack(push, 1)

// @0, 96 bytes. Fixed-size; readers seek the first block via first_block_offset
// (never sizeof), so this may grow in a minor bump without moving data.
struct FileHeader {
    uint8_t  magic[8];            // 0   FILE_MAGIC
    uint8_t  version_major;       // 8
    uint8_t  version_minor;       // 9
    uint16_t header_bytes;        // 10  = sizeof(FileHeader) = 96
    uint32_t col_flags;           // 12  ColFlag bitset (which optional columns exist file-wide)
    uint32_t codec;               // 16  CODEC_ZSTD
    uint32_t n_blocks;            // 20
    uint64_t n_records_a;         // 24  total Grain A rows
    uint64_t n_records_b;         // 32  total Grain B observations
    uint64_t first_block_offset;  // 40  = sizeof(FileHeader)
    uint16_t sample_flags;        // 48  SampleFlag bitset (library type, coord base)
    // Calibration identity: hash over the SCORE CONTRACT (model + code/spec
    // revision + params λ/d_max + score semantics/scale + training reference +
    // canonicalization), truncated. All-zero is invalid: the writer refuses to
    // emit p_damage without it. The human-readable calib_model string and the
    // canonical contract blob live in the calib section (see IndexHeader).
    uint8_t  calib_contract_id[16]; // 50
    uint8_t  reserved[30];        // 66
};

// 28 bytes, hoki-shaped. Minimal on purpose: pruning happens at the footer,
// so by the time a reader parses this it has already committed to the block.
struct BlockHeader {
    uint8_t  magic[4];            // 0   BLOCK_MAGIC
    uint8_t  version;             // 4
    uint8_t  flags;               // 5   block-local ColFlag subset
    uint16_t reserved;            // 6
    uint32_t compressed_sz;       // 8   zstd frame size (payload after this header)
    uint32_t raw_sz;              // 12  decompressed size
    uint32_t n_records_a;         // 16
    uint32_t min_target_ord;      // 20  zone map (also in trailer; kept for standalone scan)
    uint32_t max_target_ord;      // 24
};

// Raw (decompressed) payload begins with this, then the column directory,
// then columns in ColumnDirEntry order. Self-describing so minor versions
// add optional columns without changing the block reader's skip logic.
struct RawBlockHeader {
    uint32_t n_records_a;         // 0
    uint32_t n_records_b;         // 4   == sum of per-record popcount(presence)
    uint32_t n_columns;           // 8
    uint32_t column_dir_bytes;    // 12  = n_columns * sizeof(ColumnDirEntry)
};

struct ColumnDirEntry {          // 18 bytes
    uint16_t column_id;           // 0   ColumnId
    uint16_t encoding;            // 2   ColEncoding
    uint32_t elem_count;          // 4   n_records_a for Grain A cols, n_records_b for Grain B
    uint32_t byte_offset;         // 8   from start of raw payload
    uint32_t byte_size;           // 12
    uint16_t reserved;            // 16
};

// Trailer index — the only seek target for pruning. One per block. Carries
// every pushdown axis so a per-target SAAV query opens exactly one block.
struct BlockIndexEntry {
    uint64_t offset;              // 0   BlockHeader byte offset in file
    uint32_t compressed_sz;      // 8
    uint32_t raw_sz;             // 12
    uint32_t n_records_a;        // 16
    uint32_t n_records_b;        // 20  obs_count; MUST equal Σ popcount(presence)
    uint32_t min_target_ord;     // 24  zone map
    uint32_t max_target_ord;     // 28
    uint32_t min_tstart;         // 32
    uint32_t max_tend;           // 36
};

// Index section header (precedes BlockIndexEntry[]).
struct IndexHeader {
    uint8_t  magic[4];           // 0   INDEX_MAGIC
    uint8_t  version;            // 4
    uint8_t  reserved[3];        // 5
    uint32_t n_blocks;           // 8
    uint64_t qid_dict_offset;    // 12  0 if COL_QID_DICT unset
    uint64_t qid_dict_bytes;     // 20
    uint64_t tgt_dict_offset;    // 28  global target_ord -> target_id string registry
    uint64_t tgt_dict_bytes;     // 36
    uint64_t calib_offset;       // 44  -> calib section (calib_model + library_label + contract blob)
    uint64_t calib_bytes;        // 52
};

// @EOF-24. Reader seeks here first (file_size - 24).
struct Footer {
    uint64_t index_offset;       // 0   -> IndexHeader
    uint32_t index_adler32;      // 8   adler32 over IndexHeader + BlockIndexEntry[] bytes
    uint32_t reserved;           // 12
    uint8_t  magic[4];           // 16  FOOTER_MAGIC
    uint32_t footer_bytes;       // 20  = sizeof(Footer) = 24
};

#pragma pack(pop)

// ---- Column layout inside the zstd frame (SoA, column-major) ----
// Grain A (elem_count = n_records_a), in directory order:
//   QueryId       DeltaVarint ordinal   (or dict-id varint if COL_QID_DICT)
//   TargetId      DeltaVarint            target_ord; delta resets nothing —
//                 boundary where delta != 0 is exactly the tstart reset point
//   TStart        DeltaVarint            reset per target
//   Span          Varint
//   GammaQ        U16Quant  (enc_unit)
//   PosteriorQ    U16Quant
//   EmKeep        Bitset1   (LSB-first, Arrow validity-buffer order)
//   [Strand       Bitset1]  (per COL_STRAND) bit=1 => read maps to reference '-'
//                 strand; needed to test strand-symmetry of the substitution
//                 (real variants symmetric; deamination strand-biased).
//   [GammaAncientQ U16Quant] [EmMarginQ U16Quant] [EmRefNeff F32]   (per col_flags)
//   Presence      Bitset1   concatenated: ceil(span_i/8) bytes per record.
//                 Bit k => target position tstart_i + k is observed.
//                 INVARIANT: bits [span_i, 8*ceil(span_i/8)) MUST be zero;
//                 writer masks the tail byte. n_obs_i = popcount(presence_i).
// Grain B (elem_count = n_records_b = Σ n_obs_i), one entry per set bit, in
// (record, bit) order so the reader walks it against the recomputed offsets:
//   RefAa         Raw u8     reference residue at this target position (ASCII).
//   ObsAa         Raw u8     observed residue (ASCII) — the calling grain. ref==obs
//                 => no substitution; a C->T/G->A-reachable ref->obs pair is the
//                 deamination-suspect first filter. RefAa/ObsAa are MANDATORY so
//                 the consumer recovers alt_aa without decoding a codon.
//   [BestCodon    CodonU8]   (per COL_BEST_CODON) MAP (damage-corrected) codon —
//                 panbarley's "most likely codon". Needs DART's marginalized-out
//                 codon posterior, discarded for KapK, hence gated.
//   [ObsCodon     CodonU8]   (per COL_OBS_CODON) RAW observed codon (nt grain).
//                 Sliced from reads.fq + pred.gff ORF bounds at emit; lets a
//                 consumer redo the nt-level transition reasoning if wanted.
//   [BestPostQ    U16Quant]  (per COL_BEST_CODON) MAP posterior; low => tied codons.
//   PDamageQ      U16Quant
//   Dist5         Raw u16
//   [Dist3        Raw u16]   (absent if COL_DIST3_DERIVED)
//   DamageClass2  Bits2     (LSB-first)
//
// The Arrow List<Struct> offsets for Grain B are NOT stored: recompute as the
// running prefix-sum of popcount(presence_i). Decode-time check that the sum
// equals RawBlockHeader.n_records_b (== BlockIndexEntry.n_records_b) — a stray
// high bit in a tail byte otherwise desyncs every subsequent record silently.

// ---- API (implemented in dart_columnar.cpp) ----

struct DcxWriteOptions {
    int      zstd_level   = 19;
    uint32_t col_flags    = 0;                    // ColFlag bitset
    uint16_t sample_flags = 0;                    // SampleFlag bitset (library type, coord base)
    uint32_t block_raw_cap = BLOCK_TARGET_RAW_BYTES;
    // Calibration provenance (see FileHeader.calib_contract_id). calib_contract_id
    // MUST be non-zero or write_dcx throws — p_damage is meaningless without it.
    // The calib section on disk is: u16 model_len | model | u16 label_len |
    // library_label | u32 blob_len | contract_blob (the canonical bytes the id
    // hashes, so the file self-verifies).
    uint8_t     calib_contract_id[16] = {};
    std::string calib_model;                      // open namespaced string, e.g. "dfit/1.4"
    std::string library_label;                    // stable per-library id, e.g. "KapK-L03"
    std::string calib_contract_blob;              // canonical preimage of calib_contract_id
};

// Materialized writer input: Grain A columns (SoA) already sorted by
// (target_ord, tstart), the concatenated per-record presence bitmaps, and the
// concatenated Grain B observation columns. A Batch is one DART shard chunk
// (bounded, fits memory); write_dcx splits it into size-bounded blocks. The
// DART-side adapter that fills a Batch from the annotate scan is where the X
// policy lives: a covered position is a set presence bit with its RAW observed
// codon in obs_codon (never masked to X — masking is the consumer's decision to
// make from evidence, not ours to bake in); best_codon carries the MAP call for
// panbarley. A position the read spans but that has no base call is an unset
// presence bit (uncovered), distinct from an observed-but-damage-suspect one.
struct Batch {
    // ---- Grain A (one entry per kept alignment) ----
    std::vector<uint64_t> query_ord;     // read ordinal (zigzag-delta) or dict id if COL_QID_DICT
    std::vector<uint32_t> target_ord;    // global target dict ordinal (monotonic non-decreasing)
    std::vector<uint32_t> tstart;        // 0-based target start
    std::vector<uint32_t> span;          // aligned length on target; tend = tstart+span-1
    std::vector<float>    gamma;         // [0,1]
    std::vector<float>    posterior;     // [0,1]
    std::vector<uint8_t>  em_keep;       // 0/1
    std::vector<uint8_t>  strand;        // 0='+' 1='-' vs reference; present iff COL_STRAND
    std::vector<float>    gamma_ancient; // present iff COL_GAMMA_ANCIENT
    std::vector<float>    em_margin;     // present iff COL_EM_MARGIN
    std::vector<float>    em_ref_neff;   // present iff COL_EM_REF_NEFF
    // presence bitmaps concatenated in record order; ceil(span_i/8) bytes each,
    // tail bits [span_i, 8*ceil(span_i/8)) MUST be zero. n_obs_i = popcount.
    std::vector<uint8_t>  presence;

    // ---- Grain B (one entry per set presence bit, in (record,bit) order) ----
    std::vector<uint8_t>  ref_aa;        // MANDATORY: reference residue, ASCII byte. Grain-B count reference.
    std::vector<uint8_t>  obs_aa;        // MANDATORY: observed residue, ASCII byte.
    std::vector<uint8_t>  best_codon;    // 0..63, b0|b1<<2|b2<<4 (MAP); present iff COL_BEST_CODON
    std::vector<uint8_t>  obs_codon;     // 0..63 RAW observed codon; present iff COL_OBS_CODON
    std::vector<float>    best_post;     // [0,1] MAP posterior; present iff COL_BEST_CODON
    std::vector<float>    p_damage;      // [0,1]
    std::vector<uint16_t> dist5;
    std::vector<uint16_t> dist3;         // empty iff COL_DIST3_DERIVED
    std::vector<uint8_t>  damage_class;  // 0..3 (DamageClass)

    // ---- dictionaries ----
    std::vector<std::string> target_dict; // target_ord -> target_id
    std::vector<std::string> query_dict;  // used iff COL_QID_DICT

    uint32_t col_flags = 0;               // ColFlag bitset for the optional columns present
};

bool write_dcx(const char* path, const Batch& batch, const DcxWriteOptions& opt);

// Consumer-side (also usable from C++ tests). The Rust reader mirrors this:
// mmap the file, read Footer at EOF-24, adler32-check the index, prune blocks
// by BlockIndexEntry zone maps, then decode a block into borrowed column
// spans. Grain A cols -> Arrow primitive/DictionaryArray; EmKeep/Presence ->
// Arrow boolean buffers; Grain B -> child arrays of a List<Struct> whose
// offsets are the recomputed popcount prefix-sum. All zero-copy over the
// decompressed block buffer, handed to DuckDB via the Arrow C Data Interface.
struct DcxFile;
struct BlockView;

DcxFile   open_dcx_mmap(const void* base, size_t len);
BlockView read_block(const DcxFile& f, uint32_t block_ordinal);

}  // namespace dcx
