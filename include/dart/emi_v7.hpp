#pragma once
/**
 * @file emi_v7.hpp
 * @brief EMI v7 on-disk format: self-contained, random-access columnar store of
 *        mmseqs results, built as the best input for EM + damage annotation.
 *
 * Design of record: fable + sol dual-fusion (room dual-7c271f2c), with three
 * open disagreements resolved against the dart source:
 *   - damage_score is HOT: coverage_em.cpp:181/232 weights EM responsibilities
 *     by it as p_read, so it sits beside ref_idx/bit_score.
 *   - qaln/taln are KEPT (cold, ZSTD): cmd_damage_annotate.cpp:360-554 walks the
 *     full alignment to identify every AA SAAV, not just terminal codons.
 *   - obs_codon is DERIVED at annotate time from a per-READ AgdRecord embedded in
 *     the file (not duplicated per alignment); this drops the external .agd dep.
 *
 * Layout (all offsets absolute, little-endian):
 *   [ EMIv7Header (256B) ]
 *   [ Section: BLOCK_DIR      ]  BlockMeta[num_blocks]
 *   [ Section: READ_DEGREE    ]  ZSTD(varint/u16 per-read degree) -> prefix-sum => CSR at load
 *   [ Section: READ_TO_BLOCK  ]  ZSTD(u32 per read -> owning block)
 *   [ Section: BLOCK_CHUNKS   ]  per block: hot chunks then cold chunks (each a ColumnChunk)
 *   [ Section: READ_ANNOT     ]  ZSTD(EmiReadAnnot[num_reads]) — embedded AGD payload
 *   [ Section: NAME_DICT_READ ]  front-coded, ZSTD per NAME_DICT_BLOCK
 *   [ Section: NAME_DICT_REF  ]  front-coded, ZSTD per NAME_DICT_BLOCK
 *   [ Section: ID_TRAILER     ]  dense-ID -> {db uuid, shard, mmseqs key} provenance map
 *   [ Section: PROVENANCE     ]  build args, mmseqs rev, matrix, gap penalties, dbSize
 *   [ EMIv7Footer (trailing magic) ]
 *
 * Blocks are READ-BOUNDARY-ALIGNED: a read's alignments never straddle two
 * blocks (a single high-degree read yields an oversized block). Target
 * EMI_V7_BLOCK_TARGET_ROWS; grow past it only to keep a read whole. This gives
 * one decompression per read AND long monotonic runs for DELTA_ZSTD.
 */

#include "dart/columnar_index.hpp"   // reuse Codec, ColumnStats, ColumnChunk

#include <cstdint>

namespace dart {
namespace emiv7 {

// "AGPEMI02" magic shared with the columnar family; version discriminates v7.
constexpr uint64_t EMI_V7_MAGIC   = 0x3230494D45504741ULL;  // "AGPEMI02"
constexpr uint32_t EMI_V7_VERSION = 7;

// Read-aligned block target. A block ends at a read boundary at/after this many
// rows; one pathological read may exceed it. 65536 keeps DELTA_ZSTD runs long.
constexpr uint32_t EMI_V7_BLOCK_TARGET_ROWS = 65536;

// Name dictionaries are front-coded then ZSTD-compressed in fixed-count blocks.
constexpr uint32_t EMI_V7_NAME_DICT_BLOCK = 4096;

// Degree escape: a per-read degree >= this is stored in a u32 side array.
constexpr uint16_t EMI_V7_DEGREE_ESCAPE = 0xFFFF;

// Dense read ids are u32 (BlockMeta.first_read_idx, READ_ANNOT array index).
// ceiling: 2^32-1 (~4.29e9) distinct reads; alignments are u64 and uncapped.
// upgrade: widen read_idx columns + BlockMeta read fields to u64 (format bump).
constexpr uint64_t EMI_V7_MAX_READS = 0xFFFFFFFFULL;

// Hot-column codec (EMIv7Header.hot_codec). At billions of alignments the hot
// columns (ref_idx+bit_score+damage_score = 12 B/aln) may exceed RAM, so EM runs
// out-of-core: decompress a block, process, evict. That re-decompresses hot
// chunks every iteration, where LZ4 decode speed dominates ZSTD ratio. In-core
// EM (small runs) prefers ZSTD (materialize once, ratio wins). Writer picks per
// projected size; reader honors whatever is stamped here.
enum class HotCodec : uint8_t { ZSTD = 0, LZ4 = 1 };

// Feature flags (EMIv7Header.flags).
constexpr uint32_t EMI_V7_FLAG_SORTED_BY_READ    = 1u << 0;  // rows grouped/sorted by read then score
constexpr uint32_t EMI_V7_FLAG_READ_ALIGNED_BLKS = 1u << 1;  // no read straddles a block
constexpr uint32_t EMI_V7_FLAG_DEGREE_CSR        = 1u << 2;  // CSR stored as per-read degrees
constexpr uint32_t EMI_V7_FLAG_COMPRESSED_DICT   = 1u << 3;  // front-coded + ZSTD name dicts
constexpr uint32_t EMI_V7_FLAG_SELF_CONTAINED    = 1u << 4;  // READ_ANNOT present: no .agd needed
constexpr uint32_t EMI_V7_FLAG_HAS_ALN_STRINGS   = 1u << 5;  // qaln/taln cold columns present
constexpr uint32_t EMI_V7_FLAG_ID_TRAILER        = 1u << 6;  // dense-ID provenance trailer present
constexpr uint32_t EMI_V7_FLAG_HOT_OUT_OF_CORE   = 1u << 7;  // built for block-streaming EM (billions)

// Section identifiers, index into EMIv7Header.sections[].
enum class SectionID : uint32_t {
    BLOCK_DIR      = 0,
    READ_DEGREE    = 1,
    READ_TO_BLOCK  = 2,
    BLOCK_CHUNKS   = 3,
    READ_ANNOT     = 4,
    NAME_DICT_READ = 5,
    NAME_DICT_REF  = 6,
    ID_TRAILER     = 7,
    PROVENANCE     = 8,
    NUM_SECTIONS   = 9
};

// {offset, size} descriptor for one section.
struct SectionDesc {
    uint64_t offset = 0;    // absolute byte offset, 0 = absent
    uint64_t size   = 0;    // on-disk (compressed) byte size
};
static_assert(sizeof(SectionDesc) == 16, "SectionDesc must be 16 bytes");

// 256-byte file header. Byte map: 8+4+4 +8+8+8+4+4 +1+1+6 +4+4+8 + 9*16 + 40 = 256.
struct EMIv7Header {
    uint64_t magic   = EMI_V7_MAGIC;                 // +8  = 8
    uint32_t version = EMI_V7_VERSION;               // +4  = 12
    uint32_t flags   = 0;                            // +4  = 16

    uint64_t num_alignments = 0;                     // +8  = 24  (u64: billions)
    uint64_t num_reads      = 0;                     // +8  = 32
    uint64_t num_refs       = 0;                     // +8  = 40
    uint32_t num_blocks     = 0;                     // +4  = 44
    uint32_t block_target   = EMI_V7_BLOCK_TARGET_ROWS; // +4 = 48

    uint8_t  hot_codec      = 0;   // HotCodec        // +1  = 49
    uint8_t  terminal_depth = 5;   // AGD codons/end  // +1  = 50
    uint8_t  _pad0[6]       = {0};                    // +6  = 56

    // Scoring fingerprint so a reader can assert bit/e-value parity provenance.
    float    d_max          = 0.0f;                  // +4  = 60
    float    lambda         = 0.0f;                  // +4  = 64
    uint64_t scoring_fp     = 0;    // {matrix,gapO,gapE,dbSize,k,spaced} // +8 = 72

    SectionDesc sections[static_cast<size_t>(SectionID::NUM_SECTIONS)];   // +144 = 216

    uint8_t  _reserved[40] = {0};                    // +40 = 256
};
static_assert(sizeof(EMIv7Header) == 256, "EMIv7Header must be 256 bytes");

// Per-block directory entry. Chunk offsets/sizes live in the ColumnChunk array
// that follows in BLOCK_CHUNKS; here we keep the block's row span + read span so
// a reader can locate a read's block and decompress only the columns it needs.
struct BlockMeta {
    uint64_t start_row       = 0;   // global first row (inclusive)
    uint32_t num_rows        = 0;
    uint32_t first_read_idx  = 0;   // dense read id of first row
    uint32_t last_read_idx   = 0;   // dense read id of last row (inclusive)
    uint32_t _pad            = 0;   // explicit: else memberwise copy leaves the
                                    // ABI padding here uninitialized on disk ->
                                    // non-reproducible bytes (a named member is
                                    // value-initialized AND copied by copy-ctor).
    uint64_t chunks_offset   = 0;   // absolute offset of this block's ColumnChunk[NUM_COLUMNS]
};
static_assert(sizeof(BlockMeta) == 32, "BlockMeta must be 32 bytes");

// Embedded per-read AGD payload. Byte-compatible in spirit with AgdRecord but
// keyed by DENSE read id (array position) instead of a hash — no external .agd
// lookup at annotate time. id_hash is retained only to detect a stale pairing.
struct EmiReadAnnot {
    uint64_t id_hash;            // FNV-1a of query id (parity check vs read name dict)
    uint16_t seq_len;            // DNA sequence length
    uint8_t  frame_strand;       // bits[0:1]=frame(0-2), bit[7]=strand
    uint8_t  damage_pct_q;       // round(damage_pct*2), max 200
    uint8_t  p_damaged_q;        // round(p_read_damaged*255)
    uint8_t  n_5prime;           // valid 5' codons (0-5)
    uint8_t  n_3prime;           // valid 3' codons (0-5)
    uint8_t  _pad;
    uint8_t  codons_5prime[5];   // 0-63, 255=invalid
    uint8_t  codons_3prime[5];
    uint8_t  nt_5prime[3];       // 12 nt, 2 bits each
    uint8_t  nt_3prime[3];
};
static_assert(sizeof(EmiReadAnnot) == 32, "EmiReadAnnot must be 32 bytes (mirrors AgdRecord)");

// Footer: trailing magic so truncation/corruption is caught on read.
struct EMIv7Footer {
    uint64_t header_xxh = 0;     // reserved (0); not computed — integrity is trailing magic
    uint64_t magic      = EMI_V7_MAGIC;
};
static_assert(sizeof(EMIv7Footer) == 16, "EMIv7Footer must be 16 bytes");

}  // namespace emiv7
}  // namespace dart
