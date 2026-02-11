# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Standard release build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# Development build with debug symbols
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)

# Incremental rebuild (after code changes)
cmake --build build

# Clean rebuild
rm -rf build && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make -j$(nproc)
```

## Running and Testing

```bash
# Run main tool
./build/agp <input.fasta.gz> -o output.gff --fasta-aa proteins.faa -v

# Validation tools
./build/agp-validate predictions.gff reference.gff
./build/agp-validate-sequences predictions.gff reference.fasta

# Quick test with small dataset (use /scratch/tmp for temporary files)
echo ">test
ATGGCTAGCTAGCTAGCTAGCTAGCTAG" > /scratch/tmp/test.fa
./build/agp /scratch/tmp/test.fa -o /scratch/tmp/test.gff -v
```

## Project Architecture

### Core Pipeline Flow

**Two-Pass Architecture:**
1. **Pass 1** (damage_detection.cpp): Scans entire dataset to build sample-wide damage profile
   - Aggregates terminal nucleotide frequencies (T/(T+C) at 5', A/(A+G) at 3')
   - Estimates exponential decay parameters (λ) via half-life method
   - Detects library type (single-stranded vs double-stranded)
   - Outputs: `SampleDamageProfile` with calibrated damage model

2. **Pass 2** (main.cpp): Gene prediction using damage profile
   - Per-read frame selection with damage-aware scoring
   - Bayesian damage probability calculation
   - Optional Bayesian stop codon correction
   - Outputs: GFF3 + optional FASTA files

### Key Components

**Frame Selection** (frame_selector.cpp/hpp)
- `select_best_per_strand()`: Standard frame selection (6-frame translation + scoring)
- `select_best_per_strand_damage_aware()`: **Damage-aware variant** - adjusts stop codon penalties based on damage probability
  - Uses `calculate_damage_adjusted_stop_multiplier()` to weight stops by likelihood of being damage-induced
  - Prevents false negatives from C→T damage creating spurious stop codons (CAA→TAA, CGA→TGA)
  - Returns `FrameScore` with best forward/reverse frames

**Scoring System** (scoring.cpp)
- Weighted combination of 7+ features:
  - Codon usage bias (0.15 weight)
  - Stop codon penalty (0.28 weight) - **critical for damage handling**
  - Amino acid composition (0.10)
  - Dicodon/hexamer patterns (0.13)
  - Dipeptide frequencies (0.13)
  - GC3 content bonus (0.05)
  - Length bonus (0.05)
  - Damage-frame consistency (variable)
- `calculate_stop_penalty()`: Returns **multiplier** (1.0 = no stops, 0.2 = one stop, 0.05 = multiple stops)
  - **Important**: This is a multiplicative penalty, not additive
  - In damage-aware mode, replaced by `calculate_damage_adjusted_stop_multiplier()` which weights by damage probability

**Damage Model** (damage_model.cpp/hpp)
- `DamageModel`: Encapsulates ancient DNA damage patterns
  - Exponential decay from termini: δ(p) = δ_max · e^(-λp) + δ_background
  - Position-dependent C→T (5' end) and G→A (3' end) rates
  - Profile caching via `create_profile_cached()` - **must be const-qualified** (uses mutable cache)
- `stop_codon_damage_probability()`: Bayesian inference for whether a stop is damage-induced
  - Considers codon position, sequence length, and damage profile
  - Used by damage-aware frame selector
- `correct_stop_codons_bayesian()`: Post-prediction correction of likely damage stops
  - Applied AFTER frame selection (not during)
  - Stores corrected protein in `gene.corrected_protein` (metadata only)

**Sequence I/O** (sequence_io.cpp/hpp)
- `FastaReader`/`FastqReader`: Memory-efficient streaming readers
  - Support gzip via pigz or zcat (auto-detected)
  - **Security**: File paths passed to shell commands MUST be escaped via `shell_escape()`
- `GeneWriter`: GFF3 output with damage metadata (ancient_prob, damage_pct, strand_conf)
- `FastaWriter`:
  - `write_genes_nucleotide()`: Outputs DNA sequences
  - `write_genes_protein()`: Outputs **original** protein sequences (gene.protein)
    - **Never output gene.corrected_protein** - that's computational metadata, not biological sequence

**Data Flow:**
```
FASTQ/FASTA → FastaReader → Pass 1 (DamageDetector) → SampleDamageProfile
                                                              ↓
FASTQ/FASTA → FastaReader → Pass 2 → FrameSelector (damage-aware) → Gene
                                              ↓
                                      DamageModel::correct_stop_codons_bayesian()
                                              ↓
                                      Gene{protein, corrected_protein, ...}
                                              ↓
                                      GeneWriter (GFF3) + FastaWriter
```

## Critical Implementation Details

### Damage-Aware Frame Selection
The damage-aware frame selector was previously broken (produced 0 genes) due to mixing incompatible penalty types:
- `calculate_stop_penalty()` returns a **multiplier** (0.0-1.0)
- Early implementations tried to compute **additive penalties** (-1.5, -0.3, etc.)

**Current correct implementation:**
```cpp
float calculate_damage_adjusted_stop_multiplier(seq, protein, frame, damage_model) {
    // Weights each internal stop by (1 - damage_probability)
    // High damage_prob → low weight → less penalty
    // Returns multiplier compatible with scoring formula
}
```

### Performance Optimizations

**Profile Caching** (damage_model.cpp:469, 789)
- **Bug to avoid**: Never call `create_profile(seq_len)` in per-codon loops → O(N²) complexity
- **Always use**: `create_profile_cached(seq_len)` which caches profiles by sequence length
- The cached version must be `const`-qualified (uses mutable cache member)

**Shell Command Injection** (validate.cpp)
- **Critical security bug**: File paths passed to `popen()` must be escaped
- Use `shell_escape()` function with single-quote escaping:
  ```cpp
  std::string shell_escape(const std::string& arg) {
      std::string escaped = "'";
      for (char c : arg) {
          if (c == '\'') escaped += "'\\''";  // Close, escaped quote, open
          else escaped += c;
      }
      escaped += "'";
      return escaped;
  }
  ```

**SIMD Optimizations** (simd_utils.hpp)
- AVX2/AVX512 support for nucleotide operations
- Automatically detected at compile-time via CMake
- Use `USE_AVX2` / `USE_AVX512` preprocessor flags

### Gene Structure

```cpp
struct Gene {
    std::string sequence;           // DNA sequence
    std::string protein;            // Original protein translation
    std::string corrected_protein;  // Bayesian-corrected protein (metadata only)
    uint32_t start, end;            // 0-based coordinates
    bool is_forward;                // Strand orientation
    float score;                    // Coding probability (0-1)
    float ancient_prob;             // Damage pattern probability (0-1)
    float damage_score;             // Per-read damage percentage (0-100)
    float strand_conf;              // Strand confidence (0-1)
    size_t aa_corrections;          // Number of amino acid corrections
};
```

**Important**:
- `protein` = what's actually encoded in the DNA (use for FASTA output)
- `corrected_protein` = computational inference of original pre-damage sequence (metadata/debugging only)

## Version Management

Version is automatically determined from Git tags:
- Tags format: `v0.2.1` → version `0.2.1`
- Untagged commits: `0.2.1-abc1234` (version + commit hash)
- Generated in `build/include/agp/version.h` via CMake

## Common Pitfalls

1. **Don't mix penalty types**: Stop penalties are multiplicative (0-1), not additive
2. **Always escape shell arguments**: Use `shell_escape()` before passing to `popen()`
3. **Cache damage profiles**: Use `create_profile_cached()` in loops, not `create_profile()`
4. **Output original proteins**: FASTA should contain `gene.protein`, never `gene.corrected_protein`
5. **Const-qualify cached methods**: Methods using mutable cache must be marked `const`
6. **Frame selection timing**: Damage correction happens AFTER frame selection, not during

## Development Workflow

When modifying frame selection or scoring:
1. Read existing implementation in frame_selector.cpp first
2. Understand scoring weights and their interaction
3. Test with both damaged and modern DNA samples
4. Verify gene count doesn't drop to zero
5. Check that damage_prob and damage_pct metrics are reasonable

When adding new damage model features:
1. Update `DamageModel` class in damage_model.hpp
2. Ensure new methods are const-qualified if using cache
3. Add profile caching for any per-position lookups
4. Test with high-damage ancient samples (>10% terminal damage)

## Environment Notes

- Use `/scratch/tmp` for temporary files and test data
- use /projects/caeg/scratch/kbd606/agp/ for the runs
- benchmark dat in /projects/caeg/scratch/kbd606/agp/benchmark_data/

## Current Status (2026-02-01)

### Two-Channel Damage Validation (Channel B)

**Problem Solved**: AGP now distinguishes real C→T damage from compositional variation using two independent signals.

The fundamental limitation of reference-free damage detection is that elevated T/(T+C) at terminal positions can be caused by EITHER real C→T deamination OR natural compositional variation. Channel B solves this by using stop codon conversion as an independent validator.

**Two-Channel Architecture:**

- **Channel A (Nucleotide frequencies)**: Measures T/(T+C) at terminal positions
  - Can be elevated by real damage OR compositional variation
  - Computes `decay_llr` (log-likelihood ratio for exponential decay)
  - Control-adjusted via `delta_llr = decay_llr - ctrl_decay_llr`

- **Channel B (Stop codon conversion)**: Tracks CAA→TAA, CAG→TAG, CGA→TGA conversions
  - Can ONLY be elevated by real C→T damage (the "smoking gun")
  - Computes `stop_decay_llr` for position-dependent stop excess
  - If terminal has more stops than interior → real damage
  - If terminal has fewer/equal stops → compositional artifact

**Decision Logic:**
```
If Channel A fires AND Channel B fires → REAL DAMAGE → report d_max
If Channel A fires BUT Channel B flat/negative → ARTIFACT → d_max = 0
If neither fires → NO DAMAGE → d_max = 0
```

**New JSON output fields:**
```json
{
  "damage": {
    "stop_conversion_baseline": 0.42,    // Interior stop/(pre+stop) ratio
    "stop_decay_llr_5prime": 14956.82,   // Positive = real stop excess
    "stop_amplitude_5prime": 0.148,      // Fitted excess at terminal
    "channel_b_valid": true,             // Sufficient data for Channel B
    "damage_validated": true,            // Both channels agree on damage
    "damage_artifact": false             // Channel A only (false positive)
  }
}
```

**Per-Read Integration:**
When sample-level Channel B identifies no real damage:
- Per-read `ancient_prob` is capped at ~0.15
- Per-read `damage_pct` returns 0
- Prevents false positives from compositional T enrichment

**Benchmark Results:**
| Sample | metaDMG | AGP d_max | Stop LLR | Decision |
|--------|---------|-----------|----------|----------|
| 0267130b40 | 0.5% | **0.0%** | -5,460 | NO DAMAGE |
| 68825e1df0 | 8.4% | **8.7%** | +2,637 | VALIDATED |
| 521588e724 | 31.5% | **25.9%** | +653,872 | VALIDATED |
| dfb2272499 | 55.8% | **67.1%** | +1,712,774 | VALIDATED |

**Implementation Files:**
- `include/agp/frame_selector.hpp`: Channel B fields in `SampleDamageProfile`
- `src/sample_profile.cpp`: Codon tracking and joint decision logic
- `src/main.cpp`: Preservation of Pass 1 Channel B decision
- `src/damage_scoring.cpp`: Per-read probability integration
- `src/damage_probability.cpp`: Per-read damage percentage integration

### Inverted Pattern Detection

**Problem Solved**: AGP detects when reference-free damage detection is unreliable.

Some samples (particularly AT-rich organisms) show an **inverted terminal pattern** where terminal T/(T+C) is **lower** than baseline.

**How it works:**
- Compares terminal T/(T+C) at position 0 against baseline from middle of reads
- Positive gradient = normal damage pattern (detectable)
- Negative gradient = inverted pattern (Channel B will catch this)

### Position-0 Artifact Detection

**Problem Solved**: AGP detects adapter ligation artifacts at position 0.

Many ancient DNA samples show a characteristic pattern where:
- Position 0: T depleted (or near baseline) due to adapter ligation bias
- Position 1+: T enriched with exponential decay (real damage signal)

This artifact can mask or corrupt Channel A damage estimation.

**Detection criteria** (either triggers detection):
1. **Classic**: pos0 shift < -0.5% AND pos1 shift > +1%
2. **Jump pattern**: pos1 shift > +2% AND (pos1 - pos0) > +3%

**When position-0 artifact detected:**
- `terminal_shift_5prime` uses position 1 instead of position 0
- Hexamer-based inversion detection is skipped (hexamers include pos0)
- d_max preferentially uses Channel B (stop codon conversion)

**JSON output:**
```json
{
  "position_0_artifact_5prime": true,
  "terminal_shift_5prime": 0.0694
}
```

**Benchmark validation:**
| Sample | Pos0 Shift | Pos1 Shift | Artifact | AGP d_max |
|--------|------------|------------|----------|-----------|
| 68825e1df0 | -8.3% | +7.4% | YES | 6.27% |
| 521588e724 | +0.02% | +6.9% | YES | 9.06% |
| dfb2272499 | -0.3% | +1.9% | NO | 42.00% |

### AGP vs metaDMG Damage Detection Methodology

**Key Difference:**
- **metaDMG**: Reference-based - aligns reads to reference genomes and counts C→T mismatches
- **AGP**: Reference-free - uses two-channel validation (nucleotide frequencies + stop codon conversion)

**When AGP detection works:**
- Samples with sufficient convertible codon coverage (Channel B valid)
- Both damaged and undamaged samples are correctly classified

**When to use metaDMG instead:**
- Very short reads (<60bp) with insufficient interior region for baseline
- Samples where Channel B lacks sufficient convertible codon data