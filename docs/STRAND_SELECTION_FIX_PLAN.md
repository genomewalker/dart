# Strand Selection Fix Implementation Plan

## Problem Statement

**Current Status:** Strand selection accuracy is ~50% (essentially random) on synthetic CDS benchmark data.

**Root Cause Analysis:**
1. **92% of reads have stop-codon ties** - Both forward and reverse complement have 0 internal stops in their best frame
2. **Hexamer scoring doesn't break ties effectively** - The hexamer tables were trained for coding vs non-coding discrimination, not sense vs antisense strand discrimination
3. **For true CDS fragments, both strands look "coding-like"** - The antisense strand of a gene still has non-random composition

## Solution Overview

### Two-Pronged Approach

1. **Damage-Based Strand Selection** (Immediate Fix)
   - Use ancient DNA damage patterns to break ties
   - C→T at 5' and G→A at 3' of the *original sequenced molecule*
   - Works best for high-damage samples

2. **Sense/Antisense Hexamer Tables** (Long-term Improvement)
   - Train new hexamer tables that distinguish sense from antisense strands
   - Unify hexamer generation in a single program
   - Generate strand-discriminative log-likelihood ratios

---

## Part 1: Damage-Based Strand Selection

### Biological Rationale

Ancient DNA shows characteristic damage patterns:
- **C→T deamination** at the 5' end of the molecule
- **G→A deamination** at the 3' end (complement of C→T on opposite strand)

For double-stranded libraries:
- If we observe T-enrichment at position 0, that's the 5' end of the *sequenced* molecule
- The true coding strand depends on which strand was sequenced

**Key Insight:** The damage pattern is symmetric with respect to the *original molecule*, but asymmetric with respect to the *sequencing strand*. We can use this to infer strand orientation.

### Implementation Details

#### New Function: `calculate_damage_strand_score()`

```cpp
/**
 * Calculate damage-based strand preference score
 * 
 * Returns positive value if damage patterns favor forward strand,
 * negative if they favor reverse strand.
 * 
 * @param seq Original sequence (as sequenced)
 * @param sample_profile Sample-level damage statistics
 * @return Strand preference score (-1 to +1)
 */
float calculate_damage_strand_score(
    const std::string& seq,
    const SampleDamageProfile& sample_profile);
```

#### Logic:
1. Count T's at 5' terminal (positions 0-2) and A's at 3' terminal
2. Compare observed T/(T+C) and A/(A+G) to sample baseline
3. If 5' shows elevated T → this is likely the 5' end of the molecule → damage signature intact → forward strand favored
4. If 5' shows elevated A (reverse complement would have T) → reverse strand favored

#### Integration into `select_best_per_strand()`

```cpp
// After computing best_fwd and best_rev scores...

// When stop codons don't discriminate (tie)
if (std::abs(min_fwd_stops - min_rev_stops) < 1) {
    // Use damage pattern to break tie
    float damage_pref = calculate_damage_strand_score(seq, sample_profile);
    
    // Add damage preference as tie-breaker
    // Scaled to be meaningful only for close scores
    float score_gap = std::abs(best_fwd.total_score - best_rev.total_score);
    if (score_gap < 0.1f && std::abs(damage_pref) > 0.1f) {
        strand_bonus_fwd += damage_pref * 0.2f;
        strand_bonus_rev -= damage_pref * 0.2f;
    }
}
```

### Files to Modify

- `include/agp/frame_selector.hpp`: Add new function declaration
- `src/frame_selector.cpp`: Implement `calculate_damage_strand_score()`
- Modify `select_best_per_strand()` to use damage-based tie-breaking

### Expected Improvement

- For high-damage samples (>10% terminal damage): +10-20% strand accuracy
- For low-damage samples: minimal change (damage signal too weak)
- Overall: 50% → 60-70% strand accuracy on benchmark data

---

## Part 2: Sense/Antisense Hexamer Tables

### Concept

**Current hexamer tables:**
- Compare P(hexamer | coding) vs P(hexamer | non-coding)
- Both strands of a gene look "coding-like"
- Cannot distinguish sense from antisense

**New sense/antisense tables:**
- Compare P(hexamer | sense strand) vs P(hexamer | antisense strand)
- Train on known CDS sequences:
  - **Sense hexamers**: Hexamers from the coding strand (as annotated)
  - **Antisense hexamers**: Hexamers from the reverse complement of CDS
- Compute log-likelihood ratio: `log(P_sense / P_antisense)`

### Key Differences

| Feature | Sense Strand | Antisense Strand |
|---------|--------------|------------------|
| Start codons | ATG, GTG, TTG at start | CAT, CAC, CAA at start (RC) |
| Stop codons | TAA, TAG, TGA near end | TTA, CTA, TCA near end (RC) |
| Codon bias | Organism-specific | Different (non-coding bias) |
| Wobble position | Biased G/C (synonymous) | Different pattern |

### Implementation

#### Unified Hexamer Extraction Program

Consolidate all hexamer extraction into one program: `agp-train-hexamers`

```bash
./agp-train-hexamers \
    --domain gtdb \
    --cds-dir /path/to/cds/files \
    --output-dir include/agp/ \
    --extract-all
```

**Outputs:**
1. `{domain}_hexamer_table.hpp` - Overall coding vs non-coding
2. `{domain}_positional_hexamer.hpp` - START/INTERNAL/END regions
3. `{domain}_dicodon_phase.hpp` - Phase 0/1/2 for frame selection
4. `{domain}_damage_likelihood.hpp` - Damage pattern tables
5. **NEW:** `{domain}_strand_hexamer.hpp` - Sense vs antisense

#### New Header Structure: `{domain}_strand_hexamer.hpp`

```cpp
#pragma once
#include <array>
#include <cstdint>

namespace agp {
namespace strand {

// Sense strand hexamer frequencies (from CDS)
constexpr std::array<float, 4096> GTDB_SENSE_HEXAMER_FREQ = { ... };

// Antisense strand hexamer frequencies (from RC of CDS)
constexpr std::array<float, 4096> GTDB_ANTISENSE_HEXAMER_FREQ = { ... };

// Pre-computed log-likelihood ratios: log(P_sense / P_antisense)
constexpr std::array<float, 4096> GTDB_STRAND_LLR = { ... };

/**
 * Calculate strand log-likelihood ratio for a sequence
 * Positive = more likely sense strand
 * Negative = more likely antisense strand
 */
float calculate_strand_llr(const char* seq, size_t len, int frame);

} // namespace strand
} // namespace agp
```

#### Extraction Algorithm

```cpp
void process_cds_for_strand_tables(const std::string& cds_seq, ThreadState& state) {
    // Skip if not valid CDS
    if (!has_valid_start(cds_seq) || !has_valid_stop(cds_seq)) return;
    if (cds_seq.length() % 3 != 0) return;
    
    // Extract hexamers from sense strand (the CDS as provided)
    for (size_t i = 0; i + 5 < cds_seq.length(); i += 3) {  // Codon-aligned
        uint32_t code = encode_hexamer(cds_seq.c_str() + i);
        if (code != UINT32_MAX) {
            state.sense_counts[code]++;
            state.sense_total++;
        }
    }
    
    // Extract hexamers from antisense strand (reverse complement)
    std::string antisense = reverse_complement(cds_seq);
    for (size_t i = 0; i + 5 < antisense.length(); i += 3) {
        uint32_t code = encode_hexamer(antisense.c_str() + i);
        if (code != UINT32_MAX) {
            state.antisense_counts[code]++;
            state.antisense_total++;
        }
    }
}
```

### Expected Improvement

- Theoretical: Hexamers should be highly discriminative between sense/antisense
- Expected accuracy improvement: +20-30% on stop-codon-tied reads
- Overall strand accuracy: 50% → 70-80%

---

## Part 3: Combined Strategy

### Final Strand Selection Algorithm

```cpp
std::pair<FrameScore, FrameScore> select_best_per_strand_v2(
    const std::string& seq,
    const SampleDamageProfile& sample_profile,
    const DamageProfile* damage = nullptr) {
    
    // Step 1: Traditional scoring (stops, hexamers, etc.)
    auto [best_fwd, best_rev] = select_best_per_strand(seq, damage);
    
    // Step 2: Compute strand-specific signals
    float stop_diff = count_stops_fwd - count_stops_rev;  // Already computed
    float strand_llr = strand::calculate_strand_llr(seq, best_fwd.frame)
                     - strand::calculate_strand_llr(rc_seq, best_rev.frame);
    float damage_pref = calculate_damage_strand_score(seq, sample_profile);
    
    // Step 3: Weighted combination for strand decision
    float strand_decision = 0.0f;
    
    // Stop codons are definitive when they differ
    if (std::abs(stop_diff) >= 1) {
        strand_decision = stop_diff * 5.0f;  // Strong signal
    }
    
    // Strand LLR is the primary continuous signal
    strand_decision += strand_llr * 1.0f;
    
    // Damage adds additional signal for ancient DNA
    if (sample_profile.max_damage_5prime > 0.05f) {
        strand_decision += damage_pref * 0.5f;
    }
    
    // Step 4: Apply strand decision to scores
    if (strand_decision > 0) {
        best_fwd.total_score += std::min(0.3f, std::abs(strand_decision) * 0.1f);
    } else {
        best_rev.total_score += std::min(0.3f, std::abs(strand_decision) * 0.1f);
    }
    
    return {best_fwd, best_rev};
}
```

### Priority Order

1. **Stop codons** (when different) - 96% accurate
2. **Strand-specific hexamer LLR** (new tables) - 70-80% expected
3. **Damage pattern preference** (ancient DNA specific) - 60-70% on high-damage
4. **Traditional hexamer scoring** (fallback) - 50% baseline

---

## Implementation Roadmap

### Phase 1: Damage-Based Strand Selection (1-2 days)

1. Implement `calculate_damage_strand_score()` in `frame_selector.cpp`
2. Integrate into `select_best_per_strand()` as tie-breaker
3. Test on benchmark data
4. Document in README

### Phase 2: Unified Hexamer Extraction (2-3 days)

1. Refactor `agp_train_hexamers.cpp` to be the single extraction program
2. Add sense/antisense extraction capability
3. Generate new header structure with all table types
4. Test extraction on GTDB data

### Phase 3: Strand-Specific Hexamer Tables (3-5 days)

1. Run extraction for all 8 domains
2. Generate `{domain}_strand_hexamer.hpp` files
3. Implement `strand::calculate_strand_llr()`
4. Integrate into strand selection algorithm
5. Benchmark on synthetic data

### Phase 4: Integration & Testing (2-3 days)

1. Combine all strand signals in weighted scoring
2. Tune weights on validation data
3. Run full benchmark suite
4. Update documentation

---

## Success Metrics

| Metric | Current | Target | Stretch |
|--------|---------|--------|---------|
| Strand accuracy (all reads) | 50% | 70% | 80% |
| Strand accuracy (stop-codon ties) | 50% | 65% | 75% |
| Combined frame+strand accuracy | 46% | 65% | 75% |
| AA correction precision | 27% | 50% | 70% |
| AA correction recall | 22% | 40% | 60% |

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Sense/antisense hexamers not discriminative | Medium | High | Test on subset first |
| Damage signal too weak for modern DNA | High | Medium | Use as tie-breaker only |
| Computational overhead | Low | Low | Pre-compute LLRs |
| Over-fitting to benchmark data | Medium | Medium | Cross-validate |

---

## Notes

### Why This Wasn't a Problem Before

The benchmark data is **synthetic CDS fragments** where:
- Both orientations are valid reading frames (no internal stops)
- The antisense strand also looks "coding-like"

For **real metagenomic data**:
- Most sequences are NOT perfect CDS fragments
- Stop codons are more discriminative
- The ~50% strand error may be less severe

However, improving strand selection will:
1. Improve damage correction accuracy
2. Enable better protein function prediction
3. Reduce false positive rates

### Alternative Approaches Considered

1. **Reference-based strand selection**: Requires database lookup, slow
2. **Start codon detection**: ATG/GTG/TTG at expected position could indicate strand
3. **Ribosome binding site (RBS) detection**: Shine-Dalgarno motif indicates strand
4. **Neural network strand classifier**: Tested, similar performance to hexamer LLR
