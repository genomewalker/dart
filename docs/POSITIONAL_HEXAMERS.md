# Position-Specific Hexamer Tables

## Overview

Position-specific hexamer tables capture the distinct codon-pair (dicodon) frequency patterns that occur at different positions within coding sequences (CDS). Unlike uniform hexamer frequencies, these tables recognize that:

- **START region** (first 30bp): Enriched for ATG-containing hexamers and translation initiation signals
- **INTERNAL region** (middle): Represents bulk coding sequence patterns
- **END region** (last 30bp before stop): Enriched for stop codon-adjacent patterns

## Rationale

Real genes have position-specific sequence biases:

1. **Start codon context**: The Kozak/Shine-Dalgarno sequences and first few codons after ATG have distinct patterns
2. **Internal coding**: Dominated by codon usage bias for the organism
3. **Stop codon context**: Pre-stop regions often have specific patterns (e.g., rare codons for ribosome pausing)

By scoring different regions with position-appropriate frequencies, we can:
- Improve start codon detection
- Better identify truncated reads (where start/end patterns are missing)
- Distinguish genuine short ORFs from random open reading frames

## Data Sources

### GTDB (Bacteria/Archaea combined)
- **Files**: 143,614 (6,968 archaea + 136,646 bacteria)
- **Valid CDS**: 405.4M sequences
- **Start hexamers**: 3.65B
- **Internal hexamers**: 121.7B
- **End hexamers**: 3.65B
- **Source**: GTDB r220 representative genomes (protein_fna_reps)

### RefSeq Domains

| Domain | Files | Valid CDS | Internal Hexamers |
|--------|-------|-----------|-------------------|
| fungi | 665 | 6.76M | 3.0B |
| protozoa | 121 | 1.16M | 0.6B |
| viral | 15,014 | 0.68M | 0.2B |
| plant | 199 | 9.16M | 4.0B |
| invertebrate | 486 | 12.85M | 5.4B |
| vertebrate_mammalian | 250 | 12.9M | 5.5B |
| vertebrate_other | 504 | 21.26M | 9.0B |

## File Format

Each domain has a header file in `include/agp/`:

```
<domain>_positional_hexamer_table.hpp
```

Each header contains three `constexpr float[4096]` arrays:

```cpp
namespace agp {
    constexpr float START_HEXAMER_FREQ[4096] = { ... };    // First 30bp
    constexpr float INTERNAL_HEXAMER_FREQ[4096] = { ... }; // Middle region
    constexpr float END_HEXAMER_FREQ[4096] = { ... };      // Last 30bp before stop
}
```

Arrays are indexed by encoded hexamer (A=0, C=1, G=2, T=3, 12-bit encoding).

## Extraction Method

The `extract-positional-hexamers` tool processes CDS FASTA files:

1. **Validation**: Only CDS starting with ATG and ending with stop codon (TAA/TAG/TGA) are used
2. **Minimum length**: 93bp (30bp start + 30bp end + 3bp buffer + 30bp internal)
3. **Region extraction**:
   - START: First 30bp, hexamers every 3bp (codon-aligned)
   - END: Last 30bp before stop codon
   - INTERNAL: Everything between START and END regions
4. **Normalization**: Raw counts converted to frequencies

### Usage

```bash
# Single domain
./build/extract-positional-hexamers /path/to/cds/fasta_dir --output fungi --threads 32

# Multiple directories (e.g., combined archaea+bacteria for GTDB)
./build/extract-positional-hexamers /path/to/archaea /path/to/bacteria --output gtdb --threads 64
```

## Integration with Frame Selection

Positional hexamers enhance frame selection by:

1. **Position-aware scoring**: Score the first 30bp of predicted ORF against START frequencies, last 30bp against END frequencies
2. **Truncation detection**: If a read lacks START-like patterns at 5' end, it may be truncated
3. **Confidence weighting**: Reads matching position-specific patterns get higher confidence

### Scoring Formula

For a candidate ORF:

```
score = w_start * LLR(first_30bp, START_FREQ) +
        w_internal * LLR(middle, INTERNAL_FREQ) +
        w_end * LLR(last_30bp, END_FREQ)
```

Where LLR is the log-likelihood ratio against background (uniform 1/4096).

## Key Observations

### Top START hexamers (GTDB)
```
ATGAAA  Start=0.96%  Internal=0.04%  (22x enrichment)
ATGAGC  Start=0.44%  Internal=0.04%  (11x enrichment)
ATGAAG  Start=0.41%  Internal=0.05%  (8x enrichment)
ATGACC  Start=0.40%  Internal=0.05%  (8x enrichment)
AAAAAA  Start=0.38%  Internal=0.10%  (4x enrichment, Shine-Dalgarno?)
```

### Biological interpretation
- ATG-containing hexamers are 8-22x more frequent in START vs INTERNAL
- AAAAAA enrichment in START may reflect purine-rich Shine-Dalgarno sequences
- These position-specific biases are valuable discriminative features

