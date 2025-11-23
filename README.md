# Ancient Gene Predictor (AGP)

A damage-aware gene predictor for ancient DNA sequences.

## Overview

Ancient DNA carries characteristic post-mortem damage—cytosine deamination causing C→T transitions at 5' ends and G→A transitions at 3' ends. Traditional gene predictors treat these as errors, leading to false stop codons and missed predictions.

AGP treats damage patterns as positive evidence for authenticity while correcting their effects on translation. This enables reliable protein-coding identification in degraded ancient metagenomic samples.

## Quick Start

```bash
# Build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# Run
./agp reads.fastq.gz -o predictions.gff --fasta-aa proteins.faa
```

## Performance

Tested on 262M ancient DNA reads (6.7GB compressed):

| Metric | Value |
|--------|-------|
| Total time | 25 minutes |
| Pass 1 (damage scan) | 6.8 min (643K seq/s) |
| Pass 2 (prediction) | 18.6 min (236K seq/s) |
| Memory | ~3 MB constant |

## Methodology

### Two-Pass Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              INPUT (FASTQ.gz)                               │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  PASS 1: Sample-Level Damage Profiling                                      │
│  ═══════════════════════════════════════                                    │
│                                                                             │
│  For each read, count nucleotides at terminal positions:                    │
│                                                                             │
│  5' end (positions 0-14):           3' end (positions 0-14):                │
│  ┌─────────────────────┐            ┌─────────────────────┐                 │
│  │ T count → damage    │            │ A count → damage    │                 │
│  │ C count → undamaged │            │ G count → undamaged │                 │
│  └─────────────────────┘            └─────────────────────┘                 │
│                                                                             │
│  Aggregate across all reads:                                                │
│  • Position-specific T/(T+C) ratios at 5' end                               │
│  • Position-specific A/(A+G) ratios at 3' end                               │
│  • Baseline nucleotide frequencies from read centers                        │
│  • Codon-position-aware damage rates                                        │
│  • CpG vs non-CpG context damage rates                                      │
│                                                                             │
│  Output: SampleDamageProfile with sample classification                     │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  PASS 2: Gene Prediction with Damage-Aware Scoring                          │
│  ═════════════════════════════════════════════════                          │
│                                                                             │
│  For each read:                                                             │
│                                                                             │
│  1. Six-frame translation (3 forward + 3 reverse complement)                │
│     ┌───────────────────────────────────────────────────────┐               │
│     │  Frame 0: ATG GCT AGC TAG...  →  M  A  S  *           │               │
│     │  Frame 1:  TGG CTA GCT AG...  →  W  L  A              │               │
│     │  Frame 2:   GG CTA GCT AG...  →  G  L  A              │               │
│     └───────────────────────────────────────────────────────┘               │
│                                                                             │
│  2. Frame scoring (weighted combination):                                   │
│     • Codon usage bias (0.15)                                               │
│     • Stop codon penalty (0.28)                                             │
│     • Amino acid composition (0.10)                                         │
│     • Dicodon/hexamer patterns (0.13)                                       │
│     • Dipeptide frequencies (0.13)                                          │
│     • GC3 content bonus (0.05)                                              │
│     • Length bonus (0.05)                                                   │
│     • Damage-frame consistency (variable)                                   │
│                                                                             │
│  3. Per-read damage probability using sample profile calibration            │
│                                                                             │
│  4. Per-read damage percentage (log-likelihood ratio)                       │
│                                                                             │
│  Output: Best frame prediction with damage metrics                          │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         OUTPUT (GFF3 + FASTA)                               │
│                                                                             │
│  • GFF3: coordinates, strand, coding probability, damage metrics            │
│  • Nucleotide FASTA: coding sequences in predicted frame                    │
│  • Protein FASTA: translated amino acid sequences                           │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Sample Classification

Based on terminal damage rates (C→T at 5', G→A at 3'):

| Average Damage | Classification |
|----------------|----------------|
| > 10% | HIGH DAMAGE |
| 5-10% | MODERATE DAMAGE |
| 2-5% | LOW DAMAGE |
| < 2% | MINIMAL DAMAGE |

Classification also considers supporting signals:
- **CpG context enrichment**: Elevated damage in CpG dinucleotides (methylation-driven deamination)
- **Wobble position bias**: T enrichment at codon position 3 (synonymous changes preserved)

### Damage Model

Ancient DNA undergoes cytosine deamination, converting C→U (read as T). The observed damage pattern depends on library preparation:

| Library Type | 5' End | 3' End | Notes |
|--------------|--------|--------|-------|
| Double-stranded | C→T | G→A | Both strands sequenced, complementary pattern |
| Single-stranded | C→T | C→T | Only one strand, damage throughout |

AGP **automatically detects library type** from the damage pattern during Pass 1:
- **Symmetric damage** at both ends (5' C→T ≈ 3' G→A) → **double-stranded**
- **Asymmetric damage** at only one end → **single-stranded**
  - Only 5' C→T (no 3' G→A)
  - Only 3' G→A (no 5' C→T)

The detected library type is shown in verbose output (`-v`). Damage rates decay exponentially from fragment termini:

```
                    5' END                              3' END
                    C → T                               G → A

Position:    0    1    2    3    4    5          -5   -4   -3   -2   -1    0
            ─┬────┬────┬────┬────┬────┬──  ...  ──┬────┬────┬────┬────┬────┬─
Damage %:   │20% │14% │10% │ 7% │ 5% │ 4%        4% │ 5% │ 7% │10% │14% │20%│
            └────┴────┴────┴────┴────┴──        ──┴────┴────┴────┴────┴────┘
             ▲                                                           ▲
             Maximum damage                              Maximum damage
             (overhang region)                           (overhang region)
```

**Exponential decay model:**

```
δ(p) = δ_max · e^(-λp) + δ_background

Where:
  δ(p)         = damage rate at position p from terminus
  δ_max        = maximum damage rate at position 0 (10-40% for ancient samples)
  λ            = decay constant (0.3-0.5, sample-dependent)
  p            = distance from terminus (bases)
  δ_background = sequencing error rate (~1%)
```

**Position-weighted damage detection:**

AGP weights terminal positions exponentially when computing damage signals:

```cpp
weight(p) = e^(-λ × p)

Example with λ = 0.3:
  Position 0:  weight = 1.00
  Position 3:  weight = 0.41
  Position 5:  weight = 0.22
  Position 10: weight = 0.05
```

**Sample-adaptive λ estimation:**

Rather than using a fixed decay constant, AGP estimates λ from the sample during Pass 1 using the half-life method:

```
λ = ln(2) / half_position

Where half_position = position at which damage drops to 50% of maximum
```

This adapts to the sample's preservation state:
- Well-preserved samples: Lower λ (gentler decay, damage extends further into reads)
- Heavily damaged samples: Higher λ (steeper decay, damage concentrated at termini)

The estimated λ values are reported in verbose output:
```
[DEBUG] Decay λ: 5'=0.35 3'=0.28 (half-life: 2.0 / 2.5 bp)
```

### Codon-Aware Damage Detection

C→T damage at the 5' end can create false stop codons from common amino acids:

```
Original Codon    Damaged Codon    Effect
──────────────    ─────────────    ──────────────────────────────
CAA (Gln)    →    TAA              False stop (amber)
CAG (Gln)    →    TAG              False stop (ochre)
CGA (Arg)    →    TGA              False stop (opal)
CGG (Arg)    →    TGG (Trp)        Amino acid change
CCN (Pro)    →    TCN (Ser)        Amino acid change
CTN (Leu)    →    TTN (Phe/Leu)    Amino acid change (often synonymous)
```

AGP detects these **damaged stop codons** as positive evidence for authentic ancient damage rather than treating them as gene-breaking mutations.

**Wobble position enrichment:**

In true coding regions with C→T damage, the third codon position (wobble) accumulates more T's because:
1. Many wobble position changes are synonymous (don't change amino acid)
2. Non-synonymous changes at positions 1-2 are under purifying selection

```
                    Codon Position
                    1      2      3 (wobble)
                   ─┬──────┬──────┬──────┬─
Expected T/(T+C):  │ 50%  │ 50%  │ 50%  │  (no damage)
Observed ancient:  │ 55%  │ 52%  │ 65%  │  (with damage)
                   └──────┴──────┴──────┘
                                    ▲
                          Wobble enrichment = evidence of
                          coding region + authentic damage
```

AGP uses this **wobble enrichment signal** to:
1. Distinguish real coding sequences from random DNA
2. Validate that observed T's are from damage (not random composition)
3. Improve frame selection accuracy

### Per-Read Damage Scoring

Each read receives a **damage percentage** (0-100%) based on a Bayesian log-likelihood ratio:

```
               P(observed bases | ancient, damaged)
log(LR) = log ─────────────────────────────────────
               P(observed bases | modern, undamaged)

For each terminal position:
  • T at 5' end:  log(LR) += log(P_damage / P_background)
  • C at 5' end:  log(LR) += log((1 - P_damage) / (1 - P_background))
  • A at 3' end:  log(LR) += log(P_damage / P_background)
  • G at 3' end:  log(LR) += log((1 - P_damage) / (1 - P_background))

Where P_damage is calibrated from the sample profile.
```

The log-likelihood is normalized against the sample's maximum expected damage to produce a percentage:
- **0-30%**: Likely undamaged (modern contamination or well-preserved)
- **30-70%**: Moderate damage evidence
- **70-100%**: Strong damage evidence (typical ancient)

**Quality score integration:**

For FASTQ input, AGP uses base quality scores (Phred+33) to refine damage detection:

- High quality C/G at terminal positions is strong evidence *against* damage (the original base is intact)
- Low quality bases are treated as less informative (could be sequencing error or damage)

Note: Damaged bases typically appear as high-quality T's because the sequencer correctly identifies the deaminated base. Tools like [mapDamage2](https://academic.oup.com/bioinformatics/article/29/13/1682/184965) can rescale quality scores to reflect damage probability for downstream analysis.

## Installation

**Requirements:**
- C++20 compiler (GCC 10+ or Clang 12+)
- CMake 3.18+
- zlib
- OpenMP (optional)
- pigz (optional, for parallel gzip)

**Build:**
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

## Usage

```bash
# Standard run with FASTA output
agp reads.fastq.gz -o predictions.gff --fasta-nt genes.fna --fasta-aa genes.faa

# Verbose output (shows damage statistics)
agp reads.fastq.gz -o predictions.gff -v

# Adjust sensitivity
agp reads.fastq.gz -o predictions.gff --min-coding-prob 0.5

# Modern DNA (skip damage detection)
agp reads.fastq.gz -o predictions.gff --no-damage
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input` | required | Input FASTA/FASTQ (.gz supported) |
| `-o, --output` | predictions.gff | Output GFF3 file |
| `--fasta-nt` | - | Nucleotide FASTA output |
| `--fasta-aa` | - | Protein FASTA output |
| `--min-coding-prob` | 0.3 | Minimum coding probability |
| `--min-length` | 30 | Minimum sequence length |
| `--threads` | auto | Thread count |
| `--no-damage` | off | Disable damage detection |
| `--no-aggregate` | off | Disable two-pass mode |
| `--both-strands` | off | Output both strands |
| `-v, --verbose` | off | Verbose output |

## Output

### Example Run

```
Ancient Gene Predictor v0.2.0
Input: sample.fastq.gz
Output: predictions.gff
Threads: 64

[Pass 1] Scanning for damage patterns...
  Reads: 262729902 | 5' damage: 19.6% | 3' damage: 20.3% | HIGH DAMAGE | 408.6s
  [DEBUG] Decay λ: 5'=0.26 3'=0.26 (half-life: 2.7 / 2.6 bp)
  [DEBUG] Library type: double-stranded

[Pass 2] Gene prediction...
  Sequences: 262729902 | Genes: 262727878 | 1114.5s (235739 seq/s)

  Ancient Probability Distribution:
    0.9-1.0 |############        |  58.7%
    0.8-0.9 |###                 |  16.3%
    ...

  Per-Read Damage Score (mean: 57.1%):
    50-60%  |####                |  19.6%
    60-70%  |####                |  19.8%
    ...

Done.
```

The `[DEBUG]` lines appear with `-v` (verbose) flag and show:
- **Decay λ**: Estimated exponential decay constant and half-life in base pairs
- **Library type**: Auto-detected as `double-stranded` or `single-stranded`

### GFF3 Format

```
##gff-version 3
read_001  AncientGenePredictor  gene  1  150  0.85  +  .  ID=gene1;ancient_prob=0.73;damage_pct=45.2
read_002  AncientGenePredictor  gene  1  180  0.92  -  .  ID=gene1;ancient_prob=0.81;damage_pct=62.1
```

**Attributes:**
- `ancient_prob` — Probability of ancient damage patterns (0-1)
- `damage_pct` — Per-read damage score (0-100%)

### FASTA Output

Nucleotide (`--fasta-nt`):
```
>read_001 frame=0 strand=+ coding_prob=0.85 ancient_prob=0.73
ATGGCTAGCTAGCTAGCTAGCTAGCTAG...
```

Protein (`--fasta-aa`):
```
>read_001 frame=0 strand=+ coding_prob=0.85 ancient_prob=0.73
MASLLQRLLPKQWERT...
```

## Validation

```bash
# Validate against reference annotations
agp-validate predictions.gff reference.gff

# Validate sequences against reference genome
agp-validate-sequences predictions.gff reference.fasta
```

## Citation

If you use AGP in your research, please cite:

> Fernandez-Guerra, A. (2025). Ancient Gene Predictor: damage-aware gene prediction for ancient metagenomes.

## Author

**Antonio Fernandez-Guerra**
Globe Institute, University of Copenhagen

## License

MIT License
