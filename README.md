# AGP - Ancient Gene Predictor

AGP translates ancient DNA (aDNA) reads into proteins for functional annotation. It predicts protein-coding genes from short, degraded metagenomic sequences while accounting for post-mortem C→T and G→A deamination that corrupts standard translation.

## What AGP Does

Ancient DNA reads are too short and fragmented for genome assembly. AGP translates them directly into proteins that can be searched against functional databases (KEGG, CAZy, viral protein databases) to identify metabolic pathways, enzyme functions, and viral content in ancient samples.

**Primary use cases:**
- Functional profiling of ancient metagenomes (KEGG pathway analysis)
- Carbohydrate-active enzyme (CAZyme) discovery in paleofeces and dental calculus
- Viral detection in archaeological material
- Authentication of ancient proteins via damage patterns

## Features

- **Damage-aware translation**: Adjusts for C→T/G→A deamination that creates false stop codons
- **Per-protein damage scoring**: Identifies which proteins show authentic ancient damage patterns
- **Reference-free damage detection**: No alignment to reference genomes required
- **Two-channel validation**: Distinguishes real damage from natural sequence composition
- **High throughput**: ~20,000 reads/second with SIMD optimization

## Installation

Requires C++20 compiler (GCC 10+ or Clang 12+), CMake 3.18+, and zlib.

```bash
git clone https://github.com/genomewalker/agp.git
cd agp
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

## Quick Start

### Basic Gene Prediction

```bash
agp predict -i reads.fq.gz -o predictions.gff --fasta-aa proteins.faa --adaptive
```

This produces:
- `predictions.gff`: Gene coordinates with damage scores
- `proteins.faa`: Translated proteins for database search

### Functional Profiling Pipeline

Search predicted proteins against KEGG for pathway analysis:

```bash
# 1. Predict genes with damage index
agp predict -i reads.fq.gz -o out.gff \
    --fasta-aa-masked search.faa \
    --damage-index out.agd \
    --adaptive

# 2. Search against KEGG (MMseqs2 with aDNA-optimized settings)
mmseqs easy-search search.faa kegg_genes.faa hits.tsv tmp/ \
    --sub-mat VTML20.out \
    --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln"

# 3. Annotate hits with per-protein damage scores
agp damage-annotate -i hits.tsv --damage-index out.agd -o annotated_hits.tsv
```

The `--fasta-aa-masked` output replaces damage-induced stop codons with X for better database matching.

### Sample Damage Assessment

Check if a sample shows authentic ancient damage:

```bash
agp sample-damage reads.fq.gz
```

Output (JSON):
```json
{
  "d_max": 0.25,
  "lambda": 0.3,
  "damage_validated": true,
  "library_type": "double-stranded"
}
```

## Commands

### `agp predict`

Translates reads into proteins with damage-aware frame selection.

```
Usage: agp predict -i <input> -o <output> [options]

Required:
  -i, --input FILE       Input FASTQ/FASTA (gzip supported)
  -o, --output FILE      Output GFF3

Output:
  --fasta-aa FILE        Protein sequences (observed translation)
  --fasta-aa-masked FILE Protein sequences (damage stops masked as X)
  --fasta-nt FILE        Nucleotide sequences
  --damage-index FILE    Binary index for damage-annotate (.agd)
  --summary FILE         Sample statistics (JSON)

Parameters:
  --adaptive             Adaptive damage correction (recommended)
  --domain NAME          Taxonomic domain: gtdb, fungi, plant, viral,
                         vertebrate_mammalian, vertebrate_other,
                         invertebrate, protozoa (default: gtdb)
  --orf-min-aa N         Minimum protein length (default: 10)
  -t, --threads N        Thread count (default: auto)
```

### `agp sample-damage`

Estimates sample-wide damage rate without gene prediction.

```
Usage: agp sample-damage <input.fq.gz> [options]

Output fields:
  d_max              Maximum damage rate at read termini (0-1)
  lambda             Exponential decay constant
  damage_validated   True if damage confirmed by stop codon signal
  library_type       single-stranded or double-stranded
```

### `agp damage-annotate`

Annotates database hits with per-protein damage evidence. Compares observed proteins to reference proteins to identify damage-consistent amino acid substitutions (R→W, H→Y, Q→*, etc.).

```
Usage: agp damage-annotate -i <hits.tsv> --damage-index <index.agd> -o <output.tsv>

Required:
  -i, --hits FILE        MMseqs2 results (16-column format with qaln/taln)
  --damage-index FILE    AGP damage index from predict
  -o FILE                Output annotated TSV

Optional:
  --protein-summary FILE Per-protein aggregated statistics
  --corrected FILE       Reference-corrected protein FASTA
  --threshold FLOAT      Classification threshold (default: 0.7)
```

Output columns include:
- `p_read`: Per-read damage probability from prediction
- `damage_consistent`: Amino acid substitutions matching damage patterns
- `combined_score`: Weighted damage score (0.80×p_read + 0.40×has_nonsyn + 0.05×has_syn)
- `is_damaged`: Binary classification (combined_score >= threshold)

### `agp damage-profile`

Computes position-wise damage profiles for proteins with sufficient coverage.

```
Usage: agp damage-profile -i <reads.fq> --hits <hits.tsv> -o <profile.tsv.gz>

Options:
  --aggregate        Position-wise summary across all proteins
  --min-reads N      Minimum reads per protein (default: 10)
```

### `agp validate`

Compares predictions against ground truth for benchmarking.

```
Usage: agp validate <predictions.gff> <reference.gff>
```

## Output Formats

### GFF3

```
read_001  AGP  CDS  1  150  0.85  +  0  ID=gene_1;ancient_prob=0.72;damage_pct=15.3;conf=0.88
```

| Attribute | Description |
|-----------|-------------|
| ancient_prob | Probability read shows damage pattern (0-1) |
| damage_pct | Estimated damage percentage |
| conf | Prediction confidence from length and coding score |

### Damage Index (.agd)

Binary format for O(1) per-read damage lookup. Contains sequence length, frame, strand, damage probability, and terminal codon information. Required by `damage-annotate`.

### Annotated Hits TSV

Per-read output from `damage-annotate`:

| Column | Description |
|--------|-------------|
| p_read | Per-read damage probability from AGP predict |
| ct_sites | C→T substitution count at protein level |
| ga_sites | G→A substitution count at protein level |
| combined_score | Weighted damage score for classification |
| is_damaged | 1 if combined_score >= threshold |
| syn_5prime | Synonymous C→T at 5' terminus (from damage index) |
| syn_3prime | Synonymous G→A at 3' terminus (from damage index) |

## How It Works

### Two-Pass Architecture

**Pass 1** scans all reads to estimate sample-wide damage:
- Terminal nucleotide frequencies: T/(T+C) at 5', A/(A+G) at 3'
- Exponential decay fitting: δ(p) = δ_max × e^(-λp)
- Stop codon conversion rates for validation

**Pass 2** predicts genes using the calibrated damage model:
- Six-frame translation with damage-aware scoring
- Bayesian stop codon probability based on position
- Frame selection weighted by codon usage, hexamer patterns, and damage consistency

### Two-Channel Damage Validation

Reference-free damage detection has a fundamental limitation: elevated T/(T+C) at read termini could be real C→T damage OR natural sequence composition.

AGP solves this with two independent signals:

**Channel A (Nucleotide frequencies)**: Measures T/(T+C) ratio at each position from 5' terminus. Real damage shows exponential decay from ~30% at position 0 to baseline (~25%) by position 15.

**Channel B (Stop codon conversion)**: Tracks CAA→TAA, CAG→TAG, CGA→TGA conversions. These can ONLY be elevated by real C→T damage. Interior reads provide baseline; terminal excess indicates authentic damage.

Decision logic:
- Channel A fires AND Channel B fires → **Real damage** (report d_max)
- Channel A fires BUT Channel B flat → **Compositional artifact** (d_max = 0)
- Neither fires → **No damage** (d_max = 0)

### Per-Protein Damage Scoring

After database search, `damage-annotate` identifies damage-consistent amino acid substitutions by comparing query (observed) to target (reference) proteins:

| Reference | Observed | Damage Type | Codon Change |
|-----------|----------|-------------|--------------|
| R | W | C→T (5') | CGG→TGG |
| Q | * | C→T (5') | CAA→TAA |
| H | Y | C→T (5') | CAC→TAC |
| E | K | G→A (3') | GAA→AAA |
| D | N | G→A (3') | GAC→AAC |
| R | K | G→A (3') | AGA→AAA |

Positional probability weights these substitutions: sites near termini (where damage is expected) score higher than interior sites.

## Performance

Benchmarked on synthetic ancient DNA with known damage patterns (KapK dataset, 3.1M reads):

| Metric | Value |
|--------|-------|
| Read-level damage detection AUC | 0.78 |
| Per-protein damage precision | 73% (any damage), 67% (AA-changing) |
| Classification precision at threshold 0.7 | 92% |
| Classification recall at threshold 0.7 | 81% |
| Throughput | ~20,000 reads/second |

Combined scoring formula:
```
combined_score = 0.80 × p_read + 0.40 × has_nonsyn + 0.05 × has_syn
```

## Validation

Validating ancient DNA gene prediction is challenging because there is no single gold standard. Different benchmarks test different capabilities, and real archaeological samples lack complete ground truth. We therefore use a combination of synthetic data with known provenance and real samples validated against established reference-based methods.

AGP performs three distinct tasks that require independent validation: (1) gene prediction—selecting the correct reading frame and strand from six possibilities; (2) sample-wide damage estimation—quantifying the overall deamination rate for quality control and downstream filtering; and (3) per-read damage classification—identifying which individual proteins carry authentic ancient damage signatures. Each task has different success criteria and is evaluated against appropriate ground truth.

Our validation uses two complementary data sources. The **KapK synthetic community** consists of gargammel-simulated ancient DNA reads from 58 ancient and 42 modern reference genomes, with damage rates ranging from 0% to 35%. Each read has a known gene assignment, true reading frame, true strand orientation, and per-position damage annotation—enabling precise evaluation of both gene prediction accuracy and damage detection. For sample-level damage estimation, we validate against **31 real archaeological metagenomes** with reference-based damage estimates from metaDMG, the current standard for ancient DNA authentication.

### Read-Level Damage Classification

Beyond sample-level estimation, AGP assigns a per-read damage probability (p_read) reflecting the extent to which each read's terminal nucleotide patterns are consistent with ancient DNA deamination. After database search, the `damage-annotate` command combines p_read with alignment-derived evidence—specifically, amino acid substitutions characteristic of C→T and G→A damage—to classify each protein as ancient or modern. This combined score is validated against KapK ground truth where each read has a known damage status from gargammel simulation.

**Methodology**: The evaluation pipeline runs AGP predict with `--adaptive --damage-index` to create a binary damage index containing p_read for each predicted protein. Proteins are then searched against the KEGG database using MMseqs2 with the VTML20 substitution matrix (optimized for damaged sequences). The `damage-annotate` command computes a combined score:

```
combined_score = 0.80 × p_read + 0.40 × has_nonsyn + 0.05 × has_syn
```

Where `has_nonsyn` indicates whether the alignment contains non-synonymous substitutions consistent with damage (D→N, E→K, H→Y, etc.) and `has_syn` indicates synonymous damage patterns. The combined score is evaluated against binary ground truth (read contains any C→T or G→A event) using AUC-ROC.

**Per-sample results** (10 KapK synthetic samples, 3.1M reads total):

| Sample | Reads | AUC-ROC |
|--------|-------|---------|
| 119_B3...-24 | 94,407 | 0.753 |
| 119_B3...-25 | 21,737 | 0.779 |
| 119_B3...-27 | 69,152 | 0.752 |
| 119_B3...-29 | 57,252 | 0.768 |
| 69_B2...-31 | 304,896 | 0.759 |
| 69_B2...-33 | 377,497 | 0.779 |
| 69_B2...-34 | 799,059 | 0.797 |
| 69_B2...-35 | 668,142 | 0.782 |
| 69_B2...-36 | 608,086 | 0.789 |
| 69_B2...-37 | 148,617 | 0.795 |
| **Weighted mean** | **3,148,845** | **0.78** |

**Classification metrics**: At a threshold of 0.7, this achieves 92% precision and 81% recall for identifying damaged reads. Post-mapping protein-level evaluation yields precision 73% for any damage detection and 67% for AA-changing damage, with 86% recall.

**Information-theoretic constraints**: A single read carries fundamentally limited damage signal because C→T deamination produces thymine indistinguishable from natural thymine without a reference sequence. Consider a single terminal position with d=30% damage rate:

- P(T | damaged read) = 0.25 + 0.75 × 0.30 = 0.475
- P(T | undamaged read) = 0.25

The best achievable single-site AUC is approximately 0.61. The achievable multi-site AUC is bounded by the number of informative terminal positions (typically 5-10 per read end) and the fact that 68% of damaged reads contain exactly one C→T event—a single thymine that could plausibly be natural.

AGP achieves AUC 0.78 by aggregating evidence across all terminal positions on both read ends and applying Bayesian posterior inference with per-read composition baselines. This result approaches the information-theoretic ceiling for reference-free single-read classification. Higher discrimination requires either reference alignment (which metaDMG provides for sample-level estimation) or multi-read aggregation across proteins (which `damage-annotate` provides through the combined score).

**Limitations**: KapK defines "damaged" as any C→T or G→A event; a read with a single synonymous substitution in a non-terminal position is technically damaged but may be indistinguishable from an undamaged read at the protein level. Reads without database hits rely entirely on p_read (lower single-read AUC approximately 0.60). These results derive from synthetic ground truth with controlled damage parameters; real archaeological samples have variable damage base rates and may contain damage patterns not captured by the gargammel model.

### Sample-Wide Damage Validation

Authenticating ancient DNA requires estimating the sample-wide damage rate to establish that a specimen is genuinely ancient rather than modern contamination. AGP performs this estimation without reference genome alignment using two-channel validation: nucleotide decay patterns (Channel A) combined with stop codon conversion rates (Channel B). To validate these reference-free estimates, we compare AGP against metaDMG, the current standard for ancient DNA damage quantification that aligns reads to reference genomes and counts C→T mismatches at each position.

**Validation dataset**: 31 real ancient metagenomes from the KapK archaeological project, spanning a wide range of preservation states. The samples include 4 low-damage specimens (<2% metaDMG estimate), 5 moderate-damage (8-32%), and 22 high-damage samples (37-56%). metaDMG serves as ground truth: it aligns reads to reference genomes, identifies C→T mismatches, and fits an exponential decay model to estimate the maximum damage rate. AGP's reference-free approach measures T/(T+C) terminal frequencies (Channel A) and tracks CAA→TAA, CAG→TAG, CGA→TGA stop codon conversions (Channel B) to distinguish real deamination from compositional T enrichment.

**Representative samples across damage range**:

| Sample | metaDMG (%) | AGP d_max (%) | Channel B LLR | Decision |
|--------|-------------|---------------|---------------|----------|
| low_001 | 0.48 | 0.00 | -24,707 | No damage (artifact) |
| low_002 | 0.53 | 0.00 | -365 | No damage (artifact) |
| low_005 | 8.37 | 6.27 | +808 | Validated |
| low_009 | 31.47 | 9.06 | +39,631 | Validated |
| 69_34 | 39.50 | 39.00 | +5,332,371 | Validated |
| 119_48 | 44.98 | 44.00 | +3,818,288 | Validated |
| 75_205D | 55.82 | 42.00 | +1,481,034 | Validated |

**Key metrics**:

| Metric | Value |
|--------|-------|
| Pearson correlation (r) | 0.807 |
| Mean bias | +4.4% (AGP higher) |
| Mean absolute error | 8.0% |

The strong correlation (r = 0.807) demonstrates that AGP's reference-free estimates track metaDMG well across the full damage range from <1% to >55%. This is notable because the two methods measure fundamentally different signals: metaDMG counts actual C→T mismatches in aligned reads, while AGP infers damage from terminal nucleotide frequencies without any reference.

**Channel B prevents false positives**: The four low-damage samples (low_001 through low_004) illustrate Channel B's critical role. All four have metaDMG estimates below 1.3% yet show elevated terminal T/(T+C) in Channel A—exactly the pattern that would cause reference-free methods to incorrectly report damage. However, their negative Channel B LLR values (-365 to -24,707) indicate that stop codon conversion rates at termini are actually *lower* than interior baseline, proving the terminal T enrichment is compositional, not deamination. AGP correctly reports d_max = 0 for all four.

**Systematic bias exists but is expected**: AGP's +4.4% average overestimation reflects different estimands. metaDMG only analyzes reads that successfully align to reference genomes—a subset biased toward conserved, potentially better-preserved sequences. AGP analyzes all reads regardless of alignment status. This selection bias cannot be eliminated without reference genomes, but the strong correlation means AGP estimates remain useful for sample authentication and quality control.

**Limitations**: Channel B validation requires sufficient convertible codon coverage (typically >10,000 reads containing CAA, CAG, or CGA at terminal positions). Very short reads (<60 bp) may lack enough interior sequence to establish reliable baseline frequencies. AT-rich organisms can show inverted terminal patterns (lower T at termini than interior), which confounds Channel A; in these cases, Channel B becomes the primary signal.

## Mathematical Foundations

### Damage Model

Post-mortem deamination follows exponential decay from fragment termini:

```
δ(p) = δ_max × e^(-λp) + δ_background
```

Where:
- δ(p) = damage rate at position p from terminus
- δ_max = maximum damage rate (typically 0.1-0.5 for ancient samples)
- λ = decay constant (typically 0.2-0.4)
- p = position in nucleotides from terminus

### Stop Codon Rescue

For a stop codon at position p from 5' terminus, the probability it arose from damage:

```
P(damage | stop) = P(stop | damage) × P(damage) / P(stop)
```

Where:
- P(stop | damage) = probability of convertible precursor (CAA, CAG, CGA)
- P(damage) = δ(p) from the damage model
- P(stop) = observed stop frequency at that position

High P(damage | stop) triggers X-masking in `--fasta-aa-masked` output.

### Frame Scoring

Frame selection combines multiple signals:
- Codon usage bias (weight 0.15)
- Stop codon penalty (weight 0.28)
- Dicodon/hexamer patterns (weight 0.13)
- Amino acid composition (weight 0.10)
- GC3 content (weight 0.05)
- Damage-frame consistency (variable)

In damage-aware mode, stop penalties are weighted by (1 - damage_probability), reducing penalties for likely damage-induced stops.

## Citation

If you use AGP in your research, please cite:

```
[Citation pending publication]
```

## License

MIT License - see LICENSE file for details.
