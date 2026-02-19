# AGP: Damage-Aware Gene prediction from ancienT mEtagenomics Reads

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

AGP translates ancient DNA (aDNA) reads into proteins for functional annotation. It predicts protein-coding genes from short, degraded metagenomic sequences while accounting for post-mortem DNA damage that corrupts standard translation.

## Overview

<p align="center">
<img src="docs/agp_pipeline.svg" width="700" alt="AGP pipeline: three-row overview">
</p>

The pipeline runs in three stages. **Pass 1** scans all reads to measure the sample's damage level. **Pass 2** translates reads into proteins using that damage estimate to rescue genes that would otherwise be missed. **Functional annotation** searches those proteins against databases and scores each hit for authentic ancient damage.

Post-mortem deamination converts cytosines to uracil (read as T), producing C→T gradients at read termini. This breaks standard gene predictors by creating false stop codons. AGP models the damage and corrects for it during translation; a second independent signal (stop codon conversion rates) validates that the enrichment is real damage rather than genomic T-richness. For full methods and model derivations, see the [wiki](https://github.com/genomewalker/agp/wiki/Methods-and-Model).

---

## Installation

Requires C++20 compiler (GCC 10+ or Clang 12+), CMake 3.18+, zlib, and zstd.

The easiest way to get all build dependencies is via conda:

```bash
conda env create -f environment.yml
conda activate agp
```

Then build:

```bash
git clone https://github.com/genomewalker/agp.git
cd agp
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

---

## Quick start

### 5-minute tutorial run (recommended)

AGP ships a reproducible tutorial in `examples/tutorial/` with 50,000 synthetic KapK ancient reads and a matched reference protein subset:

```bash
cd examples/tutorial
AGP=../../build/agp DB=./ref_proteins.faa.gz bash run_tutorial.sh
```

The four pipeline stages:

```bash
# 1) Predict genes + AGD damage index
agp predict -i reads.fq.gz -o out/predictions.gff \
  --fasta-aa out/proteins.faa --damage-index out/predictions.agd --adaptive -t 8
```
```
Ancient Gene Predictor v0.4.0
Input: reads.fq.gz  |  Domain: gtdb  |  Threads: 8

[Pass 1] Scanning for damage patterns...
  Reads: 50000 | 5' damage: 28.6% | 3' damage: 26.3% | 4.1 s

[Pass 2] Gene prediction...
  Damage index: 192449 records
  Sequences: 50000 | Genes: 192449 | 2.0 s (25062 seq/s)

Done. Total runtime: 6.1 s
```

```bash
# 2) Search against reference proteins (VTML20 optimised for aDNA)
mmseqs easy-search out/proteins.faa ref_proteins.faa.gz out/hits.tsv out/tmp \
  --search-type 1 --min-length 12 -e 10.0 --min-seq-id 0.86 -c 0.65 --cov-mode 2 \
  --sub-mat VTML20.out --seed-sub-mat VTML20.out \
  -s 2 -k 6 --spaced-kmer-pattern 11011101 \
  --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln" \
  --threads 8
```
```
[=================================================================] 185.1K  ~10 s
Time for merging to hits.tsv: 0h 0m 0s 267ms
```

```bash
# 3) Build EMI index
agp hits2emi -i out/hits.tsv -o out/hits.emi --damage-index out/predictions.agd -t 8
```
```
Lines: 46628 (scan: 3 ms), building dictionary...
  Embedded per-read damage scores  |  d_max=25.1%  lambda=0.424
  Dictionary: 133 ms (22367 reads, 13681 refs, 8 threads)
Parsed: 46628 rows in 249 ms  |  Written: out/hits.emi  |  Total: 298 ms
```

```bash
# 4) Damage annotation
agp damage-annotate --emi out/hits.emi --damage-index out/predictions.agd \
  -o out/annotated.tsv --gene-summary out/gene_summary.tsv --auto-calibrate-spurious -t 8
```
```
Loaded damage index: 192449 records  |  d_max: 25.1%  |  lambda: 0.424
EMI: 46628 alignments, 22367 reads, 13681 refs
EM converged in 8 iterations
  Unique mappers: 13011  |  Multi-mappers: 9356  |  Reassigned: 8.1%
Reads with damage-consistent sites: 5962  |  Total sites: 6722
Auto-calibration (n=397 genes): min_positional_score=0.786
Gene summary: 591/7645 genes after filtering
Summary: assigned_reads=21819, mismatch_reads=8624, reads_with_damage_sites=5962
Runtime: 186 ms
```

Use `--fasta-aa-masked` instead of `--fasta-aa` to replace damage-induced stop codons with X so MMseqs2 can align across them.

### Sample damage assessment only

```bash
agp sample-damage reads.fq.gz
```

Output (JSON):
```json
{
  "version": "0.4.0",
  "input": "reads.fq.gz",
  "domain": "gtdb",
  "sequences": 50000,
  "damage": {
    "level": "HIGH DAMAGE",
    "d_max": 16.3,
    "lambda_5prime": 0.480,
    "lambda_3prime": 0.480,
    "damage_validated": true,
    "library_type": "double-stranded"
  }
}
```

`d_max` values are reported as percentages (0–100).

---

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
  --domain NAME          Taxonomic domain: gtdb, fungi, plant, protozoa,
                         invertebrate, viral (default: gtdb)
  --orf-min-aa N         Minimum protein length (default: 10)
  -t, --threads N        Thread count (default: auto)
```

### `agp sample-damage`

Estimates sample-wide damage rate without gene prediction.

```
Usage: agp sample-damage <input.fq> [options]

Output fields:
  damage.d_max               Maximum damage rate at read termini (percent, 0-100)
  damage.lambda_5prime / lambda_3prime
                             Exponential decay constants
  damage.damage_validated    True if damage confirmed by stop codon signal
  damage.library_type        single-stranded or double-stranded
```

### `agp hits2emi`

Converts MMseqs2 hit TSV into AGP EMI (columnar) with alignment strings required by `damage-annotate`.

```
Usage: agp hits2emi -i <hits.tsv[.gz]> -o <hits.emi> [options]

Required:
  -i, --tsv FILE         MMseqs2 16-column hits TSV (plain or gz)
  -o, --output FILE      Output EMI file

Options:
  -t, --threads N        Thread count (default: all available)
  -m, --memory SIZE      Memory budget with K/M/G/T suffix (default: 10% available RAM)
  --streaming            Force streaming mode
  --damage-index FILE    Optional AGD; embeds per-read damage probability into EMI
  --d-max FLOAT          Store sample d_max in EMI metadata
  --lambda FLOAT         Store decay lambda in EMI metadata
```

### `agp damage-annotate`

Annotates database hits with per-protein damage evidence. Compares observed proteins to reference proteins to identify damage-consistent amino acid substitutions (R→W, H→Y, Q→*, etc.).

```
Usage: agp damage-annotate --emi <hits.emi> -o <output.tsv> [options]

Required:
  --emi FILE             EMI index from `agp hits2emi` (with qaln/taln)

Optional:
  -t, --threads INT      OpenMP threads for EMI scans
  --protein-summary FILE Per-protein aggregated statistics
  --protein-filtered FILE Per-protein subset that passes mapping-pattern filters
  --combined-output FILE Per-target combined protein+gene summary
  --blast8-unique FILE   BLAST8-like export with unique mappers only
  --analysis-prefix STR  Output prefix: .reads.tsv, .protein.tsv, .blast8.unique.tsv,
                         and .categories.tsv (if --map is provided)
  --gene-summary FILE    Per-protein/gene mapping summary (coverage/diversity metrics)
  --map FILE             Gene-to-group mapping TSV (gene_id<TAB>group)
  --functional-summary FILE  Per-group stats TSV
  --anvio-ko FILE        Gene-level group abundance for anvi-estimate-metabolism
  --annotation-source STR  Source label for anvi output (default: AGP)
  --preset STR           Apply filter bundle: loose or strict
  --min-reads INT        Minimum supporting reads (default: 3)
  --min-breadth FLOAT    Minimum alignment breadth (default: 0.10)
  --min-depth FLOAT      Minimum mean depth (default: 0.5; set to 0 for ancient DNA)
  --prior-ancient FLOAT  Prior P(ancient) for Bayesian scorer (default: 0.10)
  --auto-prior-ancient   Auto-calibrate prior from mean p_read (requires AGD index)
  --min-damage-sites INT Min susceptible positions for site evidence (default: auto)
  --threshold FLOAT      Classification threshold for is_damaged (default: 0.7)
  --aln-min-identity F   Min alignment identity (default: 0)
  --aln-min-bits F       Min alignment bit score (default: 0)
  --aln-max-evalue F     Max alignment e-value (default: 1e10)
  --no-em                Disable EM reassignment (enabled by default)
  --em-iters INT         Max EM iterations (default: 100)
  --em-lambda FLOAT      EM score temperature (default: 3.0)
  --em-tol FLOAT         EM convergence tolerance (default: 1e-4)
  --em-min-prob FLOAT    Min EM responsibility to keep a read (default: 1e-6)
```

Output columns: `p_read`, `ct_sites`, `ga_sites`, `posterior`, `damage_class`, `is_damaged`, `gamma`, `em_keep`. See [Output Formats](https://github.com/genomewalker/agp/wiki/Output-Formats) for full schema.

#### Filtering presets

```bash
agp damage-annotate --emi hits.emi -o out.tsv --preset loose   # maximize recall
agp damage-annotate --emi hits.emi -o out.tsv --preset strict  # maximize precision
```

| Parameter | default | `--preset loose` | `--preset strict` |
|-----------|---------|-----------------|------------------|
| `--min-reads` | 3 | 1 | 5 |
| `--min-breadth` | 0.10 | 0 | 0.20 |
| `--min-depth` | 0.5 | 0 | 1.0 |
| `--min-positional-score` | 0 | 0 | auto-calibrated |
| `--min-terminal-ratio` | 0 | 0 | auto-calibrated |
| `--min-damage-sites` | auto | 1 | 5 |

Coverage filters (`--min-positional-score`, `--min-terminal-ratio`) and Bayesian scoring parameters are documented in the [wiki](https://github.com/genomewalker/agp/wiki/Methods-and-Model#coverage-quality-filters).

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
Usage: agp validate <fastq.gz> <aa-damage.tsv.gz> <predictions.gff> [proteins.faa] [--csv output.csv]
```

---

## Output formats

### GFF3

```
read_001  AncientGenePredictor  gene  1  150  0.85  +  .  ID=gene1;frame=0;conf=0.88;damage_signal=0.41;damage_pct=15.3;p_damaged=0.72
```

| Attribute | Description |
|-----------|-------------|
| `p_damaged` | Per-read damage posterior from prediction |
| `damage_pct` | Estimated damage percentage |
| `conf` | Prediction confidence from length and coding score |
| `frame` | Reading frame used for the prediction |

### Damage Index (.agd)

Binary format for O(1) per-read damage lookup. Contains sequence length, frame, strand, damage probability, and terminal codon information. Required by `damage-annotate`.

### EMI (.emi) and Annotated TSV

See [Output Formats](https://github.com/genomewalker/agp/wiki/Output-Formats) on the wiki for full column schemas.

---

## Benchmark

### Gene prediction

Benchmarked on 18.3 million synthetic ancient DNA reads from the KapK community (10 samples). A prediction is correct if it matches any reference protein at ≥90% sequence identity.

| Method | Recall | Precision | Avg Identity |
|--------|--------|-----------|--------------|
| AGP | 67.6% | **97.2%** | **96.2%** |
| MMseqs2 blastx [Steinegger & Söding 2017] | **68.5%** | 96.3% | 94.8% |
| FGS-rs [Van der Jeugt et al. 2022] | 19.1% | 95.3% | 94.2% |

AGP and BLASTX have equivalent functional recall (~68%), both far exceeding FGS-rs (19%), which treats damage-induced stop codons as real stops. AGP produces 1.4% higher identity hits than BLASTX and runs ~1.5× faster (~35,000 reads/s, 8 threads).

<p align="center">
<img src="docs/benchmark_comparison.png" width="700" alt="Method comparison: AGP, BLASTX, FGS-rs">
</p>

### Damage classification

Protein-level damage classification benchmarked on the same 10 KapK synthetic samples. Each protein is labelled ancient or modern based on the origin of its supporting reads (≥80% purity, ≥5 reads). Score: `assign_mean_posterior` from `damage-annotate`.

| Group | AUC-ROC |
|-------|---------|
| Overall | **0.881** |
| AT-rich samples | **0.957** |
| GC-rich samples | **0.943** |

At the Youden-optimal threshold τ = 0.70: **Precision 99.7%, Recall 89.5%, F1 = 0.943**.

<p align="center">
<img src="docs/protein_benchmark.png" width="700" alt="Protein-level damage classification: ROC curves and precision/recall">
</p>

Read-level damage classification (per-read `posterior` score):

<p align="center">
<img src="docs/read_benchmark.png" width="700" alt="Read-level damage classification benchmark">
</p>

Full benchmark details: [wiki/Benchmarks](https://github.com/genomewalker/agp/wiki/Benchmarks).

---

## Citation

If you use AGP in your research, please cite:

> **Two-million-year-old microbial communities from the Kap København Formation in North Greenland**
> Fernandez-Guerra A, Wörmer L, Borrel G, et al.
> *bioRxiv* (2025)
> DOI: [10.1101/2023.06.10.544454](https://doi.org/10.1101/2023.06.10.544454)

```bibtex
@article{fernandezguerra2025,
  title={Two-million-year-old microbial communities from the Kap København Formation in North Greenland},
  author={Fernandez-Guerra, Antonio and W{\"o}rmer, Lars and Borrel, Guillaume and Delmont, Tom O and Elberling, Bo and Elvert, Marcus and Eren, A Murat and Gribaldo, Simonetta and Henriksen, Rasmus Amund and Hinrichs, Kai-Uwe and Jochheim, Annika and Korneliussen, Thorfinn S and Krupovic, Mart and Larsen, Nicolaj K and Perez-Laso, Rafael and Pedersen, Mikkel Winther and Pedersen, Vivi K and Ruter, Anthony H and Sand, Karina K and Sikora, Martin and Steinegger, Martin and Veseli, Iva and Wang, Yucheng and Zhao, Lei and {\v{Z}}ure, Marina and Kj{\ae}r, Kurt H and Willerslev, Eske},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2023.06.10.544454}
}
```

## References

- **Briggs et al. 2007.** Patterns of damage in genomic DNA sequences from a Neandertal. *PNAS* 104(37):14616–14621. [doi:10.1073/pnas.0704665104](https://doi.org/10.1073/pnas.0704665104)
- **Dempster et al. 1977.** Maximum likelihood from incomplete data via the EM algorithm. *J Royal Statistical Society B* 39(1):1–22.
- **Jónsson et al. 2013.** mapDamage2.0. *Bioinformatics* 29(13):1682–1684. [doi:10.1093/bioinformatics/btt193](https://doi.org/10.1093/bioinformatics/btt193)
- **Lindahl 1993.** Instability and decay of the primary structure of DNA. *Nature* 362:709–715. [doi:10.1038/362709a0](https://doi.org/10.1038/362709a0)
- **Michelsen et al. 2022.** metaDMG. *bioRxiv*. [doi:10.1101/2022.12.06.519264](https://doi.org/10.1101/2022.12.06.519264)
- **Müller & Vingron 2000.** Modeling amino acid replacement. *J Computational Biology* 7(6):761–776.
- **Steinegger & Söding 2017.** MMseqs2. *Nature Biotechnology* 35(11):1026–1028. [doi:10.1038/nbt.3988](https://doi.org/10.1038/nbt.3988)
- **Van der Jeugt et al. 2022.** FragGeneScanRs. *BMC Bioinformatics* 23(1):198. [doi:10.1186/s12859-022-04736-5](https://doi.org/10.1186/s12859-022-04736-5)
- **Varadhan & Roland 2008.** SQUAREM. *Scandinavian J Statistics* 35(2):335–353.

---

## License

MIT License; see LICENSE file for details.
