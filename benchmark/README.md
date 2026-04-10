# DART Benchmark

Plotting scripts and pre-computed data for reproducing the benchmark figures in
the [DART wiki](../wiki/Benchmarks.md) and the associated paper (Fernandez-Guerra
et al. 2025).

## Directory Structure

```
benchmark/
├── scripts/
│   ├── plot_roc.py               # Read-level and protein-level damage classification ROC
│   ├── plot_gene_prediction.py   # Gene prediction method comparison (DART vs blastx vs FGS-rs)
│   └── plot_damage_comparison.R  # DART d_max vs metaDMG validation scatter
├── data/
│   ├── roc_benchmark/
│   │   ├── protein_damage_evaluation.tsv
│   │   ├── sample_damage_comparison.tsv
│   │   └── synthetic/            # NOT included — see Download below
│   ├── gene_prediction/
│   │   ├── method_comparison.tsv # Pre-computed metrics from all 18.3M reads
│   │   └── *_hits.tsv            # NOT included — see Download below
│   └── damage_comparison/
│       ├── comparison_table.tsv  # Pre-computed (10-sample subset; full 31-sample on denbi-h-micro)
│       └── *.json                # DART JSON output per sample (10-sample subset)
└── figures/                      # Pre-generated PDF and PNG output
```

## Dependencies

**Python scripts:**
```
pip install numpy pandas matplotlib scikit-learn
```

**R script:**
```r
install.packages(c("data.table", "ggplot2", "jsonlite"))
```

## Running the Scripts

### Damage validation scatter (30 seconds, data included)

```bash
Rscript scripts/plot_damage_comparison.R
# Output: figures/damage_validation_scatter.pdf
```

Note: The included `data/damage_comparison/` contains 10 of 31 samples (local subset).
Full 31-sample data is on denbi-h-micro at `/vol/cloud/agp/paper_figures/damage_comparison/`.
The pre-generated figure in `figures/` uses all 31 samples (r = 0.81).

### Gene prediction comparison (pre-computed, seconds)

```bash
python scripts/plot_gene_prediction.py
# Output: figures/benchmark_comparison.pdf
# Uses pre-computed data/gene_prediction/method_comparison.tsv
```

To recompute from raw hits files (~500 MB each, not distributed):
```bash
# Download from zenodo: [TODO: add DOI]
# Place as data/gene_prediction/agp_hits.tsv, blastx_hits.tsv, fgs_hits.tsv
python scripts/plot_gene_prediction.py
```

### ROC curves — read-level and protein-level (requires raw data)

```bash
python scripts/plot_roc.py
# Output: figures/read_benchmark.pdf, figures/protein_benchmark.pdf
```

Requires per-sample annotated.tsv and read_gt.tsv (~360 MB total):
```bash
# Download from zenodo: [TODO: add DOI]
# Place as data/roc_benchmark/synthetic/<sample>/annotated.tsv
#                                                <sample>/read_gt.tsv
```

Sample names (10 samples from KapK community, synthetic ancient reads):
- AT-rich (d_max > 0): 12-1-25, 12-1-33, 12-1-34, 12-1-35, 12-1-36, 12-1-37
- GC-rich (d_max ≈ 0): 12-1-24, 12-1-27, 12-1-29, 12-1-31

## Pre-generated Figures

The `figures/` directory contains final PDF and PNG for all four benchmark plots.
These are the figures used in the paper and DART wiki.

| File | Description | Metric |
|------|-------------|--------|
| `benchmark_comparison.pdf` | Gene prediction: DART vs blastx vs FGS-rs | Recall, Precision, Avg identity |
| `read_benchmark.pdf` | Read-level damage classification ROC | AUC = 0.844 (AT-rich) |
| `protein_benchmark.pdf` | Protein-level damage classification ROC | AUC = 0.947 (AT-rich) |
| `damage_validation_scatter.pdf` | DART d_max vs metaDMG (31 samples) | r = 0.81 |
