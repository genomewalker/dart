# DART Benchmark

Plotting scripts and pre-computed data for reproducing the benchmark figures in
the [DART wiki](https://github.com/genomewalker/dart/wiki/Benchmarks) and the associated paper (Fernandez-Guerra
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
│   │   └── synthetic/            # not included
│   ├── gene_prediction/
│   │   ├── method_comparison.tsv # Pre-computed metrics from all 18.3M reads
│   │   └── *_hits.tsv            # not included
│   └── damage_comparison/
│       ├── comparison_table.tsv  # Pre-computed (10 of the 31 samples)
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

`data/damage_comparison/` holds 10 of the 31 samples, so the plot you get differs from the
figure in `figures/`, which uses all 31.

### Gene prediction comparison (pre-computed, seconds)

```bash
python scripts/plot_gene_prediction.py
# Output: figures/benchmark_comparison.pdf
# Uses pre-computed data/gene_prediction/method_comparison.tsv
```

Recomputing the metrics needs the raw hits files (about 500 MB each), which are not
distributed yet. Place them as `data/gene_prediction/agp_hits.tsv`, `blastx_hits.tsv` and
`fgs_hits.tsv` and rerun the script.

### ROC curves — read-level and protein-level (requires raw data)

```bash
python scripts/plot_roc.py
# Output: figures/read_benchmark.pdf, figures/protein_benchmark.pdf
```

This needs the per-sample `annotated.tsv` and `read_gt.tsv` (about 360 MB in total), which are
not distributed yet. Place them as `data/roc_benchmark/synthetic/<sample>/annotated.tsv` and
`<sample>/read_gt.tsv`.

The 10 samples are synthetic ancient reads from the KapK community:
- AT-rich (d_max > 0): 12-1-25, 12-1-33, 12-1-34, 12-1-35, 12-1-36, 12-1-37
- GC-rich (d_max ≈ 0): 12-1-24, 12-1-27, 12-1-29, 12-1-31

## Pre-generated Figures

`figures/` has the PDF and PNG of all four plots, as used in the paper and the DART wiki.

| File | Description | Metric |
|------|-------------|--------|
| `benchmark_comparison.pdf` | Gene prediction: DART vs blastx vs FGS-rs | Recall, Precision, Avg identity |
| `read_benchmark.pdf` | Read-level damage classification ROC | AUC = 0.844 (AT-rich) |
| `protein_benchmark.pdf` | Protein-level damage classification ROC | AUC = 0.947 (AT-rich) |
| `damage_validation_scatter.pdf` | DART d_max vs metaDMG (31 samples) | r = 0.81 |
