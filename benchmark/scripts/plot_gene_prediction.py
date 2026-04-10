#!/usr/bin/env python3
"""
Gene prediction method comparison: DART vs MMseqs2 blastx vs FGS-rs.
Reproduces benchmark_comparison.png from the DART wiki.

Input (large files, not included in repo — see README for download):
  data/gene_prediction/agp_hits.tsv    — DART adaptive hits (mmseqs2 tsv format)
  data/gene_prediction/blastx_hits.tsv — MMseqs2 blastx hits
  data/gene_prediction/fgs_hits.tsv    — FGS-rs translated + searched hits

  Column format: query_id  target_id  fident  alnlen  read_len  bits  evalue
  Query ID encodes ground truth: <sample>___<gt_protein>---<read_id>:ancient|modern:...

Usage:
  python scripts/plot_gene_prediction.py

Output:
  figures/benchmark_comparison.pdf
  data/gene_prediction/method_comparison.tsv  (pre-computed metrics)
"""

import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

METHODS = {
    'agp':    'DART',
    'blastx': 'MMseqs2 blastx',
    'fgs':    'FGS-rs',
}
COLORS = {
    'DART':           '#2ecc71',
    'MMseqs2 blastx': '#e74c3c',
    'FGS-rs':         '#3498db',
}

BASE_DIR = os.path.join(os.path.dirname(__file__), "..")
DATA_DIR = os.path.join(BASE_DIR, "data", "gene_prediction")
OUT_DIR  = os.path.join(BASE_DIR, "figures")
MIN_FIDENT = 0.90   # ≥90% identity threshold for "correct" hit


def parse_ground_truth(query_id):
    """Extract ground truth protein from query ID.

    Query format: <sample>___<gt_protein>---<read_info>:ancient|modern:...
    Returns (gt_protein, is_ancient).
    """
    m = re.match(r'.+___(.+)---(.+):(ancient|modern):', query_id)
    if not m:
        return None, None
    return m.group(1), m.group(3) == 'ancient'


def evaluate_hits(hits_file, method):
    """Compute recall, precision for a hits TSV file."""
    if not os.path.exists(hits_file):
        print(f"  Warning: {hits_file} not found — skipping {method}", file=sys.stderr)
        return None

    best = {}   # query_id -> (fident, target_id)
    total_reads = set()

    with open(hits_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            query_id, target_id = parts[0], parts[1]
            try:
                fident = float(parts[2])
            except ValueError:
                continue
            total_reads.add(query_id)
            if query_id not in best or fident > best[query_id][0]:
                best[query_id] = (fident, target_id)

    correct_hits = 0
    any_hits     = len(best)
    for query_id, (fident, target_id) in best.items():
        gt_protein, _ = parse_ground_truth(query_id)
        if gt_protein is None:
            continue
        # Strip protein suffix (_N) from target for comparison
        target_base = re.sub(r'_\d+$', '', target_id)
        gt_base     = re.sub(r'_\d+$', '', gt_protein)
        if fident >= MIN_FIDENT and target_base == gt_base:
            correct_hits += 1

    n_reads   = len(total_reads)
    recall    = correct_hits / n_reads    if n_reads    > 0 else 0.0
    precision = correct_hits / any_hits   if any_hits   > 0 else 0.0

    return {
        'method':        method,
        'label':         METHODS[method],
        'n_reads':       n_reads,
        'correct_hits':  correct_hits,
        'any_hits':      any_hits,
        'recall':        recall,
        'precision':     precision,
    }


def load_precomputed():
    """Load pre-computed comparison table if raw data not available."""
    tsv = os.path.join(DATA_DIR, "method_comparison.tsv")
    if os.path.exists(tsv):
        return pd.read_csv(tsv, sep='\t')
    return None


def plot(metrics_df, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(8, 5))

    labels = list(metrics_df['label'])
    x      = np.arange(len(labels))
    width  = 0.35

    colors_r = [COLORS[l] for l in labels]
    colors_p = [c + 'aa' for c in colors_r]   # slightly transparent

    bars1 = ax.bar(x - width/2, metrics_df['recall'],    width,
                   label='Recall',    color=colors_r, edgecolor='black', linewidth=1.2)
    bars2 = ax.bar(x + width/2, metrics_df['precision'], width,
                   label='Precision', color=colors_r, edgecolor='black', linewidth=1.2,
                   alpha=0.65)

    for bar, val in zip(bars1, metrics_df['recall']):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.015,
                f'{val:.1%}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    for bar, val in zip(bars2, metrics_df['precision']):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.015,
                f'{val:.1%}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_ylabel('Score', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 1.15)
    ax.set_title(
        'Gene prediction method comparison\n'
        '(≥90% identity match to ground-truth proteins, 18.3M synthetic ancient reads)',
        fontsize=12
    )
    ax.legend(loc='upper right', fontsize=10)

    n_reads = metrics_df.loc[metrics_df['method'] == 'agp', 'n_reads'].values
    if len(n_reads):
        ax.text(0.98, 0.02, f'n = {n_reads[0]:,} reads',
                transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")


def main():
    # Try to compute from raw data; fall back to pre-computed table
    metrics = []
    for key in ['agp', 'blastx', 'fgs']:
        hits_file = os.path.join(DATA_DIR, f"{key}_hits.tsv")
        result = evaluate_hits(hits_file, key)
        if result:
            metrics.append(result)

    if metrics:
        df = pd.DataFrame(metrics)
        os.makedirs(DATA_DIR, exist_ok=True)
        df.to_csv(os.path.join(DATA_DIR, "method_comparison.tsv"), sep='\t', index=False)
    else:
        df = load_precomputed()
        if df is None:
            print("No data found. Download agp_hits.tsv, blastx_hits.tsv, fgs_hits.tsv "
                  "to data/gene_prediction/ — see README.", file=sys.stderr)
            sys.exit(1)
        print("Using pre-computed method_comparison.tsv")

    print(df[['label', 'recall', 'precision']].to_string(index=False))
    plot(df, os.path.join(OUT_DIR, "benchmark_comparison.pdf"))


if __name__ == "__main__":
    main()
