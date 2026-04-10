#!/usr/bin/env python3
"""
ROC curves for read-level and protein-level damage classification.
Reproduces read_benchmark.png and protein_benchmark.png from the DART wiki.

Input:
  data/roc_benchmark/synthetic/<sample>/annotated.tsv  — DART output (per-read scores)
  data/roc_benchmark/synthetic/<sample>/read_gt.tsv    — ground truth labels

Usage:
  python scripts/plot_roc.py

Output:
  figures/read_benchmark.pdf
  figures/protein_benchmark.pdf
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

SAMPLES = [
    "119_B3_116_L0_KapK-12-1-24",
    "119_B3_116_L0_KapK-12-1-25",
    "119_B3_116_L0_KapK-12-1-27",
    "119_B3_116_L0_KapK-12-1-29",
    "69_B2_97_L0_KapK-12-1-31",
    "69_B2_97_L0_KapK-12-1-33",
    "69_B2_100_L0_KapK-12-1-34",
    "69_B2_100_L0_KapK-12-1-35",
    "69_B2_100_L0_KapK-12-1-36",
    "69_B2_103_L0_KapK-12-1-37",
]

# High-GC samples (AGP detects d_max=0 due to weak terminal signal)
HIGH_GC = {
    "119_B3_116_L0_KapK-12-1-24",
    "119_B3_116_L0_KapK-12-1-27",
    "119_B3_116_L0_KapK-12-1-29",
    "69_B2_97_L0_KapK-12-1-31",
}

BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "roc_benchmark", "synthetic")
OUT_DIR  = os.path.join(os.path.dirname(__file__), "..", "figures")


def load_read_data(sample_dir):
    gt = {}
    with open(os.path.join(sample_dir, "read_gt.tsv")) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gt[parts[0]] = int(parts[1])

    preds = {}
    with open(os.path.join(sample_dir, "annotated.tsv")) as f:
        header = f.readline().strip().split('\t')
        query_idx  = header.index('query_id')
        score_idx  = header.index('combined_score')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(query_idx, score_idx):
                continue
            read_name = re.match(r'(.+)_[+-]_[012]$', parts[query_idx])
            read_name = read_name.group(1) if read_name else parts[query_idx]
            try:
                score = float(parts[score_idx])
            except ValueError:
                continue
            if read_name not in preds or score > preds[read_name]:
                preds[read_name] = score

    y_true, y_score = [], []
    for name, score in preds.items():
        if name in gt:
            y_true.append(gt[name])
            y_score.append(score)
    return np.array(y_true), np.array(y_score)


def load_protein_data(sample_dir):
    """Aggregate per-read scores to protein level (mean of top-scoring read per protein)."""
    gt_read = {}
    with open(os.path.join(sample_dir, "read_gt.tsv")) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gt_read[parts[0]] = int(parts[1])

    protein_scores = {}
    protein_gt     = {}
    with open(os.path.join(sample_dir, "annotated.tsv")) as f:
        header    = f.readline().strip().split('\t')
        query_idx = header.index('query_id')
        score_idx = header.index('p_damaged')   # protein-level posterior
        target_idx = header.index('target_id')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(query_idx, score_idx, target_idx):
                continue
            read_match = re.match(r'(.+)_[+-]_[012]$', parts[query_idx])
            read_name  = read_match.group(1) if read_match else parts[query_idx]
            protein    = parts[target_idx]
            try:
                score = float(parts[score_idx])
            except ValueError:
                continue
            if protein not in protein_scores or score > protein_scores[protein]:
                protein_scores[protein] = score
            if read_name in gt_read and protein not in protein_gt:
                protein_gt[protein] = gt_read[read_name]

    y_true, y_score = [], []
    for prot, score in protein_scores.items():
        if prot in protein_gt:
            y_true.append(protein_gt[prot])
            y_score.append(score)
    return np.array(y_true), np.array(y_score)


def plot_roc_panel(results_at, results_gc, title, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    at_y = np.concatenate([r['y_true']  for r in results_at])
    at_s = np.concatenate([r['y_score'] for r in results_at])
    gc_y = np.concatenate([r['y_true']  for r in results_gc])
    gc_s = np.concatenate([r['y_score'] for r in results_gc])
    all_y = np.concatenate([at_y, gc_y])
    all_s = np.concatenate([at_s, gc_s])

    at_auc   = roc_auc_score(at_y, at_s)
    gc_auc   = roc_auc_score(gc_y, gc_s)
    mean_auc = np.mean([r['auc'] for r in results_at + results_gc])

    # Panel 1: score distribution
    ax = axes[0]
    bins = np.linspace(0, 1.2, 61)
    ax.hist(all_s[all_y == 1], bins=bins, alpha=0.7, color='#E63946', density=True,
            label=f'Damaged (n={int((all_y==1).sum()):,})')
    ax.hist(all_s[all_y == 0], bins=bins, alpha=0.7, color='#457B9D', density=True,
            label=f'Undamaged (n={int((all_y==0).sum()):,})')
    ax.axvline(0.7, color='green', linestyle='--', linewidth=1.5, label='Threshold (0.7)')
    ax.set_xlabel('Score', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title(f'Score Distribution\n(All {len(SAMPLES)} samples)', fontsize=12)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 1.2)

    # Panel 2: ROC curves
    ax = axes[1]
    for r in results_at:
        fpr, tpr, _ = roc_curve(r['y_true'], r['y_score'])
        ax.plot(fpr, tpr, alpha=0.3, linewidth=1, color='#2A9D8F')
    for r in results_gc:
        fpr, tpr, _ = roc_curve(r['y_true'], r['y_score'])
        ax.plot(fpr, tpr, alpha=0.3, linewidth=1, color='#E76F51')
    fpr_at, tpr_at, _ = roc_curve(at_y, at_s)
    fpr_gc, tpr_gc, _ = roc_curve(gc_y, gc_s)
    ax.plot(fpr_at, tpr_at, linewidth=2.5, color='#2A9D8F',
            label=f'AT-rich (n={len(results_at)}): AUC={at_auc:.3f}')
    ax.plot(fpr_gc, tpr_gc, linewidth=2.5, color='#E76F51',
            label=f'GC-rich (n={len(results_gc)}): AUC={gc_auc:.3f}')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random')
    ax.set_xlabel('False Positive Rate', fontsize=11)
    ax.set_ylabel('True Positive Rate', fontsize=11)
    ax.set_title('ROC Curves by Sample Type', fontsize=12)
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal')

    # Panel 3: per-sample AUC bar chart
    ax = axes[2]
    all_res = sorted(results_at + results_gc, key=lambda x: x['auc'])
    samples = [r['sample'] for r in all_res]
    aucs    = [r['auc']    for r in all_res]
    colors  = ['#E76F51' if r['is_gc'] else '#2A9D8F' for r in all_res]
    bars = ax.barh(samples, aucs, color=colors, edgecolor='white')
    for bar, auc in zip(bars, aucs):
        ax.text(auc + 0.003, bar.get_y() + bar.get_height()/2,
                f'{auc:.3f}', va='center', fontsize=9)
    ax.axvline(mean_auc, color='#1D3557', linestyle='--', linewidth=1.5,
               label=f'Mean: {mean_auc:.3f}')
    ax.set_xlabel('AUC-ROC', fontsize=11)
    ax.set_title('Per-Sample AUC\n(Green=AT-rich, Orange=GC-rich)', fontsize=12)
    ax.legend(loc='lower right', fontsize=9)

    plt.suptitle(title, fontsize=13, y=1.02)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")
    print(f"  AT-rich AUC={at_auc:.3f}, GC-rich AUC={gc_auc:.3f}, mean={mean_auc:.3f}")


def main():
    read_at, read_gc = [], []
    prot_at, prot_gc = [], []

    for sample in SAMPLES:
        sample_dir = os.path.join(BASE_DIR, sample)
        if not os.path.exists(os.path.join(sample_dir, "read_gt.tsv")):
            print(f"Skipping {sample} (data not found)")
            continue

        y_true_r, y_score_r = load_read_data(sample_dir)
        y_true_p, y_score_p = load_protein_data(sample_dir)
        label = sample.split("KapK-")[1]
        is_gc = sample in HIGH_GC

        if len(y_true_r) >= 100:
            rec = {'sample': label, 'n_reads': len(y_true_r),
                   'auc': roc_auc_score(y_true_r, y_score_r),
                   'y_true': y_true_r, 'y_score': y_score_r, 'is_gc': is_gc}
            (read_gc if is_gc else read_at).append(rec)

        if len(y_true_p) >= 50:
            rec = {'sample': label, 'n_reads': len(y_true_p),
                   'auc': roc_auc_score(y_true_p, y_score_p),
                   'y_true': y_true_p, 'y_score': y_score_p, 'is_gc': is_gc}
            (prot_gc if is_gc else prot_at).append(rec)

    plot_roc_panel(read_at, read_gc,
                   "Read-level damage classification (DART damage-annotate)",
                   os.path.join(OUT_DIR, "read_benchmark.pdf"))

    plot_roc_panel(prot_at, prot_gc,
                   "Protein-level damage classification (DART damage-annotate)",
                   os.path.join(OUT_DIR, "protein_benchmark.pdf"))


if __name__ == "__main__":
    main()
