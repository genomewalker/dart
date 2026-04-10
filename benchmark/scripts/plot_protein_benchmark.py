#!/usr/bin/env python3
"""
Protein-level damage classification benchmark.
Ground truth: label each reference protein as ancient/modern based on
actual read labels (:ancient: in query_id) from emi.reads.tsv.
Score: assign_mean_posterior from emi.protein.tsv.
"""
import os
import glob
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.size': 7, 'axes.labelsize': 8, 'axes.titlesize': 9,
    'axes.linewidth': 0.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 6.5, 'legend.frameon': False,
    'lines.linewidth': 1.4, 'axes.grid': False,
    'font.family': 'sans-serif',
})

# Muted dusk palette — warm-cool contrast, publication-quality
C_AT      = '#C76B2A'   # terracotta  — AT-rich (panel A)
C_GC      = '#4A7FA5'   # slate blue  — GC-rich (panel A)
C_ANCIENT = '#C76B2A'   # terracotta  — damaged ground truth (panel B)
C_MODERN  = '#7667A1'   # muted grape — undamaged ground truth (panel B)
C_PREC    = '#953040'   # deep crimson — precision
C_REC     = '#2D5F8E'   # cobalt      — recall
C_F1      = '#2E7055'   # forest      — F1

GC_RICH = {'24', '27', '29', '31'}
AT_RICH = {'25', '33', '34', '35', '36', '37'}

MIN_READS   = 5     # minimum reads per protein for reliable ground truth
PURITY      = 0.80  # minimum fraction from one class (ancient or modern)
MIN_MODERN  = 50    # minimum modern proteins per sample (removes near-pure-ancient samples)
SCORE_COL   = 'assign_mean_posterior'


def fmt_n(n):
    return f'{n/1_000_000:.1f}M' if n >= 1_000_000 else f'{n//1000}k'


def load_sample(sample_dir):
    reads_f   = os.path.join(sample_dir, 'emi.reads.tsv')
    protein_f = os.path.join(sample_dir, 'emi.protein.tsv')
    if not os.path.exists(reads_f) or not os.path.exists(protein_f):
        return [], []

    read_counts = {}
    with open(reads_f) as fh:
        hdr = fh.readline().strip().split('\t')
        qi = hdr.index('query_id')
        ti = hdr.index('target_id')
        for ln in fh:
            p = ln.strip().split('\t')
            if len(p) <= max(qi, ti):
                continue
            target = p[ti]
            is_anc = ':ancient:' in p[qi]
            if target not in read_counts:
                read_counts[target] = [0, 0]
            if is_anc:
                read_counts[target][0] += 1
            else:
                read_counts[target][1] += 1

    dmg, nodmg = [], []
    with open(protein_f) as fh:
        hdr = fh.readline().strip().split('\t')
        pi = hdr.index('protein_id')
        si = hdr.index(SCORE_COL)
        fi = hdr.index('pass_mapping_filter') if 'pass_mapping_filter' in hdr else None
        for ln in fh:
            p = ln.strip().split('\t')
            if len(p) <= max(pi, si):
                continue
            if fi is not None and p[fi] != '1':
                continue
            protein = p[pi]
            try:
                score = float(p[si])
            except ValueError:
                continue
            counts = read_counts.get(protein, [0, 0])
            n_anc, n_mod = counts
            n_total = n_anc + n_mod
            if n_total < MIN_READS:
                continue
            frac_anc = n_anc / n_total
            if frac_anc >= PURITY:
                dmg.append(score)
            elif frac_anc <= (1 - PURITY):
                nodmg.append(score)

    return dmg, nodmg


def roc(d, n):
    if not d or not n:
        return [0, 1], [0, 1], 0.5
    ts = np.linspace(0, 1, 201)
    tpr, fpr = [], []
    for t in ts:
        tp = sum(s >= t for s in d)
        fp = sum(s >= t for s in n)
        tpr.append(tp / len(d))
        fpr.append(fp / len(n))
    auc = sum((fpr[i-1] - fpr[i]) * (tpr[i-1] + tpr[i]) / 2 for i in range(1, len(fpr)))
    return fpr, tpr, auc


def kde_curve(data, x_grid, bw=0.05):
    from scipy.stats import gaussian_kde
    if len(data) < 10 or np.std(data) < 1e-6:
        return np.zeros_like(x_grid)
    return gaussian_kde(data, bw_method=bw)(x_grid)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', required=True)
    parser.add_argument('-o', '--output', default='protein_benchmark.png')
    args = parser.parse_args()

    sample_dirs = sorted(glob.glob(os.path.join(args.data_dir, '*')))
    sample_dirs = [d for d in sample_dirs if os.path.isdir(d)
                   and os.path.exists(os.path.join(d, 'emi.protein.tsv'))]
    print(f"Found {len(sample_dirs)} samples")

    at_d, at_n, gc_d, gc_n, all_d, all_n = [], [], [], [], [], []
    for d in sample_dirs:
        sfx = os.path.basename(d).split('-')[-1]
        dmg, nodmg = load_sample(d)
        if len(nodmg) < MIN_MODERN:
            print(f"  {os.path.basename(d)}: SKIP ({len(nodmg)} modern < {MIN_MODERN})")
            continue
        print(f"  {os.path.basename(d)}: {len(dmg)} ancient, {len(nodmg)} modern proteins")
        all_d.extend(dmg); all_n.extend(nodmg)
        if sfx in AT_RICH:
            at_d.extend(dmg); at_n.extend(nodmg)
        elif sfx in GC_RICH:
            gc_d.extend(dmg); gc_n.extend(nodmg)

    fpr_at,  tpr_at,  auc_at  = roc(at_d, at_n)
    fpr_gc,  tpr_gc,  auc_gc  = roc(gc_d, gc_n)
    fpr_all, tpr_all, auc_all = roc(all_d, all_n)

    # P/R/F1 on AT-rich only. Threshold by Youden's J (TPR − FPR): robust to class imbalance.
    ts = np.linspace(0.05, 0.99, 200)
    pr = []
    for t in ts:
        tp = sum(s >= t for s in at_d)
        fp = sum(s >= t for s in at_n)
        fn = len(at_d) - tp
        tpr = tp / len(at_d) if at_d else 0
        fpr = fp / len(at_n) if at_n else 0
        p  = tp / (tp + fp) if tp + fp else 0
        r  = tpr
        pr.append((t, p, r, 2*p*r/(p+r) if p+r else 0, tpr - fpr))

    ci = np.argmax([x[4] for x in pr])
    opt_t, opt_p, opt_r, opt_f1, _ = pr[ci]

    print(f"\nAUC: Overall={auc_all:.3f}, AT-rich={auc_at:.3f}, GC-rich={auc_gc:.3f}")
    print(f"Optimal (Youden): τ={opt_t:.2f}, P={opt_p:.1%}, R={opt_r:.1%}, F1={opt_f1:.3f}")
    print(f"AT-rich proteins: {len(at_d)} ancient, {len(at_n)} modern")
    print(f"GC-rich proteins: {len(gc_d)} ancient, {len(gc_n)} modern")
    print(f"Score: {SCORE_COL}, min_reads={MIN_READS}, purity={PURITY}, min_modern={MIN_MODERN}")

    x_grid = np.linspace(0, 1, 300)
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    plt.subplots_adjust(hspace=0.38, wspace=0.32)

    # A: score by sample type
    ax = axes[0, 0]
    ax.fill_between(x_grid, kde_curve(at_d + at_n, x_grid), alpha=0.35, color=C_AT,
                    label=f'AT-rich (n={fmt_n(len(at_d+at_n))})')
    ax.fill_between(x_grid, kde_curve(gc_d + gc_n, x_grid), alpha=0.35, color=C_GC,
                    label=f'GC-rich, d_max=0 (n={fmt_n(len(gc_d+gc_n))})')
    ax.plot(x_grid, kde_curve(at_d + at_n, x_grid), color=C_AT, lw=1.2)
    ax.plot(x_grid, kde_curve(gc_d + gc_n, x_grid), color=C_GC, lw=1.2)
    ax.set_xlabel('Mean damage posterior'); ax.set_ylabel('Density')
    ax.set_title('A. Score distribution by sample type')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.99))
    ax.set_xlim(0, 1); ax.set_ylim(bottom=0)

    # B: score by ground truth (AT-rich only)
    ax = axes[0, 1]
    ax.fill_between(x_grid, kde_curve(at_n, x_grid), alpha=0.35, color=C_MODERN,
                    label=f'Undamaged (n={fmt_n(len(at_n))})')
    ax.fill_between(x_grid, kde_curve(at_d, x_grid), alpha=0.35, color=C_ANCIENT,
                    label=f'Damaged (n={fmt_n(len(at_d))})')
    ax.plot(x_grid, kde_curve(at_n, x_grid), color=C_MODERN, lw=1.2)
    ax.plot(x_grid, kde_curve(at_d, x_grid), color=C_ANCIENT, lw=1.2)
    ax.axvline(opt_t, color='#888888', ls='--', lw=0.9, alpha=0.8, label=f'τ={opt_t:.2f}')
    ax.set_xlabel('Mean damage posterior'); ax.set_ylabel('Density')
    ax.set_title('B. Score distribution by ground truth (AT-rich)')
    ax.legend(loc='upper left')
    ax.set_xlim(0, 1); ax.set_ylim(bottom=0)

    # C: ROC curves
    ax = axes[1, 0]
    ax.plot(fpr_at, tpr_at, '-', lw=1.8, color=C_AT, label=f'AT-rich (AUC={auc_at:.2f})')
    ax.plot(fpr_gc, tpr_gc, '-', lw=1.8, color=C_GC, label=f'GC-rich (AUC={auc_gc:.2f})')
    ax.plot([0, 1], [0, 1], 'k--', lw=0.6, alpha=0.3)
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
    ax.set_title('C. ROC curves')
    ax.legend(loc='lower right')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # D: P/R/F1 (AT-rich, Youden's J)
    # Empty space: threshold 0.40-0.63, score 0-0.92 (all curves at 0.95+ there)
    ax = axes[1, 1]
    ax.plot([x[0] for x in pr], [x[1] for x in pr], '-', lw=1.4, color=C_PREC, label='Precision')
    ax.plot([x[0] for x in pr], [x[2] for x in pr], '-', lw=1.4, color=C_REC,  label='Recall')
    ax.plot([x[0] for x in pr], [x[3] for x in pr], '-', lw=1.4, color=C_F1,   label='F1')
    ax.axvline(opt_t, color='#888888', ls='--', lw=0.9, alpha=0.8)
    ax.text(0.04, 0.52,
            f'\u03c4* = {opt_t:.2f}\nP  = {opt_p:.1%}\nR  = {opt_r:.1%}\nF1 = {opt_f1:.3f}',
            transform=ax.transAxes, fontsize=6.5, va='top', ha='left', color='#444444',
            linespacing=1.6)
    ax.set_xlabel('Threshold'); ax.set_ylabel('Score')
    ax.set_title('D. Precision, recall, F1 by threshold (AT-rich)')
    ax.legend(loc='lower left')
    ax.set_xlim(0.4, 1.0); ax.set_ylim(0, 1.02)

    plt.savefig(args.output, dpi=150, bbox_inches='tight')
    pdf_out = args.output.replace('.png', '.pdf')
    plt.savefig(pdf_out, bbox_inches='tight')
    print(f"Saved to: {args.output}")


if __name__ == '__main__':
    main()
