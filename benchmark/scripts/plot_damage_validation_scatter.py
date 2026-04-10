#!/usr/bin/env python3
"""
DART vs metaDMG damage validation scatter plot.
Reads damage.json + metadmg.tsv directly from each sample subdirectory.
  - Green circles for validated (Channel B+), red squares for rejected (Channel B-)
  - Single regression line through ALL samples
  - r = correlation across all n samples (bottom-right annotation)
  - Title: "DART vs metaDMG Damage Estimates (N samples)"
"""
import sys
import os
import json
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams.update({
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'axes.linewidth': 0.5,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'legend.frameon': True,
    'legend.edgecolor': '0.8',
    'axes.grid': False,
    'font.family': 'sans-serif',
})

SAMPLE_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/vol/cloud/agp/validation_workflow/results/sample_damage"
OUT_PNG = sys.argv[2] if len(sys.argv) > 2 else "damage_validation_scatter.png"
OUT_PDF = OUT_PNG.replace(".png", ".pdf")

# Load all samples that have both damage.json and metadmg.tsv
samples, metadmg_vals, dart_vals, validated = [], [], [], []

for sample_path in sorted(glob.glob(os.path.join(SAMPLE_DIR, "*"))):
    sample = os.path.basename(sample_path)
    json_path = os.path.join(sample_path, "damage.json")
    meta_path = os.path.join(sample_path, "metadmg.tsv")
    if not (os.path.exists(json_path) and os.path.exists(meta_path)):
        continue

    with open(json_path) as f:
        d = json.load(f)
    dmg = d.get("damage", d)

    with open(meta_path) as f:
        next(f)
        meta_dmax = float(f.readline().strip().split("\t")[1])

    samples.append(sample)
    metadmg_vals.append(meta_dmax * 100)   # fraction → %
    dart_vals.append(dmg["d_max"])
    validated.append(bool(dmg.get("damage_validated", dmg.get("validated", True))))

n = len(samples)
print(f"Loaded {n} samples")

metadmg_vals = np.array(metadmg_vals)
dart_vals = np.array(dart_vals)
validated = np.array(validated)

# Pearson r across ALL samples
r_all, _ = stats.pearsonr(metadmg_vals, dart_vals)

# Plot
fig, ax = plt.subplots(figsize=(5.5, 5))

ax_max = max(metadmg_vals.max(), dart_vals.max()) * 1.05
ax_max = max(ax_max, 60)

# Identity line
ax.plot([0, ax_max], [0, ax_max], '--', color='0.5', linewidth=0.8, zorder=1)

# Regression line
slope, intercept, _, _, _ = stats.linregress(metadmg_vals, dart_vals)
x_fit = np.linspace(0, ax_max, 100)
ax.plot(x_fit, slope * x_fit + intercept, '-', color='#2166AC', linewidth=1.2, zorder=2)

# Scatter
mask_v = validated
mask_r = ~validated

ax.scatter(metadmg_vals[mask_v], dart_vals[mask_v],
           c='#2CA02C', s=40, alpha=0.75, edgecolors='none', zorder=3,
           label='Validated (Channel B+)')
ax.scatter(metadmg_vals[mask_r], dart_vals[mask_r],
           c='#D62728', s=40, alpha=0.75, marker='s', edgecolors='none', zorder=3,
           label='Rejected (Channel B\u2212)')

ax.set_xlabel('metaDMG d_max (%)')
ax.set_ylabel('DART d_max (%)')
ax.set_title(f'DART vs metaDMG Damage Estimates ({n} samples)')
ax.set_xlim(-2, ax_max + 2)
ax.set_ylim(-2, ax_max + 2)
ax.set_aspect('equal')
ax.legend(loc='upper left', framealpha=0.9)

ax.text(0.97, 0.03, f'r = {r_all:.2f} (n={n})',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=7, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='0.7', alpha=0.9))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
plt.savefig(OUT_PDF, bbox_inches='tight')
print(f"Saved: {OUT_PNG}")
print(f"r = {r_all:.3f} (n={n})")
print(f"Mean bias (DART - metaDMG): {np.mean(dart_vals - metadmg_vals):.1f}%")

# Print table
print("\n| Sample | metaDMG (%) | DART d_max (%) | Validated |")
print("|--------|-------------|---------------|-----------|")
for s, m, d, v in zip(samples, metadmg_vals, dart_vals, validated):
    print(f"| {s} | {m:.2f} | {d:.2f} | {'Yes' if v else 'No'} |")
