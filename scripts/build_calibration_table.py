#!/usr/bin/env python3
"""
Build calibration table from aggregated protein hit results.

Reads GFF predictions and MMseqs2 hits from all samples,
computes FDR using target-decoy, and generates calibration header.

Usage:
    python build_calibration_table.py \
        --results-dir /projects/caeg/scratch/kbd606/agp/calibration \
        --output-header include/dart/score_calibration_protein.hpp
"""

import argparse
import gzip
import numpy as np
from pathlib import Path
from collections import defaultdict
import sys


def parse_gff_scores(gff_path: Path) -> dict:
    """Extract scores and metadata from AGP GFF output."""
    scores = {}

    opener = gzip.open if str(gff_path).endswith('.gz') else open
    with opener(gff_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            attrs = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attrs[key] = val

            seq_id = parts[0]
            strand = parts[6]
            frame = int(attrs.get('frame', 0))
            gene_id = attrs.get('ID', f"{seq_id}_{strand}_{frame}")

            # The score in GFF column 6 is the normalized score
            # We need the raw posterior which should be in coding_prob or similar
            raw_score = float(parts[5]) if parts[5] != '.' else 0.0

            scores[gene_id] = {
                'seq_id': seq_id,
                'strand': strand,
                'frame': frame,
                'score': raw_score,
                'damage_pct': float(attrs.get('damage_pct', 0.0)),
                'damage_signal': float(attrs.get('damage_signal', 0.5)),
                'length': int(parts[4]) - int(parts[3]) + 1,
            }

    return scores


def parse_mmseqs_hits(hits_path: Path) -> dict:
    """Parse MMseqs2 hits, keeping best hit per query."""
    best_hits = {}

    if not hits_path.exists():
        return best_hits

    with open(hits_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            query = parts[0]
            target = parts[1]
            identity = float(parts[2])
            alnlen = int(parts[3])
            qlen = int(parts[4])
            evalue = float(parts[6])
            bitscore = float(parts[7])

            coverage = alnlen / qlen if qlen > 0 else 0

            hit = {
                'target': target,
                'identity': identity,
                'coverage': coverage,
                'evalue': evalue,
                'bitscore': bitscore,
            }

            if query not in best_hits or bitscore > best_hits[query]['bitscore']:
                best_hits[query] = hit

    return best_hits


def compute_fdr_labels(target_hits: dict, decoy_hits: dict,
                       fdr_threshold: float = 0.01,
                       min_identity: float = 0.8,
                       min_coverage: float = 0.5) -> dict:
    """Compute FDR and assign pseudo-labels."""

    # Combine all hits with labels
    all_hits = []
    for query, hit in target_hits.items():
        all_hits.append((query, hit['bitscore'], False, hit))  # is_decoy=False
    for query, hit in decoy_hits.items():
        # Only add if better than target hit
        if query not in target_hits or hit['bitscore'] > target_hits[query]['bitscore']:
            all_hits.append((query, hit['bitscore'], True, hit))

    # Sort by bitscore descending
    all_hits.sort(key=lambda x: x[1], reverse=True)

    # Compute FDR at each score threshold
    n_target = 0
    n_decoy = 0
    labels = {}

    for query, score, is_decoy, hit in all_hits:
        if is_decoy:
            n_decoy += 1
        else:
            n_target += 1

        fdr = n_decoy / max(n_target, 1)

        # Label as correct if passes thresholds
        is_correct = (
            not is_decoy and
            fdr <= fdr_threshold and
            hit['identity'] >= min_identity and
            hit['coverage'] >= min_coverage
        )

        labels[query] = {
            'is_correct': is_correct,
            'fdr': fdr,
            'hit': hit,
        }

    return labels


def build_calibration_table(labeled_data: list,
                            score_bins: list,
                            length_bins: list,
                            damage_bins: list) -> tuple:
    """Build 3D calibration table."""

    n_score = len(score_bins) - 1
    n_length = len(length_bins) - 1
    n_damage = len(damage_bins) - 1

    counts = np.zeros((n_score, n_length, n_damage, 2))

    for data in labeled_data:
        raw_score = data['score']
        p = max(0.001, min(0.999, raw_score))
        logit_score = np.log(p / (1 - p))

        length = data['length']
        damage = data['damage_pct']

        s_bin = np.searchsorted(score_bins[1:], logit_score)
        l_bin = np.searchsorted(length_bins[1:], length)
        d_bin = np.searchsorted(damage_bins[1:], damage)

        s_bin = min(s_bin, n_score - 1)
        l_bin = min(l_bin, n_length - 1)
        d_bin = min(d_bin, n_damage - 1)

        counts[s_bin, l_bin, d_bin, 1] += 1
        if data['is_correct']:
            counts[s_bin, l_bin, d_bin, 0] += 1

    # Compute calibrated probabilities with smoothing
    calibration = np.zeros((n_score, n_length, n_damage))

    for s in range(n_score):
        for l in range(n_length):
            for d in range(n_damage):
                correct = counts[s, l, d, 0]
                total = counts[s, l, d, 1]
                calibration[s, l, d] = (correct + 1) / (total + 2)

    return calibration, counts


def generate_header(calibration: np.ndarray,
                    counts: np.ndarray,
                    score_bins: list,
                    length_bins: list,
                    damage_bins: list,
                    output_path: Path):
    """Generate C++ header file."""

    n_score, n_length, n_damage = calibration.shape

    with open(output_path, 'w') as f:
        f.write("""#pragma once
// Auto-generated calibration table from protein hit pseudo-labels
// Generated by build_calibration_table.py
//
// Uses target-decoy FDR estimation against KEGG database
// Input: 6-way frame posterior from unified_codon_scorer
// Output: P(protein-supported correct frame | score, length, damage)

#include <array>
#include <cmath>
#include <algorithm>

namespace dart {

""")

        # Score bins (logit scale)
        f.write(f"// Score bins (logit of 6-way posterior)\n")
        f.write(f"constexpr std::array<float, {len(score_bins)}> CALIB_SCORE_EDGES = {{{{\n    ")
        f.write(", ".join(f"{x:.4f}f" for x in score_bins))
        f.write("\n}};\n\n")

        # Length bins
        f.write(f"// Length bins (bp)\n")
        f.write(f"constexpr std::array<float, {len(length_bins)}> CALIB_LENGTH_EDGES = {{{{\n    ")
        f.write(", ".join(f"{x:.1f}f" for x in length_bins))
        f.write("\n}};\n\n")

        # Damage bins
        f.write(f"// Damage bins (%)\n")
        f.write(f"constexpr std::array<float, {len(damage_bins)}> CALIB_DAMAGE_EDGES = {{{{\n    ")
        f.write(", ".join(f"{x:.1f}f" for x in damage_bins))
        f.write("\n}};\n\n")

        # Calibration table
        f.write(f"// Calibration table: P(correct | score_bin, length_bin, damage_bin)\n")
        f.write(f"// Dimensions: [{n_score}][{n_length}][{n_damage}]\n")
        f.write(f"constexpr float CALIB_TABLE[{n_score}][{n_length}][{n_damage}] = {{{{\n")

        for s in range(n_score):
            f.write("    {  // score bin %d\n" % s)
            for l in range(n_length):
                f.write("        {")
                f.write(", ".join(f"{calibration[s,l,d]:.4f}f" for d in range(n_damage)))
                f.write("},\n")
            f.write("    },\n")
        f.write("}};\n\n")

        # Sample counts for reference
        f.write("// Sample counts per bin (for reference)\n")
        f.write(f"// Total samples: {int(counts[:,:,:,1].sum())}\n")
        f.write(f"// Total correct: {int(counts[:,:,:,0].sum())}\n")
        f.write(f"// Overall accuracy: {counts[:,:,:,0].sum() / max(counts[:,:,:,1].sum(), 1):.3f}\n\n")

        # Helper function
        f.write("""
inline int find_calib_bin(float value, const float* edges, int n_edges) {
    for (int i = 1; i < n_edges; ++i) {
        if (value < edges[i]) return i - 1;
    }
    return n_edges - 2;
}

// Calibrate 6-way posterior to protein-supported correctness probability
// Input: p_max = max of 6-way frame posterior (range ~0.17-0.5)
// Output: P(correct frame | p_max, length, damage)
inline float calibrate_posterior(float p_max, float length, float damage_pct) {
    // Convert posterior to logit for binning
    float p_clamped = std::clamp(p_max, 0.001f, 0.999f);
    float logit_score = std::log(p_clamped / (1.0f - p_clamped));

    int s_bin = find_calib_bin(logit_score, CALIB_SCORE_EDGES.data(), CALIB_SCORE_EDGES.size());
    int l_bin = find_calib_bin(length, CALIB_LENGTH_EDGES.data(), CALIB_LENGTH_EDGES.size());
    int d_bin = find_calib_bin(damage_pct, CALIB_DAMAGE_EDGES.data(), CALIB_DAMAGE_EDGES.size());

    s_bin = std::clamp(s_bin, 0, static_cast<int>(CALIB_SCORE_EDGES.size()) - 2);
    l_bin = std::clamp(l_bin, 0, static_cast<int>(CALIB_LENGTH_EDGES.size()) - 2);
    d_bin = std::clamp(d_bin, 0, static_cast<int>(CALIB_DAMAGE_EDGES.size()) - 2);

    return CALIB_TABLE[s_bin][l_bin][d_bin];
}

} // namespace dart
""")

    print(f"Generated header: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Build calibration table from protein hits')
    parser.add_argument('--results-dir', required=True, help='Directory with calibration results')
    parser.add_argument('--output-header', required=True, help='Output C++ header file')
    parser.add_argument('--fdr', type=float, default=0.01, help='FDR threshold')
    parser.add_argument('--min-identity', type=float, default=0.8, help='Min identity')
    parser.add_argument('--min-coverage', type=float, default=0.5, help='Min coverage')

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    predictions_dir = results_dir / 'predictions'
    hits_dir = results_dir / 'hits'

    if not predictions_dir.exists():
        print(f"ERROR: Predictions directory not found: {predictions_dir}")
        sys.exit(1)

    # Collect all samples
    samples = [d.name for d in predictions_dir.iterdir() if d.is_dir()]
    print(f"Found {len(samples)} samples")

    # Aggregate data from all samples
    all_labeled = []
    total_predictions = 0
    total_target_hits = 0
    total_decoy_hits = 0

    for sample in sorted(samples):
        gff_path = predictions_dir / sample / 'predictions.gff'
        target_hits_path = hits_dir / f'{sample}_target.tsv'
        decoy_hits_path = hits_dir / f'{sample}_decoy.tsv'

        if not gff_path.exists():
            print(f"  Skipping {sample}: no GFF")
            continue

        # Parse predictions
        scores = parse_gff_scores(gff_path)
        total_predictions += len(scores)

        # Parse hits
        target_hits = parse_mmseqs_hits(target_hits_path)
        decoy_hits = parse_mmseqs_hits(decoy_hits_path)
        total_target_hits += len(target_hits)
        total_decoy_hits += len(decoy_hits)

        # Compute FDR labels
        labels = compute_fdr_labels(target_hits, decoy_hits,
                                    args.fdr, args.min_identity, args.min_coverage)

        # Match predictions with labels
        for gene_id, score_data in scores.items():
            # Try to find matching hit (gene_id format varies)
            label_data = labels.get(gene_id)
            if label_data is None:
                # Try partial match
                for query_id in labels:
                    if gene_id in query_id or query_id in gene_id:
                        label_data = labels[query_id]
                        break

            if label_data is not None:
                all_labeled.append({
                    **score_data,
                    'is_correct': label_data['is_correct'],
                })

        print(f"  {sample}: {len(scores)} predictions, {len(target_hits)} target hits, {len(decoy_hits)} decoy hits")

    print(f"\nTotal: {total_predictions} predictions, {total_target_hits} target hits, {total_decoy_hits} decoy hits")
    print(f"Labeled: {len(all_labeled)}")

    if len(all_labeled) == 0:
        print("ERROR: No labeled data found")
        sys.exit(1)

    n_correct = sum(1 for d in all_labeled if d['is_correct'])
    print(f"Correct (FDR<{args.fdr}, id>={args.min_identity}, cov>={args.min_coverage}): {n_correct} ({100*n_correct/len(all_labeled):.1f}%)")

    # Define bins
    # Logit scale: logit(0.17)≈-1.6, logit(0.25)≈-1.1, logit(0.5)=0
    score_bins = [-4.0, -2.0, -1.5, -1.2, -1.0, -0.8, -0.5, 0.0, 0.5, 1.0, 2.0]
    length_bins = [0, 50, 75, 100, 125, 150, 200, 300, 500, 1000]
    damage_bins = [0, 2, 5, 10, 15, 20, 30, 50]

    # Build calibration table
    print("\nBuilding calibration table...")
    calibration, counts = build_calibration_table(all_labeled, score_bins, length_bins, damage_bins)

    print(f"  Score bins: {len(score_bins)-1}")
    print(f"  Length bins: {len(length_bins)-1}")
    print(f"  Damage bins: {len(damage_bins)-1}")
    print(f"  Min P(correct): {calibration.min():.3f}")
    print(f"  Max P(correct): {calibration.max():.3f}")
    print(f"  Mean P(correct): {calibration.mean():.3f}")

    # Generate header
    output_path = Path(args.output_header)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    generate_header(calibration, counts, score_bins, length_bins, damage_bins, output_path)

    print("\nDone!")


if __name__ == '__main__':
    main()
