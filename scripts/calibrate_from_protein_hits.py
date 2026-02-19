#!/usr/bin/env python3
"""
Calibrate AGP frame confidence using protein database hits as pseudo-labels.

Uses target-decoy approach to estimate FDR and generate calibration table
mapping logit(p_max) to P(correct | covariates).

Usage:
    python calibrate_from_protein_hits.py \
        --agp-output predictions.gff \
        --proteins predictions.faa \
        --target-db /path/to/uniref50.faa \
        --output-header score_calibration_protein.hpp \
        --threads 16
"""

import argparse
import subprocess
import tempfile
import gzip
import numpy as np
from pathlib import Path
from collections import defaultdict
import sys
import os


def parse_gff_scores(gff_path: str) -> dict:
    """Extract scores and metadata from AGP GFF output."""
    scores = {}

    opener = gzip.open if gff_path.endswith(".gz") else open
    with opener(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            # Parse attributes
            attrs = {}
            for attr in parts[8].split(";"):
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    attrs[key] = val

            seq_id = parts[0]
            strand = parts[6]
            frame = int(attrs.get("frame", 0))

            # Build unique gene ID
            gene_id = attrs.get("ID", f"{seq_id}_{strand}_{frame}")

            # Extract scores
            scores[gene_id] = {
                "seq_id": seq_id,
                "strand": strand,
                "frame": frame,
                "score": float(parts[5]) if parts[5] != "." else 0.0,
                "coding_prob": float(attrs.get("score", parts[5]))
                if parts[5] != "."
                else 0.0,
                "p_correct": float(attrs.get("p_correct", 0.5)),
                "damage_pct": float(attrs.get("damage_pct", 0.0)),
                "damage_signal": float(attrs.get("damage_signal", 0.5)),
                "length": int(parts[4]) - int(parts[3]) + 1,
            }

    return scores


def create_decoy_db(target_fasta: str, decoy_fasta: str):
    """Create decoy database by reversing protein sequences."""
    print(f"Creating decoy database: {decoy_fasta}")

    opener = gzip.open if target_fasta.endswith(".gz") else open
    with opener(target_fasta, "rt") as fin, open(decoy_fasta, "w") as fout:
        seq_id = None
        seq = []

        for line in fin:
            if line.startswith(">"):
                if seq_id:
                    # Write reversed sequence with DECOY_ prefix
                    fout.write(f">DECOY_{seq_id}\n")
                    rev_seq = "".join(seq)[::-1]
                    for i in range(0, len(rev_seq), 60):
                        fout.write(rev_seq[i : i + 60] + "\n")
                seq_id = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())

        # Write last sequence
        if seq_id:
            fout.write(f">DECOY_{seq_id}\n")
            rev_seq = "".join(seq)[::-1]
            for i in range(0, len(rev_seq), 60):
                fout.write(rev_seq[i : i + 60] + "\n")


def run_mmseqs_search(
    query_fasta: str, target_db: str, output_tsv: str, threads: int, tmp_dir: str
):
    """Run MMseqs2 search."""
    print(f"Running MMseqs2 search: {query_fasta} vs {target_db}")

    cmd = [
        "mmseqs",
        "easy-search",
        query_fasta,
        target_db,
        output_tsv,
        tmp_dir,
        "--format-output",
        "query,target,fident,alnlen,qlen,tlen,evalue,bits",
        "--threads",
        str(threads),
        "-e",
        "0.001",
        "--min-seq-id",
        "0.5",
        "-c",
        "0.5",
        "--cov-mode",
        "0",
        "-s",
        "5.7",
    ]

    subprocess.run(cmd, check=True)


def parse_mmseqs_hits(hits_tsv: str) -> dict:
    """Parse MMseqs2 hits, keeping best hit per query."""
    best_hits = {}

    with open(hits_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            query = parts[0]
            target = parts[1]
            identity = float(parts[2])
            alnlen = int(parts[3])
            qlen = int(parts[4])
            evalue = float(parts[6])
            bitscore = float(parts[7])

            is_decoy = target.startswith("DECOY_")
            coverage = alnlen / qlen if qlen > 0 else 0

            hit = {
                "target": target,
                "identity": identity,
                "coverage": coverage,
                "evalue": evalue,
                "bitscore": bitscore,
                "is_decoy": is_decoy,
            }

            # Keep best hit by bitscore
            if query not in best_hits or bitscore > best_hits[query]["bitscore"]:
                best_hits[query] = hit

    return best_hits


def compute_fdr(hits: dict, score_key: str = "bitscore") -> dict:
    """Compute FDR using target-decoy approach."""
    # Sort all hits by score descending
    all_hits = [(q, h) for q, h in hits.items()]
    all_hits.sort(key=lambda x: x[1][score_key], reverse=True)

    # Compute FDR at each threshold
    n_target = 0
    n_decoy = 0

    for query, hit in all_hits:
        if hit["is_decoy"]:
            n_decoy += 1
        else:
            n_target += 1

        # FDR = decoys / targets (or 2*decoys / (targets+decoys) for competitive)
        fdr = n_decoy / max(n_target, 1)
        hit["fdr"] = min(fdr, 1.0)
        hit["q_value"] = fdr  # Will be updated with monotonic q-value

    # Convert to q-values (monotonic minimum FDR)
    min_fdr = 1.0
    for query, hit in reversed(all_hits):
        min_fdr = min(min_fdr, hit["fdr"])
        hit["q_value"] = min_fdr

    return hits


def assign_pseudo_labels(
    scores: dict,
    hits: dict,
    fdr_threshold: float = 0.01,
    min_identity: float = 0.8,
    min_coverage: float = 0.5,
) -> dict:
    """Assign pseudo-labels based on protein hits."""
    labeled = {}

    for gene_id, score_data in scores.items():
        # Extract read name from gene_id (format: readname_+_0 or readname_gene1)
        # Try to match with protein FASTA header format
        read_name = score_data["seq_id"]

        # Find matching hit
        hit = None
        for query_id in hits:
            if read_name in query_id or query_id.startswith(read_name):
                hit = hits[query_id]
                break

        if hit is None:
            continue

        # Skip decoy hits
        if hit["is_decoy"]:
            continue

        # Check quality thresholds
        is_correct = (
            hit["q_value"] <= fdr_threshold
            and hit["identity"] >= min_identity
            and hit["coverage"] >= min_coverage
        )

        labeled[gene_id] = {
            **score_data,
            "hit": hit,
            "is_correct": is_correct,
        }

    return labeled


def build_calibration_table(
    labeled_data: dict, score_bins: list, length_bins: list, damage_bins: list
) -> np.ndarray:
    """Build 3D calibration table: P(correct | score_bin, length_bin, damage_bin)."""

    n_score = len(score_bins) - 1
    n_length = len(length_bins) - 1
    n_damage = len(damage_bins) - 1

    # Count correct and total per bin
    counts = np.zeros((n_score, n_length, n_damage, 2))  # [correct, total]

    for gene_id, data in labeled_data.items():
        # Compute logit(p_correct) as score; currently using coding_prob
        # since that's what's available in the output
        raw_score = data["score"]

        # Convert posterior to logit for binning
        # Clamp to avoid inf
        p = max(0.01, min(0.99, raw_score))
        logit_score = np.log(p / (1 - p))

        length = data["length"]
        damage = data["damage_pct"]

        # Find bins
        s_bin = np.searchsorted(score_bins[1:], logit_score)
        l_bin = np.searchsorted(length_bins[1:], length)
        d_bin = np.searchsorted(damage_bins[1:], damage)

        # Clamp to valid range
        s_bin = min(s_bin, n_score - 1)
        l_bin = min(l_bin, n_length - 1)
        d_bin = min(d_bin, n_damage - 1)

        counts[s_bin, l_bin, d_bin, 1] += 1
        if data["is_correct"]:
            counts[s_bin, l_bin, d_bin, 0] += 1

    # Compute calibrated probabilities with smoothing
    calibration = np.zeros((n_score, n_length, n_damage))

    for s in range(n_score):
        for l in range(n_length):
            for d in range(n_damage):
                correct = counts[s, l, d, 0]
                total = counts[s, l, d, 1]

                # Beta-binomial smoothing (prior: Beta(1,1) = uniform)
                calibration[s, l, d] = (correct + 1) / (total + 2)

    return calibration, counts


def generate_header(
    calibration: np.ndarray,
    score_bins: list,
    length_bins: list,
    damage_bins: list,
    output_path: str,
):
    """Generate C++ header file with calibration table."""

    n_score, n_length, n_damage = calibration.shape

    with open(output_path, "w") as f:
        f.write("""#pragma once
// Auto-generated calibration table from protein hit pseudo-labels
// Generated by calibrate_from_protein_hits.py
//
// Input: logit(p_max) where p_max is the 6-way frame posterior
// Output: P(protein-supported correct frame | score, length, damage)

#include <array>
#include <cmath>
#include <algorithm>

namespace agp {

""")

        # Score bins (logit scale)
        f.write(f"// Score bins (logit scale): logit(p_max)\n")
        f.write(
            f"constexpr std::array<float, {len(score_bins)}> SCORE_EDGES_LOGIT = {{{{\n    "
        )
        f.write(", ".join(f"{x:.4f}f" for x in score_bins))
        f.write("\n}};\n\n")

        # Length bins
        f.write(f"// Length bins (bp)\n")
        f.write(
            f"constexpr std::array<float, {len(length_bins)}> LENGTH_EDGES_CALIB = {{{{\n    "
        )
        f.write(", ".join(f"{x:.1f}f" for x in length_bins))
        f.write("\n}};\n\n")

        # Damage bins
        f.write(f"// Damage bins (%)\n")
        f.write(
            f"constexpr std::array<float, {len(damage_bins)}> DAMAGE_EDGES_CALIB = {{{{\n    "
        )
        f.write(", ".join(f"{x:.1f}f" for x in damage_bins))
        f.write("\n}};\n\n")

        # Calibration table
        f.write(
            f"// Calibration table: P(correct | score_bin, length_bin, damage_bin)\n"
        )
        f.write(f"// Dimensions: [{n_score}][{n_length}][{n_damage}]\n")
        f.write(
            f"constexpr float CALIBRATION_PROTEIN[{n_score}][{n_length}][{n_damage}] = {{{{\n"
        )

        for s in range(n_score):
            f.write("    {\n")
            for l in range(n_length):
                f.write("        {")
                f.write(
                    ", ".join(f"{calibration[s, l, d]:.4f}f" for d in range(n_damage))
                )
                f.write("},\n")
            f.write("    },\n")
        f.write("}};\n\n")

        # Helper function
        f.write("""
inline int find_bin_logit(float value, const float* edges, int n_edges) {
    for (int i = 1; i < n_edges; ++i) {
        if (value < edges[i]) return i - 1;
    }
    return n_edges - 2;
}

// Calibrate using logit(p_max) as input
// p_max = max posterior from 6-way frame selection
inline float calibrate_from_posterior(float p_max, float length, float damage_pct) {
    // Convert posterior to logit
    float p_clamped = std::clamp(p_max, 0.01f, 0.99f);
    float logit_score = std::log(p_clamped / (1.0f - p_clamped));

    int s_bin = find_bin_logit(logit_score, SCORE_EDGES_LOGIT.data(), SCORE_EDGES_LOGIT.size());
    int l_bin = find_bin_logit(length, LENGTH_EDGES_CALIB.data(), LENGTH_EDGES_CALIB.size());
    int d_bin = find_bin_logit(damage_pct, DAMAGE_EDGES_CALIB.data(), DAMAGE_EDGES_CALIB.size());

    s_bin = std::clamp(s_bin, 0, (int)SCORE_EDGES_LOGIT.size() - 2);
    l_bin = std::clamp(l_bin, 0, (int)LENGTH_EDGES_CALIB.size() - 2);
    d_bin = std::clamp(d_bin, 0, (int)DAMAGE_EDGES_CALIB.size() - 2);

    return CALIBRATION_PROTEIN[s_bin][l_bin][d_bin];
}

} // namespace agp
""")

    print(f"Generated header: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Calibrate AGP using protein hits")
    parser.add_argument("--agp-gff", required=True, help="AGP GFF output file")
    parser.add_argument("--proteins", required=True, help="AGP protein FASTA output")
    parser.add_argument(
        "--target-db", required=True, help="Target protein database (FASTA)"
    )
    parser.add_argument(
        "--output-header",
        default="score_calibration_protein.hpp",
        help="Output C++ header file",
    )
    parser.add_argument("--threads", type=int, default=16, help="Number of threads")
    parser.add_argument(
        "--fdr", type=float, default=0.01, help="FDR threshold for pseudo-labels"
    )
    parser.add_argument(
        "--min-identity", type=float, default=0.8, help="Min identity for correct label"
    )
    parser.add_argument(
        "--min-coverage", type=float, default=0.5, help="Min coverage for correct label"
    )
    parser.add_argument(
        "--tmp-dir", default="/scratch/tmp/calibration", help="Temp directory"
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.agp_gff):
        print(f"ERROR: GFF file not found: {args.agp_gff}")
        sys.exit(1)
    if not os.path.exists(args.proteins):
        print(f"ERROR: Protein file not found: {args.proteins}")
        sys.exit(1)
    if not os.path.exists(args.target_db):
        print(f"ERROR: Target DB not found: {args.target_db}")
        sys.exit(1)

    os.makedirs(args.tmp_dir, exist_ok=True)

    # Step 1: Parse AGP scores
    print("Step 1: Parsing AGP scores...")
    scores = parse_gff_scores(args.agp_gff)
    print(f"  Loaded {len(scores)} predictions")

    # Step 2: Create decoy database
    decoy_db = os.path.join(args.tmp_dir, "decoy.faa")
    create_decoy_db(args.target_db, decoy_db)

    # Step 3: Concatenate target + decoy
    combined_db = os.path.join(args.tmp_dir, "target_decoy.faa")
    print("Step 3: Creating combined target-decoy database...")
    with open(combined_db, "w") as fout:
        # Copy target
        opener = gzip.open if args.target_db.endswith(".gz") else open
        with opener(args.target_db, "rt") as fin:
            fout.write(fin.read())
        # Append decoy
        with open(decoy_db) as fin:
            fout.write(fin.read())

    # Step 4: Run MMseqs2
    hits_tsv = os.path.join(args.tmp_dir, "hits.tsv")
    mmseqs_tmp = os.path.join(args.tmp_dir, "mmseqs")
    os.makedirs(mmseqs_tmp, exist_ok=True)
    run_mmseqs_search(args.proteins, combined_db, hits_tsv, args.threads, mmseqs_tmp)

    # Step 5: Parse hits and compute FDR
    print("Step 5: Computing FDR...")
    hits = parse_mmseqs_hits(hits_tsv)
    hits = compute_fdr(hits)

    n_target = sum(1 for h in hits.values() if not h["is_decoy"])
    n_decoy = sum(1 for h in hits.values() if h["is_decoy"])
    print(f"  Target hits: {n_target}, Decoy hits: {n_decoy}")

    # Step 6: Assign pseudo-labels
    print("Step 6: Assigning pseudo-labels...")
    labeled = assign_pseudo_labels(
        scores, hits, args.fdr, args.min_identity, args.min_coverage
    )

    n_correct = sum(1 for d in labeled.values() if d["is_correct"])
    n_total = len(labeled)
    print(
        f"  Labeled: {n_total}, Correct: {n_correct} ({100 * n_correct / max(n_total, 1):.1f}%)"
    )

    # Step 7: Build calibration table
    print("Step 7: Building calibration table...")

    # Define bins
    # Score bins in logit scale: logit(0.17) ≈ -1.6, logit(0.5) = 0, logit(0.9) ≈ 2.2
    score_bins = [-3.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
    length_bins = [0, 50, 75, 100, 150, 200, 300, 500, 1000, 2000]
    damage_bins = [0, 2, 5, 10, 20, 30, 50, 100]

    calibration, counts = build_calibration_table(
        labeled, score_bins, length_bins, damage_bins
    )

    # Print summary
    print("\nCalibration table summary:")
    print(f"  Score bins: {len(score_bins) - 1}")
    print(f"  Length bins: {len(length_bins) - 1}")
    print(f"  Damage bins: {len(damage_bins) - 1}")
    print(f"  Total cells: {calibration.size}")
    print(f"  Min P(correct): {calibration.min():.3f}")
    print(f"  Max P(correct): {calibration.max():.3f}")
    print(f"  Mean P(correct): {calibration.mean():.3f}")

    # Step 8: Generate header
    generate_header(
        calibration, score_bins, length_bins, damage_bins, args.output_header
    )

    print("\nDone!")


if __name__ == "__main__":
    main()
