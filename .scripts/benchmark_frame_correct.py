#!/usr/bin/env python3
"""
CORRECT Frame Accuracy Benchmark for AGP

This script properly evaluates frame selection accuracy by matching
predicted proteins against the known correct protein translation from
the benchmark metadata (damaged_seq_inframe_aa).

IMPORTANT: The "true frame" cannot be determined from coordinates alone!
The only reliable method is protein sequence matching.

Previous benchmarks used:
    true_frame = start % 3  # WRONG - produces ~33% (random) accuracy

Correct method:
    true_frame = find_frame_matching_known_protein(seq, damaged_seq_inframe_aa)
    # Produces accurate measurements: ~77% strand, ~93% frame|strand

Usage:
    python3 .scripts/benchmark_frame_correct.py \\
        <aa-damage.tsv.gz> \\
        <predictions.gff>

Example:
    # Run AGP first
    ./build_release/agp -i damaged_reads.fa.gz -o predictions.gff

    # Then evaluate
    python3 .scripts/benchmark_frame_correct.py \\
        /projects/.../119_B3_116_L0_KapK-12-1-24_aa-damage.tsv.gz \\
        predictions.gff

Expected Results (on synthetic damaged data, ~25% d_max):
    - Strand accuracy: ~77%
    - Frame accuracy (given correct strand): ~93%
    - Combined (strand + frame): ~72%
"""

import gzip
import sys
from collections import defaultdict

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def translate(seq, frame=0):
    """Translate DNA to protein in given frame."""
    protein = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)

def rc(seq):
    """Reverse complement."""
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    return ''.join(comp.get(c, 'N') for c in reversed(seq.upper()))

def find_true_frame(seq, target_protein, min_match_frac=0.8):
    """
    Find which of 6 frames produces the target protein.

    Returns (strand, frame) where strand is '+' or '-' and frame is 0, 1, or 2.
    Returns (None, None) if no frame matches well enough.
    """
    if not target_protein or len(target_protein) < 3:
        return None, None

    best_match = 0
    best_frame = None
    best_strand = None

    # Try all 6 frames
    for strand, s in [('+', seq), ('-', rc(seq))]:
        for frame in [0, 1, 2]:
            prot = translate(s, frame)
            min_len = min(len(prot), len(target_protein))
            if min_len < 3:
                continue
            matches = sum(a == b for a, b in zip(prot[:min_len], target_protein[:min_len]))
            if matches > best_match:
                best_match = matches
                best_frame = frame
                best_strand = strand

    # Require at least min_match_frac of the target to match
    if best_match >= len(target_protein) * min_match_frac:
        return best_strand, best_frame
    return None, None

def load_ground_truth(tsv_path, max_reads=None):
    """
    Load ground truth from aMGSIM benchmark TSV.

    Key columns:
    - read_name: Unique read identifier
    - damaged_seq: The damaged DNA sequence
    - damaged_seq_inframe_aa: Correct protein when translated in true frame
    """
    truth = {}
    opener = gzip.open if tsv_path.endswith('.gz') else open

    with opener(tsv_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        col = {name: i for i, name in enumerate(header)}

        required_cols = ['read_name', 'damaged_seq', 'damaged_seq_inframe_aa']
        for c in required_cols:
            if c not in col:
                print(f"Error: Missing required column '{c}' in TSV")
                print(f"Available columns: {list(col.keys())}")
                sys.exit(1)

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) <= max(col.values()):
                continue

            read_name = fields[col['read_name']]
            damaged_seq = fields[col['damaged_seq']]
            true_protein = fields[col['damaged_seq_inframe_aa']]

            if not damaged_seq or not true_protein:
                continue

            # Find which frame produces the correct protein
            true_strand, true_frame = find_true_frame(damaged_seq, true_protein)

            if true_strand is not None:
                truth[read_name] = {
                    'strand': true_strand,
                    'frame': true_frame,
                    'seq': damaged_seq,
                    'protein': true_protein
                }

            if max_reads and len(truth) >= max_reads:
                break

    return truth

def load_predictions(gff_path):
    """Load AGP predictions from GFF file."""
    preds = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            read_name = fields[0]
            strand = fields[6]

            # Parse attributes
            attrs = {}
            for kv in fields[8].split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    attrs[k] = v

            frame = int(attrs.get('frame', -1))

            if read_name not in preds:
                preds[read_name] = {'strand': strand, 'frame': frame}

    return preds

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print(f"\nUsage: {sys.argv[0]} <aa-damage.tsv.gz> <predictions.gff> [max_reads]")
        sys.exit(1)

    tsv_path = sys.argv[1]
    gff_path = sys.argv[2]
    max_reads = int(sys.argv[3]) if len(sys.argv) > 3 else 10000

    print(f"Loading ground truth from {tsv_path}...")
    truth = load_ground_truth(tsv_path, max_reads)
    print(f"  {len(truth)} reads with determinable true frame")

    print(f"\nLoading predictions from {gff_path}...")
    preds = load_predictions(gff_path)
    print(f"  {len(preds)} predictions")

    # Evaluate
    strand_correct = 0
    frame_correct_given_strand = 0
    both_correct = 0
    total = 0

    strand_confusion = defaultdict(lambda: defaultdict(int))
    frame_confusion = defaultdict(lambda: defaultdict(int))

    for read_name, t in truth.items():
        if read_name not in preds:
            continue

        p = preds[read_name]
        total += 1

        # Strand accuracy
        strand_confusion[t['strand']][p['strand']] += 1

        if p['strand'] == t['strand']:
            strand_correct += 1

            # Frame accuracy (only when strand is correct)
            frame_confusion[t['frame']][p['frame']] += 1

            if p['frame'] == t['frame']:
                frame_correct_given_strand += 1
                both_correct += 1

    # Results
    print(f"\n{'='*60}")
    print("FRAME ACCURACY RESULTS (Correct Method)")
    print(f"{'='*60}")
    print(f"Total matched: {total}")
    print()
    print(f"Strand accuracy:           {100*strand_correct/total:.1f}%  ({strand_correct}/{total})")
    if strand_correct > 0:
        print(f"Frame accuracy|strand:     {100*frame_correct_given_strand/strand_correct:.1f}%  ({frame_correct_given_strand}/{strand_correct})")
    print(f"Combined (strand+frame):   {100*both_correct/total:.1f}%  ({both_correct}/{total})")

    print(f"\nStrand confusion matrix:")
    for true_s in ['+', '-']:
        row = strand_confusion[true_s]
        total_row = sum(row.values())
        if total_row > 0:
            print(f"  True {true_s}: ", end='')
            for pred_s in ['+', '-']:
                pct = 100 * row[pred_s] / total_row
                print(f"pred_{pred_s}={pct:.1f}% ", end='')
            print()

    print(f"\nFrame confusion matrix (given correct strand):")
    for true_f in [0, 1, 2]:
        row = frame_confusion[true_f]
        total_row = sum(row.values())
        if total_row > 0:
            print(f"  True {true_f}: ", end='')
            for pred_f in [0, 1, 2]:
                pct = 100 * row[pred_f] / total_row
                print(f"pred_{pred_f}={pct:.1f}% ", end='')
            print()

    # Random baseline comparison
    print(f"\n{'='*60}")
    print("COMPARISON TO RANDOM")
    print(f"{'='*60}")
    print(f"Random strand: 50.0%  | AGP: {100*strand_correct/total:.1f}%  | Improvement: {100*strand_correct/total - 50:.1f}pp")
    print(f"Random frame:  33.3%  | AGP: {100*frame_correct_given_strand/strand_correct:.1f}%  | Improvement: {100*frame_correct_given_strand/strand_correct - 33.3:.1f}pp" if strand_correct > 0 else "")
    print(f"Random both:   16.7%  | AGP: {100*both_correct/total:.1f}%  | Improvement: {100*both_correct/total - 16.7:.1f}pp")

if __name__ == '__main__':
    main()
