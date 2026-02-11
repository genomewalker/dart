#!/usr/bin/env python3
"""
Self-supervised frame classifier training.

Uses CDS sequences where frame 0 is ALWAYS correct (or synthetic coding-like sequences).
No labeled data needed - learns from the structure of coding sequences.

Outputs a small neural network as a C++ header file.

Architecture:
- Input: 30 features per frame
- Hidden: 16 neurons with ReLU
- Output: 1 (probability frame is correct)
"""

import sys
import os
import gzip
import argparse
import numpy as np
from collections import defaultdict
import random

try:
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

# Codon table
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

STOP_CODONS = {'TAA', 'TAG', 'TGA'}
RARE_CODONS = {'AGG', 'AGA', 'CGA', 'CGG', 'CTA', 'ATA'}

# Expected codon frequencies (typical prokaryote)
EXPECTED_CODON_FREQ = {
    'TTT': 0.022, 'TTC': 0.016, 'TTA': 0.014, 'TTG': 0.013,
    'TCT': 0.009, 'TCC': 0.009, 'TCA': 0.007, 'TCG': 0.009,
    'TAT': 0.016, 'TAC': 0.012, 'TAA': 0.002, 'TAG': 0.000,
    'TGT': 0.005, 'TGC': 0.006, 'TGA': 0.001, 'TGG': 0.015,
    'CTT': 0.011, 'CTC': 0.011, 'CTA': 0.004, 'CTG': 0.052,
    'CCT': 0.007, 'CCC': 0.006, 'CCA': 0.008, 'CCG': 0.023,
    'CAT': 0.013, 'CAC': 0.010, 'CAA': 0.015, 'CAG': 0.029,
    'CGT': 0.021, 'CGC': 0.022, 'CGA': 0.004, 'CGG': 0.006,
    'ATT': 0.030, 'ATC': 0.025, 'ATA': 0.005, 'ATG': 0.028,
    'ACT': 0.009, 'ACC': 0.023, 'ACA': 0.007, 'ACG': 0.015,
    'AAT': 0.018, 'AAC': 0.022, 'AAA': 0.034, 'AAG': 0.010,
    'AGT': 0.009, 'AGC': 0.016, 'AGA': 0.002, 'AGG': 0.001,
    'GTT': 0.018, 'GTC': 0.015, 'GTA': 0.011, 'GTG': 0.026,
    'GCT': 0.015, 'GCC': 0.026, 'GCA': 0.020, 'GCG': 0.033,
    'GAT': 0.032, 'GAC': 0.019, 'GAA': 0.040, 'GAG': 0.018,
    'GGT': 0.025, 'GGC': 0.030, 'GGA': 0.008, 'GGG': 0.011
}

NUM_FEATURES = 30


def reverse_complement(seq):
    """Reverse complement a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))


def translate(seq, frame):
    """Translate DNA to protein starting at given frame."""
    protein = ''
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        protein += CODON_TABLE.get(codon, 'X')
    return protein


def extract_features(seq, frame):
    """
    Extract 30 features for a single frame.
    """
    seq = seq.upper()
    features = []
    
    # Get codons for this frame
    codons = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if all(c in 'ACGT' for c in codon):
            codons.append(codon)
    
    n_codons = len(codons)
    if n_codons < 3:
        return [0.0] * NUM_FEATURES
    
    # 1. Internal stop codon count (1)
    internal_stops = sum(1 for c in codons[:-1] if c in STOP_CODONS)
    features.append(min(internal_stops, 5) / 5.0)
    
    # 2. Stop codon indicator (1)
    features.append(1.0 if internal_stops > 0 else 0.0)
    
    # 3. First stop position (1)
    first_stop = -1
    for i, c in enumerate(codons[:-1]):
        if c in STOP_CODONS:
            first_stop = i
            break
    features.append(first_stop / n_codons if first_stop >= 0 else 1.0)
    
    # 4. Rare codon frequency (1)
    rare_count = sum(1 for c in codons if c in RARE_CODONS)
    features.append(rare_count / n_codons)
    
    # 5. GC content at each codon position (3)
    gc_pos = [0, 0, 0]
    total_pos = [0, 0, 0]
    for i in range(frame, len(seq)):
        pos = (i - frame) % 3
        c = seq[i]
        if c in 'ACGT':
            total_pos[pos] += 1
            if c in 'GC':
                gc_pos[pos] += 1
    for p in range(3):
        features.append(gc_pos[p] / max(1, total_pos[p]))
    
    # 6. Nucleotide frequency at each codon position (12)
    nt_freq = [[0, 0, 0, 0] for _ in range(3)]
    nt_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(frame, len(seq)):
        pos = (i - frame) % 3
        c = seq[i]
        if c in nt_map:
            nt_freq[pos][nt_map[c]] += 1
    for pos in range(3):
        total = sum(nt_freq[pos]) or 1
        for nt in range(4):
            features.append(nt_freq[pos][nt] / total)
    
    # 7. RNY rule compliance (1)
    rny_count = 0
    for i in range(frame, len(seq) - 2, 3):
        c0, c2 = seq[i], seq[i+2]
        if c0 in 'AG' and c2 in 'CT':
            rny_count += 1
    features.append(rny_count / n_codons)
    
    # 8. Codon usage bias (1)
    codon_counts = defaultdict(int)
    for c in codons:
        codon_counts[c] += 1
    bias = 0.0
    for codon, count in codon_counts.items():
        observed = count / n_codons
        expected = EXPECTED_CODON_FREQ.get(codon, 0.01)
        bias += (observed - expected) ** 2 / (expected + 0.001)
    features.append(min(bias / 10.0, 1.0))
    
    # 9. Dipeptide entropy (1)
    protein = translate(seq, frame)
    if len(protein) >= 2:
        dipeptides = defaultdict(int)
        for i in range(len(protein) - 1):
            dp = protein[i:i+2]
            dipeptides[dp] += 1
        total_dp = sum(dipeptides.values())
        entropy = 0.0
        for count in dipeptides.values():
            p = count / total_dp
            if p > 0:
                entropy -= p * np.log2(p)
        features.append(entropy / 8.64)
    else:
        features.append(0.5)
    
    # 10. Length feature (1)
    features.append(min(1.0, n_codons / 50.0))
    
    # 11. Start codon presence (1)
    features.append(1.0 if codons and codons[0] == 'ATG' else 0.0)
    
    # 12. Terminal stop presence (1)
    features.append(1.0 if codons and codons[-1] in STOP_CODONS else 0.0)
    
    # 13. Forbidden hexamer count (1)
    forbidden = 0
    for i in range(frame, len(seq) - 5, 3):
        if seq[i:i+3] in STOP_CODONS:
            forbidden += 1
    features.append(min(forbidden / 3.0, 1.0))
    
    # 14. Pyrimidine at position 3 (1)
    pyr_at_3 = sum(1 for i in range(frame + 2, len(seq), 3) if seq[i] in 'CT')
    features.append(pyr_at_3 / max(1, n_codons))
    
    # 15-17. Dinucleotide bias (3)
    cpg_count = seq.count('CG')
    tpa_count = seq.count('TA')
    aa_tt_count = seq.count('AA') + seq.count('TT')
    seq_len = len(seq)
    features.append(cpg_count / max(1, seq_len // 2) * 4)
    features.append(tpa_count / max(1, seq_len // 2) * 4)
    features.append(aa_tt_count / max(1, seq_len // 2) * 2)
    
    assert len(features) == NUM_FEATURES, f"Expected {NUM_FEATURES}, got {len(features)}"
    return features


def extract_all_frames(seq):
    """Extract features for all 3 frames."""
    return [extract_features(seq, f) for f in range(3)]


def generate_coding_sequences(n_seqs=5000, min_codons=20, max_codons=50):
    """Generate synthetic coding sequences with realistic codon usage."""
    print(f"Generating {n_seqs} synthetic coding sequences...", file=sys.stderr)
    
    sequences = []
    codons = list(EXPECTED_CODON_FREQ.keys())
    weights = [EXPECTED_CODON_FREQ[c] + 0.001 for c in codons]
    
    for _ in range(n_seqs):
        n_codons = random.randint(min_codons, max_codons)
        seq = ''
        for _ in range(n_codons):
            codon = random.choices(codons, weights=weights)[0]
            # Avoid internal stops
            if codon in STOP_CODONS and len(seq) > 0:
                codon = 'GCT'  # Replace with alanine
            seq += codon
        sequences.append(seq)
    
    return sequences


def simulate_reads(cds_sequences, read_len=100, n_reads=50000):
    """Simulate short reads from CDS sequences."""
    print(f"Simulating {n_reads} reads...", file=sys.stderr)
    
    reads = []
    true_frames = []
    
    for _ in range(n_reads):
        cds = random.choice(cds_sequences)
        if len(cds) < read_len + 10:
            continue
        
        max_start = len(cds) - read_len
        start = random.randint(0, max_start)
        true_frame = start % 3
        read = cds[start:start + read_len]
        
        if len(read) == read_len and 'N' not in read:
            reads.append(read)
            true_frames.append(true_frame)
    
    return reads, true_frames


def create_training_data(reads, true_frames):
    """Create training dataset: correct frame = 1, wrong frames = 0."""
    X = []
    y = []
    
    for read, true_frame in zip(reads, true_frames):
        features = extract_all_frames(read)
        for frame in range(3):
            X.append(features[frame])
            y.append(1 if frame == true_frame else 0)
    
    return np.array(X, dtype=np.float32), np.array(y, dtype=np.int32)


def evaluate_frame_selection(model, scaler, reads, true_frames):
    """Evaluate frame selection accuracy (pick best of 3)."""
    correct = 0
    for read, true_frame in zip(reads, true_frames):
        features = np.array(extract_all_frames(read), dtype=np.float32)
        features_scaled = scaler.transform(features)
        probs = model.predict_proba(features_scaled)[:, 1]
        predicted = np.argmax(probs)
        if predicted == true_frame:
            correct += 1
    return correct / len(reads)


def load_ground_truth(tsv_path, max_samples=5000):
    """Load aMGSIM ground truth for validation."""
    print(f"Loading validation from {tsv_path}...", file=sys.stderr)
    
    data = []
    opener = gzip.open if tsv_path.endswith('.gz') else open
    mode = 'rt' if tsv_path.endswith('.gz') else 'r'
    
    with opener(tsv_path, mode) as f:
        header = f.readline().strip().split('\t')
        cols = {name: i for i, name in enumerate(header)}
        
        required = ['damaged_seq', 'damaged_seq_inframe_nt', 'Strand_read', 'Strand_gene']
        if not all(c in cols for c in required):
            print(f"  Missing columns", file=sys.stderr)
            return []
        
        for line in f:
            if len(data) >= max_samples:
                break
            
            parts = line.strip().split('\t')
            if len(parts) <= max(cols.values()):
                continue
            
            damaged_seq = parts[cols['damaged_seq']].upper()
            inframe_nt = parts[cols['damaged_seq_inframe_nt']].upper()
            strand_read = parts[cols['Strand_read']]
            strand_gene = parts[cols['Strand_gene']]
            
            if not damaged_seq or len(damaged_seq) < 30:
                continue
            
            is_forward = (strand_read == strand_gene)
            if not is_forward:
                damaged_seq = reverse_complement(damaged_seq)
            
            if not inframe_nt or len(inframe_nt) < 6:
                continue
            
            # Find true frame
            true_frame = None
            for frame in range(3):
                test_seq = damaged_seq[frame:]
                if inframe_nt in test_seq:
                    pos = test_seq.find(inframe_nt)
                    if pos % 3 == 0:
                        true_frame = frame
                        break
            
            if true_frame is not None:
                data.append((damaged_seq, true_frame))
    
    print(f"  Loaded {len(data)} samples", file=sys.stderr)
    return data


def export_to_cpp(model, scaler, output_path):
    """Export trained MLP to C++ header."""
    print(f"Exporting to {output_path}...", file=sys.stderr)
    
    # Extract weights from sklearn MLP
    W1 = model.coefs_[0]  # (input_size, hidden_size)
    b1 = model.intercepts_[0]  # (hidden_size,)
    W2 = model.coefs_[1]  # (hidden_size, output_size)
    b2 = model.intercepts_[1]  # (output_size,)
    
    input_size = W1.shape[0]
    hidden_size = W1.shape[1]
    
    with open(output_path, 'w') as f:
        f.write("#pragma once\n\n")
        f.write("/**\n")
        f.write(" * Auto-generated frame classifier neural network.\n")
        f.write(" * Generated by train_frame_nn.py\n")
        f.write(" * \n")
        f.write(f" * Architecture: Input({input_size}) -> Hidden({hidden_size}, ReLU) -> Output(2, Softmax)\n")
        f.write(" */\n\n")
        f.write("#include <array>\n")
        f.write("#include <cmath>\n")
        f.write("#include <algorithm>\n\n")
        f.write("namespace agp {\n")
        f.write("namespace frame_nn {\n\n")
        
        f.write(f"constexpr int INPUT_SIZE = {input_size};\n")
        f.write(f"constexpr int HIDDEN_SIZE = {hidden_size};\n\n")
        
        # Input normalization
        f.write("// Input normalization (StandardScaler)\n")
        f.write("constexpr std::array<float, INPUT_SIZE> INPUT_MEAN = {{\n    ")
        f.write(",\n    ".join(f"{x:.8f}f" for x in scaler.mean_))
        f.write("\n}};\n\n")
        
        f.write("constexpr std::array<float, INPUT_SIZE> INPUT_STD = {{\n    ")
        f.write(",\n    ".join(f"{max(x, 1e-8):.8f}f" for x in scaler.scale_))
        f.write("\n}};\n\n")
        
        # Layer 1
        f.write("// Layer 1: INPUT_SIZE -> HIDDEN_SIZE\n")
        f.write("constexpr std::array<std::array<float, HIDDEN_SIZE>, INPUT_SIZE> W1 = {{\n")
        for i, row in enumerate(W1):
            f.write("    {{" + ", ".join(f"{x:.8f}f" for x in row) + "}}")
            f.write(",\n" if i < len(W1) - 1 else "\n")
        f.write("}};\n\n")
        
        f.write("constexpr std::array<float, HIDDEN_SIZE> B1 = {{\n    ")
        f.write(", ".join(f"{x:.8f}f" for x in b1))
        f.write("\n}};\n\n")
        
        # Layer 2 (output layer has 2 outputs for binary classification)
        f.write("// Layer 2: HIDDEN_SIZE -> 2 (class probabilities)\n")
        f.write("constexpr std::array<std::array<float, 2>, HIDDEN_SIZE> W2 = {{\n")
        for i, row in enumerate(W2):
            f.write("    {{" + ", ".join(f"{x:.8f}f" for x in row) + "}}")
            f.write(",\n" if i < len(W2) - 1 else "\n")
        f.write("}};\n\n")
        
        f.write("constexpr std::array<float, 2> B2 = {{\n    ")
        f.write(", ".join(f"{x:.8f}f" for x in b2))
        f.write("\n}};\n\n")
        
        # Forward pass
        f.write("/**\n")
        f.write(" * Score a reading frame using the trained neural network.\n")
        f.write(" * @param features Array of 30 features for this frame\n")
        f.write(" * @return Probability (0-1) that this is the correct frame\n")
        f.write(" */\n")
        f.write("inline float score_frame(const std::array<float, INPUT_SIZE>& features) {\n")
        f.write("    // Normalize input\n")
        f.write("    std::array<float, INPUT_SIZE> x;\n")
        f.write("    for (int i = 0; i < INPUT_SIZE; i++) {\n")
        f.write("        x[i] = (features[i] - INPUT_MEAN[i]) / INPUT_STD[i];\n")
        f.write("    }\n\n")
        f.write("    // Hidden layer with ReLU\n")
        f.write("    std::array<float, HIDDEN_SIZE> h;\n")
        f.write("    for (int j = 0; j < HIDDEN_SIZE; j++) {\n")
        f.write("        float sum = B1[j];\n")
        f.write("        for (int i = 0; i < INPUT_SIZE; i++) {\n")
        f.write("            sum += x[i] * W1[i][j];\n")
        f.write("        }\n")
        f.write("        h[j] = std::max(0.0f, sum);  // ReLU\n")
        f.write("    }\n\n")
        f.write("    // Output layer (softmax, return prob of class 1)\n")
        f.write("    float z0 = B2[0], z1 = B2[1];\n")
        f.write("    for (int j = 0; j < HIDDEN_SIZE; j++) {\n")
        f.write("        z0 += h[j] * W2[j][0];\n")
        f.write("        z1 += h[j] * W2[j][1];\n")
        f.write("    }\n")
        f.write("    // Softmax: exp(z1) / (exp(z0) + exp(z1)) = 1 / (1 + exp(z0-z1))\n")
        f.write("    float diff = z0 - z1;\n")
        f.write("    if (diff > 20.0f) return 0.0f;\n")
        f.write("    if (diff < -20.0f) return 1.0f;\n")
        f.write("    return 1.0f / (1.0f + std::exp(diff));\n")
        f.write("}\n\n")
        
        f.write("} // namespace frame_nn\n")
        f.write("} // namespace agp\n")
    
    print(f"  Done!", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Train self-supervised frame classifier')
    parser.add_argument('--cds', nargs='*', help='CDS FASTA files (optional)')
    parser.add_argument('--validation', nargs='*', help='aMGSIM TSV files for validation')
    parser.add_argument('--output', default='frame_nn', help='Output base name')
    parser.add_argument('--n-reads', type=int, default=50000, help='Number of simulated reads')
    parser.add_argument('--read-len', type=int, default=100, help='Read length')
    parser.add_argument('--hidden-size', type=int, default=16, help='Hidden layer size')
    args = parser.parse_args()
    
    if not HAS_SKLEARN:
        print("Error: sklearn required. Install with: pip install scikit-learn", file=sys.stderr)
        sys.exit(1)
    
    # Generate synthetic CDS (frame 0 is always correct)
    cds_sequences = generate_coding_sequences(n_seqs=5000)
    
    # Simulate reads
    reads, true_frames = simulate_reads(cds_sequences, 
                                        read_len=args.read_len,
                                        n_reads=args.n_reads)
    
    # Split train/test
    n_test = min(5000, len(reads) // 5)
    test_reads, test_frames = reads[-n_test:], true_frames[-n_test:]
    train_reads, train_frames = reads[:-n_test], true_frames[:-n_test]
    
    # Create training data
    print("Extracting features...", file=sys.stderr)
    X_train, y_train = create_training_data(train_reads, train_frames)
    X_test, y_test = create_training_data(test_reads, test_frames)
    
    print(f"Training: {len(X_train)} samples ({y_train.sum()} positive)", file=sys.stderr)
    print(f"Test: {len(X_test)} samples ({y_test.sum()} positive)", file=sys.stderr)
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train MLP
    print(f"\nTraining MLP (hidden_size={args.hidden_size})...", file=sys.stderr)
    model = MLPClassifier(
        hidden_layer_sizes=(args.hidden_size,),
        activation='relu',
        max_iter=200,
        random_state=42,
        early_stopping=True,
        validation_fraction=0.1,
        verbose=True
    )
    model.fit(X_train_scaled, y_train)
    
    # Evaluate
    print("\n=== Results ===", file=sys.stderr)
    train_acc = model.score(X_train_scaled, y_train)
    test_acc = model.score(X_test_scaled, y_test)
    print(f"Per-sample accuracy (train): {train_acc:.4f}", file=sys.stderr)
    print(f"Per-sample accuracy (test): {test_acc:.4f}", file=sys.stderr)
    
    # Frame selection accuracy
    train_frame_acc = evaluate_frame_selection(model, scaler, train_reads[:2000], train_frames[:2000])
    test_frame_acc = evaluate_frame_selection(model, scaler, test_reads, test_frames)
    print(f"Frame selection accuracy (train): {train_frame_acc:.4f}", file=sys.stderr)
    print(f"Frame selection accuracy (test): {test_frame_acc:.4f}", file=sys.stderr)
    
    # Validate on real data
    if args.validation:
        print("\n=== Validation on real aDNA ===", file=sys.stderr)
        for val_path in args.validation:
            val_data = load_ground_truth(val_path)
            if val_data:
                val_reads = [d[0] for d in val_data]
                val_frames = [d[1] for d in val_data]
                val_acc = evaluate_frame_selection(model, scaler, val_reads, val_frames)
                print(f"  {os.path.basename(val_path)}: {val_acc:.4f}", file=sys.stderr)
    
    # Export
    export_to_cpp(model, scaler, f"{args.output}.hpp")
    
    print("\nDone!", file=sys.stderr)


if __name__ == '__main__':
    main()
