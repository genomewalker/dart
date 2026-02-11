#!/usr/bin/env python3
"""
Train frame classifier on real aMGSIM ground truth data.

This uses labeled ancient DNA data where we know the true reading frame.
The model learns to distinguish correct from incorrect frames on damaged sequences.

Usage:
    python train_frame_selector.py --train data1.tsv.gz data2.tsv.gz --test test.tsv.gz --output model
"""

import sys
import os
import gzip
import argparse
import numpy as np
from collections import defaultdict

try:
    from sklearn.neural_network import MLPClassifier
    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("Error: sklearn required", file=sys.stderr)
    sys.exit(1)

# Constants
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
NUM_FEATURES = 30


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))


def translate(seq, frame):
    protein = ''
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        protein += CODON_TABLE.get(codon, 'X')
    return protein


def extract_features(seq, frame):
    """Extract 30 features for frame scoring."""
    seq = seq.upper()
    features = []
    
    # Get codons
    codons = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if all(c in 'ACGT' for c in codon):
            codons.append(codon)
    
    n_codons = len(codons)
    if n_codons < 2:
        return [0.0] * NUM_FEATURES
    
    # 1. Internal stop count (1)
    internal_stops = sum(1 for c in codons[:-1] if c in STOP_CODONS)
    features.append(min(internal_stops, 5) / 5.0)
    
    # 2. Has any internal stop (1)
    features.append(1.0 if internal_stops > 0 else 0.0)
    
    # 3. First stop position normalized (1)
    first_stop = next((i for i, c in enumerate(codons[:-1]) if c in STOP_CODONS), -1)
    features.append(first_stop / n_codons if first_stop >= 0 else 1.0)
    
    # 4. GC content per codon position (3)
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
    
    # 5. Nucleotide freq per position (12)
    nt_freq = [[0]*4 for _ in range(3)]
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
    
    # 6. RNY compliance (1)
    rny = sum(1 for i in range(frame, len(seq)-2, 3) 
              if seq[i] in 'AG' and seq[i+2] in 'CT')
    features.append(rny / n_codons)
    
    # 7. Dipeptide entropy (1)
    protein = translate(seq, frame)
    if len(protein) >= 2:
        dp_counts = defaultdict(int)
        for i in range(len(protein)-1):
            dp_counts[protein[i:i+2]] += 1
        total = sum(dp_counts.values())
        entropy = -sum((c/total) * np.log2(c/total) for c in dp_counts.values() if c > 0)
        features.append(entropy / 8.64)  # Normalized
    else:
        features.append(0.5)
    
    # 8. Length normalized (1)
    features.append(min(1.0, n_codons / 50.0))
    
    # 9. Has start codon (1)
    features.append(1.0 if codons and codons[0] == 'ATG' else 0.0)
    
    # 10. Has terminal stop (1)
    features.append(1.0 if codons and codons[-1] in STOP_CODONS else 0.0)
    
    # 11. Pyrimidine at wobble position (1)
    pyr_wobble = sum(1 for i in range(frame+2, len(seq), 3) if seq[i] in 'CT')
    features.append(pyr_wobble / max(1, n_codons))
    
    # 12-14. Dinucleotide features (3)
    seq_len = len(seq)
    features.append(seq.count('CG') / max(1, seq_len//2) * 4)
    features.append(seq.count('TA') / max(1, seq_len//2) * 4)
    features.append((seq.count('AA') + seq.count('TT')) / max(1, seq_len//2) * 2)
    
    # 15. GC content overall (1)
    features.append((seq.count('G') + seq.count('C')) / len(seq))
    
    # 16. Purine at position 1 (1)
    pur_pos1 = sum(1 for i in range(frame, len(seq), 3) if seq[i] in 'AG')
    features.append(pur_pos1 / max(1, n_codons))
    
    # 17. A/T content at position 2 (1) - wobble position tends to be A/T rich
    at_pos2 = sum(1 for i in range(frame+1, len(seq), 3) if seq[i] in 'AT')
    features.append(at_pos2 / max(1, n_codons))
    
    assert len(features) == NUM_FEATURES, f"Expected {NUM_FEATURES}, got {len(features)}"
    return features


def load_ground_truth(tsv_paths, max_per_file=50000):
    """Load ground truth from aMGSIM TSV files."""
    all_data = []
    
    for tsv_path in tsv_paths:
        print(f"Loading {os.path.basename(tsv_path)}...", file=sys.stderr)
        
        opener = gzip.open if tsv_path.endswith('.gz') else open
        mode = 'rt' if tsv_path.endswith('.gz') else 'r'
        count = 0
        
        with opener(tsv_path, mode) as f:
            header = f.readline().strip().split('\t')
            cols = {name: i for i, name in enumerate(header)}
            
            required = ['damaged_seq', 'damaged_seq_inframe_nt', 'Strand_read', 'Strand_gene']
            if not all(c in cols for c in required):
                print(f"  Skipping - missing columns", file=sys.stderr)
                continue
            
            for line in f:
                if count >= max_per_file:
                    break
                
                parts = line.strip().split('\t')
                if len(parts) <= max(cols.values()):
                    continue
                
                damaged_seq = parts[cols['damaged_seq']].upper()
                inframe_nt = parts[cols['damaged_seq_inframe_nt']].upper()
                strand_read = parts[cols['Strand_read']]
                strand_gene = parts[cols['Strand_gene']]
                
                if len(damaged_seq) < 30:
                    continue
                
                # Determine orientation
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
                    all_data.append((damaged_seq, true_frame))
                    count += 1
        
        print(f"  Loaded {count} samples", file=sys.stderr)
    
    return all_data


def create_dataset(data):
    """Create X, y arrays for training."""
    X = []
    y = []
    
    for seq, true_frame in data:
        for frame in range(3):
            features = extract_features(seq, frame)
            X.append(features)
            y.append(1 if frame == true_frame else 0)
    
    return np.array(X, dtype=np.float32), np.array(y, dtype=np.int32)


def evaluate_frame_accuracy(model, scaler, data):
    """Evaluate frame selection accuracy."""
    correct = 0
    for seq, true_frame in data:
        features = np.array([extract_features(seq, f) for f in range(3)])
        features_scaled = scaler.transform(features)
        
        if hasattr(model, 'predict_proba'):
            probs = model.predict_proba(features_scaled)[:, 1]
        else:
            probs = model.decision_function(features_scaled)
        
        predicted = np.argmax(probs)
        if predicted == true_frame:
            correct += 1
    
    return correct / len(data)


def export_mlp_to_cpp(model, scaler, output_path):
    """Export MLP to C++ header."""
    W1 = model.coefs_[0]
    b1 = model.intercepts_[0]
    W2 = model.coefs_[1]
    b2 = model.intercepts_[1]
    
    input_size = W1.shape[0]
    hidden_size = W1.shape[1]
    output_size = W2.shape[1]  # Usually 1 for binary classification
    
    with open(output_path, 'w') as f:
        f.write("#pragma once\n\n")
        f.write("/**\n")
        f.write(" * Frame classifier trained on aMGSIM ground truth.\n")
        f.write(" * Generated by train_frame_selector.py\n")
        f.write(f" * Architecture: {input_size} -> {hidden_size} (ReLU) -> {output_size} (sigmoid)\n")
        f.write(" */\n\n")
        f.write("#include <array>\n")
        f.write("#include <cmath>\n")
        f.write("#include <algorithm>\n\n")
        f.write("namespace agp {\n")
        f.write("namespace frame_nn {\n\n")
        
        f.write(f"constexpr int INPUT_SIZE = {input_size};\n")
        f.write(f"constexpr int HIDDEN_SIZE = {hidden_size};\n\n")
        
        # Scaler
        f.write("constexpr std::array<float, INPUT_SIZE> INPUT_MEAN = {{\n    ")
        f.write(",\n    ".join(f"{x:.8f}f" for x in scaler.mean_))
        f.write("\n}};\n\n")
        
        f.write("constexpr std::array<float, INPUT_SIZE> INPUT_STD = {{\n    ")
        f.write(",\n    ".join(f"{max(x, 1e-8):.8f}f" for x in scaler.scale_))
        f.write("\n}};\n\n")
        
        # W1
        f.write("constexpr std::array<std::array<float, HIDDEN_SIZE>, INPUT_SIZE> W1 = {{\n")
        for i, row in enumerate(W1):
            f.write("    {{" + ", ".join(f"{x:.8f}f" for x in row) + "}}")
            f.write(",\n" if i < len(W1)-1 else "\n")
        f.write("}};\n\n")
        
        f.write("constexpr std::array<float, HIDDEN_SIZE> B1 = {{\n    ")
        f.write(", ".join(f"{x:.8f}f" for x in b1))
        f.write("\n}};\n\n")
        
        # W2 - flatten if output_size is 1
        f.write("constexpr std::array<float, HIDDEN_SIZE> W2 = {{\n    ")
        f.write(", ".join(f"{x:.8f}f" for x in W2.flatten()))
        f.write("\n}};\n\n")
        
        f.write(f"constexpr float B2 = {b2[0]:.8f}f;\n\n")
        
        # Forward function (binary sigmoid output)
        f.write("inline float score_frame(const std::array<float, INPUT_SIZE>& features) {\n")
        f.write("    std::array<float, INPUT_SIZE> x;\n")
        f.write("    for (int i = 0; i < INPUT_SIZE; i++)\n")
        f.write("        x[i] = (features[i] - INPUT_MEAN[i]) / INPUT_STD[i];\n\n")
        f.write("    std::array<float, HIDDEN_SIZE> h;\n")
        f.write("    for (int j = 0; j < HIDDEN_SIZE; j++) {\n")
        f.write("        float sum = B1[j];\n")
        f.write("        for (int i = 0; i < INPUT_SIZE; i++)\n")
        f.write("            sum += x[i] * W1[i][j];\n")
        f.write("        h[j] = std::max(0.0f, sum);\n")
        f.write("    }\n\n")
        f.write("    // Binary output with sigmoid\n")
        f.write("    float z = B2;\n")
        f.write("    for (int j = 0; j < HIDDEN_SIZE; j++)\n")
        f.write("        z += h[j] * W2[j];\n")
        f.write("    return 1.0f / (1.0f + std::exp(-std::clamp(z, -20.0f, 20.0f)));\n")
        f.write("}\n\n")
        
        f.write("} // namespace frame_nn\n")
        f.write("} // namespace agp\n")


def export_logistic_to_cpp(model, scaler, output_path):
    """Export logistic regression to C++ header."""
    coef = model.coef_[0]
    intercept = model.intercept_[0]
    
    with open(output_path, 'w') as f:
        f.write("#pragma once\n\n")
        f.write("/**\n")
        f.write(" * Frame classifier (logistic regression) trained on aMGSIM.\n")
        f.write(" * Generated by train_frame_selector.py\n")
        f.write(" */\n\n")
        f.write("#include <array>\n")
        f.write("#include <cmath>\n\n")
        f.write("namespace agp {\n")
        f.write("namespace frame_nn {\n\n")
        
        f.write(f"constexpr int INPUT_SIZE = {len(coef)};\n\n")
        
        f.write("constexpr std::array<float, INPUT_SIZE> INPUT_MEAN = {{\n    ")
        f.write(",\n    ".join(f"{x:.8f}f" for x in scaler.mean_))
        f.write("\n}};\n\n")
        
        f.write("constexpr std::array<float, INPUT_SIZE> INPUT_STD = {{\n    ")
        f.write(",\n    ".join(f"{max(x, 1e-8):.8f}f" for x in scaler.scale_))
        f.write("\n}};\n\n")
        
        f.write("constexpr std::array<float, INPUT_SIZE> WEIGHTS = {{\n    ")
        f.write(",\n    ".join(f"{x:.8f}f" for x in coef))
        f.write("\n}};\n\n")
        
        f.write(f"constexpr float INTERCEPT = {intercept:.8f}f;\n\n")
        
        f.write("inline float score_frame(const std::array<float, INPUT_SIZE>& features) {\n")
        f.write("    float z = INTERCEPT;\n")
        f.write("    for (int i = 0; i < INPUT_SIZE; i++) {\n")
        f.write("        float x = (features[i] - INPUT_MEAN[i]) / INPUT_STD[i];\n")
        f.write("        z += WEIGHTS[i] * x;\n")
        f.write("    }\n")
        f.write("    return 1.0f / (1.0f + std::exp(-std::clamp(z, -20.0f, 20.0f)));\n")
        f.write("}\n\n")
        
        f.write("} // namespace frame_nn\n")
        f.write("} // namespace agp\n")


def main():
    parser = argparse.ArgumentParser(description='Train frame classifier on aMGSIM data')
    parser.add_argument('--train', nargs='+', required=True, help='Training TSV files')
    parser.add_argument('--test', nargs='*', help='Test TSV files (optional)')
    parser.add_argument('--output', default='frame_classifier', help='Output base name')
    parser.add_argument('--model', choices=['mlp', 'logistic', 'gbdt'], default='mlp')
    parser.add_argument('--hidden-size', type=int, default=32)
    parser.add_argument('--max-samples', type=int, default=100000)
    args = parser.parse_args()
    
    # Load training data
    print("=== Loading Training Data ===", file=sys.stderr)
    train_data = load_ground_truth(args.train, max_per_file=args.max_samples)
    print(f"Total training samples: {len(train_data)}", file=sys.stderr)
    
    if len(train_data) < 1000:
        print("Error: Not enough training data", file=sys.stderr)
        sys.exit(1)
    
    # Split for validation
    np.random.shuffle(train_data)
    val_size = min(10000, len(train_data) // 10)
    val_data = train_data[-val_size:]
    train_data = train_data[:-val_size]
    
    # Create datasets
    print("\n=== Creating Datasets ===", file=sys.stderr)
    X_train, y_train = create_dataset(train_data)
    X_val, y_val = create_dataset(val_data)
    
    print(f"Training: {len(X_train)} samples", file=sys.stderr)
    print(f"Validation: {len(X_val)} samples", file=sys.stderr)
    
    # Scale
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_val_scaled = scaler.transform(X_val)
    
    # Train
    print(f"\n=== Training {args.model.upper()} ===", file=sys.stderr)
    
    if args.model == 'mlp':
        model = MLPClassifier(
            hidden_layer_sizes=(args.hidden_size,),
            activation='relu',
            max_iter=300,
            early_stopping=True,
            validation_fraction=0.1,
            random_state=42,
            verbose=True
        )
    elif args.model == 'logistic':
        model = LogisticRegression(max_iter=1000, random_state=42, verbose=1)
    else:
        model = GradientBoostingClassifier(
            n_estimators=100,
            max_depth=5,
            random_state=42,
            verbose=1
        )
    
    model.fit(X_train_scaled, y_train)
    
    # Evaluate
    print("\n=== Results ===", file=sys.stderr)
    train_acc = model.score(X_train_scaled, y_train)
    val_acc = model.score(X_val_scaled, y_val)
    print(f"Per-sample accuracy (train): {train_acc:.4f}", file=sys.stderr)
    print(f"Per-sample accuracy (val): {val_acc:.4f}", file=sys.stderr)
    
    # Frame selection accuracy
    train_frame_acc = evaluate_frame_accuracy(model, scaler, train_data[:5000])
    val_frame_acc = evaluate_frame_accuracy(model, scaler, val_data)
    print(f"Frame selection (train): {train_frame_acc:.4f}", file=sys.stderr)
    print(f"Frame selection (val): {val_frame_acc:.4f}", file=sys.stderr)
    
    # Test on held-out files
    if args.test:
        print("\n=== Test Set Results ===", file=sys.stderr)
        test_data = load_ground_truth(args.test, max_per_file=10000)
        if test_data:
            test_frame_acc = evaluate_frame_accuracy(model, scaler, test_data)
            print(f"Frame selection (test): {test_frame_acc:.4f}", file=sys.stderr)
    
    # Export
    print(f"\n=== Exporting ===", file=sys.stderr)
    if args.model == 'mlp':
        export_mlp_to_cpp(model, scaler, f"{args.output}.hpp")
    elif args.model == 'logistic':
        export_logistic_to_cpp(model, scaler, f"{args.output}.hpp")
    else:
        print("Note: GBDT export not supported, saving sklearn model", file=sys.stderr)
        import pickle
        with open(f"{args.output}.pkl", 'wb') as f:
            pickle.dump({'model': model, 'scaler': scaler}, f)
    
    print("\nDone!", file=sys.stderr)


if __name__ == '__main__':
    main()
