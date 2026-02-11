#!/usr/bin/env python3
"""
Train a gradient boosting classifier for reading frame selection.

This script trains on ground truth data from aMGSIM to learn which
features best discriminate between correct and incorrect reading frames.

Input: Ground truth TSV files from benchmark_data/
Output: 
  - Trained model (joblib format)
  - Feature importance analysis
  - C++ header with model weights (for simple models)

Features extracted per frame:
1. Codon usage bias (64 codons → compressed to PCA or direct)
2. AA composition (20 values)
3. Dipeptide frequencies (top 30)
4. Stop codon count/position
5. GC content at each codon position (GC1, GC2, GC3)
6. Nucleotide periodicity scores
7. Forbidden pattern count
8. RNY rule compliance
9. Hexamer rarity score
10. Length

Relative features (comparing between frames):
- Score differences between this frame and alternatives
- Rank of this frame among all 3
"""

import sys
import os
import gzip
import argparse
import numpy as np
from collections import defaultdict
import pickle

# Try to import sklearn, fallback to simple logistic regression if not available
try:
    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split, cross_val_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import classification_report, confusion_matrix
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("Warning: scikit-learn not available, using simple training", file=sys.stderr)

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

AA_ORDER = 'ACDEFGHIKLMNPQRSTVWY'
STOP_CODONS = {'TAA', 'TAG', 'TGA'}
RARE_CODONS = {'AGG', 'AGA', 'CGA', 'CGG', 'CTA', 'ATA'}

def reverse_complement(seq):
    """Reverse complement a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

def translate(seq, frame):
    """Translate DNA to protein."""
    protein = ''
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        protein += CODON_TABLE.get(codon, 'X')
    return protein

def count_internal_stops(protein):
    """Count internal stop codons (excluding terminal)."""
    if len(protein) <= 1:
        return 0
    return protein[:-1].count('*')

def calc_gc_by_position(seq, frame):
    """Calculate GC content at each codon position."""
    gc = [0, 0, 0]
    total = [0, 0, 0]
    
    for i in range(frame, len(seq)):
        pos = (i - frame) % 3
        c = seq[i].upper()
        if c in 'ACGT':
            total[pos] += 1
            if c in 'GC':
                gc[pos] += 1
    
    return [gc[p] / max(1, total[p]) for p in range(3)]

def count_forbidden_patterns(seq, frame):
    """Count hexamers starting with stop codons at codon boundaries."""
    count = 0
    for i in range(frame, len(seq) - 5, 3):
        if seq[i:i+3].upper() in STOP_CODONS:
            count += 1
        if seq[i+3:i+6].upper() in STOP_CODONS:
            count += 1
    return count

def calc_rny_compliance(seq, frame):
    """Calculate RNY rule compliance (R-N-Y pattern)."""
    rny = 0
    total = 0
    for i in range(frame, len(seq) - 2, 3):
        c0, c2 = seq[i].upper(), seq[i+2].upper()
        total += 1
        if c0 in 'AG' and c2 in 'CT':
            rny += 1
    return rny / max(1, total)

def calc_nucleotide_periodicity(seq, frame):
    """Calculate nucleotide periodicity features."""
    nt_counts = [[0, 0, 0, 0] for _ in range(3)]  # [pos][A,C,G,T]
    nt_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    for i, c in enumerate(seq.upper()):
        if c not in nt_map:
            continue
        pos = (i - frame + 300) % 3
        nt_counts[pos][nt_map[c]] += 1
    
    features = []
    for pos in range(3):
        total = sum(nt_counts[pos])
        if total < 3:
            features.extend([0.25, 0.25, 0.25, 0.25])
        else:
            features.extend([c / total for c in nt_counts[pos]])
    
    return features  # 12 features total

def extract_frame_features(seq, frame):
    """Extract features for a single frame."""
    features = []
    
    # Translate
    protein = translate(seq, frame)
    prot_len = len(protein)
    
    if prot_len < 2:
        return [0.0] * 50  # Return zeros for short sequences
    
    # 1. AA composition (20 features)
    aa_counts = defaultdict(int)
    for aa in protein:
        if aa in AA_ORDER:
            aa_counts[aa] += 1
    total_aa = sum(aa_counts.values()) or 1
    for aa in AA_ORDER:
        features.append(aa_counts[aa] / total_aa)
    
    # 2. Internal stop count (1 feature)
    stops = count_internal_stops(protein)
    features.append(stops)
    
    # 3. Stop position if any (1 feature) - normalized position of first stop
    first_stop = protein.find('*')
    features.append(first_stop / prot_len if first_stop >= 0 and first_stop < prot_len - 1 else 1.0)
    
    # 4. GC by position (3 features)
    gc_pos = calc_gc_by_position(seq, frame)
    features.extend(gc_pos)
    
    # 5. Forbidden patterns (1 feature)
    forbidden = count_forbidden_patterns(seq, frame)
    features.append(forbidden)
    
    # 6. RNY compliance (1 feature)
    rny = calc_rny_compliance(seq, frame)
    features.append(rny)
    
    # 7. Nucleotide periodicity (12 features)
    periodicity = calc_nucleotide_periodicity(seq, frame)
    features.extend(periodicity)
    
    # 8. Rare codon count (1 feature)
    rare_count = 0
    for i in range(frame, len(seq) - 2, 3):
        if seq[i:i+3].upper() in RARE_CODONS:
            rare_count += 1
    features.append(rare_count / max(1, prot_len))
    
    # 9. Length features (2 features)
    features.append(min(1.0, prot_len / 50.0))  # Normalized length
    features.append(prot_len % 3 == 0)  # Divisible by 3
    
    # 10. Terminal features (4 features)
    # 5' terminal T ratio (first 5 positions)
    t_5prime = sum(1 for i in range(min(5, len(seq))) if seq[i].upper() == 'T') / 5
    features.append(t_5prime)
    # 3' terminal A ratio (last 5 positions)
    a_3prime = sum(1 for i in range(max(0, len(seq)-5), len(seq)) if seq[i].upper() == 'A') / 5
    features.append(a_3prime)
    # Start codon (ATG) present
    has_start = 1.0 if seq[frame:frame+3].upper() == 'ATG' else 0.0
    features.append(has_start)
    # Near-start codon (within first 2 codons)
    near_start = 0.0
    for offset in [0, 3]:
        if frame + offset + 3 <= len(seq) and seq[frame+offset:frame+offset+3].upper() == 'ATG':
            near_start = 1.0
            break
    features.append(near_start)
    
    return features

def extract_all_features(seq):
    """Extract features for all 6 frames with relative comparisons."""
    # Get raw features for each frame
    raw_features = []
    
    # Forward strand
    for frame in range(3):
        raw_features.append(extract_frame_features(seq, frame))
    
    # Reverse strand
    rc_seq = reverse_complement(seq)
    for frame in range(3):
        raw_features.append(extract_frame_features(rc_seq, frame))
    
    # Convert to numpy for easier manipulation
    raw_features = np.array(raw_features, dtype=np.float32)
    
    # Create final feature vector for each frame
    # Include both raw features and relative features
    all_features = []
    
    for i in range(6):
        frame_features = list(raw_features[i])
        
        # Add relative features: difference from mean of other frames
        other_mean = np.mean(raw_features[[j for j in range(6) if j != i]], axis=0)
        frame_features.extend(raw_features[i] - other_mean)
        
        # Add rank features for key metrics
        # Internal stops rank (0 = fewest)
        stop_idx = 20  # Position of stop count in features
        stop_counts = raw_features[:, stop_idx]
        frame_features.append(np.sum(stop_counts < stop_counts[i]))  # Rank
        
        # Forbidden patterns rank
        forbidden_idx = 24  # Position of forbidden count
        forbidden_counts = raw_features[:, forbidden_idx]
        frame_features.append(np.sum(forbidden_counts < forbidden_counts[i]))
        
        all_features.append(frame_features)
    
    return np.array(all_features, dtype=np.float32)

def determine_true_frame(damaged_seq, inframe_nt, strand_read, strand_gene):
    """
    Determine the true reading frame and strand.
    
    Returns: (frame, is_forward) where frame is 0, 1, or 2
    """
    # Determine if we should use forward or reverse complement
    is_forward = (strand_read == strand_gene)
    
    if not is_forward:
        damaged_seq = reverse_complement(damaged_seq)
    
    damaged_seq = damaged_seq.upper()
    inframe_nt = inframe_nt.upper() if inframe_nt else ""
    
    if not inframe_nt or len(inframe_nt) < 6:
        return None, None
    
    # Find where inframe_nt aligns in damaged_seq
    # The inframe_nt should be a substring at a codon boundary
    for frame in range(3):
        # Check if damaged_seq starting at frame contains inframe_nt
        test_seq = damaged_seq[frame:]
        if inframe_nt in test_seq:
            # Verify it's at a codon boundary
            pos = test_seq.find(inframe_nt)
            if pos % 3 == 0:
                return frame, is_forward
    
    # Fallback: try to find by translating
    return None, None

def load_ground_truth(tsv_path, max_samples=None):
    """
    Load ground truth data from aMGSIM TSV file.
    
    Returns: List of (sequence, true_frame, is_forward) tuples
    """
    print(f"Loading ground truth from {tsv_path}...", file=sys.stderr)
    
    data = []
    skipped = 0
    
    opener = gzip.open if tsv_path.endswith('.gz') else open
    mode = 'rt' if tsv_path.endswith('.gz') else 'r'
    
    with opener(tsv_path, mode) as f:
        header = f.readline().strip().split('\t')
        
        # Find column indices
        cols = {name: i for i, name in enumerate(header)}
        
        required = ['damaged_seq', 'damaged_seq_inframe_nt', 'Strand_read', 'Strand_gene']
        for col in required:
            if col not in cols:
                print(f"Error: Required column '{col}' not found in TSV", file=sys.stderr)
                print(f"Available columns: {list(cols.keys())}", file=sys.stderr)
                sys.exit(1)
        
        for line_num, line in enumerate(f):
            if max_samples and len(data) >= max_samples:
                break
            
            parts = line.strip().split('\t')
            if len(parts) <= max(cols.values()):
                skipped += 1
                continue
            
            damaged_seq = parts[cols['damaged_seq']]
            inframe_nt = parts[cols['damaged_seq_inframe_nt']]
            strand_read = parts[cols['Strand_read']]
            strand_gene = parts[cols['Strand_gene']]
            
            if not damaged_seq or len(damaged_seq) < 30:
                skipped += 1
                continue
            
            frame, is_forward = determine_true_frame(damaged_seq, inframe_nt, strand_read, strand_gene)
            
            if frame is None:
                skipped += 1
                continue
            
            data.append((damaged_seq.upper(), frame, is_forward))
            
            if (len(data) % 10000) == 0:
                print(f"  Loaded {len(data)} samples...", file=sys.stderr)
    
    print(f"Loaded {len(data)} samples, skipped {skipped}", file=sys.stderr)
    return data

def create_training_data(ground_truth_data):
    """
    Create training dataset from ground truth.
    
    For each read, we create 6 examples (one per frame) with labels.
    """
    X = []
    y = []
    
    for seq, true_frame, is_forward in ground_truth_data:
        features = extract_all_features(seq)
        
        # True label: frame 0-2 for forward, 3-5 for reverse
        true_label = true_frame if is_forward else true_frame + 3
        
        # For now, we'll train a binary classifier for each frame position
        # Or we can train a single 6-way classifier
        X.append(features[true_label])  # Features for the correct frame
        y.append(1)  # Positive label
        
        # Add negative examples from wrong frames
        for i in range(6):
            if i != true_label:
                X.append(features[i])
                y.append(0)
    
    return np.array(X, dtype=np.float32), np.array(y, dtype=np.int32)

def train_binary_classifier(X, y, model_type='gradient_boosting'):
    """Train a binary classifier to score frames."""
    
    print(f"Training {model_type} classifier...", file=sys.stderr)
    print(f"  Samples: {len(X)} ({sum(y)} positive, {len(y)-sum(y)} negative)", file=sys.stderr)
    print(f"  Features: {X.shape[1]}", file=sys.stderr)
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train model
    if model_type == 'gradient_boosting':
        model = GradientBoostingClassifier(
            n_estimators=100,
            max_depth=5,
            learning_rate=0.1,
            random_state=42
        )
    elif model_type == 'random_forest':
        model = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1
        )
    else:
        model = LogisticRegression(
            max_iter=1000,
            random_state=42
        )
    
    model.fit(X_train_scaled, y_train)
    
    # Evaluate
    train_score = model.score(X_train_scaled, y_train)
    test_score = model.score(X_test_scaled, y_test)
    
    print(f"  Train accuracy: {train_score:.4f}", file=sys.stderr)
    print(f"  Test accuracy: {test_score:.4f}", file=sys.stderr)
    
    # Cross-validation
    cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=5)
    print(f"  CV accuracy: {cv_scores.mean():.4f} (+/- {cv_scores.std()*2:.4f})", file=sys.stderr)
    
    # Feature importance
    if hasattr(model, 'feature_importances_'):
        importances = model.feature_importances_
        top_features = np.argsort(importances)[-10:][::-1]
        print("\n  Top 10 features:", file=sys.stderr)
        for idx in top_features:
            print(f"    Feature {idx}: {importances[idx]:.4f}", file=sys.stderr)
    
    return model, scaler

def export_to_cpp_header(model, scaler, feature_names, output_path):
    """Export model to C++ header file."""
    
    print(f"\nExporting model to {output_path}...", file=sys.stderr)
    
    with open(output_path, 'w') as f:
        f.write("#pragma once\n\n")
        f.write("// Auto-generated frame classifier model\n")
        f.write("// Generated by train_frame_classifier.py\n\n")
        f.write("#include <array>\n")
        f.write("#include <cmath>\n\n")
        f.write("namespace agp {\n")
        f.write("namespace frame_ml {\n\n")
        
        # Export scaler parameters
        f.write("// Feature scaling parameters\n")
        f.write(f"constexpr int NUM_FEATURES = {len(scaler.mean_)};\n\n")
        
        f.write("constexpr std::array<float, NUM_FEATURES> FEATURE_MEAN = {{\n    ")
        f.write(",\n    ".join(f"{x:.8f}f" for x in scaler.mean_))
        f.write("\n}};\n\n")
        
        f.write("constexpr std::array<float, NUM_FEATURES> FEATURE_STD = {{\n    ")
        f.write(",\n    ".join(f"{max(x, 1e-8):.8f}f" for x in scaler.scale_))
        f.write("\n}};\n\n")
        
        # For logistic regression, export weights directly
        if hasattr(model, 'coef_'):
            coef = model.coef_.flatten()
            intercept = model.intercept_[0] if hasattr(model.intercept_, '__len__') else model.intercept_
            
            f.write("// Logistic regression weights\n")
            f.write("constexpr std::array<float, NUM_FEATURES> WEIGHTS = {{\n    ")
            f.write(",\n    ".join(f"{x:.8f}f" for x in coef))
            f.write("\n}};\n\n")
            
            f.write(f"constexpr float INTERCEPT = {intercept:.8f}f;\n\n")
            
            f.write("/**\n")
            f.write(" * Score a frame using the trained logistic regression model.\n")
            f.write(" * Returns probability that this is the correct frame.\n")
            f.write(" */\n")
            f.write("inline float score_frame_ml(const std::array<float, NUM_FEATURES>& features) {\n")
            f.write("    float z = INTERCEPT;\n")
            f.write("    for (int i = 0; i < NUM_FEATURES; i++) {\n")
            f.write("        float scaled = (features[i] - FEATURE_MEAN[i]) / FEATURE_STD[i];\n")
            f.write("        z += WEIGHTS[i] * scaled;\n")
            f.write("    }\n")
            f.write("    return 1.0f / (1.0f + std::exp(-z));\n")
            f.write("}\n\n")
        else:
            # For tree-based models, just export feature importance
            f.write("// This model type requires full sklearn export\n")
            f.write("// See the .pkl file for the full model\n\n")
        
        f.write("} // namespace frame_ml\n")
        f.write("} // namespace agp\n")
    
    print(f"  Exported to {output_path}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Train frame classifier')
    parser.add_argument('--tsv', required=True, nargs='+', help='Ground truth TSV.gz file(s)')
    parser.add_argument('--output', default='frame_classifier', help='Output base name')
    parser.add_argument('--max-samples', type=int, default=100000, help='Max training samples per file')
    parser.add_argument('--model-type', choices=['gradient_boosting', 'random_forest', 'logistic'],
                        default='logistic', help='Model type')
    args = parser.parse_args()
    
    if not HAS_SKLEARN:
        print("Error: scikit-learn is required for training", file=sys.stderr)
        sys.exit(1)
    
    # Load all ground truth files
    all_data = []
    for tsv_path in args.tsv:
        data = load_ground_truth(tsv_path, args.max_samples)
        all_data.extend(data)
    
    print(f"\nTotal samples: {len(all_data)}", file=sys.stderr)
    
    if len(all_data) < 1000:
        print("Error: Not enough training data", file=sys.stderr)
        sys.exit(1)
    
    # Create training data
    print("\nCreating training dataset...", file=sys.stderr)
    X, y = create_training_data(all_data)
    
    # Train model
    model, scaler = train_binary_classifier(X, y, args.model_type)
    
    # Save model
    model_path = f"{args.output}.pkl"
    with open(model_path, 'wb') as f:
        pickle.dump({'model': model, 'scaler': scaler}, f)
    print(f"\nSaved model to {model_path}", file=sys.stderr)
    
    # Export to C++ header (only for logistic regression)
    if args.model_type == 'logistic':
        header_path = f"{args.output}.hpp"
        export_to_cpp_header(model, scaler, None, header_path)
    
    print("\nDone!", file=sys.stderr)

if __name__ == '__main__':
    main()
