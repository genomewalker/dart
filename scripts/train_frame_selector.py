#!/usr/bin/env python3
"""
Train a small neural network for frame selection.

This NN takes features from all 6 frames and predicts which frame is correct.
Features per frame:
- AA composition (20 values)
- Dipeptide frequencies (top 50 dipeptides)
- Stop codon count
- GC3 content
- Sequence length

Input: 6 frames x ~75 features = ~450 features
Output: 6-way softmax (which frame is correct)
"""

import gzip
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from collections import defaultdict
import sys
import os

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
TOP_DIPEPTIDES = ['LL', 'AA', 'AL', 'LA', 'LS', 'VL', 'LG', 'EL', 'GL', 'SL',
                  'LE', 'LV', 'AG', 'VA', 'AV', 'SS', 'GG', 'LK', 'EA', 'EE',
                  'GA', 'LR', 'AE', 'RL', 'TL', 'DL', 'LD', 'KL', 'IL', 'LT']

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

def translate(seq, frame):
    protein = ''
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        protein += CODON_TABLE.get(codon, 'X')
    return protein

def extract_features(seq):
    """Extract features for all 6 frames."""
    features = []

    # Forward and reverse complement
    seqs = [seq, reverse_complement(seq)]

    for strand_idx, s in enumerate(seqs):
        for frame in range(3):
            frame_features = []

            # Translate
            protein = translate(s, frame)

            if len(protein) < 3:
                # Pad with zeros for short sequences
                frame_features = [0.0] * 53
                features.extend(frame_features)
                continue

            # AA composition (20 features)
            aa_counts = defaultdict(int)
            for aa in protein:
                if aa in AA_ORDER:
                    aa_counts[aa] += 1
            total_aa = sum(aa_counts.values()) or 1
            for aa in AA_ORDER:
                frame_features.append(aa_counts[aa] / total_aa)

            # Dipeptide frequencies (30 features)
            dp_counts = defaultdict(int)
            for i in range(len(protein) - 1):
                dp = protein[i:i+2]
                if all(c in AA_ORDER for c in dp):
                    dp_counts[dp] += 1
            total_dp = sum(dp_counts.values()) or 1
            for dp in TOP_DIPEPTIDES:
                frame_features.append(dp_counts[dp] / total_dp)

            # Stop codon count (1 feature)
            internal_stops = protein[:-1].count('*') if len(protein) > 1 else 0
            frame_features.append(internal_stops / (len(protein) or 1))

            # GC3 content (1 feature)
            gc3_count = sum(1 for i in range(frame + 2, len(s), 3) if s[i].upper() in 'GC')
            total_pos3 = len(range(frame + 2, len(s), 3)) or 1
            frame_features.append(gc3_count / total_pos3)

            # Length (1 feature, normalized)
            frame_features.append(min(1.0, len(protein) / 50.0))

            features.extend(frame_features)

    return np.array(features, dtype=np.float32)

class FrameSelectorNN(nn.Module):
    """Small MLP for frame selection."""
    def __init__(self, input_size=318, hidden_size=64):
        super().__init__()
        self.model = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_size, hidden_size // 2),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_size // 2, 6),  # 6 frames
        )

    def forward(self, x):
        return self.model(x)

def load_training_data(tsv_path, fastq_path, max_samples=100000):
    """Load training data from ground truth TSV and FASTQ."""
    print(f"Loading training data from {tsv_path}...", file=sys.stderr)

    # Read ground truth
    ground_truth = {}
    with gzip.open(tsv_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        strand_idx = header.index('Strand_read')
        read_name_idx = header.index('read_name')

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(strand_idx, read_name_idx):
                continue

            read_name = parts[read_name_idx]
            strand = parts[strand_idx]

            # Extract frame from read name (format: ...:frame_X_Y_Z)
            # The frame info is in the read name pattern
            # For now, we'll use the strand info
            ground_truth[read_name.split(':')[0]] = strand

    print(f"Loaded {len(ground_truth)} ground truth entries", file=sys.stderr)

    # Read sequences and create training data
    X = []
    y = []

    with gzip.open(fastq_path, 'rt') as f:
        count = 0
        for line in f:
            if count >= max_samples * 4:
                break

            if count % 4 == 0:
                read_id = line.strip().split()[0][1:]  # Remove @
                base_id = read_id.split(':')[0]
            elif count % 4 == 1:
                seq = line.strip()

                if base_id in ground_truth and len(seq) >= 30:
                    strand = ground_truth[base_id]

                    # Determine correct frame
                    # For simplicity: forward = frames 0,1,2; reverse = frames 3,4,5
                    # The exact frame within strand needs more info from ground truth
                    correct_frame = 0 if strand == '+' else 3  # Default to frame 0

                    features = extract_features(seq)
                    X.append(features)
                    y.append(correct_frame)

            count += 1

    print(f"Created {len(X)} training samples", file=sys.stderr)
    return np.array(X), np.array(y)

def train_model(X, y, epochs=50, batch_size=256):
    """Train the frame selector model."""
    # Convert to tensors
    X_tensor = torch.tensor(X, dtype=torch.float32)
    y_tensor = torch.tensor(y, dtype=torch.long)

    # Split train/val
    n = len(X)
    perm = np.random.permutation(n)
    train_idx = perm[:int(0.8 * n)]
    val_idx = perm[int(0.8 * n):]

    X_train, y_train = X_tensor[train_idx], y_tensor[train_idx]
    X_val, y_val = X_tensor[val_idx], y_tensor[val_idx]

    # Create model
    model = FrameSelectorNN(input_size=X.shape[1])
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    print(f"Training model with {len(X_train)} samples...", file=sys.stderr)

    for epoch in range(epochs):
        model.train()
        total_loss = 0

        # Mini-batch training
        perm = torch.randperm(len(X_train))
        for i in range(0, len(X_train), batch_size):
            idx = perm[i:i+batch_size]
            batch_x, batch_y = X_train[idx], y_train[idx]

            optimizer.zero_grad()
            outputs = model(batch_x)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        # Validation
        model.eval()
        with torch.no_grad():
            val_outputs = model(X_val)
            val_preds = val_outputs.argmax(dim=1)
            val_acc = (val_preds == y_val).float().mean().item()

            # Also check strand accuracy (frames 0-2 vs 3-5)
            pred_strand = (val_preds >= 3).long()
            true_strand = (y_val >= 3).long()
            strand_acc = (pred_strand == true_strand).float().mean().item()

        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}: loss={total_loss:.4f}, val_acc={val_acc:.4f}, strand_acc={strand_acc:.4f}",
                  file=sys.stderr)

    return model

def export_model_weights(model, output_path):
    """Export model weights to a format that can be loaded in C++."""
    weights = {}
    for name, param in model.named_parameters():
        weights[name] = param.detach().numpy()

    np.savez(output_path, **weights)
    print(f"Saved model weights to {output_path}", file=sys.stderr)

    # Also save as simple text format for C++ loading
    txt_path = output_path.replace('.npz', '.txt')
    with open(txt_path, 'w') as f:
        for name, arr in weights.items():
            f.write(f"# {name} shape: {arr.shape}\n")
            if arr.ndim == 1:
                f.write(' '.join(f'{x:.6f}' for x in arr) + '\n')
            else:
                for row in arr:
                    f.write(' '.join(f'{x:.6f}' for x in row) + '\n')
    print(f"Saved weights as text to {txt_path}", file=sys.stderr)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Train frame selector NN')
    parser.add_argument('--tsv', required=True, help='Ground truth TSV.gz file')
    parser.add_argument('--fastq', required=True, help='FASTQ.gz file')
    parser.add_argument('--output', default='frame_selector_model.npz', help='Output model path')
    parser.add_argument('--epochs', type=int, default=50, help='Training epochs')
    parser.add_argument('--max-samples', type=int, default=100000, help='Max training samples')
    args = parser.parse_args()

    # Load data
    X, y = load_training_data(args.tsv, args.fastq, args.max_samples)

    if len(X) < 100:
        print("Not enough training data!", file=sys.stderr)
        sys.exit(1)

    # Train
    model = train_model(X, y, epochs=args.epochs)

    # Export
    export_model_weights(model, args.output)

    # Save full model
    torch.save(model.state_dict(), args.output.replace('.npz', '.pt'))
    print(f"Saved PyTorch model to {args.output.replace('.npz', '.pt')}", file=sys.stderr)
