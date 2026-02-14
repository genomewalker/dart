#!/usr/bin/env python3
"""
Train a neural network for ancient DNA damage detection.

Combines:
1. Terminal nucleotide patterns
2. Codon-based mutation signals
3. Hexamer distributions (GTDB-derived)
4. Sample-level context from Pass 1

Usage:
    python train_damage_nn.py --data-dir /path/to/benchmark_data --output model.h

Author: AGP Development Team
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np

# Check for required packages
try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import Dataset, DataLoader
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score, precision_recall_curve, f1_score
except ImportError as e:
    print(f"Missing required package: {e}")
    print("Install with: pip install torch scikit-learn")
    sys.exit(1)


# Hexamer tables (Python subset; full table is in C++)


def encode_hexamer(seq: str) -> int:
    """Encode 6-mer as 12-bit integer."""
    base_map = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
    code = 0
    for c in seq[:6]:
        if c not in base_map:
            return -1
        code = (code << 2) | base_map[c]
    return code


# High-damage hexamers (TGA* patterns from C→T damage creating stops)
HIGH_DAMAGE_HEXAMERS = {
    "TGAACG",
    "TGAACC",
    "TGAAGC",
    "TGAAGG",
    "TGAATC",
    "TGAATG",
    "TGACAC",
    "TGACAG",
    "TGACCA",
    "TGACCC",
    "TGACCG",
    "TGACCT",
    "TGACGA",
    "TGACGC",
    "TGACGG",
    "TGACGT",
    "TGACTA",
    "TGACTC",
    "TGACTG",
    "TGACTT",
    "TGAGAA",
    "TGAGAC",
    "TGAGAG",
    "TGAGAT",
    "TGAGCA",
    "TGAGCC",
    "TGAGCG",
    "TGAGCT",
    "TGAGGA",
    "TGAGGC",
    "TGAGGG",
    "TGAGGT",
    "TGAGTA",
    "TGAGTC",
    "TGAGTG",
    "TGAGTT",
    "TGATAA",
    "TGATAC",
    "TGATAG",
    "TGATAT",
    "TGATCA",
    "TGATCC",
    "TGATCG",
    "TGATCT",
    "TGATGA",
    "TGATGC",
    "TGATGG",
    "TGATGT",
    "TGATTA",
    "TGATTC",
    "TGATTG",
    "TGATTT",
    # TAG* patterns
    "TAGAAA",
    "TAGAAC",
    "TAGAAG",
    "TAGAAT",
    "TAGACA",
    "TAGACC",
    # TAA* patterns
    "TAAAAA",
    "TAAAAC",
    "TAAAAG",
    "TAAAAT",
    "TAAACA",
    "TAAACC",
}


# Feature extraction


class DamageFeatureExtractor:
    """Extract features for damage detection from DNA sequences."""

    def __init__(self, sample_profile: Optional[Dict] = None):
        """
        Initialize feature extractor.

        Args:
            sample_profile: Sample-level damage statistics from Pass 1
        """
        self.sample_profile = sample_profile or {
            "damage_rate_5prime": 0.0,
            "damage_rate_3prime": 0.0,
            "decay_lambda_5": 0.69,
            "decay_lambda_3": 0.69,
            "baseline_tc": 0.5,
            "baseline_ag": 0.5,
            "is_double_stranded": True,
            "damage_level": 0,  # 0=minimal, 1=moderate, 2=high
        }

    def extract_features(self, seq: str) -> np.ndarray:
        """
        Extract all features for a sequence.

        Returns:
            numpy array of shape (56,) with all features
        """
        seq = seq.upper()
        features = []

        # Group 1: Terminal damage signals (12 features)
        features.extend(self._terminal_features(seq))

        # Group 2: Sample profile context (8 features)
        features.extend(self._sample_features())

        # Group 3: Codon-based features (18 features)
        features.extend(self._codon_features(seq))

        # Group 4: Hexamer features (12 features)
        features.extend(self._hexamer_features(seq))

        # Group 5: Sequence properties (6 features)
        features.extend(self._sequence_features(seq))

        return np.array(features, dtype=np.float32)

    def _terminal_features(self, seq: str) -> List[float]:
        """Extract terminal nucleotide pattern features."""
        features = []
        length = len(seq)

        # 5' end: T/(T+C) ratios at positions 0-4
        for i in range(5):
            if i < length:
                base = seq[i]
                if base == "T":
                    features.append(1.0)
                elif base == "C":
                    features.append(0.0)
                else:
                    features.append(0.5)  # Neither T nor C
            else:
                features.append(0.5)

        # 3' end: A/(A+G) ratios at positions -1 to -5
        for i in range(5):
            pos = length - 1 - i
            if pos >= 0:
                base = seq[pos]
                if base == "A":
                    features.append(1.0)
                elif base == "G":
                    features.append(0.0)
                else:
                    features.append(0.5)
            else:
                features.append(0.5)

        # Weighted sums with exponential decay (λ=0.69)
        lambda_decay = 0.69
        weighted_t_5 = sum(features[i] * np.exp(-lambda_decay * i) for i in range(5))
        weighted_a_3 = sum(
            features[5 + i] * np.exp(-lambda_decay * i) for i in range(5)
        )
        features.append(weighted_t_5 / 2.0)  # Normalize
        features.append(weighted_a_3 / 2.0)

        # Additional discriminative features

        # Count T's in first 3bp (strong damage signal)
        t_count_5 = sum(1 for i in range(min(3, length)) if seq[i] == "T")
        features.append(t_count_5 / 3.0)

        # Count A's in last 3bp
        a_count_3 = sum(1 for i in range(max(0, length - 3), length) if seq[i] == "A")
        features.append(a_count_3 / 3.0)

        # First base is T (very strong 5' damage indicator)
        features.append(1.0 if length > 0 and seq[0] == "T" else 0.0)

        # Last base is A (very strong 3' damage indicator)
        features.append(1.0 if length > 0 and seq[-1] == "A" else 0.0)

        # T run at 5' (TT or TTT pattern)
        t_run_5 = 0
        for i in range(min(3, length)):
            if seq[i] == "T":
                t_run_5 += 1
            else:
                break
        features.append(t_run_5 / 3.0)

        # A run at 3'
        a_run_3 = 0
        for i in range(length - 1, max(-1, length - 4), -1):
            if seq[i] == "A":
                a_run_3 += 1
            else:
                break
        features.append(a_run_3 / 3.0)

        return features  # 18 features

    def _sample_features(self) -> List[float]:
        """Extract sample-level context features."""
        p = self.sample_profile
        return [
            p.get("damage_rate_5prime", 0.0),
            p.get("damage_rate_3prime", 0.0),
            p.get("decay_lambda_5", 0.69),
            p.get("decay_lambda_3", 0.69),
            p.get("baseline_tc", 0.5),
            p.get("baseline_ag", 0.5),
            1.0 if p.get("is_double_stranded", True) else 0.0,
            p.get("damage_level", 0) / 2.0,  # Normalize to [0, 1]
        ]  # 8 features

    def _codon_features(self, seq: str) -> List[float]:
        """Extract codon-based mutation features for each frame."""
        features = []
        length = len(seq)

        for frame in range(3):
            # Count stop codons
            stop_counts = {"TAA": 0, "TAG": 0, "TGA": 0}
            total_codons = 0
            damaged_stops = 0

            for i in range(frame, length - 2, 3):
                codon = seq[i : i + 3]
                if len(codon) == 3:
                    total_codons += 1
                    if codon in stop_counts:
                        stop_counts[codon] += 1
                        # Check if stop is in damage-prone position
                        if i < 15 or i > length - 15:
                            damaged_stops += 1

            total_stops = sum(stop_counts.values())

            # Features for this frame
            features.append(total_stops / max(1, total_codons))  # Stop density
            features.append(stop_counts["TAA"] / max(1, total_codons))
            features.append(stop_counts["TAG"] / max(1, total_codons))
            features.append(stop_counts["TGA"] / max(1, total_codons))
            features.append(
                damaged_stops / max(1, total_stops) if total_stops > 0 else 0.0
            )

            # T at wobble position (3rd codon position)
            wobble_t = 0
            for i in range(frame + 2, length, 3):
                if seq[i] == "T":
                    wobble_t += 1
            features.append(wobble_t / max(1, total_codons))

        return features  # 18 features (6 per frame × 3 frames)

    def _hexamer_features(self, seq: str) -> List[float]:
        """Extract hexamer-based features."""
        features = []
        length = len(seq)

        # 5' region (first 12bp)
        hexamers_5 = [seq[i : i + 6] for i in range(min(7, length - 5))]
        high_damage_5 = sum(1 for h in hexamers_5 if h in HIGH_DAMAGE_HEXAMERS)
        features.append(len(hexamers_5) / 7.0 if hexamers_5 else 0.0)  # Coverage
        features.append(high_damage_5 / max(1, len(hexamers_5)))  # High damage ratio
        features.append(
            1.0 if any(h.startswith("TGA") for h in hexamers_5) else 0.0
        )  # TGA present

        # 3' region (last 12bp)
        start_3 = max(0, length - 12)
        hexamers_3 = [seq[i : i + 6] for i in range(start_3, length - 5)]
        high_damage_3 = sum(1 for h in hexamers_3 if h in HIGH_DAMAGE_HEXAMERS)
        features.append(len(hexamers_3) / 7.0 if hexamers_3 else 0.0)
        features.append(high_damage_3 / max(1, len(hexamers_3)))
        features.append(1.0 if any(h.startswith("TGA") for h in hexamers_3) else 0.0)

        # Middle region
        mid_start = length // 3
        mid_end = 2 * length // 3
        hexamers_mid = [
            seq[i : i + 6] for i in range(mid_start, min(mid_end, length - 5))
        ]
        high_damage_mid = sum(1 for h in hexamers_mid if h in HIGH_DAMAGE_HEXAMERS)
        features.append(len(hexamers_mid) / 20.0 if hexamers_mid else 0.0)  # Coverage
        features.append(high_damage_mid / max(1, len(hexamers_mid)))

        # Ratios
        total_5 = len(hexamers_5)
        total_mid = len(hexamers_mid)
        features.append(
            high_damage_5 / max(1, high_damage_mid)
            if high_damage_mid > 0
            else high_damage_5
        )
        features.append(
            high_damage_3 / max(1, high_damage_mid)
            if high_damage_mid > 0
            else high_damage_3
        )

        # GC content in hexamers
        gc_5 = sum(h.count("G") + h.count("C") for h in hexamers_5) / max(
            1, total_5 * 6
        )
        features.append(gc_5)

        # GC content in 3' hexamers
        total_3 = len(hexamers_3)
        gc_3 = sum(h.count("G") + h.count("C") for h in hexamers_3) / max(
            1, total_3 * 6
        )
        features.append(gc_3)

        return features  # 12 features

    def _sequence_features(self, seq: str) -> List[float]:
        """Extract general sequence properties."""
        length = len(seq)

        # Base counts
        a_count = seq.count("A")
        c_count = seq.count("C")
        g_count = seq.count("G")
        t_count = seq.count("T")

        gc_content = (g_count + c_count) / length if length > 0 else 0.5

        # Terminal GC content
        term_5 = seq[:10] if length >= 10 else seq
        gc_5 = (term_5.count("G") + term_5.count("C")) / len(term_5) if term_5 else 0.5

        term_3 = seq[-10:] if length >= 10 else seq
        gc_3 = (term_3.count("G") + term_3.count("C")) / len(term_3) if term_3 else 0.5

        # Skews
        gc_skew = (
            (g_count - c_count) / (g_count + c_count)
            if (g_count + c_count) > 0
            else 0.0
        )
        at_skew = (
            (a_count - t_count) / (a_count + t_count)
            if (a_count + t_count) > 0
            else 0.0
        )

        return [
            length / 150.0,  # Normalize to typical read length
            gc_content,
            gc_5,
            gc_3,
            (gc_skew + 1) / 2,  # Normalize to [0, 1]
            (at_skew + 1) / 2,
        ]  # 6 features


# Data loading


def parse_aa_damage_file(filepath: str) -> Dict[str, Dict]:
    """
    Parse aMGSIM aa-damage.tsv.gz file.

    Columns include:
    - read_name: the read identifier
    - damage: comma-separated positions with damage (e.g., "19,33" or "None")
    - damage_inframe_ct: C→T damage positions
    - damage_inframe_ag: G→A damage positions
    - Strand_read, Strand_gene: strand information
    - read_length: length of read

    Returns dict mapping read_name to ground truth info.
    """
    ground_truth = {}

    opener = gzip.open if filepath.endswith(".gz") else open

    with opener(filepath, "rt") as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")

            if header is None:
                header = parts
                continue

            if len(parts) < 10:
                continue

            try:
                record = dict(zip(header, parts))
                read_name = record.get("read_name", "")

                # Get damage positions from 'damage' column
                damage_str = record.get("damage", "None")
                read_length = int(record.get("read_length", 100))

                # Parse damage positions
                damage_positions = []
                if damage_str and damage_str != "None":
                    for pos_str in damage_str.split(","):
                        try:
                            pos = int(pos_str.strip())
                            damage_positions.append(pos)
                        except ValueError:
                            continue

                # Check if damage in first/last 10bp
                has_damage_5 = any(p < 10 for p in damage_positions)
                has_damage_3 = any(p >= read_length - 10 for p in damage_positions)

                # Also check C→T and G→A specific damage
                ct_damage = record.get("damage_inframe_ct", "")
                ag_damage = record.get("damage_inframe_ag", "")

                has_ct = bool(ct_damage and ct_damage != "")
                has_ag = bool(ag_damage and ag_damage != "")

                ground_truth[read_name] = {
                    "has_damage": len(damage_positions) > 0,
                    "damage_5prime": has_damage_5,
                    "damage_3prime": has_damage_3,
                    "damage_positions": damage_positions,
                    "has_ct_damage": has_ct,
                    "has_ag_damage": has_ag,
                    "strand_read": record.get("Strand_read", "+"),
                    "strand_gene": record.get("Strand_gene", "+"),
                    "read_length": read_length,
                }

            except Exception as e:
                continue

    return ground_truth


def load_fastq_sequences(filepath: str, max_reads: int = 100000) -> Dict[str, str]:
    """Load sequences from FASTQ file."""
    sequences = {}

    opener = gzip.open if filepath.endswith(".gz") else open

    with opener(filepath, "rt") as f:
        count = 0
        while count < max_reads:
            # Read 4 lines at a time
            header = f.readline().strip()
            if not header:
                break

            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()

            if header.startswith("@"):
                read_name = header[1:].split()[0]
                sequences[read_name] = seq
                count += 1

    return sequences


# Neural network model


class DamageDetectorNN(nn.Module):
    """Neural network for damage detection."""

    def __init__(self, input_dim: int = 62, hidden_dims: List[int] = [64, 32, 16]):
        super().__init__()

        layers = []
        prev_dim = input_dim

        for hidden_dim in hidden_dims:
            layers.extend(
                [
                    nn.Linear(prev_dim, hidden_dim),
                    nn.BatchNorm1d(hidden_dim),
                    nn.ReLU(),
                    nn.Dropout(0.2),
                ]
            )
            prev_dim = hidden_dim

        # Output layer
        layers.append(nn.Linear(prev_dim, 1))

        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return torch.sigmoid(self.network(x))


class DamageDataset(Dataset):
    """Dataset for damage detection training."""

    def __init__(self, features: np.ndarray, labels: np.ndarray):
        self.features = torch.FloatTensor(features)
        self.labels = torch.FloatTensor(labels)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return self.features[idx], self.labels[idx]


# Training


def train_model(
    train_loader: DataLoader,
    val_loader: DataLoader,
    model: nn.Module,
    epochs: int = 100,
    lr: float = 1e-3,
    device: str = "cuda" if torch.cuda.is_available() else "cpu",
) -> Tuple[nn.Module, Dict]:
    """Train the damage detection model."""

    model = model.to(device)
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, epochs)

    best_auc = 0.0
    best_state = None
    history = {"train_loss": [], "val_loss": [], "val_auc": []}

    for epoch in range(epochs):
        # Training
        model.train()
        train_loss = 0.0
        for features, labels in train_loader:
            features, labels = features.to(device), labels.to(device)

            optimizer.zero_grad()
            outputs = model(features).squeeze()
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()

        train_loss /= len(train_loader)

        # Validation
        model.eval()
        val_loss = 0.0
        all_preds = []
        all_labels = []

        with torch.no_grad():
            for features, labels in val_loader:
                features, labels = features.to(device), labels.to(device)
                outputs = model(features).squeeze()
                loss = criterion(outputs, labels)
                val_loss += loss.item()

                all_preds.extend(outputs.cpu().numpy())
                all_labels.extend(labels.cpu().numpy())

        val_loss /= len(val_loader)
        val_auc = roc_auc_score(all_labels, all_preds)

        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)
        history["val_auc"].append(val_auc)

        scheduler.step()

        if val_auc > best_auc:
            best_auc = val_auc
            best_state = model.state_dict().copy()

        if (epoch + 1) % 10 == 0:
            print(
                f"Epoch {epoch + 1}/{epochs}: train_loss={train_loss:.4f}, "
                f"val_loss={val_loss:.4f}, val_auc={val_auc:.4f}"
            )

    # Load best model
    if best_state is not None:
        model.load_state_dict(best_state)

    return model, history


# Export to C++


def export_to_cpp_header(model: nn.Module, filepath: str, input_dim: int = 62):
    """Export trained model weights to C++ header file."""

    with open(filepath, "w") as f:
        f.write("// Auto-generated damage detection neural network weights\n")
        f.write("// Do not edit manually\n\n")
        f.write("#pragma once\n\n")
        f.write("#include <array>\n")
        f.write("#include <cmath>\n\n")
        f.write("namespace agp {\n\n")

        f.write("class DamageNeuralNetwork {\n")
        f.write("public:\n")

        # Write inference function with correct input dimension
        f.write(
            f"    float predict(const std::array<float, {input_dim}>& features) const {{\n"
        )

        layer_idx = 0
        prev_var = "features"
        linear_count = 0

        for name, module in model.network.named_modules():
            if isinstance(module, nn.Linear):
                in_features = module.in_features
                out_features = module.out_features

                f.write(
                    f"        // Layer {linear_count}: {in_features} -> {out_features}\n"
                )
                f.write(f"        std::array<float, {out_features}> h{linear_count};\n")
                f.write(f"        for (int i = 0; i < {out_features}; ++i) {{\n")
                f.write(f"            float sum = b{linear_count}_[i];\n")
                f.write(f"            for (int j = 0; j < {in_features}; ++j) {{\n")
                f.write(
                    f"                sum += {prev_var}[j] * W{linear_count}_[j][i];\n"
                )
                f.write(f"            }}\n")

                # Check if this is the last layer
                if linear_count < 3:  # Hidden layers have activation
                    f.write(
                        f"            h{linear_count}[i] = sum > 0 ? sum : 0;  // ReLU\n"
                    )
                else:  # Output layer
                    f.write(
                        f"            return 1.0f / (1.0f + std::exp(-sum));  // Sigmoid\n"
                    )

                f.write(f"        }}\n")

                prev_var = f"h{linear_count}"
                linear_count += 1

        f.write("        return 0.5f;  // Fallback\n")
        f.write("    }\n\n")

        # Write weights
        f.write("private:\n")
        linear_count = 0

        for name, module in model.network.named_modules():
            if isinstance(module, nn.Linear):
                weights = module.weight.detach().cpu().numpy()
                biases = module.bias.detach().cpu().numpy()

                in_features, out_features = weights.shape[1], weights.shape[0]

                # Weight matrix (transposed for row-major access)
                f.write(
                    f"    static constexpr float W{linear_count}_[{in_features}][{out_features}] = {{\n"
                )
                for i in range(in_features):
                    f.write("        {")
                    f.write(
                        ", ".join(f"{weights[j, i]:.8f}f" for j in range(out_features))
                    )
                    f.write("},\n")
                f.write("    };\n\n")

                # Bias vector
                f.write(
                    f"    static constexpr float b{linear_count}_[{out_features}] = {{\n        "
                )
                f.write(", ".join(f"{b:.8f}f" for b in biases))
                f.write("\n    };\n\n")

                linear_count += 1

        f.write("};\n\n")
        f.write("} // namespace agp\n")

    print(f"Exported model to {filepath}")


# Main


def main():
    parser = argparse.ArgumentParser(description="Train damage detection NN")
    parser.add_argument(
        "--data-dir", type=str, required=True, help="Directory with benchmark data"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="damage_nn_weights.hpp",
        help="Output C++ header file",
    )
    parser.add_argument("--epochs", type=int, default=100, help="Training epochs")
    parser.add_argument("--batch-size", type=int, default=256, help="Batch size")
    parser.add_argument(
        "--max-reads", type=int, default=100000, help="Max reads per dataset"
    )
    args = parser.parse_args()

    data_dir = Path(args.data_dir)

    # Find datasets
    fastq_files = list(data_dir.glob("*_art.fq.gz"))
    damage_files = list(data_dir.glob("*_aa-damage.tsv.gz"))

    print(f"Found {len(fastq_files)} FASTQ files and {len(damage_files)} damage files")

    # Build training data
    all_features = []
    all_labels = []

    extractor = DamageFeatureExtractor()

    for fastq_path in fastq_files[:5]:  # Limit for initial training
        # Find matching damage file
        base_name = fastq_path.name.replace("_art.fq.gz", "")
        damage_path = data_dir / f"{base_name}_aa-damage.tsv.gz"

        if not damage_path.exists():
            print(f"Skipping {base_name} - no damage file")
            continue

        print(f"Processing {base_name}...")

        # Load data
        sequences = load_fastq_sequences(str(fastq_path), args.max_reads)
        ground_truth = parse_aa_damage_file(str(damage_path))

        print(
            f"  Loaded {len(sequences)} sequences, {len(ground_truth)} ground truth entries"
        )

        # Extract features for reads with ground truth
        for read_name, seq in sequences.items():
            if read_name in ground_truth:
                gt = ground_truth[read_name]

                # Use terminal damage as label (first/last 10bp)
                # This is what we actually want to detect
                label = 1.0 if (gt["damage_5prime"] or gt["damage_3prime"]) else 0.0

                features = extractor.extract_features(seq)
                all_features.append(features)
                all_labels.append(label)

    if len(all_features) == 0:
        print("No training data found!")
        sys.exit(1)

    X = np.array(all_features)
    y = np.array(all_labels)

    print(f"\nTotal samples: {len(y)}")
    print(f"Damaged: {sum(y)} ({100 * sum(y) / len(y):.1f}%)")
    print(f"Undamaged: {len(y) - sum(y)} ({100 * (1 - sum(y) / len(y)):.1f}%)")

    # Split data
    X_train, X_temp, y_train, y_temp = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    X_val, X_test, y_val, y_test = train_test_split(
        X_temp, y_temp, test_size=0.5, random_state=42, stratify=y_temp
    )

    print(f"\nTrain: {len(y_train)}, Val: {len(y_val)}, Test: {len(y_test)}")

    # Create data loaders
    train_dataset = DamageDataset(X_train, y_train)
    val_dataset = DamageDataset(X_val, y_val)
    test_dataset = DamageDataset(X_test, y_test)

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=args.batch_size)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size)

    # Create and train model
    # Calculate actual feature count
    input_dim = len(all_features[0])
    print(f"Feature dimension: {input_dim}")

    model = DamageDetectorNN(input_dim=input_dim, hidden_dims=[64, 32, 16])
    print(f"\nModel parameters: {sum(p.numel() for p in model.parameters())}")

    model, history = train_model(train_loader, val_loader, model, epochs=args.epochs)

    # Test evaluation
    model.eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)

    all_preds = []
    all_labels = []

    with torch.no_grad():
        for features, labels in test_loader:
            features = features.to(device)
            outputs = model(features).squeeze()
            all_preds.extend(outputs.cpu().numpy())
            all_labels.extend(labels.numpy())

    test_auc = roc_auc_score(all_labels, all_preds)
    print(f"\nTest AUC-ROC: {test_auc:.4f}")

    # Find best threshold
    precision, recall, thresholds = precision_recall_curve(all_labels, all_preds)
    f1_scores = 2 * precision * recall / (precision + recall + 1e-8)
    best_idx = np.argmax(f1_scores)
    best_threshold = thresholds[best_idx] if best_idx < len(thresholds) else 0.5

    print(f"Best threshold: {best_threshold:.3f}")
    print(f"Best F1: {f1_scores[best_idx]:.4f}")
    print(f"Precision: {precision[best_idx]:.4f}")
    print(f"Recall: {recall[best_idx]:.4f}")

    # Export to C++
    export_to_cpp_header(model, args.output)

    print("\nDone!")


if __name__ == "__main__":
    main()
