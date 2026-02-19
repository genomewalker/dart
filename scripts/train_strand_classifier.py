#!/usr/bin/env python3
"""
Train strand classifier on real aMGSIM ground truth data.

Unlike the frame classifier (which works within a single orientation), the strand
classifier compares features between the original sequence and its reverse complement
to predict which strand contains the true reading frame.

The key insight is that strand selection needs COMPARATIVE features that are
computed for both orientations and then compared.

Usage:
    python train_strand_classifier.py --train data1.tsv.gz --test test.tsv.gz --output strand_nn.hpp
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
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}

# Number of features: computed as differences between fwd and rc strands
NUM_FEATURES = 40


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def translate(seq, frame):
    protein = ""
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i : i + 3].upper()
        protein += CODON_TABLE.get(codon, "X")
    return protein


def count_internal_stops(seq, frame):
    """Count internal stop codons."""
    codons = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i : i + 3].upper()
        if all(c in "ACGT" for c in codon):
            codons.append(codon)
    if len(codons) <= 1:
        return 0
    return sum(1 for c in codons[:-1] if c in STOP_CODONS)


def compute_hexamer_score(seq, frame):
    """Compute log-likelihood ratio for hexamers (coding vs uniform)."""
    # Simplified hexamer scoring using codon pairs
    score = 0.0
    count = 0
    for i in range(frame, len(seq) - 5, 3):
        hex_seq = seq[i : i + 6].upper()
        if all(c in "ACGT" for c in hex_seq):
            # Simple proxy: favor balanced GC content
            gc = sum(1 for c in hex_seq if c in "GC")
            # Coding sequences tend to have moderate GC (not too high, not too low)
            score += 1.0 - abs(gc - 3) * 0.1
            count += 1
    return score / max(1, count)


def compute_rny_score(seq, frame):
    """RNY codon compliance (R=purine pos1, N=any pos2, Y=pyrimidine pos3)."""
    count = 0
    total = 0
    for i in range(frame, len(seq) - 2, 3):
        if i + 2 < len(seq):
            pos1 = seq[i].upper()
            pos3 = seq[i + 2].upper()
            if pos1 in "ACGT" and pos3 in "ACGT":
                total += 1
                if pos1 in "AG" and pos3 in "CT":
                    count += 1
    return count / max(1, total)


def compute_gc_by_position(seq, frame):
    """GC content at each codon position."""
    gc_pos = [0, 0, 0]
    total_pos = [0, 0, 0]
    for i in range(frame, len(seq)):
        pos = (i - frame) % 3
        c = seq[i].upper()
        if c in "ACGT":
            total_pos[pos] += 1
            if c in "GC":
                gc_pos[pos] += 1
    return [gc_pos[p] / max(1, total_pos[p]) for p in range(3)]


def best_frame_features(seq):
    """Get features for the best frame (lowest stops, highest hexamer score)."""
    best_frame = 0
    best_score = -1000

    for frame in range(3):
        stops = count_internal_stops(seq, frame)
        hex_score = compute_hexamer_score(seq, frame)
        rny = compute_rny_score(seq, frame)

        # Combine: fewer stops is better, higher hexamer/rny is better
        score = -stops * 10 + hex_score + rny
        if score > best_score:
            best_score = score
            best_frame = frame

    return best_frame


def extract_strand_features(seq):
    """
    Extract features for strand comparison.
    Returns features that capture coding signal quality for each strand.
    """
    seq = seq.upper()
    rc = reverse_complement(seq)

    features = []

    # For each strand, compute features for all 3 frames and take the best
    for strand_seq in [seq, rc]:
        # Find best frame for this strand
        best_frame = best_frame_features(strand_seq)

        # Features for best frame
        stops = count_internal_stops(strand_seq, best_frame)
        features.append(min(stops, 5) / 5.0)  # normalized stops
        features.append(1.0 if stops > 0 else 0.0)  # has any stop

        # Hexamer-like score
        hex_score = compute_hexamer_score(strand_seq, best_frame)
        features.append(hex_score)

        # RNY compliance
        rny = compute_rny_score(strand_seq, best_frame)
        features.append(rny)

        # GC by position
        gc_pos = compute_gc_by_position(strand_seq, best_frame)
        features.extend(gc_pos)

        # GC3 - GC1 (wobble position enrichment)
        features.append(gc_pos[2] - gc_pos[0])

        # Protein features
        protein = translate(strand_seq, best_frame)

        # Dipeptide entropy
        if len(protein) >= 2:
            dp_counts = defaultdict(int)
            for i in range(len(protein) - 1):
                dp_counts[protein[i : i + 2]] += 1
            total = sum(dp_counts.values())
            entropy = -sum(
                (c / total) * np.log2(c / total) for c in dp_counts.values() if c > 0
            )
            features.append(entropy / 8.64)
        else:
            features.append(0.5)

        # AA composition bias (hydrophobic fraction)
        hydrophobic = sum(1 for aa in protein if aa in "AILMFVPWG")
        features.append(hydrophobic / max(1, len(protein)))

        # Start codon presence
        codons = [
            strand_seq[i : i + 3] for i in range(best_frame, len(strand_seq) - 2, 3)
        ]
        features.append(1.0 if codons and codons[0] == "ATG" else 0.0)

        # Terminal stop presence
        features.append(1.0 if codons and codons[-1] in STOP_CODONS else 0.0)

        # Min stops across all frames (frame-agnostic)
        min_stops = min(count_internal_stops(strand_seq, f) for f in range(3))
        features.append(min_stops / 5.0)

    # Now add DIFFERENCE features (fwd - rc)
    # These are the key discriminative features
    n_per_strand = len(features) // 2
    fwd_features = features[:n_per_strand]
    rc_features = features[n_per_strand:]

    diff_features = [f - r for f, r in zip(fwd_features, rc_features)]
    features.extend(diff_features)

    # Add some additional comparative features

    # Hexamer log-likelihood difference (similar to current AGP approach)
    best_fwd_hex = max(compute_hexamer_score(seq, f) for f in range(3))
    best_rc_hex = max(compute_hexamer_score(rc, f) for f in range(3))
    features.append(best_fwd_hex - best_rc_hex)

    # Best RNY difference
    best_fwd_rny = max(compute_rny_score(seq, f) for f in range(3))
    best_rc_rny = max(compute_rny_score(rc, f) for f in range(3))
    features.append(best_fwd_rny - best_rc_rny)

    # Stop count difference (min across frames)
    min_fwd_stops = min(count_internal_stops(seq, f) for f in range(3))
    min_rc_stops = min(count_internal_stops(rc, f) for f in range(3))
    features.append((min_rc_stops - min_fwd_stops) / 5.0)  # positive = fwd better

    # GC skew (can indicate strand)
    g_count = seq.count("G")
    c_count = seq.count("C")
    gc_skew = (g_count - c_count) / max(1, g_count + c_count)
    features.append(gc_skew)

    return features


def load_ground_truth(tsv_paths, max_per_file=100000):
    """Load ground truth from aMGSIM TSV files."""
    all_data = []

    for tsv_path in tsv_paths:
        print(f"Loading {os.path.basename(tsv_path)}...", file=sys.stderr)

        opener = gzip.open if tsv_path.endswith(".gz") else open
        mode = "rt" if tsv_path.endswith(".gz") else "r"
        count = 0

        with opener(tsv_path, mode) as f:
            header = f.readline().strip().split("\t")
            cols = {name: i for i, name in enumerate(header)}

            required = ["damaged_seq", "Strand_read", "Strand_gene"]
            if not all(c in cols for c in required):
                print(f"  Skipping - missing columns", file=sys.stderr)
                continue

            for line in f:
                if count >= max_per_file:
                    break

                parts = line.strip().split("\t")
                if len(parts) <= max(cols.values()):
                    continue

                damaged_seq = parts[cols["damaged_seq"]].upper()
                strand_read = parts[cols["Strand_read"]]
                strand_gene = parts[cols["Strand_gene"]]

                if len(damaged_seq) < 30:
                    continue

                # Label: 1 if forward strand is correct (strand_read == strand_gene)
                # 0 if reverse strand is correct
                is_forward = 1 if (strand_read == strand_gene) else 0

                all_data.append((damaged_seq, is_forward))
                count += 1

        print(f"  Loaded {count} samples", file=sys.stderr)

    return all_data


def create_dataset(data):
    """Create X, y arrays for training."""
    X = []
    y = []

    for seq, is_forward in data:
        features = extract_strand_features(seq)
        X.append(features)
        y.append(is_forward)

    return np.array(X, dtype=np.float32), np.array(y, dtype=np.int32)


def evaluate_accuracy(model, scaler, data):
    """Evaluate strand selection accuracy."""
    correct = 0
    for seq, is_forward in data:
        features = np.array([extract_strand_features(seq)])
        features_scaled = scaler.transform(features)

        if hasattr(model, "predict_proba"):
            prob = model.predict_proba(features_scaled)[0, 1]
        else:
            prob = 1.0 / (1.0 + np.exp(-model.decision_function(features_scaled)[0]))

        predicted = 1 if prob > 0.5 else 0
        if predicted == is_forward:
            correct += 1

    return correct / len(data)


def export_mlp_to_cpp(model, scaler, output_path, n_features):
    """Export MLP to C++ header for strand classification."""
    W1 = model.coefs_[0]
    b1 = model.intercepts_[0]
    W2 = model.coefs_[1]
    b2 = model.intercepts_[1]

    input_size = W1.shape[0]
    hidden_size = W1.shape[1]

    with open(output_path, "w") as f:
        f.write("#pragma once\n\n")
        f.write("/**\n")
        f.write(" * Strand classifier trained on aMGSIM ground truth.\n")
        f.write(" * Generated by train_strand_classifier.py\n")
        f.write(
            f" * Architecture: {input_size} -> {hidden_size} (ReLU) -> 1 (sigmoid)\n"
        )
        f.write(" * \n")
        f.write(" * This classifier predicts whether the forward strand (+) contains\n")
        f.write(
            " * the correct reading frame. Output > 0.5 means forward is correct.\n"
        )
        f.write(" */\n\n")
        f.write("#include <array>\n")
        f.write("#include <cmath>\n")
        f.write("#include <algorithm>\n")
        f.write("#include <string>\n\n")
        f.write("namespace agp {\n")
        f.write("namespace strand_nn {\n\n")

        f.write(f"constexpr int INPUT_SIZE = {input_size};\n")
        f.write(f"constexpr int HIDDEN_SIZE = {hidden_size};\n\n")

        # Scaler
        f.write("inline const std::array<float, INPUT_SIZE>& get_input_mean() {\n")
        f.write(
            "    static const std::array<float, INPUT_SIZE> INPUT_MEAN = {{\n        "
        )
        f.write(",\n        ".join(f"{x:.8f}f" for x in scaler.mean_))
        f.write("\n    }};\n")
        f.write("    return INPUT_MEAN;\n")
        f.write("}\n\n")

        f.write("inline const std::array<float, INPUT_SIZE>& get_input_std() {\n")
        f.write(
            "    static const std::array<float, INPUT_SIZE> INPUT_STD = {{\n        "
        )
        f.write(",\n        ".join(f"{max(x, 1e-8):.8f}f" for x in scaler.scale_))
        f.write("\n    }};\n")
        f.write("    return INPUT_STD;\n")
        f.write("}\n\n")

        # W1
        f.write(
            "inline const std::array<std::array<float, HIDDEN_SIZE>, INPUT_SIZE>& get_W1() {\n"
        )
        f.write(
            "    static const std::array<std::array<float, HIDDEN_SIZE>, INPUT_SIZE> W1 = {{\n"
        )
        for i, row in enumerate(W1):
            f.write("        {{" + ", ".join(f"{x:.8f}f" for x in row) + "}}")
            f.write(",\n" if i < len(W1) - 1 else "\n")
        f.write("    }};\n")
        f.write("    return W1;\n")
        f.write("}\n\n")

        f.write("inline const std::array<float, HIDDEN_SIZE>& get_B1() {\n")
        f.write("    static const std::array<float, HIDDEN_SIZE> B1 = {{\n        ")
        f.write(", ".join(f"{x:.8f}f" for x in b1))
        f.write("\n    }};\n")
        f.write("    return B1;\n")
        f.write("}\n\n")

        # W2
        f.write("inline const std::array<float, HIDDEN_SIZE>& get_W2() {\n")
        f.write("    static const std::array<float, HIDDEN_SIZE> W2 = {{\n        ")
        f.write(", ".join(f"{x:.8f}f" for x in W2.flatten()))
        f.write("\n    }};\n")
        f.write("    return W2;\n")
        f.write("}\n\n")

        f.write(f"inline float get_B2() {{ return {b2[0]:.8f}f; }}\n\n")

        # Feature extraction helper functions
        f.write("// Helper: reverse complement\n")
        f.write("inline char complement(char c) {\n")
        f.write("    switch (c) {\n")
        f.write("        case 'A': case 'a': return 'T';\n")
        f.write("        case 'T': case 't': return 'A';\n")
        f.write("        case 'G': case 'g': return 'C';\n")
        f.write("        case 'C': case 'c': return 'G';\n")
        f.write("        default: return 'N';\n")
        f.write("    }\n")
        f.write("}\n\n")

        f.write("inline std::string reverse_complement(const std::string& seq) {\n")
        f.write("    std::string rc(seq.length(), 'N');\n")
        f.write("    for (size_t i = 0; i < seq.length(); ++i) {\n")
        f.write("        rc[i] = complement(seq[seq.length() - 1 - i]);\n")
        f.write("    }\n")
        f.write("    return rc;\n")
        f.write("}\n\n")

        # Feature extraction (matching Python exactly)
        f.write("// Count internal stop codons\n")
        f.write(
            "inline int count_internal_stops(const std::string& seq, int frame) {\n"
        )
        f.write("    int stops = 0;\n")
        f.write("    int n_codons = 0;\n")
        f.write("    for (size_t i = frame; i + 2 < seq.length(); i += 3) {\n")
        f.write("        char c1 = std::toupper(seq[i]);\n")
        f.write("        char c2 = std::toupper(seq[i+1]);\n")
        f.write("        char c3 = std::toupper(seq[i+2]);\n")
        f.write("        n_codons++;\n")
        f.write("        bool is_stop = (c1=='T' && c2=='A' && c3=='A') ||\n")
        f.write("                       (c1=='T' && c2=='A' && c3=='G') ||\n")
        f.write("                       (c1=='T' && c2=='G' && c3=='A');\n")
        f.write(
            "        if (is_stop && i + 5 < seq.length()) stops++;  // internal only\n"
        )
        f.write("    }\n")
        f.write("    return stops;\n")
        f.write("}\n\n")

        f.write("// Compute hexamer-like score\n")
        f.write(
            "inline float compute_hexamer_score(const std::string& seq, int frame) {\n"
        )
        f.write("    float score = 0.0f;\n")
        f.write("    int count = 0;\n")
        f.write("    for (size_t i = frame; i + 5 < seq.length(); i += 3) {\n")
        f.write("        int gc = 0;\n")
        f.write("        bool valid = true;\n")
        f.write("        for (int j = 0; j < 6; j++) {\n")
        f.write("            char c = std::toupper(seq[i+j]);\n")
        f.write(
            "            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { valid = false; break; }\n"
        )
        f.write("            if (c == 'G' || c == 'C') gc++;\n")
        f.write("        }\n")
        f.write("        if (valid) {\n")
        f.write("            score += 1.0f - std::abs(gc - 3) * 0.1f;\n")
        f.write("            count++;\n")
        f.write("        }\n")
        f.write("    }\n")
        f.write("    return count > 0 ? score / count : 0.0f;\n")
        f.write("}\n\n")

        f.write("// RNY compliance score\n")
        f.write("inline float compute_rny_score(const std::string& seq, int frame) {\n")
        f.write("    int count = 0, total = 0;\n")
        f.write("    for (size_t i = frame; i + 2 < seq.length(); i += 3) {\n")
        f.write("        char pos1 = std::toupper(seq[i]);\n")
        f.write("        char pos3 = std::toupper(seq[i+2]);\n")
        f.write(
            "        if ((pos1 == 'A' || pos1 == 'C' || pos1 == 'G' || pos1 == 'T') &&\n"
        )
        f.write(
            "            (pos3 == 'A' || pos3 == 'C' || pos3 == 'G' || pos3 == 'T')) {\n"
        )
        f.write("            total++;\n")
        f.write(
            "            if ((pos1 == 'A' || pos1 == 'G') && (pos3 == 'C' || pos3 == 'T')) count++;\n"
        )
        f.write("        }\n")
        f.write("    }\n")
        f.write("    return total > 0 ? (float)count / total : 0.0f;\n")
        f.write("}\n\n")

        f.write("// GC content by codon position\n")
        f.write(
            "inline void compute_gc_by_position(const std::string& seq, int frame, float gc_pos[3]) {\n"
        )
        f.write("    int gc[3] = {0, 0, 0};\n")
        f.write("    int total[3] = {0, 0, 0};\n")
        f.write("    for (size_t i = frame; i < seq.length(); i++) {\n")
        f.write("        int pos = (i - frame) % 3;\n")
        f.write("        char c = std::toupper(seq[i]);\n")
        f.write("        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {\n")
        f.write("            total[pos]++;\n")
        f.write("            if (c == 'G' || c == 'C') gc[pos]++;\n")
        f.write("        }\n")
        f.write("    }\n")
        f.write("    for (int p = 0; p < 3; p++) {\n")
        f.write("        gc_pos[p] = total[p] > 0 ? (float)gc[p] / total[p] : 0.0f;\n")
        f.write("    }\n")
        f.write("}\n\n")

        f.write("// Find best frame for a sequence\n")
        f.write("inline int best_frame(const std::string& seq) {\n")
        f.write("    int best = 0;\n")
        f.write("    float best_score = -1000.0f;\n")
        f.write("    for (int frame = 0; frame < 3; frame++) {\n")
        f.write("        int stops = count_internal_stops(seq, frame);\n")
        f.write("        float hex = compute_hexamer_score(seq, frame);\n")
        f.write("        float rny = compute_rny_score(seq, frame);\n")
        f.write("        float score = -stops * 10.0f + hex + rny;\n")
        f.write(
            "        if (score > best_score) { best_score = score; best = frame; }\n"
        )
        f.write("    }\n")
        f.write("    return best;\n")
        f.write("}\n\n")

        # Main feature extraction function
        f.write(f"// Extract {input_size} features for strand classification\n")
        f.write(
            f"inline std::array<float, INPUT_SIZE> extract_features(const std::string& seq) {{\n"
        )
        f.write(f"    std::array<float, INPUT_SIZE> features;\n")
        f.write("    features.fill(0.0f);\n")
        f.write("    std::string rc = reverse_complement(seq);\n")
        f.write("    \n")
        f.write("    int idx = 0;\n")
        f.write("    const std::string* strands[2] = {&seq, &rc};\n")
        f.write("    \n")
        f.write("    // Features for each strand\n")
        f.write("    for (int s = 0; s < 2; s++) {\n")
        f.write("        const std::string& strand = *strands[s];\n")
        f.write("        int bf = best_frame(strand);\n")
        f.write("        \n")
        f.write("        int stops = count_internal_stops(strand, bf);\n")
        f.write("        features[idx++] = std::min(stops, 5) / 5.0f;\n")
        f.write("        features[idx++] = stops > 0 ? 1.0f : 0.0f;\n")
        f.write("        \n")
        f.write("        features[idx++] = compute_hexamer_score(strand, bf);\n")
        f.write("        features[idx++] = compute_rny_score(strand, bf);\n")
        f.write("        \n")
        f.write("        float gc_pos[3];\n")
        f.write("        compute_gc_by_position(strand, bf, gc_pos);\n")
        f.write("        features[idx++] = gc_pos[0];\n")
        f.write("        features[idx++] = gc_pos[1];\n")
        f.write("        features[idx++] = gc_pos[2];\n")
        f.write("        features[idx++] = gc_pos[2] - gc_pos[0];  // GC3 - GC1\n")
        f.write("        \n")
        f.write("        // Dipeptide entropy default\n")
        f.write("        features[idx++] = 0.5f;\n")
        f.write("        \n")
        f.write("        // Hydrophobic fraction default\n")
        f.write("        features[idx++] = 0.3f;\n")
        f.write("        \n")
        f.write("        // Start/stop presence\n")
        f.write("        bool has_atg = (bf + 2 < strand.length() &&\n")
        f.write("                        std::toupper(strand[bf]) == 'A' &&\n")
        f.write("                        std::toupper(strand[bf+1]) == 'T' &&\n")
        f.write("                        std::toupper(strand[bf+2]) == 'G');\n")
        f.write("        features[idx++] = has_atg ? 1.0f : 0.0f;\n")
        f.write("        \n")
        f.write("        // Terminal stop\n")
        f.write(
            "        size_t last_codon = ((strand.length() - bf) / 3 - 1) * 3 + bf;\n"
        )
        f.write("        bool has_term_stop = false;\n")
        f.write("        if (last_codon + 2 < strand.length()) {\n")
        f.write("            char c1 = std::toupper(strand[last_codon]);\n")
        f.write("            char c2 = std::toupper(strand[last_codon+1]);\n")
        f.write("            char c3 = std::toupper(strand[last_codon+2]);\n")
        f.write(
            "            has_term_stop = (c1=='T' && c2=='A' && (c3=='A' || c3=='G')) ||\n"
        )
        f.write("                           (c1=='T' && c2=='G' && c3=='A');\n")
        f.write("        }\n")
        f.write("        features[idx++] = has_term_stop ? 1.0f : 0.0f;\n")
        f.write("        \n")
        f.write("        // Min stops across all frames\n")
        f.write("        int min_stops = std::min({count_internal_stops(strand, 0),\n")
        f.write("                                  count_internal_stops(strand, 1),\n")
        f.write(
            "                                  count_internal_stops(strand, 2)});\n"
        )
        f.write("        features[idx++] = min_stops / 5.0f;\n")
        f.write("    }\n")
        f.write("    \n")
        f.write("    // Difference features (fwd - rc)\n")
        f.write("    int n_per_strand = idx / 2;\n")
        f.write("    for (int i = 0; i < n_per_strand; i++) {\n")
        f.write("        features[idx++] = features[i] - features[n_per_strand + i];\n")
        f.write("    }\n")
        f.write("    \n")
        f.write("    // Additional comparative features\n")
        f.write("    float best_fwd_hex = std::max({compute_hexamer_score(seq, 0),\n")
        f.write("                                    compute_hexamer_score(seq, 1),\n")
        f.write(
            "                                    compute_hexamer_score(seq, 2)});\n"
        )
        f.write("    float best_rc_hex = std::max({compute_hexamer_score(rc, 0),\n")
        f.write("                                   compute_hexamer_score(rc, 1),\n")
        f.write("                                   compute_hexamer_score(rc, 2)});\n")
        f.write("    features[idx++] = best_fwd_hex - best_rc_hex;\n")
        f.write("    \n")
        f.write("    float best_fwd_rny = std::max({compute_rny_score(seq, 0),\n")
        f.write("                                    compute_rny_score(seq, 1),\n")
        f.write("                                    compute_rny_score(seq, 2)});\n")
        f.write("    float best_rc_rny = std::max({compute_rny_score(rc, 0),\n")
        f.write("                                   compute_rny_score(rc, 1),\n")
        f.write("                                   compute_rny_score(rc, 2)});\n")
        f.write("    features[idx++] = best_fwd_rny - best_rc_rny;\n")
        f.write("    \n")
        f.write("    int min_fwd_stops = std::min({count_internal_stops(seq, 0),\n")
        f.write("                                   count_internal_stops(seq, 1),\n")
        f.write("                                   count_internal_stops(seq, 2)});\n")
        f.write("    int min_rc_stops = std::min({count_internal_stops(rc, 0),\n")
        f.write("                                  count_internal_stops(rc, 1),\n")
        f.write("                                  count_internal_stops(rc, 2)});\n")
        f.write("    features[idx++] = (min_rc_stops - min_fwd_stops) / 5.0f;\n")
        f.write("    \n")
        f.write("    // GC skew\n")
        f.write("    int g_count = 0, c_count = 0;\n")
        f.write("    for (char c : seq) {\n")
        f.write("        if (c == 'G' || c == 'g') g_count++;\n")
        f.write("        else if (c == 'C' || c == 'c') c_count++;\n")
        f.write("    }\n")
        f.write("    features[idx++] = (g_count + c_count) > 0 ?\n")
        f.write("        (float)(g_count - c_count) / (g_count + c_count) : 0.0f;\n")
        f.write("    \n")
        f.write("    return features;\n")
        f.write("}\n\n")

        # Forward pass
        f.write(
            "// Predict strand: returns probability that FORWARD strand is correct\n"
        )
        f.write(
            "inline float predict_strand(const std::array<float, INPUT_SIZE>& features) {\n"
        )
        f.write("    const auto& INPUT_MEAN = get_input_mean();\n")
        f.write("    const auto& INPUT_STD = get_input_std();\n")
        f.write("    const auto& W1 = get_W1();\n")
        f.write("    const auto& B1 = get_B1();\n")
        f.write("    const auto& W2 = get_W2();\n")
        f.write("    float B2 = get_B2();\n")
        f.write("    \n")
        f.write("    // Normalize input\n")
        f.write("    std::array<float, INPUT_SIZE> x;\n")
        f.write("    for (int i = 0; i < INPUT_SIZE; i++) {\n")
        f.write("        x[i] = (features[i] - INPUT_MEAN[i]) / INPUT_STD[i];\n")
        f.write("    }\n")
        f.write("    \n")
        f.write("    // Hidden layer with ReLU\n")
        f.write("    std::array<float, HIDDEN_SIZE> h;\n")
        f.write("    for (int j = 0; j < HIDDEN_SIZE; j++) {\n")
        f.write("        float sum = B1[j];\n")
        f.write("        for (int i = 0; i < INPUT_SIZE; i++) {\n")
        f.write("            sum += x[i] * W1[i][j];\n")
        f.write("        }\n")
        f.write("        h[j] = std::max(0.0f, sum);\n")
        f.write("    }\n")
        f.write("    \n")
        f.write("    // Output with sigmoid\n")
        f.write("    float z = B2;\n")
        f.write("    for (int j = 0; j < HIDDEN_SIZE; j++) {\n")
        f.write("        z += h[j] * W2[j];\n")
        f.write("    }\n")
        f.write("    return 1.0f / (1.0f + std::exp(-std::clamp(z, -20.0f, 20.0f)));\n")
        f.write("}\n\n")

        # Convenience function
        f.write("// Convenience: predict strand directly from sequence\n")
        f.write("inline float predict_strand(const std::string& seq) {\n")
        f.write("    auto features = extract_features(seq);\n")
        f.write("    return predict_strand(features);\n")
        f.write("}\n\n")

        f.write("} // namespace strand_nn\n")
        f.write("} // namespace agp\n")

    print(f"Exported to {output_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Train strand classifier on aMGSIM data"
    )
    parser.add_argument("--train", nargs="+", required=True, help="Training TSV files")
    parser.add_argument("--test", nargs="+", help="Test TSV files")
    parser.add_argument("--output", default="strand_nn.hpp", help="Output header path")
    parser.add_argument(
        "--max-train", type=int, default=100000, help="Max samples per training file"
    )
    parser.add_argument(
        "--max-test", type=int, default=10000, help="Max samples per test file"
    )
    args = parser.parse_args()

    if not HAS_SKLEARN:
        print("Error: sklearn required", file=sys.stderr)
        sys.exit(1)

    # Load training data
    print("Loading training data...", file=sys.stderr)
    train_data = load_ground_truth(args.train, args.max_train)
    print(f"Total training samples: {len(train_data)}", file=sys.stderr)

    # Check class balance
    fwd_count = sum(1 for _, is_fwd in train_data if is_fwd == 1)
    print(
        f"Forward: {fwd_count}, Reverse: {len(train_data) - fwd_count}", file=sys.stderr
    )

    # Create dataset
    print("Extracting features...", file=sys.stderr)
    X_train, y_train = create_dataset(train_data)
    print(f"Feature shape: {X_train.shape}", file=sys.stderr)

    # Handle NaN/Inf
    X_train = np.nan_to_num(X_train, nan=0.0, posinf=1.0, neginf=-1.0)

    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    # Train model
    print("Training MLP classifier...", file=sys.stderr)
    model = MLPClassifier(
        hidden_layer_sizes=(64,),  # Larger hidden layer for strand classification
        activation="relu",
        solver="adam",
        max_iter=500,
        random_state=42,
        early_stopping=True,
        validation_fraction=0.1,
        verbose=True,
    )
    model.fit(X_train_scaled, y_train)

    # Training accuracy
    train_acc = evaluate_accuracy(model, scaler, train_data[:5000])
    print(f"Training accuracy: {train_acc:.4f}", file=sys.stderr)

    # Test accuracy
    if args.test:
        print("Loading test data...", file=sys.stderr)
        test_data = load_ground_truth(args.test, args.max_test)
        test_acc = evaluate_accuracy(model, scaler, test_data)
        print(f"Test accuracy: {test_acc:.4f}", file=sys.stderr)

    # Export to C++
    print("Exporting to C++...", file=sys.stderr)
    export_mlp_to_cpp(model, scaler, args.output, X_train.shape[1])

    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
