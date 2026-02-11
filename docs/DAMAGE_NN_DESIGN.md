# Neural Network Damage Detection Design

## Overview

Train a lightweight neural network to detect ancient DNA damage by combining:
1. **Terminal nucleotide patterns** (from Pass 1 sample profile)
2. **Codon-based mutation patterns** (frame-aware)
3. **Hexamer distributions** (GTDB-derived)

## Objectives

1. **Detect damaged reads**: Binary classification (damaged vs undamaged)
2. **Improve gene prediction**: Use damage probability to adjust frame selection

## Feature Engineering

### Feature Group 1: Terminal Damage Signals (12 features)

Position-weighted nucleotide frequencies at 5' and 3' ends:

```
5' end (positions 0-4):
  - T/(T+C) ratio at each position (5 values)
  - Weighted sum with exp(-λp) decay

3' end (positions -1 to -5):
  - A/(A+G) ratio at each position (5 values)
  - Weighted sum with exp(-λp) decay

Derived:
  - excess_T_5prime = T_ratio_5' - baseline_TC
  - excess_A_3prime = A_ratio_3' - baseline_AG
```

### Feature Group 2: Sample Profile Context (8 features)

From Pass 1 aggregation:
```
- sample_damage_rate_5prime (float)
- sample_damage_rate_3prime (float)
- sample_decay_lambda_5 (float)
- sample_decay_lambda_3 (float)
- sample_baseline_TC (float)
- sample_baseline_AG (float)
- is_double_stranded (bool → 0/1)
- sample_damage_level (0=minimal, 1=moderate, 2=high)
```

### Feature Group 3: Codon-Based Features (18 features)

Frame-specific mutation analysis:
```
For each frame (0, 1, 2):
  - internal_stop_count
  - stop_codon_types: TAA_count, TAG_count, TGA_count
  - damaged_stop_probability (stops in high-damage positions)
  - synonymous_T_enrichment (T at wobble positions)
```

### Feature Group 4: Hexamer Features (12 features)

Using GTDB hexamer tables:
```
5' region (first 12bp, 7 hexamers):
  - mean_hexamer_freq (from GTDB_HEXAMER_FREQ)
  - mean_damage_prob (from get_hexamer_damage_prob)
  - count_high_damage_hexamers (prob > 0.9)

3' region (last 12bp, 7 hexamers):
  - Same 3 features

Middle region:
  - mean_hexamer_freq (baseline coding signal)
  - hexamer_coding_score (log-likelihood ratio vs random)

Ratios:
  - terminal_vs_middle_freq_ratio
  - damage_hexamer_enrichment_5prime
  - damage_hexamer_enrichment_3prime
```

### Feature Group 5: Sequence Properties (6 features)

```
- sequence_length
- gc_content
- gc_content_5prime (first 10bp)
- gc_content_3prime (last 10bp)
- gc_skew = (G-C)/(G+C)
- at_skew = (A-T)/(A+T)
```

**Total: ~56 features**

## Neural Network Architecture

### Option A: Simple MLP (Recommended for C++ deployment)

```
Input (56 features)
    ↓
Dense(64, ReLU) + BatchNorm + Dropout(0.2)
    ↓
Dense(32, ReLU) + BatchNorm + Dropout(0.2)
    ↓
Dense(16, ReLU)
    ↓
Dense(1, Sigmoid) → P(damaged)
```

**Rationale**:
- Small enough to embed weights in C++ header (~6K parameters)
- Fast inference (~1µs per read)
- No external dependencies

### Option B: Feature Interaction Network

```
Terminal Features (12) → Dense(16) → Concat
Sample Features (8)    → Dense(8)  → Concat → Dense(32) → Dense(1)
Codon Features (18)    → Dense(16) → Concat
Hexamer Features (12)  → Dense(16) → Concat
Sequence Features (6)  → Pass-through → Concat
```

**Rationale**:
- Learns feature group interactions
- Better interpretability

## Training Strategy

### Dataset Construction

**Positive examples (damaged reads)**:
- From aMGSIM synthetic data: reads with damage in aa-damage.tsv.gz
- Ground truth: `damage_in_first_10bp OR damage_in_last_10bp`

**Negative examples (undamaged reads)**:
- Modern GTDB sequences (no damage simulation)
- Random subsequences from coding regions

### Training Data Sources

1. **KapK datasets** (2M years, high damage): ~60K coding reads each
   - 119_B3_116_L0_KapK-12-1-25
   - 119_B3_116_L0_KapK-12-1-24
   - etc.

2. **Mediterranean datasets** (varied ages): ~60K coding reads each
   - MED-2021-* series
   - MED-2022-* series

3. **GTDB undamaged** (negative class):
   - Sample ~500K coding sequences from gtdb_proteins_nt.tar.gz

### Training Protocol

```python
# Stratified split
train: 70%
val: 15%
test: 15%

# Class balancing
# Use class weights or SMOTE for imbalanced data

# Loss function
binary_crossentropy with label smoothing (0.1)

# Optimizer
Adam(lr=1e-3) with cosine decay

# Early stopping
patience=10, monitor='val_auc'

# Regularization
L2 weight decay (1e-4)
Dropout (0.2-0.3)
```

## C++ Integration

### Weight Export

```python
# After training, export weights to C++ header
def export_to_cpp_header(model, filename):
    with open(filename, 'w') as f:
        f.write("#pragma once\n")
        f.write("namespace agp {\n")

        for layer in model.layers:
            if hasattr(layer, 'kernel'):
                weights = layer.kernel.numpy()
                biases = layer.bias.numpy()
                # Write as constexpr arrays
                ...

        f.write("} // namespace agp\n")
```

### Inference Code

```cpp
// include/agp/damage_nn.hpp
class DamageNeuralNetwork {
public:
    float predict(const DamageFeatures& features) const {
        // Layer 1: 56 → 64
        auto h1 = relu(matmul(features, W1) + b1);

        // Layer 2: 64 → 32
        auto h2 = relu(matmul(h1, W2) + b2);

        // Layer 3: 32 → 16
        auto h3 = relu(matmul(h2, W3) + b3);

        // Output: 16 → 1
        float logit = dot(h3, W4) + b4;
        return sigmoid(logit);
    }

private:
    // Weights embedded as constexpr arrays
    static constexpr float W1[56][64] = { ... };
    static constexpr float b1[64] = { ... };
    // etc.
};
```

## Integration with Gene Prediction

### Pass 1: Sample Profile + NN Calibration

```cpp
// Compute sample-level damage profile
SampleDamageProfile profile = scan_for_damage(reads);

// NN uses sample profile as context features
DamageNeuralNetwork nn;
nn.set_sample_context(profile);
```

### Pass 2: Per-Read Damage Classification

```cpp
for (auto& read : reads) {
    // Extract features
    DamageFeatures features = extract_features(read, profile);

    // NN prediction
    float damage_prob = nn.predict(features);

    // Adjust frame selection
    if (damage_prob > 0.7) {
        // Use damage-aware frame selector
        // Reduce stop codon penalty for C→T/G→A contexts
    }

    gene.ancient_prob = damage_prob;
}
```

## Expected Performance

Based on feature analysis:

| Metric | Current | Target |
|--------|---------|--------|
| AUC-ROC | 0.39 | >0.85 |
| Precision@0.5 | 21% | >70% |
| Recall@0.5 | 21% | >80% |

## Next Steps

1. **Implement feature extraction** in Python for training
2. **Build training dataset** from synthetic data
3. **Train and validate** the NN
4. **Export weights** to C++ header
5. **Integrate** into AGP Pass 2
6. **Benchmark** on held-out datasets
