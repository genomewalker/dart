# Long Read Frame Error Analysis

## Executive Summary

Investigation of the 22.4% frame error rate revealed a **critical length-dependent failure mode**: reads longer than 130bp have near-random frame selection despite high prediction scores.

### Key Findings

| Read Length | Frame Accuracy | Cohen's d | Score Discrimination |
|-------------|---------------|-----------|---------------------|
| 30-75bp | 78% | 0.76 | Good |
| 75-100bp | 83% | 0.81 | Good |
| 100-130bp | 71% | 0.55 | Medium |
| **130-200bp** | **57%** | **0.15** | **Negligible** |

## Root Cause Analysis

### 1. Statistical Convergence Effect

As sequence length increases, the scoring signals (periodicity, RNY rule, amino acid composition, etc.) converge to similar values across all three reading frames due to the **Law of Large Numbers**:

- **Short reads**: High variance in averaged signals → clear frame discrimination
- **Long reads**: Low variance in averaged signals → frames become indistinguishable

Example: If rare codon frequency varies ±5% in short reads but only ±0.5% in long reads, the frame discrimination disappears.

### 2. Frame NN Length Saturation

The neural network's feature 20 (length normalization) saturates at 50 codons (150bp):
```cpp
features[20] = std::min(1.0f, static_cast<float>(n_codons) / 50.0f);
```

All reads >150bp have identical length features, so the NN cannot use length to adjust its predictions.

### 3. Stop Codon Dilution

For long sequences:
- The probability of having 0 internal stops in ALL 3 frames increases
- Stop penalty contributes 0% when all frames have equal stops (22% weight wasted)
- System relies entirely on weak secondary signals

### 4. Symmetric Error Pattern

Frame confusion matrix shows near-equal distribution:
- Shift +1: 49.4% of errors
- Shift +2: 46.5% of errors

This is essentially **random frame selection** among 3 options.

## Evidence

### Score Distribution Analysis

```
Length      n_correct    n_wrong  mean_corr  mean_wrng  diff    overlap%
30-75bp     1,418,587    172,151  0.330      0.124      +0.206  13.6%
75-100bp    237,387      17,167   0.404      0.180      +0.224  18.2%
100-130bp   141,360      14,681   0.422      0.263      +0.158  30.1%
130-200bp   44,882       26,810   0.427      0.383      +0.045  45.7%
```

For 130-200bp reads:
- Score difference between correct and wrong: **only 0.045** (vs 0.224 for 75-100bp)
- **45.7%** of wrong predictions have scores above the mean of correct predictions
- Cohen's d = 0.15 is classified as **negligible effect size**

## Proposed Improvements

### Option 1: Sliding Window Approach (Recommended)

Instead of averaging scores over the whole sequence, score in fixed-size windows (e.g., 60bp) and use **consistency** as the discrimination signal:

```
Intuition: The correct frame should score best in MOST windows
           Wrong frames may win some windows but inconsistently
```

**Implementation:**
1. Split read into overlapping 60bp windows
2. Score each window in all 3 frames
3. Count how often each frame wins
4. Use consistency metric: (best_frame_wins - second_best_wins) / total_windows

### Option 2: Length-Specific Scoring Weights

Different weight distributions for different length bins:

```cpp
if (length < 100) {
    // Current weights work well
} else if (length < 150) {
    // Increase dicodon weight, decrease stop weight
} else {
    // Heavy dicodon + unorthodox, minimal stop weight
}
```

### Option 3: Train Length-Aware NN

Current NN was trained on mixed-length data. Train separate models:
- Short NN: 30-100bp (current model)
- Long NN: 100-200bp (new model with length-specific features)

### Option 4: Stop Codon Position Analysis

Instead of counting stops, analyze their **relative positions**:
- Correct frame: stops tend to be at end or absent
- Wrong frames: stops distributed randomly throughout

Use first_stop_position / sequence_length as a feature with higher weight.

### Option 5: Ensemble with High-Confidence Filtering

For long reads where scoring fails:
1. Flag reads >130bp as "low confidence"
2. Optionally filter them from output
3. Or return all 3 frame predictions for downstream analysis

## Impact Assessment

Currently:
- 78,704 reads are 130-200bp (3.3% of dataset)
- 26,810 have frame errors (34% of all frame-only errors)
- Fixing this would reduce frame error rate by ~1-2 percentage points

However, the current 2D calibration already identifies these as low-confidence predictions, so users filtering by confidence are not significantly affected.

## Diagnostic Scripts Created

1. `.scripts/analyze_errors_from_csv.py` - Analyzes validation CSV by length/score bins
2. `.scripts/analyze_long_read_errors.py` - Deep dive into 130+bp error patterns
3. `.scripts/diagnose_long_read_frame_errors.py` - Frame confusion matrix analysis
4. `.scripts/analyze_score_discrimination.py` - Cohen's d analysis by length
5. `.scripts/analyze_stop_positions.py` - Stop codon distribution analysis
6. `.scripts/analyze_window_consistency.py` - Frame accuracy by length bin
7. `.scripts/analyze_long_read_failures.py` - Characteristics of failing predictions
8. `.scripts/analyze_stop_distribution.py` - Stop count analysis (correct vs wrong)
9. `.scripts/analyze_first_stop_position.py` - First stop position analysis
10. `.scripts/analyze_score_margin.py` - Score distribution analysis

## Implementation Attempts (Dec 2025)

### Sliding Window Consistency (Tested)

Implemented sliding window frame voting but found no improvement:
- Window size: 60-90bp with 30bp step
- Scores each window, counts which frame wins most
- **Result**: No improvement (56.1% → 56.1%)
- **Reason**: Same convergence problem exists within windows

### Damage-Aware Stop Scoring (Tested)

Fixed bug where `DamageProfile` was passed but unused. Now correctly extracts `delta_max`.
- **Result**: No improvement for low-damage samples (1.6% terminal damage)
- **Reason**: Damage-induced stops are rare in minimal-damage samples

### Key Diagnostic Finding: Stop Codon Presence

```
No-stop rate analysis for 130-200bp reads:
  Correct predictions: 74% have NO internal stops
  Wrong predictions (true frame): 21% have NO internal stops
  Wrong predictions (pred frame): 33% have NO internal stops
```

This reveals the core problem:
- When correct frame has NO stops (74% of cases), we get it right
- When correct frame HAS stops (26% of cases), we often pick wrong frame with fewer stops
- We can't distinguish "correct frame with artifact stops" from "wrong frame with fewer random stops"

### First Stop Position Analysis

Stop position distribution is **uniform** across all groups:
- Early (0-0.33): ~56% in all cases
- Middle (0.33-0.67): ~25% in all cases
- Late (0.67-1.0): ~17% in all cases

Stop position does NOT help discriminate correct from wrong frames.

### Score Calibration Problem

46% of wrong predictions have scores **above** the mean of correct predictions:
```
Correct: mean=0.579, range=[0.0005, 0.994]
Wrong:   mean=0.478, range=[0.0005, 0.994]

Wrong predictions with high confidence (>0.8): 24%
```

The scoring system gives false confidence to long read predictions.

## Conclusion

The long read frame error is a **fundamental statistical limitation** of averaging-based scoring signals. Multiple approaches were tested without success:

1. **Sliding window consistency**: Same convergence in windows
2. **Damage-aware scoring**: Minimal benefit for low-damage samples
3. **Stop position weighting**: Position doesn't discriminate
4. **Length-adaptive weights**: Already implemented, not sufficient

The ~56% accuracy for 130-200bp reads appears to be near the limit of what's achievable with current features. The remaining 26% of cases where the correct frame has internal stops are fundamentally difficult because we lack discriminative signals beyond stop codons.

**Recommended mitigation**: Use the 2D calibration system to flag these predictions as low-confidence, allowing downstream analysis to filter or process them appropriately. Users should be aware that frame predictions for reads >130bp have higher uncertainty.
