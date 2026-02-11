# META Mode: Multi-Domain Ensemble Scoring

META mode is AGP's default domain mode, designed for ancient metagenomes where the taxonomic origin of sequences is unknown. Rather than assuming all sequences come from a single taxonomic group, META mode uses ensemble-weighted scoring across all trained domains.

## Overview

AGP has pre-trained hexamer frequency tables for 8 taxonomic domains:

| Domain | Description | Source |
|--------|-------------|--------|
| GTDB | Bacteria + Archaea | GTDB representative genomes |
| Fungi | Fungal eukaryotes | RefSeq fungi |
| Protozoa | Single-celled eukaryotes | RefSeq protozoa |
| Invertebrate | Invertebrate animals | RefSeq invertebrate |
| Plant | Land plants | RefSeq plant |
| Vertebrate (Mammalian) | Mammals | RefSeq vertebrate_mammalian |
| Vertebrate (Other) | Non-mammalian vertebrates | RefSeq vertebrate_other |
| Viral | Viruses | RefSeq viral |

## How META Mode Works

### 1. Per-Sequence Domain Probability Estimation

For each input sequence, AGP computes domain probabilities using a softmax over hexamer log-likelihood scores:

```
For each domain d in {GTDB, Fungi, ...}:
    score[d] = Σ log₂(hexamer_freq[d] / random_freq)

Normalize via softmax with temperature scaling:
    prob[d] = exp((score[d] - max_score) * T) / Σ exp(...)
```

This produces a probability distribution like:
- GTDB: 0.65
- Fungi: 0.15
- Viral: 0.12
- Others: 0.08

### 2. Ensemble-Weighted Lookups

When ensemble mode is enabled, all hexamer-based lookups use weighted averages:

**Hexamer Frequency (coding potential)**:
```cpp
weighted_freq = gtdb_prob * GTDB_HEXAMER_FREQ[code] +
                fungi_prob * FUNGI_HEXAMER_FREQ[code] +
                protozoa_prob * PROTOZOA_HEXAMER_FREQ[code] +
                ... (all 8 domains)
```

**Positional Hexamers (START/INTERNAL/END regions)**:
```cpp
weighted_start_freq = gtdb_prob * GTDB_START_HEXAMER_FREQ[code] +
                      fungi_prob * FUNGI_START_HEXAMER_FREQ[code] +
                      ... (all 8 domains)
```

**Damage Likelihood (ancient DNA patterns)**:
```cpp
weighted_damage_llr = gtdb_prob * GTDB_DAMAGE_LLR_5PRIME[code] +
                      fungi_prob * FUNGI_DAMAGE_LLR_5PRIME[code] +
                      ... (all 8 domains)
```

### 3. Subsystems with Ensemble Support

| Subsystem | File | Ensemble Function |
|-----------|------|-------------------|
| Hexamer frequency | `hexamer_tables.hpp` | `get_ensemble_hexamer_freq()` |
| Dicodon scoring | `hexamer_tables.hpp` | `calculate_dicodon_score_weighted()` |
| Positional hexamers | `positional_hexamer.hpp` | `get_ensemble_start_hexamer_freq()`, etc. |
| Damage likelihood | `damage_likelihood.hpp` | `get_ensemble_damage_llr_5prime()`, etc. |

## Usage

META mode is the default. Explicitly select it with:

```bash
# Default (META mode)
agp input.fq.gz -o output.gff

# Explicit META mode
agp input.fq.gz -o output.gff --domain meta

# Alternative names for META mode
agp input.fq.gz -o output.gff --domain metagenome
agp input.fq.gz -o output.gff --domain all
```

For samples with known taxonomy, use a specific domain:

```bash
# Bacteria-only sample
agp input.fq.gz -o output.gff --domain gtdb

# Known fungal sample
agp input.fq.gz -o output.gff --domain fungi
```

## Implementation Details

### Thread-Local State

Ensemble mode uses thread-local storage for domain probabilities:

```cpp
// Thread-local ensemble state
inline bool& ensemble_mode_enabled() {
    static thread_local bool enabled = false;
    return enabled;
}

inline MultiDomainResult& current_domain_probs() {
    static thread_local MultiDomainResult result;
    return result;
}
```

This allows parallel processing where each thread maintains its own per-sequence domain weights.

### MultiDomainResult Structure

```cpp
struct MultiDomainResult {
    float gtdb_prob = 0.0f;
    float fungi_prob = 0.0f;
    float protozoa_prob = 0.0f;
    float invertebrate_prob = 0.0f;
    float plant_prob = 0.0f;
    float vertebrate_mammalian_prob = 0.0f;
    float vertebrate_other_prob = 0.0f;
    float viral_prob = 0.0f;
    Domain best_domain = Domain::GTDB;
    float best_score = 0.0f;
};
```

### Fallback Behavior

When ensemble mode is disabled (e.g., using `--domain gtdb`), ensemble functions fall back to single-domain lookups:

```cpp
inline float get_ensemble_hexamer_freq(uint32_t code) {
    if (!ensemble_mode_enabled()) {
        return get_hexamer_freq(code, get_active_domain());
    }
    // ... weighted calculation
}
```

## Performance Considerations

- **META mode is ~30% slower** than single-domain mode due to weighted calculations
- Domain probability estimation happens once per sequence, not per hexamer
- Pre-computed lookup tables (4096 entries each) ensure fast access

## When to Use Specific Domains

| Scenario | Recommended Mode |
|----------|-----------------|
| Ancient metagenomes (mixed community) | `meta` (default) |
| Ancient human/mammalian samples | `mammal` |
| Gut microbiome | `gtdb` |
| Environmental samples | `meta` |
| Pure bacterial culture | `gtdb` |
| Known fungal pathogen | `fungi` |

## Technical Notes

### Hexamer Encoding

All hexamers are encoded as 12-bit integers (0-4095):
```cpp
// A=0, C=1, G=2, T=3
// ATGCAA = 0*4^5 + 3*4^4 + 2*4^3 + 1*4^2 + 0*4^1 + 0*4^0 = 872
uint32_t encode_hexamer(const char* seq);
```

### Temperature Scaling

The softmax uses temperature scaling to control confidence:
- Temperature 0.1 in `score_all_domains()` (broad distribution)
- Temperature 0.5 in `score_positional_all_domains()` (sharper)

Lower temperature = more confident (winner-take-all), higher = more uniform.
