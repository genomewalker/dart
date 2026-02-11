# Hexamer Table Training Guide

This document describes how to generate hexamer frequency tables for AGP from reference CDS sequences.

## Overview

AGP uses three types of hexamer tables for gene prediction:

1. **Overall Hexamer Frequencies** (`{domain}_hexamer_table.hpp`)
   - General codon-pair (dicodon) frequencies across all CDS positions
   - Used for coding potential scoring

2. **Positional Hexamer Frequencies** (`{domain}_positional_hexamer.hpp`)
   - Position-specific frequencies for START (first 30bp), INTERNAL, and END (last 30bp) regions
   - Improves start/stop codon detection and truncation handling

3. **Damage Likelihood Ratios** (`{domain}_damage_likelihood.hpp`)
   - Log-likelihood ratios for damage-consistent hexamer patterns
   - 5' C→T damage: compares T-starting hexamers to their C-precursors
   - 3' G→A damage: compares A-ending hexamers to their G-precursors
   - Used for ancient DNA damage detection and correction

## Supported Domains

| Domain | Data Source | Description |
|--------|-------------|-------------|
| `gtdb` | GTDB r220 | Bacteria and Archaea (default for metagenomes) |
| `fungi` | RefSeq | Fungi |
| `protozoa` | RefSeq | Protozoa |
| `viral` | RefSeq | Viruses |
| `plant` | RefSeq | Plants |
| `invertebrate` | RefSeq | Invertebrates |
| `vertebrate_mammalian` | RefSeq | Mammalian vertebrates |
| `vertebrate_other` | RefSeq | Non-mammalian vertebrates |

## Prerequisites

### Software Requirements

```bash
# Build tools
cmake >= 3.18
g++ with C++20 support
zlib development headers

# Download tools
aria2c    # Parallel downloads (apt install aria2)
curl      # For fetching assembly summaries
```

### Build the Training Tools

```bash
cd /path/to/ancient-gene-prediction
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc) agp-train-hexamers
```

This builds `agp-train-hexamers`, the unified hexamer extraction tool.

## Step 1: Download Reference CDS Sequences

### Option A: GTDB (Bacteria/Archaea)

Download GTDB representative genomes from https://gtdb.ecogenomic.org/

```bash
# Create output directory
GTDB_DIR="/path/to/databases/gtdb"
mkdir -p ${GTDB_DIR}

# Download GTDB r220 protein_fna_reps (CDS nucleotide sequences)
# Bacteria (~85 GB compressed)
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/bac120_cds_reps.tar.gz
tar -xzf bac120_cds_reps.tar.gz -C ${GTDB_DIR}/bacteria/

# Archaea (~3 GB compressed)
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/ar53_cds_reps.tar.gz
tar -xzf ar53_cds_reps.tar.gz -C ${GTDB_DIR}/archaea/
```

### Option B: RefSeq (Eukaryotes and Viruses)

Use the provided download script or manual commands:

```bash
# Using the provided script
./scripts/download_refseq_cds.sh

# Or manually for a single domain:
DOMAIN="fungi"  # or: protozoa, viral, plant, invertebrate, vertebrate_mammalian, vertebrate_other
OUTDIR="/path/to/databases/refseq_cds/${DOMAIN}"
mkdir -p ${OUTDIR}

# Download assembly summary
curl -s "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/${DOMAIN}/assembly_summary.txt" \
  -o ${OUTDIR}/assembly_summary.txt

# Generate CDS download URLs
grep -v "^#" ${OUTDIR}/assembly_summary.txt | \
  awk -F'\t' '{
    ftp_path=$20;
    if (ftp_path != "na" && ftp_path != "") {
      n=split(ftp_path, parts, "/");
      asm_name=parts[n];
      print ftp_path "/" asm_name "_cds_from_genomic.fna.gz"
    }
  }' > ${OUTDIR}/urls_cds.txt

# Download with aria2 (parallel)
aria2c -x 16 -j 16 --auto-file-renaming=false --continue=true \
  --max-tries=3 -d ${OUTDIR}/cds -i ${OUTDIR}/urls_cds.txt
```

### Expected Data Sizes

| Domain | Assemblies | CDS Files | Disk Space |
|--------|------------|-----------|------------|
| GTDB (bacteria) | ~85,000 | ~85,000 | ~85 GB |
| GTDB (archaea) | ~4,000 | ~4,000 | ~3 GB |
| fungi | ~665 | ~665 | ~2 GB |
| protozoa | ~121 | ~121 | ~1 GB |
| viral | ~15,000 | ~15,000 | ~3 GB |
| plant | ~199 | ~199 | ~5 GB |
| invertebrate | ~486 | ~486 | ~8 GB |
| vertebrate_mammalian | ~250 | ~250 | ~15 GB |
| vertebrate_other | ~504 | ~504 | ~20 GB |

## Step 2: Generate Hexamer Tables

### Using agp-train-hexamers

The unified training tool generates all three hexamer table types in a single pass:

```bash
cd /path/to/ancient-gene-prediction

# Single domain (e.g., fungi)
./build/agp-train-hexamers \
  --domain fungi \
  /path/to/databases/refseq_cds/fungi/cds \
  --output-dir include/agp \
  --threads 32

# Combined domains (e.g., GTDB = bacteria + archaea)
./build/agp-train-hexamers \
  --domain gtdb \
  /path/to/databases/gtdb/bacteria \
  /path/to/databases/gtdb/archaea \
  --output-dir include/agp \
  --threads 64
```

### Output Files

For each domain, three files are generated:

```
include/agp/{domain}_hexamer_table.hpp        # Overall frequencies
include/agp/{domain}_positional_hexamer.hpp   # Position-specific (START/INTERNAL/END)
include/agp/{domain}_damage_likelihood.hpp    # Damage LLR tables
```

### Processing Details

The tool performs these steps:

1. **CDS Validation**
   - Requires valid start codon (ATG, GTG, or TTG)
   - Requires valid stop codon (TAA, TAG, or TGA)
   - Minimum length: 90bp (to have distinct START/INTERNAL/END regions)
   - Must be divisible by 3 (complete codons)

2. **Hexamer Extraction**
   - Extracts hexamers at codon boundaries (every 3bp, frame 0)
   - START region: first 30bp of CDS
   - END region: last 30bp before stop codon
   - INTERNAL region: everything between START and END

3. **Damage LLR Calculation**
   - For 5' damage: calculates log2(P(C-hexamer) / P(T-hexamer)) for T-starting hexamers
   - For 3' damage: calculates log2(P(G-hexamer) / P(A-hexamer)) for A-ending hexamers
   - Positive values indicate damage-consistent patterns
   - Only includes ratios > 1.5 (informative patterns)

4. **Frequency Normalization**
   - Converts raw counts to frequencies
   - Applies Laplace smoothing for overall frequencies
   - Generates combined weighted table (20% start + 60% internal + 20% end)

## Step 3: Rebuild AGP

After generating new hexamer tables:

```bash
cd /path/to/ancient-gene-prediction/build
cmake --build . --target agp -j$(nproc)
```

## Training Statistics

### Current Hexamer Table Counts

The following statistics are from the hexamer tables currently included with AGP:

| Domain | Files | Valid CDS | Total Hexamers |
|--------|------:|----------:|---------------:|
| `fungi` | 665 | 6.76M | 3.18B |
| `protozoa` | 121 | 1.16M | 0.60B |
| `viral` | 15,014 | 0.68M | 0.16B |
| `plant` | 199 | 9.16M | 4.24B |
| `invertebrate` | 486 | 12.85M | 8.43B |
| `vertebrate_mammalian` | 250 | 12.92M | 8.70B |
| `vertebrate_other` | 504 | 21.26M | 15.10B |
| **Total (RefSeq)** | **17,239** | **64.69M** | **40.41B** |

*Note: GTDB hexamer tables (for bacteria/archaea) are not yet included.*

### Example Output

Training output from GTDB (expected):

```
=== RESULTS ===
Files processed: 143,614
Total sequences: 487.2M
Valid CDS: 405.4M
Overall hexamers: 129.0B
START hexamers: 3.65B
INTERNAL hexamers: 121.7B
END hexamers: 3.65B

=== TERMINAL NUCLEOTIDE BIAS ===
5' end T/(T+C) by position:
  32.1% 31.8% 31.5% 31.3% 31.2% 31.1% 31.0% 31.0% 30.9% 30.9%
3' end A/(A+G) by position:
  33.2% 32.9% 32.6% 32.4% 32.3% 32.2% 32.1% 32.0% 32.0% 31.9%
```

The terminal nucleotide bias should be approximately equal (~30-35%) for undamaged reference data. Ancient samples will show elevated T/(T+C) at 5' and A/(A+G) at 3'.

## File Format Reference

### Hexamer Encoding

Hexamers are encoded as 12-bit integers (0-4095):
- A = 0, C = 1, G = 2, T = 3
- Code = (b0 << 10) | (b1 << 8) | (b2 << 6) | (b3 << 4) | (b4 << 2) | b5

Example: `ATGAAA` = (0<<10)|(3<<8)|(2<<6)|(0<<4)|(0<<2)|0 = 0 + 768 + 128 + 0 + 0 + 0 = 896

### Header File Structure

```cpp
// {domain}_hexamer_table.hpp
#pragma once
namespace agp {
    // Laplace-smoothed frequencies
    constexpr float {DOMAIN}_HEXAMER_FREQ[4096] = { ... };
}

// {domain}_positional_hexamer.hpp
#pragma once
namespace agp {
    constexpr float {DOMAIN}_START_HEXAMER_FREQ[4096] = { ... };
    constexpr float {DOMAIN}_INTERNAL_HEXAMER_FREQ[4096] = { ... };
    constexpr float {DOMAIN}_END_HEXAMER_FREQ[4096] = { ... };
    // Combined: 20% start + 60% internal + 20% end
    constexpr float {DOMAIN}_POSITIONAL_HEXAMER_FREQ[4096] = { ... };
}

// {domain}_damage_likelihood.hpp
#pragma once
namespace agp {
    // log2(P(C-precursor) / P(T-form)) for T-starting hexamers
    static constexpr float {DOMAIN}_DAMAGE_LLR_5PRIME[4096] = { ... };
    // log2(P(G-precursor) / P(A-form)) for A-ending hexamers
    static constexpr float {DOMAIN}_DAMAGE_LLR_3PRIME[4096] = { ... };
}
```

## Troubleshooting

### No hexamers extracted

Check that your CDS files:
- Are gzip-compressed FASTA format (`.fna.gz`, `.fasta.gz`, or `.fa.gz`)
- Contain valid CDS with start and stop codons
- Have sequences >= 90bp

### Missing download files

If NCBI downloads fail:
```bash
# Convert FTP to HTTPS (more reliable)
sed -i 's|ftp://ftp.ncbi.nlm.nih.gov|https://ftp.ncbi.nlm.nih.gov|' urls_cds.txt

# Retry with lower concurrency
aria2c -x 4 -j 4 --max-tries=10 --retry-wait=5 -d cds -i urls_cds.txt
```

### Out of memory

Reduce thread count:
```bash
./build/agp-train-hexamers --domain fungi /path/to/cds --threads 8
```

## Quick Reference

```bash
# Full workflow for a new domain
DOMAIN="fungi"
CDS_DIR="/path/to/databases/refseq_cds/${DOMAIN}/cds"

# 1. Download CDS data (see Step 1)

# 2. Generate hexamer tables
./build/agp-train-hexamers --domain ${DOMAIN} ${CDS_DIR} \
  --output-dir include/agp --threads 32

# 3. Rebuild AGP
cd build && cmake --build . -j$(nproc)

# 4. Test
./agp --domain ${DOMAIN} test.fq.gz -o test.gff -v
```
