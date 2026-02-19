#!/bin/bash
# AGP Tutorial: Ancient Gene Prediction + Damage Annotation
#
# Dataset: 50,000 synthetic ancient DNA reads from the KapK sediment benchmark
# (sample KapK-12-1-37, AT-rich community, ~25% 5' C→T damage, ~25% 3' G→A).
#
# The included ref_proteins.faa.gz is the subset of reference proteins
# (13,681 of 149,276) that are actually matched by these reads.
#
# Prerequisites:
#   - agp built: cmake --build build -j8 && cp build/agp /usr/local/bin/
#   - mmseqs2 installed: conda install -c bioconda mmseqs2 (includes VTML20.out)
#
# Runtime: ~5 s (predict) + ~10 s (MMseqs2 against tutorial DB) + <1 s (annotate)
#
# Usage:
#   cd examples/tutorial
#   DB=./ref_proteins.faa.gz bash run_tutorial.sh [output_dir]
#
# To use your own protein database (e.g. KEGG, UniRef90):
#   DB=/path/to/mmseqs2_db bash run_tutorial.sh [output_dir]
#   (pre-index with: mmseqs createdb ref.faa ref_db && mmseqs createindex ref_db tmp)

set -euo pipefail

AGP=${AGP:-agp}
MMSEQS=${MMSEQS:-mmseqs}
DB=${DB:?Please set DB= to a protein FASTA or mmseqs2 database (e.g. DB=./ref_proteins.faa.gz)}
THREADS=${THREADS:-8}

TUTORIAL_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="${1:-$TUTORIAL_DIR/output}"

mkdir -p "$OUT_DIR" "$OUT_DIR/tmp"

echo "=== AGP Tutorial ==="
echo "Reads:     $TUTORIAL_DIR/reads.fq.gz (50,000 reads, KapK-12-1-37)"
echo "Reference: $DB"
echo "Output:    $OUT_DIR"
echo ""

# ---------------------------------------------------------------------------
# Step 1: Predict genes and build damage index
# ---------------------------------------------------------------------------
echo "Step 1/4: Gene prediction..."
"$AGP" predict \
    -i "$TUTORIAL_DIR/reads.fq.gz" \
    -o "$OUT_DIR/predictions.gff" \
    --fasta-aa "$OUT_DIR/proteins.faa" \
    --damage-index "$OUT_DIR/predictions.agd" \
    --adaptive \
    -t "$THREADS"
echo ""

# ---------------------------------------------------------------------------
# Step 2: Search predicted proteins against reference
# ---------------------------------------------------------------------------
echo "Step 2/4: MMseqs2 search against reference proteins..."
"$MMSEQS" easy-search \
    "$OUT_DIR/proteins.faa" \
    "$DB" \
    "$OUT_DIR/hits.tsv" \
    "$OUT_DIR/tmp" \
    --search-type 1 \
    --min-length 12 \
    -e 10.0 \
    --min-seq-id 0.86 \
    -c 0.65 \
    --cov-mode 2 \
    --sub-mat VTML20.out \
    --seed-sub-mat VTML20.out \
    -s 2 -k 6 \
    --spaced-kmer-pattern 11011101 \
    --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln" \
    --threads "$THREADS"
echo ""

# ---------------------------------------------------------------------------
# Step 3: Build columnar EMI index (required by damage-annotate)
# ---------------------------------------------------------------------------
echo "Step 3/4: Building EMI index..."
"$AGP" hits2emi \
    -i "$OUT_DIR/hits.tsv" \
    -o "$OUT_DIR/hits.emi" \
    --damage-index "$OUT_DIR/predictions.agd" \
    -t "$THREADS"
echo ""

# ---------------------------------------------------------------------------
# Step 4: Damage annotation
# ---------------------------------------------------------------------------
echo "Step 4/4: Damage annotation..."
"$AGP" damage-annotate \
    --emi "$OUT_DIR/hits.emi" \
    --damage-index "$OUT_DIR/predictions.agd" \
    -o "$OUT_DIR/annotated.tsv" \
    --gene-summary "$OUT_DIR/gene_summary.tsv" \
    --auto-calibrate-spurious \
    -t "$THREADS" \
    -v
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "=== Results ==="
GENES=$(grep -vc "^#" "$OUT_DIR/predictions.gff" 2>/dev/null || echo 0)
HITS=$(wc -l < "$OUT_DIR/hits.tsv")
ANN=$(awk 'NR>1' "$OUT_DIR/annotated.tsv" | wc -l)
POSTERIOR_COL=$(
    awk -F'\t' '
        NR==1 {
            for (i=1; i<=NF; ++i) {
                if ($i == "posterior") {
                    print i;
                    exit;
                }
            }
        }
    ' "$OUT_DIR/annotated.tsv"
)
ANCIENT="n/a"
if [[ -n "${POSTERIOR_COL}" ]]; then
    ANCIENT=$(
        awk -F'\t' -v c="$POSTERIOR_COL" '
            NR>1 && $c != "NA" && $c+0 >= 0.5 { n++ }
            END { print n+0 }
        ' "$OUT_DIR/annotated.tsv"
    )
fi

echo "  Predicted genes:          $GENES"
echo "  MMseqs2 hits:             $HITS"
echo "  Damage-annotated reads:   $ANN"
echo "  Ancient (posterior≥0.5):  $ANCIENT"
echo ""
echo "Key output files:"
echo "  $OUT_DIR/predictions.gff     — gene predictions (GFF3)"
echo "  $OUT_DIR/proteins.faa        — predicted protein sequences"
echo "  $OUT_DIR/annotated.tsv       — per-read damage posteriors"
echo "  $OUT_DIR/gene_summary.tsv    — per-gene damage statistics"
echo "  $TUTORIAL_DIR/data/benchmark_summary.tsv  — benchmark summary fixture for docs/tests"
echo ""
echo "Compare against expected outputs in: $TUTORIAL_DIR/expected/"
