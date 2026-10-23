#!/bin/bash
# DART Tutorial: Ancient Gene Prediction + Damage Annotation
#
# Dataset: 50,000 synthetic ancient DNA reads from the KapK sediment benchmark
# (sample KapK-12-1-37, AT-rich community). predict prints the measured 5' C→T and 3' G→A rates.
#
# The included ref_proteins.faa.gz is the subset of reference proteins
# (13,681 of 149,276) that are actually matched by these reads.
#
# Prerequisites:
#   - dart built (cmake --build build -j8); pass its path as DART=, e.g. DART=../../build/dart
#   - mmseqs2 installed: conda install -c bioconda mmseqs2
#     (VTML20.out substitution matrix is bundled in examples/tutorial/data/)
#
# Runtime: about a minute with 8 threads, most of it in the two MMseqs2 searches
#
# Usage:
#   cd examples/tutorial
#   DB=./ref_proteins.faa.gz bash run_tutorial.sh [output_dir]
#
# To use your own protein database (e.g. KEGG, UniRef90):
#   DB=/path/to/mmseqs2_db bash run_tutorial.sh [output_dir]
#   (pre-index with: mmseqs createdb ref.faa ref_db && mmseqs createindex ref_db tmp)
#
# Standalone-run choices made explicit (damage-annotate and predict refuse to guess them):
#   --allow-self-profile  predict estimates the damage profile from these reads; production runs
#                         pass --library-profile from a library-level profile instead.
#   --decoy-every 10      predict also writes proteins.decoy.faa (reversed proteins of every
#                         10th read, ids prefixed DECOY_); the decoys are searched with the same
#                         MMseqs2 parameters and `hits2emi --decoy-hits` fits the null into
#                         hits.emi.null, which damage-annotate picks up by itself. Off by
#                         default: for a study, fit one null from a few samples and reuse it.
#   EM tape dir           an O_DIRECT-capable local dir for the EM tape; defaults to
#                         $OUT_DIR/em_tape, override with DART_EM_TAPE_DIR (not tmpfs).

set -euo pipefail

DART=${DART:-dart}
MMSEQS=${MMSEQS:-mmseqs}
DB=${DB:?Please set DB= to a protein FASTA or mmseqs2 database (e.g. DB=./ref_proteins.faa.gz)}
THREADS=${THREADS:-8}

TUTORIAL_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="${1:-$TUTORIAL_DIR/output}"
EM_TAPE_DIR="${DART_EM_TAPE_DIR:-$OUT_DIR/em_tape}"

mkdir -p "$OUT_DIR" "$OUT_DIR/tmp"

echo "=== DART Tutorial ==="
echo "Reads:     $TUTORIAL_DIR/reads.fq.gz (50,000 reads, KapK-12-1-37)"
echo "Reference: $DB"
echo "Output:    $OUT_DIR"
echo ""

# ---------------------------------------------------------------------------
# Step 1: Predict genes and build damage index
# ---------------------------------------------------------------------------
echo "Step 1/5: Gene prediction (+ decoy peptides)..."
"$DART" predict \
    -i "$TUTORIAL_DIR/reads.fq.gz" \
    -o "$OUT_DIR/predictions.gff" \
    --fasta-aa "$OUT_DIR/proteins.faa" \
    --damage-index "$OUT_DIR/predictions.agd" \
    --allow-self-profile \
    --decoy-every 10 \
    -t "$THREADS"
echo ""

# ---------------------------------------------------------------------------
# Step 2: Search predicted proteins against reference
# ---------------------------------------------------------------------------
echo "Step 2/5: MMseqs2 search against reference proteins..."
MMSEQS_ARGS=(
    --search-type 1
    --min-length 12
    -e 10.0
    --min-seq-id 0.86
    -c 0.65
    --cov-mode 2
    --sub-mat "$TUTORIAL_DIR/data/VTML20.out"
    --seed-sub-mat "$TUTORIAL_DIR/data/VTML20.out"
    -s 2 -k 6
    --spaced-kmer-pattern 11011101
    --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln"
    --threads "$THREADS"
)
"$MMSEQS" easy-search \
    "$OUT_DIR/proteins.faa" \
    "$DB" \
    "$OUT_DIR/hits.tsv" \
    "$OUT_DIR/tmp" \
    "${MMSEQS_ARGS[@]}"
echo ""

# ---------------------------------------------------------------------------
# Step 3: Search the decoy peptides with the same parameters (null model input)
# ---------------------------------------------------------------------------
echo "Step 3/5: MMseqs2 search of the decoy peptides..."
"$MMSEQS" easy-search \
    "$OUT_DIR/proteins.decoy.faa" \
    "$DB" \
    "$OUT_DIR/decoy_hits.tsv" \
    "$OUT_DIR/tmp" \
    "${MMSEQS_ARGS[@]}"
echo ""

# ---------------------------------------------------------------------------
# Step 4: Build columnar EMI index (required by damage-annotate) and fit the null
#         model from the decoy hits into hits.emi.null
# ---------------------------------------------------------------------------
echo "Step 4/5: Building EMI index and null model..."
"$DART" hits2emi \
    -i "$OUT_DIR/hits.tsv" \
    -o "$OUT_DIR/hits.emi" \
    --damage-index "$OUT_DIR/predictions.agd" \
    --decoy-hits "$OUT_DIR/decoy_hits.tsv" \
    --decoy-fasta "$OUT_DIR/proteins.decoy.faa" \
    -t "$THREADS"
echo ""

# ---------------------------------------------------------------------------
# Step 5: Damage annotation (finds hits.emi.null next to the EMI by itself)
# ---------------------------------------------------------------------------
echo "Step 5/5: Damage annotation..."
"$DART" damage-annotate \
    --emi "$OUT_DIR/hits.emi" \
    --damage-index "$OUT_DIR/predictions.agd" \
    -o "$OUT_DIR/annotated.tsv" \
    --gene-summary "$OUT_DIR/gene_summary.tsv" \
    --protein-summary "$OUT_DIR/protein_summary.tsv" \
    --auto-calibrate-spurious \
    --em-tape-dir "$EM_TAPE_DIR" \
    -t "$THREADS" \
    -v
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "=== Results ==="
GENES=$(grep -vc "^#" "$OUT_DIR/predictions.gff" 2>/dev/null || echo 0)
HITS=$(wc -l < "$OUT_DIR/hits.tsv")
DECOYS=$(grep -c "^>" "$OUT_DIR/proteins.decoy.faa")
DECOY_HITS=$(wc -l < "$OUT_DIR/decoy_hits.tsv")
ANN=$(awk 'NR>1' "$OUT_DIR/annotated.tsv" | wc -l)
DAMAGED_COL=$(
    awk -F'\t' '
        NR==1 {
            for (i=1; i<=NF; ++i) {
                if ($i == "is_damaged") {
                    print i;
                    exit;
                }
            }
        }
    ' "$OUT_DIR/annotated.tsv"
)
ANCIENT="n/a"
if [[ -n "${DAMAGED_COL}" ]]; then
    ANCIENT=$(
        awk -F'\t' -v c="$DAMAGED_COL" '
            NR>1 && $c == 1 { n++ }
            END { print n+0 }
        ' "$OUT_DIR/annotated.tsv"
    )
fi

echo "  Predicted genes:          $GENES"
echo "  MMseqs2 hits:             $HITS"
echo "  Decoy peptides / hits:    $DECOYS / $DECOY_HITS"
echo "  Damage-annotated reads:   $ANN"
echo "  Reads called damaged:     $ANCIENT (is_damaged = 1, posterior >= 0.7)"
echo ""
echo "Key output files:"
echo "  $OUT_DIR/predictions.gff     — gene predictions (GFF3)"
echo "  $OUT_DIR/proteins.faa        — predicted protein sequences"
echo "  $OUT_DIR/proteins.decoy.faa  — reversed decoy peptides (every 10th read)"
echo "  $OUT_DIR/hits.emi.null       — decoy null model (inspect with: dart damage-info hits.emi.null)"
echo "  $OUT_DIR/annotated.tsv       — per-read damage posteriors"
echo "  $OUT_DIR/gene_summary.tsv    — per-gene damage statistics"
echo "  $OUT_DIR/protein_summary.tsv — per-protein authentication (protein_auth_tier, auth_binom_call)"
echo ""
echo "Compare against expected outputs in: $TUTORIAL_DIR/expected/"
