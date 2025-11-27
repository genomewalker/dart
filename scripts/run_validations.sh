#!/bin/bash
# Validate AGP benchmark results against ground truth
# Usage: ./scripts/run_validations.sh [DATA_DIR] [OUT_DIR]

set -e

DATA_DIR="${1:-/projects/caeg/scratch/kbd606/agp/benchmark_data}"
OUT_DIR="${2:-/projects/caeg/scratch/kbd606/agp/benchmark_v3}"
VALIDATOR="$(dirname "$0")/../build/agp-validate"

# Resolve validator path if it's relative
VALIDATOR=$(cd "$(dirname "$VALIDATOR")" && pwd)/$(basename "$VALIDATOR")

# List of samples with matching reference files
SAMPLES=(
    "119_B3_116_L0_KapK-12-1-24"
    "119_B3_116_L0_KapK-12-1-25"
    "119_B3_116_L0_KapK-12-1-27"
    "119_B3_116_L0_KapK-12-1-29"
    "69_B2_100_L0_KapK-12-1-34"
    "MED-2021-16-ver15-2LFQY-210811_S15"
    "MED-2021-24-ver15-2LFQY-210811_S21"
    "MED-2021-25-ver15-2LFQY-210811_S22"
    "MED-2021-4-ver15-2LFQY-210811_S4"
    "MED-2021-5-ver15-2LFQY-210811_S5"
    "MED-2022-11_S9"
    "MED-2022-17_S2"
    "MED-2022-23_S3"
    "MED-2022-32_S4"
    "MED-2022-7_S7"
)

# Create results TSV header
RESULTS_FILE="${OUT_DIR}/validation_results.tsv"
echo -e "sample\tframe_accuracy\tstrand_accuracy\tcombined_accuracy\tseq_identity\tdamage_auc\tn_coding\tn_predicted" > "$RESULTS_FILE"

for name in "${SAMPLES[@]}"; do
    fq="${DATA_DIR}/${name}_art.fq.gz"
    ref="${DATA_DIR}/${name}_aa-damage.tsv.gz"
    gff="${OUT_DIR}/${name}.gff"
    faa="${OUT_DIR}/${name}.faa"
    log="${OUT_DIR}/${name}_validation.log"

    echo "[$(date +%H:%M:%S)] Validating: $name"

    # Run validator and capture output
    "$VALIDATOR" "$fq" "$ref" "$gff" "$faa" > "$log" 2>&1 || true

    # Extract metrics from log
    frame=$(grep "Correct frame:" "$log" | awk '{print $3}' | sed 's/.*(//' | sed 's/%).*//')
    strand=$(grep "Correct strand:" "$log" | awk '{print $3}' | sed 's/.*(//' | sed 's/%).*//')
    combined=$(grep "Both correct:" "$log" | awk '{print $3}' | sed 's/.*(//' | sed 's/%).*//')
    identity=$(grep "Average identity:" "$log" | awk '{print $3}' | sed 's/%//')
    auc=$(grep "AUC-ROC:" "$log" | awk '{print $2}')
    n_coding=$(grep "Reference coding reads:" "$log" | awk '{print $4}')
    n_predicted=$(grep "Predicted genes:" "$log" | awk '{print $3}')

    echo -e "${name}\t${frame}\t${strand}\t${combined}\t${identity}\t${auc}\t${n_coding}\t${n_predicted}" >> "$RESULTS_FILE"
    echo "[$(date +%H:%M:%S)] Completed: $name (Frame: ${frame}%, Strand: ${strand}%, AUC: ${auc})"
done

echo ""
echo "All validations completed!"
echo ""
column -t -s$'\t' "$RESULTS_FILE"
