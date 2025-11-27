#!/bin/bash
# Run AGP benchmarks on all samples
# Usage: ./scripts/run_benchmarks.sh [DATA_DIR] [OUT_DIR]

set -e

DATA_DIR="${1:-/projects/caeg/scratch/kbd606/agp/benchmark_data}"
OUT_DIR="${2:-/projects/caeg/scratch/kbd606/agp/benchmark_v3}"
AGP="$(dirname "$0")/../build/agp"

# Resolve AGP path if it's relative
AGP=$(cd "$(dirname "$AGP")" && pwd)/$(basename "$AGP")

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

mkdir -p "$OUT_DIR"

run_sample() {
    local name=$1
    local fq="${DATA_DIR}/${name}_art.fq.gz"
    local out_gff="${OUT_DIR}/${name}.gff"
    local out_faa="${OUT_DIR}/${name}.faa"
    local log="${OUT_DIR}/${name}.log"

    echo "[$(date +%H:%M:%S)] Starting: $name"
    "$AGP" "$fq" -o "$out_gff" --fasta-aa "$out_faa" --best-strand -v > "$log" 2>&1
    echo "[$(date +%H:%M:%S)] Completed: $name"
}

export -f run_sample
export DATA_DIR OUT_DIR AGP

# Run 4 samples in parallel
printf '%s\n' "${SAMPLES[@]}" | xargs -P 4 -I {} bash -c 'run_sample "$@"' _ {}

echo "All samples completed!"
