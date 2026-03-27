#!/bin/bash -l
set -e  # Exit immediately if a command exits with a non-zero status

# --- Configuration ---
TMP_ACC="/tmp/acc_files"
RAMDISK="/dev/shm/input"
SUBSET_DIR="${RAMDISK}/subset_input"
SCRATCH="/var/tmp/tmp.ak9sAoblPN"
RESULTS_DIR="${SCRATCH}/RIsearch2_results"
NEW_PIPELINE_DIR="/dev/shm/src/RIsearch_pipeline"

# Input Source
INPUT_FASTA="${RAMDISK}/huesken_on.fa"

# CLI Argument for subset size (default to 10 if not provided)
SUBSET_SIZE=${1:-10}

module load RIsearch2
module load bedtools

# --- 1. Preparation ---
echo "--- Checking Pre-requisites ---"
mkdir -p "$TMP_ACC" "$RAMDISK" "$SUBSET_DIR"

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Source file $INPUT_FASTA not found."
    exit 1
fi

if ! command -v uv &> /dev/null; then
    echo "Error: 'uv' is not installed or not in PATH."
    exit 1
fi

if ! command -v hyperfine &> /dev/null; then
    echo "Error: hyperfine is not installed."
    exit 1
fi

# --- 2. Random Subsetting ---
echo "--- Creating random subset of $SUBSET_SIZE siRNAs ---"

# Files to be generated
SUBSET_FA="${SUBSET_DIR}/subset_${SUBSET_SIZE}.fa"
SUBSET_IDS="${SUBSET_DIR}/subset_${SUBSET_SIZE}.ids"
SUBSET_RESULTS="${SUBSET_DIR}/subset_${SUBSET_SIZE}_results"

# Step A: Subset the FASTA file
# 1. 'paste - -' pairs Header+Sequence onto one line
# 2. 'shuf' selects random lines
# 3. 'tr' restores the newline for FASTA format
paste - - < "$INPUT_FASTA" | shuf -n "$SUBSET_SIZE" | tr '\t' '\n' > "$SUBSET_FA"

# Step B: Generate the matching IDs file from the subset FASTA
# We extract the headers (removing '>') to create the ID list
grep "^>" "$SUBSET_FA" | sed 's/^>//' > "$SUBSET_IDS"

# Step C: Copy matching RIsearch result files to a subset directory
# The new pipeline takes a directory; this ensures both pipelines see the same inputs
rm -rf "$SUBSET_RESULTS" && mkdir -p "$SUBSET_RESULTS"
while IFS= read -r sid; do
    src="${RESULTS_DIR}/risearch_${sid}.out.gz"
    if [ -f "$src" ]; then
        cp "$src" "$SUBSET_RESULTS/"
    else
        echo "Warning: $src not found, skipping"
    fi
done < "$SUBSET_IDS"

# Step D: Build siRNA → on-target gene ID mapping (gene ID = first field of siRNA ID)
SUBSET_ON_TARGET_MAP="${SUBSET_DIR}/subset_${SUBSET_SIZE}.on_target_map.tsv"
while IFS= read -r sid; do
    gene_id=$(echo "$sid" | cut -d'-' -f1)
    printf '%s\t%s\n' "$sid" "$gene_id"
done < "$SUBSET_IDS" > "$SUBSET_ON_TARGET_MAP"

echo "Subset created:"
echo "  FASTA:   $SUBSET_FA"
echo "  IDs:     $SUBSET_IDS"
echo "  Results: $SUBSET_RESULTS ($(ls "$SUBSET_RESULTS" | wc -l) files)"
echo "  OT map:  $SUBSET_ON_TARGET_MAP"

TIME_LOG="${RAMDISK}/benchmark_metrics_${SUBSET_SIZE}.log"

# --- 3. Hyperfine Benchmark ---
echo "--- Starting Hyperfine Benchmark ---"

# --- COMMAND DEFINITIONS ---

# 1. OLD VERSION (Conda + Parallel + Python 2.7)
CONDA_SETUP="source /home/users/lorenzo/miniconda3/etc/profile.d/conda.sh && conda activate /home/users/lorenzo/.conda/envs/pitone2"

CMD_OLD="${CONDA_SETUP} && \
mkdir -p ${SCRATCH}/old && \
cd ${SCRATCH} && \
parallel -j 32 --colsep '-' \"python2.7 /dev/shm/src/pipeline.py \
-r '${RESULTS_DIR}/risearch_{1}-{2}-{3}.out.gz' \
-o 'old/{1}-{2}-{3}.H1299.off' \
-t ${RAMDISK}/E-MTAB-2770_fixed.bed \
-alpha '0.5;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
-gamma '0.55;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
-oi {1} --sort -p ${TMP_ACC} \
-theta '0.5;0.6;0.7;0.75;0.8;0.85;0.9;0.95' \
-rx risearch2.x \
-q '{1}-{2}-{3}' \
-os ${SUBSET_FA}\" :::: ${SUBSET_IDS}"

# 2. NEW VERSION (uv + Python 3, directory mode)
CMD_NEW="cd ${NEW_PIPELINE_DIR} && \
mkdir -p ${SCRATCH}/new && \
uv run src/RIsearch_pipeline/cli.py off-targets \
-r ${SUBSET_RESULTS} \
-t ${RAMDISK}/E-MTAB-2770_fixed.bed \
-a ${TMP_ACC} \
--alpha '0.5;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
--gamma '0.55;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
--theta '0.5;0.6;0.7;0.75;0.8;0.85;0.9;0.95' \
-oi ${SUBSET_ON_TARGET_MAP} \
-o ${SCRATCH}/new/output_${SUBSET_SIZE}.tsv"

# --- EXECUTION ---
hyperfine \
    --shell bash \
    --warmup 1 \
    --runs 3 \
    --show-output \
    --export-markdown benchmark_results_${SUBSET_SIZE}.md \
    -n "Old Version (Parallel+Python2.7)" \
    "$CMD_OLD" \
    -n "New Version (uv+Python3)" \
    "$CMD_NEW"

echo "--- Benchmark Complete ---"
echo "Results saved to: benchmark_results_${SUBSET_SIZE}.md"
