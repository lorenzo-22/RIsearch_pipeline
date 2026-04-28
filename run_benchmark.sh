#!/bin/bash -l
set -e

# --- Configuration ---
TMP_ACC="/tmp/acc_files"
RAMDISK="/dev/shm/input"
SUBSET_DIR="/var/tmp/subset_input"
SCRATCH="/var/tmp/results"
RESULTS_DIR="${SCRATCH}/RIsearch2_results"
NEW_PIPELINE_DIR="/dev/shm/src/RIsearch_pipeline"

INPUT_FASTA="${RAMDISK}/huesken_on.fa"

if [ $# -gt 0 ]; then
    SUBSET_SIZES=("$@")
else
    SUBSET_SIZES=(500 1000 2000)
fi
N_RUNS=3
N_WARMUP=1

RESULTS_TSV="${NEW_PIPELINE_DIR}/benchmark_results.tsv"
PERF_LOG="${NEW_PIPELINE_DIR}/benchmark_perf.log"
LOG_DIR="${NEW_PIPELINE_DIR}/benchmark_logs"
TIMEFILE=$(mktemp)

module load RIsearch2
module load bedtools

# --- Helper: convert /usr/bin/time -v wall-clock string to seconds ---
# Accepts h:mm:ss or m:ss.ss (output format depends on GNU time version)
parse_wall_clock() {
    local raw="$1"
    if [[ "$raw" =~ ^([0-9]+):([0-9]+):([0-9.]+)$ ]]; then
        echo "$((${BASH_REMATCH[1]} * 3600 + ${BASH_REMATCH[2]} * 60))$(echo "+${BASH_REMATCH[3]}" | bc)" | bc
    elif [[ "$raw" =~ ^([0-9]+):([0-9.]+)$ ]]; then
        echo "${BASH_REMATCH[1]} * 60 + ${BASH_REMATCH[2]}" | bc
    else
        echo "0"
    fi
}

# --- Helper: run command under /usr/bin/time -v, record results to TSV ---
# Captures stdout/stderr to a log file; keeps it only on failure.
# Usage: run_and_record <subset_size> <pipeline_label> <run_number> <command_string>
run_and_record() {
    local subset_size="$1"
    local pipeline="$2"
    local run_num="$3"
    local cmd="$4"
    local log_file="${LOG_DIR}/${pipeline}_n${subset_size}_run${run_num}.log"

    echo "    [run ${run_num}] ${pipeline} n=${subset_size}..."
    set +e
    /usr/bin/time -v -o "${TIMEFILE}" bash -c "$cmd" > "$log_file" 2>&1
    local exit_code=$?
    set -e

    local wall_raw
    wall_raw=$(grep "Elapsed (wall clock)" "${TIMEFILE}" | sed 's/.*): //' | tr -d ' ')
    local wall_sec
    wall_sec=$(parse_wall_clock "$wall_raw")
    local peak_rss
    peak_rss=$(grep "Maximum resident set size" "${TIMEFILE}" | awk '{print $NF}')

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$subset_size" "$pipeline" "$run_num" "$wall_sec" "$peak_rss" "$exit_code" \
        >> "$RESULTS_TSV"

    # Extract [perf] block (new pipeline emits it at end of run)
    if grep -q '^\[perf\]' "$log_file" 2>/dev/null; then
        echo "    --- [perf] ${pipeline} n=${subset_size} run${run_num} ---"
        grep -A 999 '^\[perf\]' "$log_file" | head -30
        {
            printf '\n=== %s n=%s run%s ===\n' "$pipeline" "$subset_size" "$run_num"
            grep -A 999 '^\[perf\]' "$log_file" | head -30
        } >> "$PERF_LOG"
    fi

    if [ "$exit_code" -eq 0 ]; then
        rm -f "$log_file"
    else
        echo "    WARNING: ${pipeline} run ${run_num} exited with code ${exit_code}"
        echo "    Log kept: $log_file"
    fi
}

# --- Prerequisites ---
echo "--- Checking prerequisites ---"
mkdir -p "$TMP_ACC" "$RAMDISK" "$SUBSET_DIR" "$LOG_DIR"

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: $INPUT_FASTA not found."
    exit 1
fi

if ! command -v uv &> /dev/null; then
    echo "Error: 'uv' not in PATH."
    exit 1
fi

if ! /usr/bin/time --version &> /dev/null 2>&1; then
    echo "Error: /usr/bin/time not available (need GNU time for -v flag)."
    exit 1
fi

# --- Write TSV header only if file is empty/missing (append-friendly) ---
if [ ! -s "$RESULTS_TSV" ]; then
    printf 'subset_size\tpipeline\trun\twall_clock_s\tpeak_rss_kb\texit_code\n' > "$RESULTS_TSV"
fi
printf '# benchmark run: %s  sizes: %s\n' "$(date -Iseconds)" "${SUBSET_SIZES[*]}" >> "$PERF_LOG"
echo "Results will be appended to: $RESULTS_TSV"

# --- Conda setup string (reused in old pipeline command) ---
CONDA_SETUP="source /home/users/lorenzo/miniconda3/etc/profile.d/conda.sh && conda activate /home/users/lorenzo/.conda/envs/pitone2"

# --- Generate master FASTA at max size; derive smaller subsets as prefixes ---
# Nested subsets (500 ⊂ 1000 ⊂ 2000) ensure per-siRNA workload is identical
# across sizes, so cross-size comparisons reflect only siRNA count.
MAX_SIZE=$(printf '%s\n' "${SUBSET_SIZES[@]}" | sort -n | tail -1)
MASTER_FA="${SUBSET_DIR}/master_${MAX_SIZE}.fa"
echo "--- Generating master subset of ${MAX_SIZE} siRNAs (nested) ---"
paste - - < "$INPUT_FASTA" | shuf -n "$MAX_SIZE" | tr '\t' '\n' > "$MASTER_FA"

# --- Main loop ---
for SUBSET_SIZE in "${SUBSET_SIZES[@]}"; do
    echo ""
    echo "===== BENCHMARK SIZE: ${SUBSET_SIZE} siRNAs ====="

    SUBSET_FA="${SUBSET_DIR}/subset_${SUBSET_SIZE}.fa"
    SUBSET_IDS="${SUBSET_DIR}/subset_${SUBSET_SIZE}.ids"
    SUBSET_PARQUET="${SUBSET_DIR}/subset_${SUBSET_SIZE}_parquet"
    SUBSET_ON_TARGET_MAP="${SUBSET_DIR}/subset_${SUBSET_SIZE}.on_target_map.tsv"

    # Derive nested subset as prefix of master (each siRNA = 2 FASTA lines)
    echo "--- Deriving subset of $SUBSET_SIZE siRNAs (prefix of master) ---"
    head -n $((SUBSET_SIZE * 2)) "$MASTER_FA" > "$SUBSET_FA"
    grep "^>" "$SUBSET_FA" | sed 's/^>//' > "$SUBSET_IDS"

    while IFS= read -r sid; do
        gene_id=$(echo "$sid" | cut -d'-' -f1)
        printf '%s\t%s\n' "$sid" "$gene_id"
    done < "$SUBSET_IDS" > "$SUBSET_ON_TARGET_MAP"

    echo "  IDs: $(wc -l < "$SUBSET_IDS") siRNAs"

    # --- Old pipeline (3 runs) ---
    echo "--- Old pipeline ---"
    # Clean up any stale intersection files to ensure a fair benchmark and avoid corruption
    rm -f "${SCRATCH}"/*.targets "${SCRATCH}"/*.targets.* "${SCRATCH}"/*.temp.bed

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

    echo "  Warming up (${N_WARMUP} runs)..."
    for run in $(seq 1 $N_WARMUP); do
        rm -rf "${SCRATCH}/old"
        bash -c "$CMD_OLD" > /dev/null 2>&1 || true
    done

    for run in $(seq 1 $N_RUNS); do
        rm -rf "${SCRATCH}/old"
        run_and_record "$SUBSET_SIZE" "old" "$run" "$CMD_OLD"
    done
    rm -rf "${SCRATCH}/old"

    # --- Conversion: .out.gz → .parquet (1 run, reported separately) ---
    echo "--- Conversion (.out.gz → .parquet) ---"
    mkdir -p "$SUBSET_PARQUET"
    CMD_CONVERT="cd ${NEW_PIPELINE_DIR} && \
uv run python3 convert_risearch_to_parquet.py ${RESULTS_DIR} \
    --ids-file ${SUBSET_IDS} \
    --out-dir ${SUBSET_PARQUET} \
    --workers 32"

    run_and_record "$SUBSET_SIZE" "convert" "1" "$CMD_CONVERT"
    EXPECTED=$(wc -l < "$SUBSET_IDS")
    ACTUAL=$(ls "$SUBSET_PARQUET" 2>/dev/null | wc -l)
    echo "  Parquet files: ${ACTUAL}/${EXPECTED}"
    if [ "$ACTUAL" -lt "$EXPECTED" ]; then
        echo "  WARNING: only ${ACTUAL}/${EXPECTED} parquet files produced — possible disk pressure"
    fi

    # --- New pipeline on parquet (3 runs) ---
    echo "--- New pipeline (parquet input) ---"
    CMD_NEW="cd ${NEW_PIPELINE_DIR} && \
mkdir -p ${SCRATCH}/new && \
uv run src/RIsearch_pipeline/cli.py off-targets \
-j 32 \
-r ${SUBSET_PARQUET} \
-t ${RAMDISK}/E-MTAB-2770_fixed.bed \
--accessibility-dir ${TMP_ACC} \
--alpha '0.5;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
--gamma '0.55;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
--theta '0.5;0.6;0.7;0.75;0.8;0.85;0.9;0.95' \
-oi ${SUBSET_ON_TARGET_MAP} \
-o ${SCRATCH}/new/output_${SUBSET_SIZE}.tsv \
--summary-only"

    echo "  Warming up (${N_WARMUP} runs)..."
    for run in $(seq 1 $N_WARMUP); do
        pkill -f "spawn_main|forkserver" 2>/dev/null || true
        rm -rf "${SCRATCH}/new"
        bash -c "$CMD_NEW" > /dev/null 2>&1 || true
    done

    for run in $(seq 1 $N_RUNS); do
        pkill -f "spawn_main|forkserver" 2>/dev/null || true
        rm -rf "${SCRATCH}/new"
        run_and_record "$SUBSET_SIZE" "new" "$run" "$CMD_NEW"
    done
    rm -rf "${SCRATCH}/new"

    # --- Cleanup: delete intermediates and kill orphaned workers ---
    echo "--- Cleanup for size ${SUBSET_SIZE} ---"
    rm -rf "$SUBSET_PARQUET"
    rm -f "$SUBSET_FA" "$SUBSET_IDS" "$SUBSET_ON_TARGET_MAP"
    pkill -f "spawn_main|forkserver" 2>/dev/null || true
    sleep 2
    echo "  Done. $(free -h | awk '/^Mem:/{print "Available: " $7}')"
done

rm -f "$MASTER_FA" "$TIMEFILE"

echo ""
echo "===== ALL BENCHMARKS COMPLETE ====="
echo ""
echo "NOTE: peak_rss_kb is the peak RSS of the parent process only."
echo "For multi-process workloads (GNU parallel / ProcessPoolExecutor),"
echo "this reflects the parent shell, not the aggregate of all workers."
echo ""
echo "Results: $RESULTS_TSV"
echo ""
column -t -s $'\t' "$RESULTS_TSV"
