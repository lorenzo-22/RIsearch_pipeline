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

N_RUNS_LARGE=5
N_WARMUP=1
N_ITERATIONS=5

RESULTS_CSV="${NEW_PIPELINE_DIR}/benchmark_results.csv"
LOG_DIR="${NEW_PIPELINE_DIR}/benchmark_logs"
TIMEFILE=$(mktemp)

module load RIsearch2
module load bedtools

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

# Sets LAST_WALL_SEC and LAST_PEAK_KB globals; keeps log on failure.
run_timed() {
    local cmd="$1"
    local log_file="$2"
    set +e
    /usr/bin/time -v -o "${TIMEFILE}" bash -c "$cmd" > "$log_file" 2>&1
    local exit_code=$?
    set -e

    local wall_raw
    wall_raw=$(grep "Elapsed (wall clock)" "${TIMEFILE}" | sed 's/.*): //' | tr -d ' ')
    LAST_WALL_SEC=$(parse_wall_clock "$wall_raw")
    LAST_PEAK_KB=$(grep "Maximum resident set size" "${TIMEFILE}" | awk '{print $NF}')

    if [ "$exit_code" -ne 0 ]; then
        echo "    WARNING: command exited with code ${exit_code}. Log: $log_file"
    else
        rm -f "$log_file"
    fi
}

record_result() {
    printf '%s,%s,%s,%s,%s\n' "$1" "$2" "$3" "$4" "$5" >> "$RESULTS_CSV"
}

run_and_record_actual() {
    local pipeline="$1"
    local count="$2"
    local iter="$3"
    local cmd="$4"
    local log_file="${LOG_DIR}/${pipeline}_n${count}_iter${iter}.log"

    echo "    [iter ${iter}] ${pipeline} n=${count}..."
    run_timed "$cmd" "$log_file"
    record_result "$pipeline" "$count" "$iter" "$LAST_WALL_SEC" "$LAST_PEAK_KB"
}

# --- Prerequisites ---
echo "--- Checking prerequisites ---"
mkdir -p "$TMP_ACC" "$RAMDISK" "$SUBSET_DIR" "$LOG_DIR"
rm -f "${SCRATCH}"/*.temp.bed "${SCRATCH}"/*.targets "${SCRATCH}"/*.targets.*

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

if [ ! -s "$RESULTS_CSV" ]; then
    printf 'pipeline_version,siRNA_count,iteration_index,duration_seconds,peak_memory_kb\n' > "$RESULTS_CSV"
fi
echo "Results will be appended to: $RESULTS_CSV"

CONDA_SETUP="source /home/users/lorenzo/miniconda3/etc/profile.d/conda.sh && conda activate /home/users/lorenzo/.conda/envs/pitone2"
RISEARCH2_BIN="$(module show RIsearch2 2>&1 | awk '/prepend-path.*PATH/{print $3}')/risearch2.x"

# --- Generate master 2000 FASTA ---
MASTER_FA="${SUBSET_DIR}/master_2000.fa"
MASTER_IDS="${SUBSET_DIR}/master_2000.ids"
MASTER_ON_TARGET_MAP="${SUBSET_DIR}/master_2000.on_target_map.tsv"
MASTER_PARQUET="${SUBSET_DIR}/master_2000_parquet"

echo "--- Generating master subset of 2000 siRNAs ---"
paste - - < "$INPUT_FASTA" | shuf -n 2000 | tr '\t' '\n' > "$MASTER_FA"
grep "^>" "$MASTER_FA" | sed 's/^>//' > "$MASTER_IDS"
while IFS= read -r sid; do
    gene_id=$(echo "$sid" | cut -d'-' -f1)
    printf '%s\t%s\n' "$sid" "$gene_id"
done < "$MASTER_IDS" > "$MASTER_ON_TARGET_MAP"
echo "  IDs: $(wc -l < "$MASTER_IDS") siRNAs"

# --- Convert master 2000 to parquet (setup, not benchmarked) ---
echo "--- Converting 2000 siRNAs to parquet ---"
mkdir -p "$MASTER_PARQUET"
bash -c "cd ${NEW_PIPELINE_DIR} && \
uv run python3 convert_risearch_to_parquet.py ${RESULTS_DIR} \
    --ids-file ${MASTER_IDS} \
    --out-dir ${MASTER_PARQUET} \
    --workers 32" > /dev/null 2>&1
echo "  Parquet files: $(ls "$MASTER_PARQUET" | wc -l)/$(wc -l < "$MASTER_IDS")"

CMD_NEW_2000="cd ${NEW_PIPELINE_DIR} && \
mkdir -p ${SCRATCH}/new && \
uv run src/RIsearch_pipeline/cli.py off-targets \
-j 32 \
-r ${MASTER_PARQUET} \
-t ${RAMDISK}/E-MTAB-2770_fixed.bed \
--accessibility-dir ${TMP_ACC} \
--alpha '0.5;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
--gamma '0.55;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
--theta '0.5;0.6;0.7;0.75;0.8;0.85;0.9;0.95' \
-oi ${MASTER_ON_TARGET_MAP} \
-o ${SCRATCH}/new/output_2000.tsv \
--summary-only"

CMD_OLD_2000="${CONDA_SETUP} && \
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
-rx ${RISEARCH2_BIN} \
-q '{1}-{2}-{3}' \
-os ${MASTER_FA}\" :::: ${MASTER_IDS}"

# === TIER 1: 2000 siRNAs ===
echo ""
echo "===== BENCHMARK SIZE: 2000 siRNAs ====="

echo "--- New pipeline (2000) ---"
echo "  Warming up..."
for i in $(seq 1 $N_WARMUP); do
    pkill -f "spawn_main|forkserver" 2>/dev/null || true
    rm -rf "${SCRATCH}/new"
    bash -c "$CMD_NEW_2000" > /dev/null 2>&1 || true
    rm -rf "${SCRATCH}/new"
done

for run in $(seq 1 $N_RUNS_LARGE); do
    pkill -f "spawn_main|forkserver" 2>/dev/null || true
    run_and_record_actual "new" 2000 "$run" "$CMD_NEW_2000"
    rm -rf "${SCRATCH}/new"
done

echo "--- Old pipeline (2000) ---"
rm -f "${SCRATCH}"/*.targets "${SCRATCH}"/*.targets.* "${SCRATCH}"/*.temp.bed
echo "  Warming up..."
for i in $(seq 1 $N_WARMUP); do
    rm -rf "${SCRATCH}/old"
    bash -c "$CMD_OLD_2000" > /dev/null 2>&1 || true
    rm -rf "${SCRATCH}/old"
    rm -f "${SCRATCH}"/*.targets "${SCRATCH}"/*.targets.* "${SCRATCH}"/*.temp.bed
done

for run in $(seq 1 $N_RUNS_LARGE); do
    run_and_record_actual "old" 2000 "$run" "$CMD_OLD_2000"
    rm -rf "${SCRATCH}/old"
    rm -f "${SCRATCH}"/*.targets "${SCRATCH}"/*.targets.* "${SCRATCH}"/*.temp.bed
done

# Master parquet no longer needed after 2000 tier
rm -rf "$MASTER_PARQUET"

# === TIERS 2 & 3: 1000 and 500 siRNAs (5 random-subset iterations each) ===
for SUBSET_SIZE in 1000 500; do
    echo ""
    echo "===== BENCHMARK SIZE: ${SUBSET_SIZE} siRNAs ====="

    for iter in $(seq 1 $N_ITERATIONS); do
        echo "--- Iteration ${iter}/${N_ITERATIONS} ---"

        ITER_FA="${SUBSET_DIR}/iter_${SUBSET_SIZE}_${iter}.fa"
        ITER_IDS="${SUBSET_DIR}/iter_${SUBSET_SIZE}_${iter}.ids"
        ITER_ON_TARGET_MAP="${SUBSET_DIR}/iter_${SUBSET_SIZE}_${iter}.on_target_map.tsv"
        ITER_PARQUET="${SUBSET_DIR}/iter_${SUBSET_SIZE}_${iter}_parquet"

        paste - - < "$MASTER_FA" | shuf -n "$SUBSET_SIZE" | tr '\t' '\n' > "$ITER_FA"
        grep "^>" "$ITER_FA" | sed 's/^>//' > "$ITER_IDS"
        while IFS= read -r sid; do
            gene_id=$(echo "$sid" | cut -d'-' -f1)
            printf '%s\t%s\n' "$sid" "$gene_id"
        done < "$ITER_IDS" > "$ITER_ON_TARGET_MAP"

        mkdir -p "$ITER_PARQUET"
        bash -c "cd ${NEW_PIPELINE_DIR} && \
uv run python3 convert_risearch_to_parquet.py ${RESULTS_DIR} \
    --ids-file ${ITER_IDS} \
    --out-dir ${ITER_PARQUET} \
    --workers 32" > /dev/null 2>&1

        CMD_NEW_ITER="cd ${NEW_PIPELINE_DIR} && \
mkdir -p ${SCRATCH}/new && \
uv run src/RIsearch_pipeline/cli.py off-targets \
-j 32 \
-r ${ITER_PARQUET} \
-t ${RAMDISK}/E-MTAB-2770_fixed.bed \
--accessibility-dir ${TMP_ACC} \
--alpha '0.5;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
--gamma '0.55;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1' \
--theta '0.5;0.6;0.7;0.75;0.8;0.85;0.9;0.95' \
-oi ${ITER_ON_TARGET_MAP} \
-o ${SCRATCH}/new/output_${SUBSET_SIZE}_${iter}.tsv \
--summary-only"

        CMD_OLD_ITER="${CONDA_SETUP} && \
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
-rx ${RISEARCH2_BIN} \
-q '{1}-{2}-{3}' \
-os ${ITER_FA}\" :::: ${ITER_IDS}"

        # New pipeline
        echo "  [new] Warming up..."
        pkill -f "spawn_main|forkserver" 2>/dev/null || true
        rm -rf "${SCRATCH}/new"
        bash -c "$CMD_NEW_ITER" > /dev/null 2>&1 || true
        rm -rf "${SCRATCH}/new"

        pkill -f "spawn_main|forkserver" 2>/dev/null || true
        run_and_record_actual "new" "$SUBSET_SIZE" "$iter" "$CMD_NEW_ITER"
        rm -rf "${SCRATCH}/new"

        # Old pipeline
        rm -f "${SCRATCH}"/*.targets "${SCRATCH}"/*.targets.* "${SCRATCH}"/*.temp.bed
        echo "  [old] Warming up..."
        rm -rf "${SCRATCH}/old"
        bash -c "$CMD_OLD_ITER" > /dev/null 2>&1 || true
        rm -rf "${SCRATCH}/old"
        rm -f "${SCRATCH}"/*.targets "${SCRATCH}"/*.targets.* "${SCRATCH}"/*.temp.bed

        run_and_record_actual "old" "$SUBSET_SIZE" "$iter" "$CMD_OLD_ITER"
        rm -rf "${SCRATCH}/old"
        rm -f "${SCRATCH}"/*.targets "${SCRATCH}"/*.targets.* "${SCRATCH}"/*.temp.bed

        rm -f "$ITER_FA" "$ITER_IDS" "$ITER_ON_TARGET_MAP"
        rm -rf "$ITER_PARQUET"

        pkill -f "spawn_main|forkserver" 2>/dev/null || true
        echo "  Done. $(free -h | awk '/^Mem:/{print "Available: " $7}')"
    done
done

rm -f "$MASTER_FA" "$MASTER_IDS" "$MASTER_ON_TARGET_MAP" "$TIMEFILE"

echo ""
echo "===== ALL BENCHMARKS COMPLETE ====="
echo ""
echo "NOTE: peak_memory_kb reflects the parent process only."
echo "Multi-process workloads (GNU parallel / ProcessPoolExecutor) fork workers;"
echo "aggregate memory across all workers will exceed this figure."
echo ""
echo "Results: $RESULTS_CSV"
echo ""
column -t -s ',' "$RESULTS_CSV"
