#!/bin/bash
# run AF3 predictions using both GPUs in parallel (no SLURM).
#
# splits 296 jobs into two halves:
#   GPU 0: jobs 1-148 (odd lines from job list)
#   GPU 1: jobs 149-296 (even lines from job list)
#
# usage:
#   cd /storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination
#   mkdir -p structures/af3_logs structures/af3_outputs
#   nohup bash scripts/12_run_af3_parallel.sh > structures/af3_logs/parallel_master.log 2>&1 &

set -euo pipefail

PROJECT_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
JOBLIST="${PROJECT_DIR}/structures/af3_jobs.txt"
OUTPUT_DIR="${PROJECT_DIR}/structures/af3_outputs"
LOG_DIR="${PROJECT_DIR}/structures/af3_logs"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# activate conda
source ~/miniforge3/etc/profile.d/conda.sh
conda activate af3_mmseqs2

# AF3 environment
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95
# mmseqs2 is in base conda, not af3_mmseqs2 env
export PATH="/home/kiran/miniforge3/bin:$PATH"

TOTAL_JOBS=$(wc -l < "$JOBLIST")
HALF=$(( (TOTAL_JOBS + 1) / 2 ))

echo "total jobs: $TOTAL_JOBS, split at: $HALF"
echo "GPU 0: jobs 1-$HALF"
echo "GPU 1: jobs $((HALF+1))-$TOTAL_JOBS"
echo "started: $(date)"

# function to run one job
run_one() {
    local json_path="$1"
    local gpu="$2"
    local job_num="$3"
    local job_name=$(basename "$json_path" .json)

    # skip if done
    if [ -f "${OUTPUT_DIR}/${job_name}/${job_name}_summary_confidences.json" ]; then
        echo "[GPU${gpu}] SKIP ${job_name} (already done)"
        return 0
    fi

    echo "[GPU${gpu}] START ${job_name} (job ${job_num}/${TOTAL_JOBS}) $(date)"

    CUDA_VISIBLE_DEVICES=$gpu python /storage/kiran-stuff/AF3_mmseqs2/run_alphafold.py \
        --json_path="$json_path" \
        --model_dir=/storage/kiran-stuff/alphafold3_models \
        --db_dir=/storage/kiran-stuff/alphafold_databases \
        --output_dir="$OUTPUT_DIR" \
        --buckets='256,512,768,1024,1536,2048' \
        --flash_attention_implementation=xla \
        > "${LOG_DIR}/${job_name}.log" 2>&1

    if [ $? -eq 0 ]; then
        echo "[GPU${gpu}] DONE  ${job_name} $(date)"
    else
        echo "[GPU${gpu}] FAIL  ${job_name} $(date)"
    fi
}

# run GPU 0 batch in background
(
    i=0
    while IFS= read -r json_path; do
        i=$((i + 1))
        if [ $i -le $HALF ]; then
            run_one "$json_path" 0 "$i"
        fi
    done < "$JOBLIST"
    echo "[GPU0] ALL DONE $(date)"
) &
PID_GPU0=$!

# run GPU 1 batch in background
(
    i=0
    while IFS= read -r json_path; do
        i=$((i + 1))
        if [ $i -gt $HALF ]; then
            run_one "$json_path" 1 "$i"
        fi
    done < "$JOBLIST"
    echo "[GPU1] ALL DONE $(date)"
) &
PID_GPU1=$!

echo "GPU 0 PID: $PID_GPU0"
echo "GPU 1 PID: $PID_GPU1"
echo "monitor: tail -f ${LOG_DIR}/parallel_master.log"

# wait for both
wait $PID_GPU0
wait $PID_GPU1

echo "========================================="
echo "ALL PREDICTIONS COMPLETE $(date)"
echo "========================================="

# count results
completed=$(find "$OUTPUT_DIR" -name "*_summary_confidences.json" | wc -l)
echo "completed: $completed / $TOTAL_JOBS"
