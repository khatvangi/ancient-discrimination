#!/bin/bash
# run AF3 pilot batch: 6 jobs split across 2 GPUs.
#
# GPU 0: S4 (RNA, DNA) + HHH (RNA)
# GPU 1: HHH (DNA) + S1 (RNA, DNA)
#
# usage:
#   cd /storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination
#   nohup bash scripts/12_run_af3_pilot.sh > structures/af3_logs/pilot_master.log 2>&1 &

set -euo pipefail

PROJECT_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
JOBLIST="${PROJECT_DIR}/structures/af3_pilot_jobs.txt"
OUTPUT_DIR="${PROJECT_DIR}/structures/af3_outputs"
LOG_DIR="${PROJECT_DIR}/structures/af3_logs"
AF3_PYTHON="/home/kiran/miniforge3/envs/af3_mmseqs2/bin/python"
AF3_SCRIPT="/storage/kiran-stuff/AF3_mmseqs2/run_alphafold.py"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# ensure mmseqs2 is on PATH
export PATH="/home/kiran/miniforge3/bin:$PATH"
# AF3 environment variables
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

TOTAL_JOBS=$(wc -l < "$JOBLIST")
echo "pilot batch: $TOTAL_JOBS jobs on 2 GPUs"
echo "started: $(date)"

run_one() {
    local json_path="$1"
    local gpu="$2"
    local job_name=$(basename "$json_path" .json)

    # skip if done
    if [ -f "${OUTPUT_DIR}/${job_name}/${job_name}_summary_confidences.json" ]; then
        echo "[GPU${gpu}] SKIP ${job_name} (already done)"
        return 0
    fi

    echo "[GPU${gpu}] START ${job_name} $(date)"

    CUDA_VISIBLE_DEVICES=$gpu $AF3_PYTHON $AF3_SCRIPT \
        --json_path="$json_path" \
        --model_dir=/storage/kiran-stuff/alphafold3_models \
        --db_dir=/storage/kiran-stuff/alphafold_databases \
        --output_dir="$OUTPUT_DIR" \
        --buckets='256,512,768,1024,1536,2048' \
        --flash_attention_implementation=xla \
        > "${LOG_DIR}/${job_name}.log" 2>&1

    local status=$?
    if [ $status -eq 0 ]; then
        echo "[GPU${gpu}] DONE  ${job_name} $(date)"
    else
        echo "[GPU${gpu}] FAIL  ${job_name} (exit $status) $(date)"
    fi
    return $status
}

# GPU 0: jobs 1, 2, 3 (S4 RNA, S4 DNA, HHH RNA)
(
    i=0
    while IFS= read -r json_path; do
        i=$((i + 1))
        if [ $i -le 3 ]; then
            run_one "$json_path" 0
        fi
    done < "$JOBLIST"
    echo "[GPU0] ALL DONE $(date)"
) &
PID_GPU0=$!

# GPU 1: jobs 4, 5, 6 (HHH DNA, S1 RNA, S1 DNA)
(
    i=0
    while IFS= read -r json_path; do
        i=$((i + 1))
        if [ $i -gt 3 ]; then
            run_one "$json_path" 1
        fi
    done < "$JOBLIST"
    echo "[GPU1] ALL DONE $(date)"
) &
PID_GPU1=$!

echo "GPU 0 PID: $PID_GPU0 (S4-RNA, S4-DNA, HHH-RNA)"
echo "GPU 1 PID: $PID_GPU1 (HHH-DNA, S1-RNA, S1-DNA)"

wait $PID_GPU0
wait $PID_GPU1

echo ""
echo "========================================="
echo "PILOT BATCH COMPLETE $(date)"
echo "========================================="

# check results
completed=$(find "$OUTPUT_DIR" -name "*_summary_confidences.json" 2>/dev/null | wc -l)
echo "completed: $completed / $TOTAL_JOBS"

# show confidence scores
echo ""
echo "PILOT RESULTS:"
for conf in $(find "$OUTPUT_DIR" -name "*_summary_confidences.json" 2>/dev/null | sort); do
    job=$(basename $(dirname "$conf"))
    iptm=$($AF3_PYTHON -c "import json; d=json.load(open('$conf')); print(f'{d.get(\"iptm\", 0):.3f}')")
    ptm=$($AF3_PYTHON -c "import json; d=json.load(open('$conf')); print(f'{d.get(\"ptm\", 0):.3f}')")
    echo "  $job: ipTM=$iptm pTM=$ptm"
done
