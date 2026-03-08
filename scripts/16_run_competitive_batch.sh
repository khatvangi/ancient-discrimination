#!/bin/bash
# NOTE: Phase 3 script (planned, not yet completed).
# competitive AF3 batch — not yet completed.
# run 93 competitive AF3 jobs on 2 GPUs (proteins <= 700 aa)
# uses standard AF3 (Docker + HMMER)
# ~30 min per job × 47 jobs per GPU ≈ 24 hours total

set -uo pipefail

INPUT_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/af3_competitive_local_batch"
OUTPUT_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/af3_competitive_output"
LOG_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/af3_logs/competitive_batch"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

run_one() {
    local json_file="$1"
    local gpu="$2"
    local job_name=$(basename "$json_file" .json)
    local log_file="${LOG_DIR}/${job_name}.log"

    # skip if already completed
    if [ -f "${OUTPUT_DIR}/${job_name}/${job_name}_summary_confidences.json" ]; then
        echo "[GPU${gpu}] SKIP ${job_name} (already done)"
        return 0
    fi

    echo "[GPU${gpu}] START ${job_name} $(date)"

    docker run --rm \
        --gpus '"device='"${gpu}"'"' \
        -e CUDA_VISIBLE_DEVICES=0 \
        -e XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter" \
        -e XLA_PYTHON_CLIENT_PREALLOCATE=false \
        -v /storage/kiran-stuff/alphafold3_models:/models:ro \
        -v /storage/kiran-stuff/alphafold_databases:/databases:ro \
        -v "${INPUT_DIR}":/inputs:ro \
        -v "${OUTPUT_DIR}":/output \
        alphafold3-local \
        python /app/alphafold/run_alphafold.py \
            --json_path="/inputs/${job_name}.json" \
            --model_dir=/models \
            --db_dir=/databases \
            --output_dir=/output \
            --flash_attention_implementation=xla \
            --jackhmmer_n_cpu=8 \
            --nhmmer_n_cpu=8 \
            --buckets='256,512,768,1024' \
        > "${log_file}" 2>&1

    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "[GPU${gpu}] DONE  ${job_name} (exit 0) $(date)"
    else
        echo "[GPU${gpu}] FAIL  ${job_name} (exit ${exit_code}) $(date)"
    fi
    return $exit_code
}

# collect all job names, sorted by protein size (smallest first = fastest first)
# this way we get early results quickly
ALL_JOBS=($(ls "$INPUT_DIR"/*.json | sort))
TOTAL=${#ALL_JOBS[@]}

echo "competitive AF3 local batch: ${TOTAL} jobs on 2 GPUs"
echo "started: $(date)"
echo "input:   ${INPUT_DIR}"
echo "output:  ${OUTPUT_DIR}"
echo "logs:    ${LOG_DIR}"
echo ""

# split odd/even into two GPU queues
GPU0_JOBS=()
GPU1_JOBS=()
for i in "${!ALL_JOBS[@]}"; do
    if (( i % 2 == 0 )); then
        GPU0_JOBS+=("${ALL_JOBS[$i]}")
    else
        GPU1_JOBS+=("${ALL_JOBS[$i]}")
    fi
done

echo "GPU 0: ${#GPU0_JOBS[@]} jobs"
echo "GPU 1: ${#GPU1_JOBS[@]} jobs"
echo ""

# GPU 0
(
    completed=0
    failed=0
    for json_file in "${GPU0_JOBS[@]}"; do
        job_name=$(basename "$json_file" .json)
        run_one "$job_name" 0
        if [ $? -eq 0 ]; then ((completed++)); else ((failed++)); fi
    done
    echo "[GPU0] FINISHED: ${completed} done, ${failed} failed $(date)"
) &
GPU0_PID=$!

# GPU 1
(
    completed=0
    failed=0
    for json_file in "${GPU1_JOBS[@]}"; do
        job_name=$(basename "$json_file" .json)
        run_one "$job_name" 1
        if [ $? -eq 0 ]; then ((completed++)); else ((failed++)); fi
    done
    echo "[GPU1] FINISHED: ${completed} done, ${failed} failed $(date)"
) &
GPU1_PID=$!

echo "GPU 0 PID: ${GPU0_PID}"
echo "GPU 1 PID: ${GPU1_PID}"

wait $GPU0_PID
wait $GPU1_PID

echo ""
echo "========================================="
echo "BATCH COMPLETE $(date)"
echo "========================================="
total_done=$(find "$OUTPUT_DIR" -name "*_summary_confidences.json" -maxdepth 2 -type f 2>/dev/null | wc -l)
echo "total completed predictions: ${total_done}"
