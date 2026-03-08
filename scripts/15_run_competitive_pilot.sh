#!/bin/bash
# run competitive AF3 pilot: 9 small proteins on 2 GPUs
# uses standard AF3 (Docker + HMMER), NOT the mmseqs2 fork
#
# proteins are tiny (70-199 aa + 45 nt = 115-244 tokens)
# each job ~15-30 min on Titan RTX 24GB
# total: ~1-2 hours for all 9 jobs
#
# GPU 0: 5 jobs (CSD family + Alba_Mjann + PCBP1)
# GPU 1: 4 jobs (Alba Sso10b + Rpp25 + KH NusA pair)

set -uo pipefail

AF3_DIR="/storage/kiran-stuff/alphafold3"
INPUT_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/af3_competitive_local"
OUTPUT_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/af3_competitive_output"
LOG_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/af3_logs/competitive"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# --- helper: run one AF3 job via Docker ---
run_one() {
    local json_file="$1"
    local gpu="$2"
    local job_name=$(basename "$json_file" .json)
    local log_file="${LOG_DIR}/${job_name}.log"

    echo "[GPU${gpu}] START ${job_name} $(date)"

    ## fix: --gpus needs single+double quote wrapping for Docker
    ## fix: inside container GPU is always device 0 regardless of physical device
    docker run --rm \
        --gpus "'\"device=${gpu}\"'" \
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

echo "competitive AF3 pilot: 9 jobs on 2 GPUs"
echo "started: $(date)"
echo "input:   ${INPUT_DIR}"
echo "output:  ${OUTPUT_DIR}"
echo "logs:    ${LOG_DIR}"
echo ""

# --- GPU 0: 5 jobs ---
(
    run_one "competitive_CspA_Ecoli.json" 0
    run_one "competitive_CspB_Hfx.json" 0
    run_one "competitive_YBX1_Human.json" 0
    run_one "competitive_Alba_Mjann.json" 0
    run_one "competitive_PCBP1_KH1_Human.json" 0
    echo "[GPU0] ALL DONE $(date)"
) &
GPU0_PID=$!

# --- GPU 1: 4 jobs ---
(
    run_one "competitive_Sso10b_Ssolf.json" 1
    run_one "competitive_Rpp25_Human.json" 1
    run_one "competitive_NusA_KH_Ecoli.json" 1
    run_one "competitive_NusA_Afulg.json" 1
    echo "[GPU1] ALL DONE $(date)"
) &
GPU1_PID=$!

echo "GPU 0 PID: ${GPU0_PID} (CspA, CspB, YBX1, Alba_Mj, PCBP1)"
echo "GPU 1 PID: ${GPU1_PID} (Sso10b, Rpp25, NusA_Ec, NusA_Af)"

# wait for both
wait $GPU0_PID
wait $GPU1_PID

echo ""
echo "========================================="
echo "PILOT COMPLETE $(date)"
echo "========================================="

# count results
completed=$(find "$OUTPUT_DIR" -name "summary_confidences*.json" -type f 2>/dev/null | wc -l)
echo "completed predictions: ${completed}"
echo "check results: ls ${OUTPUT_DIR}/"
