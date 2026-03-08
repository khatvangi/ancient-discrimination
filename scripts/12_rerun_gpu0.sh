#!/bin/bash
# rerun the 3 GPU 0 pilot jobs that failed due to mmseqs race condition.
# databases now exist, so createdb will be skipped.

set -uo pipefail  # removed -e so all jobs run even if one fails

PROJECT_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
OUTPUT_DIR="${PROJECT_DIR}/structures/af3_outputs"
LOG_DIR="${PROJECT_DIR}/structures/af3_logs"
AF3_PYTHON="/home/kiran/miniforge3/envs/af3_mmseqs2/bin/python"
AF3_SCRIPT="/storage/kiran-stuff/AF3_mmseqs2/run_alphafold.py"

export PATH="/home/kiran/miniforge3/bin:$PATH"
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

JOBS=(
    "PF01479_bsu_P37557_RNA"
    "PF01479_bsu_P37557_DNA"
    "PF00633_eco_P0AB83_RNA"
)

echo "rerunning ${#JOBS[@]} GPU0 jobs (databases now exist)"
echo "started: $(date)"

for job_name in "${JOBS[@]}"; do
    json_path="${PROJECT_DIR}/structures/af3_inputs/${job_name}.json"

    # skip if already done
    if [ -f "${OUTPUT_DIR}/${job_name}/${job_name}_summary_confidences.json" ]; then
        echo "SKIP ${job_name} (already done)"
        continue
    fi

    echo "START ${job_name} $(date)"
    CUDA_VISIBLE_DEVICES=0 $AF3_PYTHON $AF3_SCRIPT \
        --json_path="$json_path" \
        --model_dir=/storage/kiran-stuff/alphafold3_models \
        --db_dir=/storage/kiran-stuff/alphafold_databases \
        --output_dir="$OUTPUT_DIR" \
        --buckets='256,512,768,1024,1536,2048' \
        --flash_attention_implementation=xla \
        > "${LOG_DIR}/${job_name}.log" 2>&1

    status=$?
    if [ $status -eq 0 ]; then
        echo "DONE  ${job_name} $(date)"
    else
        echo "FAIL  ${job_name} (exit $status) $(date)"
    fi
done

echo ""
echo "GPU0 RERUN COMPLETE $(date)"
completed=$(find "$OUTPUT_DIR" -name "*_summary_confidences.json" 2>/dev/null | wc -l)
echo "total completed: $completed / 6"
