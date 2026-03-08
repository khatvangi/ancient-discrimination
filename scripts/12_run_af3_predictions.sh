#!/bin/bash
#SBATCH --job-name=af3_ancient
#SBATCH --output=structures/af3_logs/af3_%A_%a.out
#SBATCH --error=structures/af3_logs/af3_%A_%a.err
#SBATCH --array=1-296%2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#
# runs AF3 predictions as a SLURM array job.
# %2 = max 2 concurrent jobs (one per GPU).
# each array task picks one JSON input file from the job list.
#
# usage:
#   cd /storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination
#   mkdir -p structures/af3_logs structures/af3_outputs
#   sbatch scripts/12_run_af3_predictions.sh

set -euo pipefail

PROJECT_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
JOBLIST="${PROJECT_DIR}/structures/af3_jobs.txt"
OUTPUT_DIR="${PROJECT_DIR}/structures/af3_outputs"

# get the JSON path for this array task
JSON_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")

if [ -z "$JSON_PATH" ]; then
    echo "ERROR: no job found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

JOB_NAME=$(basename "$JSON_PATH" .json)

# skip if already completed
if [ -f "${OUTPUT_DIR}/${JOB_NAME}/${JOB_NAME}_summary_confidences.json" ]; then
    echo "SKIP: ${JOB_NAME} already completed"
    exit 0
fi

echo "========================================="
echo "AF3 prediction: ${JOB_NAME}"
echo "Array task: ${SLURM_ARRAY_TASK_ID} / 296"
echo "GPU: ${CUDA_VISIBLE_DEVICES:-auto}"
echo "JSON: ${JSON_PATH}"
echo "Output: ${OUTPUT_DIR}"
echo "========================================="

# activate conda environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate af3_mmseqs2

# set AF3 environment variables
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

# run AF3
python /storage/kiran-stuff/AF3_mmseqs2/run_alphafold.py \
    --json_path="$JSON_PATH" \
    --model_dir=/storage/kiran-stuff/alphafold3_models \
    --db_dir=/storage/kiran-stuff/alphafold_databases \
    --output_dir="$OUTPUT_DIR" \
    --buckets='256,512,768,1024,1536,2048' \
    --flash_attention_implementation=xla

echo "========================================="
echo "Completed: ${JOB_NAME}"
echo "========================================="
