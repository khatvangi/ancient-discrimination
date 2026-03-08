#!/bin/bash
# NOTE: Phase 3 script (planned, not yet completed).
# AF3 pilot v2 — feasibility test only, not used in current analysis.
# AF3 pilot v2: fixed binary paths for hmmbuild, hmmsearch, etc.
#
# previous failures:
#   v1 GPU 0: mmseqs createdb race condition
#   v1/v2 GPU 0+1: hmmbuild not found — AF3 uses os.path.exists() not PATH lookup
#
# fix: pass absolute paths via --hmmbuild_binary_path etc.
#      (AF3's check_binary_exists uses os.path.exists, not shutil.which)
#
# databases (uniref90, mgnify, small_bfd, uniprot) already created from v1 runs
#
# usage:
#   cd /storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination
#   nohup bash scripts/12_run_af3_pilot_v2.sh > structures/af3_logs/pilot_v2_master.log 2>&1 &

set -uo pipefail  # no -e: continue if a job fails

PROJECT_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
OUTPUT_DIR="${PROJECT_DIR}/structures/af3_outputs"
LOG_DIR="${PROJECT_DIR}/structures/af3_logs"
AF3_PYTHON="/home/kiran/miniforge3/envs/af3_mmseqs2/bin/python"
AF3_SCRIPT="/storage/kiran-stuff/AF3_mmseqs2/run_alphafold.py"

# absolute paths to all binaries AF3 needs
# (AF3 uses os.path.exists() not PATH lookup — must be absolute)
ENV_BIN="/home/kiran/miniforge3/envs/af3_mmseqs2/bin"
HMMBUILD_PATH="${ENV_BIN}/hmmbuild"
HMMSEARCH_PATH="${ENV_BIN}/hmmsearch"
HMMALIGN_PATH="${ENV_BIN}/hmmalign"
JACKHMMER_PATH="${ENV_BIN}/jackhmmer"
NHMMER_PATH="${ENV_BIN}/nhmmer"
MMSEQS2_PATH="/home/kiran/miniforge3/bin/mmseqs"

# still need PATH for mmseqs subprocesses that AF3 spawns internally
export PATH="${ENV_BIN}:/home/kiran/miniforge3/bin:$PATH"
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# verify binaries before starting
for bin in mmseqs hmmbuild hmmsearch hmmalign; do
    if ! command -v $bin &>/dev/null; then
        echo "ERROR: $bin not found on PATH"
        exit 1
    fi
done
echo "all binaries found on PATH"

run_one() {
    local job_name="$1"
    local gpu="$2"
    local json_path="${PROJECT_DIR}/structures/af3_inputs/${job_name}.json"

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
        --hmmbuild_binary_path="$HMMBUILD_PATH" \
        --hmmsearch_binary_path="$HMMSEARCH_PATH" \
        --hmmalign_binary_path="$HMMALIGN_PATH" \
        --jackhmmer_binary_path="$JACKHMMER_PATH" \
        --nhmmer_binary_path="$NHMMER_PATH" \
        --mmseqs2_binary_path="$MMSEQS2_PATH" \
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

echo "pilot v2: 6 jobs on 2 GPUs (PATH fixed)"
echo "started: $(date)"

# GPU 0: S4-RNA, S4-DNA, HHH-RNA
(
    run_one "PF01479_bsu_P37557_RNA" 0
    run_one "PF01479_bsu_P37557_DNA" 0
    run_one "PF00633_eco_P0AB83_RNA" 0
    echo "[GPU0] ALL DONE $(date)"
) &
PID_GPU0=$!

# GPU 1: HHH-DNA, S1-RNA, S1-DNA
(
    run_one "PF00633_eco_P0AB83_DNA" 1
    run_one "PF00575_mja_Q57840_RNA" 1
    run_one "PF00575_mja_Q57840_DNA" 1
    echo "[GPU1] ALL DONE $(date)"
) &
PID_GPU1=$!

echo "GPU 0 PID: $PID_GPU0 (S4-RNA, S4-DNA, HHH-RNA)"
echo "GPU 1 PID: $PID_GPU1 (HHH-DNA, S1-RNA, S1-DNA)"

wait $PID_GPU0
wait $PID_GPU1

echo ""
echo "========================================="
echo "PILOT V2 COMPLETE $(date)"
echo "========================================="

completed=$(find "$OUTPUT_DIR" -name "*_summary_confidences.json" 2>/dev/null | wc -l)
echo "completed: $completed / 6"

# show confidence scores
echo ""
echo "PILOT RESULTS:"
for conf in $(find "$OUTPUT_DIR" -name "*_summary_confidences.json" 2>/dev/null | sort); do
    job=$(basename $(dirname "$conf"))
    iptm=$($AF3_PYTHON -c "import json; d=json.load(open('$conf')); print(f'{d.get(\"iptm\", 0):.3f}')")
    ptm=$($AF3_PYTHON -c "import json; d=json.load(open('$conf')); print(f'{d.get(\"ptm\", 0):.3f}')")
    echo "  $job: ipTM=$iptm pTM=$ptm"
done
