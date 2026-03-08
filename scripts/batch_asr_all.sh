#!/bin/bash
# batch ASR pipeline for all complete families in nasbench_full_luca.tsv
# runs 4 families in parallel with 4 threads each (16 total)
# usage: bash scripts/batch_asr_all.sh 2>&1 | tee results/asr/batch_log.txt

cd /storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination

# families: pfam_id rna_pdb dna_pdb
FAMILIES=(
  "PF00013 1EC6 1ZTG"
  "PF00136 4Q5V 3QEX"
  "PF00270 3PEY 6CRM"
  "PF00398 8CSQ 9MN5"
  "PF00575 7DID 9HVQ"
  "PF00753 5A0T 8DQ1"
  "PF01131 9GDA 1MW8"
  "PF01336 1C0A 3F2B"
  "PF01479 4LGT 9GUW"
  "PF01588 7K98 7D8T"
  "PF02171 4Z4D 6T5T"
  "PF03129 4YYE 8F69"
  "PF03372 9HDR 5HT2"
  "PF04851 9LOV 2D7D"
  "PF06733 7ML4 4A15"
  "PF07521 9BCU 9BCU"
  "PF07650 8CF1 9GUW"
  "PF10996 9BCU 9BCU"
)

MAX_PARALLEL=4
THREADS=4
RUNNING=0
PIDS=()

echo "=============================================="
echo "BATCH ASR PIPELINE: ${#FAMILIES[@]} families"
echo "parallel: $MAX_PARALLEL, threads per job: $THREADS"
echo "started: $(date)"
echo "=============================================="

for entry in "${FAMILIES[@]}"; do
  read -r pfam rna_pdb dna_pdb <<< "$entry"

  # skip if already completed (convergence summary exists)
  if [ -f "results/asr/${pfam,,}/${pfam}_convergence_summary.tsv" ]; then
    echo "[SKIP] $pfam — already completed"
    continue
  fi

  echo "[START] $pfam (RNA=$rna_pdb DNA=$dna_pdb)"

  # run in background, log per-family output
  mkdir -p results/asr/${pfam,,}
  python3 scripts/asr_run_family.py "$pfam" "$rna_pdb" "$dna_pdb" --threads $THREADS --max-seqs 200 \
    > "results/asr/${pfam,,}/${pfam}_run.log" 2>&1 &
  PIDS+=($!)
  RUNNING=$((RUNNING + 1))

  # wait if we've hit the parallel limit
  if [ $RUNNING -ge $MAX_PARALLEL ]; then
    # wait for any one to finish
    wait -n
    RUNNING=$((RUNNING - 1))
  fi
done

# wait for remaining
echo "[WAITING] for remaining jobs..."
wait

echo "=============================================="
echo "BATCH COMPLETE: $(date)"
echo "=============================================="

# summary: check which families completed
echo ""
echo "RESULTS:"
for entry in "${FAMILIES[@]}"; do
  read -r pfam _ _ <<< "$entry"
  summary="results/asr/${pfam,,}/${pfam}_convergence_summary.tsv"
  if [ -f "$summary" ]; then
    n_cols=$(tail -n +2 "$summary" | wc -l)
    echo "  [OK] $pfam — $n_cols critical columns"
  else
    log="results/asr/${pfam,,}/${pfam}_run.log"
    if [ -f "$log" ]; then
      last_line=$(tail -1 "$log")
      echo "  [FAIL] $pfam — $last_line"
    else
      echo "  [FAIL] $pfam — no output"
    fi
  fi
done
