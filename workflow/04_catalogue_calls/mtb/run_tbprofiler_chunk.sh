#!/usr/bin/env bash

set -uo pipefail

: "${SAMPLE_LIST:?set SAMPLE_LIST}"
: "${SAMPLE_WORKER:?set SAMPLE_WORKER}"
: "${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}"

SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-10}"
n_samples=$(wc -l < "$SAMPLE_LIST")
start=$(( (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
end=$(( start + SAMPLES_PER_TASK - 1 ))
(( end > n_samples )) && end=$n_samples

echo "[INFO] task=$SLURM_ARRAY_TASK_ID samples=$start-$end total=$n_samples"
failed=0
while IFS= read -r sample; do
  [[ -n "$sample" ]] || continue
  if ! SAMPLE_ID="$sample" "$SAMPLE_WORKER"; then
    echo "[WARN] failed sample=$sample" >&2
    failed=$((failed + 1))
  fi
done < <(sed -n "${start},${end}p" "$SAMPLE_LIST")

echo "[OK] processed=$((end - start + 1)) failed=$failed"
exit 0
