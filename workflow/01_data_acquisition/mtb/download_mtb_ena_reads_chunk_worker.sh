#!/usr/bin/env bash

set -uo pipefail

: "${METADATA_TSV:?set METADATA_TSV}"
: "${OUTROOT:?set OUTROOT}"
: "${WORKER:?set WORKER}"
: "${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}"

SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-25}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

n_rows=$(( $(wc -l < "$METADATA_TSV") - 1 ))
start=$(( (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
end=$(( start + SAMPLES_PER_TASK - 1 ))
(( end > n_rows )) && end=$n_rows

echo "[INFO] task=$SLURM_ARRAY_TASK_ID rows=$start-$end n_rows=$n_rows"

failed=0
for ((idx=start; idx<=end; idx++)); do
  echo "[INFO] metadata_row=$idx"
  if ! SLURM_ARRAY_TASK_ID="$idx" "$PYTHON_BIN" "$WORKER"; then
    echo "[WARN] worker failed metadata_row=$idx" >&2
    failed=$((failed + 1))
  fi
done

echo "[OK] processed=$((end - start + 1)) failed=$failed"
exit 0
