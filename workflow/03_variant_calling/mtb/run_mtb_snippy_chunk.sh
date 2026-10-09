#!/usr/bin/env bash
set -euo pipefail

start=$(( (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
end=$(( start + SAMPLES_PER_TASK - 1 ))
failed=0
while IFS= read -r sample; do
  [[ -n "$sample" ]] || continue
  if ! SAMPLE_ID="$sample" bash "$SAMPLE_WORKER"; then
    echo "[WARN] failed sample=$sample" >&2
    failed=$((failed + 1))
  fi
done < <(sed -n "${start},${end}p" "$SAMPLE_LIST")
echo "[OK] processed_range=${start}-${end} failed=$failed"
[[ "$failed" -eq 0 ]]
