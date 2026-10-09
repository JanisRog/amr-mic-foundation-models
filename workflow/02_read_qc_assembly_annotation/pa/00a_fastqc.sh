#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_shared_gwas.sh"
setup_stage_logging "00a_fastqc"

require_dir "$RAW_READS_DIR"
[[ -e "$FASTQC_IMG" ]] || { echo "[ERROR] missing FastQC container/image: $FASTQC_IMG" >&2; exit 1; }
CTR_BIN="$(resolve_container_bin)"

STAMP="$(date +%Y%m%d_%H%M%S)"
mkdir -p "$QC_FASTQC_OUTROOT"
RUN_OUT="$QC_FASTQC_OUTROOT/run_${STAMP}"
mkdir -p "$RUN_OUT"

FASTQC_BATCH_SIZE="${FASTQC_BATCH_SIZE:-200}"
mapfile -d '' FASTQS < <(find "$RAW_READS_DIR" -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' -o -name '*.fastq' -o -name '*.fq' \) -print0 | sort -z)

echo "[INFO] raw_reads_dir=$RAW_READS_DIR"
echo "[INFO] run_out=$RUN_OUT"
echo "[INFO] fastq_count=${#FASTQS[@]}"

if [[ "${#FASTQS[@]}" -eq 0 ]]; then
  touch "$RUN_OUT/SKIPPED_NO_FASTQ.txt"
  refresh_symlink "$RUN_OUT" "$QC_FASTQC_OUTROOT/latest"
  echo "[WARN] no FASTQ files found; skipping FastQC"
  exit 0
fi

unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"

batch=()
for fq in "${FASTQS[@]}"; do
  batch+=("$fq")
  if [[ "${#batch[@]}" -ge "$FASTQC_BATCH_SIZE" ]]; then
    "$CTR_BIN" exec --bind "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$FASTQC_IMG" fastqc "${batch[@]}" --outdir "$RUN_OUT" -t "${SLURM_CPUS_PER_TASK:-1}"
    batch=()
  fi
done
if [[ "${#batch[@]}" -gt 0 ]]; then
  "$CTR_BIN" exec --bind "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$FASTQC_IMG" fastqc "${batch[@]}" --outdir "$RUN_OUT" -t "${SLURM_CPUS_PER_TASK:-1}"
fi

FIRST_FASTQC_OUTPUT="$(find "$RUN_OUT" -maxdepth 1 -type f \( -name '*_fastqc.html' -o -name '*_fastqc.zip' \) -print -quit)"
[[ -n "$FIRST_FASTQC_OUTPUT" ]] || { echo "[ERROR] FastQC finished but produced no *_fastqc.html or *_fastqc.zip outputs" >&2; exit 1; }
refresh_symlink "$RUN_OUT" "$QC_FASTQC_OUTROOT/latest"
echo "[OK] fastqc_latest=$QC_FASTQC_OUTROOT/latest"
