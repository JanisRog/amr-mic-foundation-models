#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_shared_gwas.sh"
setup_stage_logging "00b_multiqc"

[[ -e "$MULTIQC_IMG" ]] || { echo "[ERROR] missing MultiQC container/image: $MULTIQC_IMG" >&2; exit 1; }
CTR_BIN="$(resolve_container_bin)"

FASTQC_DIR="$QC_FASTQC_OUTROOT/latest"
require_dir "$FASTQC_DIR"
FASTQC_SCAN_DIR="$(cd "$FASTQC_DIR" && pwd -P)"

STAMP="$(date +%Y%m%d_%H%M%S)"
mkdir -p "$QC_MULTIQC_READS_OUTROOT"
RUN_OUT="$QC_MULTIQC_READS_OUTROOT/run_${STAMP}"
mkdir -p "$RUN_OUT"

echo "[INFO] fastqc_dir=$FASTQC_DIR"
echo "[INFO] fastqc_scan_dir=$FASTQC_SCAN_DIR"
echo "[INFO] run_out=$RUN_OUT"

FIRST_FASTQC_OUTPUT="$(find "$FASTQC_SCAN_DIR" -maxdepth 1 -type f \( -name '*_fastqc.html' -o -name '*_fastqc.zip' \) -print -quit)"
if [[ -z "$FIRST_FASTQC_OUTPUT" ]]; then
  touch "$RUN_OUT/SKIPPED_NO_FASTQC_OUTPUTS.txt"
  refresh_symlink "$RUN_OUT" "$QC_MULTIQC_READS_OUTROOT/latest"
  echo "[WARN] no FastQC outputs found; skipping MultiQC"
  exit 0
fi

unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"
CMD=("$CTR_BIN" exec --bind "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$MULTIQC_IMG" multiqc "$FASTQC_SCAN_DIR" -o "$RUN_OUT")
"${CMD[@]}"

require_file "$RUN_OUT/multiqc_report.html"
refresh_symlink "$RUN_OUT" "$QC_MULTIQC_READS_OUTROOT/latest"
echo "[OK] multiqc_reads_latest=$QC_MULTIQC_READS_OUTROOT/latest"
