#!/usr/bin/env bash
set -euo pipefail

CONFIG_ENV="${CONFIG_ENV:-}"
if [[ -n "$CONFIG_ENV" && -f "$CONFIG_ENV" ]]; then
  GWAS_PIPELINE_ROOT="$(cd "$(dirname "$CONFIG_ENV")/.." && pwd)"
  SHARED_GWAS="$GWAS_PIPELINE_ROOT/wrappers/_shared_gwas.sh"
else
  SHARED_GWAS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_shared_gwas.sh"
fi
source "$SHARED_GWAS"
setup_stage_logging "01a_snippy_hybrid_worker"

require_file "$SNIPPY_SAMPLE_TSV"
require_file "$REFERENCE_FASTA"
require_file "$SNIPPY_SIF"
CTR_BIN="$(resolve_container_bin)"

TASK="${SLURM_ARRAY_TASK_ID:-}"
[[ -n "$TASK" ]] || { echo "[ERROR] SLURM_ARRAY_TASK_ID missing" >&2; exit 1; }
LINE="$(awk -F $'\t' -v n="$((TASK+1))" 'NR==n {printf "%s\037%s\037%s\037%s", $1, $2, $3, $4; exit}' "$SNIPPY_SAMPLE_TSV" || true)"
[[ -n "$LINE" ]] || { echo "[INFO] no line for task=$TASK"; exit 0; }
IFS=$'\037' read -r SAMPLE ASM R1 R2 <<< "$LINE"
SAMPLE="${SAMPLE//$'\r'/}"
ASM="${ASM//$'\r'/}"
R1="${R1//$'\r'/}"
R2="${R2//$'\r'/}"
[[ -n "$SAMPLE" ]] || { echo "[ERROR] empty sample at task=$TASK" >&2; exit 1; }

OUT="$SNIPPY_RUN_ROOT/runs/$SAMPLE"
STATUS_DIR="$SNIPPY_RUN_ROOT/status"
mkdir -p "$OUT" "$STATUS_DIR"
INPUT_SELECTED=""
FINAL_STATUS="unknown"
FINAL_RC=0
SNIPPY_FORCE_REDO="${SNIPPY_FORCE_REDO:-0}"

write_status_marker() {
  {
    printf "task\t%s\n" "$TASK"
    printf "sample\t%s\n" "$SAMPLE"
    printf "status\t%s\n" "$FINAL_STATUS"
    printf "input_selected\t%s\n" "${INPUT_SELECTED:-}"
    printf "rc\t%s\n" "${FINAL_RC:-0}"
  } > "$STATUS_DIR/${TASK}.done"
}

trap write_status_marker EXIT

if [[ "$SNIPPY_FORCE_REDO" != "1" && -s "$OUT/snps.vcf" && -s "$OUT/snps.tab" ]]; then
  echo "[SKIP] sample=$SAMPLE already complete"
  FINAL_STATUS="complete_existing"
  exit 0
fi
find "$OUT" -mindepth 1 -maxdepth 1 -exec rm -rf {} +

unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"

echo "[INFO] sample=$SAMPLE"
echo "[INFO] snippy_input_mode=$SNIPPY_INPUT_MODE"
echo "[INFO] snippy_force_redo=$SNIPPY_FORCE_REDO"
echo "[INFO] asm=$ASM"
echo "[INFO] r1=$R1"
echo "[INFO] r2=$R2"

CMD=("$CTR_BIN" exec --cleanenv -B "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$SNIPPY_SIF" snippy --force --cpus "${SLURM_CPUS_PER_TASK:-1}" --outdir "$OUT" --ref "$REFERENCE_FASTA" --minqual "$SNIPPY_MINQUAL" --mincov "$SNIPPY_MINCOV" --minfrac "$SNIPPY_MINFRAC")
INPUT_SELECTED=""
case "$SNIPPY_INPUT_MODE" in
  assemblies_only)
    if [[ -n "$ASM" && -f "$ASM" ]]; then
      CMD+=(--ctgs "$ASM")
      INPUT_SELECTED="assembly"
    else
      echo "[WARN] sample=$SAMPLE missing assembly for assemblies_only mode" >&2
      FINAL_STATUS="missing_input"
      exit 0
    fi
    ;;
  reads_only)
    if [[ -n "$R1" && -n "$R2" && -f "$R1" && -f "$R2" ]]; then
      CMD+=(--R1 "$R1" --R2 "$R2")
      INPUT_SELECTED="reads"
    else
      echo "[WARN] sample=$SAMPLE missing read pair for reads_only mode" >&2
      FINAL_STATUS="missing_input"
      exit 0
    fi
    ;;
  hybrid)
    if [[ -n "$R1" && -n "$R2" && -f "$R1" && -f "$R2" ]]; then
      CMD+=(--R1 "$R1" --R2 "$R2")
      INPUT_SELECTED="reads"
    elif [[ -n "$ASM" && -f "$ASM" ]]; then
      CMD+=(--ctgs "$ASM")
      INPUT_SELECTED="assembly_fallback"
    else
      echo "[WARN] sample=$SAMPLE missing both read pair and assembly" >&2
      FINAL_STATUS="missing_input"
      exit 0
    fi
    ;;
  *)
    echo "[ERROR] unsupported SNIPPY_INPUT_MODE=$SNIPPY_INPUT_MODE" >&2
    FINAL_STATUS="config_error"
    FINAL_RC=1
    exit 1
    ;;
esac
echo "[INFO] input_selected=$INPUT_SELECTED"

set +e
"${CMD[@]}" > "$OUT/snippy.stdout.log" 2> "$OUT/snippy.stderr.log"
rc=$?
set -e
FINAL_RC="$rc"
if (( rc != 0 )); then
  echo "[WARN] sample=$SAMPLE failed rc=$rc"
  FINAL_STATUS="failed"
  exit 0
fi
if [[ -s "$OUT/snps.vcf" && -s "$OUT/snps.tab" ]]; then
  echo "[OK] sample=$SAMPLE complete"
  FINAL_STATUS="complete"
else
  echo "[WARN] sample=$SAMPLE missing expected outputs"
  FINAL_STATUS="partial"
fi
