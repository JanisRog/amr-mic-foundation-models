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
setup_stage_logging "02b_amrfinderplus_worker"

AMRFINDERPLUS_SIF="${AMRFINDERPLUS_SIF:?set AMRFINDERPLUS_SIF in config/paths.env}"
AMRFINDERPLUS_DB="${AMRFINDERPLUS_DB:?set AMRFINDERPLUS_DB in config/paths.env}"
AMRFINDERPLUS_OUTROOT="${AMRFINDERPLUS_OUTROOT:-$RESULTS_ROOT/amrfinderplus/assemblies}"
AMRFINDERPLUS_RUN_OUT="${AMRFINDERPLUS_RUN_OUT:-}"
AMRFINDERPLUS_PLUS="${AMRFINDERPLUS_PLUS:-1}"
AMRFINDERPLUS_PRINT_NODE="${AMRFINDERPLUS_PRINT_NODE:-1}"
AMRFINDERPLUS_ORGANISM_RESOLVED="${AMRFINDERPLUS_ORGANISM_RESOLVED:-${AMRFINDERPLUS_ORGANISM:-}}"
if [[ -d "$AMRFINDERPLUS_DB/latest" ]]; then
  AMRFINDERPLUS_DB="$AMRFINDERPLUS_DB/latest"
fi

require_dir "$ASSEMBLIES_DIR"
require_file "$AMRFINDERPLUS_SIF"
require_dir "$AMRFINDERPLUS_DB"
require_file "$AMRFINDERPLUS_SAMPLE_FILE"
[[ -n "$AMRFINDERPLUS_RUN_OUT" ]] || { echo "[ERROR] AMRFINDERPLUS_RUN_OUT not set" >&2; exit 1; }
CTR_BIN="$(resolve_container_bin)"

SAMPLE="${SAMPLE_ID:-}"
if [[ -z "$SAMPLE" ]]; then
  TASK="${SLURM_ARRAY_TASK_ID:-}"
  [[ -n "$TASK" ]] || { echo "[ERROR] SLURM_ARRAY_TASK_ID missing" >&2; exit 1; }
  SAMPLE="$(sed -n "${TASK}p" "$AMRFINDERPLUS_SAMPLE_FILE" | tr -d '\r')"
fi
[[ -n "$SAMPLE" ]] || { echo "[ERROR] empty sample" >&2; exit 1; }

resolve_assembly() {
  local sample="$1"
  for ext in fasta fa fna; do
    if [[ -f "$ASSEMBLIES_DIR/${sample}.${ext}" ]]; then
      printf '%s\n' "$ASSEMBLIES_DIR/${sample}.${ext}"
      return 0
    fi
  done
  echo "[ERROR] missing assembly for $sample" >&2
  exit 1
}

ASSEMBLY="$(resolve_assembly "$SAMPLE")"
STATUS_DIR="$AMRFINDERPLUS_RUN_OUT/status"
TMP_ROOT="$AMRFINDERPLUS_RUN_OUT/tmp"
SAMPLE_OUT="$AMRFINDERPLUS_RUN_OUT/$SAMPLE"
TSV_OUT="$SAMPLE_OUT/${SAMPLE}.amrfinder.tsv"
SAMPLE_TMP="$TMP_ROOT/${SAMPLE}_${SLURM_JOB_ID:-$$}"
mkdir -p "$STATUS_DIR" "$TMP_ROOT" "$SAMPLE_OUT" "$SAMPLE_TMP"

if [[ -s "$TSV_OUT" ]]; then
  echo "[SKIP] sample=$SAMPLE already processed"
  exit 0
fi

unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"
export TMPDIR="$SAMPLE_TMP"

CMD=(
  "$CTR_BIN" exec --cleanenv -B "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$AMRFINDERPLUS_SIF" amrfinder
  -n "$ASSEMBLY"
  -o "$TSV_OUT"
  -d "$AMRFINDERPLUS_DB"
  --threads "${SLURM_CPUS_PER_TASK:-1}"
  --name "$SAMPLE"
)
[[ "$AMRFINDERPLUS_PLUS" == "1" ]] && CMD+=(--plus)
[[ "$AMRFINDERPLUS_PRINT_NODE" == "1" ]] && CMD+=(--print_node)
[[ -n "$AMRFINDERPLUS_ORGANISM_RESOLVED" ]] && CMD+=(--organism "$AMRFINDERPLUS_ORGANISM_RESOLVED")

echo "[INFO] sample=$SAMPLE"
echo "[INFO] assembly=$ASSEMBLY"
echo "[INFO] outdir=$SAMPLE_OUT"
echo "[INFO] organism=${AMRFINDERPLUS_ORGANISM_RESOLVED:-<none>}"
echo "[INFO] cmd=${CMD[*]}"
"${CMD[@]}" >"$SAMPLE_OUT/${SAMPLE}.stdout.log" 2>"$SAMPLE_OUT/${SAMPLE}.stderr.log"
require_file "$TSV_OUT"
rm -rf "$SAMPLE_TMP"
echo "[OK] sample=$SAMPLE"
