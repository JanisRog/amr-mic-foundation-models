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
setup_stage_logging "02a_bakta_worker"

require_dir "$ASSEMBLIES_DIR"
require_file "$BAKTA_SIF"
require_dir "$BAKTA_DB"
require_file "$BAKTA_SAMPLE_FILE"
CTR_BIN="$(resolve_container_bin)"
BAKTA_AMRFINDERPLUS_DB="${BAKTA_AMRFINDERPLUS_DB:-}"

SAMPLE="${SAMPLE_ID:-}"
if [[ -z "$SAMPLE" ]]; then
  TASK="${SLURM_ARRAY_TASK_ID:-}"
  [[ -n "$TASK" ]] || { echo "[ERROR] SLURM_ARRAY_TASK_ID missing" >&2; exit 1; }
  SAMPLE="$(sed -n "${TASK}p" "$BAKTA_SAMPLE_FILE" | tr -d '\r')"
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
STATUS_DIR="$BAKTA_OUTROOT/status"
TMP_ROOT="$BAKTA_OUTROOT/tmp"
SAMPLE_OUT="$BAKTA_OUTROOT/$SAMPLE"
GFF3_OUT="$SAMPLE_OUT/${SAMPLE}.gff3"
SAMPLE_TMP="$TMP_ROOT/${SAMPLE}_${SLURM_JOB_ID:-$$}"
DB_LINK_ROOT="$TMP_ROOT/db_link_${SAMPLE}_${SLURM_JOB_ID:-$$}"
mkdir -p "$STATUS_DIR" "$TMP_ROOT" "$SAMPLE_OUT" "$SAMPLE_TMP"

if [[ -s "$GFF3_OUT" ]]; then
  echo "[SKIP] sample=$SAMPLE already annotated"
  exit 0
fi

RUN_BAKTA_DB="$BAKTA_DB"
if [[ -n "$BAKTA_AMRFINDERPLUS_DB" && -d "$BAKTA_AMRFINDERPLUS_DB" ]]; then
  RUN_BAKTA_DB="$DB_LINK_ROOT"
  rm -rf "$DB_LINK_ROOT"
  mkdir -p "$DB_LINK_ROOT"
  for item in "$BAKTA_DB"/*; do
    [[ -e "$item" ]] || continue
    base="$(basename "$item")"
    [[ "$base" == "amrfinderplus-db" ]] && continue
    ln -s "$item" "$DB_LINK_ROOT/$base"
  done
  ln -s "$BAKTA_AMRFINDERPLUS_DB" "$DB_LINK_ROOT/amrfinderplus-db"
fi

RUN_FORCE=0
if [[ -d "$SAMPLE_OUT" && ! -s "$GFF3_OUT" ]]; then
  RUN_FORCE=1
fi

unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"
export TMPDIR="$SAMPLE_TMP"
CMD=("$CTR_BIN" exec "$BAKTA_SIF" bakta --db "$RUN_BAKTA_DB" --output "$SAMPLE_OUT" --prefix "$SAMPLE" --threads "${SLURM_CPUS_PER_TASK:-1}" --tmp-dir "$SAMPLE_TMP" "$ASSEMBLY")
[[ "$RUN_FORCE" == "1" ]] && CMD+=(--force)
[[ "$BAKTA_SKIP_TRNA" == "1" ]] && CMD+=(--skip-trna)
[[ "$BAKTA_DEBUG" == "1" ]] && CMD+=(--debug --verbose)

echo "[INFO] sample=$SAMPLE"
echo "[INFO] assembly=$ASSEMBLY"
echo "[INFO] outdir=$SAMPLE_OUT"
echo "[INFO] cmd=${CMD[*]}"
"${CMD[@]}"
require_file "$GFF3_OUT"
rm -rf "$SAMPLE_TMP"
echo "[OK] sample=$SAMPLE"
