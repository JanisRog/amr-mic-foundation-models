#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_shared_gwas.sh"
setup_stage_logging "02r_bakta_reference"

require_file "$REFERENCE_FASTA"
require_file "$BAKTA_SIF"
require_dir "$BAKTA_DB"
CTR_BIN="$(resolve_container_bin)"

BAKTA_REFERENCE_NAME="${BAKTA_REFERENCE_NAME:-$REFERENCE_NAME}"
BAKTA_REFERENCE_OUTROOT="${BAKTA_REFERENCE_OUTROOT:-$BAKTA_OUTROOT}"
BAKTA_AMRFINDERPLUS_DB="${BAKTA_AMRFINDERPLUS_DB:-}"
BAKTA_REFERENCE_FORCE="${BAKTA_REFERENCE_FORCE:-0}"

REF_OUT="$BAKTA_REFERENCE_OUTROOT/$BAKTA_REFERENCE_NAME"
GFF3_OUT="$REF_OUT/${BAKTA_REFERENCE_NAME}.gff3"
TMP_ROOT="$BAKTA_REFERENCE_OUTROOT/tmp"
REF_TMP="$TMP_ROOT/${BAKTA_REFERENCE_NAME}_${SLURM_JOB_ID:-$$}"
DB_LINK_ROOT="$TMP_ROOT/db_link_${BAKTA_REFERENCE_NAME}_${SLURM_JOB_ID:-$$}"

mkdir -p "$REF_OUT" "$TMP_ROOT" "$REF_TMP"

if [[ -s "$GFF3_OUT" && "$BAKTA_REFERENCE_FORCE" != "1" ]]; then
  echo "[SKIP] reference already annotated: $GFF3_OUT"
  refresh_symlink "$REF_OUT" "$BAKTA_REFERENCE_OUTROOT/reference_latest"
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

unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"
export TMPDIR="$REF_TMP"

CMD=(
  "$CTR_BIN" exec "$BAKTA_SIF" bakta
  --db "$RUN_BAKTA_DB"
  --output "$REF_OUT"
  --prefix "$BAKTA_REFERENCE_NAME"
  --threads "${SLURM_CPUS_PER_TASK:-1}"
  --tmp-dir "$REF_TMP"
  "$REFERENCE_FASTA"
)

if [[ "$BAKTA_REFERENCE_FORCE" == "1" || -d "$REF_OUT" ]]; then
  CMD+=(--force)
fi
[[ "$BAKTA_SKIP_TRNA" == "1" ]] && CMD+=(--skip-trna)
[[ "$BAKTA_DEBUG" == "1" ]] && CMD+=(--debug --verbose)

echo "[INFO] reference_fasta=$REFERENCE_FASTA"
echo "[INFO] reference_name=$BAKTA_REFERENCE_NAME"
echo "[INFO] outdir=$REF_OUT"
echo "[INFO] gff3_out=$GFF3_OUT"
echo "[INFO] cmd=${CMD[*]}"

"${CMD[@]}"

require_file "$GFF3_OUT"
rm -rf "$REF_TMP"
refresh_symlink "$REF_OUT" "$BAKTA_REFERENCE_OUTROOT/reference_latest"

echo "[OK] bakta_reference_gff3=$GFF3_OUT"
echo "[OK] bakta_reference_latest=$BAKTA_REFERENCE_OUTROOT/reference_latest"
