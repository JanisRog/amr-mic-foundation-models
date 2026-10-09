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
setup_stage_logging "00c_quast_worker"

[[ -e "$QUAST_IMG" ]] || { echo "[ERROR] missing QUAST container/image: $QUAST_IMG" >&2; exit 1; }
require_file "${QUAST_SAMPLE_TSV:-}"
require_dir "${QUAST_RUN_OUT:-}"
CTR_BIN="$(resolve_container_bin)"

TASK="${SLURM_ARRAY_TASK_ID}"
N=$(($(wc -l < "$QUAST_SAMPLE_TSV") - 1))
if (( TASK > N )); then
  echo "[INFO] task $TASK > N $N; nothing to do"
  exit 0
fi

LINE=$(sed -n "$((TASK+1))p" "$QUAST_SAMPLE_TSV")
SAMPLE=$(echo "$LINE" | cut -f1)
ASM=$(echo "$LINE" | cut -f2-)
OUT="$QUAST_RUN_OUT/quast_${SAMPLE}"

require_file "$ASM"
mkdir -p "$OUT"

if [[ -s "$OUT/report.tsv" ]]; then
  echo "[SKIP] sample=$SAMPLE already has QUAST report"
  exit 0
fi

find "$OUT" -mindepth 1 -maxdepth 1 -exec rm -rf {} +

echo "[INFO] sample=$SAMPLE"
echo "[INFO] asm=$ASM"
echo "[INFO] out=$OUT"

unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"
CMD=("$CTR_BIN" exec --bind "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$QUAST_IMG" quast --threads "${SLURM_CPUS_PER_TASK:-1}" --min-contig "$QUAST_MIN_CONTIG" -o "$OUT" "$ASM")
set +e
"${CMD[@]}" > "$OUT/quast.stdout.log" 2> "$OUT/quast.stderr.log"
rc=$?
set -e

if (( rc != 0 )); then
  echo "[WARN] sample=$SAMPLE failed rc=$rc"
  exit 0
fi

if [[ -s "$OUT/report.tsv" ]]; then
  echo "[OK] sample=$SAMPLE complete"
else
  echo "[WARN] sample=$SAMPLE finished but report.tsv missing"
fi
