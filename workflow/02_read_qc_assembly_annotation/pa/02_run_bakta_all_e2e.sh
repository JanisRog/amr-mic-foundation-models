#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_shared_gwas.sh"
setup_stage_logging "02_run_bakta_all_e2e"

require_dir "$ASSEMBLIES_DIR"
require_file "$BAKTA_SIF"
require_dir "$BAKTA_DB"
WORKER="$GWAS_PIPELINE_ROOT/wrappers/02a_bakta_worker.sh"
require_file "$WORKER"
mkdir -p "$BAKTA_OUTROOT"

BAKTA_SAMPLE_FILE="${BAKTA_SAMPLE_FILE_OVERRIDE:-}"
if [[ -n "$BAKTA_SAMPLE_FILE" ]]; then
  require_file "$BAKTA_SAMPLE_FILE"
else
  BAKTA_SAMPLE_FILE="$BAKTA_OUTROOT/bakta_all_sample_ids.txt"
  find -L "$ASSEMBLIES_DIR" -maxdepth 1 -type f \( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' \) -printf '%f\n' | sed -E 's/\.(fasta|fa|fna)$//' | sort -u > "$BAKTA_SAMPLE_FILE"
fi
N=$(wc -l < "$BAKTA_SAMPLE_FILE")
(( N >= 1 )) || { echo "[ERROR] no samples listed for bakta" >&2; exit 1; }

echo "[INFO] assemblies_dir=$ASSEMBLIES_DIR"
echo "[INFO] sample_file=$BAKTA_SAMPLE_FILE"
echo "[INFO] n_samples=$N"
echo "[INFO] worker_cpus=$BAKTA_WORKER_CPUS"
echo "[INFO] worker_mem=$BAKTA_WORKER_MEM"
echo "[INFO] worker_time=$BAKTA_WORKER_TIME"
ARRAY_JOBID=$(sbatch --parsable --array=1-"$N"%"$BAKTA_MAX_PARALLEL" \
  --cpus-per-task="$BAKTA_WORKER_CPUS" \
  --mem="$BAKTA_WORKER_MEM" \
  --time="$BAKTA_WORKER_TIME" \
  --export=ALL,CONFIG_ENV="$CONFIG_ENV",BAKTA_SAMPLE_FILE="$BAKTA_SAMPLE_FILE" \
  "$WORKER")
echo "[INFO] submitted array=$ARRAY_JOBID"
while squeue -h -j "$ARRAY_JOBID" | grep -q .; do
  sleep 30
done

echo "[INFO] bakta worker array completed"
DONE_FILE="$BAKTA_OUTROOT/bakta_done_from_outputs_all.txt"
for d in "$BAKTA_OUTROOT"/*; do
  [[ -d "$d" ]] || continue
  s="$(basename "$d")"
  [[ "$s" == "status" || "$s" == "tmp" ]] && continue
  [[ -s "$d/$s.gff3" ]] && echo "$s"
done | sort > "$DONE_FILE"
echo "[OK] done_file=$DONE_FILE"
echo "[OK] done_count=$(wc -l < "$DONE_FILE")"
