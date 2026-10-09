#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_shared_gwas.sh"
setup_stage_logging "02b_run_amrfinderplus_all_e2e"

AMRFINDERPLUS_SIF="${AMRFINDERPLUS_SIF:?set AMRFINDERPLUS_SIF in config/paths.env}"
AMRFINDERPLUS_DB="${AMRFINDERPLUS_DB:?set AMRFINDERPLUS_DB in config/paths.env}"
AMRFINDERPLUS_OUTROOT="${AMRFINDERPLUS_OUTROOT:-$RESULTS_ROOT/amrfinderplus/assemblies}"
AMRFINDERPLUS_MAX_PARALLEL="${AMRFINDERPLUS_MAX_PARALLEL:-16}"
AMRFINDERPLUS_WORKER_CPUS="${AMRFINDERPLUS_WORKER_CPUS:-8}"
AMRFINDERPLUS_WORKER_MEM="${AMRFINDERPLUS_WORKER_MEM:-16G}"
AMRFINDERPLUS_WORKER_TIME="${AMRFINDERPLUS_WORKER_TIME:-08:00:00}"
AMRFINDERPLUS_SAMPLE_FILE_OVERRIDE="${AMRFINDERPLUS_SAMPLE_FILE_OVERRIDE:-}"
AMRFINDERPLUS_ORGANISM="${AMRFINDERPLUS_ORGANISM:-Streptococcus_pneumoniae}"
AMRFINDERPLUS_PLUS="${AMRFINDERPLUS_PLUS:-1}"
AMRFINDERPLUS_PRINT_NODE="${AMRFINDERPLUS_PRINT_NODE:-1}"
if [[ -d "$AMRFINDERPLUS_DB/latest" ]]; then
  AMRFINDERPLUS_DB="$AMRFINDERPLUS_DB/latest"
fi

require_dir "$ASSEMBLIES_DIR"
require_file "$AMRFINDERPLUS_SIF"
require_dir "$AMRFINDERPLUS_DB"
WORKER="$GWAS_PIPELINE_ROOT/wrappers/02b_amrfinderplus_worker.sh"
require_file "$WORKER"
CTR_BIN="$(resolve_container_bin)"

mkdir -p "$AMRFINDERPLUS_OUTROOT"
RUN_OUT="$AMRFINDERPLUS_OUTROOT/run_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUN_OUT" "$RUN_OUT/meta"

AMRFINDERPLUS_SAMPLE_FILE="$AMRFINDERPLUS_SAMPLE_FILE_OVERRIDE"
if [[ -n "$AMRFINDERPLUS_SAMPLE_FILE" ]]; then
  require_file "$AMRFINDERPLUS_SAMPLE_FILE"
else
  AMRFINDERPLUS_SAMPLE_FILE="$RUN_OUT/meta/amrfinderplus_all_sample_ids.txt"
  find -L "$ASSEMBLIES_DIR" -maxdepth 1 -type f \( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' \) -printf '%f\n' \
    | sed -E 's/\.(fasta|fa|fna)$//' | sort -u > "$AMRFINDERPLUS_SAMPLE_FILE"
fi
N=$(wc -l < "$AMRFINDERPLUS_SAMPLE_FILE")
(( N >= 1 )) || { echo "[ERROR] no samples listed for amrfinderplus" >&2; exit 1; }

ORG_LIST="$RUN_OUT/meta/amrfinderplus_supported_organisms.txt"
ORG_LIST_ERR="$RUN_OUT/meta/amrfinderplus_supported_organisms.stderr.log"
AMRFINDERPLUS_ORGANISM_RESOLVED="$AMRFINDERPLUS_ORGANISM"
if [[ -n "$AMRFINDERPLUS_ORGANISM" ]]; then
  if "$CTR_BIN" exec --cleanenv -B "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$AMRFINDERPLUS_SIF" amrfinder --list_organisms -d "$AMRFINDERPLUS_DB" >"$ORG_LIST" 2>"$ORG_LIST_ERR"; then
    if ! grep -Fxq "$AMRFINDERPLUS_ORGANISM" "$ORG_LIST"; then
      echo "[WARN] organism '$AMRFINDERPLUS_ORGANISM' not found in --list_organisms; proceeding without --organism"
      AMRFINDERPLUS_ORGANISM_RESOLVED=""
    fi
  else
    echo "[WARN] could not query --list_organisms; proceeding with organism='$AMRFINDERPLUS_ORGANISM'"
  fi
fi

echo "[INFO] assemblies_dir=$ASSEMBLIES_DIR"
echo "[INFO] sample_file=$AMRFINDERPLUS_SAMPLE_FILE"
echo "[INFO] n_samples=$N"
echo "[INFO] amrfinderplus_sif=$AMRFINDERPLUS_SIF"
echo "[INFO] amrfinderplus_db=$AMRFINDERPLUS_DB"
echo "[INFO] amrfinderplus_outroot=$AMRFINDERPLUS_OUTROOT"
echo "[INFO] run_out=$RUN_OUT"
echo "[INFO] organism_requested=${AMRFINDERPLUS_ORGANISM:-<none>}"
echo "[INFO] organism_resolved=${AMRFINDERPLUS_ORGANISM_RESOLVED:-<none>}"
echo "[INFO] plus=$AMRFINDERPLUS_PLUS"
echo "[INFO] print_node=$AMRFINDERPLUS_PRINT_NODE"
echo "[INFO] worker_cpus=$AMRFINDERPLUS_WORKER_CPUS"
echo "[INFO] worker_mem=$AMRFINDERPLUS_WORKER_MEM"
echo "[INFO] worker_time=$AMRFINDERPLUS_WORKER_TIME"

ARRAY_JOBID=$(sbatch --parsable --array=1-"$N"%"$AMRFINDERPLUS_MAX_PARALLEL" \
  --cpus-per-task="$AMRFINDERPLUS_WORKER_CPUS" \
  --mem="$AMRFINDERPLUS_WORKER_MEM" \
  --time="$AMRFINDERPLUS_WORKER_TIME" \
  --export=ALL,CONFIG_ENV="$CONFIG_ENV",AMRFINDERPLUS_SAMPLE_FILE="$AMRFINDERPLUS_SAMPLE_FILE",AMRFINDERPLUS_RUN_OUT="$RUN_OUT",AMRFINDERPLUS_ORGANISM_RESOLVED="$AMRFINDERPLUS_ORGANISM_RESOLVED" \
  "$WORKER")
echo "[INFO] submitted array=$ARRAY_JOBID"
while squeue -h -j "$ARRAY_JOBID" | grep -q .; do
  sleep 30
done

echo "[INFO] amrfinderplus worker array completed"
DONE_FILE="$RUN_OUT/amrfinderplus_done_from_outputs_all.txt"
for d in "$RUN_OUT"/*; do
  [[ -d "$d" ]] || continue
  s="$(basename "$d")"
  [[ "$s" == "status" || "$s" == "tmp" || "$s" == "meta" ]] && continue
  [[ -s "$d/$s.amrfinder.tsv" ]] && echo "$s"
done | sort > "$DONE_FILE"

MERGED_TSV="$RUN_OUT/amrfinderplus_all_samples.tsv"
: > "$MERGED_TSV"
header_written=0
for d in "$RUN_OUT"/*; do
  [[ -d "$d" ]] || continue
  s="$(basename "$d")"
  [[ "$s" == "status" || "$s" == "tmp" || "$s" == "meta" ]] && continue
  f="$d/$s.amrfinder.tsv"
  [[ -s "$f" ]] || continue
  if [[ "$header_written" == "0" ]]; then
    awk -v OFS='\t' -v sample="$s" 'NR==1{print "sample_id",$0; next} {print sample,$0}' "$f" > "$MERGED_TSV"
    header_written=1
  else
    awk -v OFS='\t' -v sample="$s" 'NR>1{print sample,$0}' "$f" >> "$MERGED_TSV"
  fi
done

refresh_symlink "$RUN_OUT" "$AMRFINDERPLUS_OUTROOT/latest"
echo "[OK] done_file=$DONE_FILE"
echo "[OK] done_count=$(wc -l < "$DONE_FILE")"
echo "[OK] merged_tsv=$MERGED_TSV"
echo "[OK] amrfinderplus_latest=$AMRFINDERPLUS_OUTROOT/latest"
