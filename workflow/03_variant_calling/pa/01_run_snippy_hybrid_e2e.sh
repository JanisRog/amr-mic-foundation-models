#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_shared_gwas.sh"
setup_stage_logging "01_run_snippy_hybrid_e2e"

require_file "$PYTHON_BIN"
require_file "$REFERENCE_FASTA"
require_file "$SNIPPY_SIF"
require_file "$INPUT_PREP_ROOT/latest/sample_manifest.cleaned.tsv"
require_file "$SPLIT_DIR/latest/eligible_samples_manifest.tsv"
WORKER="$GWAS_PIPELINE_ROOT/wrappers/01a_snippy_hybrid_worker.sh"
require_file "$WORKER"
CTR_BIN="$(resolve_container_bin)"

STAMP="$(date +%Y%m%d_%H%M%S)"
mkdir -p "$SNIPPY_ROOT" "$SNP_ML_ROOT"
RUN_ROOT="$SNIPPY_ROOT/run_${STAMP}"
SAMPLE_TSV="$RUN_ROOT/snippy_samples.tsv"
CORE_DIR="$RUN_ROOT/core"
RUNS_DIR="$RUN_ROOT/runs"
STATUS_DIR="$RUN_ROOT/status"
mkdir -p "$RUN_ROOT" "$RUNS_DIR" "$CORE_DIR" "$STATUS_DIR"

cut -f1,3,4,5 "$SPLIT_DIR/latest/eligible_samples_manifest.tsv" > "$SAMPLE_TSV"
N=$(( $(wc -l < "$SAMPLE_TSV") - 1 ))
(( N >= 1 )) || { echo "[ERROR] no eligible samples for snippy" >&2; exit 1; }

echo "[INFO] reference_fasta=$REFERENCE_FASTA"
echo "[INFO] snippy_input_mode=$SNIPPY_INPUT_MODE"
echo "[INFO] sample_tsv=$SAMPLE_TSV"
echo "[INFO] run_root=$RUN_ROOT"
echo "[INFO] n_samples=$N"
echo "[INFO] worker_cpus=$SNIPPY_WORKER_CPUS"
echo "[INFO] worker_mem=$SNIPPY_WORKER_MEM"
echo "[INFO] worker_time=$SNIPPY_WORKER_TIME"

ARRAY_JOBID=$(sbatch --parsable --array=1-"$N"%"$SNIPPY_MAX_PARALLEL" \
  --cpus-per-task="$SNIPPY_WORKER_CPUS" \
  --mem="$SNIPPY_WORKER_MEM" \
  --time="$SNIPPY_WORKER_TIME" \
  --export=ALL,CONFIG_ENV="$CONFIG_ENV",SNIPPY_SAMPLE_TSV="$SAMPLE_TSV",SNIPPY_RUN_ROOT="$RUN_ROOT" \
  "$WORKER")
echo "[INFO] submitted array=$ARRAY_JOBID"

while true; do
  N_DONE=$(find "$STATUS_DIR" -maxdepth 1 -type f -name '*.done' | wc -l | tr -d ' ')
  if (( N_DONE >= N )); then
    break
  fi
  sleep 30
done

echo "[INFO] worker array completed"
echo "[INFO] worker_status_markers=$N_DONE/$N"
if find "$STATUS_DIR" -maxdepth 1 -type f -name '*.done' | grep -q .; then
  echo "[INFO] worker_status_summary:"
  awk -F'\t' 'FNR==3{c[$2]++} END{for(k in c) print "  " k "=" c[k]}' "$STATUS_DIR"/*.done | sort
fi
find "$RUNS_DIR" -mindepth 1 -maxdepth 1 -type d \
  -exec test -s "{}/snps.tab" \; -exec test -s "{}/snps.vcf" \; -print | sort > "$CORE_DIR/snippy_run_dirs.list"
N_OK=$(wc -l < "$CORE_DIR/snippy_run_dirs.list" | tr -d ' ')
echo "[INFO] successful_runs=$N_OK"
(( N_OK >= 2 )) || { echo "[ERROR] need at least 2 successful snippy runs" >&2; exit 1; }

unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"
"$CTR_BIN" exec --cleanenv -B "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$SNIPPY_SIF" snippy-core \
  --ref "$REFERENCE_FASTA" \
  --prefix "$CORE_DIR/core" \
  $(cat "$CORE_DIR/snippy_run_dirs.list") \
  > "$CORE_DIR/snippy-core.stdout.log" 2> "$CORE_DIR/snippy-core.stderr.log"

require_file "$CORE_DIR/core.vcf"
require_file "$CORE_DIR/core.full.aln"
refresh_symlink "$RUN_ROOT" "$SNIPPY_ROOT/latest"
refresh_symlink "$CORE_DIR" "$SNP_ML_ROOT/latest"

echo "[OK] snippy_latest=$SNIPPY_ROOT/latest"
echo "[OK] core_vcf=$CORE_DIR/core.vcf"
echo "[OK] core_full_aln=$CORE_DIR/core.full.aln"
