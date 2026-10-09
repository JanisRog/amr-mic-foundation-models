#!/usr/bin/env bash
#SBATCH --job-name=mtb_snippy_submit
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
OUTROOT="${OUTROOT:-$MTB/results/gwas/snippy_h37rv}"
DEFAULT_SAMPLE_LIST="$MTB/results/tbprofiler/reads/tbprofiler_complete_sample_ids.txt"
if [[ ! -s "$DEFAULT_SAMPLE_LIST" ]]; then
  DEFAULT_SAMPLE_LIST="$MTB/reads/paired_sample_ids.txt"
fi
SAMPLE_LIST="${SAMPLE_LIST:-$DEFAULT_SAMPLE_LIST}"
SAMPLE_WORKER="${SAMPLE_WORKER:-$MTB/scripts/run_mtb_snippy_sample.sh}"
CHUNK_WORKER="${CHUNK_WORKER:-$MTB/scripts/run_mtb_snippy_chunk.sh}"
LOGDIR="${LOGDIR:-$MTB/00_logs/snippy_h37rv/workers}"
SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-1}"
MAX_PARALLEL="${MAX_PARALLEL:-50}"
WORKER_CPUS="${WORKER_CPUS:-8}"
WORKER_MEM="${WORKER_MEM:-24G}"
WORKER_TIME="${WORKER_TIME:-01:00:00}"
FORCE="${FORCE:-0}"

[[ -s "$SAMPLE_LIST" ]] || { echo "[ERROR] missing or empty SAMPLE_LIST: $SAMPLE_LIST" >&2; exit 2; }
[[ -s "$SAMPLE_WORKER" ]] || { echo "[ERROR] missing SAMPLE_WORKER: $SAMPLE_WORKER" >&2; exit 2; }
[[ -s "$CHUNK_WORKER" ]] || { echo "[ERROR] missing CHUNK_WORKER: $CHUNK_WORKER" >&2; exit 2; }

mkdir -p "$OUTROOT/runs" "$LOGDIR"
n=$(grep -cve '^[[:space:]]*$' "$SAMPLE_LIST")
tasks=$(( (n + SAMPLES_PER_TASK - 1) / SAMPLES_PER_TASK ))
echo "[INFO] sample_list=$SAMPLE_LIST samples=$n per_task=$SAMPLES_PER_TASK tasks=$tasks"

jid=$(sbatch --parsable \
  --job-name=mtb_snippy \
  --partition=standard \
  --array="1-${tasks}%${MAX_PARALLEL}" \
  --cpus-per-task="$WORKER_CPUS" \
  --mem="$WORKER_MEM" \
  --time="$WORKER_TIME" \
  --output="$LOGDIR/mtb_snippy_%A_%a.out" \
  --error="$LOGDIR/mtb_snippy_%A_%a.err" \
  --export=ALL,MTB="$MTB",OUTROOT="$OUTROOT",SAMPLE_LIST="$SAMPLE_LIST",SAMPLE_WORKER="$SAMPLE_WORKER",SAMPLES_PER_TASK="$SAMPLES_PER_TASK",THREADS="$WORKER_CPUS",FORCE="$FORCE" \
  "$CHUNK_WORKER")

echo "[OK] snippy_array=$jid samples=$n tasks=$tasks outroot=$OUTROOT"
