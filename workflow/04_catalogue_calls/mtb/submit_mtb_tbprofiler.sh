#!/usr/bin/env bash
#SBATCH --job-name=mtb_tbprof_submit
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00

set -euo pipefail

MTB_ROOT="${MTB_ROOT:-${MTB:?set MTB in config/paths.env}}"
READS_ROOT="${READS_ROOT:-$MTB_ROOT/reads}"
OUTROOT="${OUTROOT:-$MTB_ROOT/results/tbprofiler/reads}"
TBPROFILER_SIF="${TBPROFILER_SIF:?set TBPROFILER_SIF in config/paths.env}"
SAMPLE_LIST="${SAMPLE_LIST:-$READS_ROOT/paired_sample_ids.txt}"
SAMPLE_WORKER="${SAMPLE_WORKER:-$MTB_ROOT/scripts/run_tbprofiler_sample.sh}"
CHUNK_WORKER="${CHUNK_WORKER:-$MTB_ROOT/scripts/run_tbprofiler_chunk.sh}"
LOGDIR="${LOGDIR:-$MTB_ROOT/00_logs/tbprofiler/workers}"
SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-10}"
MAX_PARALLEL="${MAX_PARALLEL:-50}"
WORKER_CPUS="${WORKER_CPUS:-4}"
WORKER_MEM="${WORKER_MEM:-12G}"
WORKER_TIME="${WORKER_TIME:-24:00:00}"
KEEP_BAM="${KEEP_BAM:-0}"
FORCE="${FORCE:-0}"
CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-}"

for path in "$TBPROFILER_SIF" "$SAMPLE_LIST" "$SAMPLE_WORKER" "$CHUNK_WORKER"; do
  [[ -e "$path" ]] || { echo "[ERROR] missing file: $path" >&2; exit 1; }
done

mkdir -p "$OUTROOT/samples" "$LOGDIR"
n_samples=$(wc -l < "$SAMPLE_LIST")
[[ "$n_samples" -gt 0 ]] || { echo "[ERROR] empty sample list: $SAMPLE_LIST" >&2; exit 1; }
n_tasks=$(( (n_samples + SAMPLES_PER_TASK - 1) / SAMPLES_PER_TASK ))

jid=$(sbatch --parsable \
  --job-name=mtb_tbprofiler \
  --partition=standard \
  --array="1-${n_tasks}%${MAX_PARALLEL}" \
  --cpus-per-task="$WORKER_CPUS" \
  --mem="$WORKER_MEM" \
  --time="$WORKER_TIME" \
  --output="$LOGDIR/mtb_tbprofiler_%A_%a.out" \
  --error="$LOGDIR/mtb_tbprofiler_%A_%a.err" \
  --export=ALL,MTB_ROOT="$MTB_ROOT",READS_ROOT="$READS_ROOT",OUTROOT="$OUTROOT",TBPROFILER_SIF="$TBPROFILER_SIF",SAMPLE_LIST="$SAMPLE_LIST",SAMPLE_WORKER="$SAMPLE_WORKER",SAMPLES_PER_TASK="$SAMPLES_PER_TASK",THREADS="$WORKER_CPUS",KEEP_BAM="$KEEP_BAM",FORCE="$FORCE",CONTAINER_RUNTIME="$CONTAINER_RUNTIME" \
  "$CHUNK_WORKER")

echo "[OK] tbprofiler_array=$jid samples=$n_samples tasks=$n_tasks"
echo "[OK] outputs=$OUTROOT/samples"
echo "[OK] logs=$LOGDIR/mtb_tbprofiler_${jid}_*.out"
