#!/usr/bin/env bash
#SBATCH --job-name=mtb_reads_submit
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail

MTB_ROOT="${MTB_ROOT:-${MTB:?set MTB in config/paths.env}}"
SOURCE_METADATA_TSV="${SOURCE_METADATA_TSV:-$MTB_ROOT/metadata/mtb_metadata_ena_samples.tsv}"
OUTROOT="${OUTROOT:-$MTB_ROOT/reads}"
TODO_TSV="${TODO_TSV:-$MTB_ROOT/metadata/mtb_metadata_ena_reads.todo.tsv}"
BUILDER="${BUILDER:-$MTB_ROOT/scripts/build_mtb_reads_missing_list.py}"
WORKER="${WORKER:-$MTB_ROOT/scripts/download_mtb_ena_reads_worker.py}"
CHUNK_WORKER="${CHUNK_WORKER:-$MTB_ROOT/scripts/download_mtb_ena_reads_chunk_worker.sh}"
LOGDIR="${LOGDIR:-$MTB_ROOT/00_logs/reads/workers}"
SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-25}"
MAX_PARALLEL="${MAX_PARALLEL:-40}"
WORKER_MEM="${WORKER_MEM:-4G}"
WORKER_TIME="${WORKER_TIME:-24:00:00}"

for path in "$SOURCE_METADATA_TSV" "$BUILDER" "$WORKER" "$CHUNK_WORKER"; do
  [[ -f "$path" ]] || { echo "[ERROR] missing file: $path" >&2; exit 1; }
done

mkdir -p "$OUTROOT/fastq_by_sample" "$OUTROOT/status" "$OUTROOT/manifests" "$LOGDIR"
python3 "$BUILDER" "$SOURCE_METADATA_TSV" "$OUTROOT" "$TODO_TSV"

n_samples=$(( $(wc -l < "$TODO_TSV") - 1 ))
if (( n_samples == 0 )); then
  echo "[OK] no missing read samples"
  exit 0
fi
n_tasks=$(( (n_samples + SAMPLES_PER_TASK - 1) / SAMPLES_PER_TASK ))

jid=$(sbatch --parsable \
  --job-name=mtb_ena_reads_missing \
  --partition=standard \
  --array="1-${n_tasks}%${MAX_PARALLEL}" \
  --cpus-per-task=1 \
  --mem="$WORKER_MEM" \
  --time="$WORKER_TIME" \
  --output="$LOGDIR/mtb_ena_reads_%A_%a.out" \
  --error="$LOGDIR/mtb_ena_reads_%A_%a.err" \
  --export=ALL,METADATA_TSV="$TODO_TSV",OUTROOT="$OUTROOT",WORKER="$WORKER",SAMPLES_PER_TASK="$SAMPLES_PER_TASK" \
  "$CHUNK_WORKER")

echo "[OK] job=$jid samples=$n_samples tasks=$n_tasks samples_per_task=$SAMPLES_PER_TASK"
echo "[OK] logs=$LOGDIR/mtb_ena_reads_${jid}_*.out"
