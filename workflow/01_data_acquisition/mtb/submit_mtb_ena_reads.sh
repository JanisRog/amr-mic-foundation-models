#!/usr/bin/env bash
#SBATCH --job-name=mtb_ena_reads_submit
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00

set -euo pipefail

MTB_ROOT="${MTB_ROOT:-${MTB:?set MTB in config/paths.env}}"
METADATA_TSV="${METADATA_TSV:-$MTB_ROOT/assemblies/missing_failed_samples.tsv}"
OUTROOT="${OUTROOT:-$MTB_ROOT/reads}"
WORKER="${WORKER:-$MTB_ROOT/scripts/download_mtb_ena_reads_worker.py}"
LOGDIR="${LOGDIR:-$MTB_ROOT/00_logs/reads}"
MAX_PARALLEL="${MAX_PARALLEL:-20}"
WORKER_MEM="${WORKER_MEM:-8G}"
WORKER_TIME="${WORKER_TIME:-06:00:00}"
FORCE="${FORCE:-0}"

[[ -f "$METADATA_TSV" ]] || { echo "[ERROR] missing metadata: $METADATA_TSV" >&2; exit 1; }
[[ -f "$WORKER" ]] || { echo "[ERROR] missing worker: $WORKER" >&2; exit 1; }
mkdir -p "$OUTROOT/fastq_by_sample" "$OUTROOT/status" "$OUTROOT/manifests" "$OUTROOT/views/fastq_by_sample" "$LOGDIR"
ln -sfn "$OUTROOT/fastq_by_sample" "$OUTROOT/views/fastq_by_sample/latest"

first=$(head -n 1 "$METADATA_TSV" | cut -f1)
N=$(wc -l < "$METADATA_TSV")
if [[ "$first" == "sample_id" || "$first" == "Sample_IDs" ]]; then
  N=$((N - 1))
fi
[[ "$N" -gt 0 ]] || { echo "[ERROR] no metadata rows" >&2; exit 1; }

jid=$(sbatch --parsable \
  --job-name=mtb_ena_reads \
  --partition=standard \
  --array=1-${N}%${MAX_PARALLEL} \
  --cpus-per-task=1 \
  --mem="$WORKER_MEM" \
  --time="$WORKER_TIME" \
  --output="$LOGDIR/mtb_ena_reads_%A_%a.out" \
  --error="$LOGDIR/mtb_ena_reads_%A_%a.err" \
  --export=ALL,METADATA_TSV="$METADATA_TSV",OUTROOT="$OUTROOT",FORCE="$FORCE" \
  "$WORKER")

echo "mtb_ena_reads=$jid"
echo "n_samples=$N"
echo "outroot=$OUTROOT"
echo "fastq_by_sample=$OUTROOT/fastq_by_sample"
echo "latest=$OUTROOT/views/fastq_by_sample/latest"
echo "logs=$LOGDIR/mtb_ena_reads_${jid}_*.out"
