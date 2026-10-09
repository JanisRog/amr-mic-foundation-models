#!/usr/bin/env bash
#SBATCH --job-name=ardap_submit
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail

PSEAER="${PSEAER:?set PSEAER in config/paths.env}"
PIPE="${PIPE:?set PIPE in config/paths.env}"
READS_ROOT="${READS_ROOT:-$PSEAER/raw_reads/fastq}"
OUTROOT="${OUTROOT:-$PSEAER/results/ardap/all_samples}"
WORKER="${WORKER:-$PIPE/scripts/run_ardap_pseaer_worker.sh}"
SAMPLE_LIST="${SAMPLE_LIST:-$OUTROOT/sample_ids.with_pairs.txt}"
LOGDIR="${LOGDIR:-$PSEAER/00_logs/ardap/workers}"
MAX_PARALLEL="${MAX_PARALLEL:-50}"
WORKER_CPUS="${WORKER_CPUS:-8}"
WORKER_MEM="${WORKER_MEM:-32G}"
WORKER_TIME="${WORKER_TIME:-12:00:00}"
ARDAP_SIZE="${ARDAP_SIZE:-1000000}"
FORCE="${FORCE:-0}"

[[ -d "$READS_ROOT" ]] || { echo "[ERROR] missing reads root: $READS_ROOT" >&2; exit 1; }
[[ -f "$WORKER" ]] || { echo "[ERROR] missing worker: $WORKER" >&2; exit 1; }
mkdir -p "$OUTROOT" "$LOGDIR"

find "$READS_ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' |
while IFS= read -r sample; do
  [[ -s "$READS_ROOT/$sample/${sample}_1.fastq.gz" ]] || continue
  [[ -s "$READS_ROOT/$sample/${sample}_2.fastq.gz" ]] || continue
  printf '%s\n' "$sample"
done | sort > "$SAMPLE_LIST"

n_samples=$(wc -l < "$SAMPLE_LIST")
[[ "$n_samples" -gt 0 ]] || { echo "[ERROR] no paired samples found" >&2; exit 1; }

jid=$(sbatch --parsable \
  --job-name=ardap_pseaer \
  --partition=standard \
  --array="1-${n_samples}%${MAX_PARALLEL}" \
  --cpus-per-task="$WORKER_CPUS" \
  --mem="$WORKER_MEM" \
  --time="$WORKER_TIME" \
  --output="$LOGDIR/ardap_%A_%a.out" \
  --error="$LOGDIR/ardap_%A_%a.err" \
  --export=ALL,PSEAER="$PSEAER",READS_ROOT="$READS_ROOT",OUTROOT="$OUTROOT",SAMPLE_LIST="$SAMPLE_LIST",ARDAP_SIZE="$ARDAP_SIZE",FORCE="$FORCE" \
  --wrap='SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST"); export SAMPLE_ID; exec "'"$WORKER"'"')

echo "[OK] ardap_array=$jid samples=$n_samples"
echo "[OK] sample_list=$SAMPLE_LIST"
echo "[OK] outputs=$OUTROOT/samples"
echo "[OK] logs=$LOGDIR/ardap_${jid}_*.out"
