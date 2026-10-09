#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_shared_gwas.sh"
setup_stage_logging "00c_quast"

require_dir "$ASSEMBLIES_DIR"
[[ -e "$QUAST_IMG" ]] || { echo "[ERROR] missing QUAST container/image: $QUAST_IMG" >&2; exit 1; }

STAMP="$(date +%Y%m%d_%H%M%S)"
mkdir -p "$QC_QUAST_OUTROOT"
RUN_OUT="$QC_QUAST_OUTROOT/run_${STAMP}"
SAMPLE_TSV="$RUN_OUT/quast_samples.tsv"
mkdir -p "$RUN_OUT"

{
  printf "sample\tassembly\n"
  find "$ASSEMBLIES_DIR" -maxdepth 1 -type f \( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' \) \
    | sort \
    | while read -r asm; do
        sample="$(basename "$asm")"
        sample="${sample%.fasta}"
        sample="${sample%.fa}"
        sample="${sample%.fna}"
        printf "%s\t%s\n" "$sample" "$asm"
      done
} > "$SAMPLE_TSV"

N=$(($(wc -l < "$SAMPLE_TSV") - 1))
if (( N < 1 )); then
  echo "[ERROR] no assemblies found in $ASSEMBLIES_DIR" >&2
  exit 1
fi

WORKER="$GWAS_PIPELINE_ROOT/wrappers/00c_quast_worker.sh"
require_file "$WORKER"

echo "[INFO] assemblies_dir=$ASSEMBLIES_DIR"
echo "[INFO] run_out=$RUN_OUT"
echo "[INFO] sample_tsv=$SAMPLE_TSV"
echo "[INFO] n_samples=$N"
echo "[INFO] worker_cpus=$QUAST_WORKER_CPUS"
echo "[INFO] worker_mem=$QUAST_WORKER_MEM"
echo "[INFO] worker_time=$QUAST_WORKER_TIME"

ARRAY_JOBID=$(sbatch --parsable \
  --array=1-"$N"%${QUAST_MAX_PARALLEL} \
  --cpus-per-task="$QUAST_WORKER_CPUS" \
  --mem="$QUAST_WORKER_MEM" \
  --time="$QUAST_WORKER_TIME" \
  --export=ALL,CONFIG_ENV="$CONFIG_ENV",QUAST_RUN_OUT="$RUN_OUT",QUAST_SAMPLE_TSV="$SAMPLE_TSV" \
  "$WORKER")
echo "[INFO] submitted array=$ARRAY_JOBID"

while squeue -h -j "$ARRAY_JOBID" | grep -q .; do
  sleep 30
done
echo "[INFO] worker array completed"

awk 'NR>1 {print $1}' "$SAMPLE_TSV" | sort > "$RUN_OUT/quast_all_samples.txt"
find "$RUN_OUT" -mindepth 1 -maxdepth 1 -type d -name 'quast_*' \
  -exec test -s "{}/report.tsv" \; -print \
  | sed 's#.*/quast_##' | sort > "$RUN_OUT/quast_done_samples.txt"
comm -23 "$RUN_OUT/quast_all_samples.txt" "$RUN_OUT/quast_done_samples.txt" > "$RUN_OUT/quast_missing_samples.txt"

echo "[INFO] successful_runs=$(wc -l < "$RUN_OUT/quast_done_samples.txt")"
echo "[INFO] missing_runs=$(wc -l < "$RUN_OUT/quast_missing_samples.txt")"

refresh_symlink "$RUN_OUT" "$QC_QUAST_OUTROOT/latest"
echo "[OK] quast_latest=$QC_QUAST_OUTROOT/latest"
