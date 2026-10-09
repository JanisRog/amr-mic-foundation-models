#!/usr/bin/env bash

set -euo pipefail

: "${SAMPLE_ID:?set SAMPLE_ID}"

PSEAER="${PSEAER:?set PSEAER in config/paths.env}"
ARDAP="${ARDAP:?set ARDAP in config/paths.env}"
ARDAP_ENV="${ARDAP_ENV:?set ARDAP_ENV in config/paths.env}"
JAVA_ENV="${JAVA_ENV:?set JAVA_ENV in config/paths.env}"
PICARD_ENV="${PICARD_ENV:?set PICARD_ENV in config/paths.env}"
READS_ROOT="${READS_ROOT:-$PSEAER/raw_reads/fastq}"
OUTROOT="${OUTROOT:-$PSEAER/results/ardap/all_samples}"
NXF_VER="${NXF_VER:-22.10.8}"
SIZE="${ARDAP_SIZE:-1000000}"
FORCE="${FORCE:-0}"

export PATH="$ARDAP/bin:$PICARD_ENV/bin:$ARDAP_ENV/bin:$PATH"
export JAVA_CMD="$JAVA_ENV/bin/java"
export NXF_OPTS="-XX:+UseSerialGC -XX:ActiveProcessorCount=${SLURM_CPUS_PER_TASK:-8} -Xms512m -Xmx4g"
export NXF_VER
export NXF_HOME="$OUTROOT/.nextflow"

sample_reads="$READS_ROOT/$SAMPLE_ID"
sample_out="$OUTROOT/samples/$SAMPLE_ID"
status_dir="$OUTROOT/status"
done_file="$status_dir/$SAMPLE_ID.done"

if [[ "$FORCE" != "1" && -s "$done_file" ]]; then
  echo "[SKIP] completed sample=$SAMPLE_ID"
  exit 0
fi
[[ -d "$sample_reads" ]] || { echo "[ERROR] missing reads directory: $sample_reads" >&2; exit 2; }

r1="$sample_reads/${SAMPLE_ID}_1.fastq.gz"
r2="$sample_reads/${SAMPLE_ID}_2.fastq.gz"
[[ -s "$r1" ]] || { echo "[ERROR] missing R1: $r1" >&2; exit 2; }
[[ -s "$r2" ]] || { echo "[ERROR] missing R2: $r2" >&2; exit 2; }

for path in "$ARDAP/main.nf" "$JAVA_CMD"; do
  [[ -e "$path" ]] || { echo "[ERROR] missing required file: $path" >&2; exit 2; }
done
for command_name in nextflow bwa samtools picard bedtools gatk mosdepth bc; do
  command -v "$command_name" >/dev/null 2>&1 || {
    echo "[ERROR] required ARDaP command not found: $command_name" >&2
    exit 2
  }
done

mkdir -p "$sample_out/input" "$sample_out/work" "$status_dir" "$NXF_HOME"
ln -sfn "$r1" "$sample_out/input/${SAMPLE_ID}_1.fastq.gz"
ln -sfn "$r2" "$sample_out/input/${SAMPLE_ID}_2.fastq.gz"

echo "[INFO] sample=$SAMPLE_ID"
echo "[INFO] sample_out=$sample_out"

cd "$sample_out"
nextflow run "$ARDAP/main.nf" \
  -work-dir "$sample_out/work" \
  --executor local \
  --species Pseudomonas_aeruginosa \
  --fastq "$sample_out/input/*_{1,2}.fastq.gz" \
  --size "$SIZE" \
  -resume

[[ -d "$sample_out/Outputs/AbR_reports" ]] || {
  echo "[ERROR] missing resistance report directory after completion" >&2
  exit 2
}

printf 'sample_id\tstatus\toutput\n%s\tcomplete\t%s\n' \
  "$SAMPLE_ID" "$sample_out/Outputs" > "$done_file"
echo "[OK] ARDaP sample=$SAMPLE_ID outputs=$sample_out/Outputs"
