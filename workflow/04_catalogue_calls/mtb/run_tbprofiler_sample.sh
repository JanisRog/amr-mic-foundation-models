#!/usr/bin/env bash

set -euo pipefail

: "${SAMPLE_ID:?set SAMPLE_ID}"
: "${MTB_ROOT:?set MTB_ROOT}"
: "${TBPROFILER_SIF:?set TBPROFILER_SIF}"

READS_ROOT="${READS_ROOT:-$MTB_ROOT/reads}"
OUTROOT="${OUTROOT:-$MTB_ROOT/results/tbprofiler/reads}"
THREADS="${THREADS:-4}"
KEEP_BAM="${KEEP_BAM:-0}"
FORCE="${FORCE:-0}"
CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-}"

resolve_container_runtime() {
  if command -v apptainer >/dev/null 2>&1; then command -v apptainer; return 0; fi
  if command -v singularity >/dev/null 2>&1; then command -v singularity; return 0; fi
  module load apptainer >/dev/null 2>&1 || module load Apptainer >/dev/null 2>&1 || module load singularityce >/dev/null 2>&1 || true
  if command -v apptainer >/dev/null 2>&1; then command -v apptainer; return 0; fi
  if command -v singularity >/dev/null 2>&1; then command -v singularity; return 0; fi
  return 1
}

if [[ -z "$CONTAINER_RUNTIME" ]]; then
  if ! CONTAINER_RUNTIME="$(resolve_container_runtime)"; then
    echo "[ERROR] neither apptainer nor singularity found; set CONTAINER_RUNTIME to its absolute path" >&2
    exit 2
  fi
fi
[[ -x "$CONTAINER_RUNTIME" ]] || {
  echo "[ERROR] container runtime is not executable: $CONTAINER_RUNTIME" >&2
  exit 2
}

manifest="$READS_ROOT/manifests/${SAMPLE_ID}.fastq_manifest.tsv"
sample_out="$OUTROOT/samples/$SAMPLE_ID"
result_json="$sample_out/results/${SAMPLE_ID}.results.json"

[[ -f "$manifest" ]] || { echo "[ERROR] missing manifest: $manifest" >&2; exit 2; }
if [[ "$FORCE" != "1" && -s "$result_json" ]]; then
  echo "[SKIP] completed sample=$SAMPLE_ID result=$result_json"
  exit 0
fi

mkdir -p "$sample_out"
tmp="${SLURM_TMPDIR:-/tmp}/tbprofiler_${SLURM_JOB_ID:-local}_${SAMPLE_ID}"
rm -rf "$tmp"
mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT

# Prevent bcftools and other container tools from inheriting a nonexistent
# login-node scratch path.
export TMPDIR="$tmp"
export TMP="$tmp"
export TEMP="$tmp"
export APPTAINERENV_TMPDIR="$tmp"
export APPTAINERENV_TMP="$tmp"
export APPTAINERENV_TEMP="$tmp"
export SINGULARITYENV_TMPDIR="$tmp"
export SINGULARITYENV_TMP="$tmp"
export SINGULARITYENV_TEMP="$tmp"

mapfile -t r1_files < <(awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++){gsub(/\r/,"",$i); h[$i]=i}; next}
  {for(i=1;i<=NF;i++) gsub(/\r/,"",$i)}
  $(h["read_index"])=="1" {print $(h["fastq"])}
' "$manifest")
mapfile -t r2_files < <(awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++){gsub(/\r/,"",$i); h[$i]=i}; next}
  {for(i=1;i<=NF;i++) gsub(/\r/,"",$i)}
  $(h["read_index"])=="2" {print $(h["fastq"])}
' "$manifest")

[[ "${#r1_files[@]}" -gt 0 && "${#r1_files[@]}" -eq "${#r2_files[@]}" ]] || {
  echo "[ERROR] invalid pairs sample=$SAMPLE_ID r1=${#r1_files[@]} r2=${#r2_files[@]}" >&2
  exit 2
}
for fastq in "${r1_files[@]}" "${r2_files[@]}"; do
  [[ -s "$fastq" ]] || { echo "[ERROR] missing/empty FASTQ: $fastq" >&2; exit 2; }
done

r1="$tmp/${SAMPLE_ID}_1.fastq.gz"
r2="$tmp/${SAMPLE_ID}_2.fastq.gz"

# Some ENA downloads are plain FASTQ despite a .gz suffix. Normalize every
# input stream so Trimmomatic always receives valid gzip-compressed FASTQs.
combine_fastqs_gz() {
  local output="$1"
  shift
  {
    local fastq magic
    for fastq in "$@"; do
      magic="$(od -An -N2 -tx1 "$fastq" | tr -d '[:space:]')"
      if [[ "$magic" == "1f8b" ]]; then
        gzip -cd -- "$fastq"
      else
        cat -- "$fastq"
      fi
    done
  } | gzip -c > "$output"
  gzip -t "$output"
}

combine_fastqs_gz "$r1" "${r1_files[@]}"
combine_fastqs_gz "$r2" "${r2_files[@]}"

echo "[INFO] sample=$SAMPLE_ID pairs=${#r1_files[@]} out=$sample_out"
cd "$sample_out"
unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"
"$CONTAINER_RUNTIME" exec --bind "${BIND_PATHS:?set BIND_PATHS in config/paths.env}","$tmp":"$tmp" "$TBPROFILER_SIF" \
  tb-profiler profile \
  -1 "$r1" \
  -2 "$r2" \
  -p "$SAMPLE_ID" \
  -t "$THREADS" \
  --txt \
  --csv

[[ -s "$result_json" ]] || { echo "[ERROR] missing result after profile: $result_json" >&2; exit 2; }

if [[ "$KEEP_BAM" != "1" ]]; then
  rm -rf "$sample_out/bam"
fi

echo "[OK] tbprofiler sample=$SAMPLE_ID result=$result_json"
