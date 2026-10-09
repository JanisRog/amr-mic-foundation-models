#!/usr/bin/env bash
set -euo pipefail

: "${SAMPLE_ID:?set SAMPLE_ID}"
: "${MTB:?set MTB}"

SNIPPY_SIF="${SNIPPY_SIF:?set SNIPPY_SIF in config/paths.env}"
REF="${REF:-$MTB/data/reference/tbprofiler_tbdb_72ef6fa/tbdb.fasta}"
OUTROOT="${OUTROOT:-$MTB/results/gwas/snippy_h37rv}"
THREADS="${THREADS:-8}"
FORCE="${FORCE:-0}"
MANIFEST="$MTB/reads/manifests/${SAMPLE_ID}.fastq_manifest.tsv"
OUT="$OUTROOT/runs/$SAMPLE_ID"

REQUIRED_OUTPUTS=(snps.vcf snps.tab snps.consensus.fa snps.aligned.fa)

missing_outputs() {
  local f out=""
  for f in "${REQUIRED_OUTPUTS[@]}"; do
    [[ -s "$OUT/$f" ]] || out="${out:+$out,}$f"
  done
  printf '%s' "$out"
}

if [[ "$FORCE" != "1" && -z "$(missing_outputs)" ]]; then
  echo "[SKIP] sample=$SAMPLE_ID"
  exit 0
fi

module load apptainer >/dev/null 2>&1 || true
RUNTIME="$(command -v apptainer || command -v singularity)"
[[ -x "$RUNTIME" ]] || { echo "[ERROR] no container runtime" >&2; exit 2; }
[[ -s "$SNIPPY_SIF" ]] || { echo "[ERROR] missing snippy sif: $SNIPPY_SIF" >&2; exit 2; }
[[ -s "$REF" ]] || { echo "[ERROR] missing ref: $REF" >&2; exit 2; }
[[ -s "$MANIFEST" ]] || { echo "[ERROR] missing manifest: $MANIFEST" >&2; exit 2; }

tmp="${SLURM_TMPDIR:-/tmp}/snippy_${SLURM_JOB_ID:-local}_${SAMPLE_ID}"
rm -rf "$tmp"
mkdir -p "$tmp" "$OUTROOT/runs"
trap 'rm -rf "$tmp"' EXIT

mapfile -t r1s < <(awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++){gsub(/\r/,"",$i); h[$i]=i}; next}
  {for(i=1;i<=NF;i++) gsub(/\r/,"",$i)}
  $(h["read_index"])=="1" {print $NF}
' "$MANIFEST")
mapfile -t r2s < <(awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++){gsub(/\r/,"",$i); h[$i]=i}; next}
  {for(i=1;i<=NF;i++) gsub(/\r/,"",$i)}
  $(h["read_index"])=="2" {print $NF}
' "$MANIFEST")
[[ "${#r1s[@]}" -gt 0 && "${#r1s[@]}" -eq "${#r2s[@]}" ]] || { echo "[ERROR] invalid pairs sample=$SAMPLE_ID" >&2; exit 2; }

combine_fastqs_gz() {
  local output="$1"; shift
  {
    local f magic
    for f in "$@"; do
      [[ -s "$f" ]] || { echo "[ERROR] missing/empty FASTQ: $f" >&2; exit 2; }
      magic="$(od -An -N2 -tx1 "$f" | tr -d '[:space:]')"
      [[ "$magic" == 1f8b ]] && gzip -cd -- "$f" || cat -- "$f"
    done
  } | gzip -c > "$output"
}
combine_fastqs_gz "$tmp/R1.fastq.gz" "${r1s[@]}"
combine_fastqs_gz "$tmp/R2.fastq.gz" "${r2s[@]}"

rm -rf "$OUT"
unset SINGULARITY_BINDPATH APPTAINER_BINDPATH || true
"$RUNTIME" exec --cleanenv --bind "${BIND_PATHS:?set BIND_PATHS in config/paths.env}","$tmp":"$tmp" "$SNIPPY_SIF" \
  snippy \
    --force \
    --cpus "$THREADS" \
    --outdir "$OUT" \
    --ref "$REF" \
    --R1 "$tmp/R1.fastq.gz" \
    --R2 "$tmp/R2.fastq.gz" \
    --tmpdir "$tmp"

miss="$(missing_outputs)"
[[ -z "$miss" ]] || { echo "[ERROR] incomplete output sample=$SAMPLE_ID missing=$miss" >&2; exit 2; }
echo "[OK] snippy sample=$SAMPLE_ID out=$OUT"
