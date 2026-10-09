#!/bin/bash
#SBATCH --job-name=pa_extract
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-2126%64

# Rebuild each P. aeruginosa Snippy VCF at allele fraction 0.1 and extract the
# 62 curated target sequences from PAO1.
#
# The AF 0.1 rebuild is the MTB panel's section 11 finding, ported: rebuilding
# the consensus at 0.1 instead of freebayes' ploidy-2 call improved all twelve
# isoniazid comparisons and moved the catalogue rule's own holdout cost from
# 0.166 to 0.139. It was never testable on assemblies -- these isolates have
# read-level VCFs with AO/DP at every site, so it is.
#
# One array task per isolate, two steps:
#   snps.raw.vcf   -> snps.af010.vcf.gz   (QUAL>=100, DP>=10, AO/DP>=0.1)
#   snps.af010.vcf -> <sample>.targets.tsv + <sample>.sequences.fasta
#
# `snps.raw.vcf` is freebayes' pre-filter output and is the only file carrying
# AO/DP at every site; `snps.vcf` has already been reduced to 1/1 genotypes,
# which is the 0.9-equivalent behaviour this is replacing. Reading the wrong one
# produces a plausible file with no minority calls in it.

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
SCRIPTS="${SCRIPTS:-$P/scripts}"
PY="${PY:?set PY in config/paths.env}"
RUNS="${RUNS:-$P/results/gwas/snippy_hybrid/latest/runs}"
REF="${REF:-$P/results/bakta/reference_latest/PAO1_NC_002516.2.fna}"
CLUMPS="${CLUMPS:-$P/targets/all_gene_target_clumps.tsv}"
TLIST="${TLIST:-$P/targets/union_targets.txt}"
OUT="${OUT:-$P/results/target_sequences/af010}"
AF="${AF:-0.1}"
SAMPLES="${SAMPLES:-$P/results/gwas/snippy_hybrid/latest/all_samples_snippy_ok.txt}"

for f in "$REF" "$CLUMPS" "$TLIST" "$SAMPLES"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done

S=$(awk -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==n{print $1}' "$SAMPLES")
[ -n "$S" ] || { echo "[ERROR] no sample at row $((SLURM_ARRAY_TASK_ID + 1))"; exit 1; }
D="$RUNS/$S"
RAW="$D/snps.raw.vcf"
[ -e "$RAW" ] || { echo "[ERROR] $S: no snps.raw.vcf"; exit 1; }

mkdir -p "$OUT"
# MTB names these af010 / af090 (AF x 100, zero-padded). Keep both organisms
# on one convention so a path can be read the same way in either project.
AFVCF="$D/snps.af$(printf '%03d' "$(echo "$AF * 100 / 1" | bc)").vcf.gz"
echo "[INFO] $S  ->  $AFVCF"

# The filter writes plain text whatever the output is called, so it writes to
# the bare .vcf name and the file is compressed afterwards. Naming a text file
# .vcf.gz fails later, in whatever reads it, rather than here.
PLAIN="${AFVCF%.gz}"
if [ ! -s "$AFVCF" ]; then
  $PY "$SCRIPTS/filter_vcf_by_allele_fraction.py" \
    --in-vcf "$RAW" --out-vcf "$PLAIN" --threshold "$AF" \
    --stats "$OUT/$S.afstats.tsv"
  gzip -f "$PLAIN"
  [ -s "$AFVCF" ] || { echo "[ERROR] $S: gzip produced no $AFVCF"; exit 1; }
else
  echo "[INFO] reusing existing $AFVCF"
fi

$PY "$SCRIPTS/extract_pseaer_target_sequences.py" \
  --reference "$REF" --clumps "$CLUMPS" --target-list "$TLIST" \
  --vcf "$AFVCF" --sample-id "$S" --out-dir "$OUT"
