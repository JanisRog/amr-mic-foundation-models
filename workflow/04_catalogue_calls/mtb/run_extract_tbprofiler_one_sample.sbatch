#!/usr/bin/env bash
#SBATCH --job-name=mtb_tbp_extract_test
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PIPE="${PIPE:?set PIPE in config/paths.env}"
SAMPLE="${SAMPLE:-site.02.subj.0001.lab.2014222001.iso.1}"
SAMPLE_DIR="${SAMPLE_DIR:-$MTB/results/tbprofiler/reads/samples/$SAMPLE}"
DB_DIR="${DB_DIR:-$MTB/data/reference/tbprofiler_tbdb_72ef6fa}"
OUTDIR="${OUTDIR:-$MTB/results/gwas/sample_windows/tbprofiler_all_targets_gene_up300/test/$SAMPLE}"
SCRIPT="${SCRIPT:-$PIPE/scripts/extract_tbprofiler_target_sequences.py}"

python3 "$SCRIPT" \
  --sample-dir "$SAMPLE_DIR" \
  --db-dir "$DB_DIR" \
  --outdir "$OUTDIR" \
  --upstream-bp 300
