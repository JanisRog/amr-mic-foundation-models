#!/usr/bin/env bash
#SBATCH --job-name=mtb_tbp_matrix
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
SCRIPT="${SCRIPT:-$MTB/scripts/build_tbprofiler_amrfinder_style_matrix.py}"
VARIANTS_CSV="${VARIANTS_CSV:-$MTB/results/tbprofiler/reads/gather/tbprofiler.variants.csv}"
SAMPLE_LIST="${SAMPLE_LIST:-$MTB/results/tbprofiler/reads/tbprofiler_complete_sample_ids.txt}"
OUTDIR="${OUTDIR:-$MTB/results/tbprofiler/reads/gather/matrix}"

python3 "$SCRIPT" \
  --variants-csv "$VARIANTS_CSV" \
  --sample-list "$SAMPLE_LIST" \
  --outdir "$OUTDIR"
