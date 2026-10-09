#!/usr/bin/env bash
#SBATCH --job-name=mtb_tbp_target_gather
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PIPE="${PIPE:?set PIPE in config/paths.env}"
OUTROOT="${OUTROOT:-$MTB/results/gwas/sample_windows/tbprofiler_all_targets_gene_up300}"
SAMPLE_LIST="${SAMPLE_LIST:-$MTB/results/tbprofiler/reads/tbprofiler_complete_sample_ids.txt}"
SCRIPT="${SCRIPT:-$PIPE/scripts/gather_tbprofiler_target_sequences.py}"

python3 "$SCRIPT" \
  --samples-root "$OUTROOT/samples" \
  --sample-list "$SAMPLE_LIST" \
  --outdir "$OUTROOT/gather" \
  --expected-rows-per-sample 145
