#!/bin/bash
#SBATCH --job-name=cat_only
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# One drug's development CV on catalogue features alone -- no embeddings.
#
# The feature-set ablation (§20) needs two arms without a sequence model:
# catalogued mutation columns alone, and the same columns with TB-Profiler's call.
# They are scored on the drug's own development fold file, so every isolate and
# partition is the one the embedding arms were scored on, and the comparison is
# paired. Same 1,000 configurations, same median aggregation, same breakpoint, same
# seed; the only difference from an embedding arm is the absence of the cache.
#
# Submitted by submit_mtb_catalogue_arm.sh with CATALOGUE_ONLY=1, which builds the
# merged matrix and passes it as MATRIX. Run directory: runs/catalogue_<TAG>.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
NJOBS="${NJOBS:-16}"

: "${DRUG:?set DRUG}" "${CVDIR:?set CVDIR}" "${MATRIX:?set MATRIX}" \
  "${SR_BOUNDARY:?set SR_BOUNDARY}" "${TAG:?set TAG}"
D="$ROOT/$CVDIR"
# Development sweeps use the full grid. A refit on other folds (the lineage test,
# §21) passes the development-selected config instead, so it is the same model
# and not a re-selection. CVDIR may be any fold folder under ROOT, e.g.
# cv_<drug>/lineage or cv_<drug>/lineage/sizematched.
CONFIGS="${CONFIGS:-$SCRIPTS/mtb_rf1000_configs.json}"
RUNNAME="${RUNNAME:-catalogue_${TAG}}"

for f in "$D/cv_folds.tsv" "$MATRIX" "$CONFIGS"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done
n_cols=$(( $(head -1 "$MATRIX" | tr '\t' '\n' | wc -l) - 1 ))
echo "[INFO] drug     = $DRUG   folds = $D/cv_folds.tsv"
echo "[INFO] matrix   = $MATRIX ($n_cols feature columns)"
echo "[INFO] boundary = $SR_BOUNDARY   tag = $TAG   configs = $CONFIGS"

$PY "$SCRIPTS/cv_run_architecture.py" \
  --feature-tsv "$MATRIX" \
  --folds-tsv "$D/cv_folds.tsv" \
  --configs "$CONFIGS" \
  --sr-boundary "$SR_BOUNDARY" \
  --label "$RUNNAME" \
  --out-dir "$D/runs/$RUNNAME" \
  --n-jobs "$NJOBS" --seed 2026

echo "[DONE] $D/runs/$RUNNAME"
