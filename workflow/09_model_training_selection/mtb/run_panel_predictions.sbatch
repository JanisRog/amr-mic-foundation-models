#!/bin/bash
#SBATCH --job-name=panel_pred
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# Re-score one drug's holdout with its declared config and keep the per-isolate
# predictions this time.
#
# The holdout has already been scored; `cv_run_architecture.py` wrote the fold
# metrics and threw the predictions away. A predicted-versus-actual plot needs
# them back, and summary metrics cannot be inverted into them.
#
# This is a re-fit, so it is only sound if it reproduces what was already
# reported. Same cache, same config file, same seed, same boundary, one config
# instead of a thousand -- nothing that enters the fit differs, so the metrics
# recomputed from the saved predictions must match the workbook row. The plotting
# script checks exactly that and refuses to draw if they do not.
#
#   DRUG        isoniazid
#   CODE        inh                (names the rule-bit file and the config)
#   MODEL       evo2               (one model -- the one declared for this drug)
#   SR_BOUNDARY -3.0               (no default: the parameter that fails silently)

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="$MTB/scripts"
NJOBS="${NJOBS:-16}"
TARGETSET="${TARGETSET:-both}"
OUTDIR="${OUTDIR:-$ROOT/predictions}"

: "${DRUG:?set DRUG, e.g. isoniazid}"
: "${CODE:?set CODE, e.g. inh}"
: "${MODEL:?set MODEL, e.g. evo2}"
: "${SR_BOUNDARY:?set SR_BOUNDARY explicitly -- it has no safe default}"

# Isoniazid's fold directory is `cv_isoniazid_final`, not `cv_isoniazid`: it was
# built before this script existed and the name was never reconciled. Overridable
# rather than symlinked, so the real directory stays visible in the log -- and
# note the sibling `cv_isoniazid_rifS_sizematched`, which is the lineage-shift
# experiment and NOT this drug's panel folds.
D="${CVDIR:-$ROOT/cv_${DRUG}}"
FIXED="$ROOT/panel_targets/targets_${DRUG}_${TARGETSET}.txt"
# Isoniazid and ethionamide predate the rule_bit_<code>_tbprofiler.tsv naming
# that started at rifampicin, and isoniazid has three candidate files of
# identical size sitting beside each other -- two reconstructions at different
# allele fractions and TB-Profiler's own call. Overridable for that reason, and
# the plotting script's metric check is what confirms the choice: the wrong rule
# bit does not reproduce the reported holdout row.
RULE="${RULE:-$ROOT/rule_bit_${CODE}_tbprofiler.tsv}"
CFG="$SCRIPTS/final_config_${CODE}_${MODEL}.json"
CACHE="${CACHE:-$D/cache_${MODEL}_${TARGETSET}_HOLDOUT}"

for f in "$FIXED" "$RULE" "$CFG"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done

# Isoniazid predates the drug-parameterised scripts and its caches are named by
# hand, so the constructed path is a guess for that drug and a certainty for the
# rest. Say what exists rather than rebuilding a cache silently -- a rebuilt
# holdout cache would be a second PCA fit and the numbers would not be the
# reported ones.
if [ ! -e "$CACHE/precompute_manifest.json" ]; then
  echo "[ERROR] no holdout cache at $CACHE"
  echo "        caches present under $D:"
  ls -1d "$D"/cache_* 2>/dev/null | sed 's/^/          /' || echo "          (none)"
  echo "        pass the right one as CACHE=<path>. Do not let this job build"
  echo "        one: a fresh PCA fit would not reproduce the reported holdout."
  exit 1
fi

TOPK=$(grep -vc '^#' "$FIXED")
echo "[INFO] drug     = $DRUG ($CODE)   model = $MODEL"
echo "[INFO] cache    = $CACHE"
echo "[INFO] config   = $CFG"
echo "[INFO] targets  = $FIXED ($TOPK)   boundary = $SR_BOUNDARY"
echo "[INFO] rule     = $RULE"

mkdir -p "$OUTDIR"
DEST="$OUTDIR/${DRUG}_${MODEL}_holdout_predictions.tsv"

$PY "$SCRIPTS/cv_run_architecture.py" \
  --cache-dir "$CACHE" \
  --extra-feature-tsv "$RULE" \
  --mode shared4 --top-k "$TOPK" \
  --configs "$CFG" \
  --sr-boundary "$SR_BOUNDARY" \
  --label "${MODEL}_PREDICTIONS" \
  --out-dir "$D/runs/${MODEL}_PREDICTIONS" \
  --save-predictions "$DEST" \
  --n-jobs "$NJOBS" --seed 2026

echo
echo "[DONE] $DEST"
echo "       both aggregations are in the file; the plot picks the one that"
echo "       reproduces the reported metrics."
