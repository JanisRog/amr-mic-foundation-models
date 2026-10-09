#!/bin/bash
#SBATCH --job-name=pa_dev
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# One P. aeruginosa drug's development CV: three models, 1,000 configs, 3x5 folds.
#
# Embeddings only -- no rule bit. AMRFinder emits gene hits with drug classes,
# not a per-drug S/R call, so there is no equivalent of TB-Profiler's own call to
# feed in. Reconstructing one would repeat the ethionamide mistake (section 11:
# substituting a reconstruction for the tool's own call cost 0.077 before anyone
# noticed), so the catalogue feature is deliberately absent until the acquired-gene
# half is designed. That makes this a pure-sequence arm, and the comparison to
# AMRFinder is correspondingly stronger if it wins and weaker if it loses.
#
# Folds are stratified on MIC class only, as for MTB. These 2,127 isolates come
# from eight studies whose resistance rates differ substantially, so study is a
# plausible confound -- it is treated as a post-hoc check (leave-one-study-out),
# not as a fold constraint.
#
#   DRUG   ciprofloxacin      CODE  cip
#   TOPK   38                 (must equal the curated target-list length)
#   SR_BOUNDARY -1.0          (no default: the parameter that fails silently)

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
EMB_SUBDIR="${EMB_SUBDIR:-ml_inputs/embeddings/curated_af010}"
NJOBS="${NJOBS:-16}"
TAG="${TAG:-final1k}"
# Optional extra feature block (acquired-gene presence). Empty = embeddings only,
# which is the arm every other arm is measured against.
EXTRA="${EXTRA:-}"

: "${DRUG:?set DRUG, e.g. ciprofloxacin}"
: "${TOPK:?set TOPK to the number of targets in the curated list}"
: "${SR_BOUNDARY:?set SR_BOUNDARY explicitly -- it has no safe default}"

# CVROOT lets a re-labelled panel live beside the original instead of
# overwriting it. The imputed-censoring rebuild uses CVROOT=$P/imputed, so the
# published cohorts stay intact and reproducible until the new ones are shown to
# be better -- and cefepime and levofloxacin, whose labels are bit-identical
# under both rules, become a regression test rather than a rerun.
D="${CVROOT:-$P}/cv_${DRUG}"
# FIXED defaults to the curated chromosomal set; arm C points it at the combined
# list. ZEROFILL names the targets where absence is meaningful rather than fatal.
FIXED="${FIXED:-$P/targets/targets_${DRUG}_curated.txt}"
ZEROFILL="${ZEROFILL:-}"
CONFIGS="${CONFIGS:-$MTB/scripts/mtb_rf1000_configs.json}"

for f in "$D/cv_folds.tsv" "$FIXED" "$CONFIGS"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done
n_fixed=$(grep -vc '^#' "$FIXED")
[ "$n_fixed" -eq "$TOPK" ] || {
  echo "[ERROR] TOPK=$TOPK but $FIXED holds $n_fixed targets"; exit 1; }

echo "[INFO] drug     = $DRUG   boundary $SR_BOUNDARY"
echo "[INFO] targets  = $FIXED ($n_fixed)"
[ -n "$ZEROFILL" ] && echo "[INFO] zero-fill = $ZEROFILL ($(grep -vc '^#' "$ZEROFILL") targets may be absent)"
echo "[INFO] embeddings subdir = $EMB_SUBDIR"
if [ -n "$EXTRA" ]; then
  [ -e "$EXTRA" ] || { echo "[ERROR] missing extra features: $EXTRA"; exit 1; }
  echo "[INFO] extra features = $EXTRA ($(( $(head -1 "$EXTRA" | tr '\t' '\n' | wc -l) - 1 )) columns)"
else
  echo "[INFO] extra features = NONE (embeddings only)"
fi

declare -A MROOT=( [evo2]="$P/results/ml/pa_evo2"
                   [carbon]="$P/results/ml/pa_carbon"
                   [dnabert]="$P/results/ml/pa_dnabert" )

for M in evo2 carbon dnabert; do
  echo
  echo "################ $M"
  # The cache name must vary with BOTH the model and the feature set. An earlier
  # version took a CACHE_PREFIX from the environment, which cannot contain $M --
  # the submitting shell has no such variable -- so all three models shared one
  # directory, and the second and third took the "reusing" branch and were scored
  # on the first model's embeddings. Three identical numbers that look like a
  # model comparison. Derived here instead, where $M exists.
  # Defaults to `curated` so arms that differ only in their extra-feature TSV
  # share one precompute -- embeddings-only and embeddings+presence have
  # identical cached features, the presence columns being appended at scoring
  # time. An arm with a different target set passes CACHE_TAG explicitly.
  CACHE="$D/cache_${M}_${CACHE_TAG:-curated}"
  if [ ! -e "$CACHE/precompute_manifest.json" ]; then
    $PY "$SCRIPTS/cv_precompute_fold_features.py" \
      --root "${MROOT[$M]}" \
      --embedding-subdir "$EMB_SUBDIR" \
      --folds-tsv "$D/cv_folds.tsv" \
      --pca-mode shared4 --pca-dim 4 \
      --var-threshold 0.99 --max-components 25 \
      --top-targets 256 --seed 2026 \
      --fixed-targets "$FIXED" \
      ${ZEROFILL:+--zero-fill-targets "$ZEROFILL"} \
      --out-dir "$CACHE"
  else
    echo "[INFO] reusing $CACHE"
  fi

  $PY "$SCRIPTS/cv_run_architecture.py" \
    --cache-dir "$CACHE" \
    ${EXTRA:+--extra-feature-tsv "$EXTRA"} \
    --mode shared4 --top-k "$TOPK" \
    --configs "$CONFIGS" \
    --sr-boundary "$SR_BOUNDARY" \
    --label "${M}_${TAG}" \
    --out-dir "$D/runs/${M}_${TAG}" \
    --n-jobs "$NJOBS" --seed 2026
done

echo
echo "[DONE] $D/runs/{evo2,carbon,dnabert}_${TAG}"
