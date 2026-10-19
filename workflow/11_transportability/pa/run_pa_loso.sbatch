#!/bin/bash
#SBATCH --job-name=pa_loso
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=0-9

# Leave-one-study-out for one P. aeruginosa drug: the §1 confound check that was
# declared and never run.
#
# Two prediction sets per drug, which is what makes this a measurement rather
# than a number:
#
#   LOSO    the reported config refit with one whole collection held out
#   RANDOM  the same config under the ordinary random CV, predictions saved
#
# Comparing them **on the same isolates** separates transport cost from study
# difficulty. A collection the model handles badly under both is simply hard; one
# it handles well under RANDOM and badly under LOSO was being predicted from its
# cohort. ARDaP needs no arm here -- it is training-free, so its per-study
# performance is already the difficulty anchor.
#
# Development isolates only. The holdout is not touched: this is a declared
# post-hoc check and must not spend the one-shot evaluation.
#
# The reported config is refit, never reselected. The question is whether the
# reported result transports, not whether a model re-tuned per study could.

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
PASCRIPTS="${PASCRIPTS:-$P/scripts}"
NJOBS="${NJOBS:-16}"

CVROOT="${CVROOT:-$P/imputed_strict}"
LABELS="${LABELS:-$P/labels_imputed_strict}"
CFGDIR="${CFGDIR:-$P/scripts/final_configs_strict}"
OUTDIR="${OUTDIR:-$P/loso_strict}"
# ARM=CR: the reported model for the eight antimicrobials ARDaP calls; the two it
# does not call exit cleanly. Outputs carry an _armCR suffix.
ARM="${ARM:-D}"
TAG="${TAG:-arm${ARM}1k}"
CACHE_TAG="${CACHE_TAG:-embseq}"
EMB_SUBDIR="${EMB_SUBDIR:-ml_inputs/embeddings/curated_plus_amrseq}"
METAROOT="${METAROOT:-$P/models/pseaer_all_antibiotics_amrfinder_groupedpca4_rf500_cv_package/data/metadata}"
MIN_ISOLATES="${MIN_ISOLATES:-40}"
MIN_PER_CLASS="${MIN_PER_CLASS:-8}"

DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)
declare -A SRB=(
  [amikacin]=4 [aztreonam]=4 [cefepime]=3 [ceftazidime]=3
  [ciprofloxacin]=-1 [imipenem]=2 [levofloxacin]=1 [meropenem]=3
  [piperacillin-tazobactam]=4 [tobramycin]=1
)
DRUG="${DRUGS[$SLURM_ARRAY_TASK_ID]}"
SR_BOUNDARY="${SRB[$DRUG]}"

D="$CVROOT/cv_${DRUG}"
FIXED="$P/targets/combined/${DRUG}_all.txt"
ZEROFILL="$P/targets/amrseq_targets.txt"
EXTRA="$LABELS/amrfinder/by_drug_tsv/${DRUG}.tsv"
SFX=""
if [ "$ARM" = CR ]; then
  case "$DRUG" in aztreonam|levofloxacin) echo "[SKIP] $DRUG: no ARDaP call"; exit 0 ;; esac
  EXTRA="$LABELS/ardap_features/${DRUG}_armCR.tsv"; SFX=_armCR
fi
META="$METAROOT/$DRUG/metadata.tsv"

for f in "$D/cv_folds_final.tsv" "$D/cv_folds.tsv" "$FIXED" "$ZEROFILL" \
         "$EXTRA" "$META" "$PASCRIPTS/build_loso_folds.py" \
         "$PASCRIPTS/pick_reported_model.py"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done

# The reported model, by the same script the figures and the SHAP jobs use, so
# all four answer the question identically. --config-dir is mandatory against a
# full run directory: folds.tsv there holds all 1,000 configs and averaging over
# the grid would pick a different model than the panel reports.
MODEL=$($PY "$PASCRIPTS/pick_reported_model.py" --run-dir "$D" --tag "$TAG" \
          --config-dir "$CFGDIR" --drug "$DRUG")
CFG="$CFGDIR/final_config_${DRUG}_${MODEL}_arm${ARM}.json"
[ -e "$CFG" ] || { echo "[ERROR] missing $CFG"; exit 1; }
TOPK=$(grep -vc '^#' "$FIXED")

echo "[INFO] drug $DRUG   model $MODEL   boundary $SR_BOUNDARY   top-k $TOPK"
echo "[INFO] config $CFG"

mkdir -p "$OUTDIR"

# ---- 1. leave-one-study-out folds ------------------------------------------
echo
echo "################ folds"
$PY "$PASCRIPTS/build_loso_folds.py" \
  --metadata "$META" \
  --folds-final "$D/cv_folds_final.tsv" \
  --drug "$DRUG" \
  --boundary "$SR_BOUNDARY" \
  --min-isolates "$MIN_ISOLATES" \
  --min-per-class "$MIN_PER_CLASS" \
  --out "$D/cv_folds_loso.tsv" \
  --summary "$OUTDIR/${DRUG}_studies.tsv"

n_fold=$($PY - "$D/cv_folds_loso.tsv" <<'PYEOF'
import csv, sys
f = {int(r["fold"]) for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
print(len([x for x in f if x >= 0]))
PYEOF
)
if [ "$n_fold" -eq 0 ]; then
  echo "[DONE] $DRUG: no scoreable study — a single-collection drug. Nothing to run."
  exit 0
fi
echo "[INFO] $n_fold scoreable study fold(s)"

declare -A MROOT=( [evo2]="$P/results/ml/pa_evo2"
                   [carbon]="$P/results/ml/pa_carbon"
                   [dnabert]="$P/results/ml/pa_dnabert" )

# ---- 2. features, refit inside every study fold ----------------------------
# Non-negotiable: a projection fitted across all studies has already seen the
# held-out collection, and the whole point of the check is that it has not.
echo
echo "################ precompute (PCA refit per study fold)"
CACHE="$D/cache_${MODEL}_${CACHE_TAG}_LOSO"
if [ ! -e "$CACHE/precompute_manifest.json" ]; then
  $PY "$SCRIPTS/cv_precompute_fold_features.py" \
    --root "${MROOT[$MODEL]}" \
    --embedding-subdir "$EMB_SUBDIR" \
    --folds-tsv "$D/cv_folds_loso.tsv" \
    --pca-mode shared4 --pca-dim 4 \
    --var-threshold 0.99 --max-components 25 \
    --top-targets 256 --seed 2026 \
    --fixed-targets "$FIXED" \
    --zero-fill-targets "$ZEROFILL" \
    --out-dir "$CACHE"
else
  echo "[INFO] reusing $CACHE"
fi

# ---- 3. LOSO predictions ---------------------------------------------------
echo
echo "################ LOSO"
$PY "$SCRIPTS/cv_run_architecture.py" \
  --cache-dir "$CACHE" \
  --extra-feature-tsv "$EXTRA" \
  --mode shared4 --top-k "$TOPK" \
  --configs "$CFG" \
  --sr-boundary "$SR_BOUNDARY" \
  --label "${MODEL}${SFX}_LOSO" \
  --save-predictions "$OUTDIR/${DRUG}_${MODEL}${SFX}_loso_predictions.tsv" \
  --out-dir "$D/runs/${MODEL}${SFX}_LOSO" \
  --n-jobs "$NJOBS" --seed 2026

# ---- 4. the matched random-CV arm ------------------------------------------
# Same config, same isolates, ordinary folds. Without this the LOSO numbers can
# only be read against a different cohort, and a study that is simply hard would
# be indistinguishable from one the model was predicting by its cohort.
echo
echo "################ matched random CV"
DEVCACHE="$D/cache_${MODEL}_${CACHE_TAG}"
[ -e "$DEVCACHE/precompute_manifest.json" ] || {
  echo "[ERROR] missing $DEVCACHE — the development cache should already exist"; exit 1; }
$PY "$SCRIPTS/cv_run_architecture.py" \
  --cache-dir "$DEVCACHE" \
  --extra-feature-tsv "$EXTRA" \
  --mode shared4 --top-k "$TOPK" \
  --configs "$CFG" \
  --sr-boundary "$SR_BOUNDARY" \
  --label "${MODEL}${SFX}_CVPRED" \
  --save-predictions "$OUTDIR/${DRUG}_${MODEL}${SFX}_cv_predictions.tsv" \
  --out-dir "$D/runs/${MODEL}${SFX}_CVPRED" \
  --n-jobs "$NJOBS" --seed 2026

echo
echo "[DONE] $DRUG ($MODEL${SFX}): $OUTDIR/${DRUG}_${MODEL}${SFX}_{loso,cv}_predictions.tsv"
