#!/bin/bash
#SBATCH --job-name=hold_arm
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-9

# Held-out evaluation of a feature arm other than the originally reported one.
#
# Written for the wider target set (gwasfiltered + rule bit, §27). It is separate
# from run_panel_holdout.sbatch on purpose: that script takes its configuration
# from final_config_<code>_<model>.json and writes runs/<model>_HOLDOUT, so reusing
# it would score the wrong configuration and overwrite the reported results.
#
# Per drug and model, in one pass:
#   configuration  the one 1-SE selection chose on development folds for this arm
#                  (slim_arms/configs/<drug>__<model>__<targetset>_<tag>.json)
#   top-k          recovered from that development run's feature count, and the
#                  refit's feature count is checked against it
#   projection     fitted on the development isolates only (--only-fold 0:0)
#   output         runs/<model>_<targetset>_<tag>_HOLDOUT, with per-isolate
#                  predictions saved in the same fit
#
# SCORED ONCE. If the run directory already holds results the task stops; there is
# no FORCE. A second score is a decision to be made by deleting the first by hand.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
EMB_SUBDIR="${EMB_SUBDIR:-ml_inputs/embeddings/panel_af010}"
NJOBS="${NJOBS:-16}"
TARGETSET="${TARGETSET:-gwasfiltered}"
TAG="${TAG:-gwfrule1k}"
MODELS="${MODELS:-evo2 carbon dnabert}"

# drug  code  boundary  fold-dir
SPEC=(
  "isoniazid    inh -3.0 cv_isoniazid_final"
  "rifampicin   rif -1.0 cv_rifampicin"
  "ethambutol   emb  2.0 cv_ethambutol"
  "ethionamide  eth  2.0 cv_ethionamide"
  "kanamycin    kan  2.0 cv_kanamycin"
  "amikacin     ami  0.0 cv_amikacin"
  "levofloxacin lev  0.0 cv_levofloxacin"
  "moxifloxacin mxf  0.0 cv_moxifloxacin"
  "bedaquiline  bdq -2.0 cv_bedaquiline"
  "delamanid    dlm -3.0 cv_delamanid"
)
read -r DRUG CODE BOUND FOLDDIR <<< "${SPEC[$SLURM_ARRAY_TASK_ID]}"
D="$ROOT/$FOLDDIR"
FOLDS_FINAL="$D/cv_folds_final.tsv"
FIXED="$ROOT/panel_targets/targets_${DRUG}_${TARGETSET}.txt"
RULE="$ROOT/rule_bit_${CODE}_tbprofiler.tsv"
PREDDIR="$ROOT/predictions"
for f in "$FOLDS_FINAL" "$FIXED" "$RULE"; do
  [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; exit 1; }
done
n_fixed=$(grep -vc '^#' "$FIXED")
N_EXTRA=$(( $(head -1 "$RULE" | awk -F'\t' '{print NF}') - 1 ))
mkdir -p "$PREDDIR"
echo "[INFO] drug=$DRUG ($CODE)  arm=${TARGETSET}_${TAG}  boundary=$BOUND"
echo "[INFO] folds=$FOLDS_FINAL  targets=$FIXED ($n_fixed)  rule=$RULE ($N_EXTRA column)"

# A case statement, not an associative array, so the script also runs under
# bash 3 (macOS), where it is tested against stand-ins before any real run.
model_root() {
  case "$1" in
    evo2)    echo "$MTB/results/ml/03_mtb_evo2" ;;
    carbon)  echo "$MTB/results/ml/03_mtb" ;;
    dnabert) echo "$MTB/results/ml/03_mtb_dnabert" ;;
    *) echo "[ERROR] unknown model $1" >&2; return 1 ;;
  esac
}

for M in $MODELS; do
  echo
  echo "################ $M ################"
  CFG="$ROOT/slim_arms/configs/${FOLDDIR#cv_}__${M}__${TARGETSET}_${TAG}.json"
  DEVF="$ROOT/slim_arms/$FOLDDIR/runs/${M}_${TARGETSET}_${TAG}/folds.tsv"
  DEV="$D/cache_${M}_${TARGETSET}"
  OUT="$D/cache_${M}_${TARGETSET}_HOLDOUT"
  RUN="$D/runs/${M}_${TARGETSET}_${TAG}_HOLDOUT"
  PREDS="$PREDDIR/${DRUG}_${M}_holdout_predictions_${TAG}.tsv"
  for f in "$CFG" "$DEVF" "$DEV/precompute_manifest.json"; do
    [ -e "$f" ] || { echo "[ERROR] $DRUG/$M: missing $f"; exit 1; }
  done
  if [ -e "$RUN/folds.tsv" ] || [ -e "$PREDS" ]; then
    echo "[ERROR] $DRUG/$M: held-out results already exist ($RUN or $PREDS)."
    echo "        The held-out set is scored once. Not rescoring."
    exit 1
  fi

  DEV_NF=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="n_features")c=i;next}
                       NR==2{print $c; exit}' "$DEVF")
  [ $(( (DEV_NF - N_EXTRA) % 4 )) -eq 0 ] || {
    echo "[ERROR] development n_features $DEV_NF minus $N_EXTRA is not a multiple of 4"; exit 1; }
  TOPK=$(( (DEV_NF - N_EXTRA) / 4 ))
  [ "$TOPK" -eq "$n_fixed" ] || {
    echo "[ERROR] development used $TOPK targets but $FIXED holds $n_fixed"; exit 1; }
  echo "[INFO] config $(basename "$CFG")  tag=$($PY -c \
    "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$CFG")  top-k=$TOPK"

  if [ ! -e "$OUT/precompute_manifest.json" ]; then
    $PY "$SCRIPTS/cv_precompute_fold_features.py" \
      --root "$(model_root "$M")" --embedding-subdir "$EMB_SUBDIR" \
      --folds-tsv "$FOLDS_FINAL" \
      --pca-mode shared4 --pca-dim 4 \
      --var-threshold 0.99 --max-components 25 \
      --top-targets 256 --seed 2026 \
      --fixed-targets "$FIXED" \
      --only-fold 0:0 \
      --out-dir "$OUT"
  else
    echo "[INFO] reusing $OUT"
  fi
  $PY - "$DEV/precompute_manifest.json" "$OUT/precompute_manifest.json" <<'PYEOF'
import json, sys
dev, hold = (json.load(open(p)) for p in sys.argv[1:3])
bad = [k for k in ("n_targets_total", "n_fixed_targets", "pca_mode", "pca_dim",
                   "var_threshold", "max_components", "top_targets", "seed")
       if dev.get(k) != hold.get(k)]
if bad:
    sys.exit("[ERROR] held-out cache differs from development on "
             + ", ".join(f"{k} {dev.get(k)!r} vs {hold.get(k)!r}" for k in bad))
print(f"[OK] geometry matches development ({hold['n_fixed_targets']} targets, "
      f"{hold['n_samples']:,} isolates)")
PYEOF

  $PY "$SCRIPTS/cv_run_architecture.py" \
    --cache-dir "$OUT" \
    --extra-feature-tsv "$RULE" \
    --mode shared4 --top-k "$TOPK" \
    --configs "$CFG" \
    --sr-boundary "$BOUND" \
    --label "${M}_${TARGETSET}_${TAG}_HOLDOUT" \
    --save-predictions "$PREDS" \
    --out-dir "$RUN" \
    --n-jobs "$NJOBS" --seed 2026

  NEW_NF=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="n_features")c=i;next}
                       NR==2{print $c; exit}' "$RUN/folds.tsv")
  [ "$NEW_NF" = "$DEV_NF" ] || {
    echo "[ERROR] $DRUG/$M: held-out fit has $NEW_NF features, development had $DEV_NF"; exit 1; }
  echo "[OK] $DRUG/$M: n_features $NEW_NF matches development; predictions $PREDS"
done
echo
echo "[DONE] $DRUG -> $D/runs/{${MODELS// /,}}_${TARGETSET}_${TAG}_HOLDOUT"
