#!/bin/bash
#SBATCH --job-name=hold_cat
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-9

# Held-out evaluation of mutation columns + TB-Profiler call, no embeddings (§27).
#
# The paper compares embeddings against the catalogue's own mutations on the MIC
# ladder. Until now that comparison existed on development folds only. This scores
# the mutation-column model once on each drug's held-out set, so the comparison
# rests on the same isolates as every other held-out result.
#
#   matrix         rebuilt over ALL isolates of the final fold file -- the
#                  development matrix has no held-out rows. The column set is fixed
#                  per drug (no carrier floor), so it is the same model; a column
#                  constant among development isolates is dropped at fit, as before.
#   configuration  the one selected on development folds
#                  (slim_arms/configs/<drug>__catalogue__snprule1k.json)
#   partition      cv_folds_final.tsv: train on development, score the held-out set
#
# SCORED ONCE: stops if results already exist.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
MATDIR="${MATDIR:-$MTB/results/tbprofiler/reads/gather/matrix}"
NJOBS="${NJOBS:-8}"
TAG="${TAG:-snprule1k}"

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
RULE="$ROOT/rule_bit_${CODE}_tbprofiler.tsv"
CFG="$ROOT/slim_arms/configs/${FOLDDIR#cv_}__catalogue__${TAG}.json"
DEVMAT="$ROOT/catalogue_combined/catalogue_${CODE}_${TAG}.tsv"
MATRIX="$ROOT/catalogue_combined/catalogue_${CODE}_${TAG}_allisolates.tsv"
RUN="$D/runs/catalogue_${TAG}_HOLDOUT"
PREDS="$ROOT/predictions/${DRUG}_catalogue_holdout_predictions_${TAG}.tsv"

for f in "$FOLDS_FINAL" "$RULE" "$CFG" "$DEVMAT" \
         "$MATDIR/tbprofiler_resistance_mutation_binary.tsv" \
         "$MATDIR/tbprofiler_amrfinder_style_long.tsv"; do
  [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; exit 1; }
done
if [ -e "$RUN/folds.tsv" ] || [ -e "$PREDS" ]; then
  echo "[ERROR] $DRUG: held-out results already exist ($RUN or $PREDS)."
  echo "        The held-out set is scored once. Not rescoring."
  exit 1
fi
mkdir -p "$ROOT/predictions"

$PY "$SCRIPTS/merge_catalogue_features.py" \
  --rule-tsv "$RULE" \
  --mutation-matrix "$MATDIR/tbprofiler_resistance_mutation_binary.tsv" \
  --long "$MATDIR/tbprofiler_amrfinder_style_long.tsv" \
  --drug "$DRUG" \
  --samples "$FOLDS_FINAL" \
  --out "$MATRIX"

# Same columns, in the same order, as the matrix the configuration was selected on.
[ "$(head -1 "$MATRIX")" = "$(head -1 "$DEVMAT")" ] || {
  echo "[ERROR] $DRUG: all-isolate matrix has different columns from the development matrix"
  exit 1; }
echo "[INFO] drug=$DRUG  config=$(basename "$CFG")  tag=$($PY -c \
  "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$CFG")"
echo "[INFO] matrix columns identical to development: $(( $(head -1 "$MATRIX" | tr '\t' '\n' | wc -l) - 1 ))"

$PY "$SCRIPTS/cv_run_architecture.py" \
  --feature-tsv "$MATRIX" \
  --folds-tsv "$FOLDS_FINAL" \
  --configs "$CFG" \
  --sr-boundary "$BOUND" \
  --label "catalogue_${TAG}_HOLDOUT" \
  --save-predictions "$PREDS" \
  --out-dir "$RUN" \
  --n-jobs "$NJOBS" --seed 2026

echo "[DONE] $DRUG -> $RUN ; predictions $PREDS"
