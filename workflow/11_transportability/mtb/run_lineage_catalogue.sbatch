#!/bin/bash
#SBATCH --job-name=lolo_cat
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-9

# §21: does the embeddings' advantage on the MIC ladder survive leaving a lineage
# out? Refits mutation columns + TB-Profiler call (no embeddings) on the lineage
# folds and size-matched folds that run_panel_lineage.sbatch already built, so it
# can be set against the embeddings + call refits that job made on the same folds.
#
# Whole-gene embeddings carry each isolate's lineage background; catalogue
# mutation columns carry almost none. An embedding advantage on the ladder that is
# partly lineage recognition would shrink here. That is the test.
#
# No selection: the configuration is the one 1-SE selection chose on development
# folds (slim_arms/configs/<drug>__catalogue__snprule1k.json), and the matrix is
# the one the development run used. Needs run_panel_lineage.sbatch to have built
# the drug's folds, and slim_mtb_runs.sh to have collected the development run.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
NJOBS="${NJOBS:-8}"
TAG="${TAG:-snprule1k}"

# drug  code  boundary  fold-dir   (same table as run_panel_lineage.sbatch)
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
L="$ROOT/$FOLDDIR/lineage"

CFG="$ROOT/slim_arms/configs/${FOLDDIR#cv_}__catalogue__${TAG}.json"
MATRIX="$ROOT/catalogue_combined/catalogue_${CODE}_${TAG}.tsv"
DEVF="$ROOT/slim_arms/$FOLDDIR/runs/catalogue_${TAG}/folds.tsv"
for f in "$CFG" "$MATRIX" "$DEVF" "$L/cv_folds.tsv" "$L/sizematched/cv_folds.tsv"; do
  [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; exit 1; }
done
echo "[INFO] drug=$DRUG  config=$(basename "$CFG")  tag=$($PY -c \
  "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$CFG")"
echo "[INFO] matrix $(basename "$MATRIX"): $(( $(head -1 "$MATRIX" | tr '\t' '\n' | wc -l) - 1 )) columns"

for split in lineage sizematched; do
  if [ "$split" = lineage ]; then FD="$L"; NAME=catalogue_lolo_snprule
  else FD="$L/sizematched"; NAME=catalogue_sm_snprule; fi
  echo "=== $split"
  $PY "$SCRIPTS/cv_run_architecture.py" \
    --feature-tsv "$MATRIX" \
    --folds-tsv "$FD/cv_folds.tsv" \
    --configs "$CFG" \
    --sr-boundary "$BOUND" \
    --label "$NAME" \
    --out-dir "$FD/runs/$NAME" \
    --n-jobs "$NJOBS" --seed 2026
done
echo "[DONE] $DRUG -> $L/runs/catalogue_lolo_snprule, $L/sizematched/runs/catalogue_sm_snprule"
