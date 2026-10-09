#!/usr/bin/env bash
# Every MTB figure, from the downloaded cluster outputs, in one command.
#
# This exists instead of a command list in the protocol. A list in prose is
# checked by nobody and drifts the first time a flag changes; this is the thing
# that actually runs, so the protocol points here rather than restating it.
#
# What it needs, and which cluster job produces it:
#
#   <drug>_evo2_holdout_predictions.tsv     run_panel_predictions.sbatch
#   isoniazid_evo2_holdout_predictions_tbrule.tsv
#                                           run_inh_redo.sbatch MODE=holdout
#   cv_<drug>/cv_folds_final.tsv            already on the cluster; tar them over
#   AMRPredict_panel_results.xlsx           build_results_workbook.py
#
# Every figure re-derives its metrics from the predictions and refuses to draw
# unless they reproduce the workbook row, so a stale download fails loudly rather
# than producing a figure that disagrees with the table beside it.

set -euo pipefail

IN="${IN:-$HOME/Downloads}"
WB="${WB:-$HOME/Desktop/AMRPredict_panel_results.xlsx}"
OUT="${OUT:-figures/mtb}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# drug + S/R boundary. The boundaries are section 16's table, itself recovered by
# solving rather than by recall; regenerate with recover_sr_boundaries.py.
DRUGS=(
  "isoniazid -3"  "rifampicin -1"  "ethambutol 2"  "ethionamide 2"
  "kanamycin 2"   "amikacin 0"     "levofloxacin 0" "moxifloxacin 0"
  "bedaquiline -2" "delamanid -3"
)

preds() {  # isoniazid's redo is named apart from the panel's
  if [ "$1" = isoniazid ]; then
    echo "$IN/isoniazid_evo2_holdout_predictions_tbrule.tsv"
  else
    echo "$IN/${1}_evo2_holdout_predictions.tsv"
  fi
}
folds() {  # and its fold directory is cv_isoniazid_final, not cv_isoniazid
  if [ "$1" = isoniazid ]; then
    echo "$IN/cv_isoniazid_final/cv_folds_final.tsv"
  else
    echo "$IN/cv_${1}/cv_folds_final.tsv"
  fi
}

for spec in "${DRUGS[@]}"; do
  set -- $spec
  for f in "$(preds "$1")" "$(folds "$1")"; do
    [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
  done
done
[ -e "$WB" ] || { echo "[ERROR] missing workbook $WB"; exit 1; }
echo "[INFO] inputs present for ${#DRUGS[@]} drugs"

echo "=== panel figures (workbook only)"
python3 "$HERE/make_panel_figures.py" --workbook "$WB" --out-dir "$OUT"

echo "=== per-drug confusion grids and class distributions"
for spec in "${DRUGS[@]}"; do
  set -- $spec
  python3 "$HERE/plot_mic_scatter.py" --predictions "$(preds "$1")" \
    --workbook "$WB" --drug "$1" --model evo2 --sr-boundary "$2" \
    --out-dir "$OUT" --stem "scatter_$1"
  python3 "$HERE/plot_class_distribution.py" --predictions "$(preds "$1")" \
    --folds "$(folds "$1")" --workbook "$WB" --drug "$1" --model evo2 \
    --sr-boundary "$2" --out-dir "$OUT" --stem "classdist_$1"
done

echo "=== SHAP beeswarms (only for drugs whose SHAP job has been run)"
for spec in "${DRUGS[@]}"; do
  set -- $spec
  [ -e "$IN/${1}_evo2_shap_long.tsv" ] || { echo "[skip] $1: no SHAP output"; continue; }
  python3 "$HERE/plot_shap_beeswarm.py" --shap-long "$IN/${1}_evo2_shap_long.tsv" \
    --predictions "$IN/${1}_evo2_predictions.tsv" --workbook "$WB" \
    --drug "$1" --model evo2 --sr-boundary "$2" \
    --out-dir "$OUT" --stem "shap_$1"
done

echo "=== model comparison"
python3 "$HERE/plot_model_comparison.py" --workbook "$WB" --panel mtb \
  --label "M. tuberculosis" --stem model_comparison --out-dir "$OUT"

echo "=== combined sheets"
python3 "$HERE/make_panel_grids.py" --pred-dir "$IN" --fold-dir "$IN" \
  --workbook "$WB" --out-dir "$OUT"

echo "=== review page"
python3 "$HERE/build_figure_page.py" --workbook "$WB" --fig-dir "$OUT" \
  --out "$OUT/index.html"

echo
echo "[DONE] $OUT — $(ls "$OUT"/*.pdf | wc -l | tr -d ' ') PDFs, and index.html"
