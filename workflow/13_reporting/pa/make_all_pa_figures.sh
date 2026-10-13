#!/usr/bin/env bash
# Every P. aeruginosa figure, from the downloaded cluster outputs, in one command.
#
# The MTB twin of this script, with four differences that matter:
#
#   * PA and MTB share drug names -- amikacin and levofloxacin are in both panels
#     -- so PA inputs live in their own directory. Downloading PA predictions
#     into the MTB directory silently overwrites two MTB drugs, and every MTB
#     figure then fails its workbook check with no obvious cause.
#   * The reported arm is D (embeddings + acquired sequence + presence), declared
#     as the complete nested arm rather than as the development winner: arm C has
#     the lower development cost on 6 of 10 drugs.
#   * Breakpoints are the PA protocol's derived table, recovered from
#     resistance_label_TMP rather than supplied from memory.
#   * The reported model is Evo 2 for every drug, the primary model designated
#     for M. tuberculosis (until 2026-10-13 it was chosen per drug on development
#     cost, which gave carbon for amikacin and meropenem). It is asked of
#     pick_reported_model.py rather than hardcoded here, so the figures, the SHAP
#     jobs and the study-held-out jobs cannot disagree.

set -euo pipefail

IN="${IN:-$HOME/Downloads/pa_strict}"
PRED="${PRED:-$IN/predictions_strict}"
SLIM="${SLIM:-$IN/slim_strict}"
WB="${WB:-$HOME/Desktop/AMRPredict_pseaer_strict_results.xlsx}"
OUT="${OUT:-figures/pseaer_strict}"
TAG="${TAG:-armD1k}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MTB="${MTB_SCRIPTS:-$HERE/../mtb_tbprofiler}"          # the drawing scripts are organism-agnostic
PICK="$HERE/pick_reported_model.py"

DRUGS=(
  "ciprofloxacin -1" "levofloxacin 1" "tobramycin 1" "imipenem 2" "meropenem 3"
  "ceftazidime 3"    "cefepime 3"     "aztreonam 4"  "amikacin 4"
  "piperacillin-tazobactam 4"
)

# ---- which model each drug reports, by mean development cost ---------------
# No associative arrays: macOS ships bash 3.2 and this script runs on the laptop.
# MAP is the single source of truth; model_of() reads it back.
[ -e "$PICK" ] || { echo "[ERROR] missing $PICK"; exit 1; }
MAP=""
for spec in "${DRUGS[@]}"; do
  set -- $spec
  m=$(python3 "$PICK" --run-dir "$SLIM/cv_$1" --tag "$TAG")
  [ -n "$m" ] || { echo "[ERROR] could not pick a model for $1"; exit 1; }
  MAP="${MAP:+$MAP,}$1=$m"
done
model_of() { printf '%s' "$MAP" | tr ',' '\n' | grep "^$1=" | cut -d= -f2; }
echo "[INFO] reported model per drug: $MAP"

for spec in "${DRUGS[@]}"; do
  set -- $spec
  M=$(model_of "$1")
  for f in "$PRED/${1}_${M}_holdout_predictions.tsv" "$SLIM/cv_${1}/cv_folds_final.tsv"; do
    [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
  done
done
[ -e "$WB" ] || { echo "[ERROR] missing workbook $WB"; exit 1; }
echo "[INFO] inputs present for ${#DRUGS[@]} drugs"

# ---- per-drug figures ------------------------------------------------------
for spec in "${DRUGS[@]}"; do
  set -- $spec
  M=$(model_of "$1")
  python3 "$MTB/plot_mic_scatter.py" --predictions "$PRED/${1}_${M}_holdout_predictions.tsv" \
    --workbook "$WB" --drug "$1" --model "$M" --sr-boundary "$2" \
    --out-dir "$OUT" --stem "scatter_$1"
  python3 "$MTB/plot_class_distribution.py" --predictions "$PRED/${1}_${M}_holdout_predictions.tsv" \
    --folds "$SLIM/cv_${1}/cv_folds_final.tsv" --workbook "$WB" --drug "$1" --model "$M" \
    --sr-boundary "$2" --out-dir "$OUT" --stem "classdist_$1"
done

# ---- per-drug SHAP beeswarms ----------------------------------------------
# Only for drugs whose SHAP run exists: SHAP comes from a separate cluster job,
# so a panel is routinely half-attributed and a missing beeswarm is a state, not
# an error. The file is named for the reported model, so a drug whose model
# changed picks up its new run rather than silently drawing the old one.
n_shap=0
for spec in "${DRUGS[@]}"; do
  set -- $spec
  M=$(model_of "$1")
  L="$PRED/${1}_${M}_shap_long.tsv"
  [ -e "$L" ] || { echo "[skip] no SHAP for $1/$M"; continue; }
  python3 "$MTB/plot_shap_beeswarm.py" --shap-long "$L" \
    --predictions "$PRED/${1}_${M}_holdout_predictions.tsv" \
    --workbook "$WB" --drug "$1" --model "$M" --sr-boundary "$2" \
    --out-dir "$OUT" --stem "shap_$1"
  n_shap=$((n_shap + 1))
done
echo "[INFO] $n_shap of ${#DRUGS[@]} drugs have SHAP beeswarms"

# ---- F1 and one-tier accuracy as grouped bars ------------------------------
python3 "$MTB/plot_metric_bars.py" --workbook "$WB" --panel pseaer --out-dir "$OUT"

# ---- the three models side by side ----------------------------------------
python3 "$MTB/plot_model_comparison.py" --workbook "$WB" --panel pseaer \
  --label "P. aeruginosa" --stem pa_model_comparison --out-dir "$OUT"

# ---- cross-panel: F1 against each organism's own baseline ------------------
# Needs both organisms, so it lives here rather than in either panel's folder and
# is skipped when the MTB side is not present. Generated from the PA driver
# because the PA numbers are the ones that move.
MTB_WB="${MTB_WB:-$HOME/Desktop/AMRPredict_panel_results.xlsx}"
CROSS_OUT="${CROSS_OUT:-figures/cross}"
if [ -e "$MTB_WB" ] && [ -d "$IN/cv_amikacin" ]; then
  python3 "$MTB/plot_f1_boxes.py" --mtb-workbook "$MTB_WB" \
    --pa-workbook "$WB" --pa-baseline-dir "$IN" --out-dir "$CROSS_OUT"
  python3 "$MTB/plot_metric_bars.py" --panel pooled --workbook "$MTB_WB" \
    --pa-workbook "$WB" --out-dir "$CROSS_OUT"
else
  echo "[skip] cross-panel F1 figure: need $MTB_WB and $IN/cv_<drug>/baseline_rows.tsv"
fi

# ---- combined sheets -------------------------------------------------------
# One figure each for the scatters, the class distributions and the beeswarms.
# --model-map carries the same per-drug choice, so a sheet cannot plot a
# different model than the single-drug figure beside it in the paper.
python3 "$MTB/make_panel_grids.py" \
  --pred-dir "$PRED" --fold-dir "$SLIM" --panel pseaer \
  --model-map "$MAP" --workbook "$WB" \
  --stem-prefix pa_ --out-dir "$OUT"

echo
echo "[DONE] $OUT — $(ls "$OUT"/*.pdf 2>/dev/null | wc -l | tr -d ' ') PDFs"
