#!/usr/bin/env bash
# The P. aeruginosa results workbook, from downloaded cluster outputs.
#
# Never scripted before -- the previous panel's workbook was assembled by hand,
# one --model line at a time, which is how a comparison table stops being built
# and starts being retyped. Every figure script checks its numbers against this
# file, so a hand-assembled workbook makes the checks circular.
#
# Two passes per drug, deliberately:
#   1. the comparator alone, so ARDaP's rows exist before any model row does
#   2. the three models appended
# --append replaces a drug's rows rather than duplicating them, so re-running is
# safe. ARDaP covers 8 of 10 drugs; the two it does not are model-only.

set -euo pipefail

IN="${IN:-$HOME/Downloads/pa_strict}"
PRED="${PRED:-$IN/predictions_strict}"
SLIM="${SLIM:-$IN/slim_strict}"
BASE="${BASE:-$IN}"                       # cv_<drug>/baseline_rows.tsv
OUT="${OUT:-$HOME/Desktop/AMRPredict_pseaer_strict_results.xlsx}"
TAG="${TAG:-armD1k}"
HTAG="${HTAG:-armD_HOLDOUT}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS="${TOOLS:-${MTB_SCRIPTS:-$HERE/../tools}}"
NOTE="${NOTE:-1-SE on cost=2*vMER+MER, ties broken on bal_OT, median aggregation; arm D declared in advance}"

DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)
MODELS=(evo2 carbon dnabert)

rm -f "$OUT"
for DRUG in "${DRUGS[@]}"; do
  B="$BASE/cv_${DRUG}/baseline_rows.tsv"
  # Pass 1: the comparator, before the model it judges.
  if [ -e "$B" ]; then
    python3 "$TOOLS/build_results_workbook.py" --drug "$DRUG" \
      --baseline "$B" --append --out "$OUT"
  else
    echo "[INFO] $DRUG: no ARDaP baseline (not covered)"
  fi

  # Pass 2: the models. The config tag is read from the selected-config JSON so
  # the workbook and the holdout cannot disagree about which config was scored.
  args=()
  for M in "${MODELS[@]}"; do
    dev="$SLIM/cv_${DRUG}/runs/${M}_${TAG}"
    hold="$SLIM/cv_${DRUG}/runs/${M}_${HTAG}"
    cfg="$IN/scripts/final_configs_strict/final_config_${DRUG}_${M}_armD.json"
    for f in "$dev/folds.tsv" "$hold/folds.tsv" "$cfg"; do
      [ -e "$f" ] || { echo "[ERROR] $DRUG/$M: missing $f"; exit 1; }
    done
    tag=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$cfg")
    args+=(--model "${M}=${dev}:${hold}:${tag}")
  done
  python3 "$TOOLS/build_results_workbook.py" --drug "$DRUG" \
    "${args[@]}" --selection-note "$NOTE" --append --out "$OUT"
done

echo
echo "[DONE] $OUT"
