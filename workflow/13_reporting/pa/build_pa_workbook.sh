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
# The reported model mixes arms (CR for eight antimicrobials, D for two), so it is built
# from the tree stage_reported_pa.py writes: TAG=reported1k HTAG=reported_HOLDOUT CFGARM=reported.
CFGARM="${CFGARM:-armD}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS="${TOOLS:-${MTB_SCRIPTS:-$HERE/../tools}}"
NOTE="${NOTE:-1-SE on cost=2*vMER+MER, ties broken on bal_OT, median aggregation; arm D declared in advance}"

DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)
read -ra MODELS <<< "${MODELS:-evo2 carbon dnabert}"

rm -f "$OUT"
n_base=0
for DRUG in "${DRUGS[@]}"; do
  B="$BASE/cv_${DRUG}/baseline_rows.tsv"
  # Models and comparator go in ONE call. build_results_workbook.py deletes every
  # existing row of a drug when that drug is written again, so the earlier
  # two-pass form (baseline, then models) silently dropped every ARDaP row: the
  # strict workbook was built without its comparator.
  base_args=()
  if [ -e "$B" ]; then
    base_args=(--baseline "$B"); n_base=$((n_base + 1))
  else
    echo "[INFO] $DRUG: no ARDaP baseline (not covered)"
  fi

  # The config tag is read from the selected-config JSON so the workbook and the
  # holdout cannot disagree about which config was scored.
  args=()
  for M in "${MODELS[@]}"; do
    dev="$SLIM/cv_${DRUG}/runs/${M}_${TAG}"
    hold="$SLIM/cv_${DRUG}/runs/${M}_${HTAG}"
    cfg="$IN/scripts/final_configs_strict/final_config_${DRUG}_${M}_${CFGARM}.json"
    for f in "$dev/folds.tsv" "$hold/folds.tsv" "$cfg"; do
      [ -e "$f" ] || { echo "[ERROR] $DRUG/$M: missing $f"; exit 1; }
    done
    tag=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$cfg")
    args+=(--model "${M}=${dev}:${hold}:${tag}")
  done
  python3 "$TOOLS/build_results_workbook.py" --drug "$DRUG" \
    "${args[@]}" "${base_args[@]+"${base_args[@]}"}" \
    --selection-note "$NOTE" --append --out "$OUT"
done

# The comparator's absence is silent downstream, so it is checked here.
python3 - "$OUT" "$n_base" <<'CHECKEOF'
import sys
from openpyxl import load_workbook
rows = list(load_workbook(sys.argv[1], read_only=True)["summary"].iter_rows(values_only=True))[1:]
base = {r[0] for r in rows if r[1] not in ("evo2", "carbon", "dnabert")}
if len(base) != int(sys.argv[2]):
    sys.exit(f"[ERROR] baseline rows for {len(base)} drugs in the workbook, "
             f"{sys.argv[2]} baseline files were read")
print(f"[OK] comparator rows present for {len(base)} drugs ({len(rows)} summary rows)")
CHECKEOF

echo
echo "[DONE] $OUT"
