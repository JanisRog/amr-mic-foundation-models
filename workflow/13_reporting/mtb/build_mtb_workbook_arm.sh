#!/usr/bin/env bash
# The M. tuberculosis results workbook for a feature arm, from downloaded outputs.
#
# Written for the adopted model (wider target set + rule bit, §27). Development is
# the slimmed selected configuration, held-out is that configuration scored once,
# and the incumbent is TB-Profiler from each drug's baseline_rows.tsv. The earlier
# model's workbook (AMRPredict_panel_results.xlsx) is a separate file and is not
# touched.
#
#   SLIM   slim_arms download (development, selected configs)
#   HOLD   download holding cv_*/runs/*_HOLDOUT
#   OUT    workbook to write (replaced if present)

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS="${TOOLS:-${MTB_SCRIPTS:-$HERE/../tools}}"
SLIM="${SLIM:?set SLIM to the slim_arms directory}"
HOLD="${HOLD:?set HOLD to the held-out download}"
OUT="${OUT:?set OUT to the workbook path}"
ARM="${ARM:-gwasfiltered_gwfrule1k}"
NOTE="${NOTE:-1-SE on cost=2*vMER+MER, ties broken on bal_OT, median aggregation; wider target set (TB-Profiler + CRyPTIC GWAS loci) + rule bit, adopted in protocol section 27}"

# drug | fold directory | summaries directory
SPECS=(
  "isoniazid cv_isoniazid_final cv_final_summaries"
  "rifampicin cv_rifampicin cv_rif_summaries"
  "ethambutol cv_ethambutol cv_emb_summaries"
  "ethionamide cv_ethionamide cv_eth_summaries"
  "kanamycin cv_kanamycin cv_kan_summaries"
  "amikacin cv_amikacin cv_ami_summaries"
  "levofloxacin cv_levofloxacin cv_lev_summaries"
  "moxifloxacin cv_moxifloxacin cv_mxf_summaries"
  "bedaquiline cv_bedaquiline cv_bdq_summaries"
  "delamanid cv_delamanid cv_dlm_summaries"
)

[ ! -e "$OUT" ] || { mv "$OUT" "$OUT.previous"; echo "[INFO] existing workbook kept as $OUT.previous"; }
for spec in "${SPECS[@]}"; do
  set -- $spec
  DRUG="$1"; FD="$2"; SD="$3"
  B="$HERE/$SD/baseline_rows.tsv"
  [ -e "$B" ] || { echo "[ERROR] $DRUG: missing $B"; exit 1; }
  # One call per drug, models and baseline together. build_results_workbook.py
  # drops every existing row of a drug when that drug is written again, so a
  # baseline written in a first call is deleted by a second call with the models.
  args=()
  for M in evo2 carbon dnabert; do
    dev="$SLIM/$FD/runs/${M}_${ARM}"
    hold="$HOLD/$FD/runs/${M}_${ARM}_HOLDOUT"
    cfg="$SLIM/configs/${FD#cv_}__${M}__${ARM}.json"
    for f in "$dev/folds.tsv" "$hold/folds.tsv" "$cfg"; do
      [ -e "$f" ] || { echo "[ERROR] $DRUG/$M: missing $f"; exit 1; }
    done
    tag=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$cfg")
    args+=(--model "${M}=${dev}:${hold}:${tag}")
  done
  python3 "$TOOLS/build_results_workbook.py" --drug "$DRUG" "${args[@]}" \
    --baseline "$B" --selection-note "$NOTE" --append --out "$OUT"
done

# The comparator must be in the result: its absence is silent everywhere downstream.
python3 - "$OUT" "${#SPECS[@]}" <<'CHECKEOF'
import sys
from openpyxl import load_workbook
rows = list(load_workbook(sys.argv[1], read_only=True)["summary"].iter_rows(values_only=True))[1:]
enc = ("evo2", "carbon", "dnabert")
base = {r[0] for r in rows if r[1] not in enc}
models = {r[0] for r in rows if r[1] in enc}
want = int(sys.argv[2])
if len(base) != want or len(models) != want:
    sys.exit(f"[ERROR] workbook has model rows for {len(models)} drugs and baseline rows "
             f"for {len(base)}; expected {want} of each")
print(f"[OK] {want} drugs, each with model rows and baseline rows ({len(rows)} summary rows)")
CHECKEOF
echo
echo "[DONE] $OUT"
