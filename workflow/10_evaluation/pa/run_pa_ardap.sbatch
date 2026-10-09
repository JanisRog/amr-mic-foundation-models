#!/bin/bash
#SBATCH --job-name=pa_ardap
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00

# ARDaP's baseline rows on the operator-aware cohorts.
#
# The comparator has to move whenever the cohorts move. ARDaP's calls are fixed
# per isolate -- they are a property of the genome, not of our labelling -- but
# which isolates are held out is not, and a comparator scored on a different
# holdout than the model is not a comparison. This is cheap (a table join) and
# is rerun on every relabelling for that reason.
#
# It is also the fourth independent check on the rebuild. ARDaP's call for an
# isolate cannot change; only cohort membership can. So its S/R truth counts must
# match the label files drug for drug, and any drift means a cohort changed when
# it should not have.

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
PY="${PY:?set PY in config/paths.env}"
PASCRIPTS="${PASCRIPTS:-$P/scripts}"
CVROOT="${CVROOT:-$P/imputed_strict}"
CALLS="${CALLS:-$P/labels/ardap_calls.tsv}"

declare -A SRB=(
  [amikacin]=4 [aztreonam]=4 [cefepime]=3 [ceftazidime]=3
  [ciprofloxacin]=-1 [imipenem]=2 [levofloxacin]=1 [meropenem]=3
  [piperacillin-tazobactam]=4 [tobramycin]=1
)
DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)

# The calls file's location was never recorded, which is the same gap that cost
# an afternoon on the metadata package. Rather than fail with a path, say where
# the file actually is so the default above can be corrected once and for all.
if [ ! -e "$CALLS" ]; then
  echo "[ERROR] no calls file at $CALLS"
  echo "        candidates under $P:"
  find "$P" -maxdepth 4 -name "ardap_calls*.tsv" 2>/dev/null | head -10
  echo "        (pass CALLS=<path> to override)"
  exit 1
fi
[ -e "$PASCRIPTS/ardap_baseline.py" ] || {
  echo "[ERROR] missing $PASCRIPTS/ardap_baseline.py"; exit 1; }

echo "[INFO] calls  = $CALLS ($(( $(wc -l < "$CALLS") - 1 )) isolates)"
echo "[INFO] cvroot = $CVROOT"

# ARDaP does not cover the whole panel, and a drug it has no column for is not a
# failure -- it is the reason the comparison is reported as "n of 8" rather than
# "n of 10". Which drugs those are is derived from the table's own header, so the
# count in the write-up cannot drift from the file it came from.
HDR=$(head -1 "$CALLS")
COVERED=(); UNCOVERED=()
for DRUG in "${DRUGS[@]}"; do
  if printf '%s' "$HDR" | tr '\t' '\n' | grep -qxF "$DRUG"; then
    COVERED+=("$DRUG")
  else
    UNCOVERED+=("$DRUG")
  fi
done
echo "[INFO] ARDaP covers ${#COVERED[@]} of ${#DRUGS[@]}: ${COVERED[*]}"
[ "${#UNCOVERED[@]}" -eq 0 ] || echo "[INFO] no ARDaP column: ${UNCOVERED[*]}"
[ "${#COVERED[@]}" -gt 0 ] || { echo "[ERROR] no drug has a column in $CALLS"; exit 1; }

for DRUG in "${COVERED[@]}"; do
  D="$CVROOT/cv_${DRUG}"
  FINAL="$D/cv_folds_final.tsv"
  [ -e "$FINAL" ] || { echo "[ERROR] missing $FINAL"; exit 1; }
  echo
  echo "################ $DRUG (boundary ${SRB[$DRUG]}) ################"
  $PY "$PASCRIPTS/ardap_baseline.py" \
    --folds-final "$FINAL" \
    --calls "$CALLS" \
    --drug "$DRUG" \
    --boundary "${SRB[$DRUG]}" \
    --out "$D/baseline_rows.tsv"
done

n=$(ls -1 "$CVROOT"/cv_*/baseline_rows.tsv 2>/dev/null | wc -l)
echo
echo "[DONE] $n baseline_rows.tsv under $CVROOT"
[ "$n" -eq "${#COVERED[@]}" ] || {
  echo "[ERROR] expected ${#COVERED[@]}, one per covered drug"; exit 1; }
