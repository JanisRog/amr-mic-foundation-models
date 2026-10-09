#!/usr/bin/env bash
# TreeSHAP for every P. aeruginosa drug's reported model.
#
# The model is not evo2 throughout: it is chosen per drug on development cost,
# and picking it here with the same script the figures use means the beeswarm and
# the scatter beside it in the paper cannot be different models. --config-dir is
# not optional on the cluster -- a full run directory's folds.tsv holds all 1,000
# configs, and averaging over the grid instead of over the selected config would
# silently pick a different model than the panel reports.
#
#   DRYRUN=1 ./submit_pa_shap.sh   print the chosen models without submitting

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
PY="${PY:?set PY in config/paths.env}"
CVROOT="${CVROOT:-$P/imputed_strict}"
CFGDIR="${CFGDIR:-$P/scripts/final_configs_strict}"
OUTDIR="${OUTDIR:-$P/shap_strict}"
TAG="${TAG:-armD1k}"
SCRIPT="${SCRIPT:-$P/scripts/run_pa_shap.sbatch}"
PICK="${PICK:-$P/scripts/pick_reported_model.py}"
DRYRUN="${DRYRUN:-}"

declare -A SRB=(
  [amikacin]=4 [aztreonam]=4 [cefepime]=3 [ceftazidime]=3
  [ciprofloxacin]=-1 [imipenem]=2 [levofloxacin]=1 [meropenem]=3
  [piperacillin-tazobactam]=4 [tobramycin]=1
)
DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)

for f in "$SCRIPT" "$PICK"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done
mkdir -p "$OUTDIR"

for DRUG in "${DRUGS[@]}"; do
  M=$($PY "$PICK" --run-dir "$CVROOT/cv_${DRUG}" --tag "$TAG" \
        --config-dir "$CFGDIR" --drug "$DRUG")
  [ -n "$M" ] || { echo "[ERROR] no model for $DRUG"; exit 1; }

  CACHE="$CVROOT/cv_${DRUG}/cache_${M}_embseq_HOLDOUT"
  [ -e "$CACHE/precompute_manifest.json" ] || {
    echo "[ERROR] $DRUG/$M: missing $CACHE/precompute_manifest.json"; exit 1; }

  if [ -n "$DRYRUN" ]; then
    printf '%-24s model=%-8s boundary=%s\n' "$DRUG" "$M" "${SRB[$DRUG]}"
    continue
  fi
  # `< /dev/null`: sbatch reads stdin, and a loop that feeds it the remaining
  # drug names loses every other job to a 0-second FAILED with no error message.
  jid=$(sbatch --parsable --job-name="pas_${DRUG:0:6}" \
        --export="ALL,DRUG=$DRUG,MODEL=$M,SR_BOUNDARY=${SRB[$DRUG]},CVROOT=$CVROOT,CFGDIR=$CFGDIR,OUTDIR=$OUTDIR" \
        "$SCRIPT" < /dev/null)
  echo "${DRUG}_shap=$jid"
done
