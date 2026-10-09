#!/bin/bash
#SBATCH --job-name=inh_redo
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Isoniazid's development sweep, re-run on TB-Profiler's own call.
#
# Section 11: the isoniazid arms were fitted with `rule_calls/tbprofiler_calls_
# af010.tsv` -- a reconstruction of the catalogue call as "carries a catalogued
# variant at AF >= 0.1" -- while every drug from rifampicin onward was fitted
# with TB-Profiler's own graded call. The reconstruction is the weaker feature
# (development cost 0.107 against 0.102), so isoniazid has been competing with a
# handicap the other nine did not carry. This removes the inconsistency.
#
# Everything else is held fixed, so the movement is attributable to the feature:
# the same folds, the same caches, the same 1,000 configs, the same seed.
#
# A dedicated script rather than run_panel_dev.sbatch with five overrides,
# because isoniazid matches that script on none of its conventions. It predates
# all of them:
#
#   fold directory   cv_isoniazid_final          not cv_isoniazid
#   cache name       cache_<model>_shared4_af010 not cache_<model>_both
#   fixed targets    targets_mechfirst.txt (24)  not panel_targets/targets_isoniazid_both.txt (23)
#   top-k            12                          -> 12*4 + 1 rule bit = 49 features
#
# **The caches are reused, never rebuilt.** They carry the PCA fitted for the
# reported run; a rebuild would refit it and the comparison to the published
# number would no longer be a comparison of one thing.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="$MTB/scripts"
NJOBS="${NJOBS:-16}"
# TAG names the output directory and must differ between the two modes: both
# write to runs/<model>_<TAG>, so a shared default would have the holdout
# overwrite the sweep's summary.tsv -- the file the config was selected from.
TAG="${TAG:-}"
# MODE=dev sweeps the 1,000 configs on the development folds; MODE=holdout scores
# the one reselected config, once, and keeps its per-isolate predictions. Both in
# one file because isoniazid's conventions are what this script exists to hold,
# and a second copy of them is a second place for them to drift.
MODE="${MODE:-dev}"
TOPK="${TOPK:-12}"
SR_BOUNDARY="${SR_BOUNDARY:--3.0}"

CV="$ROOT/cv_isoniazid_final"
RULE="${RULE:-$ROOT/rule_bit_inh_tbprofiler.tsv}"
case "$MODE" in
  dev)     SUFFIX=""         ; TAG="${TAG:-tbrule1k}"
           CONFIGS="${CONFIGS:-$SCRIPTS/mtb_rf1000_configs.json}" ;;
  holdout) SUFFIX="_HOLDOUT" ; TAG="${TAG:-tbrule_HOLDOUT}"
           CONFIGS="${CONFIGS:-}" ;;
  *) echo "[ERROR] MODE must be dev or holdout, not '$MODE'"; exit 1 ;;
esac

for f in "$RULE" "$CV/cv_folds.tsv"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done

# The whole point of the job is which rule bit is being used, so it is named in
# the log with its columns and its resistant count before anything is fitted.
echo "[INFO] rule bit = $RULE"
echo "[INFO] columns  = $(head -1 "$RULE")"
echo "[INFO] calls    = $(awk -F'\t' 'NR>1 && $2==1' "$RULE" | wc -l) resistant of $(( $(wc -l < "$RULE") - 1 ))"
echo "[INFO] mode     = $MODE   top-k = $TOPK   boundary = $SR_BOUNDARY"
case "$(basename "$RULE")" in
  *af0*) echo "[WARN] this looks like a reconstruction, not TB-Profiler's own call."
         echo "       That is the feature this job exists to replace (section 11).";;
esac

for M in evo2 carbon dnabert; do
  echo
  echo "################ $M ################"
  CACHE="$CV/cache_${M}_shared4_af010${SUFFIX}"
  # In holdout mode each model carries its own reselected config, so the file is
  # resolved per model rather than shared.
  CFG="${CONFIGS:-$CV/final_config_${M}_tbrule.json}"
  [ -e "$CFG" ] || { echo "[ERROR] missing config $CFG"; exit 1; }
  if [ ! -e "$CACHE/precompute_manifest.json" ]; then
    echo "[ERROR] no cache at $CACHE"
    ls -1d "$CV"/cache_* 2>/dev/null | sed 's/^/          /'
    exit 1
  fi
  n_fixed=$($PY -c "import json,sys;print(json.load(open(sys.argv[1]))['n_fixed_targets'])" \
              "$CACHE/precompute_manifest.json")
  echo "[INFO] cache $CACHE ($n_fixed fixed targets, slicing the first $TOPK)"

  # The holdout is scored once and its predictions are kept in the same pass --
  # a second fit to recover them is a second chance for something to differ.
  PREDS=""
  if [ "$MODE" = "holdout" ]; then
    mkdir -p "$ROOT/predictions"
    PREDS="$ROOT/predictions/isoniazid_${M}_holdout_predictions_tbrule.tsv"
  fi

  $PY "$SCRIPTS/cv_run_architecture.py" \
    --cache-dir "$CACHE" \
    --extra-feature-tsv "$RULE" \
    --mode shared4 --top-k "$TOPK" \
    --configs "$CFG" \
    --sr-boundary "$SR_BOUNDARY" \
    --label "${M}_${TAG}" \
    ${PREDS:+--save-predictions "$PREDS"} \
    --out-dir "$CV/runs/${M}_${TAG}" \
    --n-jobs "$NJOBS" --seed 2026
done

echo
echo "[DONE] $CV/runs/{evo2,carbon,dnabert}_${TAG}"
if [ "$MODE" = "dev" ]; then
  echo "       next: select_final_config.py on each summary.tsv, then MODE=holdout."
else
  echo "       predictions in $ROOT/predictions/"
fi
