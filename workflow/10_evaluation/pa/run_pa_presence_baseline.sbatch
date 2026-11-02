#!/bin/bash
#SBATCH --job-name=pa_pres
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-9

# The catalogue-only baseline for P. aeruginosa: AMRFinderPlus presence columns
# (acquired genes and curated point mutations) alone, no embeddings, on the
# operator-aware cohort. The analogue of the M. tuberculosis mutation-column model
# (run_holdout_catalogue.sbatch): same learner, same 1,000 configs, same folds as
# the reported arm D, so arm D minus this baseline is what the embeddings add on
# top of the catalogue.
#
# The earlier variant-table baseline's `presence` mode is the same model but points at the
# drop-rule cohort and has no held-out stage; this script replaces it for that mode.
#
#   STAGE=dev      development sweep on cv_folds.tsv -> runs/base_presence
#   STAGE=holdout  select one config from that sweep by the panel's declared rule
#                  (cost = 2*vMER + MER, 1 SE, tie-break bal_OT, median), then score
#                  it once on cv_folds_final.tsv -> runs/base_presence_HOLDOUT and
#                  predictions/<drug>_presence_holdout_predictions.tsv
#
# SCORED ONCE: the held-out stage stops if results already exist.
#
# CALL=1: the same baseline with ARDaP's own call added, i.e. AMRFinderPlus presence + ARDaP call,
# no embeddings -- the P. aeruginosa counterpart of the M. tuberculosis model on mutation columns +
# TB-Profiler call. Its features are arm DR's table (presence columns, then the call), already built by
# submit_pa_dev.sh ARM=DR, so the call is encoded exactly as in the reported model. Writes runs/
# base_prescall and predictions/<drug>_prescall_holdout_predictions.tsv; aztreonam and levofloxacin,
# which ARDaP does not call, exit without running.

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
CONFIGS="${CONFIGS:-$MTB/scripts/mtb_rf1000_configs.json}"
CVROOT="${CVROOT:-$P/imputed_strict}"
PRESDIR="${PRESDIR:-$P/labels_imputed_strict/amrfinder/by_drug_tsv}"
CFGDIR="${CFGDIR:-$P/scripts/final_configs_strict}"
OUTDIR="${OUTDIR:-$P/predictions_strict}"
NJOBS="${NJOBS:-8}"
STAGE="${STAGE:?set STAGE=dev or STAGE=holdout}"
CALL="${CALL:-0}"
LABELS="${LABELS:-$P/labels_imputed_strict}"
if [ "$CALL" = 1 ]; then LAB="base_prescall"; PSTEM="prescall"; else LAB="base_presence"; PSTEM="presence"; fi

DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)
# Same S/R boundaries as submit_pa_dev.sh, in the same order as DRUGS.
BOUNDS=(4 4 3 3 -1 2 1 3 4 1)

DRUG="${DRUGS[$SLURM_ARRAY_TASK_ID]}"
B="${BOUNDS[$SLURM_ARRAY_TASK_ID]}"
D="$CVROOT/cv_${DRUG}"
TSV="$PRESDIR/${DRUG}.tsv"
if [ "$CALL" = 1 ]; then
  case "$DRUG" in aztreonam|levofloxacin) echo "[INFO] $DRUG: no ARDaP call, nothing to run"; exit 0 ;; esac
  TSV="$LABELS/ardap_features/${DRUG}_armDR.tsv"
  [ -e "$TSV" ] || { echo "[ERROR] $DRUG: missing $TSV (built by submit_pa_dev.sh ARM=DR)"; exit 1; }
  head -1 "$TSV" | tr '\t' '\n' | tail -1 | grep -q "^ardap_call_" || {
    echo "[ERROR] $DRUG: last column of $TSV is not the ARDaP call"; exit 1; }
fi
DEV="$D/runs/${LAB}"

echo "[INFO] drug=$DRUG  boundary=$B  stage=$STAGE"
echo "[INFO] features = $TSV ($(( $(head -1 "$TSV" | tr '\t' '\n' | wc -l) - 1 )) columns)"

case "$STAGE" in
  dev)
    for f in "$D/cv_folds.tsv" "$TSV" "$CONFIGS"; do
      [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; exit 1; }
    done
    $PY "$SCRIPTS/cv_run_architecture.py" \
      --feature-tsv "$TSV" \
      --folds-tsv "$D/cv_folds.tsv" \
      --configs "$CONFIGS" \
      --sr-boundary "$B" \
      --label "${DRUG}_${LAB}" \
      --out-dir "$DEV" \
      --n-jobs "$NJOBS" --seed 2026
    echo "[DONE] $DEV"
    ;;
  holdout)
    CFG="$CFGDIR/final_config_${DRUG}_${PSTEM}.json"
    RUN="$D/runs/${LAB}_HOLDOUT"
    PREDS="$OUTDIR/${DRUG}_${PSTEM}_holdout_predictions.tsv"
    for f in "$D/cv_folds_final.tsv" "$TSV" "$DEV/summary.tsv"; do
      [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; exit 1; }
    done
    if [ -e "$RUN/folds.tsv" ] || [ -e "$PREDS" ]; then
      echo "[ERROR] $DRUG: held-out results already exist ($RUN or $PREDS)."
      echo "        The held-out set is scored once. Not rescoring."
      exit 1
    fi
    # The cluster has held a stale copy of the selector without --objective
    # (see run_pa_select.sbatch); stop before the held-out set is touched.
    $PY "$SCRIPTS/select_final_config.py" --help 2>&1 | grep -q -- "--objective" || {
      echo "[ERROR] $SCRIPTS/select_final_config.py predates --objective"; exit 1; }
    mkdir -p "$CFGDIR" "$OUTDIR"
    $PY "$SCRIPTS/select_final_config.py" \
      --run-dir "$DEV" \
      --configs "$CONFIGS" \
      --aggregation median --family rf \
      --w-vmer 2.0 --w-mer 1.0 --se-multiple 1.0 \
      --objective cost --tie-break bal_OT --tie-metric bal_OT \
      --out "$CFG"
    echo "[INFO] config = $($PY -c \
      "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$CFG")"
    $PY "$SCRIPTS/cv_run_architecture.py" \
      --feature-tsv "$TSV" \
      --folds-tsv "$D/cv_folds_final.tsv" \
      --configs "$CFG" \
      --sr-boundary "$B" \
      --label "${DRUG}_${LAB}_HOLDOUT" \
      --save-predictions "$PREDS" \
      --out-dir "$RUN" \
      --n-jobs "$NJOBS" --seed 2026
    echo "[DONE] $RUN ; predictions $PREDS"
    ;;
  *) echo "[ERROR] STAGE must be dev or holdout, got '$STAGE'"; exit 1 ;;
esac
