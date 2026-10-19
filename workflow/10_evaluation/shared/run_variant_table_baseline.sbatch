#!/bin/bash
#SBATCH --job-name=varonly
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-9

# "Variants only": the comparator for "embeddings only", same loci, same calls.
#
# Every variant in the reported model's target loci (the AF 0.1 Snippy calls the
# consensus sequences were embedded from) as one binary column, no catalogue call,
# no curation -- a variant is a column if >= 10 and <= 95% of development isolates
# carry it. Same 1,000 random-forest configs, same folds, same S/R boundaries as
# the embeddings-only models, so the two differ in representation only.
#
#   ORG=mtb  targets panel_targets/targets_<drug>_gwasfiltered.txt (the reported set)
#            compared with <model>_gwasfiltered_embonly1k
#   ORG=pa   chromosomal targets of the reported arm (targets/combined/<drug>_all.txt;
#            acquired genes have no coordinates and are carried instead as
#            AMRFinderPlus gene presence, curated point mutations dropped)
#            compared with arm C (armC1k)
#
#   STAGE=dev      build the table over all isolates of cv_folds_final.tsv (floor
#                  counted on development rows only), then the development sweep
#   STAGE=holdout  select one config by the panel's rule (cost = 2*vMER + MER, 1 SE,
#                  tie-break bal_OT, median) and score it once on the held-out set
#
# The table has 10-20x the columns the grid was built for and `max_features` is a
# fraction, so a sweep can take much longer than the embedding runs did. The log
# prints a timestamp at start and finish; the first drugs to finish give the rest.
#
# SCORED ONCE: the held-out stage stops if results already exist.

set -euo pipefail

: "${ORG:?set ORG=mtb or ORG=pa}"
: "${STAGE:?set STAGE=dev or STAGE=holdout}"
PY="${PY:?set PY in config/paths.env}"
MTB="${MTB:?set MTB in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
CONFIGS="${CONFIGS:-$MTB/scripts/mtb_rf1000_configs.json}"
NJOBS="${NJOBS:-16}"
MINCAR="${MINCAR:-10}"
MAXFRAC="${MAXFRAC:-0.95}"
I="$SLURM_ARRAY_TASK_ID"
EXTRA_ARGS=()

case "$ORG" in
  mtb)
    ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
    SPEC=("isoniazid -3.0 cv_isoniazid_final" "rifampicin -1.0 cv_rifampicin"
          "ethambutol 2.0 cv_ethambutol" "ethionamide 2.0 cv_ethionamide"
          "kanamycin 2.0 cv_kanamycin" "amikacin 0.0 cv_amikacin"
          "levofloxacin 0.0 cv_levofloxacin" "moxifloxacin 0.0 cv_moxifloxacin"
          "bedaquiline -2.0 cv_bedaquiline" "delamanid -3.0 cv_delamanid")
    read -r DRUG B FOLDDIR <<< "${SPEC[$I]}"
    D="$ROOT/$FOLDDIR"
    TARGETS="$ROOT/panel_targets/targets_${DRUG}_gwasfiltered.txt"
    CLUMPS="${CLUMPS:-$MTB/results/gwas/selections/mtb_all_gene_target_clumps/all_gene_target_clumps.tsv}"
    RUNS="${RUNS:-$MTB/results/gwas/snippy_h37rv/runs}"
    TABLE="$ROOT/variant_tables/${DRUG}_gwasfiltered_variants.tsv"
    LAB="variants_gwasfiltered_varonly1k"
    CFG="$ROOT/variant_tables/configs/final_config_${DRUG}_varonly.json"
    PREDS="$ROOT/predictions/${DRUG}_varonly_holdout_predictions_gwasfiltered.tsv"
    ;;
  pa)
    P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
    LABELS="${LABELS:-$P/labels_imputed_strict}"
    DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
           levofloxacin meropenem piperacillin-tazobactam tobramycin)
    BOUNDS=(4 4 3 3 -1 2 1 3 4 1)
    DRUG="${DRUGS[$I]}"; B="${BOUNDS[$I]}"
    D="${CVROOT:-$P/imputed_strict}/cv_${DRUG}"
    TARGETS="$P/targets/combined/${DRUG}_all.txt"
    CLUMPS="${CLUMPS:-$P/targets/all_gene_target_clumps.tsv}"
    RUNS="${RUNS:-$P/results/gwas/snippy_hybrid/latest/runs}"
    PRES="$LABELS/amrfinder/by_drug_tsv/${DRUG}.tsv"
    # amrfinder_columns.tsv is written next to by_drug_tsv/ by build_amrfinder_presence.py
    POINTMETA="${POINTMETA:-$LABELS/amrfinder/amrfinder_columns.tsv}"
    # PRESENCE=0: chromosomal variants alone (the partner of arm A), no acquired genes.
    if [ "${PRESENCE:-1}" = 0 ]; then
      LAB="varchromo"
    else
      LAB="varonly"
      EXTRA_ARGS=(--extra-tsv "$PRES" --drop-point-from "$POINTMETA")
    fi
    # The full table keeps the name its completed runs were built with.
    if [ "$LAB" = varonly ]; then TABLE="$LABELS/variant_tables/${DRUG}_variants.tsv"
    else TABLE="$LABELS/variant_tables/${DRUG}_${LAB}.tsv"; fi
    CFG="${CFGDIR:-$P/scripts/final_configs_strict}/final_config_${DRUG}_${LAB}.json"
    PREDS="${OUTDIR:-$P/predictions_strict}/${DRUG}_${LAB}_holdout_predictions.tsv"
    ;;
  *) echo "[ERROR] ORG must be mtb or pa, got '$ORG'"; exit 1 ;;
esac

DEV="$D/runs/${LAB}"
echo "[INFO] org=$ORG  drug=$DRUG  boundary=$B  stage=$STAGE"
echo "[INFO] targets = $TARGETS ($(grep -vc '^#' "$TARGETS"))"
echo "[INFO] table   = $TABLE"

case "$STAGE" in
  dev)
    NEED=("$D/cv_folds.tsv" "$D/cv_folds_final.tsv" "$TARGETS" "$CLUMPS" "$CONFIGS")
    if [ "$ORG" = pa ] && [ "${PRESENCE:-1}" != 0 ]; then NEED+=("$PRES" "$POINTMETA"); fi
    for f in "${NEED[@]}"; do
      [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; exit 1; }
    done
    [ -d "$RUNS" ] || { echo "[ERROR] no such Snippy runs directory: $RUNS"; exit 1; }
    probe=$(awk -F'\t' 'NR==2{print $1}' "$D/cv_folds_final.tsv")
    [ -e "$RUNS/$probe/snps.af010.vcf.gz" ] || {
      echo "[ERROR] $RUNS/$probe/snps.af010.vcf.gz does not exist -- wrong layout or AF"; exit 1; }
    if [ -e "$DEV/summary.tsv" ]; then
      echo "[ERROR] $DEV already holds results; not overwriting"; exit 1; fi
    echo "[INFO] build started  = $(date '+%F %T')"
    $PY "$SCRIPTS/build_locus_variant_table.py" \
      --snippy-runs "$RUNS" --vcf-name snps.af010.vcf.gz \
      --clumps "$CLUMPS" --targets "$TARGETS" \
      --rows "$D/cv_folds_final.tsv" --floor-fold -1 \
      --min-carriers "$MINCAR" --max-carrier-frac "$MAXFRAC" \
      "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}" \
      --out "$TABLE"
    echo "[INFO] sweep started  = $(date '+%F %T')"
    $PY "$SCRIPTS/cv_run_architecture.py" \
      --feature-tsv "$TABLE" \
      --folds-tsv "$D/cv_folds.tsv" \
      --configs "$CONFIGS" \
      --sr-boundary "$B" \
      --label "${DRUG}_${LAB}" \
      --out-dir "$DEV" \
      --n-jobs "$NJOBS" --seed 2026
    echo "[INFO] finished       = $(date '+%F %T')"
    echo "[DONE] $DEV"
    ;;
  holdout)
    RUN="$D/runs/${LAB}_HOLDOUT"
    for f in "$D/cv_folds_final.tsv" "$TABLE" "$DEV/summary.tsv"; do
      [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; exit 1; }
    done
    if [ -e "$RUN/folds.tsv" ] || [ -e "$PREDS" ]; then
      echo "[ERROR] $DRUG: held-out results already exist ($RUN or $PREDS)."
      echo "        The held-out set is scored once. Not rescoring."
      exit 1
    fi
    $PY "$SCRIPTS/select_final_config.py" --help 2>&1 | grep -q -- "--objective" || {
      echo "[ERROR] $SCRIPTS/select_final_config.py predates --objective"; exit 1; }
    mkdir -p "$(dirname "$CFG")" "$(dirname "$PREDS")"
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
      --feature-tsv "$TABLE" \
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
