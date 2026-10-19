#!/usr/bin/env bash
# Submit arm D's holdout for all ten drugs. Companion to submit_pa_dev.sh, same
# breakpoint table, same reason for existing.
#
# Scored once, on the config selected from development and nothing else. The
# holdout is not a place to compare arms, objectives or aggregations -- each of
# those comparisons is a max over the held-out isolates, which is the bias the
# whole selection procedure exists to avoid.
#
# ARM=C scores the embeddings-only arm for the model comparison only; arm D
# stays the reported arm. Its predictions get an _armC suffix so they never
# overwrite arm D's.
#
#   DEPEND=<jobid> ./submit_pa_holdout.sh   wait for the selection job
#   DRYRUN=1 ./submit_pa_holdout.sh         print without submitting

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
CVROOT="${CVROOT:-$P/imputed_strict}"
LABELS="${LABELS:-$P/labels_imputed_strict}"
CFGDIR="${CFGDIR:-$P/scripts/final_configs_strict}"
OUTDIR="${OUTDIR:-$P/predictions_strict}"
ARM="${ARM:-D}"
# ARM=CR / DR (exploratory): arm C / arm D plus ARDaP's own call as one column
# (make_ardap_feature_tsv.py, S=0 I=1 R=2). ARDaP covers 8 drugs; Evo 2 only.
case "$ARM" in A|C|D|CR|DR) ;; *) echo "[ERROR] ARM=$ARM: A, C, D, CR or DR"; exit 1 ;; esac
if [ -z "${MODELS:-}" ]; then
  if [ "$ARM" = CR ] || [ "$ARM" = DR ]; then MODELS=evo2; else MODELS="evo2 carbon dnabert"; fi
fi
TAG="${TAG:-arm${ARM}_HOLDOUT}"
SCRIPT="${SCRIPT:-$P/scripts/run_pa_holdout.sbatch}"
DEPEND="${DEPEND:-}"
DRYRUN="${DRYRUN:-}"

declare -A SRB=(
  [amikacin]=4 [aztreonam]=4 [cefepime]=3 [ceftazidime]=3
  [ciprofloxacin]=-1 [imipenem]=2 [levofloxacin]=1 [meropenem]=3
  [piperacillin-tazobactam]=4 [tobramycin]=1
)
DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)

[ -e "$SCRIPT" ] || { echo "[ERROR] missing $SCRIPT"; exit 1; }
mkdir -p "$OUTDIR"

for DRUG in "${DRUGS[@]}"; do
  FIXED="$P/targets/combined/${DRUG}_all.txt"
  EXTRA="$LABELS/amrfinder/by_drug_tsv/${DRUG}.tsv"
  [ "$ARM" = D ] || EXTRA=none
  ZF="$P/targets/amrseq_targets.txt"; CACHE_TAG=embseq
  # Arm A: the chromosomal list written by submit_pa_dev.sh, own cache. The zero-fill
  # list stays (absent acquired genes must be declared; they are not in the list).
  if [ "$ARM" = A ]; then FIXED="$P/targets/chromosomal/${DRUG}.txt"; CACHE_TAG=chromo; fi
  if [ "$ARM" = CR ] || [ "$ARM" = DR ]; then
    case "$DRUG" in aztreonam|levofloxacin) continue ;; esac
    EXTRA="$LABELS/ardap_features/${DRUG}_arm${ARM}.tsv"     # written by submit_pa_dev.sh
  fi
  FINAL="$CVROOT/cv_${DRUG}/cv_folds_final.tsv"

  miss=0
  for f in "$FIXED" "$FINAL" $([ "$EXTRA" = none ] || echo "$EXTRA"); do
    [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; miss=1; }
  done
  # The configs may not exist yet when this is chained behind the selection job,
  # so their absence is only fatal when running unchained.
  if [ -z "$DEPEND" ]; then
    for M in $MODELS; do
      f="$CFGDIR/final_config_${DRUG}_${M}_arm${ARM}.json"
      [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; miss=1; }
    done
  fi
  [ "$miss" -eq 0 ] || { echo "[ERROR] nothing submitted"; exit 1; }

  EXPORTS="ALL,DRUG=$DRUG,SR_BOUNDARY=${SRB[$DRUG]},CVROOT=$CVROOT"
  EXPORTS="$EXPORTS,TAG=$TAG,ARM=$ARM,CACHE_TAG=$CACHE_TAG,CFGDIR=$CFGDIR,OUTDIR=$OUTDIR,MODELS=$MODELS"
  EXPORTS="$EXPORTS,EMB_SUBDIR=ml_inputs/embeddings/curated_plus_amrseq"
  EXPORTS="$EXPORTS,FIXED=$FIXED,ZEROFILL=$ZF,EXTRA=$EXTRA"

  if [ -n "$DRYRUN" ]; then
    printf '%-24s boundary=%-4s %s\n' "$DRUG" "${SRB[$DRUG]}" "$EXTRA"
    continue
  fi
  # `< /dev/null`: sbatch reads stdin, and a loop that feeds it the remaining
  # drug names loses every other job to a 0-second FAILED with no error message.
  jid=$(sbatch --parsable --job-name="ph${ARM}_${DRUG:0:6}" \
        ${DEPEND:+--dependency=afterok:$DEPEND} \
        --export="$EXPORTS" "$SCRIPT" < /dev/null)
  echo "${DRUG}_holdout=$jid"
done
