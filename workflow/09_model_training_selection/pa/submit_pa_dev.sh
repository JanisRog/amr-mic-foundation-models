#!/usr/bin/env bash
# Submit arm D's development sweep for all ten P. aeruginosa drugs.
#
# Typed into a submit block the first time, which is why §8 can say what the
# 2025-06 pipeline's hyperparameters were but not what cohort or feature set they
# belonged to. The per-drug parameters that cannot be derived -- the S/R
# breakpoint above all -- live here, in version control, next to the arm they
# describe.
#
# Arm D = curated chromosomal targets + acquired-gene sequence embeddings
# (FIXED/EMB_SUBDIR) + the acquired-gene presence block (EXTRA). It is the
# complete nested arm; see run_pa_holdout.sbatch on why it is declared rather
# than chosen.
#
# ARM=C runs the same embeddings without the presence block (chromosomal +
# acquired-gene sequence embeddings only): the embeddings-only model compared
# across the three language models in Figure 1. It reuses arm D's cache.
#
#   ./submit_pa_dev.sh            submit
#   DRYRUN=1 ./submit_pa_dev.sh   print the submissions and the derived TOPK

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
CVROOT="${CVROOT:-$P/imputed_strict}"
LABELS="${LABELS:-$P/labels_imputed_strict}"
ARM="${ARM:-D}"
case "$ARM" in C|D) ;; *) echo "[ERROR] ARM=$ARM: C or D"; exit 1 ;; esac
TAG="${TAG:-arm${ARM}1k}"
SCRIPT="${SCRIPT:-$P/scripts/run_pa_dev.sbatch}"
DRYRUN="${DRYRUN:-}"

# log2 mg/L above which an isolate is resistant. run_pa_dev.sbatch refuses to
# start without this, because it is the parameter that fails silently: a wrong
# breakpoint produces vMER and MER that look entirely plausible. Derived from the
# R/S labels in the source metadata and cross-checked against the cipro package's
# hardcoded -1.0; see §1.
declare -A SRB=(
  [amikacin]=4 [aztreonam]=4 [cefepime]=3 [ceftazidime]=3
  [ciprofloxacin]=-1 [imipenem]=2 [levofloxacin]=1 [meropenem]=3
  [piperacillin-tazobactam]=4 [tobramycin]=1
)
DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)

[ -e "$SCRIPT" ] || { echo "[ERROR] missing $SCRIPT"; exit 1; }

for DRUG in "${DRUGS[@]}"; do
  FIXED="$P/targets/combined/${DRUG}_all.txt"
  ZEROFILL="$P/targets/amrseq_targets.txt"
  EXTRA="$LABELS/amrfinder/by_drug_tsv/${DRUG}.tsv"
  [ "$ARM" = D ] || EXTRA=""
  FOLDS="$CVROOT/cv_${DRUG}/cv_folds.tsv"

  miss=0
  for f in "$FIXED" "$ZEROFILL" ${EXTRA:+"$EXTRA"} "$FOLDS"; do
    [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; miss=1; }
  done
  [ "$miss" -eq 0 ] || { echo "[ERROR] nothing submitted"; exit 1; }

  # Derived from the target file rather than typed. TOPK slices cached PCA
  # columns, so a value below the target count silently trains on a prefix of the
  # loci; deriving it from the same file the cache is built from is the only
  # version that cannot drift.
  TOPK=$(grep -vc '^#' "$FIXED")

  EXPORTS="ALL,DRUG=$DRUG,TOPK=$TOPK,SR_BOUNDARY=${SRB[$DRUG]},CVROOT=$CVROOT"
  EXPORTS="$EXPORTS,TAG=$TAG,CACHE_TAG=embseq"
  EXPORTS="$EXPORTS,EMB_SUBDIR=ml_inputs/embeddings/curated_plus_amrseq"
  EXPORTS="$EXPORTS,FIXED=$FIXED,ZEROFILL=$ZEROFILL,EXTRA=$EXTRA"

  if [ -n "$DRYRUN" ]; then
    printf '%-24s topk=%-4s boundary=%-4s %s\n' "$DRUG" "$TOPK" "${SRB[$DRUG]}" "$EXTRA"
    continue
  fi
  # `< /dev/null` is not optional inside a loop: sbatch reads stdin, and a loop
  # that feeds it the remaining drug names loses every other job to a 0-second
  # FAILED with no error message. That cost five drugs once already.
  jid=$(sbatch --parsable --job-name="pa${ARM}_${DRUG:0:6}" \
        --export="$EXPORTS" "$SCRIPT" < /dev/null)
  echo "${DRUG}_submit=$jid"
done
