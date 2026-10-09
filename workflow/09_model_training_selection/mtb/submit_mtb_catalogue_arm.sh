#!/usr/bin/env bash
# The catalogue-variant arm, declared in §17 before it was run.
#
#   embeddings + that drug's catalogued mutation columns, and NO rule bit
#
# The point of dropping the rule bit is that the features then carry the
# catalogue's observations and none of its verdict: the model does not need
# TB-Profiler to have already decided which isolates are resistant. §11 measured
# this arm on isoniazid at the best bal_OT of any arm (0.762) and a slightly worse
# cost than embeddings + rule, which is the expected shape.
#
# WITH_RULE=1 adds the rule bit back (the 76-feature variant §17 keeps available).
#
# Development only. §17 fixes what happens next: the development comparison is
# reported for all ten drugs whichever way it goes, and the reported model does
# not change on a development win alone.
#
# The embedding cache is shared with the existing arm -- the extra columns are
# appended at scoring time, not baked into the precompute -- so this reuses
# `cache_<model>_both` and costs a sweep, not a re-embedding.
#
#   DRYRUN=1 ./submit_mtb_arm76.sh   print the merged widths without submitting

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
# Where the gathered TB-Profiler matrices actually live. Recorded here because
# it was written down nowhere, which is the same gap that cost an afternoon on
# the P. aeruginosa metadata package and another on ardap_calls.tsv.
MATDIR="${MATDIR:-$MTB/results/tbprofiler/reads/gather/matrix}"
OUTDIR="${OUTDIR:-$ROOT/catalogue_combined}"
WITH_RULE="${WITH_RULE:-0}"
if [ "$WITH_RULE" = 1 ]; then TAG="${TAG:-arm76}"; else TAG="${TAG:-arm75}"; fi
TARGETSET="${TARGETSET:-both}"
SCRIPT="${SCRIPT:-$SCRIPTS/run_panel_dev.sbatch}"
DRYRUN="${DRYRUN:-}"
# Space-separated drug names to submit instead of all ten. For re-running a
# single failure without resubmitting nine jobs that already succeeded.
ONLY="${ONLY:-}"

# drug | code | S/R boundary | fold directory.
#
# The codes are the ones the rule-bit files on disk use -- `ls rule_bit_*` is the
# source of truth, not the drug name's first three letters. Moxifloxacin is `mxf`,
# not `mox`, and getting that wrong fails only when a rule bit is actually read,
# so it survived one script that never reads one.. The code, the boundary and the
# fold directory cannot be derived; TOPK can be, and is, from the target file the
# cache was built from.
#
# **Isoniazid matches the panel conventions on none of its parameters** (§15): its
# folds are in cv_isoniazid_final, its cache is cache_<model>_shared4_af010, its
# targets are targets_mechfirst.txt and its top-k is 12 rather than the target
# count. Those conventions live in run_inh_redo.sbatch, and this script routes
# isoniazid there rather than restating them -- that script's own comment is that
# a second copy of them is a second place for them to drift.
SPECS=(
  "isoniazid inh -3.0 cv_isoniazid_final"
  "rifampicin rif -1.0 cv_rifampicin"
  "ethambutol emb 2.0 cv_ethambutol"
  "ethionamide eth 2.0 cv_ethionamide"
  "kanamycin kan 2.0 cv_kanamycin"
  "amikacin ami 0.0 cv_amikacin"
  "levofloxacin lev 0.0 cv_levofloxacin"
  "moxifloxacin mxf 0.0 cv_moxifloxacin"
  "bedaquiline bdq -2.0 cv_bedaquiline"
  "delamanid dlm -3.0 cv_delamanid"
)
INH_SCRIPT="${INH_SCRIPT:-$SCRIPTS/run_inh_redo.sbatch}"

MATRIX="$MATDIR/tbprofiler_resistance_mutation_binary.tsv"
LONG="$MATDIR/tbprofiler_amrfinder_style_long.tsv"
for f in "$MATRIX" "$LONG"; do
  [ -e "$f" ] || {
    echo "[ERROR] missing $f"
    echo "        Build the panel-wide matrices once first:"
    echo "          \$PY $SCRIPTS/build_tbprofiler_amrfinder_style_matrix.py \\"
    echo "            --variants-csv <tbprofiler variants csv> \\"
    echo "            --sample-list $ROOT/embedded_samples.txt \\"
    echo "            --outdir $MATDIR"
    echo "        The gathered variants table is normally at:"
    echo "          $MTB/results/tbprofiler/reads/gather/tbprofiler.variants.csv"
    exit 1; }
done
[ -e "$SCRIPTS/merge_catalogue_features.py" ] || {
  echo "[ERROR] missing $SCRIPTS/merge_catalogue_features.py"; exit 1; }

mkdir -p "$OUTDIR"

# Every drug is checked before any is submitted, so one dry run reports every
# problem rather than the first one.
problems=0
for spec in "${SPECS[@]}"; do
  set -- $spec
  DRUG="$1"; CODE="$2"; D="$ROOT/$4"
  if [ -n "$ONLY" ] && ! printf '%s\n' $ONLY | grep -qx "$DRUG"; then continue; fi
  NEED=("$D/cv_folds.tsv")
  [ "$DRUG" = isoniazid ] || NEED+=("$ROOT/panel_targets/targets_${DRUG}_${TARGETSET}.txt")
  [ "$WITH_RULE" = 1 ] && NEED+=("$ROOT/rule_bit_${CODE}_tbprofiler.tsv")
  for f in "${NEED[@]}"; do
    [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; problems=$((problems + 1)); }
  done
done
[ "$problems" -eq 0 ] || {
  echo "[ERROR] $problems missing input(s); nothing submitted"; exit 1; }

for spec in "${SPECS[@]}"; do
  set -- $spec
  DRUG="$1"; CODE="$2"; SRB="$3"; D="$ROOT/$4"
  if [ -n "$ONLY" ] && ! printf '%s\n' $ONLY | grep -qx "$DRUG"; then continue; fi
  RULE="$ROOT/rule_bit_${CODE}_tbprofiler.tsv"
  FIXED="$ROOT/panel_targets/targets_${DRUG}_${TARGETSET}.txt"
  MERGED="$OUTDIR/catalogue_${CODE}_combined.tsv"

  # Built as an array rather than inline expansion: a `[ ... ] && ...` that
  # yields 1 is the last command of its line, and under `set -e` that ends the
  # script with no message.
  RULE_ARGS=()
  if [ "$WITH_RULE" = 1 ]; then RULE_ARGS=(--rule-tsv "$RULE"); fi

  $PY "$SCRIPTS/merge_catalogue_features.py" \
    "${RULE_ARGS[@]+"${RULE_ARGS[@]}"}" \
    --mutation-matrix "$MATRIX" \
    --long "$LONG" \
    --drug "$DRUG" \
    --samples "$D/cv_folds.tsv" \
    --out "$MERGED"

  # Isoniazid's top-k is 12 by its own convention, not the target count.
  if [ "$DRUG" = isoniazid ]; then TOPK=12; else TOPK=$(grep -vc '^#' "$FIXED"); fi
  if [ -n "$DRYRUN" ]; then
    n=$(( $(head -1 "$MERGED" | tr '\t' '\n' | wc -l) - 1 ))
    printf '%-14s code=%-4s topk=%-4s boundary=%-5s extra=%-5s rule=%-2s %s\n' \
      "$DRUG" "$CODE" "$TOPK" "$SRB" "$n" "$WITH_RULE" "$(basename "$D")"
    continue
  fi
  # `< /dev/null`: sbatch reads stdin, and a loop that feeds it the remaining
  # drug names loses every other job to a 0-second FAILED with no error message.
  if [ "$DRUG" = isoniazid ]; then
    jid=$(sbatch --parsable --job-name="${TAG}_${CODE}" \
          --export="ALL,MODE=dev,TAG=$TAG,RULE=$MERGED,ROOT=$ROOT,TOPK=$TOPK,\
SR_BOUNDARY=$SRB" \
          "$INH_SCRIPT" < /dev/null)
  else
    jid=$(sbatch --parsable --job-name="${TAG}_${CODE}" \
          --export="ALL,DRUG=$DRUG,CODE=$CODE,SR_BOUNDARY=$SRB,TOPK=$TOPK,\
TARGETSET=$TARGETSET,TAG=$TAG,RULE=$MERGED,ROOT=$ROOT" \
          "$SCRIPT" < /dev/null)
  fi
  echo "${DRUG}_${TAG}=$jid"
done
