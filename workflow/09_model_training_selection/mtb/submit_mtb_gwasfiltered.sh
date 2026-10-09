#!/usr/bin/env bash
# The co-resistance-filtered target sets, across the ten panel drugs.
#
#   targets = catalogue genes for the drug
#           + published GWAS hits for the drug
#           - any locus the catalogue annotates as ANOTHER drug's mechanism
#
# The exclusion is the point. The published GWAS reports rpoB under isoniazid with
# six significant signals, because 94.3% of rifampicin-resistant isolates are also
# isoniazid-resistant; for isoniazid the filter drops exactly rpoB, rpsL, embB and
# gid, which is the set §10 excluded by hand after leave-one-lineage-out. The rule
# reproduces that decision and generalises it.
#
# Needs no sequence extraction: all 435 union targets are already embedded at
# AF 0.10 for every panel drug (§10). It costs a precompute and a sweep per drug.
#
# Isoniazid goes through run_panel_dev.sbatch here, not run_inh_redo.sbatch. That
# script exists to hold isoniazid's cache name and target file, and this arm
# replaces both, so the only convention left to carry is its fold directory.
#
#   DRYRUN=1 ./submit_mtb_gwasfiltered.sh   print the widths without submitting

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
SRC="${SRC:-$ROOT/panel_targets_gwas}"
TARGETSET="${TARGETSET:-gwasfiltered}"

# Both arms embed the same gwasfiltered targets and differ only in which
# catalogue feature is appended at scoring time, so the comparison isolates that
# choice rather than confounding it with the target set.
#
#   CATALOGUE=rule  the rule bit -- TB-Profiler's own resistance call, one column
#   CATALOGUE=snps  its detected resistance mutations, one column each
#
# §11 measured these against each other on isoniazid's old target set: the
# columns win bal_OT by +0.012 in 14/15 folds, the rule wins cost in 11-14/15.
# They are not substitutes -- the bit carries the catalogue's grading, the columns
# carry which allele -- so running both over the new targets is the only way to
# know whether the wider set changes that balance.
CATALOGUE="${CATALOGUE:-rule}"
case "$CATALOGUE" in
  rule) TAG="${TAG:-gwfrule1k}" ;;
  snps) TAG="${TAG:-gwfsnp1k}"  ;;
  *) echo "[ERROR] CATALOGUE must be rule or snps, not '$CATALOGUE'"; exit 1 ;;
esac
SNPDIR="${SNPDIR:-$ROOT/catalogue_combined}"
SCRIPT="${SCRIPT:-$SCRIPTS/run_panel_dev.sbatch}"
DRYRUN="${DRYRUN:-}"
# Space-separated drug names to submit instead of all ten.
ONLY="${ONLY:-}"

# drug | code | S/R boundary | fold directory.
#
# The codes are the ones the rule-bit files on disk use -- `ls rule_bit_*` is the
# source of truth, not the drug name's first three letters. Moxifloxacin is `mxf`,
# not `mox`, and getting that wrong fails only when a rule bit is actually read,
# so it survived one script that never reads one.
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

# run_panel_dev.sbatch builds FIXED from $ROOT/panel_targets/, so the filtered
# sets are placed there under their own TARGETSET name rather than the script
# being taught a second directory. The names do not collide with the existing
# `both` and `catgenes` files.
mkdir -p "$ROOT/panel_targets"
for spec in "${SPECS[@]}"; do
  set -- $spec
  src="$SRC/targets_${1}_${TARGETSET}.txt"
  [ -e "$src" ] || { echo "[ERROR] missing $src"; exit 1; }
  cp -f "$src" "$ROOT/panel_targets/targets_${1}_${TARGETSET}.txt"
done

catalogue_file() {   # $1 = code
  if [ "$CATALOGUE" = rule ]; then
    echo "$ROOT/rule_bit_${1}_tbprofiler.tsv"
  else
    echo "$SNPDIR/catalogue_${1}_combined.tsv"
  fi
}

problems=0
for spec in "${SPECS[@]}"; do
  set -- $spec
  if [ -n "$ONLY" ] && ! printf '%s\n' $ONLY | grep -qx "$1"; then continue; fi
  for f in "$ROOT/$4/cv_folds.tsv" "$(catalogue_file "$2")" \
           "$ROOT/panel_targets/targets_${1}_${TARGETSET}.txt"; do
    [ -e "$f" ] || { echo "[ERROR] $1: missing $f"; problems=$((problems + 1)); }
  done
done
[ "$problems" -eq 0 ] || { echo "[ERROR] nothing submitted"; exit 1; }

for spec in "${SPECS[@]}"; do
  set -- $spec
  DRUG="$1"; CODE="$2"; SRB="$3"; CVDIR="$ROOT/$4"
  if [ -n "$ONLY" ] && ! printf '%s\n' $ONLY | grep -qx "$DRUG"; then continue; fi
  FIXED="$ROOT/panel_targets/targets_${DRUG}_${TARGETSET}.txt"
  TOPK=$(grep -vc '^#' "$FIXED")

  CAT=$(catalogue_file "$CODE")
  if [ -n "$DRYRUN" ]; then
    n=$(( $(head -1 "$CAT" | tr '\t' '\n' | wc -l) - 1 ))
    printf '%-14s code=%-4s targets=%-5s %s=%-5s boundary=%-5s %s\n' \
      "$DRUG" "$CODE" "$TOPK" "$CATALOGUE" "$n" "$SRB" "$(basename "$CVDIR")"
    continue
  fi
  # `< /dev/null`: sbatch reads stdin, and a loop that feeds it the remaining
  # drug names loses every other job to a 0-second FAILED with no error message.
  jid=$(sbatch --parsable --job-name="${TAG}_${CODE}" \
        --export="ALL,DRUG=$DRUG,CODE=$CODE,SR_BOUNDARY=$SRB,TOPK=$TOPK,\
TARGETSET=$TARGETSET,TAG=$TAG,ROOT=$ROOT,CVDIR=$CVDIR,RULE=$CAT" \
        "$SCRIPT" < /dev/null)
  echo "${DRUG}_${CATALOGUE}=$jid"
done
