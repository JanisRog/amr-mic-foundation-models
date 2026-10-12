#!/usr/bin/env bash
# TreeSHAP for the adopted M. tuberculosis model: the wider target set + rule bit.
#
# run_shap.sbatch defaults describe the model first reported (target set `both`,
# final_config_<code>_<model>.json, output shap/<drug>_<model>). This passes the
# wider arm's held-out cache, its development-selected configuration and its
# target count, and writes to shap_<tag>/ so the earlier SHAP files stay intact.
#
#   DRYRUN=1 bash submit_mtb_shap_arm.sh     print what would be submitted
#   ONLY="rifampicin ethambutol" bash ...    a subset

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
MODEL="${MODEL:-evo2}"
TARGETSET="${TARGETSET:-gwasfiltered}"
TAG="${TAG:-gwfrule1k}"
DRYRUN="${DRYRUN:-}"
ONLY="${ONLY:-}"

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

for spec in "${SPECS[@]}"; do
  set -- $spec
  DRUG="$1"; CODE="$2"; SRB="$3"; FOLDDIR="$4"
  if [ -n "$ONLY" ] && ! printf '%s\n' $ONLY | grep -qx "$DRUG"; then continue; fi
  CV="$ROOT/$FOLDDIR"
  CACHE="$CV/cache_${MODEL}_${TARGETSET}_HOLDOUT"
  CFG="$ROOT/slim_arms/configs/${FOLDDIR#cv_}__${MODEL}__${TARGETSET}_${TAG}.json"
  FIXED="$ROOT/panel_targets/targets_${DRUG}_${TARGETSET}.txt"
  OUTPREFIX="$ROOT/shap_${TAG}/${DRUG}_${MODEL}"
  for f in "$CACHE/precompute_manifest.json" "$CFG" "$FIXED"; do
    [ -e "$f" ] || { echo "[ERROR] $DRUG: missing $f"; exit 1; }
  done
  TOPK=$(grep -vc '^#' "$FIXED")
  if [ -n "$DRYRUN" ]; then
    printf '%-13s topk=%-3s boundary=%-5s cache=%s out=%s\n' "$DRUG" "$TOPK" "$SRB" \
      "$(basename "$CACHE")" "${OUTPREFIX#$ROOT/}"
    continue
  fi
  jid=$(sbatch --parsable --job-name="shap_${TAG}_${CODE}" \
        --export="ALL,DRUG=$DRUG,CODE=$CODE,SR_BOUNDARY=$SRB,MODEL=$MODEL,CVDIR=$CV,\
CACHE=$CACHE,CFG=$CFG,TOPK=$TOPK,OUTPREFIX=$OUTPREFIX,ROOT=$ROOT" \
        "$SCRIPTS/run_shap.sbatch" < /dev/null)
  echo "${DRUG}_shap=$jid"
done
