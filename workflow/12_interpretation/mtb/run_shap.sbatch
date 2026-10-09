#!/bin/bash
#SBATCH --job-name=mtb_shap
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00

# Per-isolate SHAP for one drug's reported model, summed per target.
#
# Isoniazid keeps its own conventions (section 15): fold directory
# cv_isoniazid_final, cache cache_<model>_shared4_af010_HOLDOUT, 24 fixed targets
# sliced to the first 12, and since the redo its rule bit is TB-Profiler's own
# call. Every other drug follows the panel layout, so both are parameterised
# rather than branched on.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="$MTB/scripts"
MODEL="${MODEL:-evo2}"
NJOBS="${NJOBS:-16}"

: "${DRUG:?set DRUG, e.g. isoniazid}"
: "${SR_BOUNDARY:?set SR_BOUNDARY explicitly -- it has no safe default}"

if [ "$DRUG" = "isoniazid" ]; then
  CV="${CVDIR:-$ROOT/cv_isoniazid_final}"
  CACHE="${CACHE:-$CV/cache_${MODEL}_shared4_af010_HOLDOUT}"
  CFG="${CFG:-$CV/final_config_${MODEL}_tbrule.json}"
  RULE="${RULE:-$ROOT/rule_bit_inh_tbprofiler.tsv}"
  TOPK="${TOPK:-12}"
else
  : "${CODE:?set CODE for a panel drug, e.g. rif}"
  CV="${CVDIR:-$ROOT/cv_${DRUG}}"
  CACHE="${CACHE:-$CV/cache_${MODEL}_both_HOLDOUT}"
  CFG="${CFG:-$SCRIPTS/final_config_${CODE}_${MODEL}.json}"
  RULE="${RULE:-$ROOT/rule_bit_${CODE}_tbprofiler.tsv}"
  TOPK="${TOPK:-$(grep -vc '^#' "$ROOT/panel_targets/targets_${DRUG}_both.txt")}"
fi
# Gene symbols are cosmetic -- the analysis is keyed on locus tags -- so the
# path is resolved if it exists and dropped if it does not, rather than being
# allowed to fail a job whose expensive step has already finished.
CLUMPS="${CLUMPS:-}"
if [ -z "$CLUMPS" ]; then
  # Real paths, found rather than guessed. panel_target_clumps.tsv sits beside
  # panel_targets/ and covers the panel; the gwas selections copy is the
  # genome-wide superset and is the fallback.
  for c in "$ROOT/panel_target_clumps.tsv" \
           "$MTB/results/gwas/selections/mtb_all_gene_target_clumps/all_gene_target_clumps.tsv" \
           "$ROOT/embedset_target_clumps.tsv"; do
    [ -e "$c" ] && { CLUMPS="$c"; break; }
  done
fi
[ -n "$CLUMPS" ] && echo "[INFO] clumps $CLUMPS" || echo "[INFO] no clumps file found; locus tags only"

for f in "$CACHE/precompute_manifest.json" "$CFG" "$RULE"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done
echo "[INFO] drug $DRUG  model $MODEL  top-k $TOPK  boundary $SR_BOUNDARY"
echo "[INFO] cache  $CACHE"
echo "[INFO] config $CFG"
echo "[INFO] rule   $RULE"

mkdir -p "$ROOT/shap"
$PY "$SCRIPTS/shap_target_values.py" \
  --cache-dir "$CACHE" \
  --configs "$CFG" \
  --extra-feature-tsv "$RULE" \
  --mode shared4 --top-k "$TOPK" \
  --sr-boundary "$SR_BOUNDARY" \
  ${CLUMPS:+--clumps "$CLUMPS"} \
  --n-jobs "$NJOBS" --seed 2026 \
  --out-prefix "$ROOT/shap/${DRUG}_${MODEL}"

echo
echo "[DONE] $ROOT/shap/${DRUG}_${MODEL}_shap_long.tsv"
