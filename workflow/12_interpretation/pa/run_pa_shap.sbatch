#!/bin/bash
#SBATCH --job-name=pa_shap
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00

# Per-isolate SHAP for one P. aeruginosa drug's reported arm D model.
#
# PA is the more interesting panel for this than MTB: the feature set contains
# acquired-gene sequence embeddings and a presence block, so the attribution can
# say whether blaVIM or mexZ is doing the work rather than only chromosomal loci.
# The three kinds are separate categories in the plot.
#
#   DRUG        meropenem
#   SR_BOUNDARY 3           (no default: the parameter that fails silently)

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
MODEL="${MODEL:-evo2}"
NJOBS="${NJOBS:-16}"
CVROOT="${CVROOT:-$P/imputed_strict}"
CFGDIR="${CFGDIR:-$P/scripts/final_configs_strict}"
CACHE_TAG="${CACHE_TAG:-embseq}"

: "${DRUG:?set DRUG, e.g. meropenem}"
: "${SR_BOUNDARY:?set SR_BOUNDARY explicitly -- it has no safe default}"

D="$CVROOT/cv_${DRUG}"
CACHE="${CACHE:-$D/cache_${MODEL}_${CACHE_TAG}_HOLDOUT}"
CFG="${CFG:-$CFGDIR/final_config_${DRUG}_${MODEL}_armD.json}"
EXTRA="${EXTRA:-$P/labels_imputed_strict/amrfinder/by_drug_tsv/${DRUG}.tsv}"
FIXED="${FIXED:-$P/targets/combined/${DRUG}_all.txt}"
CLUMPS="${CLUMPS:-$P/targets/all_gene_target_clumps.tsv}"
OUTDIR="${OUTDIR:-$P/shap_strict}"

for f in "$CACHE/precompute_manifest.json" "$CFG" "$EXTRA" "$FIXED"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done
TOPK="${TOPK:-$(grep -vc '^#' "$FIXED")}"
[ -e "$CLUMPS" ] || { echo "[INFO] no clumps file; locus tags only"; CLUMPS=""; }

echo "[INFO] drug $DRUG  model $MODEL  top-k $TOPK  boundary $SR_BOUNDARY"
echo "[INFO] cache  $CACHE"
echo "[INFO] config $CFG"
echo "[INFO] extra  $EXTRA ($(( $(head -1 "$EXTRA" | tr '\t' '\n' | wc -l) - 1 )) presence columns)"

mkdir -p "$OUTDIR"
$PY "$SCRIPTS/shap_target_values.py" \
  --cache-dir "$CACHE" \
  --configs "$CFG" \
  --extra-feature-tsv "$EXTRA" \
  --mode shared4 --top-k "$TOPK" \
  --sr-boundary "$SR_BOUNDARY" \
  ${CLUMPS:+--clumps "$CLUMPS"} \
  --n-jobs "$NJOBS" --seed 2026 \
  --out-prefix "$OUTDIR/${DRUG}_${MODEL}"

echo
echo "[DONE] $OUTDIR/${DRUG}_${MODEL}_shap_long.tsv"
