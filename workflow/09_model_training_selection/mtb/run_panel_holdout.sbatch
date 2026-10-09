#!/bin/bash
#SBATCH --job-name=panel_holdout
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00

# One drug's holdout: three models, one config each, scored once.
#
# Parameterised because ethionamide's version was a copy of isoniazid's and
# rifampicin's would have been a copy of that; three more drugs follow, and a
# per-drug copy is how SR_BOUNDARY ends up at the previous drug's value.
# Everything drug-specific arrives through the environment and is echoed before
# any model is fitted:
#
#   DRUG        rifampicin            (column name in TB-Profiler's matrix)
#   CODE        rif                   (names the rule-bit file)
#   TARGETSET   both                  (names the target list)
#   TOPK        8                     (targets in that list)
#   SR_BOUNDARY -1.0                  (MUST be set per drug -- see below)
#
# SR_BOUNDARY has no safe default. It is the one parameter that fails silently:
# a drug whose MIC classes straddle another drug's boundary scores against the
# wrong threshold and still produces a plausible number. There is no default here
# for that reason, and the job refuses to start without it.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="$MTB/scripts"
EMB_SUBDIR="${EMB_SUBDIR:-ml_inputs/embeddings/panel_af010}"
NJOBS="${NJOBS:-16}"

: "${DRUG:?set DRUG, e.g. rifampicin}"
: "${CODE:?set CODE, e.g. rif}"
: "${TARGETSET:=both}"
: "${TOPK:?set TOPK to the number of targets in the target list}"
: "${SR_BOUNDARY:?set SR_BOUNDARY explicitly -- it has no safe default}"

# Isoniazid's fold directory is `cv_isoniazid_final`, not `cv_isoniazid`: it was
# built before this script existed and the name was never reconciled. Overridable
# rather than symlinked, so the real directory stays visible in the log -- and
# note the sibling `cv_isoniazid_rifS_sizematched`, which is the lineage-shift
# experiment and NOT this drug's panel folds.
D="${CVDIR:-$ROOT/cv_${DRUG}}"
FOLDS_FINAL="$D/cv_folds_final.tsv"
FIXED="$ROOT/panel_targets/targets_${DRUG}_${TARGETSET}.txt"
RULE="$ROOT/rule_bit_${CODE}_tbprofiler.tsv"

for f in "$FOLDS_FINAL" "$FIXED" "$RULE"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done

n_fixed=$(grep -vc '^#' "$FIXED")
[ "$n_fixed" -eq "$TOPK" ] || {
  echo "[ERROR] TOPK=$TOPK but $FIXED holds $n_fixed targets."
  echo "        top-k slices the cached PCA columns; a mismatch silently uses"
  echo "        the wrong number of targets rather than failing."
  exit 1; }

echo "[INFO] drug     = $DRUG ($CODE), target set '$TARGETSET'"
echo "[INFO] folds    = $FOLDS_FINAL"
echo "[INFO] fixed    = $FIXED ($n_fixed targets)"
echo "[INFO] rule     = $RULE"
echo "[INFO] rule columns: $(head -1 "$RULE")"
echo "[INFO] boundary = $SR_BOUNDARY   embeddings = $EMB_SUBDIR"

declare -A MROOT=( [evo2]="$MTB/results/ml/03_mtb_evo2"
                   [carbon]="$MTB/results/ml/03_mtb"
                   [dnabert]="$MTB/results/ml/03_mtb_dnabert" )

for M in evo2 carbon dnabert; do
  echo
  echo "################ $M ################"
  DEV="$D/cache_${M}_${TARGETSET}"
  OUT="$D/cache_${M}_${TARGETSET}_HOLDOUT"
  CFG="$SCRIPTS/final_config_${CODE}_${M}.json"
  [ -e "$CFG" ] || { echo "[ERROR] missing config $CFG"; exit 1; }

  if [ ! -e "$OUT/precompute_manifest.json" ]; then
    $PY "$SCRIPTS/cv_precompute_fold_features.py" \
      --root "${MROOT[$M]}" \
      --embedding-subdir "$EMB_SUBDIR" \
      --folds-tsv "$FOLDS_FINAL" \
      --pca-mode shared4 --pca-dim 4 \
      --var-threshold 0.99 --max-components 25 \
      --top-targets 256 --seed 2026 \
      --fixed-targets "$FIXED" \
      --only-fold 0:0 \
      --out-dir "$OUT"
  else
    echo "[INFO] reusing existing $OUT"
  fi

  # the holdout cache must match the development cache in everything but folds
  $PY - "$DEV/precompute_manifest.json" "$OUT/precompute_manifest.json" <<'PYEOF'
import json, sys
dev, hold = (json.load(open(p)) for p in sys.argv[1:3])
bad = [k for k in ("n_targets_total", "n_fixed_targets", "pca_mode", "pca_dim",
                   "var_threshold", "max_components", "top_targets", "seed")
       if dev.get(k) != hold.get(k)]
if bad:
    sys.exit(f"[ERROR] holdout cache differs from development on {bad}: "
             + ", ".join(f"{k} {dev.get(k)!r} vs {hold.get(k)!r}" for k in bad))
print(f"[OK] geometry matches development "
      f"({hold['n_fixed_targets']} targets, {hold['n_samples']:,} isolates)")
PYEOF

  $PY "$SCRIPTS/cv_run_architecture.py" \
    --cache-dir "$OUT" \
    --extra-feature-tsv "$RULE" \
    --mode shared4 --top-k "$TOPK" \
    --configs "$CFG" \
    --sr-boundary "$SR_BOUNDARY" \
    --label "${M}_HOLDOUT" \
    --out-dir "$D/runs/${M}_HOLDOUT" \
    --n-jobs "$NJOBS" --seed 2026
done

echo
echo "[DONE] $D/runs/{evo2,carbon,dnabert}_HOLDOUT"
