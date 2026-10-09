#!/bin/bash
#SBATCH --job-name=panel_dev
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# One drug's development CV: three models, 1,000 configs each, 3x5 folds.
#
# Drug-specific settings arrive through the environment and are echoed before any
# cache is built. SR_BOUNDARY has no default on purpose -- it is the parameter
# that fails silently, scoring against another drug's threshold while producing a
# plausible number.
#
#   DRUG        ethambutol      CODE  emb
#   TARGETSET   catgenes        TOPK  12      (must equal the target-list length)
#   SR_BOUNDARY 2.0
#
# The cohort is restricted to isolates that have panel embeddings. 27 of
# rifampicin's 8,390 did not, and the precompute failed with a cell count
# (27 x 403 targets) that described the map rather than the drug's target set.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="$MTB/scripts"
EMB_SUBDIR="${EMB_SUBDIR:-ml_inputs/embeddings/panel_af010}"
NJOBS="${NJOBS:-16}"
TAG="${TAG:-final1k}"

: "${DRUG:?set DRUG, e.g. ethambutol}"
: "${CODE:?set CODE, e.g. emb}"
: "${TARGETSET:=both}"
: "${TOPK:?set TOPK to the number of targets in the target list}"
: "${SR_BOUNDARY:?set SR_BOUNDARY explicitly -- it has no safe default}"

# Isoniazid's fold directory is `cv_isoniazid_final`, not `cv_isoniazid`: it was
# built before this script existed and the name was never reconciled. Overridable
# rather than symlinked, so the real directory stays visible in the log -- and
# note the sibling `cv_isoniazid_rifS_sizematched`, which is the lineage-shift
# experiment and NOT this drug's panel folds.
D="${CVDIR:-$ROOT/cv_${DRUG}}"
MICTSV="$ROOT/${DRUG}_log2_mic.tsv"
# Overridable so the combined catalogue arm (§17) can pass rule bit plus
# mutation columns as one file without a second copy of this script.
RULE="${RULE:-$ROOT/rule_bit_${CODE}_tbprofiler.tsv}"
FIXED="$ROOT/panel_targets/targets_${DRUG}_${TARGETSET}.txt"
CONFIGS="$SCRIPTS/mtb_rf1000_configs.json"
EMBEDDED="$ROOT/embedded_samples.txt"

for f in "$MICTSV" "$RULE" "$FIXED" "$CONFIGS" "$EMBEDDED"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done

n_fixed=$(grep -vc '^#' "$FIXED")
[ "$n_fixed" -eq "$TOPK" ] || {
  echo "[ERROR] TOPK=$TOPK but $FIXED holds $n_fixed targets."
  echo "        top-k slices cached PCA columns; a mismatch silently uses the"
  echo "        wrong number of targets rather than failing."
  exit 1; }

echo "[INFO] drug     = $DRUG ($CODE), target set '$TARGETSET'"
echo "[INFO] mic      = $MICTSV"
echo "[INFO] rule     = $RULE"
echo "[INFO] rule columns: $(head -1 "$RULE")"
echo "[INFO] fixed    = $FIXED ($n_fixed targets)"
echo "[INFO] boundary = $SR_BOUNDARY   embeddings = $EMB_SUBDIR"

# ---- folds -----------------------------------------------------------------
if [ ! -e "$D/cv_folds.tsv" ]; then
  echo "=== building folds (restricted to embedded isolates)"
  $PY "$SCRIPTS/cv_build_folds.py" \
    --mic-from "$MICTSV" \
    --restrict-samples "$EMBEDDED" \
    --holdout-frac 0.2 --folds 5 --repeats 3 --seed 2026 \
    --out-dir "$D"
else
  echo "[INFO] reusing $D/cv_folds.tsv"
fi
for f in "$D/cv_folds.tsv" "$D/cv_folds_final.tsv"; do
  [ -e "$f" ] || { echo "[ERROR] fold builder did not write $f"; exit 1; }
done
$PY - "$D/cv_folds.tsv" "$D/cv_folds_final.tsv" <<'PYEOF'
import csv, sys, collections
for p in sys.argv[1:3]:
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    s = {r["sample_id"] for r in rows}
    hold = {r["sample_id"] for r in rows if int(r["fold"]) < 0}
    rep = {r["repeat"] for r in rows}
    print(f"[INFO] {p.split('/')[-1]}: {len(s):,} isolates, "
          f"{len(hold):,} at fold<0, {len(rep)} repeat(s)")
PYEOF

declare -A MROOT=( [evo2]="$MTB/results/ml/03_mtb_evo2"
                   [carbon]="$MTB/results/ml/03_mtb"
                   [dnabert]="$MTB/results/ml/03_mtb_dnabert" )

for M in evo2 carbon dnabert; do
  echo
  echo "################ $M ################"
  CACHE="$D/cache_${M}_${TARGETSET}"
  if [ ! -e "$CACHE/precompute_manifest.json" ]; then
    $PY "$SCRIPTS/cv_precompute_fold_features.py" \
      --root "${MROOT[$M]}" \
      --embedding-subdir "$EMB_SUBDIR" \
      --folds-tsv "$D/cv_folds.tsv" \
      --pca-mode shared4 --pca-dim 4 \
      --var-threshold 0.99 --max-components 25 \
      --top-targets 256 --seed 2026 \
      --fixed-targets "$FIXED" \
      --out-dir "$CACHE"
  else
    echo "[INFO] reusing $CACHE"
  fi

  $PY "$SCRIPTS/cv_run_architecture.py" \
    --cache-dir "$CACHE" \
    --extra-feature-tsv "$RULE" \
    --mode shared4 --top-k "$TOPK" \
    --configs "$CONFIGS" \
    --sr-boundary "$SR_BOUNDARY" \
    --label "${M}_${TARGETSET}_${TAG}" \
    --out-dir "$D/runs/${M}_${TARGETSET}_${TAG}" \
    --n-jobs "$NJOBS" --seed 2026
done

echo
echo "[DONE] $D/runs/{evo2,carbon,dnabert}_${TARGETSET}_${TAG}"
