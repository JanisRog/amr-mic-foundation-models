#!/bin/bash
#SBATCH --job-name=rif_dev
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Rifampicin development CV: three models, 1,000 configs each, 3x5 folds.
#
# Fixed before submission and not to be revisited afterwards:
#   target set   `both` -- rpoB, rpoC, rpoA, nusG as gene + 300bp upstream
#                (8 targets, 33 features with the rule bit)
#   embeddings   AF 0.1 consensus, ml_inputs/embeddings/panel_af010
#   rule bit     TB-Profiler's own rifampicin call, NOT a reconstruction
#   boundary     -1.0, the 0.5 mg/L ECOFF; phenotype and catalogue agree on
#                96.4% of isolates there
#   PCA          shared4, fitted inside each fold's training rows
#   objective    minimise 2*vMER+MER, ties within 1 SE broken on bal_OT
#   primary      evo2, declared in advance
#
# Reported ordinal metric is **bal_OT_common**, declared here because MIC classes
# 0, 1 and 2 hold 34, 45 and 40 isolates and would be 7-9 in the holdout -- three
# of ten classes, and bal_OT weights classes equally. bal_OT is reported beside
# it. Cost is unaffected, being an S/R quantity.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
R="$ROOT/cv_rifampicin"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="$MTB/scripts"
EMB_SUBDIR="${EMB_SUBDIR:-ml_inputs/embeddings/panel_af010}"
SR_BOUNDARY="${SR_BOUNDARY:--1.0}"
TOPK="${TOPK:-8}"
NJOBS="${NJOBS:-16}"
TAG="${TAG:-final1k}"

MICTSV="$ROOT/rifampicin_log2_mic.tsv"
RULE="$ROOT/rule_bit_rif_tbprofiler.tsv"
FIXED="$ROOT/panel_targets/targets_rifampicin_both.txt"
CONFIGS="$SCRIPTS/mtb_rf1000_configs.json"

for f in "$MICTSV" "$RULE" "$FIXED" "$CONFIGS"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done
echo "[INFO] mic     = $MICTSV"
echo "[INFO] rule    = $RULE"
echo "[INFO] rule columns: $(head -1 "$RULE")"
echo "[INFO] fixed   = $FIXED ($(grep -vc '^#' "$FIXED") targets)"
echo "[INFO] configs = $CONFIGS"
echo "[INFO] boundary = $SR_BOUNDARY"

# ---- folds -----------------------------------------------------------------
if [ ! -e "$R/cv_folds.tsv" ]; then
  echo "=== building folds"
  $PY "$SCRIPTS/cv_build_folds.py" \
    --mic-from "$MICTSV" \
    --holdout-frac 0.2 --folds 5 --repeats 3 --seed 2026 \
    --out-dir "$R"
else
  echo "[INFO] reusing $R/cv_folds.tsv"
fi
for f in "$R/cv_folds.tsv" "$R/cv_folds_final.tsv"; do
  [ -e "$f" ] || { echo "[ERROR] fold builder did not write $f"; exit 1; }
done
$PY - "$R/cv_folds.tsv" "$R/cv_folds_final.tsv" <<'PYEOF'
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
  CACHE="$R/cache_${M}_both"
  if [ ! -e "$CACHE/precompute_manifest.json" ]; then
    $PY "$SCRIPTS/cv_precompute_fold_features.py" \
      --root "${MROOT[$M]}" \
      --embedding-subdir "$EMB_SUBDIR" \
      --folds-tsv "$R/cv_folds.tsv" \
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
    --label "${M}_${TAG}" \
    --out-dir "$R/runs/${M}_${TAG}" \
    --n-jobs "$NJOBS" --seed 2026
done

echo
echo "[DONE] $R/runs/{evo2,carbon,dnabert}_${TAG}"
