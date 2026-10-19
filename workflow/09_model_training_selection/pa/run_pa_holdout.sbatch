#!/bin/bash
#SBATCH --job-name=pa_holdout
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00

# One P. aeruginosa drug's holdout: three models, one config each, scored once,
# with every isolate's prediction kept.
#
# **Arm D by declaration, not because it won development.** Under the panel's
# own rule -- minimise 2*vMER+MER, ties within 1 SE broken on bal_OT -- arm C
# has the lower development cost on 6 of 10 drugs and D is worse by more than
# 1 SE on four (aztreonam, ciprofloxacin, meropenem, piperacillin-tazobactam).
# None of those four separates at p<0.05, and the single significant difference
# favours D (amikacin, p=0.018). D is carried forward as the complete nested arm
# -- it contains B and C -- and because it is never distinguishably worse. The
# write-up must say that, not "arm D won".
#
# Scoring both C and D here and reporting the better one would be section 3's
# max-over-configs bias reappearing at the last possible step, which is why the
# arm is fixed before this job runs.
#
#   DRUG        ciprofloxacin
#   SR_BOUNDARY -1            (no default: the parameter that fails silently)

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
NJOBS="${NJOBS:-16}"
# ARM=C: the embeddings-only arm, scored for the model comparison only, with
# EXTRA=none and predictions suffixed _armC. Arm D remains the reported arm.
ARM="${ARM:-D}"
TAG="${TAG:-arm${ARM}_HOLDOUT}"

# Arm D = the acquired-sequence cache with the presence block appended at
# scoring time. Arms differing only in their extra-feature TSV share a
# precompute, which is why the cache is named for the embeddings and not the arm.
CACHE_TAG="${CACHE_TAG:-embseq}"
EMB_SUBDIR="${EMB_SUBDIR:-ml_inputs/embeddings/curated_plus_amrseq}"

: "${DRUG:?set DRUG, e.g. ciprofloxacin}"
: "${SR_BOUNDARY:?set SR_BOUNDARY explicitly -- it has no safe default}"

# CVROOT lets a re-labelled panel live beside the original instead of
# overwriting it. The imputed-censoring rebuild uses CVROOT=$P/imputed, so the
# published cohorts stay intact and reproducible until the new ones are shown to
# be better -- and cefepime and levofloxacin, whose labels are bit-identical
# under both rules, become a regression test rather than a rerun.
D="${CVROOT:-$P}/cv_${DRUG}"
FOLDS_FINAL="$D/cv_folds_final.tsv"
FIXED="${FIXED:-$P/targets/combined/${DRUG}_all.txt}"
ZEROFILL="${ZEROFILL:-$P/targets/amrseq_targets.txt}"
EXTRA="${EXTRA:-$P/labels/amrfinder/by_drug_tsv/${DRUG}.tsv}"
OUTDIR="${OUTDIR:-$P/predictions}"
CFGDIR="${CFGDIR:-$P/scripts/final_configs}"

# ZEROFILL=none: no acquired targets in the list (arm A), so nothing may be absent.
NEED=("$FOLDS_FINAL" "$FIXED")
ZF_ARGS=()
if [ "$ZEROFILL" != none ]; then NEED+=("$ZEROFILL"); ZF_ARGS=(--zero-fill-targets "$ZEROFILL"); fi
EXTRA_ARGS=()
if [ "$EXTRA" != none ]; then NEED+=("$EXTRA"); EXTRA_ARGS=(--extra-feature-tsv "$EXTRA"); fi
for f in "${NEED[@]}"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done
PRED_SUFFIX=""
[ "$ARM" = D ] || PRED_SUFFIX="_arm${ARM}"
TOPK="${TOPK:-$(grep -vc '^#' "$FIXED")}"

echo "[INFO] drug      = $DRUG   boundary $SR_BOUNDARY   top-k $TOPK"
echo "[INFO] folds     = $FOLDS_FINAL ($(( $(wc -l < "$FOLDS_FINAL") - 1 )) rows)"
echo "[INFO] targets   = $FIXED"
if [ "$ZEROFILL" = none ]; then echo "[INFO] zero-fill = none"
else echo "[INFO] zero-fill = $ZEROFILL ($(grep -vc '^#' "$ZEROFILL") acquired targets)"; fi
if [ "$EXTRA" = none ]; then echo "[INFO] presence  = NONE (arm $ARM, embeddings only)"
else echo "[INFO] presence  = $EXTRA ($(( $(head -1 "$EXTRA" | tr '\t' '\n' | wc -l) - 1 )) columns)"; fi
echo "[INFO] embeddings= $EMB_SUBDIR"

declare -A MROOT=( [evo2]="$P/results/ml/pa_evo2"
                   [carbon]="$P/results/ml/pa_carbon"
                   [dnabert]="$P/results/ml/pa_dnabert" )

mkdir -p "$OUTDIR"
for M in ${MODELS:-evo2 carbon dnabert}; do
  echo
  echo "################ $M ################"
  DEV="$D/cache_${M}_${CACHE_TAG}"
  OUT="$D/cache_${M}_${CACHE_TAG}_HOLDOUT"
  # CFGDIR rather than CFG: the loop runs three models, so a single CFG from the
  # environment would be applied to the first and then unset, silently sending
  # the other two to the default path.
  CFG="$CFGDIR/final_config_${DRUG}_${M}_arm${ARM}.json"
  for f in "$DEV/precompute_manifest.json" "$CFG"; do
    [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
  done

  if [ ! -e "$OUT/precompute_manifest.json" ]; then
    $PY "$SCRIPTS/cv_precompute_fold_features.py" \
      --root "${MROOT[$M]}" \
      --embedding-subdir "$EMB_SUBDIR" \
      --folds-tsv "$FOLDS_FINAL" \
      --pca-mode shared4 --pca-dim 4 \
      --var-threshold 0.99 --max-components 25 \
      --top-targets 256 --seed 2026 \
      --fixed-targets "$FIXED" \
      "${ZF_ARGS[@]+"${ZF_ARGS[@]}"}" \
      --only-fold 0:0 \
      --out-dir "$OUT"
  else
    echo "[INFO] reusing existing $OUT"
  fi

  # The holdout cache must match the development cache in everything but folds.
  # The manifest does not record --embedding-subdir (see MTB section 14), so a
  # wrong subdirectory would pass this check -- which is why it is written out
  # above and echoed, not remembered.
  $PY - "$DEV/precompute_manifest.json" "$OUT/precompute_manifest.json" <<'PYEOF'
import json, sys
dev, hold = (json.load(open(p)) for p in sys.argv[1:3])
# n_targets_total is deliberately NOT compared. It counts the targets present in
# the *cohort's* map, and the holdout cache covers a superset of the development
# isolates, so it is >= development's by construction -- ceftazidime has one
# acquired-gene target whose only carriers are held out. What has to match is the
# feature geometry: the fixed target list, its length, the PCA settings and the
# seed. Each fold refits the PCA on its own training rows, so the cached matrices
# are not expected to agree numerically and comparing them would be wrong.
bad = [k for k in ("n_fixed_targets", "pca_mode", "pca_dim", "var_threshold",
                   "max_components", "top_targets", "seed", "fixed_targets",
                   "root")
       if dev.get(k) != hold.get(k)]
if hold.get("n_targets_total", 0) < dev.get("n_targets_total", 0):
    bad.append("n_targets_total")          # fewer than development is impossible
if bad:
    sys.exit("[ERROR] holdout cache differs from development on "
             + ", ".join(f"{k}: {dev.get(k)!r} vs {hold.get(k)!r}" for k in bad))
print(f"[OK] geometry matches development "
      f"({hold['n_fixed_targets']} targets, {hold['n_samples']:,} isolates)")
PYEOF

  $PY "$SCRIPTS/cv_run_architecture.py" \
    --cache-dir "$OUT" \
    "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}" \
    --mode shared4 --top-k "$TOPK" \
    --configs "$CFG" \
    --sr-boundary "$SR_BOUNDARY" \
    --label "${M}_${TAG}" \
    --save-predictions "$OUTDIR/${DRUG}_${M}_holdout_predictions${PRED_SUFFIX}.tsv" \
    --out-dir "$D/runs/${M}_${TAG}" \
    --n-jobs "$NJOBS" --seed 2026
done

echo
echo "[DONE] $D/runs/{evo2,carbon,dnabert}_${TAG}"
echo "       predictions in $OUTDIR/"
