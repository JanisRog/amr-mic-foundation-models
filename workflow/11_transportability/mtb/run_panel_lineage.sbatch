#!/bin/bash
#SBATCH --job-name=panel_lolo
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-9

# Leave-one-lineage-out for the baseline and gwasfiltered arms.
#
# WHY THIS RUN DECIDES SOMETHING. On development folds gwfrule is cost-neutral
# against baseline (median +0.0017, 2 of 8, p=0.46) and better on every ladder
# metric (bal_OT_common +0.0269 7/1 p=0.023, MAE -0.0079 7/1 p=0.023, spearman
# 8/0 p=0.0078). That is the exact signature a co-resistance marker produces: it
# narrows an isolate's MIC distribution without encoding mechanism, so the ladder
# improves and the S/R call does not. §10 recorded a marker doing precisely this
# and giving the gain back across lineages, and §2 measured the ladder as the
# lineage-fragile quantity -- evo2 bal_OT 0.742 random, 0.683 lineage folds,
# 0.599 on held-out L4, against 0.756-0.781 for size-matched random folds.
#
# The degradation §2 measured (-0.098 to -0.133) is four to five times the gain
# being claimed (+0.027). If any part of the gain is lineage-borne this run
# erases it. If gwfrule's advantage survives, the wider target set is carrying
# mechanism the catalogue set does not.
#
# PRE-REGISTERED, before the run:
#   PASS   gwfrule keeps a positive median Δbal_OT_common against baseline on
#          lineage folds, and the across-lineage drop is no worse than
#          baseline's. The arm is then carried forward.
#   FAIL   the Δ collapses toward zero or reverses, or gwfrule drops further
#          across lineages than baseline. The development gain was lineage
#          signal and baseline stays the reported arm.
#   The size-matched control must reproduce the random-fold numbers; if it does
#   not, the lineage folds differ in training size too and nothing is testable.
#
# NO NEW SELECTION. Each (drug, model, arm) reuses the configuration 1-SE
# selection already chose on development folds, read from the configs/ JSON that
# slim_mtb_runs.sh wrote. Re-sweeping 1,000 configs here would be choosing a
# configuration on the test that is meant to judge it.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="$MTB/scripts"
EMB_SUBDIR="${EMB_SUBDIR:-ml_inputs/embeddings/panel_af010}"
LINEAGE="${LINEAGE:-$ROOT/sample_lineage.tsv}"
NJOBS="${NJOBS:-16}"
MODELS="${MODELS:-evo2}"          # evo2 is the reported basis; override to widen
# FOLDS is not configurable: it is derived below from the number of lineages
# actually present in the drug's cohort. A fixed 4 is what packed six groups
# into four folds and left rifampicin and bedaquiline with an empty one.
SIZEMATCH="${SIZEMATCH:-1}"       # also build the random size-matched control

# drug  code  boundary  fold-dir
SPEC=(
  "isoniazid    inh -3.0 cv_isoniazid_final"
  "rifampicin   rif -1.0 cv_rifampicin"
  "ethambutol   emb  2.0 cv_ethambutol"
  # Task 3 fails cleanly until ethionamide's development runs land and
  # slim_mtb_runs.sh writes its selected config -- run_arm refuses rather than
  # inventing one. Resubmit just this task with --array=3 when it is ready.
  "ethionamide  eth  2.0 cv_ethionamide"
  "kanamycin    kan  2.0 cv_kanamycin"
  "amikacin     ami  0.0 cv_amikacin"
  "levofloxacin lev  0.0 cv_levofloxacin"
  "moxifloxacin mxf  0.0 cv_moxifloxacin"
  "bedaquiline  bdq -2.0 cv_bedaquiline"
  "delamanid    dlm -3.0 cv_delamanid"
)
read -r DRUG CODE BOUND FOLDDIR <<< "${SPEC[$SLURM_ARRAY_TASK_ID]}"

R="$ROOT/$FOLDDIR"
L="$R/lineage"
mkdir -p "$L"
echo "[INFO] drug=$DRUG code=$CODE boundary=$BOUND dir=$FOLDDIR"

MICTSV="$ROOT/${DRUG}_log2_mic.tsv"
RULE="$ROOT/rule_bit_${CODE}_tbprofiler.tsv"
for f in "$MICTSV" "$RULE" "$LINEAGE"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done
echo "[INFO] lineage columns: $(head -1 "$LINEAGE")"
# sample_lineage.tsv carries THREE candidate columns -- main_lineage, lineage
# and sub_lineage. main_lineage is the right one: §2's test is four folds (L4,
# L2, L1, L3) and the panel is L4 51.7% / L2 34.8% / L3 7.7% / L1 5.1%, which is
# main lineage. Grouping on sub_lineage would make dozens of groups, pack them
# into four folds, and put relatives from one lineage into different folds --
# destroying the test while still leaving fold sizes unequal enough to pass a
# naive shape check. The census below is printed so the choice is visible, and
# the group count is bounded so a sub-lineage column cannot be used by mistake.
LINEAGE_COL="${LINEAGE_COL:-main_lineage}"
head -1 "$LINEAGE" | tr '\t' '\n' | grep -qxF "$LINEAGE_COL" || {
  echo "[ERROR] no column '$LINEAGE_COL' in $LINEAGE"
  echo "        columns are: $(head -1 "$LINEAGE" | tr '\t' ' ')"
  echo "        set LINEAGE_COL=<name> to pick the grouping column"
  exit 1; }
$PY - "$LINEAGE" "$LINEAGE_COL" <<'PYEOF'
import collections, csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
col = sys.argv[2]
c = collections.Counter((r.get(col) or "").strip() or "(blank)" for r in rows)
print(f"[INFO] grouping column '{col}': {len(c)} distinct value(s) over "
      f"{len(rows):,} isolates")
for k, n in c.most_common(12):
    print(f"[INFO]   {k:<24}{n:>7,}  {n / len(rows):>6.1%}")
if not 2 <= len(c) <= 12:
    sys.exit(f"[ERROR] {len(c)} groups is not a leave-one-lineage-out "
             f"partition. Expected a handful of main lineages; this looks like "
             f"a sub-lineage column. Set LINEAGE_COL to the main-lineage one.")
PYEOF

# Target sets. `gwasfiltered` is a SUPERSET of `both`: catalogue loci | curated
# mechanism loci | unannotated GWAS loci. Both files must exist and the second
# must not be smaller, or the arms are not what they are named.
T_BASE="$ROOT/panel_targets/targets_${DRUG}_both.txt"
T_GWF="$ROOT/panel_targets/targets_${DRUG}_gwasfiltered.txt"
# Isoniazid's baseline predates the panel conventions: it ranks 12 of the 24
# targets in targets_mechfirst.txt, not 12 of the 23 in targets_isoniazid_both.txt.
# The feature-count check cannot tell those apart (both give 49), so the list is
# pinned here and its length checked. Override with T_BASE_INH=<path>.
if [ "$DRUG" = "isoniazid" ]; then
  T_BASE=""
  for c in "${T_BASE_INH:-}" "$ROOT/targets_mechfirst.txt" \
           "$R/targets_mechfirst.txt" "$ROOT/panel_targets/targets_mechfirst.txt"; do
    if [ -n "$c" ] && [ -e "$c" ]; then T_BASE="$c"; break; fi
  done
  [ -n "$T_BASE" ] || {
    echo "[ERROR] targets_mechfirst.txt not found; locate it with"
    echo "        find $ROOT -name targets_mechfirst.txt"
    echo "        and resubmit with T_BASE_INH=<path>"; exit 1; }
  n_mf=$(grep -vc '^#' "$T_BASE")
  [ "$n_mf" -eq 24 ] || { echo "[ERROR] $T_BASE has $n_mf targets, expected 24"; exit 1; }
  echo "[INFO] isoniazid baseline candidate list: $T_BASE (24 targets)"
fi
for f in "$T_BASE" "$T_GWF"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; exit 1; }
done
n_base=$(grep -vc '^#' "$T_BASE"); n_gwf=$(grep -vc '^#' "$T_GWF")
echo "[INFO] targets: both=$n_base gwasfiltered=$n_gwf"
[ "$n_gwf" -ge "$n_base" ] || {
  echo "[ERROR] gwasfiltered ($n_gwf) smaller than both ($n_base) -- not a superset"; exit 1; }

# ---- lineage folds, and the size-matched control -------------------------
# main_lineage carries six values, not four: the four lineages plus `mixed`
# (45) and `unknown` (26). Those two are not lineages, and holding `unknown` out
# is not a test of lineage generalisation -- but leaving them in costs more than
# tidiness. cv_build_folds.py only gives each group its own fold when
# len(units) <= --folds; with six groups and four folds it packs them, which
# merged lineage1 into lineage2's fold for isoniazid and left an empty fold for
# rifampicin and bedaquiline. Restricting to the real lineages restores §2's
# design exactly: one lineage per fold, L4 / L2 / L3 / L1.
KEEP_GROUPS="${KEEP_GROUPS:-lineage1 lineage2 lineage3 lineage4}"
KEPT="$L/lineage_kept.tsv"
$PY - "$LINEAGE" "$LINEAGE_COL" "$KEPT" $KEEP_GROUPS <<'PYEOF'
import collections, csv, sys
src, col, dst, keep = sys.argv[1], sys.argv[2], sys.argv[3], set(sys.argv[4:])
rows = list(csv.DictReader(open(src), delimiter="\t"))
kept = [r for r in rows if (r.get(col) or "").strip() in keep]
seen = collections.Counter((r.get(col) or "").strip() for r in kept)
missing = keep - set(seen)
if missing:
    sys.exit(f"[ERROR] requested group(s) absent from {col}: {sorted(missing)}")
with open(dst, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample_id", col])
    for r in kept:
        w.writerow([r["sample_id"], (r.get(col) or "").strip()])
print(f"[INFO] kept {len(kept):,} of {len(rows):,} isolates in "
      f"{len(seen)} group(s): "
      + ", ".join(f"{k} {n:,}" for k, n in seen.most_common()))
print(f"[INFO] dropped {len(rows) - len(kept):,} without a lineage assignment")
PYEOF
[ -s "$KEPT" ] || { echo "[ERROR] no lineage-assigned isolates"; exit 1; }

# The development cohort as bare sample ids. In cv_folds.tsv fold<0 is the
# holdout and fold>=0 the development folds, so this is the set the arms were
# selected on -- the lineage folds must cover exactly it and nothing else.
DEV="$L/dev_samples.txt"
awk -F'\t' '
  NR==1 { for (i = 1; i <= NF; i++) { if ($i == "sample_id") si = i
                                      if ($i == "fold") fi = i }
          if (!si || !fi) { print "[ERROR] cv_folds.tsv lacks sample_id/fold" \
                                  > "/dev/stderr"; exit 1 }
          next }
  $fi >= 0 { print $si }' "$R/cv_folds.tsv" | sort -u > "$DEV.all"
# Only isolates with a real lineage can enter a leave-one-lineage-out fold.
awk -F'\t' 'NR>1 {print $1}' "$KEPT" | sort -u > "$DEV.lin"
comm -12 "$DEV.all" "$DEV.lin" > "$DEV"
n_all=$(wc -l < "$DEV.all"); n_dev=$(wc -l < "$DEV")
echo "[INFO] development cohort: $n_dev of $n_all isolates have a lineage -> $DEV"
[ "$n_dev" -gt 0 ] || { echo "[ERROR] no development isolates extracted"; exit 1; }

# One fold per lineage present in this drug's cohort. Hardcoding 4 would pack
# or leave empty folds whenever a drug is missing a lineage, which is exactly
# how rifampicin and bedaquiline failed.
FOLDS=$($PY - "$KEPT" "$LINEAGE_COL" "$DEV" <<'PYEOF'
import csv, sys
keep = {l.strip() for l in open(sys.argv[3]) if l.strip()}
g = {(r.get(sys.argv[2]) or "").strip() for r in
     csv.DictReader(open(sys.argv[1]), delimiter="\t")
     if r["sample_id"] in keep}
print(len(g))
PYEOF
)
echo "[INFO] $FOLDS lineage(s) present in this cohort -> $FOLDS folds"
[ "$FOLDS" -ge 3 ] || { echo "[ERROR] only $FOLDS lineage(s); not a partition"; exit 1; }

if [ ! -e "$L/cv_folds.tsv" ]; then
  echo "=== building $FOLDS lineage folds"
  $PY "$SCRIPTS/cv_build_folds.py" \
    --mic-from "$MICTSV" --restrict-samples "$DEV" \
    --allele-tsv "${KEPT}#${LINEAGE_COL}" \
    --holdout-frac 0.0 --folds "$FOLDS" --repeats 1 --seed 2026 \
    --out-dir "$L"
else
  echo "[INFO] reusing $L/cv_folds.tsv"
fi
# The grouping must actually have produced lineage-shaped folds. Four roughly
# equal folds mean the column was read but carried one value per isolate, which
# would be a random split wearing a lineage label.
$PY - "$L/cv_folds.tsv" <<'PYEOF'
import collections, csv, sys
c = collections.Counter(r["fold"] for r in
                        csv.DictReader(open(sys.argv[1]), delimiter="\t")
                        if int(r["fold"]) >= 0)
sizes = sorted(c.values(), reverse=True)
print(f"[INFO] lineage fold sizes {sizes}")
if len(sizes) < 2 or sizes[0] / sizes[-1] < 3:
    sys.exit(f"[ERROR] folds are nearly equal (max/min "
             f"{sizes[0] / sizes[-1]:.1f}); the grouping did not take")
PYEOF
SIZES=$($PY - "$L/cv_folds.tsv" "$FOLDS" <<'PYEOF'
import collections, csv, sys
want = int(sys.argv[2])
c = collections.Counter(int(r["fold"]) for r in
                        csv.DictReader(open(sys.argv[1]), delimiter="\t")
                        if int(r["fold"]) >= 0)
# Every fold must be non-empty and counted. Dropping empty folds here while
# --folds stays at FOLDS is what made cv_build_folds.py reject --fold-sizes
# with "gives 3 sizes for 4 folds".
sizes = [c.get(i, 0) for i in range(want)]
if any(n == 0 for n in sizes):
    sys.exit(f"[ERROR] fold(s) {[i for i, n in enumerate(sizes) if not n]} are "
             f"empty; the grouping did not produce one lineage per fold")
print(",".join(str(n) for n in sorted(sizes, reverse=True)))
PYEOF
)
echo "[INFO] lineage fold sizes: $SIZES"

if [ "$SIZEMATCH" = "1" ] && [ ! -e "$L/sizematched/cv_folds.tsv" ]; then
  echo "=== building size-matched random folds ($SIZES)"
  $PY "$SCRIPTS/cv_build_folds.py" \
    --mic-from "$MICTSV" --restrict-samples "$DEV" \
    --holdout-frac 0.0 --folds "$FOLDS" --repeats 1 --seed 2026 \
    --fold-sizes "$SIZES" --out-dir "$L/sizematched"
fi

declare -A MROOT=( [evo2]="$MTB/results/ml/03_mtb_evo2"
                   [carbon]="$MTB/results/ml/03_mtb"
                   [dnabert]="$MTB/results/ml/03_mtb_dnabert" )

# arm -> target file, extra-feature file, the suffix slim wrote its config under
run_arm () {
  local M="$1" ARM="$2" FIXED="$3" EXTRA="$4" SUFFIXES="$5" FOLDSDIR="$6" LABEL="$7"

  # The already-selected configuration, not a new sweep.
  # `[ -e ] && { ...; break 2; }` as the last command in a loop body returns
  # non-zero when the file is absent, and `set -e` then aborts the whole job on
  # what is only a miss in the search. if/fi keeps the test out of $?.
  local CFG="" SUF=""
  for s in $SUFFIXES; do
    for c in "$ROOT/slim_arms/configs/${FOLDDIR#cv_}__${M}__${s}.json" \
             "$ROOT/slim_arms/configs/${DRUG}__${M}__${s}.json"; do
      if [ -e "$c" ]; then CFG="$c"; SUF="$s"; break 2; fi
    done
  done
  [ -n "$CFG" ] || { echo "[ERROR] no selected config for $M/$ARM (tried $SUFFIXES)"; return 1; }

  # top-k is part of the model, not a run setting. The development runs used a
  # different value per drug and arm -- isoniazid's baseline ranks 12 of 23
  # targets within each fold, every other arm takes its whole target set, so
  # gwasfiltered runs at 35-72 targets. A fixed --top-k 8 here silently cut every
  # arm to the same 8 targets, so the first lineage run compared two 8-target
  # subsets rather than two target sets, and the moxifloxacin control (4 targets,
  # under the cap) could not see it. The value is now recovered from the
  # development run's own feature count and the result is checked against it.
  local DEVF="$ROOT/slim_arms/$FOLDDIR/runs/${M}_${SUF}/folds.tsv"
  [ -e "$DEVF" ] || { echo "[ERROR] no development folds at $DEVF"; return 1; }
  local DEV_NF N_EXTRA TOPK_ARM
  DEV_NF=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="n_features")c=i;next}
                       NR==2{print $c; exit}' "$DEVF")
  N_EXTRA=0
  [ -n "$EXTRA" ] && N_EXTRA=$(( $(head -1 "$EXTRA" | awk -F'\t' '{print NF}') - 1 ))
  [ $(( (DEV_NF - N_EXTRA) % 4 )) -eq 0 ] || {
    echo "[ERROR] $M/$ARM: development n_features $DEV_NF minus $N_EXTRA extra"
    echo "        column(s) is not a multiple of 4 PCs; cannot recover top-k"; return 1; }
  TOPK_ARM=$(( (DEV_NF - N_EXTRA) / 4 ))
  echo "    development n_features=$DEV_NF  extra=$N_EXTRA  -> top-k $TOPK_ARM"
  echo "--- $M $ARM  config=$CFG  tag=$($PY -c \
      "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$CFG")"

  local CACHE="$FOLDSDIR/cache_${M}_${ARM}"
  if [ ! -e "$CACHE/precompute_manifest.json" ]; then
    $PY "$SCRIPTS/cv_precompute_fold_features.py" \
      --root "${MROOT[$M]}" --embedding-subdir "$EMB_SUBDIR" \
      --folds-tsv "$FOLDSDIR/cv_folds.tsv" \
      --pca-mode shared4 --pca-dim 4 \
      --var-threshold 0.99 --max-components 25 \
      --top-targets 256 --seed 2026 \
      --fixed-targets "$FIXED" --out-dir "$CACHE"
  fi
  $PY "$SCRIPTS/cv_run_architecture.py" \
    --cache-dir "$CACHE" ${EXTRA:+--extra-feature-tsv "$EXTRA"} \
    --mode shared4 --top-k "$TOPK_ARM" \
    --configs "$CFG" --sr-boundary "$BOUND" \
    --label "${M}_${LABEL}" \
    --out-dir "$FOLDSDIR/runs/${M}_${LABEL}" \
    --n-jobs "$NJOBS" --seed 2026

  # The refit must be the same model as the one selected, feature for feature.
  local NEW_NF
  NEW_NF=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="n_features")c=i;next}
                       NR==2{print $c; exit}' "$FOLDSDIR/runs/${M}_${LABEL}/folds.tsv")
  [ "$NEW_NF" = "$DEV_NF" ] || {
    echo "[ERROR] $M/$ARM: refit has $NEW_NF features, development had $DEV_NF"
    return 1; }
  echo "    n_features $NEW_NF matches development"
}

for M in $MODELS; do
  run_arm "$M" baseline "$T_BASE" "$RULE" \
          "both_final1k tbrule1k final1k" "$L" lolo_baseline
  run_arm "$M" gwfrule  "$T_GWF"  "$RULE" \
          "gwasfiltered_gwfrule1k"        "$L" lolo_gwfrule
  if [ "$SIZEMATCH" = "1" ]; then
    run_arm "$M" baseline "$T_BASE" "$RULE" \
            "both_final1k tbrule1k final1k" "$L/sizematched" sm_baseline
    run_arm "$M" gwfrule  "$T_GWF"  "$RULE" \
            "gwasfiltered_gwfrule1k"        "$L/sizematched" sm_gwfrule
  fi
done

echo
echo "[DONE] $DRUG -> $L/runs and $L/sizematched/runs"
