#!/bin/bash
#SBATCH --job-name=pa_cohorts
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-9

# Labels, folds and the presence block for one P. aeruginosa drug -- the three
# steps that stand between the source metadata and a CV sweep.
#
# These were run by hand for the first imputed rebuild and written down only as a
# numbered list afterwards, which is exactly the failure this file removes: the
# list said *what* was run, not with which flags, against which paths, in which
# order. A relabelling now costs one submission instead of an afternoon of
# reconstruction.
#
# The rule being applied here (section 1): clamp to [LO, HI]; a STRICT bound
# ('>', '<') moves one dilution outward because it excludes the value it names;
# an INCLUSIVE bound ('>=', '<=') is transcribed as-is because it admits the
# value it names, and moving it would record a dilution the lab never claimed.
# Rebelo et al. excluded. Nothing is dropped.
#
# Order within a drug is forced: the presence block's carrier floor is counted on
# the development isolates, so it cannot be built until the folds exist.

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
PASCRIPTS="${PASCRIPTS:-$P/scripts}"

CENSORING="${CENSORING:-impute}"
LABELS="${LABELS:-$P/labels_imputed_strict}"
CVROOT="${CVROOT:-$P/imputed_strict}"
MINCAR="${MINCAR:-10}"
SEED="${SEED:-2026}"

METAROOT="${METAROOT:-$P/models/pseaer_all_antibiotics_amrfinder_groupedpca4_rf500_cv_package/data/metadata}"
AMRDIR="${AMRDIR:-$P/labels/amrfinder}"
EMBEDDED="${EMBEDDED:-$P/embedded_samples.txt}"
# The panel this one is measured against, for the regression check below. Empty
# skips the check.
PREV="${PREV:-$P/imputed}"
PREVLABELS="${PREVLABELS:-$P/labels_imputed}"

DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)
DRUG="${DRUGS[$SLURM_ARRAY_TASK_ID]}"

META="$METAROOT/$DRUG/metadata.tsv"
COLS="$AMRDIR/amr_columns_by_drug/${DRUG}.txt"
PRESENCE="$AMRDIR/amrfinder_presence.tsv"

# Every input named before anything is written. Four submissions were lost on the
# MTB panel to a path that did not exist, discovered one stage in.
fail=0
for f in "$META" "$COLS" "$PRESENCE" "$EMBEDDED" \
         "$PASCRIPTS/prepare_pseaer_drug_inputs.py" \
         "$PASCRIPTS/make_drug_presence_tsv.py" \
         "$SCRIPTS/cv_build_folds.py"; do
  [ -e "$f" ] || { echo "[ERROR] missing $f"; fail=1; }
done
[ "$fail" -eq 0 ] || { echo "[ERROR] nothing was written"; exit 1; }

echo "[INFO] drug      = $DRUG"
echo "[INFO] censoring = $CENSORING"
echo "[INFO] labels    -> $LABELS"
echo "[INFO] folds     -> $CVROOT/cv_$DRUG"

mkdir -p "$LABELS" "$LABELS/amrfinder/by_drug_tsv" "$CVROOT/cv_$DRUG"

# ---- 1. labels -------------------------------------------------------------
echo
echo "################ labels"
$PY "$PASCRIPTS/prepare_pseaer_drug_inputs.py" \
  --metadata "$META" \
  --drug "$DRUG" \
  --censoring "$CENSORING" \
  --out-dir "$LABELS"

MICTSV="$LABELS/${DRUG}_log2_mic.tsv"
[ -e "$MICTSV" ] || { echo "[ERROR] no $MICTSV"; exit 1; }

# ---- 2. folds --------------------------------------------------------------
echo
echo "################ folds"
D="$CVROOT/cv_$DRUG"
if [ ! -e "$D/cv_folds.tsv" ]; then
  $PY "$SCRIPTS/cv_build_folds.py" \
    --mic-from "$MICTSV" \
    --restrict-samples "$EMBEDDED" \
    --holdout-frac 0.2 --folds 5 --repeats 3 --seed "$SEED" \
    --out-dir "$D"
else
  echo "[INFO] reusing $D/cv_folds.tsv"
fi
for f in "$D/cv_folds.tsv" "$D/cv_folds_final.tsv"; do
  [ -e "$f" ] || { echo "[ERROR] fold builder did not write $f"; exit 1; }
done

# The regression check, free and immediate. Three drugs (cefepime, imipenem,
# levofloxacin) carry no inclusive bounds at all, so their labels are unchanged
# by the rule and their folds must come out byte-identical. If they do not, the
# fold builder is not deterministic under this seed and every cross-panel
# comparison below is meaningless -- better to learn that here than after a
# 24-hour sweep.
if [ -n "$PREV" ] && [ -e "$PREV/cv_$DRUG/cv_folds.tsv" ]; then
  if cmp -s "$D/cv_folds.tsv" "$PREV/cv_$DRUG/cv_folds.tsv"; then
    echo "[INFO] folds identical to $PREV/cv_$DRUG -- labels unchanged by the rule"
  else
    n=$($PY - "$MICTSV" "$PREVLABELS/${DRUG}_log2_mic.tsv" <<'PYEOF'
import csv, sys
def load(p):
    return {r["sample_id"]: r["log2_mic"]
            for r in csv.DictReader(open(p), delimiter="\t")}
try:
    a, b = load(sys.argv[1]), load(sys.argv[2])
except OSError:
    print("?"); raise SystemExit
print(sum(1 for k in a if k in b and a[k] != b[k]))
PYEOF
)
    echo "[INFO] folds differ from $PREV/cv_$DRUG ($n label(s) moved)"
  fi
fi

# ---- 3. presence block -----------------------------------------------------
# --samples is development only, so the carrier floor never sees a holdout
# isolate; --emit-samples is the full cohort, so holdout scoring has rows for
# every isolate it is asked about. Widening --samples instead would leak the
# holdout into feature selection; leaving --emit-samples off aborts the holdout.
echo
echo "################ presence"
$PY "$PASCRIPTS/make_drug_presence_tsv.py" \
  --presence "$PRESENCE" \
  --columns "$COLS" \
  --samples "$D/cv_folds.tsv" \
  --emit-samples "$D/cv_folds_final.tsv" \
  --min-carriers "$MINCAR" \
  --out "$LABELS/amrfinder/by_drug_tsv/${DRUG}.tsv"

echo
echo "[DONE] $DRUG: $MICTSV, $D/cv_folds{,_final}.tsv, presence block"
