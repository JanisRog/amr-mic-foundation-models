#!/usr/bin/env bash
# Select each arm's config by the declared rule, then keep only that config's rows.
#
# A sweep directory is 1,000 configs x 15 folds x 2 aggregations; the comparison
# uses one config per arm per drug per model. Shipping the rest is gigabytes of
# which ~0.1% is read, which is why the PA panel has the same step.
#
# Selection is delegated to select_final_config.py rather than reimplemented. The
# rule -- minimise cost = 2*vMER + MER, break ties within 1 SE on bal_OT, median
# aggregation, rf family -- has to be identical to the one every other reported
# number used, and a second copy of it is a second place for it to drift.
#
# The config column is located BY NAME. cv_run_architecture writes the metrics
# first and appends config/aggregation/fold afterwards, so $1 is `exact`, not
# `config`; filtering on position silently keeps nothing but the header, which is
# what it did the first time on the PA panel.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
SCRIPTS="${SCRIPTS:-$MTB/scripts}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb/metadata}"
CONFIGS="${CONFIGS:-$SCRIPTS/mtb_rf1000_configs.json}"
OUT="${OUT:-$ROOT/slim_arms}"
CFGOUT="${CFGOUT:-$OUT/configs}"

# Run-directory suffixes to collect. Three naming conventions are in play,
# because two drugs have bespoke dev scripts that take no target-set variable:
#
#   run_panel_dev.sbatch   <model>_<targetset>_<tag>   both_final1k, both_arm75
#   run_inh_redo.sbatch    <model>_<tag>               tbrule1k  (isoniazid)
#   run_rif_dev.sbatch     <model>_<tag>               final1k   (rifampicin)
#
# All of them are listed so the glob finds whichever exists for a given drug.
# Dropping the bare `final1k` silently costs rifampicin its baseline arm and
# drops the panel to n=8 without any error being raised.
RUNS="${RUNS:-both_final1k tbrule1k final1k both_arm75 arm75 gwasfiltered_gwfrule1k gwasfiltered_gwfsnp1k}"

DRUGS=(isoniazid_final rifampicin ethambutol ethionamide kanamycin amikacin
       levofloxacin moxifloxacin bedaquiline delamanid)
MODELS=(evo2 carbon dnabert)

[ -e "$CONFIGS" ] || { echo "[ERROR] missing $CONFIGS"; exit 1; }
rm -rf "$OUT"; mkdir -p "$OUT" "$CFGOUT"

n_ok=0; n_miss=0
for d in "${DRUGS[@]}"; do
  for m in "${MODELS[@]}"; do
    for suffix in $RUNS; do
      src="$ROOT/cv_$d/runs/${m}_${suffix}"
      [ -e "$src/summary.tsv" ] || { n_miss=$((n_miss + 1)); continue; }

      cfg="$CFGOUT/${d}__${m}__${suffix}.json"
      if ! $PY "$SCRIPTS/select_final_config.py" \
             --run-dir "$src" --configs "$CONFIGS" \
             --aggregation median --family rf \
             --w-vmer 2.0 --w-mer 1.0 --se-multiple 1.0 \
             --objective cost --tie-break bal_OT --tie-metric bal_OT \
             --out "$cfg" > "$CFGOUT/${d}__${m}__${suffix}.log" 2>&1; then
        echo "[WARN] selection failed for $d/$m/$suffix"
        continue
      fi
      tag=$($PY -c "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$cfg")

      dst="$OUT/cv_$d/runs/${m}_${suffix}"; mkdir -p "$dst"
      for f in folds.tsv per_class.tsv; do
        [ -e "$src/$f" ] || continue
        awk -F'\t' -v t="$tag" '
          NR==1 { for (i = 1; i <= NF; i++) if ($i == "config") c = i
                  if (!c) { print "[ERROR] no config column in " FILENAME > "/dev/stderr"; exit 1 }
                  print; next }
          $c == t' "$src/$f" > "$dst/$f"
        # A header-only result means the tag was not found; that is the
        # position-instead-of-name failure, and it must not pass silently.
        [ "$(wc -l < "$dst/$f")" -gt 1 ] || {
          echo "[ERROR] $d/$m/$suffix: $f kept no rows for config '$tag'"; exit 1; }
      done
      n_ok=$((n_ok + 1))
    done
  done
done

cp -r "$CFGOUT" "$OUT/" 2>/dev/null || true

# Per-drug arm census. A drug that silently arrives with fewer arms than the
# others shrinks every paired comparison downstream without raising anything --
# this is exactly how rifampicin's `final1k` baseline went missing and reduced
# the panel to n=8. Print the census and name the gaps.
echo
echo "=== arms collected per drug ==="
gaps=0
for d in "${DRUGS[@]}"; do
  found=""
  for suffix in $RUNS; do
    for m in "${MODELS[@]}"; do
      [ -e "$OUT/cv_$d/runs/${m}_${suffix}/folds.tsv" ] && { found="$found $suffix"; break; }
    done
  done
  n=$(echo $found | wc -w | tr -d ' ')
  printf '%-20s %s\n' "$d" "${found:- (none)}"
  [ "$n" -ge 4 ] || gaps=$((gaps + 1))
done
[ "$gaps" -eq 0 ] || echo "[WARN] $gaps drug(s) carry fewer than 4 arms -- paired \
comparisons will drop those drugs. Check the run-directory suffix before reporting an n."

echo
echo "[DONE] $n_ok run(s) slimmed into $OUT  ($n_miss absent)"
du -sh "$OUT"
