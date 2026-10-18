#!/usr/bin/env bash
# Filter each run's folds.tsv / per_class.tsv down to the chosen config.
#
# The workbook reads only the chosen config's rows, so shipping 1,000 configs x
# 15 folds x 3 models x 10 drugs is ~95 MB of which ~0.1% is used. This keeps the
# header and the matching rows.
#
# The config column is located **by name**, not by position: cv_run_architecture
# writes the metrics first and appends config/aggregation/fold afterwards, so
# `$1` is `exact`, not `config`. Filtering on $1 silently keeps nothing but the
# header -- which is exactly what it did the first time, and is the same
# position-instead-of-name mistake this project has made twice before.

set -euo pipefail

P="${P:-$PWD}"
CVROOT="${CVROOT:-$P/imputed}"
CFGDIR="${CFGDIR:-$P/scripts/final_configs_imputed}"
OUT="${OUT:-$P/slim}"
PY="${PY:?set PY in config/paths.env}"
# The run directories to slim. They must match the CFGDIR being used: running
# with the ladder configs against the cost-selected holdout looks for a tag that
# is genuinely not in that file, and the script stops -- correctly, but for a
# reason that looks like corrupt data rather than a mismatched argument.
# Each run takes the config of its own arm (armC1k / armC_HOLDOUT -> *_armC.json).
RUNS="${RUNS:-armD1k armD_HOLDOUT armC1k armC_HOLDOUT}"
BASE_RUNS="${BASE_RUNS:-base_presence base_presence_HOLDOUT varonly varonly_HOLDOUT}"
DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)

rm -rf "$OUT"; mkdir -p "$OUT"
for d in "${DRUGS[@]}"; do
  for m in evo2 carbon dnabert; do
    for suffix in $RUNS; do
      run="${m}_${suffix}"
      src="$CVROOT/cv_$d/runs/$run"; dst="$OUT/cv_$d/runs/$run"
      [ -d "$src" ] || continue
      arm=$(printf '%s' "$suffix" | sed -E 's/^(arm[A-Z]).*/\1/')
      cfg="$CFGDIR/final_config_${d}_${m}_${arm}.json"
      [ -e "$cfg" ] || { echo "[WARN] no config for $d/$m/$arm"; continue; }
      tag=$($PY -c "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$cfg")
      mkdir -p "$dst"
      for f in folds.tsv per_class.tsv; do
        [ -e "$src/$f" ] || continue
        awk -F'\t' -v t="$tag" '
          NR==1 { for (i = 1; i <= NF; i++) if ($i == "config") c = i
                  if (!c) { print "[ERROR] no config column in " FILENAME > "/dev/stderr"; exit 1 }
                  print; next }
          $c == t' "$src/$f" > "$dst/$f"
        n=$(( $(wc -l < "$dst/$f") - 1 ))
        [ "$n" -gt 0 ] || { echo "[ERROR] $d/$m $f: no rows matched $tag"; exit 1; }
      done
      cp "$src/summary.tsv" "$dst/" 2>/dev/null || true
    done
  done
  # Model-free baselines, one run per drug: base_presence (run_pa_presence_baseline
  # .sbatch, config final_config_<drug>_presence.json) and varonly
  # (run_variant_table_baseline.sbatch, final_config_<drug>_varonly.json).
  for run in $BASE_RUNS; do
    src="$CVROOT/cv_$d/runs/$run"; dst="$OUT/cv_$d/runs/$run"
    [ -d "$src" ] || continue
    cname="${run%_HOLDOUT}"; cname="${cname#base_}"
    cfg="$CFGDIR/final_config_${d}_${cname}.json"
    [ -e "$cfg" ] || { echo "[WARN] no config for $d/$run"; continue; }
    tag=$($PY -c "import json,sys;print(json.load(open(sys.argv[1]))[0]['tag'])" "$cfg")
    mkdir -p "$dst"
    for f in folds.tsv per_class.tsv; do
      [ -e "$src/$f" ] || continue
      awk -F'\t' -v t="$tag" '
        NR==1 { for (i = 1; i <= NF; i++) if ($i == "config") c = i
                if (!c) { print "[ERROR] no config column in " FILENAME > "/dev/stderr"; exit 1 }
                print; next }
        $c == t' "$src/$f" > "$dst/$f"
      n=$(( $(wc -l < "$dst/$f") - 1 ))
      [ "$n" -gt 0 ] || { echo "[ERROR] $d/$run $f: no rows matched $tag"; exit 1; }
    done
    cp "$src/summary.tsv" "$dst/" 2>/dev/null || true
  done
  cp "$CVROOT/cv_$d/cv_folds_final.tsv" "$OUT/cv_$d/" 2>/dev/null || true
  cp "$CVROOT/cv_$d/baseline_rows.tsv"  "$OUT/cv_$d/" 2>/dev/null || true
done
echo "[OK] $OUT — $(find "$OUT" -name folds.tsv | wc -l | tr -d ' ') run(s)"
