#!/bin/bash
#SBATCH --job-name=pa_select
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00

# Pick one config per drug per model from the development sweep, by the panel's
# declared rule and nothing else:
#
#   minimise cost = 2*vMER + MER
#   among configs within 1 SE of the minimum, take the highest bal_OT
#   median aggregation, rf family
#
# Those are select_final_config.py's defaults, passed explicitly anyway, because
# a selection rule that lives in a tool's defaults is a selection rule nobody can
# audit from the run log.
#
# --objective bal_OT is a SEPARATE declared analysis, not an alternative headline
# for this one, and is deliberately not run here. Choosing between two objectives
# after seeing both holdouts is the max-over-configs bias arriving at the last
# possible step.

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
MTB="${MTB:?set MTB in config/paths.env}"
PY="${PY:?set PY in config/paths.env}"
TOOLS="${TOOLS:-$MTB/scripts}"
CVROOT="${CVROOT:-$P/imputed_strict}"
CFGDIR="${CFGDIR:-$P/scripts/final_configs_strict}"
CONFIGS="${CONFIGS:-$MTB/scripts/mtb_rf1000_configs.json}"
# ARM names the config files (final_config_<drug>_<model>_arm<ARM>.json).
ARM="${ARM:-D}"
TAG="${TAG:-arm${ARM}1k}"

DRUGS=(amikacin aztreonam cefepime ceftazidime ciprofloxacin imipenem
       levofloxacin meropenem piperacillin-tazobactam tobramycin)
MODELS=(evo2 carbon dnabert)

[ -e "$CONFIGS" ] || { echo "[ERROR] missing $CONFIGS"; exit 1; }
[ -e "$TOOLS/select_final_config.py" ] || {
  echo "[ERROR] missing $TOOLS/select_final_config.py"; exit 1; }

# Two copies of this tool live on the cluster and they are not the same age.
# --objective was added when the ladder analysis was run and only one copy got
# it, so the stale one dies on the first drug with an argparse usage dump that
# does not say which file it came from.
$PY "$TOOLS/select_final_config.py" --help 2>&1 | grep -q -- "--objective" || {
  echo "[ERROR] $TOOLS/select_final_config.py predates --objective."
  echo "        scp the version in git over it, or point TOOLS at a copy that has it."
  exit 1; }

mkdir -p "$CFGDIR"

for DRUG in "${DRUGS[@]}"; do
  for M in "${MODELS[@]}"; do
    RUN="$CVROOT/cv_${DRUG}/runs/${M}_${TAG}"
    [ -e "$RUN/summary.tsv" ] || { echo "[ERROR] missing $RUN/summary.tsv"; exit 1; }
    echo
    echo "################ $DRUG / $M ################"
    $PY "$TOOLS/select_final_config.py" \
      --run-dir "$RUN" \
      --configs "$CONFIGS" \
      --aggregation median --family rf \
      --w-vmer 2.0 --w-mer 1.0 --se-multiple 1.0 \
      --objective cost --tie-break bal_OT --tie-metric bal_OT \
      --out "$CFGDIR/final_config_${DRUG}_${M}_arm${ARM}.json"
  done
done

n=$(ls -1 "$CFGDIR"/final_config_*_arm${ARM}.json 2>/dev/null | wc -l)
echo
echo "[DONE] $n config(s) in $CFGDIR"
[ "$n" -eq 30 ] || { echo "[ERROR] expected 30"; exit 1; }
