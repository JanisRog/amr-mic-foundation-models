#!/usr/bin/env bash
#SBATCH --job-name=mtb_calls_af
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=8:00:00

# Every isolate's catalogue S/R call at each allele fraction.
#
# It reads 10,170 raw VCFs end to end -- the same pass sweep_consensus_threshold
# makes -- which is a login-node job only if you want to be unpopular. Memory is
# small: one VCF at a time plus a call matrix of samples x drugs.
#
# The output is not only diagnostic. The rule bit in the model's feature set is
# one of these calls, currently taken from TB-Profiler at its default 0.1 while
# the embeddings come from a ~0.9 consensus, so emitting the call per threshold
# is what lets the feature blocks agree on what an isolate carries.

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PYTHON_BIN="${PYTHON_BIN:-${PY:?set PY in config/paths.env}}"
SCRIPT="${SCRIPT:-$MTB/scripts/tbprofiler_calls_by_threshold.py}"

SNIPPY_RUNS="${SNIPPY_RUNS:-$MTB/results/gwas/snippy_h37rv/runs}"
: "${CATALOGUE:?set CATALOGUE to call_sources_by_variant.tsv}"
: "${OUT_PREFIX:?set OUT_PREFIX}"
MIC_TSV="${MIC_TSV:-}"
MIC_SUFFIX="${MIC_SUFFIX:-_int}"
THRESHOLDS="${THRESHOLDS:-0.05,0.1,0.2,0.5,0.9}"
REFERENCE_THRESHOLD="${REFERENCE_THRESHOLD:-0.9}"
SR_BOUNDARY="${SR_BOUNDARY:-INH=-3.0}"
MIN_QUAL="${MIN_QUAL:-100}"
MIN_DP="${MIN_DP:-10}"
DRUGS="${DRUGS:-}"
SCORE_SAMPLES="${SCORE_SAMPLES:-}"

[[ -s "$SCRIPT" ]]    || { echo "[ERROR] missing $SCRIPT" >&2; exit 2; }
[[ -s "$CATALOGUE" ]] || { echo "[ERROR] missing $CATALOGUE" >&2; exit 2; }
[[ -d "$SNIPPY_RUNS" ]] || { echo "[ERROR] missing $SNIPPY_RUNS" >&2; exit 2; }
[[ -z "$MIC_TSV" || -s "$MIC_TSV" ]] || { echo "[ERROR] missing $MIC_TSV" >&2; exit 2; }

export PYTHONUNBUFFERED=1

echo "[INFO] runs       = $SNIPPY_RUNS"
echo "[INFO] catalogue  = $CATALOGUE"
echo "[INFO] thresholds = $THRESHOLDS (reference $REFERENCE_THRESHOLD)"
echo "[INFO] floors     = QUAL>=$MIN_QUAL DP>=$MIN_DP"
echo "[INFO] out        = ${OUT_PREFIX}_af*.tsv"
[[ -n "$MIC_TSV" ]] && echo "[INFO] mic        = $MIC_TSV ($SR_BOUNDARY)"
[[ -n "$DRUGS"   ]] && echo "[INFO] drugs      = $DRUGS"
[[ -n "$SCORE_SAMPLES" ]] && echo "[INFO] scored on  = $SCORE_SAMPLES"

EXTRA=()
[[ -n "$MIC_TSV" ]] && EXTRA+=(--mic-tsv "$MIC_TSV" --mic-suffix "$MIC_SUFFIX")
[[ -n "$DRUGS"   ]] && EXTRA+=(--drugs "$DRUGS")
[[ -n "$SCORE_SAMPLES" ]] && EXTRA+=(--score-samples "$SCORE_SAMPLES")

"$PYTHON_BIN" "$SCRIPT" \
  --snippy-runs "$SNIPPY_RUNS" \
  --catalogue "$CATALOGUE" \
  --thresholds "$THRESHOLDS" \
  --reference-threshold "$REFERENCE_THRESHOLD" \
  --sr-boundary "$SR_BOUNDARY" \
  --min-qual "$MIN_QUAL" \
  --min-dp "$MIN_DP" \
  "${EXTRA[@]}" \
  --out-prefix "$OUT_PREFIX"

echo "[OK] ${OUT_PREFIX}_flips.tsv"
