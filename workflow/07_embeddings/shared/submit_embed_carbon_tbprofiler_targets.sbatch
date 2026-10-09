#!/usr/bin/env bash
#SBATCH --job-name=mtb_carbon_embed
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PIPE="${PIPE:?set PIPE in config/paths.env}"
PYTHON_BIN="${PYTHON_BIN:-${PY:?set PY in config/paths.env}}"
# The three that decide WHAT gets embedded and WHERE it lands use ${VAR-default},
# not ${VAR:-default}: an exported-but-empty value (a placeholder that never got
# filled in) must fail the checks below, not silently fall back to the defaults
# and embed the wrong fasta into the right-looking output directory.
CARBON_DIR="${CARBON_DIR-${CARBON_MODEL_DIR:?set CARBON_MODEL_DIR in config/paths.env}}"
IN_FASTA="${IN_FASTA-$MTB/results/gwas/sample_windows/tbprofiler_all_targets_gene_up300/gather/unique_window_sequences.fasta}"
OUTROOT="${OUTROOT-$MTB/results/gwas/carbon_embeddings/carbon_8b/tbprofiler_all_targets_gene_up300}"
BATCH_SIZE="${BATCH_SIZE:-32}"
POOLING="${POOLING:-mean}"
LAYER="${LAYER:--1}"
AMBIGUOUS_POLICY="${AMBIGUOUS_POLICY:-keep}"
MAX_RECORDS="${MAX_RECORDS:-0}"
DNA_TAG="${DNA_TAG:-both}"
CONFIG_OVERRIDE="${CONFIG_OVERRIDE:-}"
DTYPE="${DTYPE:-bf16}"
MAX_LENGTH="${MAX_LENGTH:-0}"   # 0 = infer from model config
# where this job's logs live; set it per model so carbon/evo2/dnabert2 stay apart
LOG_DIR="${LOG_DIR:-$MTB/00_logs/carbon}"

[[ -n "$IN_FASTA" ]]  || { echo "[ERROR] IN_FASTA is empty" >&2; exit 2; }
[[ -n "$OUTROOT" ]]   || { echo "[ERROR] OUTROOT is empty" >&2; exit 2; }
[[ -s "$IN_FASTA" ]] || { echo "[ERROR] missing input fasta: $IN_FASTA" >&2; exit 2; }
[[ -s "$CARBON_DIR/config.json" ]] || { echo "[ERROR] missing model: $CARBON_DIR" >&2; exit 2; }

mkdir -p "$OUTROOT" "$LOG_DIR"
STAMP="$(date -u +%Y%m%d_%H%M%S)"
RUN_OUT="$OUTROOT/run_${STAMP}"
mkdir -p "$RUN_OUT"

echo "[INFO] input=$IN_FASTA"
echo "[INFO] sequences=$(grep -c '^>' "$IN_FASTA")"
echo "[INFO] model=$CARBON_DIR"
echo "[INFO] out=$RUN_OUT"

"$PYTHON_BIN" "$PIPE/scripts/09c_embed_with_carbon_hidden.py" \
  --input-fasta "$IN_FASTA" \
  --outdir "$RUN_OUT" \
  --model-dir "$CARBON_DIR" \
  --device auto \
  --dtype "$DTYPE" \
  --max-length "$MAX_LENGTH" \
  --batch-size "$BATCH_SIZE" \
  --pooling "$POOLING" \
  --layer "$LAYER" \
  --ambiguous-policy "$AMBIGUOUS_POLICY" \
  --dna-tag "$DNA_TAG" \
  --trust-remote-code \
  ${CONFIG_OVERRIDE:+--config-override $CONFIG_OVERRIDE} \
  --max-records "$MAX_RECORDS"

ln -sfn "$RUN_OUT" "$OUTROOT/latest"
echo "[OK] model=$(basename "$CARBON_DIR")  max_length=$MAX_LENGTH  dtype=$DTYPE"
echo "[OK] embed_latest=$OUTROOT/latest"
echo "[OK] embeddings=$RUN_OUT/embeddings.npy"
echo "[OK] embedding_index=$RUN_OUT/embedding_index.tsv"
