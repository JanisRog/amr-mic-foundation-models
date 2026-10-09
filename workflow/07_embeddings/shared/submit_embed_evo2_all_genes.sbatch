#!/usr/bin/env bash
#SBATCH --job-name=mtb_evo2_allgene
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PIPE="${PIPE:?set PIPE in config/paths.env}"
PYTHON_BIN="${PYTHON_BIN:-${PY:?set PY in config/paths.env}}"

CHAIN="${CHAIN:-$MTB/results/gwas/sample_windows/mtb_all_gene_consensus_chain/latest}"
IN_FASTA="${IN_FASTA:-$CHAIN/unique_window_sequences.fasta}"
OUTROOT="${OUTROOT:-$MTB/results/gwas/evo2_embeddings/evo2_7b/mtb_all_gene_consensus_chain_finalnorm_09b}"
MODEL_ID="${MODEL_ID:-arcinstitute/evo2_7b}"
SPLIT_THRESHOLD="${SPLIT_THRESHOLD:-8192}"
SPLIT_PARTS="${SPLIT_PARTS:-2}"
EXPECTED_DIM="${EXPECTED_DIM:-4096}"
HOOK_LAYER="${HOOK_LAYER:-norm}"
MAX_RECORDS="${MAX_RECORDS:-0}"
# Set DELTA_MAP to run 09f (per-position deltas from a per-target anchor)
# instead of 09b. It needs the sample_target_seqid_map to group alleles by
# target and pick each target's modal allele as the anchor. All the CUDA /
# transformer_engine workarounds below are why this branches here rather than
# living in its own submit script.
DELTA_MAP="${DELTA_MAP:-}"
DELTA_TOPK="${DELTA_TOPK:-8}"
DELTA_REDUCTIONS="${DELTA_REDUCTIONS:-sum,maxabs,topk}"

[[ -s "$IN_FASTA" ]] || { echo "[ERROR] missing input fasta: $IN_FASTA" >&2; exit 2; }

module load cuda >/dev/null 2>&1 || true

# Compute nodes reach huggingface.co only through a proxy that is frequently
# refusing connections. The evo2 weights are already in the shared HF cache, so
# resolve from it rather than re-checking the revision online. Set
# HF_HUB_OFFLINE=0 explicitly if the model genuinely needs downloading.
export HF_HUB_OFFLINE="${HF_HUB_OFFLINE:-1}"
export TRANSFORMERS_OFFLINE="${TRANSFORMERS_OFFLINE:-1}"
export HF_HOME="${HF_HOME:-$HOME/.cache/huggingface}"
echo "[INFO] HF_HOME              = $HF_HOME"
echo "[INFO] HF_HUB_OFFLINE       = $HF_HUB_OFFLINE"

# transformer_engine's _load_nvrtc() runs `ldconfig -p | grep libnvrtc` and dies
# when grep finds nothing, which happens on nodes where libnvrtc lives only in
# the conda env or a module rather than the system ld.so cache. Shim ldconfig -p
# so it also reports the libnvrtc copies we can actually see.
CONDA_LIB="$(dirname "$(dirname "$PYTHON_BIN")")/lib"
NVRTC_FOUND=()
for d in "$CONDA_LIB" "$CONDA_LIB"/python3.12/site-packages/nvidia/cuda_nvrtc/lib \
         "${CUDA_HOME:-/usr/local/cuda}/lib64" /usr/local/cuda*/lib64; do
  for f in "$d"/libnvrtc.so*; do
    [[ -e "$f" ]] && NVRTC_FOUND+=("$f")
  done
done

if [[ ${#NVRTC_FOUND[@]} -gt 0 ]]; then
  SHIM="$(mktemp -d)"
  {
    echo '#!/usr/bin/env bash'
    echo 'if [[ " $* " == *" -p "* || "$1" == "-p" ]]; then'
    echo '  /sbin/ldconfig -p 2>/dev/null || true'
    for f in "${NVRTC_FOUND[@]}"; do
      echo "  echo '        $(basename "$f") (libc6,x86-64) => $f'"
    done
    echo '  exit 0'
    echo 'fi'
    echo 'exec /sbin/ldconfig "$@"'
  } > "$SHIM/ldconfig"
  chmod +x "$SHIM/ldconfig"
  export PATH="$SHIM:$PATH"
  export LD_LIBRARY_PATH="$(dirname "${NVRTC_FOUND[0]}"):${LD_LIBRARY_PATH:-}"
  echo "[INFO] ldconfig shim active; libnvrtc: ${NVRTC_FOUND[*]}"
else
  echo "[WARN] no libnvrtc found; transformer_engine import may fail" >&2
fi

# transformer_engine also calls _nvidia_cudart_include_dir(), which does
# Path(nvidia.__file__) -- but `nvidia` is a namespace package so __file__ is
# None and it raises TypeError. Setting NVTE_CUDA_INCLUDE_DIR short-circuits the
# `not os.getenv(...) and _nvidia_cudart_include_dir()` guard before it runs.
if [[ -z "${NVTE_CUDA_INCLUDE_DIR:-}" ]]; then
  for d in "${CUDA_HOME:-}/include" \
           "$CONDA_LIB"/python3.12/site-packages/nvidia/cuda_runtime/include \
           /usr/local/cuda*/include; do
    if [[ -f "$d/cuda_runtime.h" ]]; then
      export NVTE_CUDA_INCLUDE_DIR="$d"
      break
    fi
  done
fi
if [[ -z "${NVTE_CUDA_INCLUDE_DIR:-}" ]]; then
  # no real include dir found; still set it so the buggy probe is skipped
  export NVTE_CUDA_INCLUDE_DIR="${CUDA_HOME:-/usr/local/cuda}/include"
  echo "[WARN] no cuda_runtime.h found; NVTE_CUDA_INCLUDE_DIR set to a placeholder" >&2
fi
echo "[INFO] NVTE_CUDA_INCLUDE_DIR = $NVTE_CUDA_INCLUDE_DIR"

# `module load cuda` above puts CUDA 13.1 on PATH, but this Triton (shipped with
# torch 2.7.x) rejects any CUDA newer than its version table -- "Triton only
# support CUDA 10.0 or higher, but got CUDA version: 13.1". Pin ptxas to the
# conda copy, which is version-matched to the installed Triton.
CONDA_BIN="$(dirname "$PYTHON_BIN")"
if [[ -z "${TRITON_PTXAS_PATH:-}" && -x "$CONDA_BIN/ptxas" ]]; then
  export TRITON_PTXAS_PATH="$CONDA_BIN/ptxas"
fi
echo "[INFO] TRITON_PTXAS_PATH    = ${TRITON_PTXAS_PATH:-<unset>}"
if [[ -n "${TRITON_PTXAS_PATH:-}" ]]; then
  echo "[INFO] ptxas version        = $("$TRITON_PTXAS_PATH" --version 2>/dev/null | tail -2 | tr '\n' ' ')"
fi

mkdir -p "$OUTROOT" "$MTB/00_logs/evo2"
STAMP="$(date -u +%Y%m%d_%H%M%S)"
RUN_OUT="$OUTROOT/run_${STAMP}"
mkdir -p "$RUN_OUT"

# evo2 tokenizes per nucleotide, so long windows must be halved first. The map
# is what 13c/add_embedding_rows uses to stitch the halves back together.
SPLIT_FASTA="$RUN_OUT/unique_window_sequences.split_long_gt${SPLIT_THRESHOLD}_halves.fasta"
SPLIT_MAP="$RUN_OUT/unique_window_sequences.split_long_gt${SPLIT_THRESHOLD}_halves.map.tsv"

echo "[INFO] input_fasta   = $IN_FASTA ($(grep -c '^>' "$IN_FASTA") sequences)"
echo "[INFO] split_above   = ${SPLIT_THRESHOLD} bp into ${SPLIT_PARTS} parts"

"$PYTHON_BIN" "$MTB/scripts/split_long_fasta_records.py" \
  --input-fasta "$IN_FASTA" \
  --output-fasta "$SPLIT_FASTA" \
  --map-tsv "$SPLIT_MAP" \
  --threshold "$SPLIT_THRESHOLD" \
  --parts "$SPLIT_PARTS"

echo "[INFO] split_fasta   = $SPLIT_FASTA ($(grep -c '^>' "$SPLIT_FASTA") sequences)"
echo "[INFO] model         = $MODEL_ID"
echo "[INFO] run_out       = $RUN_OUT"

EXTRA=()
[[ "$MAX_RECORDS" != "0" ]] && EXTRA+=(--max-records "$MAX_RECORDS")

if [[ -n "$DELTA_MAP" ]]; then
  [[ -s "$DELTA_MAP" ]] || { echo "[ERROR] DELTA_MAP not readable: $DELTA_MAP" >&2; exit 2; }
  echo "[INFO] mode          = positional delta (09f)"
  echo "[INFO] delta_map     = $DELTA_MAP"
  echo "[INFO] reductions    = $DELTA_REDUCTIONS  topk=$DELTA_TOPK"

  "$PYTHON_BIN" "$PIPE/scripts/09f_evo2_positional_delta_embeddings.py" \
    --input-fasta "$SPLIT_FASTA" \
    --map-tsv "$DELTA_MAP" \
    --out-root "$OUTROOT" \
    --script-09b "$PIPE/scripts/09b_embed_with_evo2_hidden.py" \
    --model-id "$MODEL_ID" \
    --device auto \
    --hook-layer "$HOOK_LAYER" \
    --fp8-mode disable \
    --topk "$DELTA_TOPK" \
    --reductions "$DELTA_REDUCTIONS" \
    "${EXTRA[@]}"

  # 09f writes one variant directory per reduction, each with its own latest/,
  # so the single-run symlink below does not apply.
  echo "[OK] model=$MODEL_ID  hook=$HOOK_LAYER  mode=delta"
  echo "[OK] variants under $OUTROOT/delta_*"
  exit 0
fi

"$PYTHON_BIN" "$PIPE/scripts/09b_embed_with_evo2_hidden.py" \
  --outdir "$RUN_OUT" \
  --device auto \
  --fp8-mode disable \
  --embedding-source residual \
  --expected-embedding-dim "$EXPECTED_DIM" \
  --input-fasta "$SPLIT_FASTA" \
  --hook-layer "$HOOK_LAYER" \
  --model-id "$MODEL_ID" \
  "${EXTRA[@]}"

ln -sfn "$RUN_OUT" "$OUTROOT/latest"
echo "[OK] model=$MODEL_ID  hook=$HOOK_LAYER  split_threshold=$SPLIT_THRESHOLD"
echo "[OK] embed_latest=$OUTROOT/latest"
echo "[OK] embeddings=$RUN_OUT/embeddings.npy"
echo "[OK] embedding_index=$RUN_OUT/embedding_index.tsv"
echo "[OK] split_map=$SPLIT_MAP"
