#!/usr/bin/env bash
# Assemble the per-model ML roots the MTB precompute expects.
#
# `cv_precompute_fold_features.py` reads three files from
# <root>/<embedding-subdir>/ : embeddings.npy, embedding_index.tsv and
# sample_target_seqid_map.tsv. The embeddings land in per-model run directories
# under results/gwas/, and the map is shared by all three because every model
# embedded the same FASTA. This links them into one place per model.
#
# Symlinks, not copies: the arrays are 40 MB - 230 MB each and there is no
# reason to have two of them. The run directory stays the provenance record.
#
# The precompute's --embedding-subdir default is `ml_inputs/embeddings/tb_profiler`,
# which for MTB is the all-gene set and here does not exist. Every P. aeruginosa
# job must pass EMB_SUBDIR explicitly -- the MTB panel lost an afternoon to a
# cache whose embedding source could not be recovered afterwards because the
# manifest does not record it.

set -euo pipefail

P="${P:-${PSEAER:?set PSEAER in config/paths.env}}"
SET="${SET:-curated_af010}"
MAP="${MAP:-$P/ml_inputs/embeddings/$SET/sample_target_seqid_map.tsv}"

[ -s "$MAP" ] || { echo "[ERROR] no map at $MAP"; exit 1; }
n_map=$(( $(wc -l < "$MAP") - 1 ))
echo "[INFO] map: $n_map sample-target rows"

link_one () {  # model  run_glob  ml_root
  local model="$1" glob="$2" root="$3"
  local run
  run=$(ls -dt $glob 2>/dev/null | head -1 || true)
  if [ -z "$run" ] || [ ! -s "$run/embeddings.npy" ]; then
    echo "[SKIP] $model: no completed run under $glob"
    return 0
  fi
  local dst="$root/ml_inputs/embeddings/$SET"
  mkdir -p "$dst"
  ln -sfn "$run/embeddings.npy"       "$dst/embeddings.npy"
  ln -sfn "$run/embedding_index.tsv"  "$dst/embedding_index.tsv"
  ln -sfn "$MAP"                      "$dst/sample_target_seqid_map.tsv"
  local rows dim
  rows=$(( $(wc -l < "$run/embedding_index.tsv") - 1 ))
  dim=$(python3 -c "import numpy as np;print(np.load('$run/embeddings.npy',mmap_mode='r').shape[1])" 2>/dev/null || echo "?")
  echo "[OK]   $model: $rows sequences, dim $dim"
  echo "       $dst"
}

link_one evo2 \
  "$P/results/gwas/evo2_embeddings/evo2_7b/${SET}_finalnorm_09b/run_*" \
  "$P/results/ml/pa_evo2"
link_one carbon \
  "$P/results/gwas/carbon_embeddings/carbon_8b/$SET/run_*" \
  "$P/results/ml/pa_carbon"
link_one dnabert \
  "$P/results/gwas/dnabert_embeddings/dnabert2_117m/$SET/run_*" \
  "$P/results/ml/pa_dnabert"

echo
echo "Every sequence is under 8,192 bp, so no split map is needed and all three"
echo "models see identical input -- unlike MTB, where the longest 1% reached evo2"
echo "and DNABERT-2 as the mean of two halves."
echo
echo "Pass these to the precompute explicitly:"
echo "  --root $P/results/ml/pa_<model>  --embedding-subdir ml_inputs/embeddings/$SET"
