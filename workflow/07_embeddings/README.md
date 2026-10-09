# 07 — Embeddings

Frozen models, no fine-tuning. One vector per unique sequence, mean-pooled over
positions.

| model | script | settings | dim |
|---|---|---|---|
| Evo 2 7B | `09b_embed_with_evo2_hidden.py` via `submit_embed_evo2_all_genes.sbatch` | final normalisation layer, FP8 off, 8,192 nt context (split input) | 4,096 |
| Carbon 8B | `09c_embed_with_carbon_hidden.py` via `submit_embed_carbon_tbprofiler_targets.sbatch` | layer 32 (`LAYER=-1`), `DNA_TAG=both`, bf16, whole sequences | 4,096 |
| DNABERT-2 117M | the same script, with `CARBON_DIR` set to the DNABERT-2 directory | split input | 768 |

**Carbon needs `DNA_TAG`.** Without the `<dna>…</dna>` wrapper its tokenizer reads
nucleotides as English text and the embeddings are meaningless.

**ML roots.** `mtb/setup_cluster_ml_root.sh` and `pa/setup_pseaer_ml_roots.sh` build
one analysis root per model, linking the embedding index and the
sample-to-sequence map in the layout stage 09 reads.
`pa/merge_embedding_blocks.py` fuses the chromosomal and acquired-gene embedding
blocks; `mtb/build_multimodel_embedding_block.py` concatenates several models'
embeddings for the same sequences.

`alternatives_not_adopted/` holds the delta-embedding variants (each allele minus
its target's anchor allele; per-position nonlinear reductions). They were
evaluated and not used for any reported model; `submit_embed_evo2_all_genes.sbatch`
runs 09f only when `DELTA_MAP` is set.
