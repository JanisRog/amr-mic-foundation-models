# 06 — Sequence extraction

One sequence per isolate per target, then de-duplicated: each unique sequence is
embedded once.

## *M. tuberculosis*

`submit_mtb_13l_consensus_chain.sbatch` runs
`shared/13l_build_snp_gene_consensus_chain_sequences.py`. It rebuilds each isolate's
consensus at allele fraction 0.1 and projects the reference target intervals onto
it through the reference-to-consensus chain, preserving isolate-specific variants
within a fixed target definition. `REGENERATE=0` reuses existing consensus files.
`build_sample_target_seqid_map.py` joins the per-isolate sequences to their
sequence ids.

`split_long_fasta_records.py` halves sequences longer than 8,192 bp. Evo 2 and
DNABERT-2 embed the halves and average them; Carbon reads every sequence whole.
This affects 1.04% of unique sequences, all in the longest genes.

## *P. aeruginosa*

`run_pseaer_extract.sbatch` rebuilds each VCF at allele fraction 0.1, then
`extract_pseaer_target_sequences.py` applies the variants **inside each target
window only** — no coordinate chain — so an indel shifts the rest of its own
window and cannot move another window's boundaries. Contig names are matched to
the reference by verified length; a mismatch is an error, not a silently
unmodified reference. `gather_pseaer_sequences.py` collapses the result into the
embedding layout. `extract_amr_gene_sequences.py` takes acquired-gene sequences
from each isolate's own assembly, at the AMRFinderPlus hit coordinates.
No *P. aeruginosa* sequence exceeds 8,192 bp.
