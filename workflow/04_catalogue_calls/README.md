# 04 — Catalogue calls

The comparators, and the catalogue-derived features.

## *M. tuberculosis* — TB-Profiler 6.3.0, database tbdb commit 72ef6fa

1. `submit_mtb_tbprofiler.sh` → `run_tbprofiler_chunk.sh` → `run_tbprofiler_sample.sh`:
   TB-Profiler on each isolate's reads. `gather_mtb_tbprofiler.sh` collects the results.
2. `run_extract_tbprofiler_one_sample.sbatch` → `extract_tbprofiler_target_sequences.py`,
   once per isolate (`SAMPLE=<id>`, `OUTDIR=$MTB/results/gwas/sample_windows/tbprofiler_all_targets_gene_up300/samples/<id>`):
   per-isolate tables including TB-Profiler's own per-drug resistance call. The loop
   that submitted it over all isolates was not recorded.
   `submit_gather_tbprofiler_target_sequences.sbatch` → `gather_tbprofiler_target_sequences.py`
   merges them into `tbprofiler_drug_resistance_binary.tsv`, **the source of the
   rule bit** (stage 08).
3. `submit_build_tbprofiler_amrfinder_style_matrix.sbatch` →
   `build_tbprofiler_amrfinder_style_matrix.py`: one binary column per catalogued
   resistance mutation (the mutation-column feature sets).
4. `submit_tbprofiler_calls_by_threshold.sbatch` → `tbprofiler_calls_by_threshold.py`:
   the catalogue's call at allele fractions 0.10 and 0.90, reported as additional
   baselines.
5. `extract_tbprofiler_lineage.py`: per-isolate lineage (`sample_lineage.tsv`), used by
   the leave-one-lineage-out test.

`compare_call_sources_vs_mic.py` checks TB-Profiler's calls against Snippy's
consensus at catalogue positions, scored against the measured MIC.

## *P. aeruginosa* — ARDaP

`submit_ardap_pseaer_all.sh` → `run_ardap_pseaer_worker.sh` per isolate;
`gather_ardap_calls.py` collects the per-drug S/R/I calls. ARDaP covers eight of
the ten drugs; aztreonam and levofloxacin have no comparator.
