# 03 — Variant calling

**_M. tuberculosis_.** `submit_mtb_snippy_all.sbatch` runs Snippy 4.6.0 on each
isolate's raw paired reads against H37Rv (`run_mtb_snippy_chunk.sh` →
`run_mtb_snippy_sample.sh`), with Snippy's default filters.

**_P. aeruginosa_.** `submit_pseaer_01_run_snippy_hybrid_e2e.sh` → wrappers
`01_run_snippy_hybrid_e2e.sh` / `01a_snippy_hybrid_worker.sh`, against PAO1
(NC_002516.2), with thresholds from `pseaer.env`.

**Allele fraction.** Sequences for embedding are rebuilt from calls at a minimum
allele fraction of **0.1** rather than Snippy's consensus, so minority resistant
populations are represented:

- `filter_vcf_by_allele_fraction.py` — rebuild a Snippy VCF at a chosen fraction.
  It writes plain text whatever the file name says.
- `sweep_consensus_threshold.py` — how well each threshold represents the measured
  MIC, scored on development isolates only. 0.1 was chosen this way.
