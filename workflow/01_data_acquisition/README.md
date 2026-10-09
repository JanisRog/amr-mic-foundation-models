# 01 — Data acquisition

*M. tuberculosis* only. *P. aeruginosa* starts at stage 02 from reads and polished
assemblies obtained from the source collection.

**Reads.** `submit_mtb_ena_reads.sh` downloads paired reads from ENA in chunks
(`download_mtb_ena_reads_chunk_worker.sh` → `download_mtb_ena_reads_worker.py`).
`build_mtb_reads_missing_list.py` and `submit_mtb_ena_reads_missing_chunked.sh`
retry what failed. `build_paired_sample_list.py` writes the list of isolates with
both mates, which variant calling and TB-Profiler read.

**Phenotypes.** From the CRyPTIC supplementary data (PLoS Biol 20:e3001755, 2022):

- `build_all_drug_mic_table.py` — one log2 MIC table for all thirteen drugs.
- `extract_cryptic_mics.py` — per-drug log2 MIC tables from S2 Data.

**Published association hits.** `cryptic_gwas_to_entries.py` turns S1 Data, the
curated interpretation of the CRyPTIC genome-wide association study, into
per-drug locus entries for stage 05. In S1 Data `:` names an intergenic region
between two genes; `/` separates alternatives and is **not** intergenic.
