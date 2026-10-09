# 02 — Read QC, assembly QC and annotation

*P. aeruginosa* only; *M. tuberculosis* is called from raw reads (stage 03) and
needs no assembly.

Each `submit_pseaer_*.sh` sources `pa/pseaer.env` and launches the matching wrapper,
which fans out per-isolate workers. Run in this order:

| step | submit | tool |
|---|---|---|
| read QC | `submit_pseaer_00a_fastqc.sh`, `submit_pseaer_00b_multiqc_reads.sh` | FastQC 0.11.9, MultiQC 1.33 |
| assembly QC | `submit_pseaer_00c_run_quast_all.sh` | QUAST 5.3.0 |
| annotation | `submit_pseaer_02r_bakta_reference.sh` (PAO1), `submit_pseaer_02_run_bakta_all_e2e.sh` | Bakta 1.12.0, db 6.0 |
| acquired determinants | `wrappers: 02b_run_amrfinderplus_all_e2e.sh` | AMRFinderPlus 4.0.23, db 2025-06-03.1 |

`shared/_shared_gwas.sh` holds the helpers the wrappers share (`require_file`,
`require_dir`, logging). Every tool path is checked before use, so an unset image
fails immediately with "missing file".

`pseaer.env` is the shared configuration of the read-processing pipeline. It also
defines settings for association-analysis stages that are not part of this
release; their tool paths are left empty.

**Not included:** read trimming (Trimmomatic), assembly (Unicycler) and polishing
(Pilon). Stage 02 reads the polished assemblies from
`$PSEAER/assemblies/unicycler_pilon`.
