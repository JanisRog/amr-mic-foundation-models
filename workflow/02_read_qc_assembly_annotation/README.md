# 02 — Read processing, assembly, QC and annotation

*P. aeruginosa* only; *M. tuberculosis* is called from raw reads (stage 03) and
needs no assembly.

Each `submit_pseaer_*.sh` sources `pa/pseaer.env` and launches the matching wrapper,
which fans out per-isolate workers. Run in this order:

| step | submit | tool |
|---|---|---|
| trimming, assembly, polishing | `03_assemble_pseaer_reads_unicycler.sh` (Slurm array, one isolate per task) | Trimmomatic 0.39, Unicycler 0.5.0, Pilon 1.24 (BWA 0.7.17, samtools 1.15.1) |
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

`03_assemble_pseaer_reads_unicycler.sh` is self-contained and does not source
`pseaer.env`. It trims with Nextera adapter clipping, a four-base sliding window at
mean quality 12 and a minimum length of 100 bases, assembles with Unicycler and
polishes with Pilon, and writes to `$PSEAER/assemblies/unicycler_pilon`, which the
later steps read. It needs `SING_ROOT` (the directory holding the five tool images)
and `ILLUMINACLIP` (the adapter FASTA) from `config/paths.env`. The assemblies are
used only for AMRFinderPlus and acquired-gene sequences; variants are called from
the reads (stage 03).
