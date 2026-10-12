# MIC prediction from DNA foundation-model embeddings

Code for predicting minimum inhibitory concentrations (MICs) of antimicrobials
from frozen embeddings of resistance-associated loci, produced by three DNA
foundation models (Evo 2 7B, Carbon 8B, DNABERT-2 117M), for ten antimicrobials
each in *Mycobacterium tuberculosis* and *Pseudomonas aeruginosa*.

This repository contains the code that produced the published results, and only
that code. Exploratory analyses and approaches that were evaluated and not adopted
are described in the Methods, not shipped here; the two embedding variants kept
under `workflow/07_embeddings/alternatives_not_adopted/` are the exception, as they
are referenced from the embedding scripts.

> **Status:** pre-publication. A licence will be added before release.

## Reported models

- ***M. tuberculosis*:** embeddings of each drug's target loci (the loci the
  TB-Profiler database annotates for it plus those reported for it in the CRyPTIC
  genome-wide association study, excluding loci annotated only for other drugs;
  target set `gwasfiltered`), four principal components per target, plus
  TB-Profiler's own resistance call. Evo 2 is the primary encoder.
- ***P. aeruginosa*:** embeddings of curated chromosomal loci and of acquired
  genes, plus acquired-gene presence; no catalogue call exists for this species.

Scripts for the alternatives that were compared against these (the earlier target
set, mutation columns in place of the call, embeddings alone, catalogue features
alone) are included, because the paper reports those comparisons.

## Layout

The code is organised by pipeline stage, and within a stage by organism
(`mtb/`, `pa/`) or `shared/` when both use it. Each stage folder has a README
giving run order, inputs, outputs and parameters.

| stage | what it does |
|---|---|
| [`01_data_acquisition`](workflow/01_data_acquisition) | *M. tuberculosis* reads from ENA; MIC tables and published association hits from the CRyPTIC supplementary data |
| [`02_read_qc_assembly_annotation`](workflow/02_read_qc_assembly_annotation) | *P. aeruginosa* read and assembly QC, Bakta annotation, AMRFinderPlus |
| [`03_variant_calling`](workflow/03_variant_calling) | Snippy against H37Rv and PAO1; rebuilding calls at a chosen allele fraction |
| [`04_catalogue_calls`](workflow/04_catalogue_calls) | TB-Profiler calls, mutation matrix and lineage; ARDaP calls |
| [`05_target_sets`](workflow/05_target_sets) | Per-drug candidate loci, each as gene body plus 300 bp upstream |
| [`06_sequence_extraction`](workflow/06_sequence_extraction) | Per-isolate target sequences from the consensus, de-duplicated |
| [`07_embeddings`](workflow/07_embeddings) | Evo 2, Carbon and DNABERT-2 embeddings, one per unique sequence |
| [`08_labels_features_splits`](workflow/08_labels_features_splits) | MIC labels, censoring, catalogue features, held-out split and folds |
| [`09_model_training_selection`](workflow/09_model_training_selection) | Fold-wise PCA, random-forest sweep, configuration selection, held-out scoring |
| [`10_evaluation`](workflow/10_evaluation) | Comparator baselines and paired comparisons between feature sets |
| [`11_transportability`](workflow/11_transportability) | Leave-one-lineage-out (*M. tuberculosis*), leave-one-study-out (*P. aeruginosa*) |
| [`12_interpretation`](workflow/12_interpretation) | Per-isolate TreeSHAP, summed per target |
| [`13_reporting`](workflow/13_reporting) | Results workbooks, tables and figures |

## Requirements

- A Slurm cluster with Apptainer/Singularity, and a GPU node for the embedding
  stage (one H100 was used for Evo 2 7B and Carbon 8B).
- One Python environment (Python 3.12 was used). The embedding stage needs PyTorch, `transformers`,
  `tokenizers`, `huggingface_hub`, Biopython and the Evo 2 runtime (`evo2`,
  `vortex`, Transformer Engine, FlashAttention). Everything else needs `numpy`,
  `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `shap` and, for the workbooks,
  `openpyxl`.
- Containers for the tools in [`config/paths.env.example`](config/paths.env.example),
  which also records the versions used.

## Configuration and installation

```bash
cp config/paths.env.example config/paths.env   # then fill in your paths
source config/paths.env
./install.sh --dry-run                          # show what goes where
./install.sh
```

The scripts call one another by fixed locations — `$MTB/scripts/…`,
`$PIPE/wrappers/…` — as they did when the results were produced. `install.sh`
copies each file to the location(s) recorded in [`MANIFEST.tsv`](MANIFEST.tsv)
rather than rewriting the calls. **Run the workflow from the installed tree,
not from these stage folders**; several scripts import siblings from their own
directory. An installed file that differs from the repository is reported and
left alone unless `--force` is given.

Every script reads site paths from the environment and stops with
`set <VAR> in config/paths.env` if one it needs is missing.

## Provenance

The files here are copies of the working scripts that produced the reported
results. Only site-specific paths were changed: hard-coded cluster paths became
required configuration variables, and a few lines that found helper scripts through
a sibling folder of the original repository now use `$MTB_SCRIPTS` when it is set
(two Python files gain `import os` for this). Of 143 files, 76 are byte-identical
to the working copies and 67 differ in 154 lines, every one listed with its
category in [`docs/CHANGES_FROM_WORKING_COPY.tsv`](docs/CHANGES_FROM_WORKING_COPY.tsv).
`MANIFEST.tsv` gives each file's original path and checksums of both versions.

## Not included

- **Data.** *M. tuberculosis* reads are under the ENA accessions in the CRyPTIC
  consortium's supplementary data (PLoS Biol 20:e3001755, 2022), as are the MICs
  and the curated association hits (S1 and S2 Data). *P. aeruginosa* accessions are
  listed in the paper's Data availability statement.
- **_P. aeruginosa_ read trimming, assembly and polishing** (Trimmomatic,
  Unicycler, Pilon). These steps were run, but their scripts are not part of this
  release yet; stage 02 starts from the polished assemblies.
- **The array submission for per-isolate TB-Profiler extraction.** The per-isolate
  worker is included (`04_catalogue_calls/mtb/run_extract_tbprofiler_one_sample.sbatch`);
  the loop that ran it over all isolates was not recorded.
- **The TB-Profiler database** (tbdb commit 72ef6fa) and the reference genomes,
  which are obtained from their sources.

## Citation

To be added on publication.
