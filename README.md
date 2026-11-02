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

> **Licence:** MIT (see [`LICENSE`](LICENSE)). To cite this code, see [`CITATION.cff`](CITATION.cff).

## Reported models

- ***M. tuberculosis*:** embeddings of each drug's target loci (the loci the
  TB-Profiler database annotates for it plus those reported for it in the CRyPTIC
  genome-wide association study, excluding loci annotated only for other drugs;
  target set `gwasfiltered`), four principal components per target, plus
  TB-Profiler's own resistance call. Evo 2 is the primary encoder.
- ***P. aeruginosa*:** embeddings of curated chromosomal loci and of acquired
  genes plus ARDaP's own resistance call (S/I/R as one ordered column; arm `CR`)
  for the eight antimicrobials ARDaP calls; for levofloxacin and aztreonam, which
  it does not, acquired-gene presence in place of the call (arm `D`, the
  originally declared model). Evo 2 is the reported encoder here as well.

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
| [`02_read_qc_assembly_annotation`](workflow/02_read_qc_assembly_annotation) | *P. aeruginosa* read trimming, assembly and polishing, read and assembly QC, Bakta annotation, AMRFinderPlus |
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

## Reproducing the manuscript figures and tables

The scripts that make the figures and tables of the paper are in
[`13_reporting/shared`](workflow/13_reporting/shared) (installed to `$MTB/scripts`).
They read intermediate results of the cluster runs (results workbooks, per-isolate
held-out predictions, fold files, SHAP outputs), which are not part of this
repository. Put them in one folder and point the scripts to it:

```bash
export AMR_INPUT_DIR=/path/to/staged_inputs     # default: ./inputs
export AMR_OUTPUT_DIR=/path/to/figure_output    # default: ./figures_final or ./figures_submission
build_final_figures.sh          # figures 1-5, Extended Data 1-3, Supplementary 1-6, with one-line titles
build_submission_figures.sh     # the same without titles, at the journal's size and font rules
```

`AMR_INPUT_DIR` only sets the defaults: every script also takes each input as a
command-line option (`--help` lists them), and the table builders write to the
current folder. The names below are those of our local folders; rename them or pass
your own paths.

| input (under `AMR_INPUT_DIR`) | what it is | made by |
|---|---|---|
| `AMRPredict_panel_results_wider.xlsx` | *M. tuberculosis* results workbook, reported model and comparators | `build_results_workbook.py` |
| `AMRPredict_pseaer_reported_results.xlsx`, `AMRPredict_pseaer_strict_results.xlsx` | *P. aeruginosa* workbooks: reported model; original specification | `build_pa_workbook.sh` |
| `AMRPredict_isolates_used.xlsx` | isolates analysed, with sets, folds and predictions (Supplementary Table 1) | `build_isolate_table.py` |
| `mtb_wider_figure_inputs`, `mtb_holdout_20261012`, `mtb_slim_20261014`, `mtb_slim_20261016`, `mtb_slim_final_20261012`, `mtb_verdict_20261012` | *M. tuberculosis* per-isolate predictions, fold files, development runs, feature-set and lineage results | `slim_mtb_runs.sh` and the stage 09–11 jobs |
| `pa_strict`, `pa_strict_20261020`, `pa_reported3_20261019`, `pa_prescall_20261019` | *P. aeruginosa* predictions and folds: original specification, reported model, presence + ARDaP call baseline | `pa/slim_runs.sh`, `pa/stage_reported_pa.py` |
| `mtb_shap_gwfrule1k`, `pa_reported_20261021` | TreeSHAP outputs | stage 12 jobs |
| `tbprofiler.variants.csv`, `pseaer_metadata_flat.tsv` | TB-Profiler variants; *P. aeruginosa* isolate metadata (study of origin) | stage 04 and 08 |

Panel b of Figure 1 is a hand-drawn illustration that is not part of this
repository: `make_figure1.py` leaves the panel blank, and says so, when the image
is not found. Some figure scripts also contain header text from earlier layouts
(for example "Figure 3 (draft)"); it is printed only when neither `FIG_FINAL=1` nor
`FIG_SUBMISSION=1` is set, which is what the two build scripts set.

## Provenance

The files here are copies of the working scripts that produced the reported
results. Only site-specific paths were changed: hard-coded cluster paths became
required configuration variables, and a few lines that found helper scripts through
a sibling folder of the original repository now use `$MTB_SCRIPTS` when it is set
(two Python files gain `import os` for this). Of 170 files, 99 are byte-identical
to the working copies and 71 differ in 170 lines, every one listed with its
category in [`docs/CHANGES_FROM_WORKING_COPY.tsv`](docs/CHANGES_FROM_WORKING_COPY.tsv).
`MANIFEST.tsv` gives each file's original path and checksums of both versions.

## Not included

- **Data.** *M. tuberculosis* reads are under the ENA accessions in the CRyPTIC
  consortium's supplementary data (PLoS Biol 20:e3001755, 2022), as are the MICs
  and the curated association hits (S1 and S2 Data). *P. aeruginosa* accessions are
  listed in the paper's Data availability statement.
- **The array submission for per-isolate TB-Profiler extraction.** The per-isolate
  worker is included (`04_catalogue_calls/mtb/run_extract_tbprofiler_one_sample.sbatch`);
  the loop that ran it over all isolates was not recorded.
- **The TB-Profiler database** (tbdb commit 72ef6fa) and the reference genomes,
  which are obtained from their sources.

## Citation

See [`CITATION.cff`](CITATION.cff). The archived release and its DOI will be added at publication.
