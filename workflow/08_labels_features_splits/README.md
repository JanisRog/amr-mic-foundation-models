# 08 — Labels, features and splits

## Partitioning (both organisms)

`shared/cv_build_folds.py` holds out 20% of isolates stratified on MIC class
(seed 2026), then builds five-fold cross-validation repeated three times on the
remainder. In `cv_folds.tsv`, `fold < 0` is the held-out set. With `--allele-tsv
<table>#<column>` whole groups are dealt to folds intact (used for lineage folds);
`--fold-sizes` builds random folds of given sizes (the size-matched control).
`--restrict-samples` takes **bare sample ids, one per line**.

## *M. tuberculosis*

- `prepare_drug_inputs.py` — per drug, `<drug>_log2_mic.tsv` and the rule bit
  `rule_bit_<code>_tbprofiler.tsv`: **TB-Profiler's own call** from stage 04, not a
  reconstruction. Codes are the rule-bit file codes (moxifloxacin is `mxf`).
- `recover_sr_boundaries.py` — each drug's S/R boundary, solved from the supplied
  interpretations rather than typed.
- `merge_catalogue_features.py` — a drug's mutation columns, optionally with the
  rule bit, for the mutation-column feature sets.

## *P. aeruginosa*

- `prepare_pseaer_drug_inputs.py` — per-drug labels: clamped to the range of exact
  measurements; strict bounds (`>`, `<`) moved one dilution outward, inclusive
  bounds (`≥`, `≤`) kept; the Rebelo et al. collection excluded; two isolates
  excluded from one drug each for off-series MICs.
- `audit_censoring.py` — what the censoring rule changes, per drug.
- `build_amrfinder_presence.py` → `make_drug_presence_tsv.py` — acquired-gene
  presence, families collapsed by explicit list, carrier floor of ten within each
  drug's development cohort.
- `build_pa_cohorts.sbatch` — labels, folds and presence for one drug in one job.
