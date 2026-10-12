# 10 — Evaluation

Metrics: very major error rate (vMER, resistant called susceptible, over resistant
isolates), major error rate (MER, over susceptible isolates), cost = 2·vMER + MER
(the declared objective), F1 for the resistant class, and on the MIC ladder:
one-tier accuracy, balanced one-tier accuracy, the same over classes with at least
ten isolates, exact-class accuracy and mean absolute error in dilutions. A
catalogue has no ladder position, so ladder metrics are reported for models only.

## *M. tuberculosis*

- `catalogue_baseline_rows.py` — TB-Profiler's score on the same isolates, on both
  splits, plus the catalogue at allele fractions 0.10 and 0.90.
- `compare_arms.py` — paired comparison of feature sets, fold by fold and across
  drugs (Wilcoxon signed-rank).
- `build_best_config_workbook.py` — every feature set's selected configuration,
  every metric, TB-Profiler beside it. Needs `openpyxl`.
- `build_feature_ablation_table.py` — six feature sets on identical development
  folds (embeddings, mutation columns and the call, in every combination that was
  run), with the paired comparison against embeddings + call.
- `holdout_arm_comparison.py` — held-out comparison of the reported model, the
  earlier target set and mutation columns + call, recomputed from per-isolate
  predictions; it reports only after reproducing the reference values.
- `holdout_subject_sensitivity.py` — held-out scores with and without isolates
  whose subject also has an isolate in the development set.

## *P. aeruginosa*

`run_pa_ardap.sbatch` → `ardap_baseline.py`: ARDaP on the same isolates, under both
conventions for the intermediate category; the one more favourable to ARDaP is
reported.
