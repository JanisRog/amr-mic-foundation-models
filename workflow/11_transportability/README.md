# 11 — Transportability

Post hoc tests on the development partition; the held-out isolates are not used.

## *M. tuberculosis* — leave one lineage out

`run_panel_lineage.sbatch` (one array task per drug):

- groups isolates by `main_lineage`, restricted to lineages 1–4 (isolates marked
  `mixed` or `unknown` are excluded), one lineage per fold, fold count derived
  from the lineages present;
- refits each feature set with its **development-selected configuration and
  top-k**, and fails unless the refit's feature count matches development;
- pairs every lineage split with random folds of identical sizes, since lineages
  differ about tenfold in size and holding one out also shrinks training.

`check_lineage_validity.py` checks feature parity, fold shape, the size control and
the moxifloxacin identity control, and prints no comparison between feature sets.
`analyse_lineage_arms.py` computes the comparison.

## *P. aeruginosa* — leave one study out

`build_loso_folds.py` holds out each collection with at least 40 isolates and eight
per class, leaving at least 100 training isolates with both classes.
`run_pa_loso.sbatch` refits the reported configuration and, on the same isolates,
the same configuration under ordinary cross-validation. `analyse_loso.py` reports
the difference; `plot_loso.py` relates it to training data removed and to the
collection's resistance-rate gap.
