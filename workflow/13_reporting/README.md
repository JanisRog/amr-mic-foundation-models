# 13 — Reporting

Run locally, from downloaded cluster outputs, after `install.sh`. The plotting
modules import one another (`make_panel_figures.py` holds the shared style), so
run them from the installed `$MTB/scripts`; the *P. aeruginosa* drivers find them
through `$MTB_SCRIPTS`.

- `build_results_workbook.py` — development CV, held-out scores and the comparator,
  one workbook per organism. `holdout_table.py` builds the model-against-comparator
  table from it.
- `mtb/build_mtb_workbook_arm.sh`, `pa/build_pa_workbook.sh` — the workbook per
  organism, one call per drug with models and comparator together.
- `mtb/make_all_figures.sh`, `pa/make_all_pa_figures.sh` — every figure in one
  command.
- `pa/stage_reported_pa.py` — assembles the reported *P. aeruginosa* model (ARDaP's call
  for the eight antimicrobials ARDaP calls, AMRFinderPlus presence for levofloxacin and
  aztreonam) under the names the workbook and figure scripts read, with a manifest of
  every source file; run it before `pa/build_pa_workbook.sh` and `pa/make_all_pa_figures.sh`.
- `make_panel_figures.py`, `plot_mic_scatter.py`, `plot_class_distribution.py`,
  `plot_model_comparison.py`, `plot_f1_boxes.py`, `plot_metric_bars.py`, `make_panel_grids.py` — the
  individual figures; each recomputes its metrics from per-isolate predictions and
  fails if they disagree with the workbook by more than 5 × 10⁻⁴.
- `build_figure_page.py` — all figures and their numbers on one review page.

The manuscript figures and tables (inputs: see the main README, "Reproducing the
manuscript figures and tables"):

- `shared/make_figure1.py` … `make_figure5.py` (`make_figure3_combined.py` is Figure 3) —
  main Figures 1–5. `make_extended_figure1.py`, `make_extended_figure2.py`,
  `make_extended_figure3.py` — Extended Data Figures 1–3. `make_supp_scatters.py` —
  Supplementary Figures 1–2; Supplementary Figures 3–6 come from `mtb/make_all_figures.sh`
  and `pa/make_all_pa_figures.sh`. `make_figure3.py`, `make_figure3_alt.py` and
  `make_provisional_figures.py` hold earlier layouts and shared helpers.
- `shared/build_final_figures.sh`, `shared/build_submission_figures.sh` — all of the above
  in one command, with one-line titles or without titles at the journal's size.
- `shared/build_isolate_table.py`, `shared/build_target_loci_table.py`,
  `shared/build_pa_target_loci_table.py`, `shared/build_comparison_metrics_table.py`,
  `shared/build_design_alternatives_table.py` — Supplementary Tables 1–4.
  `shared/fda_criteria.py` scores the held-out error counts against the FDA acceptance
  criteria for susceptibility test systems (sheet `FDA_criteria` of Supplementary Table 3).

