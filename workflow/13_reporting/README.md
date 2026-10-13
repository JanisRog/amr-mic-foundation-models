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
- `make_panel_figures.py`, `plot_mic_scatter.py`, `plot_class_distribution.py`,
  `plot_model_comparison.py`, `plot_f1_boxes.py`, `plot_metric_bars.py`, `make_panel_grids.py` — the
  individual figures; each recomputes its metrics from per-isolate predictions and
  fails if they disagree with the workbook by more than 5 × 10⁻⁴.
- `build_figure_page.py` — all figures and their numbers on one review page.
