# 12 — Interpretation

`shared/shap_target_values.py` computes TreeSHAP values per isolate for the reported
model on held-out isolates, against the continuous log2 MIC output. SHAP is
additive, so a target's contribution is the exact sum over its four principal
components. Run through `mtb/submit_mtb_shap_arm.sh` → `mtb/run_shap.sbatch` for the reported
*M. tuberculosis* model (the submitter passes its held-out cache, selected
configuration and target count) and `pa/submit_pa_shap.sh` →
`pa/run_pa_shap.sbatch`. `shared/plot_shap_beeswarm.py` draws target-level beeswarms.
