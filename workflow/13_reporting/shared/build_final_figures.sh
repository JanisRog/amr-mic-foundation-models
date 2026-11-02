#!/usr/bin/env bash
# The final figures (main, Extended Data, Supplementary) in $AMR_OUTPUT_DIR (default figures_final/): one-line bold title
# "Fig. N: title", no description. Figures 1, 2 and Extended Data 1 are built at 9.4 in and scaled to 7.8 in (vector),
# so every main and Extended Data figure is 7.8 in wide and every title is 11 pt. Titles: FINAL_TITLES in make_panel_figures.py.
set -euo pipefail

export AMR_INPUT_DIR="${AMR_INPUT_DIR:-inputs}"   # where the staged inputs live (README)
export FIG_FINAL=1
M="${MTB_SCRIPTS:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"      # the plotting scripts (installed: $MTB/scripts)
PA_DRIVER="${PSEAER:+$PSEAER/scripts}"; PA_DRIVER="${PA_DRIVER:-"$M"/../pseaer}/make_all_pa_figures.sh"
R="${R:-$AMR_INPUT_DIR/pa_reported3_20261019}"
B="${B:-$(mktemp -d)}"
D="${D:-${AMR_OUTPUT_DIR:-figures_final}}"
mkdir -p "$D"
python3 "$M"/make_figure1.py --pa-workbook "$AMR_INPUT_DIR/AMRPredict_pseaer_reported_results.xlsx" \
  --pa-pred-dir "$R/predictions_strict" --pa-fold-dir "$R/slim_strict" --out-dir "$B" --stem Figure_1 >/dev/null
python3 "$M"/make_figure2.py --out-dir "$B" --stem Figure_2 >/dev/null
python3 "$M"/make_figure3_combined.py --out-dir "$B" --stem Figure_3 >/dev/null
python3 "$M"/make_figure4.py --variant shares_quality --out-dir "$B" --stem Figure_4 >/dev/null
python3 "$M"/make_figure5.py --out-dir "$B" --stem Figure_5 >/dev/null
python3 "$M"/make_extended_figure1.py --out-dir "$B" --stem Extended_Data_Figure_1 >/dev/null
python3 "$M"/make_extended_figure2.py --out-dir "$B" --stem Extended_Data_Figure_2 >/dev/null
python3 "$M"/make_extended_figure3.py --out-dir "$B" --stem Extended_Data_Figure_3 >/dev/null
python3 "$M"/make_supp_scatters.py --out-dir "$B" >/dev/null
(FIG_TITLES='{"What the model used":"Supplementary Fig. 5: What the reported model used, M. tuberculosis, every antimicrobial","MIC class distribution":"Supplementary Fig. 3: MIC class distribution and essential agreement, M. tuberculosis"}' \
  IN="$AMR_INPUT_DIR/mtb_wider_figure_inputs" WB="$AMR_INPUT_DIR/AMRPredict_panel_results_wider.xlsx" \
  OUT="$B/mtb_tmp" bash "$M"/make_all_figures.sh >/dev/null)
(FIG_TITLES='{"What the model used":"Supplementary Fig. 6: What the reported model used, P. aeruginosa, every antimicrobial","MIC class distribution":"Supplementary Fig. 4: MIC class distribution and essential agreement, P. aeruginosa"}' \
  IN="$R" TAG=reported1k WB="$AMR_INPUT_DIR/AMRPredict_pseaer_reported_results.xlsx" \
  MTB_WB="$AMR_INPUT_DIR/AMRPredict_panel_results_wider.xlsx" OUT="$B/pa_tmp" CROSS_OUT="$B/cross_tmp" \
  bash "$PA_DRIVER" >/dev/null)
python3 - "$B" "$D" <<'PY'
import sys, pathlib, fitz
from PIL import Image
B, D = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
W = 7.8
for n in ["Figure_1", "Figure_2", "Figure_3", "Figure_4", "Figure_5", "Extended_Data_Figure_1",
          "Extended_Data_Figure_2", "Extended_Data_Figure_3"]:
    src = fitz.open(B / f"{n}.pdf"); r = src[0].rect; k = W * 72 / r.width
    out = fitz.open(); pg = out.new_page(width=r.width * k, height=r.height * k); pg.show_pdf_page(pg.rect, src, 0)
    out.save(D / f"{n}.pdf", garbage=3, deflate=True)
    im = Image.open(B / f"{n}.png").convert("RGB"); w, h = im.size; nw = round(W * 400)
    im.resize((nw, round(h * nw / w)), Image.LANCZOS).save(D / f"{n}.png", dpi=(400, 400))
import shutil
for e in ("pdf", "png"):
    shutil.copy(B / f"supplementary_figure1_scatters_mtb_draft.{e}", D / f"Supplementary_Figure_1.{e}")
    shutil.copy(B / f"supplementary_figure2_scatters_pa_draft.{e}", D / f"Supplementary_Figure_2.{e}")
    shutil.copy(B / f"mtb_tmp/sheet_shap.{e}", D / f"Supplementary_Figure_5.{e}")
    shutil.copy(B / f"pa_tmp/pa_sheet_shap.{e}", D / f"Supplementary_Figure_6.{e}")
    shutil.copy(B / f"mtb_tmp/sheet_classdist.{e}", D / f"Supplementary_Figure_3.{e}")
    shutil.copy(B / f"pa_tmp/pa_sheet_classdist.{e}", D / f"Supplementary_Figure_4.{e}")
print("[OK]", D)
PY
