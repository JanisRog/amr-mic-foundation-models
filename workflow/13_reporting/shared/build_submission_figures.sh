#!/usr/bin/env bash
# Submission figures (Nature specifications): no title (title and description go in the legend), 180 mm wide, at most
# 170 mm tall, all text 5-7 pt sans-serif (Arial), panel letters 8 pt bold lowercase, vector PDF with embedded fonts.
# Main and Extended Data figures: figures_submission/. The supplementary sheets keep their own size.
set -euo pipefail

export AMR_INPUT_DIR="${AMR_INPUT_DIR:-inputs}"   # where the staged inputs live (README)
export FIG_SUBMISSION=1
M="${MTB_SCRIPTS:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"      # the plotting scripts (installed: $MTB/scripts)
PA_DRIVER="${PSEAER:+$PSEAER/scripts}"; PA_DRIVER="${PA_DRIVER:-"$M"/../pseaer}/make_all_pa_figures.sh"
R="${R:-$AMR_INPUT_DIR/pa_reported3_20261019}"
B="${B:-$(mktemp -d)}"
D="${D:-${AMR_OUTPUT_DIR:-figures_submission}}"
W=7.0866                       # 180 mm
mkdir -p "$D"
# Figure 1 in the original tall proportions (for Nature Communications) and a compact version, 170 mm deep (for Nature)
FIG_TALL=1 python3 "$M"/make_figure1.py --width $W --pa-workbook "$AMR_INPUT_DIR/AMRPredict_pseaer_reported_results.xlsx" \
  --pa-pred-dir "$R/predictions_strict" --pa-fold-dir "$R/slim_strict" --out-dir "$B" --stem Figure_1 >/dev/null
FIG_TALL=0 python3 "$M"/make_figure1.py --width $W --pa-workbook "$AMR_INPUT_DIR/AMRPredict_pseaer_reported_results.xlsx" \
  --pa-pred-dir "$R/predictions_strict" --pa-fold-dir "$R/slim_strict" --out-dir "$B" --stem Figure_1_compact >/dev/null
python3 "$M"/make_figure2.py --width $W --out-dir "$B" --stem Figure_2 >/dev/null
python3 "$M"/make_figure3_combined.py --width $W --out-dir "$B" --stem Figure_3 >/dev/null
python3 "$M"/make_figure4.py --width $W --variant shares_quality --out-dir "$B" --stem Figure_4 >/dev/null
python3 "$M"/make_figure5.py --width $W --out-dir "$B" --stem Figure_5 >/dev/null
python3 "$M"/make_extended_figure1.py --width $W --out-dir "$B" --stem Extended_Data_Figure_1 >/dev/null
python3 "$M"/make_extended_figure2.py --width $W --out-dir "$B" --stem Extended_Data_Figure_2 >/dev/null
python3 "$M"/make_extended_figure3.py --width $W --out-dir "$B" --stem Extended_Data_Figure_3 >/dev/null
for n in Figure_1 Figure_1_compact Figure_2 Figure_3 Figure_4 Figure_5 Extended_Data_Figure_1 Extended_Data_Figure_2 Extended_Data_Figure_3; do
  cp "$B/$n.pdf" "$B/$n.png" "$D/"
done
echo "[OK] $D"

# Supplementary sheets (Supplementary Information; Nature wants additional figures as Extended Data, see the notes in
# figures_submission/README.txt). Same rules except the size: the sheets keep their own width.
python3 "$M"/make_supp_scatters.py --out-dir "$B" >/dev/null
(IN="$AMR_INPUT_DIR/mtb_wider_figure_inputs" WB="$AMR_INPUT_DIR/AMRPredict_panel_results_wider.xlsx" \
  OUT="$B/mtb_tmp" bash "$M"/make_all_figures.sh >/dev/null)
(IN="$R" TAG=reported1k WB="$AMR_INPUT_DIR/AMRPredict_pseaer_reported_results.xlsx" \
  MTB_WB="$AMR_INPUT_DIR/AMRPredict_panel_results_wider.xlsx" OUT="$B/pa_tmp" CROSS_OUT="$B/cross_tmp" \
  bash "$PA_DRIVER" >/dev/null)
for e in pdf png; do
  cp "$B/supplementary_figure1_scatters_mtb_draft.$e" "$D/Supplementary_Figure_1.$e"
  cp "$B/supplementary_figure2_scatters_pa_draft.$e" "$D/Supplementary_Figure_2.$e"
  cp "$B/mtb_tmp/sheet_shap.$e" "$D/Supplementary_Figure_5.$e"
  cp "$B/pa_tmp/pa_sheet_shap.$e" "$D/Supplementary_Figure_6.$e"
  cp "$B/mtb_tmp/sheet_classdist.$e" "$D/Supplementary_Figure_3.$e"
  cp "$B/pa_tmp/pa_sheet_classdist.$e" "$D/Supplementary_Figure_4.$e"
done
echo "[OK] supplementary sheets"
