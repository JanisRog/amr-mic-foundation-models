#!/usr/bin/env python3
"""Supplementary Figures 1 and 2 (draft): predicted against measured MIC for every antimicrobial.

One sheet per organism, the reported model on the held-out set, drawn with the same functions as
Figure 1d,e (plot_mic_scatter.draw_grid). Every panel is drawn from its own saved predictions and only
after the metrics recomputed from them reproduce the antimicrobial's row in the results workbook
(make_panel_grids.load_drug). Rates and accuracies are given in per cent, as in the main figures.

  Supplementary Fig. 1  M. tuberculosis: Evo 2 embeddings of the target loci + TB-Profiler call
  Supplementary Fig. 2  P. aeruginosa: Evo 2 embeddings + ARDaP call (levofloxacin and aztreonam,
                        which ARDaP does not call: + AMRFinderPlus presence)
"""

from __future__ import annotations

import argparse
import os
import pathlib
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_figure2 import MTB_BOUND, PA_BOUND  # noqa: E402
from make_panel_figures import INK, ITALIC_BOLD, MUTED, header, load, save  # noqa: E402
from make_panel_grids import load_drug  # noqa: E402
from make_provisional_figures import NAME  # noqa: E402
from plot_mic_scatter import CATS, EDGE, FACE, LABEL, SIZE_BINS, draw_grid  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
DAGGER = {"levofloxacin", "aztreonam"}
SHEETS = (
    ("mtb", 1, "M. tuberculosis", MTB_BOUND, "AMRPredict_panel_results_wider.xlsx",
     "mtb_wider_figure_inputs", "mtb_wider_figure_inputs",
     "Evo 2 embeddings of the target loci with the TB-Profiler call as one further feature"),
    ("pseaer", 2, "P. aeruginosa", PA_BOUND, "AMRPredict_pseaer_reported_results.xlsx",
     "pa_reported_20261021/predictions_strict", "pa_reported_20261021/slim_strict",
     "Evo 2 embeddings of the target loci with the ARDaP call as one further feature\n(†, levofloxacin and "
     "aztreonam, which ARDaP does not call: with AMRFinderPlus gene presence)"),
)


def sheet(k, num, org, bounds, workbook, pred_dir, fold_dir, features, out):
    summary, _ = load(HOME / workbook)
    fig = plt.figure(figsize=(11.4, 15.2))
    gs = fig.add_gridspec(4, 3, hspace=0.42, wspace=0.28)
    for i, (drug, b) in enumerate(bounds.items()):
        y, pred, cat, m, _ = load_drug(HOME / pred_dir, HOME / fold_dir, drug, "evo2", b, summary, 5e-4)
        ax = fig.add_subplot(gs[i // 3, i % 3])
        draw_grid(ax, y, pred, cat, b, counts=False, fs=0.62)
        ax.tick_params(labelsize=6.5)
        name = NAME.get(drug, drug.capitalize()).replace("Pip.–tazobactam", "Piperacillin–tazobactam")
        ax.set_title(name + (" †" if k == "pseaer" and drug in DAGGER else ""), loc="left", fontsize=9.5,
                     color=INK, pad=24, fontweight="semibold")
        ax.annotate(f"n = {m['n']:,} ({m['n_resistant']} resistant)  ·  essential agreement {100 * m['OT']:.1f}%  ·  "
                    f"exact {100 * m['exact']:.1f}%\nvMER {100 * m['vMER']:.1f}%  ·  MER {100 * m['MER']:.1f}%  ·  "
                    f"F1 {100 * m['SR_f1']:.1f}%",
                    xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 4), textcoords="offset points",
                    fontsize=6.6, color=MUTED, linespacing=1.4)
        ax.set_xlabel("Measured MIC (log₂)", fontsize=8)
        ax.set_ylabel("Predicted MIC (log₂)", fontsize=8)
        print(f"[INFO] {org} {drug}: essential agreement {100 * m['OT']:.1f}%, n = {m['n']}")

    key = fig.add_subplot(gs[3, 1:])
    key.axis("off")
    handles = [Line2D([], [], marker="o", ls="none", ms=8, mfc=FACE[c], mec=EDGE[c], mew=0.9, label=LABEL[c])
               for c, _, _, _ in CATS]
    handles += [Patch(facecolor="#cde2fb", alpha=0.35, edgecolor="#9ec5f4", lw=0.8, label="± one dilution zone"),
                Line2D([], [], color="#d03b3b", lw=1.1, ls=(0, (1.6, 1.6)), label="S/R breakpoint")]
    first = key.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.02, 1.0), ncol=1, fontsize=9,
                       handletextpad=0.8, labelspacing=1.0, borderpad=0.0)
    key.add_artist(first)
    key.legend(handles=[Line2D([], [], marker="o", ls="none", ms=(sz * 0.62 ** 2) ** 0.5 * 0.85,
                               mfc="#c9c8c1", mec="#8a8980", mew=0.7, label=lab) for _, _, lab, sz in SIZE_BINS],
               loc="upper left", bbox_to_anchor=(0.45, 1.0), ncol=1, fontsize=9, handletextpad=1.4,
               labelspacing=1.35, borderpad=0.0, title="isolates in cell", title_fontsize=8)

    top = header(fig, f"Supplementary Figure {num} (draft) — predicted against measured MIC",
                 f"Held-out sets, every antimicrobial. {features}.\nBubble size is binned by the number "
                 "of isolates in the cell. vMER, very major error rate (resistant called susceptible); MER, major "
                 "error rate (susceptible called resistant); F1 for the resistant class.",
                 pad_in=0.30 + (0.28 if os.environ.get("FIG_SUBMISSION") == "1" else 0.0))
    fig.subplots_adjust(top=top - 0.010, bottom=0.035, left=0.065, right=0.985)
    if fig.texts:
        title = fig.texts[0]                 # the organism, in italics, at the right end of the title line
        fig.text(0.985, title.get_position()[1], org, ha="right", va=title.get_va(), fontsize=title.get_fontsize(),
                 fontproperties=ITALIC_BOLD, color=INK)
    else:                                    # submission figures have no title line: the organism heads the sheet
        fig.text(0.985, 1 - 0.12 / fig.get_figheight(), org, ha="right", va="top", fontsize=8.5,
                 fontproperties=ITALIC_BOLD, color=INK)
    save(fig, out, f"supplementary_figure{num}_scatters_{'mtb' if k == 'mtb' else 'pa'}_draft")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    for s in SHEETS:
        sheet(*s, pathlib.Path(a.out_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
