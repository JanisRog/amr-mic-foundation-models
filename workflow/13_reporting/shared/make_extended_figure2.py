#!/usr/bin/env python3
"""Extended Data Figure 2 (draft): accuracy of the reported model by MIC class, both organisms.

For each antimicrobial and each true log2 MIC class, the share of the held-out isolates in that class whose MIC the
reported model predicts within one doubling dilution. The number of held-out isolates of the class is printed under the value. A cross marks a class with fewer than 30 isolates, blank cells
are classes with no isolates. This is the per-class view behind the headline essential agreement: the headline is carried by the
common MIC classes, and accuracy falls in the sparse classes and in the gap between the two peaks of a bimodal
distribution, which is why balanced essential agreement (each class weighted equally) is the primary measure.

Reads the per-class sheets of the results workbooks (M. tuberculosis: wider target set; P. aeruginosa: reported
model). Writes extended_figure2_values.tsv beside the figure.
"""

from __future__ import annotations

import argparse
import csv
import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import SUBMISSION, fig_h, INK, INK2, ITALIC_BOLD, MUTED, ORG_COLOUR, SURFACE, header, load, save  # noqa: E402
from plot_metric_bars import PANELS  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
ORG = {"mtb": "M. tuberculosis", "pseaer": "P. aeruginosa"}
NAME = {"piperacillin-tazobactam": "Pip.–tazobactam"}
MIN_N = 30


def ramp(colour):
    """White to the organism colour: low accuracy pale, high accuracy dark."""
    return LinearSegmentedColormap.from_list("org", ["#f6f7f9", colour])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-workbook", default=str(HOME / "AMRPredict_panel_results_wider.xlsx"))
    ap.add_argument("--pa-workbook", default=str(HOME / "AMRPredict_pseaer_reported_results.xlsx"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="extended_figure2_draft")
    ap.add_argument("--width", type=float, default=7.8)
    a = ap.parse_args()
    books = {"mtb": load(pathlib.Path(a.mtb_workbook))[1], "pseaer": load(pathlib.Path(a.pa_workbook))[1]}

    cells, ns, ks_all = {}, {}, {}
    for k in ("mtb", "pseaer"):
        for r in books[k]:
            if r["model"] == "evo2" and r["split"] == "holdout":
                key = (r["drug"], int(float(r["log2_mic"])))
                cells[(k, *key)] = 100 * float(r["OT"])
                ns[(k, *key)] = float(r["n"])
        ks_all[k] = sorted({q for (kk, _, q) in cells if kk == k})
    fig = plt.figure(figsize=(a.width, fig_h(7.6)))
    gs = fig.add_gridspec(2, 1, hspace=0.50)
    axes, ims = {}, {}
    for r_, k in enumerate(("mtb", "pseaer")):
        ax = fig.add_subplot(gs[r_, 0])
        drugs = PANELS[k]["drugs"]
        ks = list(range(ks_all[k][0], ks_all[k][-1] + 1))
        grid = np.full((len(drugs), len(ks)), np.nan)
        for i, d in enumerate(drugs):
            for j, q in enumerate(ks):
                if (k, d, q) in cells:
                    grid[i, j] = cells[(k, d, q)]
        im = ax.imshow(grid, cmap=ramp(ORG_COLOUR[k]), vmin=0, vmax=100, aspect="auto")
        for i, d in enumerate(drugs):
            for j, q in enumerate(ks):
                if (k, d, q) not in cells:
                    continue
                n = ns[(k, d, q)]
                v = cells[(k, d, q)]
                col_ = SURFACE if v > 62 else INK
                ax.text(j, i - 0.17, f"{v:.0f}", ha="center", va="center", fontsize=6.0, color=col_, zorder=3)
                ax.text(j, i + 0.24, f"{int(n):,}", ha="center", va="center", fontsize=5.0, color=col_, alpha=0.8,
                        zorder=3)
                if n < MIN_N:                       # a sparse class: the cell is outlined in grey
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="#b5b3a9", lw=0.9,
                                               ls=(0, (2, 2)), zorder=4))
        ax.set_xticks(range(len(ks)))
        ax.set_xticklabels([str(q) for q in ks], fontsize=6.6)
        ax.set_yticks(range(len(drugs)))
        ax.set_yticklabels([NAME.get(d, d.capitalize()) for d in drugs], fontsize=7.0, color=INK)
        ax.set_xlabel("true log₂ MIC class", fontsize=7.0)
        ax.tick_params(length=0)
        for s_ in ax.spines.values():
            s_.set_visible(False)
        axes[k], ims[k] = ax, im
    top = header(fig, "Extended Data Figure 2 (draft) — accuracy of the reported model by MIC class",
                 "Held-out sets. Share of the isolates of each true log₂ MIC class (column) whose MIC is predicted within one\n"
                 "doubling dilution, per antimicrobial (row); number in the cell, per cent, and below it the number of isolates. Dashed outline: a class with fewer\n"
                 "than 30 isolates; blank: no isolates in the class.", pad_in=0.30)
    fig.subplots_adjust(top=top - 0.025, bottom=0.07, left=0.17, right=0.93)
    fig.canvas.draw()
    for k, lt in (("mtb", "a"), ("pseaer", "b")):
        pos = axes[k].get_position()
        fig.text((pos.x0 + pos.x1) / 2, pos.y1 + 0.012, ORG[k], ha="center", va="bottom", fontsize=9.6,
                 fontproperties=ITALIC_BOLD, color=INK)
        fig.text(0.02, pos.y1 + 0.012, lt, fontsize=10, fontweight="bold", color=INK, ha="left", va="bottom")
        cax = fig.add_axes([pos.x1 + 0.012, pos.y0 + 0.2 * pos.height, 0.011, 0.6 * pos.height])
        cb = fig.colorbar(ims[k], cax=cax)
        cb.set_ticks([0, 25, 50, 75, 100])
        cb.set_ticklabels(["0", "25", "50", "75", "100"])
        cb.ax.tick_params(labelsize=6.0, length=0)
        cb.outline.set_visible(False)
        cb.set_label("essential agreement (%)", fontsize=6.4)
    out = pathlib.Path(a.out_dir)
    save(fig, out, a.stem)
    with open(out / "extended_figure2_values.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["organism", "antimicrobial", "log2_mic_class", "isolates", "within_one_dilution_pct"])
        for (k, d, q), v in sorted(cells.items()):
            w.writerow([ORG[k], d, q, int(ns[(k, d, q)]), round(v, 1)])
    sparse = sum(1 for key, n in ns.items() if n < MIN_N)
    print(f"[INFO] {len(cells)} drug-class cells, {sparse} with fewer than {MIN_N} isolates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
