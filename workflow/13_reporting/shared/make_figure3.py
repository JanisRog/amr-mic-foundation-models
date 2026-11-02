#!/usr/bin/env python3
"""Figure 3 (draft): the resistance call of the reported model against the established tool.

Two rows, one per organism, three held-out measures each: very major error rate (resistant
called susceptible), major error rate (susceptible called resistant) and F1 for the resistant
class. One row per antimicrobial, two marks joined by a line, the difference at the right
(model minus tool, positive = better for the model):

  M. tuberculosis  the reported model (Evo 2 embeddings + TB-Profiler's call) against TB-Profiler alone
  P. aeruginosa    the reported model (Evo 2 embeddings + ARDaP's call, or + AMRFinderPlus presence
                   where ARDaP gives no call) against ARDaP alone, intermediate counted as
                   resistant, for the eight antimicrobials ARDaP calls

The reported model contains the tool's own call as a feature, so the comparison asks whether the
embeddings improve on the tool, not whether they replace it. Values are read from the results
workbooks (the P. aeruginosa one is built from the staged reported model) and written to
figure3_values.tsv beside the figure.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import INK, INK2, ITALIC_BOLD, header, load, pick, save  # noqa: E402
from make_provisional_figures import A_BASE, A_EMB, dumbbell, paired_summary  # noqa: E402
from plot_metric_bars import PANELS  # noqa: E402

HOME = pathlib.Path.home()
# F1 first: the overall verdict, then the two error types it is made of
PANEL_SPEC = (("SR_f1", "F1, resistant class", "%; higher is better", False),
              ("vMER", "Very major error rate", "%; resistant called susceptible; lower is better", True),
              ("MER", "Major error rate", "%; susceptible called resistant; lower is better", True))
ORG = {"mtb": "M. tuberculosis", "pseaer": "P. aeruginosa"}


def pval(p: float) -> str:
    return f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-workbook", default=str(HOME / "Desktop/AMRPredict_panel_results_wider.xlsx"))
    ap.add_argument("--pa-workbook", default=str(HOME / "Desktop/AMRPredict_pseaer_reported_results.xlsx"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="figure3_draft")
    ap.add_argument("--width", type=float, default=7.8)
    a = ap.parse_args()
    books = {"mtb": load(pathlib.Path(a.mtb_workbook))[0], "pseaer": load(pathlib.Path(a.pa_workbook))[0]}

    rows = {}
    for k in ("mtb", "pseaer"):
        P = PANELS[k]
        drugs = [d for d in P["drugs"] if pick(books[k], d, P["baseline"], "holdout", "SR_f1") is not None]
        rows[k] = (drugs, [d for d in drugs if d in P["underpowered"]])

    fig = plt.figure(figsize=(a.width, 6.6))
    gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.50,
                          height_ratios=[len(rows["mtb"][0]), len(rows["pseaer"][0])])
    values, first = [], {}
    for r, k in enumerate(("mtb", "pseaer")):
        P, (drugs, faint) = PANELS[k], rows[k]
        for n, (key, title, xl, low) in enumerate(PANEL_SPEC):
            ax = fig.add_subplot(gs[r, n])
            u = 100.0                                          # error rates and F1 in per cent
            x = [pick(books[k], d, "evo2", "holdout", key) * u for d in drugs]
            y = [pick(books[k], d, P["baseline"], "holdout", key) * u for d in drugs]
            dumbbell(ax, drugs, x, y, a_style=A_EMB, b_style=A_BASE, title=title, xlabel=xl,
                     better_low=low, faint=faint, true_delta=True, unit=u)
            wins, ties, m, med, p = paired_summary(drugs, x, y, low, exclude=faint, true_delta=True, unit=u)
            worse = m - wins - ties
            ax.annotate(f"model better {wins}, worse {worse}" + (f", equal {ties}" if ties else "")
                        + f" of {m} · P = {pval(p)}",
                        xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points",
                        ha="left", va="bottom", fontsize=5.4, color=INK2)
            ax.set_title(title, loc="left", color=INK, pad=13, fontweight="semibold", fontsize=7.8)
            ax.set_xlabel(xl, fontsize=6.0)
            if n:
                ax.set_yticklabels([])
            else:
                first[k] = ax
            ax.annotate("abcdef"[r * 3 + n], xy=(-0.06 if n else -0.62, 1.22), xycoords="axes fraction",
                        fontsize=10, fontweight="bold", color=INK, ha="right", va="bottom",
                        annotation_clip=False)
            print(f"[INFO] {ORG[k]} {title}: model better {wins}, worse {worse}, equal {ties} of {m}, "
                  f"median {med:+.4f}, P = {p:.4f}")
            for d, xv, yv in zip(drugs, x, y):
                values.append([ORG[k], d, key, round(xv / u, 4), round(yv / u, 4), d in faint])
    top = header(fig, "Figure 3 (draft) — the resistance call against the established tools",
                 "Held-out sets. The reported model, which contains the tool's own call as a feature, against the tool "
                 "alone: TB-Profiler in\nthe first row, ARDaP in the second (intermediate counted as resistant; the eight "
                 "antimicrobials it calls). Numbers at the right:\nmodel minus tool, in percentage points "
                 "(for the error rates lower is better, so a negative\ndifference favours the model); dark type favours the model, grey type the tool. "
                 "Counts and two-sided Wilcoxon tests are over the\nantimicrobials compared (first row: eight, bedaquiline and delamanid faded and "
                 "excluded; second row: eight).", pad_in=0.30)
    fig.subplots_adjust(top=top - 0.050, bottom=0.095, left=0.165, right=0.915)
    for k in first:
        pos = first[k].get_position()
        fig.text(0.022, (pos.y0 + pos.y1) / 2, ORG[k], rotation=90, fontsize=9, fontproperties=ITALIC_BOLD,
                 color=INK, ha="center", va="center")
    fig.legend(handles=[Line2D([], [], label="Reported model (embeddings + the tool's call)", **A_EMB),
                        Line2D([], [], label="Tool alone (TB-Profiler, ARDaP)", **A_BASE)],
               loc="lower center", ncol=2, fontsize=7.0, handletextpad=0.5, columnspacing=2.2,
               bbox_to_anchor=(0.5, 0.0))
    out = pathlib.Path(a.out_dir)
    save(fig, out, a.stem)
    with open(out / "figure3_values.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["organism", "antimicrobial", "measure", "reported_model", "tool_alone", "excluded_from_tests"])
        w.writerows(values)
    return 0


if __name__ == "__main__":
    sys.exit(main())
