#!/usr/bin/env python3
"""The three foundation models side by side, on one figure per organism.

Two panels sharing an x axis of drugs: essential agreement above, S/R F1 below.
Three points per drug, one per model, with a connector spanning the spread.

**Dots, not bars.** Every value sits between roughly 0.5 and 1.0, so a bar
anchored at zero spends nine tenths of its length on range no reader is asking
about and makes a 0.01 difference invisible. The y axis here starts near the
data, which is legitimate for points and would be misleading for bars.

**The message is usually the spread, not the winner.** MTB section 9 found the
three models indistinguishable at matched capacity; the connector is what carries
that, so the subtitle states the median spread as a measured number rather than
leaving the reader to eyeball it.

Values come from the results workbook, never from a hand-typed table, so a figure
that disagrees with the protocol means the workbook moved.
"""

from __future__ import annotations

import argparse
import pathlib
import statistics
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import (AXIS, GRID, INK, INK2, MODEL_LABEL, MODELS,
                                SERIES, SURFACE, header, load, pick, save, tidy)

# Same drug order as every other sheet in the set, so a reader moving between
# figures is not re-learning the x axis.
PANELS = {
    "mtb": ["isoniazid", "rifampicin", "ethambutol", "ethionamide", "kanamycin",
            "amikacin", "levofloxacin", "moxifloxacin", "bedaquiline", "delamanid"],
    "pseaer": ["ciprofloxacin", "levofloxacin", "tobramycin", "imipenem",
               "meropenem", "ceftazidime", "cefepime", "aztreonam", "amikacin",
               "piperacillin-tazobactam"],
}
METRICS = [("OT", "Essential agreement (%)",
            "within one doubling dilution of the true MIC"),
           ("SR_f1", "S/R F1 (%)",
            "F1 on the resistant class, at the clinical breakpoint")]


def series(summary, drugs, field, split="holdout"):
    """{model: [value per drug]}, with None where a run is absent."""
    return {m: [None if (v := pick(summary, d, m, split, field)) is None else 100 * v for d in drugs]
            for m in MODELS}


def draw_panel(ax, drugs, vals, *, title, note, show_xticks):
    tidy(ax, x=False, y=True)
    xs = np.arange(len(drugs))

    for i, d in enumerate(drugs):
        got = [vals[m][i] for m in MODELS if vals[m][i] is not None]
        if len(got) > 1:
            # The spread connector is the point of the figure: it says at a
            # glance whether the models differ on this drug at all.
            ax.plot([i, i], [min(got), max(got)], color=GRID, lw=2.0,
                    solid_capstyle="round", zorder=1)
    for m in MODELS:
        ys = vals[m]
        keep = [(x, y) for x, y in zip(xs, ys) if y is not None]
        if not keep:
            continue
        kx, ky = zip(*keep)
        # A surface-coloured ring so coincident points stay countable.
        ax.scatter(kx, ky, s=62, color=SERIES[m], zorder=3,
                   edgecolors=SURFACE, linewidths=1.4)

    ax.set_xlim(-0.6, len(drugs) - 0.4)
    ax.set_xticks(xs)
    if show_xticks:
        ax.set_xticklabels([d.replace("piperacillin-tazobactam", "pip-tazo")
                            for d in drugs], rotation=30, ha="right", fontsize=8.5)
    else:
        ax.set_xticklabels([])
    ax.tick_params(axis="y", labelsize=8.5)
    # Title and note stacked in axes coordinates rather than via set_title: the
    # title's own pad is measured from the axes, so a note placed just above the
    # axes lands underneath it.
    ax.text(0, 1.085, title, transform=ax.transAxes, fontsize=10, color=INK,
            va="bottom", ha="left", fontweight="semibold")
    ax.text(0, 1.015, note, transform=ax.transAxes, fontsize=7.5, color=INK2,
            va="bottom", ha="left")
    for s in ("left",):
        ax.spines[s].set_visible(True)
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(0.8)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workbook", required=True)
    ap.add_argument("--panel", choices=sorted(PANELS), required=True)
    ap.add_argument("--split", default="holdout")
    ap.add_argument("--label", default="", help="organism name for the title")
    ap.add_argument("--stem", default="model_comparison")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    summary, _ = load(pathlib.Path(args.workbook))
    drugs = PANELS[args.panel]

    missing = [(d, m) for d in drugs for m in MODELS
               if pick(summary, d, m, args.split, "OT") is None]
    if missing:
        print(f"[WARN] {len(missing)} drug/model pair(s) absent from the "
              f"{args.split} split: {missing[:6]}")

    # A dedicated row for the legend. Placed inside a panel it sits on the data
    # for some drug orderings and not others, which is a bug that only appears
    # for one organism.
    fig = plt.figure(figsize=(9.4, 7.0))
    gs = fig.add_gridspec(2, 1, hspace=0.42)
    axes = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
    spreads = []
    for ax, (field, title, note) in zip(axes, METRICS):
        vals = series(summary, drugs, field, args.split)
        draw_panel(ax, drugs, vals, title=title, note=note,
                   show_xticks=ax is axes[-1])
        for i in range(len(drugs)):
            got = [vals[m][i] for m in MODELS if vals[m][i] is not None]
            if len(got) > 1:
                spreads.append(max(got) - min(got))

    # Figure-level and below everything: a legend row inside the gridspec sits
    # on top of the rotated drug labels, which occupy space the gridspec does
    # not know about.
    fig.legend(handles=[Line2D([], [], marker="o", ls="none", ms=8,
                               color=SERIES[m], markeredgecolor=SURFACE,
                               markeredgewidth=1.2, label=MODEL_LABEL[m])
                        for m in MODELS],
               loc="lower center", bbox_to_anchor=(0.5, 0.012), ncol=3,
               fontsize=9, handletextpad=0.6, columnspacing=2.6,
               borderpad=0.0, frameon=False)

    med = statistics.median(spreads) if spreads else float("nan")
    worst = max(spreads) if spreads else float("nan")
    org = args.label or args.panel
    top = header(fig, f"Three foundation models, {org}",
                 f"Held-out isolates, one point per model per drug, the reported "
                 f"config in each case. The bar behind each triple is the spread "
                 f"between\nbest and worst model: median {med:.1f} percentage points across the "
                 f"{len(spreads)} drug-metric pairs, largest {worst:.1f}. Values "
                 f"are read from the results workbook.", pad_in=0.18)
    fig.subplots_adjust(top=top, bottom=0.175, left=0.085, right=0.985)
    save(fig, pathlib.Path(args.out_dir), args.stem)
    print(f"[INFO] spread median {med:.3f} points  max {worst:.3f}  "
          f"over {len(spreads)} drug-metric pairs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
