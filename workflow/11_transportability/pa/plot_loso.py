#!/usr/bin/env python3
"""What the leave-one-study-out transport cost is made of.

Two panels against the same y axis: the cost of holding out a collection, plotted
against how much of the training set that collection was, and against how far its
resistance rate sits from the rest of the cohort.

The figure exists to separate two explanations that produce the same median. If
the model were predicting cohort membership and applying that collection's base
rate, the right panel would slope **up** -- the most dissimilar collections would
be the ones it fails on. It slopes down. The left panel is where the cost
actually comes from.

Points are not coloured by drug. Eight categories exceeds what a categorical
palette can separate honestly, and the claim is about the relationship rather
than about which drug sits where; the few points worth naming are labelled
directly instead.
"""

from __future__ import annotations

import argparse
import pathlib
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.environ.get("MTB_SCRIPTS") or str(pathlib.Path(__file__).resolve().parent.parent / "mtb_tbprofiler"))
from make_panel_figures import (AXIS, GRID, INK, INK2, MUTED, SURFACE, header,
                                save, tidy)

ACCENT = "#2a78d6"
WARM = "#eb6834"


def fit_line(ax, x, y):
    """Least-squares line with Pearson r and p, or nothing if scipy is absent."""
    try:
        from scipy import stats
        res = stats.linregress(x, y)
        r, p = res.rvalue, res.pvalue
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, res.intercept + res.slope * xs, color=INK2, lw=1.6,
                ls=(0, (5, 2)), zorder=2)
        return f"r = {r:+.3f},  p = {p:.3f}"
    except Exception:
        return ""


def panel(ax, d, xcol, *, title, xlabel, note, label_extremes=0, pct=True):
    tidy(ax, x=False, y=True)
    x = d[xcol].to_numpy(float)
    y = d["d_cost"].to_numpy(float)

    ax.axhline(0, color=AXIS, lw=1.0, ls=(0, (3, 3)), zorder=1)
    stat = fit_line(ax, x, y)

    # Area proportional to the isolates held out, so the large collections --
    # the ones driving the sample-size explanation -- are visible as large.
    size = 18 + 120 * (d["n"] / d["n"].max())
    ax.scatter(x, y, s=size, color=ACCENT, alpha=0.62, lw=0.8,
               edgecolors=SURFACE, zorder=3)

    if label_extremes:
        lo, hi = x.min(), x.max()
        for _, r in d.reindex(d.d_cost.abs().sort_values(ascending=False).index
                              ).head(label_extremes).iterrows():
            # Flip the label inward near the right edge, where it would otherwise
            # run off the axes and be silently clipped.
            right = (r[xcol] - lo) / max(1e-9, hi - lo) > 0.72
            ax.annotate(f"{r.drug.split('-')[0][:11]} · "
                        f"{r.study.split(' et')[0].split(' (')[0][:9]}",
                        (r[xcol], r.d_cost), textcoords="offset points",
                        xytext=(-7 if right else 7, 4),
                        ha="right" if right else "left",
                        fontsize=6.6, color=MUTED)

    ax.set_xlabel(xlabel, fontsize=8.5, color=INK2)
    ax.set_title(title, loc="left", fontsize=10, color=INK, pad=20,
                 fontweight="semibold")
    ax.text(0, 1.055, note, transform=ax.transAxes, fontsize=7.5, color=INK2,
            va="bottom", ha="left")
    if stat:
        ax.text(0.98, 0.04, stat, transform=ax.transAxes, fontsize=8,
                color=INK2, ha="right", va="bottom")
    if pct:
        ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    ax.tick_params(labelsize=8.5)
    for sp in ("left",):
        ax.spines[sp].set_visible(True)
        ax.spines[sp].set_color(AXIS)
        ax.spines[sp].set_linewidth(0.8)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="output of analyse_loso.py")
    ap.add_argument("--stem", default="loso_transport_cost")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    d = pd.read_csv(args.table, sep="\t")
    need = {"d_cost", "frac_removed", "rate_gap", "n", "drug", "study"}
    missing = need - set(d.columns)
    if missing:
        raise SystemExit(f"[ERROR] {args.table} lacks {sorted(missing)} — "
                         f"regenerate it with the current analyse_loso.py")

    big = d[d.frac_removed >= 0.30].d_cost.median()
    small = d[d.frac_removed < 0.30].d_cost.median()
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.6), sharey=True)
    # The linear trend on the left is weak (p ~ 0.12); what carries the claim is
    # the grouped comparison, so that is the number on the panel rather than the
    # correlation alone.
    try:
        from scipy import stats
        mw = stats.mannwhitneyu(d[d.frac_removed >= 0.30].d_cost,
                                d[d.frac_removed < 0.30].d_cost,
                                alternative="two-sided").pvalue
        left_note = (f"collections worth ≥30% of training cost "
                     f"{big:+.3f} against {small:+.3f} (p = {mw:.3f})")
    except Exception:
        left_note = "larger collections cost more to lose"
    panel(axes[0], d, "frac_removed",
          title="Cost rises with the data removed",
          xlabel="held-out collection as a share of the training set",
          note=left_note, label_extremes=3)
    panel(axes[1], d, "rate_gap",
          title="Not how different the collection is",
          xlabel="gap between the collection's resistance rate and the rest",
          note="the cohort-prediction hypothesis predicts the opposite slope",
          label_extremes=3)
    axes[0].set_ylabel("Δ cost  (leave-one-study-out − matched CV)",
                       fontsize=8.5, color=INK2)

    top = header(fig, "What the transport cost is made of",
                 f"Each point is one drug and one held-out collection "
                 f"({len(d)} pairs, {d.drug.nunique()} drugs); area is the number "
                 f"of isolates held out.\nAbove the dashed line the model is worse "
                 f"without that collection in training. Median Δcost "
                 f"{d.d_cost.median():+.3f} overall,\n{big:+.3f} for collections "
                 f"worth ≥30% of training against {small:+.3f} for the rest.",
                 pad_in=0.62)
    fig.subplots_adjust(top=top, bottom=0.135, left=0.088, right=0.985, wspace=0.09)
    save(fig, pathlib.Path(args.out_dir), args.stem)
    return 0


if __name__ == "__main__":
    sys.exit(main())
