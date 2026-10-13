#!/usr/bin/env python3
"""Predicted against actual MIC, as a bubble grid, for one drug and one model.

Both axes are discrete -- the truth is a doubling-dilution class and the
prediction is snapped to one -- so a point cloud would stack thousands of
isolates on a few dozen positions and show only the scatter's outline. One
bubble per (actual, predicted) cell, area proportional to count, count printed
inside: the same information a confusion matrix carries, laid out so the
diagonal and the distance from it are the things you see.

**The plot re-derives every metric from the predictions and checks them against
the results workbook before drawing.** The predictions come from a re-fit, and a
re-fit is only usable if it reproduces what was reported; the file also holds
both tree aggregations, and which of the two was declared is not recorded in the
workbook. Both problems have the same answer -- score each aggregation and keep
the one that matches the reported row. If neither matches, something in the
re-fit differed and the plot does not get drawn.

Colour is the error category, which is ordered, not nominal: exact, within one
dilution, then the two ways of being wrong that matter clinically. A cell that
crosses the S/R breakpoint is called out in status colour even when the MIC is
nearly right, because a one-dilution error at the breakpoint is a different
event from a one-dilution error in the middle of the susceptible mode.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import (AXIS, GRID, INK, INK2, MUTED, SURFACE, header,
                                load, save, tidy)

# ordered error categories, light -> dark on correctness, status hues for the
# two that cross the breakpoint
CATS = [
    ("exact",      "exact class",                     "#256abf", "#184f95"),
    ("ot",         "within one dilution",             "#9ec5f4", "#5598e7"),
    ("far_ok",     "wrong MIC, correct S/R call",     "#fab219", "#b07c00"),
    ("ot_bad",     "within one dilution, wrong S/R",  "#ec835a", "#b4522c"),
    ("far_bad",    "wrong MIC and wrong S/R",         "#d03b3b", "#8f2020"),
]
FACE = {k: f for k, _, f, _ in CATS}
EDGE = {k: e for k, _, _, e in CATS}
LABEL = {k: l for k, l, _, _ in CATS}

MODEL_LABEL = {"evo2": "Evo 2 7B", "carbon": "Carbon 8B", "dnabert": "DNABERT-2"}

# Six size bins rather than area proportional to count. A continuous area scale
# is dominated by the two susceptible modes -- 467 isolates in one cell against
# 1 in another is a 21-fold difference in radius, which renders the single-isolate
# cells as dots and hides exactly the corner errors worth looking at. Binning
# compresses that range while keeping the order readable. The cost is that a
# bubble no longer encodes its count precisely, so the key below is not optional.
SIZE_BINS = [(1, 1, "1", 30), (2, 5, "2–5", 105), (6, 20, "6–20", 190),
             (21, 50, "21–50", 340), (51, 200, "51–200", 560),
             (201, 10**9, ">200", 880)]


def bubble_size(n: int) -> float:
    for lo, hi, _, s in SIZE_BINS:
        if lo <= n <= hi:
            return s
    raise SystemExit(f"[ERROR] count {n} falls in no size bin")
CHECK = [("exact", "exact"), ("OT", "OT"), ("bal_OT", "bal_OT"),
         ("MAE", "MAE"), ("vMER", "vMER"), ("MER", "MER"), ("SR_f1", "F1")]


def balanced_recall(y, hit, min_n=0):
    vals, out = np.unique(y), []
    for c in vals:
        m = y == c
        if m.sum() >= min_n:
            out.append(float(np.mean(hit[m])))
    return float(np.mean(out)) if out else float("nan")


def metrics(y, pred, boundary):
    err = np.abs(y - pred)
    true_r, pred_r = y > boundary, pred > boundary
    n_r, n_s = int(true_r.sum()), int((~true_r).sum())
    tp = int((true_r & pred_r).sum())
    fp = int((~true_r & pred_r).sum())
    fn = int((true_r & ~pred_r).sum())
    denom = 2 * tp + fp + fn
    return {
        "exact": float(np.mean(err == 0)),
        "OT": float(np.mean(err <= 1)),
        "bal_OT": balanced_recall(y, err <= 1),
        "MAE": float(np.mean(err)),
        "vMER": float(fn / n_r) if n_r else float("nan"),
        "MER": float(fp / n_s) if n_s else float("nan"),
        # F1 on the resistant class, matching the workbook's SR_f1
        "SR_f1": float(2 * tp / denom) if denom else 0.0,
        "n": int(len(y)), "n_resistant": n_r,
    }


def categorise(y, pred, boundary):
    err = np.abs(y - pred)
    sr_ok = (y > boundary) == (pred > boundary)
    out = np.full(len(y), "far_bad", dtype=object)
    out[(err > 1) & sr_ok] = "far_ok"
    out[(err <= 1) & ~sr_ok] = "ot_bad"
    out[(err <= 1) & sr_ok] = "ot"
    out[(err == 0) & sr_ok] = "exact"
    return out


def choose_aggregation(df, boundary, want, tol):
    """Keep the aggregation that reproduces the reported row; report both."""
    rows = []
    for agg, g in df.groupby("aggregation"):
        m = metrics(g["log2_mic"].to_numpy(), g["pred"].to_numpy(), boundary)
        worst = max(abs(m[k] - want[k]) for k, _ in CHECK if want.get(k) is not None)
        rows.append((agg, m, worst))
        print(f"[INFO] {agg:<7} " + "  ".join(
            f"{lbl} {m[k]:.4f} (reported {want[k]:.4f})" for k, lbl in CHECK
            if want.get(k) is not None))
    rows.sort(key=lambda r: r[2])
    agg, m, worst = rows[0]
    if worst > tol:
        raise SystemExit(
            f"[ERROR] no aggregation reproduces the reported holdout row; "
            f"closest is {agg!r}, off by {worst:.4f} (tolerance {tol}).\n"
            f"        The re-fit differed from the scored run — check the cache, "
            f"the config file, the seed and --sr-boundary before trusting any "
            f"plot made from these predictions.")
    print(f"[OK] '{agg}' reproduces the reported row to {worst:.5f}")
    return agg, m


def draw_grid(ax, y, pred, cat, boundary, *, counts=True, fs=1.0):
    """The confusion grid itself, drawn into a given axes.

    Shared by the single-drug figure and the all-drug sheet, so the two cannot
    drift apart. `fs` scales the marks, since a panel a quarter the width needs
    smaller ones in the same proportions; `counts` is off on the sheet, where a
    bubble is too small to hold a number.
    """
    cells = {}
    for a, p_, c in zip(y, pred, cat):
        cells[(a, p_)] = (cells.get((a, p_), (0, c))[0] + 1, c)
    ticks = np.arange(min(y.min(), pred.min()), max(y.max(), pred.max()) + 1)
    lo, hi = ticks[0] - 0.7, ticks[-1] + 0.7

    tidy(ax, x=True, y=True)
    # A confusion grid reads as a matrix, so it keeps its frame -- unlike the
    # dumbbell and trade-off figures, where a box around the data adds nothing.
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_visible(True)
        ax.spines[sp].set_color(AXIS)
        ax.spines[sp].set_linewidth(0.8)
    ax.plot([lo, hi], [lo, hi], color=AXIS, lw=0.9, ls=(0, (3, 3)), zorder=1)
    ax.fill_between([lo, hi], [lo - 1, hi - 1], [lo + 1, hi + 1],
                    color="#cde2fb", alpha=0.35, lw=0, zorder=0)
    # An isolate is resistant when its MIC is *above* the boundary, so with
    # integer classes and a boundary of -3.0 the -3 class is susceptible and the
    # divider belongs between -3 and -2. Drawing it on -3.0, as the line's value
    # invites, puts a whole class on the wrong side of it.
    div = boundary + 0.5 if abs(boundary - round(boundary)) < 1e-9 else boundary
    for f in (ax.axvline, ax.axhline):
        f(div, color="#d03b3b", lw=0.8, ls=(0, (1.6, 1.6)), zorder=1)

    for (a, p_), (n, c) in sorted(cells.items(), key=lambda kv: -kv[1][0]):
        r = bubble_size(n) * fs ** 2
        ax.scatter(a, p_, s=r, facecolor=FACE[c], edgecolor=EDGE[c],
                   lw=0.8 * fs, zorder=3)
        if counts and r >= 105:
            t = 5.4 if r < 190 else (6.4 if r < 340 else (7.0 if r < 560 else 7.6))
            ax.annotate(f"{n}", (a, p_), ha="center", va="center", zorder=4,
                        fontsize=t,
                        color="#ffffff" if c in ("exact", "far_bad") else INK)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xticks(ticks); ax.set_yticks(ticks)
    ax.set_aspect("equal")
    return ticks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--workbook", default=str(pathlib.Path.home()
                                              / "Desktop/AMRPredict_panel_results.xlsx"))
    ap.add_argument("--drug", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--sr-boundary", type=float, required=True,
                    help="no default; the parameter that fails silently")
    ap.add_argument("--tolerance", type=float, default=5e-4)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="")
    ap.add_argument("--features-label", default="shared PCA4 + TB-Profiler rule bit",
                    help="what the model was fitted on, for the subtitle. The default is the "
                         "M. tuberculosis feature set; P. aeruginosa has no catalogue call and "
                         "must pass its own")
    args = ap.parse_args()

    df = pd.read_csv(args.predictions, sep="\t")
    if df["config"].nunique() != 1:
        raise SystemExit(f"[ERROR] {df['config'].nunique()} configs in the "
                         f"prediction file; pass one config's predictions")
    config = df["config"].iloc[0]

    summary, _ = load(pathlib.Path(args.workbook))
    want = next((r for r in summary if r["drug"] == args.drug
                 and r["model"] == args.model and r["split"] == "holdout"), None)
    if want is None:
        raise SystemExit(f"[ERROR] no holdout row for {args.drug}/{args.model} "
                         f"in the workbook")
    want = {k: (None if want.get(k) in ("", None) else float(want[k]))
            for k, _ in CHECK}

    agg, m = choose_aggregation(df, args.sr_boundary, want, args.tolerance)
    g = df[df["aggregation"] == agg]
    y = g["log2_mic"].to_numpy(float)
    pred = g["pred"].to_numpy(float)
    cat = categorise(y, pred, args.sr_boundary)

    # Two columns: the grid gets the left two thirds, the right third carries
    # the key above and the metrics below. Both used to sit inside the axes,
    # where they covered cells -- and on a confusion grid the covered cells are
    # the sparse corners, which is exactly what the figure is for.
    # Sized so the square grid fills its column: the axes are aspect-equal, so
    # too wide a figure leaves dead space between the plot and the key rather
    # than a bigger plot.
    fig = plt.figure(figsize=(8.8, 7.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[2, 1], wspace=0.03)
    ax = fig.add_subplot(gs[0])
    key_ax = fig.add_subplot(gs[1])
    key_ax.axis("off")
    draw_grid(ax, y, pred, cat, args.sr_boundary)
    bp_label = "S/R breakpoint"

    ax.set_xlabel("Actual  log₂ MIC")
    ax.set_ylabel("Predicted  log₂ MIC")
    ax.set_aspect("equal")

    box = "\n".join([f"n = {m['n']:,}   ({m['n_resistant']} resistant)",
                     f"exact = {m['exact']:.3f}", f"OT = {m['OT']:.3f}",
                     f"bal_OT = {m['bal_OT']:.3f}", f"MAE = {m['MAE']:.2f}",
                     f"vMER = {m['vMER']:.3f}", f"MER = {m['MER']:.3f}",
                     f"F1 = {m['SR_f1']:.3f}"])
    key_ax.annotate(box, xy=(0.02, 0.03), xycoords="axes fraction", va="bottom",
                    ha="left", fontsize=11.0, color=INK2, linespacing=1.95,
                    bbox=dict(boxstyle="round,pad=0.75", fc=SURFACE, ec=GRID,
                              lw=0.8))
    # Outside the axes: the corners of a confusion grid are exactly where the
    # interesting failures sit, and the first version put the legend on top of
    # ten isolates predicted six dilutions low.
    key_ax.legend(
        handles=[Line2D([], [], marker="o", ls="none", ms=8, mfc=FACE[k],
                        mec=EDGE[k], mew=0.9, label=LABEL[k])
                 for k, _, _, _ in CATS]
        # The band swatch carries the same fill and alpha as the band itself,
        # so it composites over the white figure to exactly the shade on the
        # plot; the hairline edge is only so a very pale patch is findable.
        + [Patch(facecolor="#cde2fb", alpha=0.35, edgecolor="#9ec5f4", lw=0.8,
                 label="± one dilution zone"),
           Line2D([], [], color="#d03b3b", lw=1.1, ls=(0, (1.6, 1.6)),
                  label=bp_label)],
        loc="upper left", bbox_to_anchor=(-0.02, 1.0), ncol=1, fontsize=10.0,
        handletextpad=0.8, labelspacing=1.25, borderpad=0.0)
    # No size key: every cell but the singletons prints its own count, which is
    # the exact number a key could only bracket. The one thing a reader cannot
    # read off the marks is what an unlabelled bubble means, so the subtitle says
    # it rather than a legend entry.

    top = header(fig, f"{args.drug.capitalize()} — predicted against actual MIC",
                 f"{MODEL_LABEL.get(args.model, args.model)}, held-out set, "
                 f"{args.features_label}. Bubble size is binned by "
                 f"the cell's isolate count, printed\ninside it; the smallest, "
                 f"unlabelled bubble is one isolate. Shaded band is the "
                 f"one-dilution zone.\nConfig {config} ({agg} of the trees).")
    fig.subplots_adjust(top=top, bottom=0.095, left=0.075, right=0.995)

    stem = args.stem or f"fig_scatter_{args.drug}_{args.model}"
    save(fig, pathlib.Path(args.out_dir), stem)
    return 0


if __name__ == "__main__":
    sys.exit(main())
