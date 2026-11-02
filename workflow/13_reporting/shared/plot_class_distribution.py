#!/usr/bin/env python3
"""MIC class distribution with holdout accuracy, one drug per figure.

The scatter says where a prediction lands; this says how many isolates the class
had to begin with, and how the holdout ones were scored. Both are needed because
the panel's headline metrics are averages over very unequal classes -- isoniazid
runs from 2,126 isolates at log2 -4 to 102 at 4 -- and `bal_OT` weights those
equally while `OT` does not. A class where the model looks poor and a class with
forty isolates are different objects, and only one of them is a modelling
problem.

Three encodings share an x-axis:

  **grey bar**   development isolates in the class. Context, not a result.
  **stacked bar** the holdout, split by the same five error categories the
                  scatter uses -- same colours, so the two figures read as one
                  system.
  **red line**    essential agreement for the class, on a right-hand axis.

**The second axis is a deliberate exception.** Counts and accuracy have no common
scale, and a dual axis invites reading a correlation between two things that were
never plotted against each other. It earns its place here because the question
the figure answers is exactly "does accuracy track class size", and the
alternative -- a second panel below -- separates the two encodings the reader is
meant to compare at a glance. The line is drawn in a colour used for nothing else
in the figure, and its axis is labelled in that colour, so which mark belongs to
which scale is never ambiguous.

Every metric is recomputed from the predictions and checked against the workbook
row before anything is drawn, as in `plot_mic_scatter.py`.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import AXIS, GRID, INK, INK2, MUTED, SURFACE, header, load, save, tidy
from plot_mic_scatter import (CATS, EDGE, FACE, LABEL, MODEL_LABEL, categorise,
                              choose_aggregation, CHECK)

TRAIN_FACE, TRAIN_EDGE = "#dedcd4", "#b5b3a9"
OT_LINE = "#d03b3b"
STACK = ["exact", "ot", "ot_bad", "far_ok", "far_bad"]


def draw_classdist(ax, ks, tr, stacks, hold, ot, boundary, *, annotate=True,
                   fs=1.0):
    """Grey development bars, stacked holdout bars and the accuracy line.

    Shared by the single-drug figure and the all-drug sheet. Returns the twin
    axes so the caller decides whether to label it -- on a sheet only the
    outermost panels should carry axis text.
    """
    x = np.arange(len(ks))
    tidy(ax, x=False, y=True)
    w = 0.38
    ax.bar(x - w / 2, tr, w, color=TRAIN_FACE, edgecolor=TRAIN_EDGE,
           lw=0.7 * fs, zorder=2)
    bottom = np.zeros(len(ks))
    for c in STACK:
        v = np.array(stacks[c], dtype=float)
        ax.bar(x + w / 2, v, w, bottom=bottom, color=FACE[c], edgecolor=EDGE[c],
               lw=0.6 * fs, zorder=2)
        bottom += v
    if annotate:
        for xi, v in zip(x - w / 2, tr):
            if v:
                ax.annotate(f"{v:,}", (xi, v), xytext=(0, 3), ha="center",
                            textcoords="offset points", fontsize=6.6, color=MUTED)
        for xi, v in zip(x + w / 2, hold):
            if v:
                ax.annotate(f"{v:,}", (xi, v), xytext=(0, 3), ha="center",
                            textcoords="offset points", fontsize=6.6, color=INK2)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{k:g}" for k in ks])
    ax.set_ylim(0, max(max(tr), max(hold)) * 1.16)
    ax.set_xlim(-0.75, len(ks) - 0.25)

    ax2 = ax.twinx()
    ax2.plot(x, [100 * v for v in ot], color=OT_LINE, lw=1.6 * fs, marker="o", ms=5 * fs,
             mec=SURFACE, mew=0.8 * fs, zorder=4)
    ax2.set_ylim(-2, 105)
    ax2.set_yticks([0, 25, 50, 75, 100])
    ax2.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
    ax2.tick_params(axis="y", colors=OT_LINE, length=0, labelsize=7.5)
    for sp in ax2.spines.values():
        sp.set_visible(False)

    # The divider sits between the last susceptible class and the first
    # resistant one, for the reason the scatter's does: resistance is MIC
    # *above* the boundary, so the boundary's own class is susceptible.
    div = next((i for i, k in enumerate(ks) if k > boundary), None)
    if div is not None and 0 < div < len(ks):
        ax.axvline(div - 0.5, color=AXIS, lw=0.8, ls=(0, (3, 3)), zorder=1)
        ax.annotate("S  |  R", xy=(div - 0.5, 0.985),
                    xycoords=("data", "axes fraction"), xytext=(0, 0),
                    textcoords="offset points", ha="center", va="top",
                    fontsize=7, color=MUTED)
    return ax2


def class_counts(folds_path, y, cat, boundary):
    """Development counts per class, and the holdout stack beside them."""
    train = {}
    with open(folds_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if int(r["fold"]) < 0:
                k = float(r["log2_mic"])
                train[k] = train.get(k, 0) + 1
    # A holdout isolate in a class with no development isolates is not a quirk --
    # it means the fold file and the predictions describe different cohorts.
    missing = sorted(set(y) - set(train))
    if missing:
        print(f"[WARN] MIC classes in the holdout with no development "
              f"isolates: {missing}")
    ks = sorted(set(train) | set(y))
    tr = [train.get(k, 0) for k in ks]
    stacks = {c: [int(((y == k) & (cat == c)).sum()) for k in ks] for c in STACK}
    hold = [sum(stacks[c][i] for c in STACK) for i in range(len(ks))]
    ot = [(sum(stacks[c][i] for c in ("exact", "ot", "ot_bad")) / hold[i])
          if hold[i] else np.nan for i in range(len(ks))]
    return ks, tr, stacks, hold, ot, sum(train.values())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--folds", required=True,
                    help="cv_folds_final.tsv: fold < 0 development, >= 0 holdout")
    ap.add_argument("--workbook", default=str(pathlib.Path.home()
                                              / "Desktop/AMRPredict_panel_results.xlsx"))
    ap.add_argument("--drug", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--sr-boundary", type=float, required=True)
    ap.add_argument("--tolerance", type=float, default=5e-4)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="")
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
        raise SystemExit(f"[ERROR] no holdout row for {args.drug}/{args.model}")
    want = {k: (None if want.get(k) in ("", None) else float(want[k]))
            for k, _ in CHECK}
    agg, m = choose_aggregation(df, args.sr_boundary, want, args.tolerance)

    g = df[df["aggregation"] == agg]
    y = g["log2_mic"].to_numpy(float)
    pred = g["pred"].to_numpy(float)
    cat = categorise(y, pred, args.sr_boundary)

    ks, tr, stacks, hold, ot, n_train = class_counts(
        args.folds, y, cat, args.sr_boundary)

    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    ax2 = draw_classdist(ax, ks, tr, stacks, hold, ot, args.sr_boundary)
    ax.set_xlabel("log₂ MIC class")
    ax.set_ylabel("isolates")
    ax.set_ylim(0, max(max(tr), max(hold)) * 1.16)
    ax.set_xlim(-0.75, len(ks) - 0.25)

    ax2.set_ylabel("essential agreement (%)", color=OT_LINE)

    handles = [Patch(facecolor=TRAIN_FACE, edgecolor=TRAIN_EDGE, lw=0.7,
                     label=f"development  (n = {n_train:,})")]
    handles += [Patch(facecolor=FACE[c], edgecolor=EDGE[c], lw=0.6, label=LABEL[c])
                for c in STACK]
    handles += [Line2D([], [], color=OT_LINE, lw=1.6, marker="o", ms=5,
                       label="essential agreement by class (right axis)")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.075, 1.0),
              fontsize=8, handletextpad=0.7, labelspacing=0.75, borderpad=0.0)

    top = header(fig, f"{args.drug.capitalize()} — MIC class distribution and "
                      f"holdout accuracy",
                 f"{MODEL_LABEL.get(args.model, args.model)}. Grey is the "
                 f"development cohort, the stack is the held-out set split by "
                 f"the same\ncategories as the confusion grid. The line reads "
                 f"against the right-hand axis. Config {config} ({agg}).")
    fig.subplots_adjust(top=top, bottom=0.135, left=0.065, right=0.735)

    stem = args.stem or f"classdist_{args.drug}"
    save(fig, pathlib.Path(args.out_dir), stem)
    return 0


if __name__ == "__main__":
    sys.exit(main())
