#!/usr/bin/env python3
"""S/R F1, reported model against its baseline, both organisms on one figure.

A box over eight drugs hides more than it shows, so the box is the summary and
the **paired points are the data**: one line per drug joining its model F1 to its
baseline F1. The pairing is the whole design — the same drug, the same isolates,
the same breakpoint — and a pair of unpaired boxes would throw that away and
invite a two-sample test that does not apply.

Drugs the baseline does not cover are absent rather than imputed: PA's ARDaP has
no aztreonam or levofloxacin column, and a missing comparator is not a score of
zero. Bedaquiline and delamanid are drawn but excluded from the box and the test,
which is what the MTB panel table does — too few resistant isolates for a verdict,
and two underpowered drugs should not move a summary over eight.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import (AXIS, GRID, INK, INK2, MODEL_LABEL, MODELS,
                                MUTED, SERIES, SURFACE, header, load, pick,
                                save, tidy)
from plot_model_comparison import PANELS

MTB_BASELINE = "TB-Profiler 6.3.0"
UNDERPOWERED = {"bedaquiline", "delamanid"}


def cost(v, m):
    return None if v is None or m is None else 2 * v + m


def reported_model(summary, drug, declared):
    if declared:
        return declared
    best, bc = None, float("inf")
    for m in MODELS:
        c = cost(pick(summary, drug, m, "development_cv", "vMER"),
                 pick(summary, drug, m, "development_cv", "MER"))
        if c is not None and c < bc:
            best, bc = m, c
    return best


def ardap_f1(base_dir: pathlib.Path, drug: str, label: str):
    f = base_dir / f"cv_{drug}" / "baseline_rows.tsv"
    if not f.exists():
        return None
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["split"] == "holdout" and r["label"] == label:
            v = r.get("SR_f1")
            return float(v) if v not in ("", None) else None
    return None


def collect(workbook, panel, declared, base_dir, base_label):
    """[(drug, model, model_f1, baseline_f1)] for drugs the baseline covers."""
    summary, _ = load(pathlib.Path(workbook))
    out = []
    for drug in PANELS[panel]:
        m = reported_model(summary, drug, declared)
        if m is None:
            continue
        mf = pick(summary, drug, m, "holdout", "SR_f1")
        bf = (ardap_f1(pathlib.Path(base_dir), drug, base_label) if base_dir
              else pick(summary, drug, MTB_BASELINE, "holdout", "SR_f1"))
        if mf is None or bf is None:
            continue
        out.append((drug, m, mf, bf))
    return out


def draw(ax, rows, *, title, baseline_name):
    counted = [r for r in rows if r[0] not in UNDERPOWERED]
    tidy(ax, x=False, y=True)
    xs = (0.0, 1.0)

    # Box over the counted drugs only; the excluded pair is drawn, not summarised.
    for x, vals in zip(xs, ([r[2] for r in counted], [r[3] for r in counted])):
        bp = ax.boxplot([vals], positions=[x], widths=0.42, showfliers=False,
                        medianprops=dict(color=INK, lw=1.8),
                        boxprops=dict(color=AXIS, lw=1.0),
                        whiskerprops=dict(color=AXIS, lw=1.0),
                        capprops=dict(color=AXIS, lw=1.0))
        for patch in bp["boxes"]:
            patch.set_zorder(1)

    rng = np.random.default_rng(7)
    for drug, m, mf, bf in rows:
        weak = drug in UNDERPOWERED
        j = rng.uniform(-0.055, 0.055)
        ax.plot([xs[0] + j, xs[1] + j], [mf, bf], color=GRID if not weak else "#efeee8",
                lw=1.4, ls="-" if not weak else (0, (2, 2)), zorder=2)
        ax.plot(xs[0] + j, mf, "o", ms=7, mfc=SERIES[m], mec=SURFACE, mew=1.1,
                alpha=1.0 if not weak else 0.45, zorder=3)
        ax.plot(xs[1] + j, bf, "o", ms=7, mfc="none", mec=INK2, mew=1.4,
                alpha=1.0 if not weak else 0.45, zorder=3)

    ax.set_xticks(xs)
    ax.set_xticklabels(["reported model", baseline_name], fontsize=9)
    ax.set_xlim(-0.55, 1.55)
    ax.tick_params(axis="y", labelsize=8.5)
    ax.set_title(title, loc="left", fontsize=10.5, color=INK, pad=26,
                 fontweight="semibold")
    for sp in ("left",):
        ax.spines[sp].set_visible(True)
        ax.spines[sp].set_color(AXIS)
        ax.spines[sp].set_linewidth(0.8)
    return counted


def describe(counted):
    """Two short lines; one long one collides with the neighbouring panel."""
    mf = [r[2] for r in counted]
    bf = [r[3] for r in counted]
    # A tie must not count as a win: moxifloxacin's model and TB-Profiler make the
    # same calls, and their F1 differ only by how the workbook rounded each (5e-7).
    higher = sum(1 for a, b in zip(mf, bf) if a - b > 1e-4)
    equal = sum(1 for a, b in zip(mf, bf) if abs(a - b) <= 1e-4)
    med = statistics.median(a - b for a, b in zip(mf, bf))
    tail = ""
    try:
        from scipy.stats import wilcoxon
        if len(counted) >= 6:
            tail = f", Wilcoxon p = {wilcoxon(mf, bf).pvalue:.3f}"
    except Exception:
        pass
    return (f"{len(counted)} drugs, model higher on {higher}"
            + (f", equal on {equal}" if equal else "") + "\n"
            f"median difference {med:+.3f}{tail}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-workbook", required=True)
    ap.add_argument("--pa-workbook", required=True)
    ap.add_argument("--pa-baseline-dir", required=True)
    ap.add_argument("--pa-baseline-label", default="ARDaP (I=R)")
    ap.add_argument("--stem", default="f1_model_vs_baseline")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    mtb = collect(args.mtb_workbook, "mtb", "evo2", "", "")
    # Evo 2 for both organisms: one designated model throughout (2026-10-13).
    pa = collect(args.pa_workbook, "pseaer", "evo2", args.pa_baseline_dir,
                 args.pa_baseline_label)

    # Shared y: the subtitle invites a comparison between the panels, and two
    # independent axes would make that comparison a reading error.
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 5.0), sharey=True)
    c1 = draw(axes[0], mtb, title="M. tuberculosis",
              baseline_name="TB-Profiler")
    c2 = draw(axes[1], pa, title="P. aeruginosa", baseline_name="ARDaP")
    for ax in axes:
        ax.set_ylabel("")
    axes[0].set_ylabel("S/R F1 (resistant class)", fontsize=9, color=INK2)

    for ax, c in zip(axes, (c1, c2)):
        ax.text(0.0, 1.005, describe(c), transform=ax.transAxes, ha="left",
                va="bottom", fontsize=7.6, color=INK2)

    # Only the models actually reported, and one neutral entry for the baseline:
    # the panels have different comparators, and naming just one of them in a
    # shared legend would label ARDaP's points "TB-Profiler".
    from matplotlib.lines import Line2D
    used = [m for m in MODELS if any(r[1] == m for r in mtb + pa)]
    handles = [Line2D([], [], marker="o", ls="none", ms=6, mfc=SERIES[m],
                      mec=SURFACE, mew=1.0, label=MODEL_LABEL[m]) for m in used]
    handles.append(Line2D([], [], marker="o", ls="none", ms=6.5, mfc="none",
                          mec=INK2, mew=1.4,
                          label="baseline (TB-Profiler / ARDaP)"))
    fig.legend(handles=handles, ncol=len(handles), loc="lower center",
               fontsize=8, handletextpad=0.5, columnspacing=1.8,
               bbox_to_anchor=(0.5, 0.012), frameon=False)

    # pad_in has to clear the panel titles AND the two-line stat block above each
    # axes, neither of which header() knows about.
    top = header(fig, "S/R F1: the model against its baseline",
                 "Held-out isolates, one line per drug joining its reported model "
                 "to its baseline. The comparison is paired, so the lines\nare the "
                 "data and the box only summarises them. Uncovered drugs are "
                 "absent rather than zeroed; bedaquiline and\ndelamanid are drawn "
                 "faintly and excluded from both summaries for want of resistant "
                 "isolates.", pad_in=0.66)
    fig.subplots_adjust(top=top, bottom=0.135, left=0.095, right=0.985, wspace=0.10)
    save(fig, pathlib.Path(args.out_dir), args.stem)
    print(f"[INFO] MTB  {describe(c1)}")
    print(f"[INFO] PA   {describe(c2)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
