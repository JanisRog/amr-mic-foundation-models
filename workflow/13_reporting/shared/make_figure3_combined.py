#!/usr/bin/env python3
"""Figure 3, combined-plane variant (draft): both organisms on one error plane.

  a      very major (x) against major (y) error rate, square-root axes, with the commonly used targets for
         susceptibility testing as a green box; one arrow per antimicrobial from the tool (open) to the reported
         model (filled); circles M. tuberculosis (TB-Profiler), triangles P. aeruginosa (ARDaP, intermediate
         counted as resistant). Arrows blue when the model is closer to the targets, orange when farther, grey
         at the same distance (distance = excess over the targets, percentage points).
  b, c   F1 for the resistant class, model minus tool, M. tuberculosis and P. aeruginosa.

Same numbers and helpers as make_figure3_alt.py; bedaquiline and delamanid are left out.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator
from scipy.stats import wilcoxon

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_figure3_alt import (CLOSER, CODE, EMB, FARTHER, NAME, T_ME, T_VME, pval, spread_labels,  # noqa: E402
                              to_box)
from make_panel_figures import FINAL, SUBMISSION, fit_fonts, INK, INK2, ITALIC, ITALIC_BOLD, MUTED, SURFACE, header, load, pick, save, scale_letters, tidy  # noqa: E402
from make_provisional_figures import dumbbell, paired_summary  # noqa: E402
from plot_metric_bars import PANELS  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
ORG = {"mtb": "M. tuberculosis", "pseaer": "P. aeruginosa"}
SHAPE = {"mtb": "o", "pseaer": "D"}     # diamonds: a triangle reads as an arrowhead
TOOL = "#8a8880"
TICKS = [0, 1, 2, 5, 10, 20, 40, 80]
PREFER = {("mtb", "isoniazid"): (-1.0, 0.0), ("pseaer", "cefepime"): (-1.0, 0.0), ("pseaer", "meropenem"): (0.0, 1.0),
          ("pseaer", "tobramycin"): (-1.0, 0.0), ("mtb", "amikacin"): (0.0, 1.0),
          ("mtb", "moxifloxacin"): (0.0, 1.0), ("mtb", "kanamycin"): (1.0, 0.0)}     # label directions (x = major, y = very major)
from make_panel_figures import ORG_COLOUR as ORGCOL, ORG_PALE as ORGPALE  # noqa: E402  navy and teal
NOARROW = {("mtb", "levofloxacin")}                 # tool and model practically identical


HEAD_PT = 7.0        # length of the solid head section, points


def dashed_arrow(ax, x0, y0, x1, y1, color):
    """Dashed shaft, solid filled head whose tip touches the model's marker.

    Call after the layout is final: the head is a fixed length in points, so the split between shaft and
    head is computed in display coordinates. (One patch with a dashed line style dashes the head's outline.)
    """
    to_d, to_x = ax.transData.transform, ax.transData.inverted().transform
    (a0, b0), (a1, b1) = to_d((x0, y0)), to_d((x1, y1))
    n = ((a1 - a0) ** 2 + (b1 - b0) ** 2) ** 0.5
    if n < 1e-9:
        return
    k = (HEAD_PT + 2.2) * ax.figure.dpi / 72.0            # head plus the gap to the marker
    k = min(k, 0.7 * n)                                  # very short arrows: head takes most of the length
    hx, hy = to_x((a1 - (a1 - a0) * k / n, b1 - (b1 - b0) * k / n))
    ax.add_patch(FancyArrowPatch((x0, y0), (hx, hy), arrowstyle="-", color=color, lw=1.1, ls=(0, (3.0, 1.8)),
                                 shrinkA=3.2, shrinkB=0, zorder=2))
    ax.add_patch(FancyArrowPatch((hx, hy), (x1, y1), arrowstyle="-|>", mutation_scale=7.5, color=color, lw=1.1,
                                 shrinkA=0, shrinkB=1.6, zorder=2))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-workbook", default=str(HOME / "AMRPredict_panel_results_wider.xlsx"))
    ap.add_argument("--pa-workbook", default=str(HOME / "AMRPredict_pseaer_reported_results.xlsx"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="figure3_combined_draft")
    ap.add_argument("--width", type=float, default=7.8)
    ap.add_argument("--scale", choices=["linear", "sqrt", "sqrt5", "lin20", "power", "split", "steps"], default="sqrt",
                    help="linear: equal distances are equal differences; sqrt: spreads the low-error corner; "
                         "power: x**0.65; split: linear to 10%%, compressed (x0.3) above")
    ap.add_argument("--inset", type=float, default=0.0,
                    help="linear scale only: add an inset magnifying 0..INSET %% on both axes (e.g. 20)")
    ap.add_argument("--upper-step", type=float, default=12.0,
                    help="lin20: length of the 20-40 and 40-80 steps, in units of 1%% below 20 (20 = linear to 40)")
    ap.add_argument("--summary", action="store_true",
                    help="print the closer/farther counts above panel a (off: they go in the text)")
    ap.add_argument("--colour", choices=["distance", "organism"], default="organism",
                    help="arrows coloured by closer/farther to the targets, or arrows and points by organism")
    a = ap.parse_args()
    books = {"mtb": load(pathlib.Path(a.mtb_workbook))[0], "pseaer": load(pathlib.Path(a.pa_workbook))[0]}

    fig = plt.figure(figsize=(a.width, 7.0))
    # F1 first (a, b: where the framework beats the tools cleanly), the error plane after it (c)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.05], hspace=0.62, wspace=0.42)
    ax = fig.add_subplot(gs[:, 1])
    tidy(ax, x=True, y=True)
    BREAK, SLOPE = 10.0, 0.3        # split scale: linear up to BREAK, SLOPE as much length per point above it
    # steps: every tick interval the same length -- 0, 1, 2, 3, 4, 5, then 10, then every 10
    STEP_TICKS = [0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
    # lin20: each percentage point the same length up to 20%; above, 20->40 and 40->80 are equal, shorter steps
    U = a.upper_step
    L20 = ([0.0, 20.0, 40.0, 80.0, 160.0], [0.0, 20.0, 20.0 + U, 20.0 + 2 * U, 20.0 + 3 * U])
    if a.scale == "lin20":
        tf = lambda v: np.interp(v, *L20)  # noqa: E731
        inv = lambda u: float(np.interp(u, L20[1], L20[0]))  # noqa: E731
    elif a.scale == "sqrt5":
        # square root above 5%; below it linear, so 0-1, 1-2, ... 4-5 are equal, and 5% sits where it does on
        # the square-root scale (the two parts meet there)
        K = np.sqrt(5.0) / 5.0
        tf = lambda v: np.where(np.asarray(v) <= 5, np.asarray(v) * K, np.sqrt(np.maximum(v, 0)))  # noqa: E731
        inv = lambda u: u / K if u <= np.sqrt(5.0) else u * u  # noqa: E731
    elif a.scale == "steps":
        tf = lambda v: np.interp(v, STEP_TICKS, range(len(STEP_TICKS)))  # noqa: E731
        inv = lambda u: float(np.interp(u, range(len(STEP_TICKS)), STEP_TICKS))  # noqa: E731
    elif a.scale == "split":
        tf = lambda v: np.where(np.asarray(v) <= BREAK, v, BREAK + (np.asarray(v) - BREAK) * SLOPE)  # noqa: E731
        inv = lambda u: u if u <= BREAK else BREAK + (u - BREAK) / SLOPE  # noqa: E731
    else:
        POW = {"sqrt": 0.5, "power": 0.65, "linear": 1.0}[a.scale]
        tf = lambda v: np.power(v, POW)  # noqa: E731
        inv = lambda u: u ** (1 / POW)  # noqa: E731
    # the very major axis skips the empty stretch between Y_LO and Y_HI (only cefepime's tool point is above it);
    # the skip is drawn as a narrow band, Y_GAP units tall, marked with three dots on the axis
    Y_LO, Y_HI, Y_GAP = 58.0, 76.0, 0.30
    skip = False            # tried 2026-10-20 (skipped band 58-76%); not used

    def ty(v):
        u = tf(v)
        if not skip:
            return u
        lo, hi = float(tf(Y_LO)), float(tf(Y_HI))
        return np.where(np.asarray(v) >= Y_HI, u - (hi - lo) + Y_GAP, np.minimum(u, lo + 0 * u))

    def ty_inv(u):
        lo = float(tf(Y_LO))
        return inv(u) if (not skip or u <= lo) else inv(u + float(tf(Y_HI)) - lo - Y_GAP)

    if a.inset and a.scale != "linear":
        raise SystemExit("[ERROR] --inset goes with --scale linear")
    x1, y1 = 40.0, 86.0      # x ends at 40: the panel edge is the 40% line      # major error rate on x, very major on y: the long range on the tall side
    ax.add_patch(Rectangle((0, 0), tf(T_ME), ty(T_VME), facecolor="#e3f1e7", edgecolor="#8fc19d", lw=0.7, zorder=0))
    if not a.inset:          # with the inset the box is labelled there, where it is large enough
        ax.annotate("target\nregion", (tf(T_ME) / 2, ty(T_VME) / 2), ha="center", va="center", fontsize=5.4,
                    color="#4f8a5f", linespacing=1.1, zorder=1)
    summary, recs = {}, []
    for k in ("mtb", "pseaer"):
        P = PANELS[k]
        drugs = [d for d in P["drugs"] if d not in P["underpowered"]
                 and pick(books[k], d, P["baseline"], "holdout", "SR_f1") is not None]
        get = lambda d, m, f: 100 * pick(books[k], d, m, "holdout", f)  # noqa: E731
        moves = []
        for d in drugs:
            tv, tm = get(d, P["baseline"], "vMER"), get(d, P["baseline"], "MER")
            mv, mm = get(d, "evo2", "vMER"), get(d, "evo2", "MER")
            dd = to_box(mv, mm) - to_box(tv, tm)
            moves.append(dd)
            if a.colour == "organism":
                col = tcol = mcol = ORGCOL[k]
            else:
                col, tcol, mcol = CLOSER if dd < -0.05 else (FARTHER if dd > 0.05 else MUTED), TOOL, EMB
            recs.append((k, d, tv, tm, mv, mm, col, tcol, mcol))
        nc, nf = sum(x < -0.05 for x in moves), sum(x > 0.05 for x in moves)
        summary[k] = (nc, nf, len(moves), wilcoxon(moves).pvalue)
        print(f"[INFO] {ORG[k]} distance to targets: closer {nc}, farther {nf} of {len(moves)}, "
              f"P = {summary[k][3]:.4f}")
    def draw(axx, keep, lab):
        """Markers now; arrows, labels and their obstacles returned, to be placed once the layout is final."""
        arrows_, labels_, obst_ = [], [], []
        for k, d, tv, tm, mv, mm, col, tcol, mcol in recs:
            if not keep(tv, tm, mv, mm):
                continue
            axx.plot(tf(tm), ty(tv), SHAPE[k], ms=5.0 if k == "mtb" else 4.4, mfc=SURFACE, mec=tcol, mew=1.0,
                     zorder=3)
            if abs(mv - tv) + abs(mm - tm) > 0.05 and (k, d) not in NOARROW:
                arrows_.append((tf(tm), ty(tv), tf(mm), ty(mv), col))
            axx.plot(tf(mm), ty(mv), SHAPE[k], ms=5.4 if k == "mtb" else 4.8, color=mcol, mew=0, zorder=4)
            dx, dy = tf(mm) - tf(tm), ty(mv) - ty(tv)
            n = (dx * dx + dy * dy) ** 0.5
            obst_.append((tf(mm), ty(mv)))
            for t_ in np.linspace(0.12, 0.88, 9):
                obst_.append((tf(tm) + t_ * dx, ty(tv) + t_ * dy))
            if lab(tv, tm):
                ux, uy = PREFER.get((k, d), (-dx / n, -dy / n) if n > 1e-9 else (0.7, 0.7))
                labels_.append((tf(tm), ty(tv), CODE[d], ux, uy))
        return arrows_, labels_, obst_

    Z = a.inset
    inside = lambda v, m: v <= Z and m <= Z  # noqa: E731
    arrows, labels, obstacles = draw(ax, lambda *_: True, (lambda v, m: not inside(v, m)) if Z else (lambda *_: True))
    if a.scale == "lin20":
        tk = [0, 5, 10, 15, 20, 40, 80]
        ax.xaxis.set_major_locator(FixedLocator([float(tf(t)) for t in tk if t <= x1]))
        ax.yaxis.set_major_locator(FixedLocator([float(ty(t)) for t in tk if t <= y1]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{inv(v):.0f}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{ty_inv(v):.0f}"))
        ax.axvline(20, color="#b9b7ae", lw=0.8, ls=(0, (1.5, 1.5)), zorder=1)
        ax.axhline(20, color="#b9b7ae", lw=0.8, ls=(0, (1.5, 1.5)), zorder=1)
    elif a.scale == "sqrt5":
        tk = [0, 1, 2, 3, 5, 10, 20, 40, 80]
        ax.xaxis.set_major_locator(FixedLocator([float(tf(t)) for t in tk if t <= x1]))
        ax.yaxis.set_major_locator(FixedLocator([float(ty(t)) for t in tk if t <= y1]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{inv(v):.0f}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{ty_inv(v):.0f}"))
        if skip:                                       # the skipped stretch: a blank band, three dots on the axis
            g0 = float(tf(Y_LO))
            ax.axhspan(g0, g0 + Y_GAP, color=SURFACE, zorder=1.5, lw=0)
            for j in (-1, 0, 1):                       # three dots on the axis, centred on the band
                ax.plot([0.0], [g0 + Y_GAP / 2 + j * Y_GAP * 0.28], "o", ms=1.6, color=INK2, clip_on=False,
                        zorder=6)
    elif a.scale == "steps":
        ax.xaxis.set_major_locator(FixedLocator([float(tf(t)) for t in STEP_TICKS if t <= x1]))
        ax.yaxis.set_major_locator(FixedLocator([float(tf(t)) for t in STEP_TICKS if t <= y1]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{inv(v):.0f}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{inv(v):.0f}"))
        for b_ in (5, 10):                             # where the step size changes
            ax.axvline(float(tf(b_)), color="#b9b7ae", lw=0.8, ls=(0, (1.5, 1.5)), zorder=1)
            ax.axhline(float(tf(b_)), color="#b9b7ae", lw=0.8, ls=(0, (1.5, 1.5)), zorder=1)
    elif a.scale == "split":
        tk = [0, 2, 4, 6, 8, 10, 20, 30, 40, 60, 80]
        ax.xaxis.set_major_locator(FixedLocator([float(tf(t)) for t in tk if t <= x1]))
        ax.yaxis.set_major_locator(FixedLocator([float(tf(t)) for t in tk if t <= y1]))
        sqfmt = FuncFormatter(lambda v, _: f"{inv(v):.0f}")
        ax.xaxis.set_major_formatter(sqfmt)
        ax.yaxis.set_major_formatter(sqfmt)
        for f_ in (ax.axvline, ax.axhline):            # where the scale changes
            f_(BREAK, color="#b9b7ae", lw=0.8, ls=(0, (1.5, 1.5)), zorder=1)
    elif a.scale != "linear":
        ax.xaxis.set_major_locator(FixedLocator([tf(t) for t in TICKS if t <= x1]))   # the same ticks on both axes
        ax.yaxis.set_major_locator(FixedLocator([tf(t) for t in TICKS if t <= y1]))
        sqfmt = FuncFormatter(lambda v, _: f"{inv(v):.0f}")
        ax.xaxis.set_major_formatter(sqfmt)
        ax.yaxis.set_major_formatter(sqfmt)
    else:
        ax.xaxis.set_major_locator(MultipleLocator(10))                                # same spacing on both axes
        ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.set_xlim(0, tf(x1))
    ax.set_ylim(0, float(ty(y1)))
    ax.set_aspect("equal", adjustable="box", anchor="NE")   # one unit the same length on both axes
    ins = None
    if Z:
        ins = ax.inset_axes([0.56, 0.47, 0.42, 0.42])
        tidy(ins, x=True, y=True)
        ins.add_patch(Rectangle((0, 0), T_ME, T_VME, facecolor="#e3f1e7", edgecolor="#8fc19d", lw=0.7, zorder=0))
        ins.annotate("target", (T_ME / 2, T_VME / 2), ha="center", va="center", fontsize=4.8, color="#4f8a5f",
                     zorder=1)
        i_arrows, i_labels, i_obst = draw(ins, lambda tv, tm, mv, mm: inside(tv, tm) or inside(mv, mm),
                                         lambda v, m: inside(v, m))
        ins.set_xlim(0, Z)
        ins.set_ylim(0, Z)
        ins.xaxis.set_major_locator(MultipleLocator(5))
        ins.yaxis.set_major_locator(MultipleLocator(5))
        ins.tick_params(labelsize=5.6)
        ins.set_facecolor(SURFACE)
        for sp in ins.spines.values():
            sp.set_visible(True)
            sp.set_color("#8a8880")
            sp.set_linewidth(0.7)
        ax.indicate_inset_zoom(ins, edgecolor="#8a8880", lw=0.7, alpha=1.0)
        ins.set_title(f"0–{Z:g}%, enlarged", loc="left", fontsize=5.8, color=INK2, pad=3)
    ax.tick_params(labelsize=6.6)
    ax.set_xlabel("Major error rate (%)", fontsize=6.8)
    ax.set_ylabel("Very major error rate (%)", fontsize=6.8)
    ax.set_title("Error rates", loc="left", color=INK, pad=26 if a.summary else 13, fontweight="semibold",
                 fontsize=7.6)
    if not a.summary:
        ax.annotate("framework vs tool, both organisms", xy=(0.0, 1.0),
                    xycoords="axes fraction", xytext=(0, 3), textcoords="offset points", ha="left", va="bottom",
                    fontsize=5.4, color=INK2)
    lines = [(ORG[k], f": closer to the targets {nc}, farther {nf}"
              + (f", unchanged {m - nc - nf}" if m - nc - nf else "") + f" of {m} · P = {pval(p)}")
             for k, (nc, nf, m, p) in summary.items()]
    letters = [("c", ax)]
    sub = []

    # b, c: F1 gain per organism
    for r, k in enumerate(("mtb", "pseaer")):
        P = PANELS[k]
        drugs = [d for d in P["drugs"] if d not in P["underpowered"]
                 and pick(books[k], d, P["baseline"], "holdout", "SR_f1") is not None]
        get = lambda d, m, f: 100 * pick(books[k], d, m, "holdout", f)  # noqa: E731
        bx = fig.add_subplot(gs[r, 0])
        x_ = [get(d, "evo2", "SR_f1") for d in drugs]
        y_ = [get(d, P["baseline"], "SR_f1") for d in drugs]
        col_ = ORGCOL[k] if a.colour == "organism" else EMB
        mshape = SHAPE[k]
        dumbbell(bx, drugs, x_, y_,
                 a_style=dict(marker=mshape, ms=5.2 if k == "mtb" else 4.6, color=col_, mew=0, ls="none"),
                 b_style=dict(marker=mshape, ms=4.8 if k == "mtb" else 4.2, mfc=SURFACE, mec=col_, mew=1.0, ls="none"),
                 title="", xlabel="", better_low=False, true_delta=True, unit=100.0)
        bx.tick_params(axis="x", labelsize=6.4)
        bx.set_xlabel("F1 (%)", fontsize=6.4)
        wins, ties, m, med, p = paired_summary(drugs, x_, y_, False, true_delta=True, unit=100.0)
        bx.set_title("F1", loc="left", color=INK, pad=13, fontweight="semibold", fontsize=7.6)
        sub.append((bx, ORG[k], f" · model better {wins}" + (f", equal {ties}" if ties else "")
                    + f" of {m} · median {med:+.1f}% · P = {pval(p)}"))
        letters.append(("ab"[r], bx))

    top = header(fig, "Figure 3, combined plane (draft) — the resistance call against the established tools",
                 "Held-out sets. The reported model, which contains the tool's own call as a feature, against the tool "
                 "alone (circles: TB-Profiler; diamonds: ARDaP,\nintermediate counted as resistant). a, b: F1 for the "
                 "resistant class, tool (open) and model (filled); numbers at the right, model minus tool\nin percentage "
                 "points. c: one arrow per antimicrobial from the tool (open) to the model (filled), labelled at the "
                 "tool's point; axes\n"
                 + {"lin20": "linear up to 20% (each point the same length);\nabove it, 20-40% and 40-80% are "
                             "equal steps (dotted lines); the same on both axes",
                    "sqrt5": "linear up to 5% and on a square-root scale above it, the same on both axes",
                    "steps": "with equal spacing per 1% up to 5%, then 5-10%, then per 10% (dotted lines)",
                    "split": "linear up to 10%; above it each point takes 0.3 of the length (dotted lines)",
                    "power": "compressed (power 0.65) to give the low error rates room", "sqrt": "on a square-root scale",
                    "linear": "linear" + (f", inset: 0–{a.inset:g}% enlarged" if a.inset else "")}[a.scale]
                 + "; green, target region (very major ≤1.5%, major ≤3%; strict point targets derived from the FDA criteria). Colours: organism (key). "
                 "Bedaquiline and\ndelamanid are left out.", pad_in=0.30)
    fig.subplots_adjust(top=top - 0.025, bottom=0.12, left=0.09 + (0.04 if FINAL or SUBMISSION else 0.0), right=0.975 if FINAL or SUBMISSION else 0.93)
    fig.canvas.draw()
    # b and c share the height of a: equal panels, b's top on a's top, c's bottom on a's bottom
    pa_, bpos = ax.get_position(), letters[1][1].get_position()
    gap = 0.115
    h_ = (pa_.y1 - pa_.y0 - gap) / 2
    letters[1][1].set_position([bpos.x0, pa_.y1 - h_, bpos.width, h_])
    letters[2][1].set_position([bpos.x0, pa_.y0, bpos.width, h_])
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for bx, org, rest in sub:                        # "M. tuberculosis · model better ..." with the name in italics
        t = bx.annotate(org, xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points",
                        ha="left", va="bottom", fontproperties=ITALIC, fontsize=5.4, color=INK2)
        w = t.get_window_extent(rend).width * 72.0 / fig.dpi
        bx.annotate(rest, xy=(0.0, 1.0), xycoords="axes fraction", xytext=(w, 3), textcoords="offset points",
                    ha="left", va="bottom", fontsize=5.4, color=INK2)
    # panel letters as in Figure 2: left of everything in the panel, level with the panel title
    inv_ = fig.transFigure.inverted()
    colx = {}                                        # one x per column: the leftmost element of that column
    for lt, axx in letters:
        c_ = axx.get_subplotspec().colspan.start
        colx[c_] = min(colx.get(c_, 1.0), inv_.transform((axx.get_tightbbox(rend).x0, 0))[0])
    for lt, axx in letters:
        ti = axx.title.get_window_extent(rend)     # baseline of the letter on the baseline of the title
        fig.text(colx[axx.get_subplotspec().colspan.start] - 0.004, inv_.transform((0, ti.y0))[1], lt,
                 fontsize=10, fontweight="bold", color=INK, ha="right", va="bottom")
    for i, (org, rest) in enumerate(lines if a.summary else []):   # species italic, the rest upright
        y = 3 + (len(lines) - 1 - i) * 9.0
        t = ax.annotate(org, xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, y), textcoords="offset points",
                        ha="left", va="bottom", fontproperties=ITALIC, fontsize=5.6, color=INK2)
        w = t.get_window_extent(rend).width * 72.0 / fig.dpi
        ax.annotate(rest, xy=(0.0, 1.0), xycoords="axes fraction", xytext=(w, y), textcoords="offset points",
                    ha="left", va="bottom", fontsize=5.6, color=INK2)
    for arr in arrows:
        dashed_arrow(ax, *arr)
    spread_labels(fig, ax, labels, fs=6.0, msz=8.0, extra=obstacles, esz=4.5)
    if ins is not None:
        for arr in i_arrows:
            dashed_arrow(ins, *arr)
        spread_labels(fig, ins, i_labels, fs=5.6, msz=8.0, extra=i_obst, esz=4.5)
    box = Rectangle((0, 0), 1, 1, facecolor="#e3f1e7", edgecolor="#8fc19d", lw=0.7,
                    label="Targets (very major ≤1.5%, major ≤3%)")
    if a.colour == "organism":
        handles = [Line2D([], [], marker="o", ls="none", ms=5.0, mfc=SURFACE, mec=INK2, mew=1.0, label="tool alone"),
                   Line2D([], [], marker="o", ls="none", ms=5.4, color=INK2, mec=SURFACE, mew=0.7,
                          label="reported model"),
                   Line2D([], [], marker="o", ls=(0, (3.2, 2.0)), ms=5.0, color=ORGCOL["mtb"], mec=SURFACE, lw=1.1,
                          label="M. tuberculosis"),
                   Line2D([], [], marker="D", ls=(0, (3.2, 2.0)), ms=4.4, color=ORGCOL["pseaer"], mec=SURFACE,
                          lw=1.1, label="P. aeruginosa"), box]
    else:
        handles = [Line2D([], [], marker="o", ls="none", ms=5.0, mfc=SURFACE, mec=TOOL, mew=1.0, label="tool alone"),
                   Line2D([], [], marker="o", ls="none", ms=5.4, color=EMB, mec=SURFACE, mew=0.7,
                          label="reported model"),
                   Line2D([], [], marker="o", ls="none", ms=5.0, mfc="none", mec=INK2, mew=1.0, label="M. tuberculosis"),
                   Line2D([], [], marker="D", ls="none", ms=4.4, mfc="none", mec=INK2, mew=1.0, label="P. aeruginosa"),
                   Line2D([], [], color=CLOSER, lw=1.6, ls=(0, (3.2, 2.0)), label="closer to the targets"),
                   Line2D([], [], color=FARTHER, lw=1.6, ls=(0, (3.2, 2.0)), label="farther"),
                   Line2D([], [], color=MUTED, lw=1.6, ls=(0, (3.2, 2.0)), label="same distance"), box]
    fig.legend(handles=handles,
               loc="upper center", ncol=len(handles), fontsize=6.4, handletextpad=0.5, columnspacing=1.4,
               bbox_to_anchor=(0.5, fig.transFigure.inverted().transform((0, ax.get_tightbbox(rend).y0))[1] - 0.012))
    for t in fig.legends[0].get_texts():          # species names in italics
        if t.get_text() in ORG.values():
            t.set_fontproperties(ITALIC)
            t.set_fontsize(6.4)
    stem = a.stem if a.colour == "organism" else a.stem.replace("_draft", "_distance_draft")
    stem = stem.replace("_draft", "_inset_draft") if Z else stem
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    scale_letters(fig)
    fit_fonts(fig)
    if FINAL or SUBMISSION:                         # full width (uniform with the other figures), trimmed below only
        from matplotlib.transforms import Bbox
        tb = fig.get_tightbbox(fig.canvas.get_renderer())
        crop = Bbox([[0.0, max(0.0, tb.y0 - 0.06)], [fig.get_figwidth(), fig.get_figheight()]])
    for ext in ("pdf", "png"):                      # trimmed: the legend sits right under the panels
        if FINAL or SUBMISSION:
            fig.savefig(out / f"{stem}.{ext}", format=ext, bbox_inches=crop)
        else:
            fig.savefig(out / f"{stem}.{ext}", format=ext, bbox_inches="tight", pad_inches=0.15)
    print(f"[OK] {out / stem}.pdf / .png")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
