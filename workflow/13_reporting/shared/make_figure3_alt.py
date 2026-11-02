#!/usr/bin/env python3
"""Figure 3, alternative layout: the resistance call as an error trade-off and an F1 gain.

The same held-out numbers as make_figure3.py (the reported model, which contains the tool's own
call, against the tool alone), drawn differently so that the figure does not repeat the dumbbell
grid of Figure 2.

  top row     one plane per organism, very major error rate (resistant called susceptible) on x,
              major error rate (susceptible called resistant) on y; one arrow per antimicrobial from
              the tool (open circle) to the reported model (filled). Lower left is better on both.
  bottom row  F1 for the resistant class, model minus tool, one bar per antimicrobial, sorted.

M. tuberculosis: eight antimicrobials (bedaquiline and delamanid, with 12 and 15 resistant
held-out isolates, are left out); P. aeruginosa: the eight that ARDaP calls.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator
import numpy as np
from scipy.stats import wilcoxon

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import INK, INK2, ITALIC_BOLD, MUTED, SERIES, SURFACE, header, load, pick, save, tidy  # noqa: E402
from make_provisional_figures import paired_summary  # noqa: E402
from plot_metric_bars import PANELS  # noqa: E402

HOME = pathlib.Path.home()
EMB = SERIES["evo2"]
ORG = {"mtb": "M. tuberculosis", "pseaer": "P. aeruginosa"}
CODE = {"isoniazid": "INH", "rifampicin": "RIF", "ethambutol": "EMB", "ethionamide": "ETH", "kanamycin": "KAN",
        "amikacin": "AMK", "levofloxacin": "LEV", "moxifloxacin": "MXF", "ciprofloxacin": "CIP",
        "tobramycin": "TOB", "imipenem": "IPM", "meropenem": "MEM", "ceftazidime": "CAZ", "cefepime": "FEP",
        "aztreonam": "ATM", "piperacillin-tazobactam": "TZP", "bedaquiline": "BDQ", "delamanid": "DLM"}
NAME = {"piperacillin-tazobactam": "Pip.–tazobactam"}
LIMS = {"mtb": ((-0.006, 0.19), (-0.016, 0.21)), "pseaer": ((-0.02, 0.86), (-0.020, 0.36))}
LIMS_PCT = {"mtb": ((-0.8, 19.0), (-1.2, 21.0)), "pseaer": ((-2.0, 86.0), (-1.5, 36.0))}   # per cent
STEP = {"mtb": (5, 5), "pseaer": (20, 10)}                                                 # gridline spacing
TICKS = {"mtb": [0, 1, 2, 5, 10, 15, 20], "pseaer": [0, 1, 5, 10, 20, 30, 40, 80]}             # square-root axes
T_VME, T_ME = 1.5, 3.0          # commonly used targets for susceptibility testing (very major, major; %)
CLOSER, FARTHER = "#2a78d6", "#e08a3c"
PREFER = {("mtb", "rifampicin"): (-1.0, 0.0),
          ("pseaer", "imipenem"): (-1.0, 0.3), ("pseaer", "piperacillin-tazobactam"): (-0.7, 1.0),
          ("pseaer", "meropenem"): (1.0, 0.6), ("pseaer", "ceftazidime"): (-0.3, -1.0),
          ("pseaer", "cefepime"): (0.6, 1.0)}      # label direction where the default crowds


def to_box(v, m):
    """Distance from (very major, major) to the target box, in percentage points; 0 inside it."""
    return (max(0.0, v - T_VME) ** 2 + max(0.0, m - T_ME) ** 2) ** 0.5


def spread_labels(fig, ax, items, fs=5.8, msz=6.5, extra=(), esz=4.0):
    """Place the antimicrobial codes beside their model points without overlapping each other.

    items: (x, y, text, ux, uy) in data units with the preferred direction (ux, uy). Boxes are
    approximated from the font size; overlapping pairs are pushed apart, then each code is drawn at
    its final offset in points, with a thin leader line when it ended up away from its point.
    """
    k = 72.0 / fig.dpi
    w = [len(t) * 0.60 * fs + 1.5 for _, _, t, _, _ in items]
    h = 1.25 * fs
    base = [tuple(c * k for c in ax.transData.transform((x, y))) for x, y, *_ in items]
    # further obstacles (data coordinates): e.g. the tool's points and points sampled along the arrows
    ext = [tuple(c * k for c in ax.transData.transform((x, y))) for x, y in extra]
    pos = [[bx + 8.5 * ux, by + 8.5 * uy] for (bx, by), (_, _, _, ux, uy) in zip(base, items)]
    for _ in range(300):
        moved = False
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                dx, dy = pos[j][0] - pos[i][0], pos[j][1] - pos[i][1]
                ox, oy = (w[i] + w[j]) / 2 - abs(dx), h - abs(dy)
                if ox > 0 and oy > 0:
                    moved = True
                    if oy < ox:
                        sgn = 1 if dy >= 0 else -1
                        pos[i][1] -= sgn * oy / 2 * 0.6; pos[j][1] += sgn * oy / 2 * 0.6
                    else:
                        sgn = 1 if dx >= 0 else -1
                        pos[i][0] -= sgn * ox / 2 * 0.6; pos[j][0] += sgn * ox / 2 * 0.6
        # the markers are obstacles too: a code must not sit on another antimicrobial's point
        for i in range(len(pos)):
            for j in range(len(base)):
                dx, dy = pos[i][0] - base[j][0], pos[i][1] - base[j][1]
                ox, oy = (w[i] + msz) / 2 - abs(dx), (h + msz) / 2 - abs(dy)
                if ox > 0 and oy > 0:
                    moved = True
                    if oy < ox:
                        pos[i][1] += (1 if dy >= 0 else -1) * oy * 0.6
                    else:
                        pos[i][0] += (1 if dx >= 0 else -1) * ox * 0.6
            for ex, ey in ext:
                dx, dy = pos[i][0] - ex, pos[i][1] - ey
                ox, oy = (w[i] + esz) / 2 - abs(dx), (h + esz) / 2 - abs(dy)
                if ox > 0 and oy > 0:
                    moved = True
                    if oy < ox:
                        pos[i][1] += (1 if dy >= 0 else -1) * oy * 0.6
                    else:
                        pos[i][0] += (1 if dx >= 0 else -1) * ox * 0.6
        if not moved:
            break
    for (x, y, t, _, _), (bx, by), (px, py) in zip(items, base, pos):
        off = (px - bx, py - by)
        far = (off[0] ** 2 + off[1] ** 2) ** 0.5 > 11
        ax.annotate(t, (x, y), xytext=off, textcoords="offset points", ha="center", va="center",
                    fontsize=fs, color=INK, zorder=5,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#b5b3a9", shrinkA=0, shrinkB=2.5)
                    if far else None)


def pval(p: float) -> str:
    return f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-workbook", default=str(HOME / "Desktop/AMRPredict_panel_results_wider.xlsx"))
    ap.add_argument("--pa-workbook", default=str(HOME / "Desktop/AMRPredict_pseaer_reported_results.xlsx"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="figure3_alt_draft")
    ap.add_argument("--width", type=float, default=7.8)
    a = ap.parse_args()
    books = {"mtb": load(pathlib.Path(a.mtb_workbook))[0], "pseaer": load(pathlib.Path(a.pa_workbook))[0]}

    # one row per organism: the error plane, then the F1 gain
    fig = plt.figure(figsize=(a.width, 7.4))
    gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.42, width_ratios=[1.05, 1.0])
    planes, labels, rowax, obstacles = {}, {}, {}, {}
    ARROW, TOOL = "#b9b7ae", "#8a8880"
    for r, k in enumerate(("mtb", "pseaer")):
        P = PANELS[k]
        drugs = [d for d in P["drugs"] if d not in P["underpowered"]
                 and pick(books[k], d, P["baseline"], "holdout", "SR_f1") is not None]
        get = lambda d, m, f: 100 * pick(books[k], d, m, "holdout", f)  # noqa: E731  per cent

        # ---- error plane on square-root axes, with the target box; arrows coloured by whether the
        # model moves closer to the box (distance = excess over the targets, in percentage points)
        ax = fig.add_subplot(gs[r, 0])
        tidy(ax, x=True, y=True)
        tf = np.sqrt
        (_, x1), (_, y1) = LIMS_PCT[k]
        ax.add_patch(Rectangle((0, 0), tf(T_VME), tf(T_ME), facecolor="#e3f1e7", edgecolor="#8fc19d", lw=0.7,
                               zorder=0))
        ax.annotate("targets", (tf(T_VME) / 2, 0), xytext=(0, 3), textcoords="offset points", ha="center",
                    va="bottom", fontsize=5.6, color="#4f8a5f", zorder=1)
        closer = []
        for d in drugs:
            tv, tm = get(d, P["baseline"], "vMER"), get(d, P["baseline"], "MER")
            mv, mm = get(d, "evo2", "vMER"), get(d, "evo2", "MER")
            dd = to_box(mv, mm) - to_box(tv, tm)
            closer.append(dd)
            col = CLOSER if dd < -0.05 else (FARTHER if dd > 0.05 else MUTED)
            ax.plot(tf(tv), tf(tm), "o", ms=5.0, mfc=SURFACE, mec=TOOL, mew=1.0, zorder=3)
            if abs(mv - tv) + abs(mm - tm) > 0.05:
                ax.annotate("", xy=(tf(mv), tf(mm)), xytext=(tf(tv), tf(tm)), zorder=2,
                            arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3, shrinkA=3.0, shrinkB=3.2,
                                            mutation_scale=8))
            ax.plot(tf(mv), tf(mm), "o", ms=5.2, color=EMB, mec=SURFACE, mew=0.7, zorder=4)
            dx, dy = tf(mv) - tf(tv), tf(mm) - tf(tm)
            norm = (dx * dx + dy * dy) ** 0.5
            ux, uy = (dx / norm, dy / norm) if norm > 1e-9 else (0.7, 0.7)
            ux, uy = PREFER.get((k, d), (ux, uy))
            labels.setdefault(k, []).append((tf(mv), tf(mm), CODE[d], ux, uy))
            # obstacles for the labels: the tool's point and points along the arrow
            obstacles.setdefault(k, []).append((tf(tv), tf(tm)))
            for t_ in np.linspace(0.12, 0.88, 9):
                obstacles[k].append((tf(tv) + t_ * dx, tf(tm) + t_ * dy))
        ticks = [t for t in TICKS[k] if t <= x1]
        ax.xaxis.set_major_locator(FixedLocator([tf(t) for t in ticks]))
        ax.yaxis.set_major_locator(FixedLocator([tf(t) for t in TICKS[k] if t <= y1]))
        sqfmt = FuncFormatter(lambda v, _: f"{v * v:.0f}")
        ax.xaxis.set_major_formatter(sqfmt)
        ax.yaxis.set_major_formatter(sqfmt)
        ax.set_xlim(0, tf(x1))
        ax.set_ylim(0, tf(y1))
        ax.tick_params(labelsize=6.6)
        ax.set_xlabel("Very major error rate (%)", fontsize=6.6)
        ax.set_ylabel("Major error rate (%)", fontsize=6.6)
        ax.set_title("Error trade-off", loc="left", color=INK, pad=13, fontweight="semibold", fontsize=8.0)
        nc, nf = sum(x < -0.05 for x in closer), sum(x > 0.05 for x in closer)
        pc = wilcoxon(closer).pvalue if any(abs(x) > 0 for x in closer) else float("nan")
        ax.annotate(f"closer to the targets {nc}, farther {nf}" + (f", unchanged {len(closer) - nc - nf}"
                    if len(closer) - nc - nf else "") + f" of {len(closer)} · P = {pval(pc)}",
                    xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points",
                    ha="left", va="bottom", fontsize=5.6, color=INK2)
        print(f"[INFO] {ORG[k]} distance to targets: closer {nc}, farther {nf} of {len(closer)}, P = {pc:.4f}")
        ax.annotate("ac"[r], xy=(-0.30, 1.06), xycoords="axes fraction", fontsize=10, fontweight="bold",
                    color=INK, ha="right", va="bottom", annotation_clip=False)
        planes[k] = ax
        rowax[k] = ax

        # ---- F1 gain, percentage points
        bx = fig.add_subplot(gs[r, 1])
        tidy(bx, x=True, y=False)
        gain = sorted(((get(d, "evo2", "SR_f1") - get(d, P["baseline"], "SR_f1"), d) for d in drugs), reverse=True)
        for i, (g, d) in enumerate(gain):
            bx.barh(i, g, 0.60, color=EMB if g > 0.05 else ("#b7c7e3" if g < -0.05 else MUTED),
                    edgecolor=SURFACE, linewidth=0.5, zorder=3)
            bx.annotate(f"{g:+.1f}" if abs(g) > 0.05 else "0.0", (g, i),
                        xytext=(3 if g >= 0 else -3, 0), textcoords="offset points",
                        ha="left" if g >= 0 else "right", va="center", fontsize=6.2, color=INK2)
        bx.axvline(0, color=INK2, lw=0.8, zorder=4)
        bx.set_yticks(range(len(gain)))
        bx.set_yticklabels([NAME.get(d, d.capitalize()) for _, d in gain], fontsize=7.0, color=INK)
        bx.invert_yaxis()
        lo, hi = min(g for g, _ in gain), max(g for g, _ in gain)
        pad = 0.20 * (hi - min(lo, 0))
        bx.set_xlim(min(lo, 0) - pad, hi + pad)
        bx.tick_params(axis="x", labelsize=6.6)
        bx.set_xlabel("F1 of the model minus F1 of the tool\n(percentage points)", fontsize=6.6, linespacing=1.3)
        wins, ties, m, med, p = paired_summary(drugs, [get(d, "evo2", "SR_f1") for d in drugs],
                                               [get(d, P["baseline"], "SR_f1") for d in drugs], False, unit=100.0)
        bx.set_title("F1 gain, resistant class", loc="left", color=INK, pad=13, fontweight="semibold", fontsize=8.0)
        bx.annotate(f"model better {wins}, worse {m - wins - ties}" + (f", equal {ties}" if ties else "")
                    + f" of {m} · median {med:+.1f} points · P = {pval(p)}",
                    xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points",
                    ha="left", va="bottom", fontsize=5.6, color=INK2)
        bx.annotate("bd"[r], xy=(-0.42, 1.06), xycoords="axes fraction", fontsize=10, fontweight="bold",
                    color=INK, ha="right", va="bottom", annotation_clip=False)
        print(f"[INFO] {ORG[k]} F1 gain: better {wins}, equal {ties} of {m}, median {med:+.2f}, P = {p:.4f}")
    top = header(fig, "Figure 3, alternative (draft) — the resistance call against the established tools",
                 "Held-out sets. The reported model, which contains the tool's own call as a feature, against the tool "
                 "alone: TB-Profiler in the first row,\nARDaP in the second (intermediate counted as resistant; the eight "
                 "antimicrobials it calls). a, c: one arrow per antimicrobial from the tool (open)\nto the model (filled), "
                 "square-root axes; very major error, resistant isolate called susceptible; major error, susceptible\n"
                 "isolate called resistant. Green: commonly used targets for susceptibility testing; arrows blue when the "
                 "model is\ncloser to them (distance = excess over the targets), orange when farther, grey at the same "
                 "distance. b, d: F1 for the\nresistant class, model minus tool. Bedaquiline and delamanid (12 and 15 "
                 "resistant held-out isolates) are left out.", pad_in=0.30)
    fig.subplots_adjust(top=top - 0.005, bottom=0.135, left=0.155, right=0.965)
    for k, ax in planes.items():
        spread_labels(fig, ax, labels[k], fs=6.0, msz=8.0, extra=obstacles[k], esz=4.5)
    for k, ax in rowax.items():
        pos = ax.get_position()
        fig.text(0.022, (pos.y0 + pos.y1) / 2, ORG[k], rotation=90, ha="center", va="center", fontsize=9,
                 fontproperties=ITALIC_BOLD, color=INK)
    fig.legend(handles=[Line2D([], [], marker="o", ls="none", ms=5.0, mfc=SURFACE, mec=TOOL, mew=1.0,
                               label="Tool alone (TB-Profiler, ARDaP)"),
                        Line2D([], [], marker="o", ls="none", ms=5.2, color=EMB, mec=SURFACE, mew=0.7,
                               label="Reported model (embeddings + the tool's call)"),
                        Line2D([], [], color=CLOSER, lw=1.6, label="closer to the targets"),
                        Line2D([], [], color=FARTHER, lw=1.6, label="farther from the targets"),
                        Line2D([], [], color=MUTED, lw=1.6, label="same distance"),
                        Rectangle((0, 0), 1, 1, facecolor="#e3f1e7", edgecolor="#8fc19d", lw=0.7,
                                  label="targets (very major ≤1.5%, major ≤3%)")],
               loc="lower center", ncol=3, fontsize=6.6, handletextpad=0.5, columnspacing=1.8,
               bbox_to_anchor=(0.5, 0.0))
    save(fig, pathlib.Path(a.out_dir), a.stem)
    return 0


if __name__ == "__main__":
    sys.exit(main())
