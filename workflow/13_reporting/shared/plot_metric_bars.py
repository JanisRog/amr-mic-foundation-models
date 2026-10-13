#!/usr/bin/env python3
"""F1 and one-tier accuracy as grouped bars: one group per drug, one bar per model.

Two panels over a shared x axis. The F1 panel also carries the organism's
comparator (TB-Profiler or ARDaP), which makes a resistance call and so has an F1;
it has no MIC prediction and therefore no bar in the one-tier panel. Held-out
sets, read from the results workbook, never recomputed.

The comparator is not a fourth hue. As in the other figures it is the incumbent
everything is read against, so it takes an open bar in neutral ink.

A second figure (bars_median_<panel>) reduces the first to one bar per model: the
median over the antimicrobials, with each antimicrobial's value drawn as a point
on its bar so the spread behind the median stays visible. Antimicrobials with too
few resistant held-out isolates are left out of the medians, as in every panel
comparison; the comparator's median is over the antimicrobials it covers.

Bars start at zero. On a truncated axis the differences between models would look
several times larger than they are, and the point of this figure is how close the
three models sit.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import statistics

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import (GRID, INK, INK2, MODEL_LABEL, MODELS, MUTED, SERIES,  # noqa: E402
                                SURFACE, header, load, pick, save, tidy)

PANELS = {
    "mtb": dict(drugs=["isoniazid", "rifampicin", "ethambutol", "ethionamide", "kanamycin",
                       "amikacin", "levofloxacin", "moxifloxacin", "bedaquiline", "delamanid"],
                baseline="TB-Profiler 6.3.0", baseline_label="TB-Profiler",
                underpowered=["bedaquiline", "delamanid"], organism="M. tuberculosis"),
    "pseaer": dict(drugs=["ciprofloxacin", "levofloxacin", "tobramycin", "amikacin", "imipenem",
                          "meropenem", "ceftazidime", "cefepime", "aztreonam",
                          "piperacillin-tazobactam"],
                   baseline="ARDaP (I=R)", baseline_label="ARDaP",
                   underpowered=[], organism="P. aeruginosa"),
}
SHORT = {"piperacillin-tazobactam": "Piperacillin–\ntazobactam"}
W = 0.19            # bar width; four bars and a gap fill one drug's slot


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workbook", required=True)
    ap.add_argument("--panel", choices=sorted(PANELS), required=True)
    ap.add_argument("--stem", default="")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    P = PANELS[a.panel]
    summary, _ = load(pathlib.Path(a.workbook))
    drugs = P["drugs"]

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.3), sharex=True,
                             gridspec_kw={"hspace": 0.30})
    n_base = 0
    for ax, field, title in ((axes[0], "SR_f1", "F1, resistant class"),
                             (axes[1], "OT", "One-tier accuracy  (predicted MIC within one "
                                             "doubling dilution)")):
        tidy(ax, x=False, y=True)
        with_base = field == "SR_f1"
        for i, d in enumerate(drugs):
            faint = d in P["underpowered"]
            series = [(m, pick(summary, d, m, "holdout", field)) for m in MODELS]
            if with_base:
                series.append(("base", pick(summary, d, P["baseline"], "holdout", field)))
            k = len(series)
            for j, (m, v) in enumerate(series):
                if v is None:
                    continue
                x = i + (j - (k - 1) / 2) * (W + 0.012)
                if m == "base":
                    n_base += 1
                    ax.bar(x, v, W, facecolor=SURFACE, edgecolor=INK2, linewidth=0.9,
                           alpha=0.55 if faint else 1.0, zorder=3)
                else:
                    ax.bar(x, v, W, color=SERIES[m], edgecolor=SURFACE, linewidth=0.5,
                           alpha=0.45 if faint else 1.0, zorder=3)
                if m == "evo2":      # one selective number per group: the reported model
                    ax.annotate(f"{v:.2f}", (x, v), xytext=(0, 2.2), textcoords="offset points",
                                ha="center", va="bottom", fontsize=5.8, color=INK2)
        ax.set_ylim(0, 1.08)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlim(-0.6, len(drugs) - 0.4)
        ax.set_title(title, loc="left", color=INK, pad=6, fontweight="medium")
        ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    axes[1].set_xticks(range(len(drugs)))
    axes[1].set_xticklabels([SHORT.get(d, d.capitalize()) for d in drugs], fontsize=7.4,
                            color=INK, rotation=32, ha="right", rotation_mode="anchor")
    if P["underpowered"]:
        i0 = drugs.index(P["underpowered"][0])
        for ax in axes:
            ax.axvline(i0 - 0.5, color=GRID, lw=0.9, ls=(0, (3, 3)), zorder=1)
        axes[0].annotate("few resistant\nheld-out isolates", (i0 - 0.42, 1.07), fontsize=6.2,
                         color=MUTED, ha="left", va="top", linespacing=1.3)

    handles = [Patch(facecolor=SERIES[m], edgecolor=SURFACE, label=MODEL_LABEL[m]) for m in MODELS]
    handles.append(Patch(facecolor=SURFACE, edgecolor=INK2, linewidth=0.9,
                         label=f"{P['baseline_label']} (F1 only)"))
    fig.legend(handles=handles, ncol=4, loc="lower center", fontsize=7.5,
               handlelength=1.2, handleheight=0.9, columnspacing=1.8,
               bbox_to_anchor=(0.5, 0.0))
    missing = [d for d in drugs if pick(summary, d, P["baseline"], "holdout", "SR_f1") is None]
    note = (f" {P['baseline_label']} makes no call for "
            f"{' and '.join(missing)}." if missing else "")
    top = header(fig, f"Resistance calls and MIC placement by model — {P['organism']}",
                 f"Held-out sets. Numbers above the bars are Evo 2, the reported model. "
                 f"{P['baseline_label']} returns a resistance call and no MIC,\nso it appears "
                 f"in the upper panel only.{note}", pad_in=0.26)
    fig.subplots_adjust(top=top, bottom=0.185, left=0.065, right=0.99)
    save(fig, pathlib.Path(a.out_dir), a.stem or f"bars_f1_onetier_{a.panel}")
    print(f"[INFO] {a.panel}: {len(drugs)} drugs, comparator bars {n_base}")
    median_figure(summary, P, a)
    return 0


def median_figure(summary, P, a):
    """One bar per model: the median over antimicrobials, each antimicrobial a point."""
    drugs = [d for d in P["drugs"] if d not in P["underpowered"]]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.9), gridspec_kw={"wspace": 0.22})
    covered = None
    for ax, field, title in ((axes[0], "SR_f1", "F1, resistant class"),
                             (axes[1], "OT", "One-tier accuracy")):
        tidy(ax, x=False, y=True)
        series = [(MODEL_LABEL[m].replace(" 7B", "").replace(" 8B", ""), SERIES[m],
                   [pick(summary, d, m, "holdout", field) for d in drugs]) for m in MODELS]
        if field == "SR_f1":
            base = [pick(summary, d, P["baseline"], "holdout", field) for d in drugs]
            covered = sum(v is not None for v in base)
            series.append((P["baseline_label"], None, base))
        for j, (label, color, vals) in enumerate(series):
            vals = [v for v in vals if v is not None]
            med = statistics.median(vals)
            if color is None:
                ax.bar(j, med, 0.56, facecolor=SURFACE, edgecolor=INK2, linewidth=0.9, zorder=3)
            else:
                ax.bar(j, med, 0.56, color=color, edgecolor=SURFACE, linewidth=0.5, zorder=3)
            # every antimicrobial, spread sideways in a fixed order so that a point
            # sits at the same offset on every bar
            n = len(vals)
            xs = [j + (k - (n - 1) / 2) * (0.40 / max(n - 1, 1)) for k in range(n)]
            ax.plot(xs, vals, "o", ms=2.6, mfc=SURFACE, mec=INK, mew=0.6, zorder=5)
            ax.annotate(f"{med:.2f}", (j, 1.035), ha="center", va="bottom", fontsize=7.2,
                        color=INK, annotation_clip=False)
        ax.set_xticks(range(len(series)))
        ax.set_xticklabels([s[0] for s in series], fontsize=7.8, color=INK)
        ax.set_ylim(0, 1.10)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlim(-0.6, len(series) - 0.4)
        ax.set_title(title, loc="left", color=INK, pad=6, fontweight="medium")
        ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    names = " and ".join(P["underpowered"])
    left_out = (f" {names[:1].upper() + names[1:]} are left out "
                f"(few resistant held-out isolates)." if P["underpowered"] else "")
    cover = ("" if covered == len(drugs) else
             f" {P['baseline_label']} covers {covered} of them.")
    top = header(fig, f"Median over antimicrobials by model — {P['organism']}",
                 f"Held-out sets. Bar and number: median over {len(drugs)} antimicrobials; "
                 f"points: each antimicrobial.{cover}\n{P['baseline_label']} returns a "
                 f"resistance call and no MIC, so it has no one-tier accuracy.{left_out}",
                 pad_in=0.26)
    fig.subplots_adjust(top=top, bottom=0.09, left=0.065, right=0.99)
    save(fig, pathlib.Path(a.out_dir), f"bars_median_{a.panel}")


if __name__ == "__main__":
    sys.exit(main())
