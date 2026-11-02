#!/usr/bin/env python3
"""Provisional Figures 2-4: one message per figure, for discussion, not for print.

  2  embeddings add information beyond known mutations (M. tuberculosis)
  3  what the model relies on, and the co-resistance it has to avoid
  4  transfer to P. aeruginosa
  extended: feature sets and the call against TB-Profiler; the lineage-held-out
            result, the genome-wide ranking and accuracy by MIC class

Everything is read from the supplementary tables and saved predictions or SHAP
values; nothing is refit. Each figure carries a title and a sentence saying what
it is meant to show, which a journal figure would move to its legend.

Two panels are marked in the figure as not final:
  extended (limits) b  the genome-wide ranking is copied from
      EVALUATION_PROTOCOL.md section 1 and was computed on an early isoniazid
      partition; it has to be recomputed on the current development set
  4b  attribution by feature kind stands in for the nested feature-set runs,
      which are on the cluster and not among the downloaded results
"""

from __future__ import annotations

import argparse
import bisect
import collections
import csv
import pathlib
import re
import statistics
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from openpyxl import load_workbook

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import (AXIS, GRID, INK, INK2, MUTED, SERIES, SURFACE, header,  # noqa: E402
                                load, save, tidy)
from make_panel_grids import load_drug  # noqa: E402
from plot_class_distribution import (OT_LINE, TRAIN_EDGE, TRAIN_FACE, class_counts,  # noqa: E402
                                     draw_classdist)
from plot_metric_bars import PANELS, SHORT  # noqa: E402
from plot_mic_scatter import CATS, EDGE, FACE, LABEL  # noqa: E402
from plot_shap_beeswarm import REGION, draw_beeswarm, locus_label, normalise, rank_loci  # noqa: E402

HOME = pathlib.Path.home()
EMB = SERIES["evo2"]            # the reported model keeps its colour from every other figure
MUT = "#8a5a2b"                 # mutation columns + call: a model without embeddings
SEQ = LinearSegmentedColormap.from_list(
    "blue", ["#f3f7fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95"])
MTB8 = [d for d in PANELS["mtb"]["drugs"] if d not in PANELS["mtb"]["underpowered"]]
NAME = {"piperacillin-tazobactam": "Pip.–tazobactam"}
ABBR = {"isoniazid": "INH", "rifampicin": "RIF", "ethambutol": "EMB", "ethionamide": "ETH",
        "kanamycin": "KAN", "amikacin": "AMI", "levofloxacin": "LEV", "moxifloxacin": "MXF",
        "bedaquiline": "BDQ", "delamanid": "DLM"}
COL = {"cost": "Cost (2 x vMER + MER)", "vMER": "vMER", "MER": "MER", "F1": "F1 (resistant)",
       "OT": "Essential agreement",
       "balc": "Balanced essential agreement (classes with >= 10 isolates)",
       "MAE": "Mean absolute error (doubling dilutions)"}


def sheet(path, name):
    rows = list(load_workbook(path, read_only=True)[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def letter(ax, s, dx=-0.02, dy=1.05):
    ax.annotate(s, xy=(dx, dy), xycoords="axes fraction", fontsize=10, fontweight="bold",
                color=INK, ha="right", va="bottom", annotation_clip=False)


def paired_summary(labels, a, b, better_low, exclude=(), true_delta=False, unit=1.0):
    """How often and by how much the first series beats the second, across antimicrobials.

    The median is oriented so that positive is better, unless true_delta: then it is the median of
    first minus second, so for a lower-is-better measure a negative value favours the first series.
    """
    from scipy.stats import wilcoxon
    d = [((y - x) if better_low else (x - y)) for la, x, y in zip(labels, a, b) if la not in exclude]
    raw = [(x - y) for la, x, y in zip(labels, a, b) if la not in exclude]
    wins = sum(1 for v in d if v > 5e-4 * unit)
    ties = sum(1 for v in d if abs(v) <= 5e-4 * unit)
    pval = float(wilcoxon(d).pvalue) if any(abs(v) > 0 for v in d) else float("nan")
    return wins, ties, len(d), statistics.median(raw if true_delta else d), pval


def dumbbell(ax, labels, a, b, *, a_style, b_style, title, xlabel, better_low, faint=(), true_delta=False, unit=1.0):
    """One row per antimicrobial, two marks joined by a line, the difference at the right."""
    tidy(ax)
    lim = None
    if faint:
        # Scale the axis to the antimicrobials that are compared. A faded row whose
        # marks fall outside is named instead of squeezing the other eight.
        core = [v for la, x, y in zip(labels, a, b) if la not in faint for v in (x, y) if v is not None]
        pad = 0.14 * (max(core) - min(core))
        lim = (min(core) - pad, max(core) + pad)
    for i, (la, x, y) in enumerate(zip(labels, a, b)):
        if x is None or y is None:
            continue
        al = 0.4 if la in faint else 1.0
        if lim and la in faint and (min(x, y) < lim[0] or max(x, y) > lim[1]):
            ax.annotate(f"off scale: {y:.2f} and {x:.2f}" if unit == 1.0 else f"off scale: {y:.1f} and {x:.1f}", xy=(0.5, i), xycoords=("axes fraction", "data"),
                        ha="center", va="center", fontsize=5.6, color=MUTED)
        else:
            ax.plot([x, y], [i, i], color=GRID, lw=2.4, solid_capstyle="round", zorder=1, alpha=al)
            ax.plot(y, i, zorder=3, alpha=al, **b_style)
            ax.plot(x, i, zorder=4, alpha=al, **a_style)
        tol = 5e-4 * unit
        d = (y - x) if better_low else (x - y)          # positive = better, sets the colour
        shown = (x - y) if true_delta else d              # the number printed: first minus second, or oriented
        ax.annotate((f"{shown:+.3f}" if unit == 1.0 else f"{shown:+.1f}") if abs(shown) > tol else ("0.000" if unit == 1.0 else "0.0"), xy=(1.03, i),
                    xycoords=("axes fraction", "data"), fontsize=6.4,
                    color=(INK2 if d > tol else MUTED), ha="left", va="center", annotation_clip=False)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([NAME.get(la, la.capitalize()) for la in labels], fontsize=7.4, color=INK)
    ax.set_ylim(-0.7, len(labels) - 0.3)
    ax.invert_yaxis()
    ax.margins(x=0.10)
    if lim:
        ax.set_xlim(*lim)
    ax.tick_params(axis="x", labelsize=6.8)
    ax.set_title(title, loc="left", color=INK, pad=6, fontweight="medium", fontsize=8.4)
    ax.set_xlabel(xlabel, fontsize=7)


A_EMB = dict(marker="o", ms=5, color=EMB, mec=SURFACE, mew=0.8, ls="none")
A_MUT = dict(marker="s", ms=4.6, color=MUT, mec=SURFACE, mew=0.8, ls="none")
A_BASE = dict(marker="o", ms=5.6, mfc="none", mec=INK2, mew=1.3, ls="none")


def results(path, name):
    out = {}
    for r in sheet(path, name):
        out[(r["Antimicrobial"], r["Model"], r["Set"])] = r
    return out


# ── Figure 2 ──────────────────────────────────────────────────────────────────
def figure2(a, out):
    """The central comparison: the same model with embeddings or with mutation columns.

    Four measures on the held-out sets, one row per antimicrobial. Each panel states
    how many antimicrobials the embeddings win and the paired test, counted over the
    eight antimicrobials compared (bedaquiline and delamanid are drawn faded).
    """
    R = results(a.metrics, "MTB_results")
    emb, mut = "Evo 2 (reported)", "Mutation columns + TB-Profiler call (no embeddings)"
    drugs = PANELS["mtb"]["drugs"]
    faint = PANELS["mtb"]["underpowered"]

    def col(model, key):
        return [R[(d, model, "held-out")][COL[key]] for d in drugs]

    fig = plt.figure(figsize=(7.2, 7.0))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.34)
    panels = (("OT", "Essential agreement", "(higher is better)", False),
              ("balc", "Balanced essential agreement", "(higher is better)", False),
              ("MAE", "Mean absolute error", "doubling dilutions  (lower is better)", True),
              ("F1", "F1, resistant class", "(higher is better)", False))
    for n, (key, title, xl, low) in enumerate(panels):
        ax = fig.add_subplot(gs[n // 2, n % 2])
        x, y = col(emb, key), col(mut, key)
        dumbbell(ax, drugs, x, y, a_style=A_EMB, b_style=A_MUT, title=title, xlabel=xl,
                 better_low=low, faint=faint)
        wins, ties, k, med, pval = paired_summary(drugs, x, y, low, exclude=faint)
        note = f"embeddings better for {wins} of {k}" + (f", equal for {ties}" if ties else "")
        ax.annotate(f"{note}  ·  median {med:+.3f}  ·  P = {pval:.4f}", xy=(0.0, 1.0),
                    xycoords="axes fraction", xytext=(0, 3.5), textcoords="offset points",
                    ha="left", va="bottom", fontsize=6.3, color=INK2)
        ax.set_title(title, loc="left", color=INK, pad=15, fontweight="semibold", fontsize=8.8)
        if n % 2:
            ax.set_yticklabels([])
        letter(ax, "abcd"[n], dx=(-0.40 if n % 2 == 0 else -0.04), dy=1.10)
        print(f"[INFO] figure 2 {title}: {note}, median {med:+.4f}, P = {pval:.4f}")
    fig.legend(handles=[Line2D([], [], label="Embeddings + TB-Profiler call (reported model)", **A_EMB),
                        Line2D([], [], label="Mutation columns + TB-Profiler call", **A_MUT)],
               loc="lower center", ncol=2, fontsize=7.2, handletextpad=0.5, columnspacing=2.2,
               bbox_to_anchor=(0.5, 0.0))
    top = header(fig, "Figure 2 (provisional) — embeddings add information beyond known mutations",
                 "M. tuberculosis, held-out sets. The same random forest, given either sequence embeddings of "
                 "the target loci or the catalogued\nmutations as columns; both also receive the TB-Profiler "
                 "call. Numbers at the right: embeddings against mutation columns,\npositive = better. Counts, "
                 "medians and two-sided Wilcoxon tests are over the eight antimicrobials compared; bedaquiline\n"
                 "and delamanid (faded) have too few resistant isolates.", pad_in=0.30)
    fig.subplots_adjust(top=top, bottom=0.095, left=0.165, right=0.925)
    save(fig, out, "figure2_provisional")


def figure2_radar(a, out, rmin=0.0, stem="figure2_radar_provisional"):
    """The same comparison as figure 2, as one star plot per antimicrobial.

    Six spokes, each oriented so that further out is better: essential agreement,
    balanced essential agreement, F1, 1 - vMER, 1 - MER and 1 - MAE. The last is not
    a rate: MAE is in doubling dilutions and has no upper bound, so 1 - MAE only
    sits on this scale because every value here is below one dilution.

    `rmin` moves the centre of the plot away from zero. It stays at zero: with a
    raised centre, 1 - MAE (0.2 to 0.7 here) falls below it and collapses to the
    middle, which draws a notch that is not in the data.
    """
    R = results(a.metrics, "MTB_results")
    emb, mut = "Evo 2 (reported)", "Mutation columns + TB-Profiler call (no embeddings)"
    drugs = PANELS["mtb"]["drugs"]
    faint = PANELS["mtb"]["underpowered"]
    spokes = [("EA", "OT", False), ("Balanced\nEA", "balc", False), ("F1", "F1", False),
              ("1 − vMER", "vMER", True), ("1 − MER", "MER", True), ("1 − MAE", "MAE", True)]
    ang = np.linspace(0, 2 * np.pi, len(spokes), endpoint=False)

    def poly(model, d):
        v = [R[(d, model, "held-out")][COL[k]] for _, k, _ in spokes]
        return [(1 - x) if inv else x for x, (_, _, inv) in zip(v, spokes)]

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 4, hspace=0.50, wspace=0.62)
    for n, d in enumerate(drugs):
        ax = fig.add_subplot(gs[n // 4, n % 4], projection="polar")
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        al = 0.45 if d in faint else 1.0
        for model, c, fill, lw, z in ((mut, MUT, False, 1.1, 3), (emb, EMB, True, 1.3, 4)):
            v = [max(x, rmin) for x in poly(model, d)]
            th = np.r_[ang, ang[0]]
            vv = np.r_[v, v[0]]
            ax.plot(th, vv, color=c, lw=lw, zorder=z, alpha=al, solid_joinstyle="round")
            if fill:
                ax.fill(th, vv, color=c, alpha=0.16 * al, zorder=2, lw=0)
        ax.set_ylim(rmin, 1.0)
        ticks = [t for t in (0.2, 0.4, 0.6, 0.8, 1.0) if t > rmin + 1e-9]
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{t:.1f}" if t < 1 else "" for t in ticks], fontsize=4.2, color=MUTED)
        ax.set_rlabel_position(0)
        ax.set_xticks(ang)
        ax.set_xticklabels([sp[0] for sp in spokes], fontsize=5.2, color=INK2, linespacing=0.95)
        ax.tick_params(axis="x", pad=-2.5)
        ax.grid(color=GRID, lw=0.5)
        ax.spines["polar"].set_color(AXIS)
        ax.spines["polar"].set_linewidth(0.7)
        ax.set_title(d.capitalize() + (" *" if d in faint else ""), fontsize=8.0, fontweight="semibold",
                     color=INK, pad=11)
    key = fig.add_subplot(gs[2, 2:4])
    key.axis("off")
    key.legend(handles=[Line2D([], [], color=EMB, lw=1.6, label="Embeddings + TB-Profiler call\n(reported model)"),
                        Line2D([], [], color=MUT, lw=1.4, label="Mutation columns + TB-Profiler call")],
               loc="upper left", bbox_to_anchor=(0.02, 0.98), fontsize=7.0, labelspacing=0.9, handlelength=1.6)
    key.annotate("Every spoke points outward for better.\n1 − vMER and 1 − MER: share of resistant and of\n"
                 "susceptible isolates called correctly.\n1 − MAE: one minus the mean error in doubling\n"
                 "dilutions (not a rate).\n* few resistant held-out isolates."
                 + (f"\nThe centre of each plot is {rmin:.1f}, not zero." if rmin else ""),
                 xy=(0.02, 0.50), xycoords="axes fraction", ha="left", va="top", fontsize=6.0, color=INK2,
                 linespacing=1.45)
    top = header(fig, "Figure 2, alternative (provisional) — embeddings against mutation columns, per antimicrobial",
                 "M. tuberculosis, held-out sets. One star per antimicrobial with six measures; the same random "
                 "forest given embeddings of the\ntarget loci (blue) or the catalogued mutations as columns "
                 "(brown), both with the TB-Profiler call."
                 + (f" Radial axis from {rmin:.1f} to 1." if rmin else " Radial axis from 0 to 1."), pad_in=0.34)
    fig.subplots_adjust(top=top, bottom=0.03, left=0.07, right=0.95)
    save(fig, out, stem)


def figure2_extended(a, out):
    R = results(a.metrics, "MTB_results")
    emb, mut, tbp = "Evo 2 (reported)", "Mutation columns + TB-Profiler call (no embeddings)", "TB-Profiler"
    drugs = PANELS["mtb"]["drugs"]
    faint = PANELS["mtb"]["underpowered"]

    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 6, hspace=0.66, wspace=1.5, height_ratios=[0.86, 1.0])

    # b: the six feature sets on development folds, one point per antimicrobial
    F = collections.defaultdict(dict)
    for r in sheet(a.metrics, "MTB_feature_sets"):
        F[r["Feature set"]][r["Antimicrobial"]] = r
    order = ["Embeddings + TB-Profiler call", "Embeddings + mutation columns + TB-Profiler call",
             "Mutation columns + TB-Profiler call", "Embeddings + mutation columns",
             "Mutation columns only", "Embeddings only"]
    for j, (key, title, xl) in enumerate((("cost", "Cost", "2·vMER + MER  (lower is better)"),
                                          ("balc", "Balanced essential agreement", "(higher is better)"))):
        ax = fig.add_subplot(gs[0, 0:3] if j == 0 else gs[0, 3:6])
        tidy(ax)
        for i, fs in enumerate(order):
            v = [F[fs][d][COL[key]] for d in MTB8 if F[fs].get(d) and F[fs][d][COL[key]] is not None]
            c = EMB if i == 0 else INK2
            ax.plot(v, [i] * len(v), "o", ms=2.8, mfc=SURFACE, mec=c, mew=0.6, zorder=3, alpha=0.9)
            m = statistics.median(v)
            ax.plot([m, m], [i - 0.32, i + 0.32], color=c, lw=2.0, zorder=4, solid_capstyle="round")
            ax.annotate(f"{m:.3f}", xy=(1.03, i), xycoords=("axes fraction", "data"), fontsize=6.4,
                        color=INK2, ha="left", va="center", annotation_clip=False)
        ax.set_yticks(range(len(order)))
        short = ["Embeddings + call", "Emb. + mutations + call", "Mutations + call",
                 "Emb. + mutations", "Mutations only", "Embeddings only"]
        ax.set_yticklabels(short if j == 0 else [], fontsize=6.8, color=INK)
        ax.set_ylim(-0.7, len(order) - 0.3)
        ax.invert_yaxis()
        ax.tick_params(axis="x", labelsize=6.8)
        ax.set_title(title + ", development folds", loc="left", color=INK, pad=6, fontweight="medium",
                     fontsize=8.4)
        ax.set_xlabel(xl, fontsize=7)
        if j == 0:
            letter(ax, "a", dx=-0.44)

    # c: the resistance call against TB-Profiler
    for j, (key, title, xl, low) in enumerate((("vMER", "Very major error rate", "resistant called susceptible", True),
                                               ("MER", "Major error rate", "susceptible called resistant", True),
                                               ("F1", "F1, resistant class", "(higher is better)", False))):
        ax = fig.add_subplot(gs[1, 2 * j:2 * j + 2])
        dumbbell(ax, MTB8, [R[(d, emb, "held-out")][COL[key]] for d in MTB8],
                 [R[(d, tbp, "held-out")][COL[key]] for d in MTB8], a_style=A_EMB, b_style=A_BASE,
                 title=title, xlabel=xl, better_low=low)
        if j:
            ax.set_yticklabels([])
        else:
            letter(ax, "b", dx=-0.55)
    fig.legend(handles=[Line2D([], [], label="Embeddings + TB-Profiler call (reported)", **A_EMB),
                        Line2D([], [], label="TB-Profiler alone", **A_BASE),
                        Line2D([], [], marker="o", ls="none", ms=2.8, mfc=SURFACE, mec=INK2, mew=0.6,
                               label="one antimicrobial (a)"),
                        Line2D([], [], color=INK2, lw=2.0, label="median (a)")],
               loc="lower center", ncol=3, fontsize=6.8, handletextpad=0.5, columnspacing=1.6,
               bbox_to_anchor=(0.5, 0.0))
    top = header(fig, "Extended figure (provisional) — which inputs the model needs, and the resistance call",
                 "M. tuberculosis. a, development folds: six feature sets; every one without the TB-Profiler call "
                 "is worse on cost.\nb, held-out sets: the reported model against TB-Profiler alone; the gain is "
                 "fewer major errors. Numbers at the right of b:\nreported model against TB-Profiler, positive = "
                 "better.", pad_in=0.26)
    fig.subplots_adjust(top=top, bottom=0.085, left=0.185, right=0.935)
    save(fig, out, "figure2_extended_provisional")


# ── panels shared by figure 3 and the extended figure ─────────────────────────
def auroc(pos, neg):
    neg = sorted(neg)
    s = 0.0
    for p in pos:
        lo, hi = bisect.bisect_left(neg, p), bisect.bisect_right(neg, p)
        s += lo + 0.5 * (hi - lo)
    return s / (len(pos) * len(neg))


def mtb_phenotypes(a):
    iso = sheet(a.isolates, "MTB_isolates")
    cut = {r["antimicrobial"]: r["resistant above (log2 MIC)"] for r in sheet(a.isolates, "MTB_counts")}
    drugs = PANELS["mtb"]["drugs"]
    mic = {d: {r["isolate"]: r[f"{d} log2 MIC"] for r in iso if r[f"{d} log2 MIC"] is not None} for d in drugs}
    res = {d: {i: v > cut[d] for i, v in mic[d].items()} for d in drugs}
    return drugs, mic, res


def p_shap(ax, a, drug="isoniazid", top=12):
    long = pd.read_csv(pathlib.Path(a.mtb_shap) / f"{drug}_evo2_shap_long.tsv", sep="\t").fillna({"gene": ""})
    long = normalise(long)
    order, meta = rank_loci(long, top=top)
    draw_beeswarm(ax, long[long["locus"].isin(order)], order, s=2.4, alpha=0.4)
    ax.set_yticklabels([locus_label(meta, lo).replace("TB-Profiler rule bit", "TB-Profiler call")
                        for lo in order], fontsize=6.6)
    ax.tick_params(axis="x", labelsize=6.8)
    ax.set_xlabel("SHAP value  (log₂ MIC)", fontsize=7.2)
    ax.set_title(f"{drug.capitalize()}: the {top} highest-ranked loci", loc="left", fontsize=8.4, color=INK,
                 pad=6, fontweight="medium")
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=4.2, color=REGION[r][0], label=lab)
                       for r, lab in (("gene", "gene body"), ("upstream", "300 bp upstream"),
                                      ("other", "TB-Profiler call"))],
              loc="upper center", bbox_to_anchor=(0.36, -0.13), ncol=3, fontsize=6.2, handletextpad=0.1,
              columnspacing=0.9)


SHARE = [("TB-Profiler call", "#8a8880"), ("loci in the TB-Profiler database", "#6b3fa0"),
         ("loci from the association study", "#e8912b")]


def p_share(ax, a):
    """Where the attribution sits, per antimicrobial."""
    tidy(ax)
    drugs = PANELS["mtb"]["drugs"]
    for i, d in enumerate(drugs):
        rows = list(csv.DictReader(open(pathlib.Path(a.mtb_shap) / f"{d}_evo2_target_summary.tsv"),
                                   delimiter="\t"))
        early = {r["target_id"] for r in csv.DictReader(open(pathlib.Path(a.mtb_earlier) /
                                                             f"{d}_evo2_target_summary.tsv"), delimiter="\t")}
        tot = sum(float(r["mean_abs_shap"]) for r in rows)
        share = [sum(float(r["mean_abs_shap"]) for r in rows if r["target_id"].startswith("rule")),
                 sum(float(r["mean_abs_shap"]) for r in rows if not r["target_id"].startswith("rule")
                     and r["target_id"] in early),
                 sum(float(r["mean_abs_shap"]) for r in rows if not r["target_id"].startswith("rule")
                     and r["target_id"] not in early)]
        x = 0.0
        for (lab, c), v in zip(SHARE, share):
            ax.barh(i, v / tot, left=x, height=0.62, color=c, edgecolor=SURFACE, linewidth=0.8, zorder=3,
                    alpha=0.5 if d in PANELS["mtb"]["underpowered"] else 1.0)
            if v / tot > 0.09:
                ax.annotate(f"{v / tot:.0%}", (x + v / tot / 2, i), ha="center", va="center", fontsize=5.8,
                            color=SURFACE if c != "#e8912b" else INK, zorder=4)
            x += v / tot
    ax.set_yticks(range(len(drugs)))
    ax.set_yticklabels([d.capitalize() + (" *" if d in PANELS["mtb"]["underpowered"] else "") for d in drugs],
                       fontsize=7.2, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25%", "50%", "75%", "100%"], fontsize=6.8)
    ax.set_xlabel("share of mean |SHAP|", fontsize=7.2)
    ax.set_title("Where the attribution sits", loc="left", fontsize=8.4, color=INK, pad=6, fontweight="medium")
    ax.legend(handles=[Patch(facecolor=c, edgecolor=SURFACE, label=lab) for lab, c in SHARE],
              loc="upper center", bbox_to_anchor=(0.42, -0.13), ncol=1, fontsize=6.4, handlelength=1.1,
              labelspacing=0.3)


def p_heat(fig, ax, a, ph):
    """Resistance to one antimicrobial predicted from the MIC of another."""
    drugs, mic, res = ph
    M = np.full((len(drugs), len(drugs)), np.nan)
    for i, da in enumerate(drugs):
        for j, db in enumerate(drugs):
            if da == db:
                continue
            ids = [x for x in res[da] if x in mic[db]]
            pos = [mic[db][x] for x in ids if res[da][x]]
            neg = [mic[db][x] for x in ids if not res[da][x]]
            if len(pos) >= 10 and neg:
                M[i, j] = auroc(pos, neg)
    im = ax.imshow(M, cmap=SEQ, vmin=0.5, vmax=1.0)
    for i in range(len(drugs)):
        for j in range(len(drugs)):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{M[i, j]:.2f}"[1:], ha="center", va="center", fontsize=5.2,
                        color=SURFACE if M[i, j] > 0.82 else INK)
    ax.set_xticks(range(len(drugs)))
    ax.set_xticklabels([ABBR[d] for d in drugs], fontsize=6.2, color=INK)
    ax.set_yticks(range(len(drugs)))
    ax.set_yticklabels([d.capitalize() for d in drugs], fontsize=6.8, color=INK)
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.set_xlabel("…from the MIC of", fontsize=7.2)
    ax.set_title("Resistance predicted from another\nantimicrobial's MIC (AUROC)", loc="left", fontsize=8.4,
                 color=INK, pad=6, fontweight="medium")
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.ax.tick_params(labelsize=6.2, length=0)
    cb.outline.set_visible(False)


def p_sites(ax, ph):
    """How often isoniazid resistance comes with rifampicin resistance, by collection site."""
    _, _, res = ph
    tidy(ax, x=False, y=True)
    site = collections.defaultdict(list)
    for i, r in res["isoniazid"].items():
        if r and i in res["rifampicin"]:
            site[re.match(r"site\.([^.]+)", i).group(1)].append(res["rifampicin"][i])
    big = sorted(((sum(v) / len(v), s_, len(v)) for s_, v in site.items() if len(v) >= 50))
    for x, (p, s_, n) in enumerate(big):
        ax.bar(x, p, 0.64, color=INK2, edgecolor=SURFACE, linewidth=0.5, zorder=3)
        ax.annotate(f"{p:.0%}", (x, p), xytext=(0, 2), textcoords="offset points", ha="center", va="bottom",
                    fontsize=6.2, color=INK2)
    ax.set_xticks(range(len(big)))
    ax.set_xticklabels([f"{s_}\nn={n}" for _, s_, n in big], fontsize=6.0, color=INK, linespacing=1.2)
    ax.set_ylim(0, 1.08)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"], fontsize=6.8)
    ax.set_xlabel("collection site (isoniazid-resistant isolates)", fontsize=7.2)
    ax.set_ylabel("also rifampicin-resistant", fontsize=7.2)
    ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    ax.set_title("Co-resistance differs between sites", loc="left", fontsize=8.4, color=INK, pad=6,
                 fontweight="medium")


def p_rank(ax):
    """Genome-wide relevance for isoniazid, all isolates against rifampicin-susceptible."""
    tidy(ax, x=False, y=True)
    loci = [("rpoB", 0.616, 0.148, MUT), ("rpsL", 0.417, 0.262, MUT), ("rpoC", 0.341, 0.077, MUT),
            ("katG", 0.153, 0.182, EMB), ("inhA promoter", 0.096, 0.115, EMB)]
    for name, v0, v1, c in loci:
        ax.plot([0, 1], [v0, v1], color=c, lw=1.4, marker="o", ms=4.5, mec=SURFACE, mew=0.7, zorder=3)
        ax.annotate(name, (0, v0), xytext=(-6, 0), textcoords="offset points", ha="right", va="center",
                    fontsize=6.8, color=c, fontstyle="italic")
        ax.annotate(f"{v1:.2f}", (1, v1), xytext=(6, 0), textcoords="offset points", ha="left", va="center",
                    fontsize=6.4, color=c)
    ax.set_xlim(-0.55, 1.35)
    ax.set_ylim(0, 0.68)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["all isolates", "rifampicin-\nsusceptible only"], fontsize=7, color=INK)
    ax.set_ylabel("|r| with isoniazid log₂ MIC", fontsize=7.2)
    ax.tick_params(axis="y", labelsize=6.8)
    ax.set_title("Genome-wide selection, isoniazid", loc="left", fontsize=8.4, color=INK, pad=6,
                 fontweight="medium")
    ax.annotate("NOT FINAL: from an early partition,\nto be recomputed", xy=(0.98, 0.97), xycoords="axes fraction",
                ha="right", va="top", fontsize=5.8, color="#d03b3b")
    ax.legend(handles=[Line2D([], [], color=MUT, lw=1.4, marker="o", ms=4, label="locus of another antimicrobial"),
                       Line2D([], [], color=EMB, lw=1.4, marker="o", ms=4, label="isoniazid locus")],
              loc="upper center", bbox_to_anchor=(0.5, -0.20), fontsize=6.4, ncol=1, labelspacing=0.3)


def p_lineage(ax, a):
    """The advantage of embeddings over mutation columns, random against lineage folds."""
    L = [r for r in sheet(a.metrics, "MTB_lineage")
         if r["Comparison"].startswith("Reported-type") and r["Measure"] == COL["balc"]]
    names = [r["Antimicrobial"] for r in L]
    adv_r = [r["First model, size-matched random folds"] - r["Second model, size-matched random folds"] for r in L]
    adv_l = [r["First model, lineage folds"] - r["Second model, lineage folds"] for r in L]
    tidy(ax)
    for i, (x0, x1) in enumerate(zip(adv_r, adv_l)):
        ax.plot([x0, x1], [i, i], color=GRID, lw=2.4, solid_capstyle="round", zorder=1)
        ax.plot(x0, i, "o", ms=5, mfc=SURFACE, mec=EMB, mew=1.3, zorder=3)
        ax.plot(x1, i, "o", ms=5, color=EMB, mec=SURFACE, mew=0.8, zorder=4)
    ax.axvline(0, color=INK2, lw=0.8, zorder=2)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([n.capitalize() for n in names], fontsize=7.2, color=INK)
    ax.set_ylim(-0.7, len(names) - 0.3)
    ax.invert_yaxis()
    ax.tick_params(axis="x", labelsize=6.8)
    ax.set_xlabel("balanced essential agreement, embeddings minus mutation columns\n"
                  "left of zero: mutation columns better      right of zero: embeddings better",
                  fontsize=6.6, linespacing=1.5)
    ax.set_title("A lineage held out", loc="left", fontsize=8.4, color=INK, pad=6, fontweight="medium")
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=5, mfc=SURFACE, mec=EMB, mew=1.3,
                              label="random folds of the same sizes"),
                       Line2D([], [], marker="o", ls="none", ms=5, color=EMB, mec=SURFACE, mew=0.8,
                              label="one lineage held out per fold")],
              loc="upper center", bbox_to_anchor=(0.5, -0.27), fontsize=6.4, ncol=1, labelspacing=0.3)


def p_classdist(ax, a, msum, drug, bound):
    pred_dir = pathlib.Path(a.mtb_pred)
    y, pred, cat, m, folds_p = load_drug(pred_dir, pred_dir, drug, "evo2", bound, msum, 5e-4)
    ks, tr, stacks, hold, ot, n_train = class_counts(folds_p, y, cat, bound)
    ax2 = draw_classdist(ax, ks, tr, stacks, hold, ot, bound, annotate=False, fs=0.72)
    ax.tick_params(labelsize=6.6)
    ax2.tick_params(labelsize=6.6)
    ax.set_xlabel("log₂ MIC class", fontsize=7.2)
    ax.set_ylabel("Isolates", fontsize=7.2)
    ax2.set_ylabel("Within one dilution", color=OT_LINE, fontsize=7.2)
    ax.set_title(f"{drug.capitalize()}: accuracy by MIC class", loc="left", fontsize=8.4, color=INK,
                 pad=6, fontweight="medium")


# ── Figure 3 ──────────────────────────────────────────────────────────────────
def figure3(a, out):
    """What the model relies on, and what else the data would let it rely on."""
    ph = mtb_phenotypes(a)
    fig = plt.figure(figsize=(7.2, 8.9))
    gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.52, height_ratios=[1.18, 1.0])
    ax = fig.add_subplot(gs[0, 0])
    p_shap(ax, a)
    letter(ax, "a", dx=-0.56)
    ax = fig.add_subplot(gs[0, 1])
    p_share(ax, a)
    letter(ax, "b", dx=-0.34)
    ax = fig.add_subplot(gs[1, 0])
    p_heat(fig, ax, a, ph)
    letter(ax, "c", dx=-0.44)
    ax = fig.add_subplot(gs[1, 1])
    p_sites(ax, ph)
    letter(ax, "d", dx=-0.22)
    top = header(fig, "Figure 3 (provisional) — what the model relies on, and the co-resistance it has to avoid",
                 "M. tuberculosis. a, TreeSHAP for isoniazid on held-out isolates: one point per isolate; in each "
                 "row the gene body is above and the\nupstream region below. b, share of the attribution carried "
                 "by the catalogue call, by loci the catalogue annotates and by\nassociation-derived loci "
                 "(* few resistant isolates). c, resistance to one antimicrobial is largely predictable from the "
                 "MIC of another.\nd, how often isoniazid resistance comes with rifampicin resistance differs "
                 "between collection sites.", pad_in=0.30)
    fig.subplots_adjust(top=top, bottom=0.075, left=0.205, right=0.95)
    save(fig, out, "figure3_provisional")


def figure_extended_limits(a, out):
    """For the supplement for now: the lineage result, the genome-wide ranking, accuracy by MIC class."""
    msum, _ = load(pathlib.Path(a.mtb_workbook))
    fig = plt.figure(figsize=(7.2, 8.6))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.50, height_ratios=[1.0, 0.95])
    ax = fig.add_subplot(gs[0, 0])
    p_lineage(ax, a)
    letter(ax, "a", dx=-0.36)
    ax = fig.add_subplot(gs[0, 1])
    p_rank(ax)
    letter(ax, "b", dx=-0.30)
    for j, (drug, bound) in enumerate((("isoniazid", -3.0), ("ethambutol", 2.0))):
        ax = fig.add_subplot(gs[1, j])
        p_classdist(ax, a, msum, drug, bound)
        if j == 0:
            letter(ax, "c", dx=-0.18)
    fig.legend(handles=[Patch(facecolor=TRAIN_FACE, edgecolor=TRAIN_EDGE, lw=0.6, label="development set (c)")]
               + [Patch(facecolor=FACE[k], edgecolor=EDGE[k], lw=0.5, label=LABEL[k] + " (c)")
                  for k, _, _, _ in CATS[:2]]
               + [Line2D([], [], color=OT_LINE, lw=1.3, marker="o", ms=3.5, label="class accuracy (c, right axis)")],
               loc="lower center", ncol=4, fontsize=6.4, handletextpad=0.5, columnspacing=1.4,
               bbox_to_anchor=(0.5, 0.0))
    top = header(fig, "Extended figure (provisional) — limits: a lineage held out, and genome-wide selection",
                 "M. tuberculosis. a, development set: the advantage of embeddings over mutation columns turns "
                 "negative when a whole lineage is\nheld out. b, selecting loci genome-wide picks markers of "
                 "resistance to other antimicrobials; they fall away among rifampicin-\nsusceptible isolates. "
                 "c, isolates per MIC class and the share predicted within one dilution, held-out sets.", pad_in=0.30)
    fig.subplots_adjust(top=top, bottom=0.07, left=0.16, right=0.93)
    save(fig, out, "figure_extended_limits_provisional")


# ── Figure 4 ──────────────────────────────────────────────────────────────────
def figure4(a, out):
    R = results(a.metrics, "PA_results")
    emb, base = "Evo 2 (reported)", "ARDaP (intermediate counted as resistant)"
    drugs = PANELS["pseaer"]["drugs"]
    cov = [d for d in drugs if (d, base, "held-out") in R]
    fig = plt.figure(figsize=(7.2, 9.0))
    gs = fig.add_gridspec(3, 6, hspace=0.62, wspace=1.6, height_ratios=[0.86, 1.15, 0.95])
    for j, (key, title, xl, low) in enumerate((("vMER", "Very major error rate", "resistant called susceptible", True),
                                               ("MER", "Major error rate", "susceptible called resistant", True),
                                               ("F1", "F1, resistant class", "(higher is better)", False))):
        ax = fig.add_subplot(gs[0, 2 * j:2 * j + 2])
        dumbbell(ax, cov, [R[(d, emb, "held-out")][COL[key]] for d in cov],
                 [R[(d, base, "held-out")][COL[key]] for d in cov], a_style=A_EMB, b_style=A_BASE,
                 title=title, xlabel=xl, better_low=low)
        if j:
            ax.set_yticklabels([])
        else:
            letter(ax, "a", dx=-0.62)

    # b: attribution by kind of feature
    shap_dir = pathlib.Path(a.pa_shap)
    kinds = ["gene", "upstream", "acquired", "mutation", "presence"]
    ax = fig.add_subplot(gs[1, 0:3])
    tidy(ax)
    for i, d in enumerate(drugs):
        long = normalise(pd.read_csv(shap_dir / f"{d}_evo2_shap_long.tsv", sep="\t").fillna({"gene": ""}))
        per = long.assign(a=long["shap"].abs()).groupby(["target_id", "region"], as_index=False)["a"].mean()
        tot = per["a"].sum()
        x = 0.0
        for k in kinds:
            v = per.loc[per["region"] == k, "a"].sum() / tot
            ax.barh(i, v, left=x, height=0.62, color=REGION[k][0], edgecolor=SURFACE, linewidth=0.8, zorder=3)
            if v > 0.10:
                ax.annotate(f"{v:.0%}", (x + v / 2, i), ha="center", va="center", fontsize=5.8,
                            color=SURFACE if k not in ("upstream", "presence") else INK, zorder=4)
            x += v
    ax.set_yticks(range(len(drugs)))
    ax.set_yticklabels([NAME.get(d, d.capitalize()) for d in drugs], fontsize=7.0, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25%", "50%", "75%", "100%"], fontsize=6.8)
    ax.set_xlabel("share of mean |SHAP|", fontsize=7.2)
    ax.set_title("Attribution by kind of feature", loc="left", fontsize=8.4, color=INK, pad=6,
                 fontweight="medium")
    ax.annotate("NOT FINAL: stands in for the nested feature-set comparison", xy=(0.5, 0.0),
                xycoords="axes fraction", xytext=(0, -30), textcoords="offset points", ha="center", va="top",
                fontsize=5.6, color="#d03b3b")
    letter(ax, "b", dx=-0.36)

    # c: ciprofloxacin, what the model used
    long = normalise(pd.read_csv(shap_dir / "ciprofloxacin_evo2_shap_long.tsv", sep="\t").fillna({"gene": ""}))
    order, meta = rank_loci(long, top=10)
    ax = fig.add_subplot(gs[1, 3:6])
    draw_beeswarm(ax, long[long["locus"].isin(order)], order, s=2.4, alpha=0.45)
    ax.set_yticklabels([locus_label(meta, lo).split(" (")[0] for lo in order], fontsize=6.6)
    ax.tick_params(axis="x", labelsize=6.8)
    ax.set_xlabel("SHAP value  (log₂ MIC)", fontsize=7.2)
    ax.set_title("Ciprofloxacin, ten highest-ranked", loc="left", fontsize=8.4, color=INK, pad=6,
                 fontweight="medium")
    letter(ax, "c", dx=-0.26)

    # d: each study held out against cross-validation on the same isolates
    S = collections.defaultdict(dict)
    for r in sheet(a.metrics, "PA_study_heldout"):
        S[(r["Antimicrobial"], r["Study"])][r["Evaluation"]] = r
    for j, (key, title) in enumerate((("OT", "Essential agreement"), ("cost", "Cost"))):
        ax = fig.add_subplot(gs[2, 0:3] if j == 0 else gs[2, 3:6])
        tidy(ax, x=True, y=True)
        xs = [v["Standard cross-validation, same isolates"][COL[key]] for v in S.values()]
        ys = [v["Study held out"][COL[key]] for v in S.values()]
        ok = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
        lo, hi = (0.2, 1.0) if key == "OT" else (0.0, max(max(p) for p in ok) * 1.06)
        ax.plot([lo, hi], [lo, hi], color=AXIS, lw=0.9, ls=(0, (3, 3)), zorder=1)
        ax.plot([p[0] for p in ok], [p[1] for p in ok], "o", ms=3.6, color=EMB, mec=SURFACE, mew=0.6,
                alpha=0.85, zorder=3)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=6.8)
        ax.set_xlabel("standard cross-validation, same isolates", fontsize=7.2)
        ax.set_ylabel("study held out", fontsize=7.2)
        worse = sum(1 for x, y in ok if (y < x - 1e-9 if key == "OT" else y > x + 1e-9))
        ax.set_title(f"{title}: worse when held out for {worse} of {len(ok)}", loc="left", fontsize=8.0,
                     color=INK, pad=6, fontweight="medium")
        if j == 0:
            letter(ax, "d", dx=-0.30)
    present = [k for k in kinds]
    fig.legend(handles=[Line2D([], [], label="Reported model (Evo 2)", **A_EMB),
                        Line2D([], [], label="ARDaP", **A_BASE)]
               + [Patch(facecolor=REGION[k][0], edgecolor=SURFACE, label=REGION[k][1]) for k in present]
               + [Line2D([], [], marker="o", ls="none", ms=3.6, color=EMB, mec=SURFACE, mew=0.6,
                         label="one study and antimicrobial (d)")],
               loc="lower center", ncol=4, fontsize=6.4, handletextpad=0.5, columnspacing=1.4,
               bbox_to_anchor=(0.5, 0.0))
    top = header(fig, "Figure 4 (provisional) — transfer to P. aeruginosa",
                 "Held-out sets unless stated. a, the resistance call against ARDaP for the eight "
                 "antimicrobials it covers; numbers: reported model\nagainst ARDaP, positive = better. b, share "
                 "of the attribution by kind of feature. c, TreeSHAP for ciprofloxacin. d, development set:\n"
                 "each study held out in turn against cross-validation on the same isolates; below the diagonal "
                 "(accuracy) or above it (cost) is worse.", pad_in=0.26)
    fig.subplots_adjust(top=top, bottom=0.095, left=0.17, right=0.925)
    save(fig, out, "figure4_provisional")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", default=str(HOME / "Desktop/AMRPredict_comparison_metrics.xlsx"))
    ap.add_argument("--isolates", default=str(HOME / "Desktop/AMRPredict_isolates_used.xlsx"))
    ap.add_argument("--mtb-workbook", default=str(HOME / "Desktop/AMRPredict_panel_results_wider.xlsx"))
    ap.add_argument("--mtb-pred", default=str(HOME / "Downloads/mtb_wider_figure_inputs"))
    ap.add_argument("--mtb-shap", default=str(HOME / "Downloads/mtb_shap_gwfrule1k/shap_gwfrule1k"))
    ap.add_argument("--mtb-earlier", default=str(HOME / "Downloads"),
                    help="target summaries of the catalogue-derived model, to tell the two sources apart")
    ap.add_argument("--pa-shap", default=str(HOME / "Downloads/pa_strict/predictions_strict"))
    ap.add_argument("--only", default="", help="comma-separated figure numbers")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out = pathlib.Path(a.out_dir)
    todo = {int(x) for x in a.only.split(",")} if a.only else {2, 3, 4}
    for n, fn in ((2, figure2), (3, figure3), (4, figure4)):
        if n in todo:
            fn(a, out)
            if n == 3:
                figure_extended_limits(a, out)
            if n == 2:
                figure2_extended(a, out)      # figure2_radar(a, out) is kept but not drawn:
                                              # at this effect size the two outlines coincide
    return 0


if __name__ == "__main__":
    sys.exit(main())
