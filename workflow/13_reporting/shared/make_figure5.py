#!/usr/bin/env python3
"""Figure 5 (draft): co-resistance, why the loci are chosen for each antimicrobial.

  a, b   one heatmap per organism. Each cell is the area under the ROC curve for telling the isolates
         resistant to the antimicrobial of the row from those susceptible to it, using only the log2 MIC
         of the antimicrobial of the column, among isolates with both MICs and at least ten resistant
         to the row's antimicrobial. Pairs inside one drug class (outlined) share a mechanism, so a high
         value there is expected; a high value outside the outlines is co-resistance in the isolates.
  c      M. tuberculosis: among isoniazid-resistant isolates, the fraction with a rifampicin-associated rpoB
         variant called by TB-Profiler (bars) and the fraction that are rifampicin-resistant by MIC
         (circles), by collection site (sites with at least 50 isoniazid-resistant isolates).

Every value is computed from the isolate workbook (Supplementary Table 1): the MICs of every isolate and
the resistance cut-off of each antimicrobial. Nothing is refit. The matrices are written to
figure5_values.tsv beside the figure.
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import pathlib
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle
from openpyxl import load_workbook
from scipy.stats import rankdata

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_figure3_alt import CODE  # noqa: E402
from make_panel_figures import (SUBMISSION, fig_h, AXIS, INK, INK2, ITALIC, ITALIC_BOLD, MUTED, ORG_COLOUR, SURFACE, header,  # noqa: E402
                                save, tidy)
from make_provisional_figures import SEQ  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# the P. aeruginosa heatmap in the colour of the organism: white to a dark coral-pink
SEQ_PA = LinearSegmentedColormap.from_list("coral", ["#fdf3f6", "#f4b8c9", "#d9577f", "#a02a55", "#6b1838"])
from plot_metric_bars import PANELS  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
DARK_RED = "#8c1d2c"
NAME = {"piperacillin-tazobactam": "Pip.–tazobactam"}
# drug classes with more than one member: pairs inside them share a mechanism of resistance
CLASSES = {"mtb": [("kanamycin", "amikacin"), ("levofloxacin", "moxifloxacin")],
           "pseaer": [("ciprofloxacin", "levofloxacin"), ("tobramycin", "amikacin"), ("imipenem", "meropenem"),
                      ("ceftazidime", "cefepime", "aztreonam", "piperacillin-tazobactam")]}


def sheet(wb, name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def phenotypes(wb, k):
    """MIC and resistance of every isolate, per antimicrobial."""
    iso_sheet, cnt_sheet, idc = (("MTB_isolates", "MTB_counts", "isolate") if k == "mtb"
                                 else ("PA_isolates", "PA_counts", "run accession"))
    cut = {r["antimicrobial"]: r["resistant above (log2 MIC)"] for r in sheet(wb, cnt_sheet)}
    drugs = PANELS[k]["drugs"]
    iso = sheet(wb, iso_sheet)
    mic = {d: {r[idc]: r[f"{d} log2 MIC"] for r in iso if r.get(f"{d} log2 MIC") is not None} for d in drugs}
    res = {d: {i: v > cut[d] for i, v in mic[d].items()} for d in drugs}
    return drugs, mic, res


def rpob_carriers(path):
    """Isolates with a TB-Profiler variant in rpoB that the catalogue links to rifampicin."""
    with open(path, newline="") as fh:
        return {r["sample"] for r in csv.DictReader(fh)
                if r["gene_name"] == "rpoB" and "rifampicin" in re.split(r"[;,]", r["drugs"])}


def auroc(pos, neg):
    """Probability that a random resistant isolate has a higher MIC than a random susceptible one."""
    x = np.asarray(list(pos) + list(neg), float)
    r = rankdata(x)
    n1, n0 = len(pos), len(neg)
    return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def matrix(drugs, mic, res):
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
    return M


def heat(ax, drugs, M, classes, cmap=SEQ):
    im = ax.imshow(M, cmap=cmap, vmin=0.5, vmax=1.0)
    for i in range(len(drugs)):
        for j in range(len(drugs)):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{100 * M[i, j]:.1f}", ha="center", va="center", fontsize=4.5,
                        color=SURFACE if M[i, j] > 0.82 else INK)
    for cls in classes:                                    # outline the pairs inside one drug class
        idx = [drugs.index(d) for d in cls]
        lo, n = min(idx), len(idx)
        ax.add_patch(Rectangle((lo - 0.5, lo - 0.5), n, n, fill=False, edgecolor=INK, lw=1.1, zorder=5))
    ax.set_xticks(range(len(drugs)))
    ax.set_xticklabels([CODE[d] for d in drugs], fontsize=6.0, color=INK)
    ax.set_yticks(range(len(drugs)))
    ax.set_yticklabels([NAME.get(d, d.capitalize()) for d in drugs], fontsize=6.6, color=INK)
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.set_xlabel("… from the MIC of", fontsize=7.0)
    ax.set_ylabel("Resistance to …", fontsize=7.0)
    return im


def sites(ax, mic, res, rpob):
    """Isoniazid-resistant isolates with a rifampicin-linked rpoB variant, by collection site."""
    tidy(ax, x=False, y=True)
    site = collections.defaultdict(list)
    for i, r in res["isoniazid"].items():
        if r and i in res["rifampicin"]:
            site[re.match(r"site\.([^.]+)", i).group(1)].append((i in rpob, res["rifampicin"][i]))
    stat = {s_: (sum(g for g, _ in v) / len(v), sum(q for _, q in v) / len(v), len(v)) for s_, v in site.items()}
    big = sorted((g, s_, n, ph) for s_, (g, ph, n) in stat.items() if n >= 50)
    allv = [x for v in site.values() for x in v]
    allp = sum(g for g, _ in allv) / len(allv)
    for x, (g, s_, n, ph) in enumerate(big):
        ax.bar(x, g, 0.64, color=DARK_RED, edgecolor=SURFACE, linewidth=0.5, zorder=3)
        ax.plot([x], [ph], marker="o", ms=4.6, mfc=SURFACE, mec=INK, mew=1.0, ls="none", zorder=5)
        ax.annotate(f"{g:.0%}", (x, max(g, ph)), xytext=(0, 5), textcoords="offset points", ha="center",
                    va="bottom", fontsize=6.2, color=INK2)
    ax.axhline(allp, color=AXIS, lw=0.9, ls=(0, (3, 3)), zorder=2)
    ax.annotate(f"all sites {allp:.0%}", xy=(-0.45, allp), xytext=(0, 3), textcoords="offset points",
                ha="left", va="bottom", fontsize=5.8, color=MUTED)
    ax.set_xticks(range(len(big)))
    ax.set_xticklabels([f"{s_}\nn={n}" for _, s_, n, _ in big], fontsize=6.0, color=INK, linespacing=1.2)
    ax.set_ylim(0, 1.12)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"], fontsize=6.8)
    ax.set_xlabel("collection site (isoniazid-resistant isolates)", fontsize=7.0)
    ax.set_ylabel("with a rifampicin-linked rpoB variant", fontsize=7.0)
    ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=DARK_RED, label="rpoB variant\n(TB-Profiler)"),
                       Line2D([], [], marker="o", ms=4.6, mfc=SURFACE, mec=INK, mew=1.0, ls="none",
                              label="rifampicin-resistant\nby MIC")],
              loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=6.4, frameon=False, handletextpad=0.6,
              borderaxespad=0.0, labelspacing=1.1)
    return big, allp


FIRST = ("isoniazid", "rifampicin", "ethambutol")
SECOND = ("ethionamide", "kanamycin", "amikacin", "levofloxacin", "moxifloxacin")
THIRD = ("bedaquiline", "delamanid")
TIER = (("first-line only", "#bcc5cf"), ("reached a second-line drug", "#667f96"),
        ("reached bedaquiline or delamanid", ORG_COLOUR["mtb"]))
MIN_TESTED = 10       # isolates tested for at least this many of the ten antimicrobials (--min-tested); 10 = all


def count_panel(ax, ph):
    """How many of the ten antimicrobials each isolate is resistant to, M. tuberculosis and P. aeruginosa.

    Isolates tested for at least MIN_TESTED of the ten antimicrobials. For M. tuberculosis each bar is split by the
    highest tier of antimicrobials the isolate is resistant to (first-line: isoniazid, rifampicin, ethambutol;
    second-line: ethionamide, aminoglycosides, fluoroquinolones; then bedaquiline and delamanid), which shows the
    staircase; P. aeruginosa has no such order and accumulates.
    """
    tidy(ax, x=False, y=True)
    stats = {}
    for k, off in (("mtb", -0.21), ("pseaer", 0.21)):
        drugs, mic, res = ph[k]
        ids = [i for i in {i for d in drugs for i in res[d]} if sum(i in res[d] for d in drugs) >= MIN_TESTED]
        n = len(ids)
        R = lambda i, ds: any(res[d].get(i, False) for d in ds)  # noqa: E731
        cnt = collections.Counter(sum(res[d].get(i, False) for d in drugs) for i in ids)
        for c in range(len(drugs) + 1):
            if k == "pseaer":
                ax.bar(c + off, 100 * cnt[c] / n, 0.40, color=ORG_COLOUR[k], edgecolor=SURFACE, linewidth=0.4, zorder=3)
                continue
            if c == 0:
                ax.bar(c + off, 100 * cnt[0] / n, 0.40, color="#d9dce3", edgecolor=SURFACE, linewidth=0.4, zorder=3)
                continue
            grp = [i for i in ids if sum(res[d].get(i, False) for d in drugs) == c]
            tiers = [sum(1 for i in grp if not R(i, SECOND + THIRD)),
                     sum(1 for i in grp if R(i, SECOND) and not R(i, THIRD)),
                     sum(1 for i in grp if R(i, THIRD))]
            bottom = 0.0
            for t, (lab, col) in zip(tiers, TIER):
                ax.bar(c + off, 100 * t / n, 0.40, bottom=bottom, color=col, edgecolor=SURFACE, linewidth=0.4,
                       zorder=3)
                bottom += 100 * t / n
        stats[k] = dict(n=n, cnt=cnt, ge8=100 * sum(v for c, v in cnt.items() if c >= 8) / n,
                        allr=100 * np.mean([all(res[d].get(i, False) for d in drugs if i in res[d]) for i in ids]))
        if k == "mtb":
            second_any = [i for i in ids if R(i, SECOND)]
            first_any = [i for i in ids if R(i, ("isoniazid", "rifampicin"))]
            both_s = [i for i in ids if not R(i, ("isoniazid", "rifampicin")) and
                      any(i in res[d] for d in ("isoniazid", "rifampicin"))]
            stats[k].update(new_n=sum(1 for i in ids if R(i, THIRD)), new_with_other=100 * np.mean(
                                [R(i, FIRST + SECOND) for i in ids if R(i, THIRD)]),
                            second_given_first=100 * np.mean([R(i, SECOND) for i in first_any]),
                            second_given_s=100 * np.mean([R(i, SECOND) for i in both_s]),
                            first_given_second=100 * np.mean([R(i, ("isoniazid", "rifampicin")) for i in second_any]),
                            second_no_first=100 * np.mean([R(i, SECOND) and not R(i, FIRST) for i in ids]))
    ax.set_xticks(range(11))
    ax.set_xticklabels([str(i) for i in range(11)], fontsize=6.8)
    ax.set_xlim(-0.6, 10.6)
    ax.set_ylim(0, 54)
    ax.set_yticks([0, 10, 20, 30, 40, 50])
    ax.set_yticklabels(["0", "10%", "20%", "30%", "40%", "50%"], fontsize=6.8)
    ax.set_xlabel("antimicrobials the isolate is resistant to (of ten)", fontsize=7.0)
    ax.set_ylabel("isolates", fontsize=7.0)
    ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    head = lambda t: Patch(facecolor="none", edgecolor="none", label=t)  # noqa: E731
    handles = [head("M. tuberculosis"), Patch(facecolor="#d9dce3", label="susceptible to all")] + \
              [Patch(facecolor=col, label=lab) for lab, col in TIER] + \
              [head("P. aeruginosa"), Patch(facecolor=ORG_COLOUR["pseaer"], label="every resistance pattern")]
    leg = ax.legend(handles=handles, loc="upper right", fontsize=5.8, frameon=False, handlelength=1.0,
                    borderaxespad=0.2, labelspacing=0.45)
    for t in leg.get_texts():
        if t.get_text() in ("M. tuberculosis", "P. aeruginosa"):
            t.set_fontproperties(ITALIC)
            t.set_fontsize(5.8)
    return stats


# column order of the profile matrix: by treatment tier (M. tuberculosis) or by class (P. aeruginosa)
PROFILE_COLS = {"mtb": [("first-line", ["isoniazid", "rifampicin", "ethambutol"]),
                        ("second-line", ["ethionamide", "levofloxacin", "moxifloxacin", "kanamycin", "amikacin"]),
                        ("new", ["bedaquiline", "delamanid"])],
                "pseaer": [("fluoroquin.", ["ciprofloxacin", "levofloxacin"]),
                           ("aminoglyc.", ["tobramycin", "amikacin"]),
                           ("carbapenems", ["imipenem", "meropenem"]),
                           ("other β-lactams", ["ceftazidime", "cefepime", "aztreonam", "piperacillin-tazobactam"])]}
NAME_GROUP = {"new": "new drugs"}


def profile_panel(axes, ph, sort="count", col_order="group"):
    """One row per isolate (tested for at least MIN_TESTED antimicrobials), one column per antimicrobial; filled where
    the isolate is resistant. Rows sorted by the number of resistances, then by profile, so an ordered accumulation
    shows as a staircase and an unordered one as scatter."""
    out = {}
    for ax, k in zip(axes, ("mtb", "pseaer")):
        drugs, mic, res = ph[k]
        groups = PROFILE_COLS[k]
        cols = [d for _, ds in groups for d in ds]
        ids = [i for i in {i for d in drugs for i in res[d]} if sum(i in res[d] for d in drugs) >= MIN_TESTED]
        if col_order == "prevalence":              # columns by the share of resistant isolates, most common first
            prev = {d: np.mean([res[d][i] for i in ids if i in res[d]]) for d in cols}
            cols = sorted(cols, key=lambda d: -prev[d])
            groups = [("", cols)]
        X = np.full((len(ids), len(cols)), np.nan)
        for r_, i in enumerate(ids):
            for c_, d in enumerate(cols):
                if i in res[d]:
                    X[r_, c_] = 1.0 if res[d][i] else 0.0
        def first_key(row):
            tested = int(np.sum(~np.isnan(row)))
            res_ = int(np.nansum(row))
            return round(res_ / tested, 6) if sort == "pct" else res_      # share of tested antimicrobials, or the count
        key = [(first_key(row), tuple(-np.nan_to_num(row, nan=0.0).astype(int))) for row in X]
        order = sorted(range(len(ids)), key=lambda r_: key[r_])
        X = X[order]
        cmap = plt.matplotlib.colors.ListedColormap(["#e9ebf0", ORG_COLOUR[k]])
        cmap.set_bad("#ffffff")
        ax.imshow(np.ma.masked_invalid(X), cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
        ax.set_yticks([])
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([CODE[d] for d in cols], fontsize=6.0, color=INK, rotation=90)
        ax.tick_params(length=0, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        # group brackets above the columns
        start = 0
        for lab, ds in (groups if col_order == "group" else []):
            ax.annotate("", xy=(start - 0.4, 1.012), xytext=(start + len(ds) - 0.6, 1.012),
                        xycoords=("data", "axes fraction"), textcoords=("data", "axes fraction"),
                        arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8, shrinkA=0, shrinkB=0), annotation_clip=False)
            ax.annotate(NAME_GROUP.get(lab, lab), xy=(start + (len(ds) - 1) / 2, 1.025), xycoords=("data", "axes fraction"),
                        ha="center", va="bottom", fontsize=5.8, color=INK2, annotation_clip=False)
            if start:
                ax.axvline(start - 0.5, color="#ffffff", lw=1.4)
            start += len(ds)
        # prevalence under the codes
        for c_, d in enumerate(cols):
            v = X[:, c_]
            ax.annotate(f"{100 * np.nanmean(v):.0f}%", xy=(c_, 0.0), xycoords=("data", "axes fraction"),
                        xytext=(0, -26), textcoords="offset points", ha="center", va="top", fontsize=5.4, color=MUTED)
        n = len(ids)
        none = int(sum(1 for r_ in range(n) if np.nansum(X[r_]) == 0))
        every = int(sum(1 for r_ in range(n) if np.nansum(X[r_]) > X.shape[1] / 2))      # resistant to more than half
        out[k] = dict(n=n, none=none, every=every)
        keys = [Patch(facecolor=ORG_COLOUR[k], label="resistant"), Patch(facecolor="#e9ebf0", label="susceptible")]
        if MIN_TESTED < 10:
            keys.append(Patch(facecolor="#ffffff", edgecolor="#b5b3a9", linewidth=0.6, label="not tested"))
        if SUBMISSION:                    # offsets in points: the shorter panels of the 170 mm figure
            from matplotlib.transforms import offset_copy
            ax.legend(handles=keys, loc="upper left", bbox_to_anchor=(0.0, 0.0),
                      bbox_transform=offset_copy(ax.transAxes, fig=ax.figure, x=0, y=-38, units="points"),
                      ncol=3, fontsize=6.0, frameon=False, handlelength=1.1, handletextpad=0.5, columnspacing=1.6,
                      borderaxespad=0.0)
        else:
            ax.legend(handles=keys, loc="upper left", bbox_to_anchor=(0.0, -0.215), ncol=3, fontsize=6.0, frameon=False,
                      handlelength=1.1, handletextpad=0.5, columnspacing=1.6, borderaxespad=0.0)
        # descriptive counts of each collection; the two collections were assembled differently, so they are not
        # compared between organisms
        ax.annotate(f"n = {n:,}" + ("" if SUBMISSION else " isolates") + f"  ·  susceptible to all: {none:,} ({100 * none / n:.0f}%)  ·  "
                    f"resistant to more than half: {every:,} ({100 * every / n:.0f}%)", xy=(0.0, 0.0), xycoords="axes fraction",
                    xytext=(0, -56 if SUBMISSION else -52), textcoords="offset points", ha="left", va="top", fontsize=5.8, color=MUTED)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--isolates", default=str(HOME / "AMRPredict_isolates_used.xlsx"))
    ap.add_argument("--tbprofiler-variants", default=str(HOME / "tbprofiler.variants.csv"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="figure5_draft")
    ap.add_argument("--width", type=float, default=7.8)
    ap.add_argument("--min-tested", type=int, default=10,
                    help="c, d: isolates tested for at least this many of the ten antimicrobials (10: all)")
    ap.add_argument("--col-order", choices=["group", "prevalence"], default="group",
                    help="columns of c, d: by treatment tier / drug class, or by the share of resistant isolates")
    ap.add_argument("--sort", choices=["count", "pct"], default="count",
                    help="rows of c, d: by the number of resistances, or by the share of the tested antimicrobials")
    ap.add_argument("--sites", action="store_true", help="also write the by-site figure (supplementary)")
    a = ap.parse_args()
    global MIN_TESTED
    MIN_TESTED = a.min_tested
    wb = load_workbook(a.isolates, read_only=True)
    ph = {k: phenotypes(wb, k) for k in ("mtb", "pseaer")}
    M = {k: matrix(*ph[k]) for k in ph}

    fig = plt.figure(figsize=(a.width, fig_h(9.1)))
    gs = fig.add_gridspec(2, 2, hspace=0.50 if SUBMISSION else 0.58, wspace=0.62,
                          height_ratios=[0.86, 0.60] if SUBMISSION else [0.80, 0.78])
    ax_a = fig.add_subplot(gs[0, 0])
    im = heat(ax_a, ph["mtb"][0], M["mtb"], CLASSES["mtb"])
    ax_b = fig.add_subplot(gs[0, 1])
    im_b = heat(ax_b, ph["pseaer"][0], M["pseaer"], CLASSES["pseaer"], cmap=SEQ_PA)
    for ax in (ax_a, ax_b):
        ax.set_anchor("N")                 # square cells: keep the heatmap at the top of its row
    gs2 = gs[1, :].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.30)
    ax_c = fig.add_subplot(gs2[0, 0])
    ax_c2 = fig.add_subplot(gs2[0, 1])
    profile_panel((ax_c, ax_c2), ph, sort=a.sort, col_order=a.col_order)
    _f, _ax = plt.subplots()
    cstats = count_panel(_ax, ph)            # the numbers quoted in the text
    plt.close(_f)
    top = header(fig, "Figure 5 (draft) — resistance to one antimicrobial is largely predictable from the MIC of another",
                 "a, b: area under the ROC curve for telling isolates resistant to the antimicrobial in the row "
                 "from susceptible ones by the log₂ MIC of the antimicrobial\nin the column (isolates with both MICs; at "
                 "least ten resistant). Outlined: pairs within one drug class, which share a mechanism.\nc, d: the same pattern isolate by isolate: "
                 "one row per isolate tested for all ten antimicrobials, sorted by the number of\n"
                 "resistances and then by profile; columns grouped by treatment tier (M. tuberculosis) or drug class "
                 "(P. aeruginosa);\npercentages under the codes: share of the isolates resistant to the antimicrobial; counts below "
                 "describe each collection and are not compared between organisms.", pad_in=0.30)
    fig.subplots_adjust(top=top - 0.045, bottom=0.165 if SUBMISSION else 0.105, left=0.17, right=0.965)
    for ax in (ax_a, ax_b):
        ax.apply_aspect()                  # final positions, after the square-cell adjustment
    pa_, pb_ = ax_a.get_position(), ax_b.get_position()
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv_ = fig.transFigure.inverted()
    lab_x0 = {}
    for axx, hm, ref in ((ax_c, ax_a, pa_), (ax_c2, ax_b, pb_)):  # c under a, d under b: the width of the heatmap with its names
        left = min(t.get_window_extent(rend).x0 for t in hm.get_yticklabels())
        x0 = inv_.transform((left, 0))[0]
        lab_x0[axx] = x0
        q = axx.get_position()
        axx.set_position([x0, q.y0, ref.x1 - x0, q.height])
    for img, ref in ((im, pa_), (im_b, pb_)):          # one colour bar under each heatmap, same scale
        w_ = 0.50 * ref.width
        cax = fig.add_axes([ref.x0 + (ref.width - w_) / 2, min(pa_.y0, pb_.y0) - 0.060, w_, 0.009])
        cb = fig.colorbar(img, cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=6.2, length=0)
        cb.outline.set_visible(False)
        cb.set_ticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        cb.set_ticklabels(["50", "60", "70", "80", "90", "100"])
        cb.set_label("AUROC (%; 50 = no information)", fontsize=6.4)
    for ax, letter, dx, dy in ((ax_a, "a", -0.40, 1.07), (ax_b, "b", -0.40, 1.07)):
        ax.annotate(letter, xy=(dx, dy), xycoords="axes fraction", fontsize=10, fontweight="bold", color=INK,
                    ha="right", va="bottom", annotation_clip=False)
    for axx, ref, letter in ((ax_c, pa_, "c"), (ax_c2, pb_, "d")):   # in the letter column of a and b
        fig.text(ref.x0 - 0.40 * ref.width, axx.get_position().y1 + 0.030, letter, fontsize=10, fontweight="bold",
                 color=INK, ha="right", va="bottom")
    # one header per organism over its column, above the panel letters (no rule)
    for ax, axc, name in ((ax_a, ax_c, "M. tuberculosis"), (ax_b, ax_c2, "P. aeruginosa")):
        pos = ax.get_position()
        x0, x1 = lab_x0[axc], pos.x1
        fig.text((x0 + x1) / 2, pos.y1 + 0.046, name, ha="center", va="bottom", fontsize=10.6,
                 fontproperties=ITALIC_BOLD, color=INK)

    out = pathlib.Path(a.out_dir)
    save(fig, out, (a.stem if a.sort == "count" else a.stem.replace("_draft", "_pctsort_draft")).replace("_draft", "_draft" if a.col_order == "group" else "_prevorder_draft").replace(
        "_draft", "_draft" if a.min_tested == 10 else "_min%d_draft" % a.min_tested))
    with open(out / "figure5_values.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["organism", "resistance_to", "mic_of", "auroc"])
        for k, org in (("mtb", "M. tuberculosis"), ("pseaer", "P. aeruginosa")):
            drugs = ph[k][0]
            for i, da in enumerate(drugs):
                for j, db in enumerate(drugs):
                    if i != j and not np.isnan(M[k][i, j]):
                        w.writerow([org, da, db, round(float(M[k][i, j]), 4)])
    # numbers for the text
    for k, org in (("mtb", "M. tuberculosis"), ("pseaer", "P. aeruginosa")):
        drugs = ph[k][0]
        inside = {(drugs.index(x), drugs.index(y)) for cls in CLASSES[k] for x in cls for y in cls if x != y}
        off = [(M[k][i, j], drugs[i], drugs[j]) for i in range(len(drugs)) for j in range(len(drugs))
               if i != j and (i, j) not in inside and not np.isnan(M[k][i, j])]
        ins = [M[k][i, j] for i, j in inside if not np.isnan(M[k][i, j])]
        off.sort(reverse=True)
        print(f"[INFO] {org}: within-class AUROC median {np.median(ins):.2f} (n={len(ins)}); between-class median "
              f"{np.median([o[0] for o in off]):.2f} (n={len(off)}), {sum(o[0] >= 0.8 for o in off)} cells >= 0.80")
        print("       highest between-class: " + "; ".join(f"{CODE[x]}<-{CODE[y]} {v:.2f}" for v, x, y in off[:5]))
    m_, p_ = cstats["mtb"], cstats["pseaer"]
    print(f"[INFO] M. tuberculosis (n={m_['n']}): resistant to >=8 {m_['ge8']:.1f}%, to every tested {m_['allr']:.1f}%; "
          f"second-line given INH or RIF resistance {m_['second_given_first']:.0f}%, given susceptible to both "
          f"{m_['second_given_s']:.1f}%; INH or RIF resistant among second-line resistant {m_['first_given_second']:.0f}%; "
          f"second-line without first-line {m_['second_no_first']:.1f}%; bedaquiline or delamanid resistant {m_['new_n']} "
          f"({m_['new_with_other']:.0f}% also resistant to a first- or second-line drug)")
    print(f"[INFO] P. aeruginosa (n={p_['n']}): resistant to >=8 {p_['ge8']:.1f}%, to every tested {p_['allr']:.1f}%; "
          f"susceptible to all {100 * p_['cnt'][0] / p_['n']:.0f}%")
    if a.sites:                                   # the by-site view, as a supplementary figure
        fig2, ax2 = plt.subplots(figsize=(6.2, 3.6))
        big, allp = sites(ax2, ph["mtb"][1], ph["mtb"][2], rpob_carriers(a.tbprofiler_variants))
        ax2.set_title("Rifampicin resistance among isoniazid-resistant isolates, by site", loc="left", fontsize=7.6,
                      color=INK, fontweight="semibold", pad=13)
        fig2.subplots_adjust(left=0.13, right=0.80, bottom=0.2, top=0.8)
        save(fig2, out, "supplementary_figure5_sites_draft")
        print("[INFO] sites: " + "; ".join(f"{s} rpoB {p:.0%}, MIC {q:.0%} (n={n})" for p, s, n, q in big)
              + f"; all sites {allp:.0%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
