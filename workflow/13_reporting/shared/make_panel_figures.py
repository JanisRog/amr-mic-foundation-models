#!/usr/bin/env python3
"""Publication figures for the MTB panel, straight from the results workbook.

One source of truth: `AMRPredict_panel_results.xlsx`, the same file section 14's
numbers were read out of. Nothing is retyped here -- a figure that disagrees with
the protocol should mean the workbook changed, not that a number was copied wrong
into a plotting script.

  1   the panel      cost and S/R F1 side by side, holdout, all three models
                     against TB-Profiler. Side by side because they disagree: the
                     model shifts toward specificity on almost every drug, which
                     F1 rewards and a 2:1 very-major-error weighting punishes.
                     Leading with either alone misrepresents the panel.
  2   the trade-off  MER against vMER with iso-cost diagonals, so that
                     disagreement is visible as a direction of movement rather
                     than asserted in a caption.
  3   the ladder     ordinal MIC accuracy, which the comparator has no equivalent
                     of, over the per-class breakdown that says where it fails.
  4   the two new    bedaquiline and delamanid, on their own axes.
      drugs
  S1  stability      development CV against holdout, as the overfitting check.

**Bedaquiline and delamanid are not in figure 1.** They have 12 and 15 resistant
isolates in the holdout, where one isolate moves bedaquiline's cost by 0.167 --
larger than the margin that decided every other drug -- and their costs are 4x
the rest of the panel's, so including them would both imply a precision they do
not have and compress the eight real results into the left fifth of the axis.
They get figure 4, at their own scale, where the finding is the catalogue's
miscalibration rather than any model's margin.

Colors are the validated categorical slots 1-3 for the three models. TB-Profiler
is deliberately *not* a fourth hue: it is the incumbent everything is read
against, so it takes neutral ink and an open marker. The shape, not the color,
distinguishes it, which also survives grayscale printing.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

MODELS = ["evo2", "carbon", "dnabert"]
MODEL_LABEL = {"evo2": "Evo 2 7B", "carbon": "Carbon 8B", "dnabert": "DNABERT-2"}
# the three foundation models in shades of purple: purple means embeddings (Figures 2, 4); lightness tells them apart
SERIES = {"evo2": "#4b2a7a", "carbon": "#8a63c2", "dnabert": "#c7b3e5"}
# the two organisms: navy and coral-pink (chosen by the author 2026-10-20, with the purple of the embeddings as one family);
# separation 35 normal, 31 deutan, 38 tritan, 19 protan -- for red-blind readers shape (circle/diamond) and lightness
# (36 / 61) carry the difference as well
ORG_COLOUR = {"mtb": "#1B3F5E", "pseaer": "#D9577F"}
ORG_PALE = {"mtb": "#bcc5cf", "pseaer": "#f0c3d1"}
BASELINE = "TB-Profiler 6.3.0"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#ffffff"
PALE = "#9ec5f4"
SEQ = LinearSegmentedColormap.from_list(
    "blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])

UNDERPOWERED = ["bedaquiline", "delamanid"]

def _italic_font(bold: bool):
    """A real italic face for organism names.

    The body font (Helvetica Neue) lives in a .ttc collection of which matplotlib loads only the
    first face, so fontstyle="italic" silently renders upright. Use an italic file directly:
    Arial where present, matplotlib's bundled DejaVu Sans Oblique otherwise.
    """
    from matplotlib.font_manager import FontProperties
    name = "Arial Bold Italic.ttf" if bold else "Arial Italic.ttf"
    for d in ("/System/Library/Fonts/Supplemental", "/Library/Fonts", "/usr/share/fonts/truetype/msttcorefonts",
              str(pathlib.Path.home() / "Library/Fonts")):
        if (pathlib.Path(d) / name).exists():
            return FontProperties(fname=str(pathlib.Path(d) / name))
    return FontProperties(fname=str(pathlib.Path(mpl.get_data_path()) / "fonts/ttf"
                                    / ("DejaVuSans-BoldOblique.ttf" if bold else "DejaVuSans-Oblique.ttf")))


# Colours of the two kinds of input, the same in every figure that separates them (Figures 2 and 4):
# embeddings of the target loci in purple, the non-embedding input (catalogue call, variant columns) in grey.
# The reported model, embeddings together with the catalogue call, stays blue (Figures 1 and 3).
PURPLE = "#6b3fa0"
GREY = "#8a8880"

ITALIC = _italic_font(False)          # use as fontproperties=ITALIC for species names
ITALIC_BOLD = _italic_font(True)


def _bold_font():
    """Upright bold from the same Arial family as the italics: the first sans-serif face (Helvetica Neue) has no
    registered bold here, so fontweight="bold" renders regular."""
    from matplotlib.font_manager import FontProperties
    for d in ("/System/Library/Fonts/Supplemental", "/Library/Fonts", "/usr/share/fonts/truetype/msttcorefonts",
              str(pathlib.Path.home() / "Library/Fonts")):
        if (pathlib.Path(d) / "Arial Bold.ttf").exists():
            return FontProperties(fname=str(pathlib.Path(d) / "Arial Bold.ttf"))
    return FontProperties(family="sans-serif", weight="bold")


BOLD = _bold_font()

mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 400,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "text.color": INK, "axes.labelcolor": INK2,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
if os.environ.get("FIG_SUBMISSION") == "1":      # one family for regular, bold and italic (the bold and italic faces are Arial)
    mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]   # fallback per glyph: Arial has no subscript two


# ── data ──────────────────────────────────────────────────────────────────────
def load(xlsx: pathlib.Path):
    import openpyxl
    wb = openpyxl.load_workbook(xlsx, data_only=True)

    def sheet(name):
        rows = list(wb[name].iter_rows(values_only=True))
        head = [str(c) for c in rows[0]]
        return [dict(zip(head, r)) for r in rows[1:] if r[0] is not None]

    return sheet("summary"), sheet("per_class")


def pick(summary, drug, model, split, field):
    for r in summary:
        if r["drug"] == drug and r["model"] == model and r["split"] == split:
            v = r[field]
            return None if v in ("", None) else float(v)
    return None


def pickp(*a, **k):
    """pick() in per cent, for the figures that print rates and accuracies."""
    v = pick(*a, **k)
    return None if v is None else 100 * v


# ── chrome ────────────────────────────────────────────────────────────────────
def tidy(ax, *, x=True, y=False):
    """Recessive chrome: gridlines behind the marks, no box around the data."""
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    if x:
        ax.xaxis.grid(True, color=GRID, lw=0.6, zorder=0)
    if y:
        ax.yaxis.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)


# Final figures: FIG_FINAL=1 draws a one-line title ("Fig. N: title") in place of the draft title and the description
# under it; FIG_TITLE overrides the title, for figures whose title the script cannot know (the per-organism sheets).
FINAL = os.environ.get("FIG_FINAL") == "1"
# Submission figures (FIG_SUBMISSION=1): no title, 180 mm wide, all text clamped to 5-7 pt and panel letters 8 pt bold,
# as the Nature figure specifications ask; the title and the description go in the legend.
SUBMISSION = os.environ.get("FIG_SUBMISSION") == "1"
SUB_WIDTH = 180 / 25.4
REF_WIDTH, TITLE_PT = 7.8, 11.0     # the narrow figures; wider ones scale title and panel letters by width / 7.8,
                                    # so that at one print width the title and letters are the same size in all figures


def title_pt(fig):
    """11 pt in the final files: figures built at 9.4 in are scaled to 7.8 in afterwards, so they are built at 11 x 9.4/7.8;
    the supplementary sheets (wider) are not scaled and keep 11 pt."""
    w = fig.get_figwidth()
    return TITLE_PT * (w / REF_WIDTH if w < 10 else 1.0)


def fig_h(h, cap=6.55):
    """Submission figures: at most 170 mm (the full page depth in the Nature specifications)."""
    return min(h, cap) if SUBMISSION else h


def fit_fonts(fig, lo=5.0, hi=7.0, letter=8.0):
    """Submission figures: text between lo and hi pt, single-letter panel labels letter pt bold upright."""
    if not SUBMISSION:
        return
    arts = list(fig.texts)
    for ax in fig.axes:
        arts += list(ax.texts) + [ax.title, ax._left_title, ax._right_title, ax.xaxis.label, ax.yaxis.label]
        arts += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        leg = ax.get_legend()
        if leg is not None:
            arts += list(leg.get_texts()) + ([leg.get_title()] if leg.get_title() else [])
    for leg in fig.legends:
        arts += list(leg.get_texts()) + ([leg.get_title()] if leg.get_title() else [])
    for t in arts:
        s_ = t.get_text()
        if not s_.strip():
            continue
        if "\u2082" in s_:                 # Arial has no subscript two: a mathtext subscript instead
            t.set_text(s_.replace("\u2082", "$_2$"))
            s_ = t.get_text()
        if len(s_) == 1 and s_.isalpha() and t.get_fontsize() >= 9.5:
            t.set_fontproperties(BOLD)
            t.set_fontsize(letter)
        else:
            t.set_fontsize(min(hi, max(lo, t.get_fontsize())))


def scale_letters(fig):
    """Final figures: the single-letter panel labels grow with the figure width, like the title."""
    k = fig.get_figwidth() / REF_WIDTH
    if not FINAL or abs(k - 1) < 0.01:
        return
    for t in list(fig.texts) + [t for ax in fig.axes for t in ax.texts]:
        if len(t.get_text()) == 1 and t.get_text().isalpha() and abs(t.get_fontsize() - 10) < 0.01:
            t.set_fontsize(10 * k)
FINAL_TITLES = [
    ("Figure 1", "Fig. 1: From genotype to phenotype: A framework for MIC prediction from\nDNA foundation-model embeddings"),
    ("Figure 2", "Fig. 2: Embeddings against the variants from the same loci"),
    ("Figure 3", "Fig. 3: Classification performance evaluation of the framework against established tools"),
    ("Figure 4", "Fig. 4: Attribution analysis of the framework's predictions"),
    ("Figure 5", "Fig. 5: Co-resistance analysis"),
    ("Extended Data Figure 1", "Extended Data Fig. 1: The three foundation models, embeddings alone"),
    ("Extended Data Figure 2", "Extended Data Fig. 2: Essential agreement of the reported model by MIC class"),
    ("Extended Data Figure 3", "Extended Data Fig. 3: How well one antimicrobial predicts another differs between centres"),
    ("Supplementary Figure 1", "Supplementary Fig. 1: Predicted against measured MIC"),
    ("Supplementary Figure 2", "Supplementary Fig. 2: Predicted against measured MIC"),
]


def final_title(title):
    """FIG_TITLES: JSON {start of the draft title: final title} for the sheets of one run (organism-specific titles);
    FIG_TITLE overrides every title of the run."""
    import json
    for k, v in json.loads(os.environ.get("FIG_TITLES", "{}")).items():
        if title.startswith(k):
            return v
    return os.environ.get("FIG_TITLE") or next((v for k, v in FINAL_TITLES if title.startswith(k)), title)


def header(fig, title, subtitle, pad_in=0.0):
    """Title block placed in inches from the top, so it cannot collide.

    The first version positioned both in figure fractions, which is a different
    distance on a 3.8in figure than on a 5.6in one -- the title and the sentence
    under it overlapped on the two short figures and not on the tall ones.
    Returns the fraction to hand `top=` so the axes start below the block.

    `pad_in` reserves the height of an axes title on top of that: `top=` places
    the axes *frame*, and a panel title is drawn above the frame, so figures with
    per-panel titles need the extra room or the title lands back in the subtitle.
    """
    h = fig.get_figheight()
    if SUBMISSION:
        return 1 - (0.10 + pad_in) / h
    if FINAL:
        fig.text(0.008, 1 - 0.26 / h, final_title(title), fontsize=title_pt(fig),
                 fontproperties=BOLD, wrap=True,
                 color=INK, va="top", ha="left")
        return 1 - (0.50 + pad_in) / h
    n = subtitle.count("\n") + 1
    fig.text(0.008, 1 - 0.26 / h, title, fontsize=11, fontweight="semibold",
             color=INK, va="top", ha="left")
    fig.text(0.008, 1 - 0.50 / h, subtitle, fontsize=7.5, color=INK2,
             va="top", ha="left", linespacing=1.55)
    return 1 - (0.62 + 0.145 * n + pad_in) / h


def legend_models(fig, *, baseline=True, y=-0.005):
    h = [Line2D([], [], marker="o", ls="none", ms=5, mfc=SERIES[m], mec=SERIES[m],
                label=MODEL_LABEL[m]) for m in MODELS]
    if baseline:
        h.append(Line2D([], [], marker="o", ls="none", ms=5.5, mfc="none",
                        mec=INK2, mew=1.2, label=BASELINE))
    fig.legend(handles=h, ncol=len(h), loc="lower center", fontsize=7.5,
               handletextpad=0.4, columnspacing=1.6, bbox_to_anchor=(0.5, y))


def save(fig, out: pathlib.Path, stem: str):
    out.mkdir(parents=True, exist_ok=True)
    scale_letters(fig)
    fit_fonts(fig)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", format=ext)
    print(f"[OK] {out/stem}.pdf / .png")
    plt.close(fig)


# ── figure 1 ──────────────────────────────────────────────────────────────────
def fig_panel(summary, out: pathlib.Path):
    """vMER, MER and F1 against the catalogue, one row per drug.

    The two error rates are shown separately rather than folded into the cost,
    because a single weighted number hides which of the two moved. For the
    reported model it is the major error rate that falls (fewer susceptible
    isolates called resistant); the very major error rate is level or higher.
    """
    drugs = [d for d in sorted({r["drug"] for r in summary}) if d not in UNDERPOWERED]
    margin = {d: pick(summary, d, BASELINE, "holdout", "cost")
                 - pick(summary, d, "evo2", "holdout", "cost") for d in drugs}
    order = sorted(drugs, key=lambda d: margin[d], reverse=True)   # as before: by cost

    def same_calls(d):
        return all(abs(pick(summary, d, "evo2", "holdout", f)
                       - pick(summary, d, BASELINE, "holdout", f)) < 1e-4
                   for f in ("vMER", "MER"))

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.7), sharey=True,
                             gridspec_kw={"wspace": 0.44})
    panels = ((axes[0], "vMER", "Very major error rate (%)",
               "resistant called susceptible\n(lower is better)", -1),
              (axes[1], "MER", "Major error rate (%)",
               "susceptible called resistant\n(lower is better)", -1),
              (axes[2], "SR_f1", "F1, resistant class (%)", "\n(higher is better)", +1))
    for ax, field, title, xlab, sign in panels:
        tidy(ax)
        for i, d in enumerate(order):
            b = pickp(summary, d, BASELINE, "holdout", field)
            vals = [pickp(summary, d, m, "holdout", field) for m in MODELS]
            ax.plot([min(vals + [b]), max(vals + [b])], [i, i], color=GRID, lw=2.4,
                    solid_capstyle="round", zorder=1)
            # The reported model is drawn last so that it is never hidden under the
            # other two where the three coincide.
            for m, v in reversed(list(zip(MODELS, vals))):
                ax.plot(v, i, "o", ms=4.6, color=SERIES[m], mec=SURFACE, mew=0.8,
                        zorder=3)
            ax.plot(b, i, "o", ms=5.6, mfc="none", mec=INK2, mew=1.3, zorder=4)
        ax.set_yticks(range(len(order)))
        # A drug whose model makes the catalogue's calls exactly is not a near-tie;
        # the marks coincide in every panel, and the row label says why.
        ax.set_yticklabels([d.capitalize() + (" *" if same_calls(d) else "")
                            for d in order], fontsize=8, color=INK)
        ax.set_title(title, loc="left", color=INK, pad=7, fontweight="medium",
                     fontsize=8.6)
        ax.set_xlabel(xlab, fontsize=7, linespacing=1.35)
        ax.tick_params(axis="x", labelsize=7)
        ax.margins(x=0.10)
        # one selective number per row: Evo 2 against the catalogue
        ax.annotate("Δ Evo 2", xy=(1.03, 1.0), xycoords="axes fraction",
                    fontsize=6.2, color=MUTED, ha="left", va="bottom")
        for i, d in enumerate(order):
            delta = sign * (pickp(summary, d, "evo2", "holdout", field)
                            - pickp(summary, d, BASELINE, "holdout", field))
            ax.annotate(f"{delta:+.1f}" if abs(delta) > 1e-2 else "0.0",
                        xy=(1.03, i), xycoords=("axes fraction", "data"),
                        fontsize=6.6, color=INK2 if delta > 1e-2 else MUTED,
                        ha="left", va="center", annotation_clip=False)
    axes[0].set_ylim(-0.7, len(order) - 0.3)
    axes[0].invert_yaxis()          # once -- the axes are shared

    # Counted from the rows drawn, never typed.
    def better(field, sign):
        return sum(1 for d in drugs
                   if sign * (pick(summary, d, "evo2", "holdout", field)
                              - pick(summary, d, BASELINE, "holdout", field)) > 1e-4)
    n_v, n_m, n_f = better("vMER", -1), better("MER", -1), better("SR_f1", +1)
    n = len(drugs)
    star = ("  * same calls as TB-Profiler." if any(same_calls(d) for d in drugs) else "")
    top = header(fig, "Sequence embeddings against the catalogue",
                 f"Held-out set, one evaluation per drug. Evo 2 has the lower very major "
                 f"error rate on {n_v} of {n}, the lower major error rate on {n_m} of {n}\n"
                 f"and the higher F1 on {n_f} of {n}. Δ: Evo 2 against TB-Profiler, "
                 f"positive = better.{star} Bedaquiline and delamanid are in figure 4.",
                 pad_in=0.28)
    legend_models(fig, y=-0.01)
    fig.subplots_adjust(top=top, bottom=0.25, left=0.125, right=0.935)
    save(fig, out, "fig1_panel")


# ── figure 2 ──────────────────────────────────────────────────────────────────
NUDGE = {          # label offsets in points, hand-set to stop collisions
    "amikacin": (8, -9, "left"), "ethambutol": (0, 9, "center"),
    "ethionamide": (0, 9, "center"), "isoniazid": (8, -3, "left"),
    "kanamycin": (8, 3, "left"), "levofloxacin": (8, -8, "left"),
    "moxifloxacin": (8, -3, "left"), "rifampicin": (-8, -6, "right"),
}


def fig_tradeoff(summary, out: pathlib.Path):
    drugs = [d for d in sorted({r["drug"] for r in summary}) if d not in UNDERPOWERED]
    fig, ax = plt.subplots(figsize=(5.8, 5.0))
    tidy(ax, x=True, y=True)

    xmax = max(pickp(summary, d, m, "holdout", "MER")
               for d in drugs for m in (BASELINE, "evo2")) * 1.22
    ymax = max(pickp(summary, d, m, "holdout", "vMER")
               for d in drugs for m in (BASELINE, "evo2")) * 1.18

    for c in np.arange(5, 75, 5):                # 2·vMER + MER = c (per cent)
        x_end = min(xmax, c)
        ax.plot([0, x_end], [c / 2, (c - x_end) / 2], color=GRID, lw=0.6, zorder=0)
    ax.annotate("iso-cost lines, 5 points apart", xy=(0.99, 0.015),
                xycoords="axes fraction", fontsize=6.6, color=MUTED, ha="right")
    # the glyph, not the phrasing: Helvetica Neue has no U+2199, and a missing
    # glyph renders as a blank box in the PDF without failing the run
    ax.annotate("toward the origin = better on both axes", xy=(0.012, 0.015),
                xycoords="axes fraction", fontsize=6.6, color=MUTED, ha="left")

    for d in drugs:
        x0 = pickp(summary, d, BASELINE, "holdout", "MER")
        y0 = pickp(summary, d, BASELINE, "holdout", "vMER")
        x1 = pickp(summary, d, "evo2", "holdout", "MER")
        y1 = pickp(summary, d, "evo2", "holdout", "vMER")
        moved = abs(x1 - x0) > 1e-2 or abs(y1 - y0) > 1e-2
        if moved:
            ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                        arrowprops=dict(
                            arrowstyle="-|>,head_width=0.13,head_length=0.30",
                            color=INK2, lw=0.9, shrinkA=4, shrinkB=4))
        ax.plot(x0, y0, "o", ms=6, mfc="none", mec=INK2, mew=1.3, zorder=3)
        ax.plot(x1, y1, "o", ms=5, color=SERIES["evo2"], mec=SURFACE, mew=0.8,
                zorder=4)
        dx, dy, ha = NUDGE[d]
        ax.annotate(d.capitalize() + ("" if moved else "  (identical)"),
                    xy=(x1, y1), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, fontsize=7, color=INK)

    ax.set_xlim(0, xmax); ax.set_ylim(0, ymax)
    ax.set_xlabel("Major error rate (%)   (1 − specificity)")
    ax.set_ylabel("Very major error rate (%)   (1 − sensitivity)")
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="none", ms=5.5, mfc="none", mec=INK2,
               mew=1.2, label="TB-Profiler 6.3.0"),
        Line2D([], [], marker="o", ls="none", ms=5, color=SERIES["evo2"],
               label="Evo 2 7B")],
        loc="upper right", fontsize=7.5, handletextpad=0.4)

    top = header(fig, "The panel trades sensitivity for specificity",
                 "Arrows run from TB-Profiler to Evo 2, held-out set. Cost weights "
                 "the vertical axis twice, so a drug can\ncross iso-cost lines in "
                 "the wrong direction while gaining F1 — ethionamide and kanamycin do.")
    fig.subplots_adjust(top=top, bottom=0.095, left=0.115, right=0.985)
    save(fig, out, "fig2_tradeoff")


# ── figure 3 ──────────────────────────────────────────────────────────────────
def fig_ladder(summary, per_class, out: pathlib.Path):
    drugs = sorted({r["drug"] for r in summary})
    order = sorted(drugs, key=lambda d: pick(summary, d, "evo2", "holdout", "OT"))

    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.5], hspace=0.40)

    ax = fig.add_subplot(gs[0])
    tidy(ax)
    for i, d in enumerate(order):
        ot = pickp(summary, d, "evo2", "holdout", "OT")
        ex = pickp(summary, d, "evo2", "holdout", "exact")
        ax.plot([ex, ot], [i, i], color=GRID, lw=2.4, solid_capstyle="round", zorder=1)
        ax.plot(ex, i, "o", ms=5, color=PALE, mec=SURFACE, mew=0.8, zorder=3)
        ax.plot(ot, i, "o", ms=5, color=SERIES["evo2"], mec=SURFACE, mew=0.8, zorder=3)
        if d in (order[0], order[-1]):
            ax.annotate(f"{ex:.1f}", (ex, i), xytext=(-6, 0), ha="right",
                        textcoords="offset points", fontsize=7, color=INK2,
                        va="center")
            ax.annotate(f"{ot:.1f}", (ot, i), xytext=(6, 0), ha="left",
                        textcoords="offset points", fontsize=7, color=INK2,
                        va="center")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([d.capitalize() for d in order], fontsize=8, color=INK)
    ax.set_xlim(28, 106)
    ax.set_xlabel("held-out isolates (%)")
    ax.set_title("MIC placement, Evo 2 — the comparator emits one bit and has no "
                 "entry here", loc="left", color=INK, pad=7, fontweight="medium")
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="none", ms=5, color=PALE, label="exact class"),
        Line2D([], [], marker="o", ls="none", ms=5, color=SERIES["evo2"],
               label="within one dilution")],
        loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=7.5,
        handletextpad=0.4, columnspacing=1.4)

    ax2 = fig.add_subplot(gs[1])
    cells, ns = {}, {}
    for r in per_class:
        if r["model"] == "evo2" and r["split"] == "holdout":
            k = (r["drug"], int(float(r["log2_mic"])))
            cells[k] = 100 * float(r["OT"]); ns[k] = float(r["n"])
    ks = list(range(min(k for _, k in cells), max(k for _, k in cells) + 1))
    # imshow draws row 0 at the top; the dumbbells above put the best drug at the
    # top, so the rows are reversed to keep one reading order across the figure.
    rows = list(reversed(order))
    grid = np.full((len(rows), len(ks)), np.nan)
    for i, d in enumerate(rows):
        for j, k in enumerate(ks):
            if (d, k) in cells:
                grid[i, j] = cells[(d, k)]
    im = ax2.imshow(grid, cmap=SEQ, vmin=0, vmax=100, aspect="auto")
    for i, d in enumerate(rows):
        for j, k in enumerate(ks):
            n = ns.get((d, k))
            if n is not None and n < 30:
                ax2.plot(j, i, marker="x", ms=3.4, mew=0.9, color="#d03b3b")
    ax2.set_xticks(range(len(ks)))
    ax2.set_xticklabels([str(k) for k in ks], fontsize=7)
    ax2.set_yticks(range(len(rows)))
    ax2.set_yticklabels([d.capitalize() for d in rows], fontsize=8, color=INK)
    ax2.set_xlabel("log₂ MIC class")
    for s in ax2.spines.values():
        s.set_visible(False)
    ax2.tick_params(length=0)
    ax2.set_title("…by MIC class: the headline is carried by the dense classes",
                  loc="left", color=INK, pad=7, fontweight="medium")
    cb = fig.colorbar(im, ax=ax2, fraction=0.020, pad=0.012, shrink=0.88)
    cb.set_label("within one dilution (%)", fontsize=7, color=INK2)
    cb.ax.tick_params(labelsize=7, length=0)
    cb.outline.set_visible(False)
    ax2.annotate("×  class holds fewer than 30 isolates; blank = no isolates",
                 xy=(0, 0), xycoords="axes fraction", xytext=(0, -32),
                 textcoords="offset points", fontsize=6.8, color=MUTED)

    top = header(fig, "The ordinal ladder is the contribution no baseline competes with",
                 "Evo 2, held-out set. Where a drug's MIC distribution is bimodal "
                 "the ladder holds within each mode and\nfails in the gap between "
                 "them — rifampicin's −2 and −1 classes are the clearest case.",
                 pad_in=0.28)
    fig.subplots_adjust(top=top, bottom=0.115, left=0.145, right=0.955)
    save(fig, out, "fig3_ladder")


# ── figure 4 ──────────────────────────────────────────────────────────────────
def fig_newdrugs(summary, out: pathlib.Path):
    """Bedaquiline and delamanid, at their own scale, reported as what they are."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2),
                             gridspec_kw={"wspace": 0.34})
    for ax, field, title, xlab in (
            (axes[0], "vMER", "Very major error rate",
             "%   1 − sensitivity   (lower is better)"),
            (axes[1], "SR_f1", "Susceptible/resistant F1", "%   (higher is better)")):
        tidy(ax)
        for i, d in enumerate(UNDERPOWERED):
            b = pickp(summary, d, BASELINE, "holdout", field)
            vals = [pickp(summary, d, m, "holdout", field) for m in MODELS]
            ax.plot([min(vals + [b]), max(vals + [b])], [i, i], color=GRID, lw=2.4,
                    solid_capstyle="round", zorder=1)
            for k, (m, v) in enumerate(zip(MODELS, vals)):
                ax.plot(v, i + (k - 1) * 0.17, "o", ms=5, color=SERIES[m],
                        mec=SURFACE, mew=0.8, zorder=3)
            ax.plot(b, i, "o", ms=6, mfc="none", mec=INK2, mew=1.3, zorder=4)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["Bedaquiline\n12 resistant", "Delamanid\n15 resistant"],
                           fontsize=8, color=INK)
        ax.set_ylim(-0.75, 1.75); ax.invert_yaxis()
        ax.set_title(title, loc="left", color=INK, pad=7, fontweight="medium")
        ax.set_xlabel(xlab)

    top = header(fig, "The two newest drugs carry no verdict",
                 "Held-out set. Bedaquiline's entire model-versus-catalogue "
                 "difference is one isolate — larger than the margin that\ndecided "
                 "every drug in figure 1. Delamanid's classifiers miss the same "
                 "eleven of fifteen. The models are offset\nvertically only so that "
                 "identical values stay visible; each row shares one x scale.",
                 pad_in=0.28)
    legend_models(fig, y=0.005)
    fig.subplots_adjust(top=top, bottom=0.27, left=0.125, right=0.985)
    save(fig, out, "fig4_newdrugs")


# ── figure S1 ─────────────────────────────────────────────────────────────────
def fig_stability(summary, out: pathlib.Path):
    drugs = sorted({r["drug"] for r in summary})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7),
                             gridspec_kw={"wspace": 0.26})
    for ax, field, name in ((axes[0], "cost", "Cost  2·vMER + MER"),
                            (axes[1], "SR_f1", "Susceptible/resistant F1")):
        tidy(ax, x=True, y=True)
        pk = pick if field == "cost" else pickp          # cost is an index and stays as it is
        pts = [(pk(summary, d, m, "development_cv", field),
                pk(summary, d, m, "holdout", field), m, d)
               for d in drugs for m in MODELS]
        lo = min(min(p[0], p[1]) for p in pts)
        hi = max(max(p[0], p[1]) for p in pts)
        pad = (hi - lo) * 0.09
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=AXIS, lw=0.8,
                ls=(0, (3, 3)), zorder=1)
        for x, y, m, d in pts:
            ax.plot(x, y, "o", ms=4.5, color=SERIES[m], mec=SURFACE, mew=0.7,
                    zorder=3)
        for d, off in zip(UNDERPOWERED, ((7, -9), (7, 5))):
            ax.annotate(d.capitalize(),
                        (pk(summary, d, "evo2", "development_cv", field),
                         pk(summary, d, "evo2", "holdout", field)),
                        xytext=off, textcoords="offset points", fontsize=6.8,
                        color=MUTED)
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlabel("development, 3×5-fold CV" + ("" if field == "cost" else " (%)"))
        ax.set_ylabel("held-out set" + ("" if field == "cost" else " (%)"))
        ax.set_title(name, loc="left", color=INK, pad=7, fontweight="medium")

    top = header(fig, "Development and holdout agree",
                 "Each point is one drug × one model, 30 per panel. Dashed line is "
                 "equality. The only departures are the\ntwo underpowered drugs, "
                 "where the holdout has twelve and fifteen resistant isolates.",
                 pad_in=0.28)
    legend_models(fig, baseline=False, y=-0.01)
    fig.subplots_adjust(top=top, bottom=0.22, left=0.085, right=0.985)
    save(fig, out, "figS1_stability")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workbook",
                    default=str(pathlib.Path.home()
                                / "Desktop/AMRPredict_panel_results.xlsx"))
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    summary, per_class = load(pathlib.Path(args.workbook))
    print(f"[INFO] {len(summary)} summary rows, {len(per_class)} per-class rows")
    out = pathlib.Path(args.out_dir)
    fig_panel(summary, out)
    fig_tradeoff(summary, out)
    fig_ladder(summary, per_class, out)
    fig_newdrugs(summary, out)
    fig_stability(summary, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
