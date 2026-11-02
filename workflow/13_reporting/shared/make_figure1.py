#!/usr/bin/env python3
"""Composite Figure 1 (draft layout): what the approach achieves, and how it works.

  a     essential agreement of the reported model (Evo 2) on the held-out set of
        every antimicrobial, one sub-panel per organism
  b.1   the three embedding models compared, median essential agreement over the
        organism-antimicrobial pairs, linked by an arrow to step 3 of the pipeline
  b.2   the pipeline illustration (an image supplied by the author, placed as it is)
  c     predicted against measured MIC, isoniazid (M. tuberculosis), held-out isolates
  d     the same for ciprofloxacin (P. aeruginosa)

Every panel reads the results workbooks or saved held-out predictions; nothing is
refit. Panels c and d are drawn by the same function as the single-drug figures
and refuse to render if the predictions do not reproduce the workbook row. They
show the antimicrobial with the highest F1 in each organism, i.e. a favourable
case, and the legend should say so.

--variant mtb draws the same figure for M. tuberculosis alone: panel a shows
essential agreement and F1 side by side, b.1 compares the three models on embeddings
with no other feature (development folds, catalogue-derived target loci; median
over the eight antimicrobials compared), c is the isoniazid scatter and d the
isoniazid MIC class distribution.
"""

from __future__ import annotations

import argparse
import csv
import os
import pathlib
import statistics
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, FancyBboxPatch, Patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import (BOLD, FINAL, SUBMISSION, GRID, title_pt, final_title, INK, INK2, ITALIC, ITALIC_BOLD, MODEL_LABEL, MODELS, ORG_COLOUR,  # noqa: E402
                                MUTED, SERIES, SURFACE, load, pick, save, tidy)
from make_panel_grids import load_drug  # noqa: E402
import holdout_arm_comparison as H  # noqa: E402
from plot_class_distribution import (OT_LINE, TRAIN_EDGE, TRAIN_FACE, class_counts,  # noqa: E402
                                     draw_classdist)
from plot_metric_bars import PANELS, SHORT  # noqa: E402
from make_figure3_alt import CODE  # noqa: E402
TALL = SUBMISSION and os.environ.get("FIG_TALL") == "1"          # original tall proportions, Nature text rules
COMPACT = SUBMISSION and not TALL                                  # the 170 mm deep version
from plot_mic_scatter import CATS, EDGE, FACE, LABEL, draw_grid  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
ABBR = {"mtb": "M. tuberculosis", "pseaer": "P. aeruginosa"}
STEP3 = "#7030a0"          # the purple of step 3 in the pipeline illustration


def letter(ax, s, dx=-0.01, dy=1.04):
    ax.annotate(s, xy=(dx, dy), xycoords="axes fraction", fontsize=10, fontweight="bold",
                color=INK, ha="right", va="bottom", annotation_clip=False)


def panel_a(ax, summary, k, ylabel, field="OT", title=None, footnote=None, dagger=(), dagger_note=""):
    """Evo 2 on the held-out sets, every antimicrobial of one organism."""
    tidy(ax, x=False, y=True)
    P = PANELS[k]
    # essential agreement is an MIC measure and rests on every held-out isolate (~1,500 for bedaquiline
    # and delamanid); only the resistance call is thin for the two, so only F1 marks them
    under = [] if field == "OT" else P["underpowered"]
    S = 100                                         # accuracies and F1 in per cent
    for x, d in enumerate(P["drugs"]):
        v = pick(summary, d, "evo2", "holdout", field) * S
        faint = d in under
        ax.bar(x, v, 0.66, color=ORG_COLOUR[k], edgecolor=SURFACE, linewidth=0.5,
               alpha=0.45 if faint else 1.0, zorder=3)
        ax.annotate(f"{v:.1f}", (x, v), xytext=(0, 2), textcoords="offset points",
                    ha="center", va="bottom", fontsize=5.6, color=INK2, zorder=5,
                    bbox=dict(fc=SURFACE, ec="none", pad=0.6) if field == "OT" else None)
    ax.set_xticks(range(len(P["drugs"])))
    if COMPACT:                                     # the codes of Figures 3-5; the names are in the legend
        ax.set_xticklabels([CODE[d] + (" *" if d in under else "") + ("†" if d in dagger else "") for d in P["drugs"]],
                           fontsize=6.0, color=INK)
    else:
        ax.set_xticklabels([SHORT.get(d, d.capitalize()) + (" *" if d in under else "")
                            + (" †" if d in dagger else "") for d in P["drugs"]], fontsize=6.6, color=INK, rotation=38,
                           ha="right", rotation_mode="anchor", linespacing=0.95)
    ax.set_ylim(0, 1.10 * S)
    ax.set_yticks([0, 0.25 * S, 0.5 * S, 0.75 * S, 1.0 * S])
    ax.set_xlim(-0.7, len(P["drugs"]) - 0.3)
    if field == "OT":                               # the FDA criterion for essential agreement (p. 18 of the guidance)
        ax.axhline(90, color="#4f8a5f", lw=0.9, ls=(0, (3, 2)), zorder=2)
        if k == "pseaer":
            ax.annotate("FDA criterion for essential agreement, 90%", xy=(1.0, 90), xycoords=("axes fraction", "data"),
                        xytext=(0, 2.5), textcoords="offset points", ha="right", va="bottom", fontsize=5.8,
                        color="#4f8a5f")
    if ylabel:
        ax.set_ylabel("Essential agreement (%)", fontsize=7.5)
    ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    if title:
        ax.set_title(title, loc="center", fontsize=8.6, fontweight="semibold", color=INK, pad=7)
    else:
        ax.set_title(ABBR[k], loc="center", fontsize=8.6, fontproperties=ITALIC_BOLD,
                     color=INK, pad=7)
    if dagger and dagger_note and not SUBMISSION:
        ax.annotate(dagger_note, xy=(0.0, -0.40), xycoords="axes fraction", fontsize=5.8, color=MUTED,
                    ha="left", va="top", annotation_clip=False)
    if footnote is None:
        footnote = "left" if ylabel else ""
    if under and footnote:
        ax.annotate("* few resistant held-out isolates", xy=(0.0 if footnote == "left" else 1.0, -0.40),
                    xycoords="axes fraction", fontsize=5.8, color=MUTED, ha=footnote, va="top",
                    annotation_clip=False)


def panel_models(ax, books):
    """Median essential agreement per embedding model over the organism-antimicrobial pairs."""
    tidy(ax, x=False, y=True)
    pairs = [(k, d) for k in ("mtb", "pseaer") for d in PANELS[k]["drugs"]
             if d not in PANELS[k]["underpowered"]]
    marker = {"mtb": "o", "pseaer": "^"}
    for j, m in enumerate(MODELS):
        vals = [(k, 100 * pick(books[k], d, m, "holdout", "OT")) for k, d in pairs]
        med = statistics.median(v for _, v in vals)
        ax.bar(j, med, 0.62, color=SERIES[m], edgecolor=SURFACE, linewidth=0.5, zorder=3)
        n = len(vals)
        for i, (k, v) in enumerate(vals):
            ax.plot(j + (i - (n - 1) / 2) * (0.46 / (n - 1)), v, marker[k], ms=2.3,
                    mfc=SURFACE, mec=INK, mew=0.5, zorder=5)
        ax.annotate(f"{med:.1f}", (j, 103), xytext=(-2, 0), textcoords="offset points", rotation=65, rotation_mode="anchor",
                    ha="left", va="bottom", fontsize=6.4, color=INK)
    ax.set_xticks(range(len(MODELS)))
    ax.set_xticklabels([MODEL_LABEL[m].replace(" 7B", "").replace(" 8B", "") for m in MODELS],
                       fontsize=6.2, color=INK, rotation=65, ha="right", rotation_mode="anchor")
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(-0.6, len(MODELS) - 0.4)
    ax.set_ylabel("Essential agreement (%), median", fontsize=7.0)
    ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    return len(pairs)


MTB_BOUND = {"isoniazid": -3.0, "rifampicin": -1.0, "ethambutol": 2.0, "ethionamide": 2.0, "kanamycin": 2.0,
             "amikacin": 0.0, "levofloxacin": 0.0, "moxifloxacin": 0.0, "bedaquiline": -2.0, "delamanid": -3.0}
PA_BOUND = {"ciprofloxacin": -1.0, "levofloxacin": 1.0, "tobramycin": 1.0, "amikacin": 4.0, "imipenem": 2.0,
            "meropenem": 3.0, "ceftazidime": 3.0, "cefepime": 3.0, "aztreonam": 4.0,
            "piperacillin-tazobactam": 4.0}


SP = 1.65                                 # spacing of the bars of one group in b.1


def panel_models_embonly_heldout(ax, mtb_dir, pa_dir):
    """The three models on embeddings alone, held-out essential agreement, one group per organism.

    Computed from the saved held-out predictions (M. tuberculosis: embeddings of the
    reported target loci, no TB-Profiler call, all ten antimicrobials;
    P. aeruginosa: arm C, embeddings of the chromosomal and acquired-gene targets, no
    catalogue input, all ten). One point per antimicrobial; the bar is the median
    within the organism. Pooling the two organisms is avoided on purpose: their levels
    differ (about 0.92 and 0.74) and so does the spread between the models.
    """
    tidy(ax, x=False, y=True)
    src = [("M. tuberculosis", pathlib.Path(mtb_dir), "holdout_predictions_embonly1k", MTB_BOUND, 0.0),
           ("P. aeruginosa", pathlib.Path(pa_dir), "holdout_predictions_armC", PA_BOUND, 5.4)]
    xticks, xlabels, n_pts = [], [], 0
    for name, d_, stem, bounds, x0 in src:
        meds = {}
        for j, m in enumerate(MODELS):
            vals = []
            for drug, b in bounds.items():
                p = H.read_preds(d_ / f"{drug}_{m}_{stem}.tsv", "median")
                vals.append(100 * H.score([p[i] for i in sorted(p)], b)["OT"])
            med = statistics.median(vals)
            meds[m] = med
            x = x0 + SP * j
            ax.bar(x, med, 0.92, color=SERIES[m], edgecolor=SURFACE, linewidth=0.5, zorder=3)
            n = len(vals)
            for i, v in enumerate(vals):
                ax.plot(x + (i - (n - 1) / 2) * (0.70 / (n - 1)), v, "o", ms=2.3, mfc=SURFACE, mec=INK,
                        mew=0.5, zorder=5)
            # just above the highest point of this bar, so every label sits the same distance above its own data
            ax.annotate(f"{med:.1f}", (x, max(vals)), xytext=(3.2, 5), textcoords="offset points", rotation=65,
                        rotation_mode="anchor", ha="left", va="bottom", fontsize=6.2, color=INK)
            xticks.append(x)
            xlabels.append(MODEL_LABEL[m].replace(" 7B", "").replace(" 8B", ""))
        n_pts += n
        # a heading with a rule above each group, as the headings of Figure 2
        ax.plot([x0 - 0.85, x0 + 2 * SP + 0.85], [126, 126], color=INK2, lw=0.9, clip_on=False, zorder=6)
        ax.text(x0 + SP, 128, name, ha="center", va="bottom", fontsize=6.4, fontproperties=ITALIC_BOLD, color=INK)
        print(f"[INFO] b.1 {name}: " + "; ".join(f"{m} {v:.1f}" for m, v in meds.items()))
    ax.set_xticks(xticks)
    if SUBMISSION:                           # vertical names: a shorter box under the panel
        ax.set_xticklabels(xlabels, fontsize=5.5, color=INK, rotation=90, ha="center")
    else:
        ax.set_xticklabels(xlabels, fontsize=5.8, color=INK, rotation=65, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", pad=7)          # names clear of the axis, room for the box above them
    ax.set_ylim(0, 138)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(-0.6, 5.4 + 2 * SP + 0.6)
    ax.set_ylabel("Essential agreement (%), median\n(embeddings only, held-out)", fontsize=6.6, linespacing=1.2)
    ax.yaxis.label.set_y(50 / 138)            # centred on the 0-100 scale, not on the headroom above it
    ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    return n_pts


def embonly_dev(slim_dir, drug, model, field):
    """Mean over the fifteen development folds of the selected configuration."""
    d = "isoniazid_final" if drug == "isoniazid" else drug
    run = f"{model}_embonly1k" if drug == "isoniazid" else f"{model}_both_embonly1k"
    f = pathlib.Path(slim_dir) / f"cv_{d}" / "runs" / run / "folds.tsv"
    rows = [r for r in csv.DictReader(open(f), delimiter="\t") if r["aggregation"] == "median"]
    if len({r["config"] for r in rows}) != 1 or len(rows) != 15:
        raise SystemExit(f"[ERROR] {f}: expected the selected configuration on 15 folds")
    if field == "cost":
        return statistics.mean(2 * float(r["vMER"]) + float(r["MER"]) for r in rows)
    return statistics.mean(float(r[field]) for r in rows)


def panel_models_embonly(ax, slim_dir):
    """The three models on embeddings alone -- no catalogue call, no mutation columns.

    Development folds, catalogue-derived target loci, the configuration chosen by
    the same rule as every other model. Median over the eight antimicrobials
    compared; each antimicrobial a point.
    """
    tidy(ax, x=False, y=True)
    drugs = [d for d in PANELS["mtb"]["drugs"] if d not in PANELS["mtb"]["underpowered"]]
    for j, m in enumerate(MODELS):
        vals = [100 * embonly_dev(slim_dir, d, m, "OT") for d in drugs]
        med = statistics.median(vals)
        ax.bar(j, med, 0.62, color=SERIES[m], edgecolor=SURFACE, linewidth=0.5, zorder=3)
        n = len(vals)
        for i, x in enumerate(vals):
            ax.plot(j + (i - (n - 1) / 2) * (0.46 / (n - 1)), x, "o", ms=2.3, mfc=SURFACE, mec=INK,
                    mew=0.5, zorder=5)
        ax.annotate(f"{med:.1f}", (j, 103), xytext=(-2, 0), textcoords="offset points", rotation=65, rotation_mode="anchor",
                    ha="left", va="bottom", fontsize=6.4, color=INK)
    ax.set_xticks(range(len(MODELS)))
    ax.set_xticklabels([MODEL_LABEL[m].replace(" 7B", "").replace(" 8B", "") for m in MODELS],
                       fontsize=6.2, color=INK, rotation=65, ha="right", rotation_mode="anchor")
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(-0.6, len(MODELS) - 0.4)
    ax.set_ylabel("Essential agreement (%), median\n(embeddings only, development)", fontsize=6.6, linespacing=1.2)
    ax.axhline(0, color=INK2, lw=0.7, zorder=4)


def panel_pipeline(ax):
    """The pipeline as a row of steps; the embedding models carry the bar colours."""
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def box(x, y, w, h, title, body, face="#f4f3ee"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6",
                                    facecolor=face, edgecolor="#b5b3a9", linewidth=0.7, zorder=2))
        ax.text(x + w / 2, y + h - 5.5, title, ha="center", va="top", fontsize=6.5,
                fontweight="semibold", color=INK, zorder=3)
        ax.text(x + w / 2, y + h - 17.5, body, ha="center", va="top", fontsize=5.5, color=INK2,
                linespacing=1.3, zorder=3)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=1,
                    arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.8, mutation_scale=7,
                                    shrinkA=0, shrinkB=0))

    top, h, w, gap = 56, 40, 21.5, 4.6
    xs = [0 + i * (w + gap) for i in range(4)]
    box(xs[0], top, w, h, "Sequencing reads", "Illumina\nshort reads")
    box(xs[1], top, w, h, "Variant calling", "Snippy; consensus\nat allele fraction\n≥ 0.1")
    box(xs[2], top, w, h, "Target loci", "gene body and 300 bp\nupstream; 4–72\nper antimicrobial")
    box(xs[3], top, w, h, "Foundation model", "frozen DNA model;\nmean over positions", face=SURFACE)
    for i in range(3):
        arrow(xs[i] + w, top + h / 2, xs[i + 1], top + h / 2)
    # the three embedding models, in the colours of the bars in the panel beside this one
    cx, cy = xs[3] + 0.9, top + 3.2
    for m, short, cw in (("evo2", "Evo 2", 4.9), ("carbon", "Carbon", 5.6), ("dnabert", "DNABERT-2", 8.4)):
        ax.add_patch(FancyBboxPatch((cx, cy), cw, 8.4,
                                    boxstyle="round,pad=0,rounding_size=1.0", facecolor=SERIES[m],
                                    edgecolor="none", zorder=3))
        ax.text(cx + cw / 2, cy + 4.2, short, ha="center", va="center",
                fontsize=4.2, color=SURFACE, fontweight="semibold", zorder=4)
        cx += cw + 0.45

    bot = 4
    xb = [xs[3], xs[2], xs[1], xs[0]]              # second row runs right to left
    box(xb[0], bot, w, h, "Shared projection", "PCA over all targets;\n4 components\nper target")
    box(xb[1], bot, w, h, "Catalogue call", "TB-Profiler:\nresistant or\nsusceptible")
    box(xb[2], bot, w, h, "Random forest", "regression on\nlog₂ MIC; median\nof the trees")
    box(xb[3], bot, w, h, "Prediction", "MIC, and resistant\nor susceptible\nat the cut-off", face="#e9f1fc")
    arrow(xs[3] + w / 2, top, xs[3] + w / 2, bot + h)
    arrow(xb[0], bot + h / 2, xb[1] + w, bot + h / 2)
    arrow(xb[1], bot + h / 2, xb[2] + w, bot + h / 2)
    arrow(xb[2], bot + h / 2, xb[3] + w, bot + h / 2)
    ax.text(xb[1] + w / 2, bot + h + 2.2, "+ one binary feature", ha="center", va="bottom",
            fontsize=5.2, color=MUTED)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-workbook", default=str(HOME / "AMRPredict_panel_results_wider.xlsx"))
    ap.add_argument("--pa-workbook", default=str(HOME / "AMRPredict_pseaer_strict_results.xlsx"))
    ap.add_argument("--mtb-pred-dir", default=str(HOME / "mtb_wider_figure_inputs"))
    ap.add_argument("--pipeline-image",
                    default=str(pathlib.Path(__file__).resolve().parent.parent
                                / "manuscript/pipeline_illustration.png"))
    ap.add_argument("--pa-pred-dir", default=str(HOME / "pa_strict/predictions_strict"))
    ap.add_argument("--pa-fold-dir", default=str(HOME / "pa_strict/slim_strict"))
    ap.add_argument("--mtb-drug", default="isoniazid:-3.0", help="drug:boundary for panel d")
    ap.add_argument("--pa-drug", default="ciprofloxacin:-1.0", help="drug:boundary for panel e")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="")
    ap.add_argument("--variant", choices=["both", "mtb"], default="both",
                    help="'mtb': M. tuberculosis alone (see the module docstring)")
    ap.add_argument("--mtb-slim", default=str(HOME / "mtb_slim_20261014/slim_arms"),
                    help="development runs, for the embeddings-only model comparison")
    ap.add_argument("--mtb-embonly-pred", default=str(HOME / "mtb_slim_20261016/predictions"),
                    help="held-out predictions of the embeddings-only M. tuberculosis models")
    ap.add_argument("--pa-embonly-pred", default=str(HOME / "pa_strict_20261020/predictions_strict"),
                    help="held-out predictions of the P. aeruginosa arm C (embeddings only) models")
    ap.add_argument("--width", type=float, default=9.4,
                    help="figure width in inches. A journal page is 7.2 wide, so anything "
                         "larger is scaled down in print and its text shrinks by the same factor")
    a = ap.parse_args()

    msum, _ = load(pathlib.Path(a.mtb_workbook))
    books = {"mtb": msum, "pseaer": load(pathlib.Path(a.pa_workbook))[0]}
    pred_dir = pathlib.Path(a.mtb_pred_dir)
    scat = []
    specs = [("mtb", a.mtb_drug, pred_dir, pred_dir)]
    if a.variant != "mtb":
        specs.append(("pseaer", a.pa_drug, pathlib.Path(a.pa_pred_dir), pathlib.Path(a.pa_fold_dir)))
    folds_of = {}
    for k, spec, pdir, fdir in specs:
        drug, bound = spec.split(":")
        y, pred, cat, m, folds_p = load_drug(pdir, fdir, drug, "evo2", float(bound), books[k], 5e-4)
        scat.append((k, drug, float(bound), y, pred, cat, m))
        folds_of[drug] = folds_p

    if os.path.exists(a.pipeline_image):
        img = plt.imread(a.pipeline_image)
    else:                                          # the artwork of panel b is not part of the repository
        img = np.ones((1462, 4271, 3))             # blank, same proportions: the layout is unchanged
        print(f"[WARN] {a.pipeline_image} not found: panel b (pipeline illustration) is left blank")
    W = a.width
    # the class-distribution panel of the M. tuberculosis variant has a right-hand axis
    left, right = 0.6 / W, 1 - (0.62 if a.variant == "mtb" else 0.12) / W
    pipe_frac = 0.70 if COMPACT else (0.72 if TALL else 0.81)           # share of the row given to the image
    h_pipe = pipe_frac * 0.985 * (right - left) * W * img.shape[0] / img.shape[1]
    if TALL:                          # original proportions; spacer rows keep the labels clear of the next row
        heights = [2.0, 0.95, h_pipe, 1.0, 0.455 * (right - left) * W]
        hspace, m_top, m_bot = 0.0, 0.34, 1.25
    elif COMPACT:                     # 170 mm: spacer rows size each gap for what lives in it
        heights = [0.95, 0.40, h_pipe, 0.95, 1.5]
        hspace, m_top, m_bot = 0.0, 0.30, 1.0
    else:
        heights = [2.0, h_pipe, 0.455 * (right - left) * W]
        hspace, m_top, m_bot = 0.44, 0.34 + (0.56 if FINAL else 0.0), 1.10
    ROW = (0, 2, 4) if SUBMISSION else (0, 1, 2)
    # One unit of height_ratios is one inch, so the pipeline row is exactly as tall
    # as the image at its width and the image is not shrunk.
    H = sum(heights) + (0.0 if SUBMISSION else 2 * hspace * (sum(heights) / 3)) + m_top + m_bot
    fig = plt.figure(figsize=(W, H))
    outer = fig.add_gridspec(len(heights), 1, height_ratios=heights, hspace=hspace,
                             left=left, right=right, top=1 - m_top / H, bottom=m_bot / H)

    # a: one sub-panel per organism
    top_row = outer[ROW[0]].subgridspec(1, 2, wspace=0.10)
    ax_mtb = fig.add_subplot(top_row[0])
    ax_pa = fig.add_subplot(top_row[1], sharey=ax_mtb)
    if a.variant == "mtb":
        # the footnote goes under the right-hand panel, clear of the note above b.1
        panel_a(ax_mtb, books["mtb"], "mtb", True, "OT", "Essential agreement", footnote="")
        ax_mtb.set_ylabel("")
        panel_a(ax_pa, books["mtb"], "mtb", False, "SR_f1", "F1, resistant class", footnote="right")
    else:
        panel_a(ax_mtb, books["mtb"], "mtb", True)
        panel_a(ax_pa, books["pseaer"], "pseaer", False, dagger={"levofloxacin", "aztreonam"},
                dagger_note="† with AMRFinderPlus presence (ARDaP does not call it)")
    plt.setp(ax_pa.get_yticklabels(), visible=False)
    letter(ax_mtb, "a", dx=-0.10, dy=1.07)
    if FINAL:
        fig.text(0.008, 1 - 0.16 / H, final_title("Figure 1"), fontsize=title_pt(fig), fontproperties=BOLD, wrap=True, color=INK,
                 va="top", ha="left")

    # b the pipeline, c the models, and the arrow that ties the models to step 3 of the pipeline
    mid = outer[ROW[1]].subgridspec(1, 2, width_ratios=[pipe_frac, 1 - pipe_frac], wspace=0.15)
    cell = mid[1].subgridspec(1, 2, width_ratios=[0.92, 0.08])      # room for two groups of bars
    ax_b = fig.add_subplot(cell[0])
    if a.variant == "mtb":
        panel_models_embonly(ax_b, a.mtb_slim)
        n_pairs = 0
    else:
        n_pairs = panel_models_embonly_heldout(ax_b, a.mtb_embonly_pred, a.pa_embonly_pred)
    ax_p = fig.add_subplot(mid[0])
    ax_p.imshow(img, interpolation="lanczos")
    ax_p.axis("off")
    letter(ax_p, "b", dx=0.0, dy=1.07)
    # shift b to the left so its letter lines up with the letters of a and d; the wider gap goes to c
    ax_p.apply_aspect()
    pm, pp = ax_mtb.get_position(), ax_p.get_position()
    shift = pp.x0 - (pm.x0 - 0.10 * pm.width)
    ax_p.set_position([pp.x0 - shift, pp.y0, pp.width, pp.height])

    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    # the letter c to the left of everything in panel c, its y-axis label included
    lab = ax_b.yaxis.label.get_window_extent(rend)
    letter(ax_b, "c", dx=ax_b.transAxes.inverted().transform((lab.x0 - 6, 0))[0], dy=1.03)
    boxes = [t.get_window_extent(rend) for t in ax_b.get_xticklabels()]
    x0, x1 = min(b.x0 for b in boxes) - 17, max(b.x1 for b in boxes) + 13
    y0, y1 = min(b.y0 for b in boxes) - 19, max(b.y1 for b in boxes) + 11
    inv = fig.transFigure.inverted()
    (fx0, fy0), (fx1, fy1) = inv.transform((x0, y0)), inv.transform((x1, y1))
    fig.add_artist(FancyBboxPatch((fx0, fy0), fx1 - fx0, fy1 - fy0, transform=fig.transFigure,
                                  boxstyle="round,pad=0,rounding_size=0.006", facecolor="none",
                                  edgecolor=STEP3, linewidth=1.1, zorder=6))
    # step 3 of the illustration: its box spans about 40.3-58.9% of the image width
    target = (0.496 * img.shape[1], 0.992 * img.shape[0])
    fig.add_artist(ConnectionPatch(xyA=(fx0, (fy0 + fy1) / 2), coordsA=fig.transFigure,
                                   xyB=target, coordsB=ax_p.transData,
                                   connectionstyle="arc3,rad=-0.32", arrowstyle="-|>",
                                   mutation_scale=9, color=STEP3, linewidth=1.1, zorder=6))

    # c, d: the single-drug scatter figures, reproduced at this size. Their key and
    # their statistics box sit beside the grid there; here the statistics go under
    # each grid and one shared key goes at the bottom of the figure.
    SINGLE_SIDE = 5.34                    # inches, the grid in plot_mic_scatter's own figure
    # In the M. tuberculosis variant the scatter is narrower, and panel d is the same
    # antimicrobial's MIC class distribution at its natural, wider proportions.
    solo = a.variant == "mtb"
    low = outer[ROW[2]].subgridspec(1, 2, width_ratios=[0.40, 0.60] if solo else [0.5, 0.5],
                               wspace=0.30 if solo else 0.20)
    for col, (lab, (k, drug, bound, y, pred, cat, m)) in enumerate(zip("de", scat)):
        ax = fig.add_subplot(low[col])
        box = ax.get_position()
        side = min(box.width * W, box.height * H)
        draw_grid(ax, y, pred, cat, bound, counts=True, fs=side / SINGLE_SIDE, scale_counts=True)
        ax.tick_params(labelsize=6.6)
        ax.set_xlabel("Measured MIC (log₂)", fontsize=7.6)
        ax.set_ylabel("Predicted MIC (log₂)", fontsize=7.6)
        ax.set_title(f"{drug.capitalize()}", loc="left", fontsize=8.6, color=INK, pad=6,
                     fontweight="semibold")
        ax.annotate(ABBR[k], xy=(1.0, 1.0), xycoords="axes fraction", xytext=(0, 6),
                    textcoords="offset points", ha="right", va="bottom", fontsize=7.4,
                    fontproperties=ITALIC, color=INK2)
        stats = (f"n = {m['n']:,}  ({m['n_resistant']} resistant)    exact = {100 * m['exact']:.1f}%    "
                 f"EA = {100 * m['OT']:.1f}%    bEA = {100 * m['bal_OT']:.1f}%\n"
                 f"MAE = {m['MAE']:.2f}    vMER = {100 * m['vMER']:.1f}%    MER = {100 * m['MER']:.1f}%    "
                 f"F1 = {100 * m['SR_f1']:.1f}%")
        ax.annotate(stats, xy=(0.5, 0.0), xycoords="axes fraction", xytext=(0, -27),
                    textcoords="offset points", ha="center", va="top", fontsize=6.2, color=INK2,
                    linespacing=1.55, annotation_clip=False)
        letter(ax, lab, dx=-0.10, dy=1.03)

    if solo:
        k, drug, bound, y, pred, cat, m = scat[0]
        ax = fig.add_subplot(low[1])
        ks, tr, stacks, hold, ot, n_train = class_counts(folds_of[drug], y, cat, bound)
        ax2 = draw_classdist(ax, ks, tr, stacks, hold, ot, bound, annotate=False, fs=0.95)
        ax.tick_params(labelsize=6.8)
        ax2.tick_params(labelsize=6.8)
        ax.set_xlabel("MIC class (log₂)", fontsize=7.6)
        ax.set_ylabel("Isolates", fontsize=7.6)
        ax2.set_ylabel("Essential agreement (%)", color=OT_LINE, fontsize=7.6)
        ax.set_title(f"{drug.capitalize()}, MIC classes", loc="left", fontsize=8.6, color=INK, pad=6,
                     fontweight="semibold")
        ax.annotate(f"development n = {n_train:,}", xy=(1.0, 1.0), xycoords="axes fraction", xytext=(0, 6),
                    textcoords="offset points", ha="right", va="bottom", fontsize=6.6, color=MUTED)
        letter(ax, "e", dx=-0.10, dy=1.03)

    handles = [Line2D([], [], marker="o", ls="none", ms=5, mfc=FACE[k], mec=EDGE[k], mew=0.7,
                      label=LABEL[k]) for k, _, _, _ in CATS]
    handles += [Patch(facecolor="#cde2fb", alpha=0.45, edgecolor="#9ec5f4", lw=0.7,
                      label="± one dilution zone"),
                Line2D([], [], color="#d03b3b", lw=1.0, ls=(0, (1.6, 1.6)), label="S/R breakpoint")]
    if a.variant == "mtb":
        handles += [Patch(facecolor=TRAIN_FACE, edgecolor=TRAIN_EDGE, lw=0.6, label="development set (e)"),
                    Line2D([], [], color=OT_LINE, lw=1.3, marker="o", ms=3.5,
                           label="class accuracy (e, right axis)"),
                    Line2D([], [], marker="o", ls="none", ms=3.0, mfc=SURFACE, mec=INK, mew=0.5,
                           label="one antimicrobial (c)")]
    else:
        handles += [Line2D([], [], marker="o", ls="none", ms=3.0, mfc=SURFACE, mec=INK, mew=0.5,
                           label="one antimicrobial (c)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=6.2, handletextpad=0.5,
               columnspacing=1.6, labelspacing=0.55, bbox_to_anchor=(0.53, 0.0))
    save(fig, pathlib.Path(a.out_dir), a.stem or ("figure1_mtb_example" if a.variant == "mtb"
                                                      else "figure1_example"))
    print(f"[INFO] panel b: median over {n_pairs} organism-antimicrobial pairs; panels c, d: "
          + "; ".join(f"{drug} ({len(y):,} held-out isolates)" for _, drug, _, y, *_ in scat))
    return 0


if __name__ == "__main__":
    sys.exit(main())
