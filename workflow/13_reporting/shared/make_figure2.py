#!/usr/bin/env python3
"""Figure 2 (draft): what the embeddings add over an explicit table of the variants in the same loci.

Two rows, one per organism, four held-out measures each: essential agreement, balanced essential
agreement (classes with at least ten isolates), mean absolute error, and F1 for the resistant
class as the only binary measure. One row per antimicrobial, two marks joined by a line:

  M. tuberculosis  Evo 2 embeddings of the reported target loci, no catalogue call
                   (<drug>_evo2_holdout_predictions_embonly1k)
                   against every variant in the same loci, no embeddings, no catalogue call
                   (<drug>_varonly_holdout_predictions_gwasfiltered)
  P. aeruginosa    Evo 2 embeddings of the chromosomal targets only (arm A)
                   against every variant in the same chromosomal loci (varchromo): no acquired
                   genes and no catalogue input on either side

Metrics are computed from the saved per-isolate held-out predictions with
cv_run_architecture.metrics_for, the function that scores every run. The numbers behind every
mark are written to figure2_values.tsv beside the figure.
"""

from __future__ import annotations

import argparse
import csv
import os
import pathlib
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from cv_run_architecture import metrics_for  # noqa: E402
from make_panel_figures import SUBMISSION, fig_h, AXIS, GREY, INK, INK2, ITALIC_BOLD, PURPLE, SURFACE, header, save  # noqa: E402
from make_provisional_figures import dumbbell, paired_summary  # noqa: E402

# embeddings only in purple, variants only in grey: the colours of Figure 4, and not the blue of the
# reported model (embeddings with the catalogue call) in Figures 1 and 3
A_EMB = dict(marker="o", ms=5, color=PURPLE, mec=SURFACE, mew=0.8, ls="none")
A_MUT = dict(marker="s", ms=4.6, color=GREY, mec=SURFACE, mew=0.8, ls="none")
from plot_metric_bars import PANELS  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
MTB_BOUND = {"isoniazid": -3.0, "rifampicin": -1.0, "ethambutol": 2.0, "ethionamide": 2.0, "kanamycin": 2.0,
             "amikacin": 0.0, "levofloxacin": 0.0, "moxifloxacin": 0.0, "bedaquiline": -2.0, "delamanid": -3.0}
PA_BOUND = {"ciprofloxacin": -1.0, "levofloxacin": 1.0, "tobramycin": 1.0, "amikacin": 4.0, "imipenem": 2.0,
            "meropenem": 3.0, "ceftazidime": 3.0, "cefepime": 3.0, "aztreonam": 4.0,
            "piperacillin-tazobactam": 4.0}
PANEL_SPEC = (("OT", "Essential agreement", "%; higher is better", False),
              ("bal_OT_common", "Balanced essential agreement", "%; higher is better", False),
              ("MAE", "Mean absolute error", "doubling dilutions; lower is better", True),
              ("SR_f1", "F1, resistant class", "%; higher is better", False))
PCT = ("OT", "bal_OT_common", "SR_f1")
ROWS = (("mtb", "M. tuberculosis", "{drug}_evo2_holdout_predictions_embonly1k",
         "{drug}_varonly_holdout_predictions_gwasfiltered", MTB_BOUND),
        ("pseaer", "P. aeruginosa", "{drug}_evo2_holdout_predictions_armA",
         "{drug}_varchromo_holdout_predictions", PA_BOUND))


def heldout(pred_dir: pathlib.Path, name: str, bound: float) -> dict:
    df = pd.read_csv(pred_dir / f"{name}.tsv", sep="\t")
    df = df[df["aggregation"] == "median"]
    if df["config"].nunique() != 1 or not df["sample_id"].is_unique:
        raise SystemExit(f"[ERROR] {name}: expected one configuration and one row per isolate")
    return metrics_for(df["log2_mic"].to_numpy(), df["pred"].to_numpy(), df["raw"].to_numpy(),
                       bound, min_class_n=10)


def pval(p: float) -> str:
    return f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-pred", default=str(HOME / "mtb_slim_20261016/predictions"))
    ap.add_argument("--pa-pred", default=str(HOME / "pa_strict_20261020/predictions_strict"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="figure2_draft")
    ap.add_argument("--width", type=float, default=9.4)
    a = ap.parse_args()
    pred = {"mtb": pathlib.Path(a.mtb_pred), "pseaer": pathlib.Path(a.pa_pred)}

    fig = plt.figure(figsize=(a.width, fig_h(6.9)))
    # four measures in five columns: the empty one is the gap between MIC prediction and
    # the resistance call
    gs = fig.add_gridspec(2, 5, hspace=0.72, wspace=0.46, width_ratios=[1, 1, 1, 0.16, 1])
    COL = (0, 1, 2, 4)
    values, first = [], {}
    for r, (k, org, emb_name, var_name, bounds) in enumerate(ROWS):
        drugs, faint = list(bounds), PANELS[k]["underpowered"]
        emb = {d: heldout(pred[k], emb_name.format(drug=d), b) for d, b in bounds.items()}
        var = {d: heldout(pred[k], var_name.format(drug=d), b) for d, b in bounds.items()}
        for n, (key, title, xl, low) in enumerate(PANEL_SPEC):
            ax = fig.add_subplot(gs[r, COL[n]])
            u = 100.0 if key in PCT else 1.0       # accuracies in per cent, F1 as 0-1
            x, y = [emb[d][key] * u for d in drugs], [var[d][key] * u for d in drugs]
            # bedaquiline and delamanid have ~1,500 held-out isolates each, so the MIC measures use all of
            # them; only the resistance call (F1) rests on their 12 and 15 resistant isolates
            f_n = faint if key == "SR_f1" else ()
            dumbbell(ax, drugs, x, y, a_style=A_EMB, b_style=A_MUT, title=title, xlabel=xl,
                     better_low=low, faint=f_n, true_delta=True, unit=u)
            wins, ties, m, med, p = paired_summary(drugs, x, y, low, exclude=f_n, true_delta=True, unit=u)
            ax.annotate(f"better {wins}/{m}" + (f", equal {ties}" if ties else "")
                        + ("\n" if SUBMISSION else " · ") + (f"median {med:+.3f}" if u == 1.0 else f"median {med:+.1f}%") + f" · P = {pval(p)}",
                        xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points",
                        ha="left", va="bottom", fontsize=5.4, color=INK2, linespacing=1.3)
            ax.set_title(title.replace("Balanced essential ", "Balanced essential\n") if SUBMISSION else title, loc="left",
                         color=INK, pad=13 + (14 if SUBMISSION else 0), fontweight="semibold", fontsize=7.6)
            ax.set_xlabel(xl, fontsize=6.2)
            if n:
                ax.set_yticklabels([])
            else:
                first[k] = ax
            ax.annotate("abcdefgh"[r * 4 + n], xy=(-0.06 if n else -0.56, 1.0 if SUBMISSION else 1.22), xycoords="axes fraction",
                        xytext=(0, 36 if SUBMISSION else 0), textcoords="offset points",
                        fontsize=10, fontweight="bold", color=INK, ha="right", va="bottom",
                        annotation_clip=False)
            print(f"[INFO] {org} {title}: better {wins}/{m}, equal {ties}, median {med:+.4f}, P = {p:.4f}")
            for d, xv, yv in zip(drugs, x, y):
                values.append([org, d, key, round(xv / u, 4), round(yv / u, 4), d in faint])
    top = header(fig, "Figure 2 (draft) — embeddings against the variants in the same loci",
                 "Held-out sets. The same random forest given either Evo 2 embeddings of the target loci or every "
                 "variant in those loci as binary columns; neither has a\ncatalogue call. Second row: chromosomal "
                 "loci only on both sides. Numbers at the right: embeddings minus variants (accuracy in percentage\npoints; for the "
                 "mean absolute error lower is better, so a negative difference favours the embeddings); dark type favours the "
                 "embeddings,\ngrey type the variants. Counts, medians and two-sided Wilcoxon tests cover the antimicrobials compared "
                 "(first row ten; F1 eight, bedaquiline and delamanid,\n12 and 15 resistant isolates, faded and excluded; "
                 "second row ten).", pad_in=0.30)
    fig.subplots_adjust(top=top - 0.085, bottom=0.085, left=0.125, right=0.925)
    # group headers, a divider and a tint tell the MIC measures from the resistance call
    fig.canvas.draw()
    rows = [(first[k].get_position(), first[k]) for k in first]
    top_ax = [a_ for a_ in fig.axes if a_.get_subplotspec().rowspan.start == 0]
    bot_ax = [a_ for a_ in fig.axes if a_.get_subplotspec().rowspan.start == 1]
    col_ax = {c: [a_ for a_ in fig.axes if a_.get_subplotspec().colspan.start == c] for c in COL}
    y_hdr = top_ax[0].get_position().y1 + 0.100
    for cols, label in (((0, 2), "MIC prediction"), ((4, 4), "Resistance call")):
        x0 = col_ax[cols[0]][0].get_position().x0
        x1 = col_ax[cols[1]][0].get_position().x1 + 0.035
        fig.add_artist(Line2D([x0, x1], [y_hdr, y_hdr], transform=fig.transFigure, color=INK2, lw=0.9))
        fig.text((x0 + x1) / 2, y_hdr + 0.007, label, ha="center", va="bottom", fontsize=8.6,
                 fontweight="semibold", color=INK)
    f_ax = col_ax[4]
    bx0 = col_ax[2][0].get_position().x1 + 0.056      # clear of panel c's difference labels
    bx1 = min(f_ax[0].get_position().x1 + 0.052, 0.998)
    by0 = min(a_.get_position().y0 for a_ in f_ax) - 0.050
    by1 = y_hdr + 0.004
    fig.add_artist(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, transform=fig.transFigure,
                             facecolor="#f3f5fa", edgecolor="none", zorder=-10))
    for a_ in f_ax:
        a_.patch.set_alpha(0)
    xs = bx0
    fig.add_artist(Line2D([xs, xs], [by0, by1], transform=fig.transFigure, color=AXIS, lw=0.8,
                          ls=(0, (3, 3))))
    for k, (_, org, *_rest) in zip(first, ROWS):
        pos = first[k].get_position()
        fig.text(0.017, (pos.y0 + pos.y1) / 2, org, rotation=90, fontsize=9, fontproperties=ITALIC_BOLD,
                 color=INK, ha="center", va="center")
    fig.legend(handles=[Line2D([], [], label="Embeddings only (Evo 2)", **A_EMB),
                        Line2D([], [], label="Variants only, same loci", **A_MUT)],
               loc="lower center", ncol=2, fontsize=7.2, handletextpad=0.5, columnspacing=2.2,
               bbox_to_anchor=(0.5, 0.0))
    out = pathlib.Path(a.out_dir)
    save(fig, out, a.stem)
    with open(out / "figure2_values.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["organism", "antimicrobial", "measure", "embeddings_only", "variants_only", "excluded_from_tests"])
        w.writerows(values)
    return 0


if __name__ == "__main__":
    sys.exit(main())
