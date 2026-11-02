#!/usr/bin/env python3
"""Extended Data Figure 1 (draft): the three foundation models compared, embeddings alone, held-out sets.

The full comparison behind Figure 1c. Two rows, one per organism, four held-out measures each: essential
agreement, balanced essential agreement (classes with at least ten isolates), mean absolute error, and F1 for
the resistant class. One row per antimicrobial, one mark per model, joined by a line from the lowest to
the highest. Above each panel, Evo 2 against each of the other two: the antimicrobials on which Evo 2 is
better and a two-sided Wilcoxon signed-rank test across antimicrobials.

  M. tuberculosis  embeddings of the reported target loci, no catalogue call
                   (<drug>_<model>_holdout_predictions_embonly1k)
  P. aeruginosa    arm C: embeddings of the chromosomal and acquired-gene targets, no catalogue input
                   (<drug>_<model>_holdout_predictions_armC)

Metrics are computed from the saved held-out predictions with cv_run_architecture.metrics_for, as in
Figure 2; the counts and P values are those of Supplementary Table 3 (Tests). The numbers behind every
mark are written to extended_figure1_values.tsv beside the figure.
"""

from __future__ import annotations

import argparse
import csv
import os
import pathlib
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_figure2 import MTB_BOUND, PA_BOUND, PANEL_SPEC, PCT, heldout, pval  # noqa: E402
from make_panel_figures import (SUBMISSION, fig_h, AXIS, GRID, INK, INK2, ITALIC_BOLD, MODEL_LABEL, MODELS, MUTED,  # noqa: E402
                                SERIES, SURFACE, header, save, tidy)
from make_provisional_figures import NAME, paired_summary  # noqa: E402
from plot_metric_bars import PANELS  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
MARK = {"evo2": dict(marker="o", ms=5.0), "carbon": dict(marker="s", ms=4.4), "dnabert": dict(marker="^", ms=5.0)}
SHORT = {m: MODEL_LABEL[m].replace(" 7B", "").replace(" 8B", "") for m in MODELS}
ROWS = (("mtb", "M. tuberculosis", "{drug}_{m}_holdout_predictions_embonly1k", MTB_BOUND),
        ("pseaer", "P. aeruginosa", "{drug}_{m}_holdout_predictions_armC", PA_BOUND))


def triple(ax, drugs, vals, faint, unit):
    """One row per antimicrobial: the three models as marks, joined from the lowest to the highest."""
    tidy(ax)
    core = [vals[m][i] for m in MODELS for i, d in enumerate(drugs) if d not in faint]
    pad = 0.14 * (max(core) - min(core))
    lim = (min(core) - pad, max(core) + pad)
    for i, d in enumerate(drugs):
        row = [vals[m][i] for m in MODELS]
        al = 0.4 if d in faint else 1.0
        if d in faint and (min(row) < lim[0] or max(row) > lim[1]):
            fmt = "{:.1f}" if unit > 1 else "{:.2f}"
            ax.annotate("off scale: " + ", ".join(fmt.format(v) for v in row), xy=(0.5, i),
                        xycoords=("axes fraction", "data"), ha="center", va="center", fontsize=5.4, color=MUTED)
            continue
        ax.plot([min(row), max(row)], [i, i], color=GRID, lw=2.4, solid_capstyle="round", zorder=1, alpha=al)
        for z, m in enumerate(reversed(MODELS)):          # Evo 2 drawn last, on top
            ax.plot(vals[m][i], i, ls="none", color=SERIES[m],
                    mec=SURFACE, mew=0.7, alpha=al, zorder=3 + z, **MARK[m])
    ax.set_yticks(range(len(drugs)))
    ax.set_yticklabels([NAME.get(d, d.capitalize()) for d in drugs], fontsize=7.2, color=INK)
    ax.set_ylim(-0.7, len(drugs) - 0.3)
    ax.invert_yaxis()
    ax.set_xlim(*lim)
    ax.tick_params(axis="x", labelsize=6.6)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-pred", default=str(HOME / "mtb_slim_20261016/predictions"))
    ap.add_argument("--pa-pred", default=str(HOME / "pa_strict_20261020/predictions_strict"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="extended_figure1_draft")
    ap.add_argument("--width", type=float, default=9.4)
    a = ap.parse_args()
    pred = {"mtb": pathlib.Path(a.mtb_pred), "pseaer": pathlib.Path(a.pa_pred)}

    fig = plt.figure(figsize=(a.width, fig_h(7.2)))
    gs = fig.add_gridspec(2, 5, hspace=0.92, wspace=0.46, width_ratios=[1, 1, 1, 0.16, 1])
    COL = (0, 1, 2, 4)
    values, first = [], {}
    for r, (k, org, stem, bounds) in enumerate(ROWS):
        drugs, faint = list(bounds), PANELS[k]["underpowered"]
        met = {m: {d: heldout(pred[k], stem.format(drug=d, m=m), b) for d, b in bounds.items()} for m in MODELS}
        for n, (key, title, xl, low) in enumerate(PANEL_SPEC):
            ax = fig.add_subplot(gs[r, COL[n]])
            u = 100.0 if key in PCT else 1.0
            f_n = faint if key == "SR_f1" else ()
            vals = {m: [met[m][d][key] * u for d in drugs] for m in MODELS}
            triple(ax, drugs, vals, f_n, u)
            lines = []
            for f_, o in (("evo2", "carbon"), ("evo2", "dnabert"), ("carbon", "dnabert")):
                wins, ties, m_, med, p = paired_summary(drugs, vals[f_], vals[o], low, exclude=f_n,
                                                        true_delta=True, unit=u)
                lines.append(f"{SHORT[f_][0]}–{SHORT[o][0]}: {wins}/{m_}" + (f", {ties} equal" if ties else "")
                             + f", P = {pval(p)}" if SUBMISSION else
                             f"{SHORT[f_]} vs {SHORT[o]}: better {wins}/{m_}" + (f", equal {ties}" if ties else "")
                             + f" · P = {pval(p)}")
                print(f"[INFO] {org} {title}: {SHORT[f_]} vs {SHORT[o]} better {wins}/{m_}, equal {ties}, "
                      f"median {med:+.4f}, P = {p:.4f}")
            ax.annotate("\n".join(lines), xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 3),
                        textcoords="offset points", ha="left", va="bottom", fontsize=5.3, color=INK2,
                        linespacing=1.35)
            ax.set_title(title.replace("Balanced essential ", "Balanced essential\n") if SUBMISSION else title, loc="left",
                         color=INK, pad=29, fontweight="semibold", fontsize=7.6)
            ax.set_xlabel(xl, fontsize=6.2)
            if n:
                ax.set_yticklabels([])
            else:
                first[k] = ax
            ax.annotate("abcdefgh"[r * 4 + n], xy=(-0.06 if n else -0.56, 1.0 if SUBMISSION else 1.38), xycoords="axes fraction",
                        xytext=(0, 38 if SUBMISSION else 0), textcoords="offset points",
                        fontsize=10, fontweight="bold", color=INK, ha="right", va="bottom", annotation_clip=False)
            for i, d in enumerate(drugs):
                values.append([org, d, key] + [round(vals[m][i] / u, 4) for m in MODELS] + [d in f_n])
    top = header(fig, "Extended Data Figure 1 (draft) — the three foundation models, embeddings alone",
                 "Held-out sets. The same random forest given the embeddings of the target loci from Evo 2, Carbon "
                 "or DNABERT-2, with no catalogue call (second row: chromosomal\nand acquired-gene targets). One mark "
                 "per model; the line spans the lowest to the highest. Above each panel: for each pair of "
                 "models, the antimicrobials on which the first is better,\nand two-sided Wilcoxon signed-rank tests across antimicrobials "
                 "(first row ten, but eight for F1, where bedaquiline and delamanid, with 12 and 15 resistant\nisolates, "
                 "are faded and excluded; second row ten). Accuracies and F1 in per cent.", pad_in=0.30)
    fig.subplots_adjust(top=top - 0.110, bottom=0.085, left=0.125, right=0.950)
    # group headers, a divider and a tint tell the MIC measures from the resistance call, as in Figure 2
    fig.canvas.draw()
    top_ax = [x for x in fig.axes if x.get_subplotspec().rowspan.start == 0]
    col_ax = {c: [x for x in fig.axes if x.get_subplotspec().colspan.start == c] for c in COL}
    y_hdr = top_ax[0].get_position().y1 + 0.120
    for cols, label in (((0, 2), "MIC prediction"), ((4, 4), "Resistance call")):
        x0 = col_ax[cols[0]][0].get_position().x0
        x1 = col_ax[cols[1]][0].get_position().x1
        fig.add_artist(Line2D([x0, x1], [y_hdr, y_hdr], transform=fig.transFigure, color=INK2, lw=0.9))
        fig.text((x0 + x1) / 2, y_hdr + 0.007, label, ha="center", va="bottom", fontsize=8.6,
                 fontweight="semibold", color=INK)
    f_ax = col_ax[4]
    bx0 = col_ax[2][0].get_position().x1 + 0.030
    bx1 = min(f_ax[0].get_position().x1 + 0.015, 0.998)
    by0 = min(x.get_position().y0 for x in f_ax) - 0.050
    by1 = y_hdr + 0.004
    fig.add_artist(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, transform=fig.transFigure,
                             facecolor="#f3f5fa", edgecolor="none", zorder=-10))
    for x in f_ax:
        x.patch.set_alpha(0)
    fig.add_artist(Line2D([bx0, bx0], [by0, by1], transform=fig.transFigure, color=AXIS, lw=0.8, ls=(0, (3, 3))))
    for k, (_, org, *_rest) in zip(first, ROWS):
        pos = first[k].get_position()
        fig.text(0.017, (pos.y0 + pos.y1) / 2, org, rotation=90, fontsize=9, fontproperties=ITALIC_BOLD,
                 color=INK, ha="center", va="center")
    fig.legend(handles=[Line2D([], [], ls="none", color=SERIES[m], mec=SURFACE, mew=0.7, label=SHORT[m], **MARK[m])
                        for m in MODELS],
               loc="lower center", ncol=3, fontsize=7.2, handletextpad=0.4, columnspacing=2.2,
               bbox_to_anchor=(0.5, 0.0))
    out = pathlib.Path(a.out_dir)
    save(fig, out, a.stem)
    with open(out / "extended_figure1_values.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["organism", "antimicrobial", "measure"] + [SHORT[m] for m in MODELS] + ["excluded_from_tests"])
        w.writerows(values)
    return 0


if __name__ == "__main__":
    sys.exit(main())
