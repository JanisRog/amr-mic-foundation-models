#!/usr/bin/env python3
"""F1 and one-tier accuracy as grouped bars: one group per drug, one bar per model.

Three panels over a shared x axis: F1, one-tier accuracy, and balanced one-tier
accuracy over the MIC classes with at least ten isolates. The F1 panel also
carries the organism's comparator (TB-Profiler or ARDaP), which makes a resistance
call and so has an F1; it has no MIC prediction and therefore no bar in the other
two. Held-out sets, read from the results workbook, never recomputed.

The balanced panel is there because plain one-tier accuracy is forgiving: most
isolates sit on two or three neighbouring MIC values, so predicting the commonest
one already scores about 0.75 for M. tuberculosis. Averaging over MIC classes
removes that, and is the declared primary measure of MIC prediction.

The comparator is not a fourth hue. As in the other figures it is the incumbent
everything is read against, so it takes an open bar in neutral ink.

Bars start at zero. On a truncated axis the differences between models would look
several times larger than they are, and the point of this figure is how close the
three models sit.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

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

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.3), sharex=True,
                             gridspec_kw={"hspace": 0.30})
    n_base = 0
    for ax, field, title in ((axes[0], "SR_f1", "F1, resistant class"),
                             (axes[1], "OT", "One-tier accuracy  (predicted MIC within one "
                                             "doubling dilution)"),
                             (axes[2], "bal_OT_common", "Balanced one-tier accuracy  (averaged "
                                                        "over MIC classes with at least ten "
                                                        "isolates)")):
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
    axes[-1].set_xticks(range(len(drugs)))
    axes[-1].set_xticklabels([SHORT.get(d, d.capitalize()) for d in drugs], fontsize=7.4,
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
                 f"in the top panel only.{note}", pad_in=0.26)
    fig.subplots_adjust(top=top, bottom=0.135, left=0.065, right=0.99)
    save(fig, pathlib.Path(a.out_dir), a.stem or f"bars_f1_onetier_{a.panel}")
    print(f"[INFO] {a.panel}: {len(drugs)} drugs, comparator bars {n_base}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
