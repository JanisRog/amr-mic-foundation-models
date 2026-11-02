#!/usr/bin/env python3
"""The per-drug figures as two single sheets, one row of the panel per drug.

Ten separate files are right for looking at one drug; they are wrong for seeing
the panel as a panel, which is the claim the paper makes. These are the same
figures, drawn by the same functions -- `draw_grid` and `draw_classdist` are
imported, not reimplemented, so a change to a single-drug figure appears here
without anybody remembering to make it twice.

Two things are dropped at sheet size rather than shrunk:

  **counts inside the bubbles**, which no longer fit, so the sheet carries the
  size key the single figures do not need.
  **per-panel metrics boxes**, replaced by the numbers in each panel title --
  eight lines of text per panel would be most of the sheet.

Every panel is still verified against its workbook row before it is drawn. Ten
verifications per sheet, same tolerance, same refusal.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import INK, INK2, MUTED, SURFACE, header, load, save, tidy
from plot_class_distribution import (OT_LINE, STACK, TRAIN_EDGE, TRAIN_FACE,
                                     class_counts, draw_classdist)
from plot_mic_scatter import (CATS, CHECK, EDGE, FACE, LABEL, SIZE_BINS,
                              categorise, choose_aggregation, draw_grid)
from plot_shap_beeswarm import (REGION, draw_beeswarm, locus_label,
                                normalise, rank_loci)

# drug, boundary; order is the panel's, wide-range drugs first. --panel picks a
# set, --drug-list overrides it -- PA and MTB share two drug names with different
# breakpoints, so a hardcoded table silently mislabels whichever panel it is not.
PANELS = {
    "mtb": [("isoniazid", -3.0), ("rifampicin", -1.0), ("ethambutol", 2.0),
            ("ethionamide", 2.0), ("kanamycin", 2.0), ("amikacin", 0.0),
            ("levofloxacin", 0.0), ("moxifloxacin", 0.0), ("bedaquiline", -2.0),
            ("delamanid", -3.0)],
    "pseaer": [("ciprofloxacin", -1.0), ("levofloxacin", 1.0),
               ("tobramycin", 1.0), ("imipenem", 2.0), ("meropenem", 3.0),
               ("ceftazidime", 3.0), ("cefepime", 3.0), ("aztreonam", 4.0),
               ("amikacin", 4.0), ("piperacillin-tazobactam", 4.0)],
}
DRUGS = PANELS["mtb"]
# What each panel's models were fitted on, for the sheet subtitles. P. aeruginosa
# has no catalogue call; printing "TB-Profiler rule bit" under its figures was wrong.
FEATURES = {"mtb": "shared PCA4 + TB-Profiler rule bit",
            "pseaer": "shared PCA4 + ARDaP call (levofloxacin and aztreonam: AMRFinderPlus presence columns)"}
# what region "other" is in each panel's SHAP key: MTB's call, or P. aeruginosa's ARDaP call
CALL_LABEL = {"mtb": "TB-Profiler rule bit", "pseaer": "ARDaP call"}
PANEL_CALL_LABEL = CALL_LABEL["mtb"]
PANEL = "mtb"
# one sentence about the leading feature, true of each panel's reported model
LEAD_NOTE = {"mtb": "Every drug's leading locus is the one its catalogue names. ",
             "pseaer": "\nARDaP's call is shown as its own row; levofloxacin and aztreonam, which it does not "
                       "call, carry the AMRFinderPlus columns instead.\n"}
FEATURES_LABEL = FEATURES["mtb"]
PREFIX = ""
# drug -> the model whose SHAP run belongs to the reported figure; filled from
# --model-map so the beeswarm sheet cannot quietly plot a different model than
# the scatter sheet beside it.
SHAP_MODEL: dict[str, str] = {}
# How many loci each drug shows on the combined sheet. The single-drug figures
# still show every target -- that is where the full attribution belongs. Ten
# panels at forty-five rows each is a page nobody can read, and the tail is
# near-zero by construction since the rows are ordered by mean |SHAP|.
SHAP_TOP = 10


def load_drug(pred_dir: pathlib.Path, fold_dir: pathlib.Path, drug: str,
              model: str, boundary: float, summary, tol: float):
    """Predictions for one drug, verified against the workbook before use."""
    pred_p = pred_dir / f"{drug}_{model}_holdout_predictions.tsv"
    if not pred_p.exists():                       # isoniazid's redo is named apart
        alt = pred_dir / f"{drug}_{model}_holdout_predictions_tbrule.tsv"
        if not alt.exists():
            raise SystemExit(f"[ERROR] no predictions for {drug}: {pred_p}")
        pred_p = alt
    folds_p = fold_dir / f"cv_{drug}" / "cv_folds_final.tsv"
    if not folds_p.exists():
        folds_p = fold_dir / f"cv_{drug}_final" / "cv_folds_final.tsv"
    if not folds_p.exists():
        raise SystemExit(f"[ERROR] no fold file for {drug} under {fold_dir}")

    df = pd.read_csv(pred_p, sep="\t")
    want = next((r for r in summary if r["drug"] == drug
                 and r["model"] == model and r["split"] == "holdout"), None)
    if want is None:
        raise SystemExit(f"[ERROR] no holdout row for {drug}/{model}")
    want = {k: (None if want.get(k) in ("", None) else float(want[k]))
            for k, _ in CHECK}
    print(f"--- {drug}")
    agg, m = choose_aggregation(df, boundary, want, tol)
    g = df[df["aggregation"] == agg]
    y = g["log2_mic"].to_numpy(float)
    pred = g["pred"].to_numpy(float)
    return y, pred, categorise(y, pred, boundary), m, folds_p


def panel_title(ax, drug, m, fields):
    ax.set_title(drug.capitalize(), loc="left", fontsize=9.5, color=INK,
                 pad=13, fontweight="semibold")
    ax.annotate("  ·  ".join(f"{lab} {100 * m[k]:.1f}%" for k, lab in fields),
                xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 4),
                textcoords="offset points", fontsize=7, color=MUTED)


def sheet_scatters(data, out: pathlib.Path):
    fig = plt.figure(figsize=(11.4, 15.2))
    gs = fig.add_gridspec(4, 3, hspace=0.40, wspace=0.28)
    for i, (drug, _) in enumerate(DRUGS):
        y, pred, cat, m, _ = data[drug]
        ax = fig.add_subplot(gs[i // 3, i % 3])
        draw_grid(ax, y, pred, cat, dict(DRUGS)[drug], counts=False, fs=0.62)
        ax.tick_params(labelsize=6.5)
        panel_title(ax, drug, m, (("exact", "exact"), ("OT", "OT"),
                                  ("vMER", "vMER"), ("SR_f1", "F1")))
        ax.set_xlabel("Measured MIC (log₂)", fontsize=8)
        ax.set_ylabel("Predicted MIC (log₂)", fontsize=8)

    key = fig.add_subplot(gs[3, 1:])
    key.axis("off")
    handles = [Line2D([], [], marker="o", ls="none", ms=8, mfc=FACE[k],
                      mec=EDGE[k], mew=0.9, label=LABEL[k])
               for k, _, _, _ in CATS]
    handles += [Patch(facecolor="#cde2fb", alpha=0.35, edgecolor="#9ec5f4",
                      lw=0.8, label="± one dilution zone"),
                Line2D([], [], color="#d03b3b", lw=1.1, ls=(0, (1.6, 1.6)),
                       label="S/R breakpoint")]
    first = key.legend(handles=handles, loc="upper left",
                       bbox_to_anchor=(0.02, 1.0), ncol=1, fontsize=9,
                       handletextpad=0.8, labelspacing=1.0, borderpad=0.0)
    key.add_artist(first)
    # The sheet needs the size key the single figures dropped: at this scale a
    # bubble cannot hold its own count.
    key.legend(handles=[Line2D([], [], marker="o", ls="none",
                               ms=(sz * 0.62 ** 2) ** 0.5 * 0.85,
                               mfc="#c9c8c1", mec="#8a8980", mew=0.7, label=lab)
                        for _, _, lab, sz in SIZE_BINS],
               loc="upper left", bbox_to_anchor=(0.45, 1.0), ncol=1, fontsize=9,
               handletextpad=1.4, labelspacing=1.35, borderpad=0.0,
               title="isolates in cell", title_fontsize=8)

    top = header(fig, "Predicted against measured MIC, every drug",
                 f"Evo 2 7B on the held-out set, {FEATURES_LABEL}. "
                 "Bubble size is binned by the number of isolates in the "
                 "cell.\nEach panel is drawn from its own saved predictions and "
                 "is not rendered unless the metrics recomputed from them "
                 "reproduce that drug's reported row.", pad_in=0.30)
    fig.subplots_adjust(top=top, bottom=0.035, left=0.065, right=0.985)
    save(fig, out, PREFIX + "sheet_scatters")


def sheet_classdist(data, out: pathlib.Path):
    fig = plt.figure(figsize=(11.4, 14.6))
    gs = fig.add_gridspec(6, 2, hspace=0.62, wspace=0.22,
                          height_ratios=[1, 1, 1, 1, 1, 0.22])
    for i, (drug, boundary) in enumerate(DRUGS):
        y, pred, cat, m, folds_p = data[drug]
        ks, tr, stacks, hold, ot, n_train = class_counts(folds_p, y, cat, boundary)
        ax = fig.add_subplot(gs[i // 2, i % 2])
        ax2 = draw_classdist(ax, ks, tr, stacks, hold, ot, boundary,
                             annotate=False, fs=0.8)
        ax.tick_params(labelsize=7)
        ax.set_xlabel("log₂ MIC class", fontsize=8)
        if i % 2 == 0:
            ax.set_ylabel("isolates", fontsize=8)
        else:
            ax2.set_ylabel("essential agreement", color=OT_LINE, fontsize=8)
        panel_title(ax, drug, m, (("bal_OT", "bEA"), ("OT", "EA")))
        ax.annotate(f"development n = {n_train:,}", xy=(1.0, 1.0),
                    xycoords="axes fraction", xytext=(0, 4), ha="right",
                    textcoords="offset points", fontsize=7, color=MUTED)

    key = fig.add_subplot(gs[5, :])
    key.axis("off")
    handles = [Patch(facecolor=TRAIN_FACE, edgecolor=TRAIN_EDGE, lw=0.7,
                     label="development cohort")]
    handles += [Patch(facecolor=FACE[c], edgecolor=EDGE[c], lw=0.6,
                      label=LABEL[c]) for c in STACK]
    handles += [Line2D([], [], color=OT_LINE, lw=1.6, marker="o", ms=5,
                       label="essential agreement by class (right axis)")]
    key.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.6),
               ncol=4, fontsize=9, handletextpad=0.8, columnspacing=2.0,
               labelspacing=0.9, borderpad=0.0)

    top = header(fig, "MIC class distribution and holdout accuracy, every drug",
                 "Grey is the development cohort, the stack is the held-out set "
                 "split by the same categories as the confusion grids, and the "
                 "line\nreads against the right-hand axis. The averages the "
                 "panel reports are taken over these very unequal classes.",
                 pad_in=0.30)
    fig.subplots_adjust(top=top, bottom=0.035, left=0.065, right=0.945)
    save(fig, out, PREFIX + "sheet_classdist")


def sheet_shap(pred_dir: pathlib.Path, out: pathlib.Path, top: int = 0):
    """The beeswarms as one sheet. Drugs without a SHAP run are skipped.

    SHAP comes from its own cluster job, separate from the predictions, so a
    partially-run panel is the normal state rather than an error. The sheet says
    which drugs are on it instead of failing or, worse, silently drawing nine.
    """
    have = [(d, b) for d, b in DRUGS
            if (pred_dir / f"{d}_{SHAP_MODEL.get(d, 'evo2')}_shap_long.tsv").exists()]
    if not have:
        print("[skip] no SHAP outputs found")
        return
    if len(have) < len(DRUGS):
        missing = [d for d, _ in DRUGS if (d, dict(DRUGS)[d]) not in have]
        print(f"[WARN] SHAP sheet covers {len(have)} of {len(DRUGS)} drugs; "
              f"missing {missing}")

    rows = (len(have) + 1) // 2
    fig = plt.figure(figsize=(11.4, 2.15 * rows + 1.9))
    gs = fig.add_gridspec(rows + 1, 2, hspace=0.75, wspace=0.42,
                          height_ratios=[1] * rows + [0.22])
    for i, (drug, _) in enumerate(have):
        long = pd.read_csv(pred_dir / f"{drug}_{SHAP_MODEL.get(drug, 'evo2')}_shap_long.tsv",
                           sep="\t").fillna({"gene": ""})
        long = normalise(long)
        order, meta = rank_loci(long, top=top)
        long = long[long["locus"].isin(order)]
        ax = fig.add_subplot(gs[i // 2, i % 2])
        draw_beeswarm(ax, long, order, s=2.6, alpha=0.4)
        ax.set_yticklabels([locus_label(meta, lo) for lo in order], fontsize=7)
        ax.tick_params(axis="x", labelsize=7)
        ax.set_title(drug.capitalize(), loc="left", fontsize=9.5, color=INK,
                     pad=6, fontweight="semibold")
        ax.set_xlabel("SHAP value", fontsize=7.5)

    # Driven by what is actually drawn, not a fixed triple: PA adds acquired
    # sequence, point mutations and presence columns, and a hardcoded list left
    # three of the five categories on the page with no key.
    drawn = set()
    for drug, _ in DRUGS:
        f = pred_dir / f"{drug}_{SHAP_MODEL.get(drug, 'evo2')}_shap_long.tsv"
        if f.exists():
            drawn |= set(normalise(pd.read_csv(f, sep="\t", usecols=["target_id", "region"])
                                   .assign(locus="", gene=""))["region"])
    seen = [r for r in REGION if r in drawn]
    key = fig.add_subplot(gs[rows, :])
    key.axis("off")
    key.legend(handles=[Line2D([], [], marker="o", ls="none", ms=6,
                               color=REGION[r][0],
                               label=PANEL_CALL_LABEL if r == "other"
                                     else REGION[r][1])
                        for r in seen],
               loc="upper center", bbox_to_anchor=(0.5, 1.7), ncol=min(len(seen), 3),
               fontsize=9, handletextpad=0.8, columnspacing=2.4, borderpad=0.0)

    shown = (f"the {top} highest-ranked loci per drug" if top
             else "every locus")
    top_y = header(fig, "What the model used, every drug",
                   "TreeSHAP on each reported model, one point per held-out "
                   "isolate. One row per locus, gene body above and 300 bp "
                   f"upstream\nbelow, ordered by mean |SHAP|; {shown}. "
                   + LEAD_NOTE[PANEL] +
                   "The single-drug figures show the full attribution.", pad_in=0.22)
    # Left margin sized for the longest label actually drawn. At 0.115 the PA
    # sheet clipped "mexZ_G195E (NFOBNJ_02085)" to "exZ_G195E" -- locus tags make
    # these labels far longer than MTB's gene names, and a clipped gene name on a
    # figure about which gene the model used is the one thing it cannot afford.
    fig.subplots_adjust(top=top_y, bottom=0.03, left=0.155, right=0.985)
    save(fig, out, PREFIX + "sheet_shap")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", required=True)
    ap.add_argument("--fold-dir", required=True)
    ap.add_argument("--panel", choices=sorted(PANELS), default="mtb")
    ap.add_argument("--drug-list", default="",
                    help="'drug:boundary,...' to override the panel table")
    ap.add_argument("--stem-prefix", default="",
                    help="prefix for the output stems, e.g. 'pa_'")
    ap.add_argument("--workbook", default=str(pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs")) / "AMRPredict_panel_results.xlsx"))
    ap.add_argument("--model", default="evo2")
    ap.add_argument("--model-map", default="",
                    help="drug=model[,drug=model...]. The reported model is "
                         "chosen per drug on development cost, so a panel is not "
                         "necessarily one model throughout -- PA reports carbon "
                         "for amikacin and meropenem and evo2 for the other "
                         "eight. Drugs not named here use --model")
    ap.add_argument("--tolerance", type=float, default=5e-4)
    ap.add_argument("--shap-top", type=int, default=SHAP_TOP,
                    help="loci per drug on the combined SHAP sheet; 0 shows "
                         "every one. The single-drug beeswarms are unaffected")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    global DRUGS
    DRUGS = PANELS[args.panel]
    if args.drug_list:
        DRUGS = [(d, float(b)) for d, b in
                 (x.split(":") for x in args.drug_list.split(","))]
    global PREFIX
    PREFIX = args.stem_prefix
    global FEATURES_LABEL, PANEL_CALL_LABEL, PANEL
    FEATURES_LABEL = FEATURES[args.panel]
    PANEL_CALL_LABEL = CALL_LABEL[args.panel]
    PANEL = args.panel
    print(f"[INFO] panel {args.panel}: {len(DRUGS)} drugs")

    mmap = dict(x.split("=") for x in args.model_map.split(",")) \
        if args.model_map else {}
    bad = [d for d in mmap if d not in {x for x, _ in DRUGS}]
    if bad:
        raise SystemExit(f"[ERROR] --model-map names drugs not in the panel: {bad}")
    if mmap:
        print("[INFO] per-drug models: " + ", ".join(
            f"{d}={mmap.get(d, args.model)}" for d, _ in DRUGS))

    global SHAP_MODEL
    SHAP_MODEL = {d: mmap.get(d, args.model) for d, _ in DRUGS}
    summary, _ = load(pathlib.Path(args.workbook))
    data = {d: load_drug(pathlib.Path(args.pred_dir), pathlib.Path(args.fold_dir),
                         d, mmap.get(d, args.model), b, summary, args.tolerance)
            for d, b in DRUGS}
    print(f"\n[OK] all {len(data)} drugs reproduce their reported row\n")
    out = pathlib.Path(args.out_dir)
    sheet_scatters(data, out)
    sheet_classdist(data, out)
    sheet_shap(pathlib.Path(args.pred_dir), out, top=args.shap_top)
    return 0


if __name__ == "__main__":
    sys.exit(main())
