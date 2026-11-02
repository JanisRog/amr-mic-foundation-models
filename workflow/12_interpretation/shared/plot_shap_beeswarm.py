#!/usr/bin/env python3
"""Target-level SHAP beeswarm: what the model used, per isolate.

One row per target, one point per held-out isolate, positioned by that target's
contribution to that isolate's predicted MIC. A bar chart of mean |SHAP| would
compress each row to a single number and lose the thing worth seeing -- whether
a target pushes every isolate a little or a few isolates a great deal. On
isoniazid those are different mechanisms: katG has a bimodal row (variant or
not), while the promoter windows have a long tail.

Colour separates the gene body from the 300 bp upstream window, because they are
separate features and the distinction is the finding: for isoniazid the promoter
windows carry signal their gene bodies do not.

**The model is verified before it is explained.** The SHAP job writes predictions
beside the values; those are checked against the workbook's holdout row here, on
the same terms as every other figure. An attribution for a model that is not the
reported one would be worse than no attribution.

The swarm is deterministic: points are binned in x and offset symmetrically
within the bin, so the figure is identical on every redraw. Jitter from a random
number generator would move points between renders and invite reading structure
into what is really just noise in the offsets.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import AXIS, GRID, INK, INK2, MUTED, SURFACE, header, load, save, tidy
from plot_mic_scatter import CHECK, MODEL_LABEL, choose_aggregation

# One row per locus, with the gene body in its upper half and the 300 bp
# upstream window in its lower half. A row per *target* put katG in two places
# and separated a locus from its own promoter by six rows; the pairing is the
# comparison the figure is for.
# PA adds two kinds MTB has no equivalent of: acquired-gene sequence targets and
# a presence block of many columns. Both get a category so the figure can say
# which *kind* of feature carries a drug, which is the arm question.
REGION = {"gene": ("#6b3fa0", "gene body", -0.19),
          "upstream": ("#c3a5de", "300 bp upstream", +0.19),
          "acquired": ("#0079f2", "acquired gene, sequence", -0.19),
          "mutation": ("#cf2a5a", "point mutation / truncation", 0.0),
          "presence": ("#e8912b", "acquired gene, present", 0.0),
          "other": ("#8a8880", "catalogue feature", 0.0)}
BAND = "#f1f0ea"

# Separation checked in OKLab across deuteranopia, protanopia and tritanopia:
# every pair clears dE 8.6 or better. `gene` and `upstream` are exempt and
# deliberately the same hue -- they are one locus seen twice, separated by
# POSITION within the row (body above, upstream below), so identity is never
# carried by colour alone there.

# An AMRFinder presence column naming a substitution (gyrA_T83I) or a truncation
# (oprD_LOF) rather than a whole gene (sul1, blaVIM).
MUTATION_RE = re.compile(r"(?:_[A-Z]\d+[A-Z*]$)|(?:_LOF$)")


def jitter(n, half, seed):
    """Vertical scatter within a half-row, from a seeded generator.

    Binning and fanning out produced visible vertical stripes at 1,778 points --
    structure in the offsets that a reader could mistake for structure in the
    data. Uniform jitter reads as the cloud it is. Seeded per row, so the figure
    is identical on every redraw.
    """
    return np.random.default_rng(seed).uniform(-half, half, n)


def normalise(long):
    """Give every AMRFinder presence column its own row, and a real name.

    `shap_target_values.py` writes `locus='rule'` for the whole catalogue block,
    which is right for MTB -- one rule bit, one row -- and wrong for PA, where
    the block is a presence column per AMR family. Grouping by locus then
    collapsed eleven distinct features into a single row labelled with MTB's
    rule-bit name, so ciprofloxacin's gyrA_T83I and parC_S87L were being drawn
    on top of each other.

    Done here rather than in the SHAP job: it is a presentation decision, and
    redoing it upstream would mean re-running ten cluster jobs to rename a row.
    """
    long = long.copy()
    # ARDaP's call (the reported P. aeruginosa model for the antimicrobials it calls) is
    # one feature, like MTB's TB-Profiler call, and carries its own name.
    ard = long["target_id"].str.startswith("rule::ardap_call_")
    long.loc[ard, "locus"] = "ardap_call"
    long.loc[ard, "gene"] = "ARDaP call"
    m = (long["region"] == "other") & long["target_id"].str.startswith("rule::amr_")
    if not m.any():
        return long                      # MTB or ARDaP call: one catalogue feature, nothing to split
    name = long.loc[m, "target_id"].str.removeprefix("rule::amr_")
    long.loc[m, "locus"] = name
    long.loc[m, "gene"] = name
    long.loc[m, "region"] = np.where(name.str.contains(MUTATION_RE),
                                     "mutation", "presence")
    return long


def rank_loci(long, top=0):
    """Loci ordered by mean |SHAP|, with the metadata needed to label them."""
    rank = (long.assign(a=long["shap"].abs())
                .groupby("locus", as_index=False)["a"].mean()
                .sort_values("a", ascending=False))
    order = list(rank["locus"])
    if top:
        order = order[:top]
    return order, long.drop_duplicates("locus").set_index("locus")


def locus_label(meta, lo):
    reg = meta.loc[lo, "region"]
    g = str(meta.loc[lo, "gene"]).strip()
    if reg in ("mutation", "presence"):
        # already the feature's own name, set by normalise()
        return g or str(lo)
    if reg == "other":
        return "TB-Profiler rule bit" if lo == "rule" else (g or str(lo))
    if reg == "acquired":
        return (g or str(lo).removeprefix("amrseq_")) + "  · acquired"
    return f"{g} ({lo})" if g else str(lo)


def draw_beeswarm(ax, long, order, *, s=4.5, alpha=0.45, band=True):
    """The swarm itself, into a given axes.

    Shared by the single-drug figure and the all-drug sheet so the two cannot
    drift apart. Rows are loci; within a row the gene body sits above and the
    300 bp upstream window below.
    """
    tidy(ax, x=True, y=False)
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_visible(True)
        ax.spines[sp].set_color(AXIS)
        ax.spines[sp].set_linewidth(0.8)
    if band:
        for i in range(len(order)):
            if i % 2:
                ax.axhspan(i - 0.5, i + 0.5, color=BAND, lw=0, zorder=0)
    ax.axvline(0, color=AXIS, lw=0.9, ls=(0, (3, 3)), zorder=1)

    for i, lo in enumerate(order):
        rows = long[long["locus"] == lo]
        for reg, g in rows.groupby("region"):
            col, _, shift = REGION[reg]
            v = g["shap"].to_numpy(float)
            ax.scatter(v, i + shift + jitter(len(v), 0.155, seed=i * 7 + len(reg)),
                       s=s, color=col, lw=0, alpha=alpha, zorder=3)
    ax.set_yticks(range(len(order)))
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.invert_yaxis()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shap-long", required=True)
    ap.add_argument("--predictions", required=True,
                    help="written by the same job; used to verify the model")
    ap.add_argument("--workbook", default=str(pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs")) / "AMRPredict_panel_results.xlsx"))
    ap.add_argument("--drug", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--sr-boundary", type=float, required=True)
    ap.add_argument("--tolerance", type=float, default=5e-4)
    ap.add_argument("--top", type=int, default=0, help="0 shows every target")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="")
    args = ap.parse_args()

    pred = pd.read_csv(args.predictions, sep="\t")
    summary, _ = load(pathlib.Path(args.workbook))
    want = next((r for r in summary if r["drug"] == args.drug
                 and r["model"] == args.model and r["split"] == "holdout"), None)
    if want is None:
        raise SystemExit(f"[ERROR] no holdout row for {args.drug}/{args.model}")
    want = {k: (None if want.get(k) in ("", None) else float(want[k]))
            for k, _ in CHECK}
    agg, m = choose_aggregation(pred, args.sr_boundary, want, args.tolerance)
    config = pred["config"].iloc[0]

    long = pd.read_csv(args.shap_long, sep="\t")
    long = normalise(long).fillna({"gene": ""})
    n_iso = long["sample_id"].nunique()
    if n_iso != m["n"]:
        raise SystemExit(f"[ERROR] {n_iso:,} isolates in the SHAP file but "
                         f"{m['n']:,} in the predictions — different runs")

    order, meta = rank_loci(long, args.top)

    fig, ax = plt.subplots(figsize=(8.8, 0.60 * len(order) + 2.4))
    draw_beeswarm(ax, long, order)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([locus_label(meta, lo) for lo in order],
                       fontsize=8.5, color=INK)
    ax.set_xlabel("SHAP value  —  contribution to the predicted log₂ MIC, "
                  "summed over the target's four PCs")
    ax.tick_params(axis="x", labelsize=8)

    present = set(long.loc[long["locus"].isin(order), "region"])
    seen = [r for r in REGION if r in present]

    def legend_text(r):
        # After normalise() the PA presence block is already split into its own
        # categories, so "other" can only be MTB's single rule bit.
        if r == "other":
            ard = long["target_id"].str.startswith("rule::ardap_call_").any()
            return "ARDaP call" if ard else "TB-Profiler rule bit"
        return REGION[r][1]

    ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=6,
                              color=REGION[r][0], label=legend_text(r))
                       for r in seen],
              loc="lower right", fontsize=8, handletextpad=0.6,
              labelspacing=0.6, borderpad=0.7)

    top = header(fig, f"{args.drug.capitalize()} — what the model used",
                 f"{MODEL_LABEL.get(args.model, args.model)}, "
                 f"{n_iso:,} held-out isolates, one point each. One row per "
                 f"locus, gene body above and 300 bp upstream window below,\n"
                 f"ordered by mean |SHAP|. A region's contribution is the sum of "
                 f"its four principal components — exact, since SHAP is additive "
                 f"in the\nfeatures. Config {config} ({agg}).")
    fig.subplots_adjust(top=top, bottom=0.95 / (0.60 * len(order) + 2.4),
                        left=0.215, right=0.985)

    stem = args.stem or f"shap_{args.drug}_{args.model}"
    save(fig, pathlib.Path(args.out_dir), stem)
    return 0


if __name__ == "__main__":
    sys.exit(main())
