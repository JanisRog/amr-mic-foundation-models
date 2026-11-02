#!/usr/bin/env python3
"""Extended Data Figure 3 (draft): how well one antimicrobial's MIC predicts resistance to another differs between centres.

The same AUROC as Figure 5a,b (the log2 MIC of the antimicrobial in the column separating the isolates resistant to
the antimicrobial in the row from the susceptible ones), computed within each collection site (M. tuberculosis) or
study (P. aeruginosa), for the between-class pairs whose pooled AUROC is at least 80%. A group is shown for a pair when
it has at least 20 resistant and at least 20 susceptible isolates, and a pair when at least three groups qualify.
One dot per site or study, the pooled AUROC as a black bar, the range as a line.

M. tuberculosis sites: the isolate IDs (site.NN); all sites used the same plate and protocol (CRyPTIC).
P. aeruginosa: the source studies, which also differ in how the MICs were measured.

Writes extended_figure3_values.tsv (every dot) beside the figure.
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
from openpyxl import load_workbook

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_figure3_alt import CODE  # noqa: E402
from make_figure5 import CLASSES, auroc, phenotypes, sheet  # noqa: E402
from make_panel_figures import SUBMISSION, fig_h, INK, INK2, ITALIC, ITALIC_BOLD, MUTED, ORG_COLOUR, SURFACE, header, save, tidy  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
ORG = {"mtb": "M. tuberculosis", "pseaer": "P. aeruginosa"}
MIN_EACH, MIN_GROUPS, MIN_POOLED = 20, 3, 0.80


def groups_of(wb, k):
    if k == "mtb":
        return lambda i: re.match(r"site\.([^.]+)", i).group(1)
    st = {r["run accession"]: r["study"] for r in sheet(wb, "PA_isolates")}
    return lambda i: st[i]


def pair_table(wb, k, rng):
    drugs, mic, res = phenotypes(wb, k)
    grp = groups_of(wb, k)
    inside = {(x, y) for c in CLASSES[k] for x in c for y in c if x != y}
    out, widths = [], []
    for da in drugs:
        for db in drugs:
            if da == db or (da, db) in inside:
                continue
            ids = [x for x in res[da] if x in mic[db]]
            pos = [mic[db][x] for x in ids if res[da][x]]
            neg = [mic[db][x] for x in ids if not res[da][x]]
            if len(pos) < 10 or not neg:
                continue
            pooled = auroc(pos, neg)
            if pooled < MIN_POOLED:
                continue
            per = {}
            for g in sorted({grp(x) for x in ids}):
                gi = [x for x in ids if grp(x) == g]
                p = np.array([mic[db][x] for x in gi if res[da][x]])
                n = np.array([mic[db][x] for x in gi if not res[da][x]])
                if len(p) >= MIN_EACH and len(n) >= MIN_EACH:
                    b = [auroc(rng.choice(p, len(p)), rng.choice(n, len(n))) for _ in range(200)]
                    lo, hi = np.percentile(b, [2.5, 97.5])
                    per[g] = (auroc(p, n), len(p), len(n))
                    widths.append(100 * (hi - lo))
            if len(per) >= MIN_GROUPS:
                out.append((da, db, pooled, per))
    return sorted(out, key=lambda t: -t[2]), widths


def panel(ax, rows, k, group_name):
    tidy(ax, x=True, y=False)
    col = ORG_COLOUR[k]
    for r, (da, db, pooled, per) in enumerate(rows):
        v = [100 * a for a, *_ in per.values()]
        ax.plot([min(v), max(v)], [r, r], color="#cfcdc5", lw=2.2, solid_capstyle="round", zorder=1)
        ax.plot(v, [r] * len(v), "o", ms=3.6, color=col, mec=SURFACE, mew=0.4, alpha=0.85, zorder=3)
        ax.plot([100 * pooled] * 2, [r - 0.32, r + 0.32], color=INK, lw=1.5, zorder=4)
        ax.annotate(f"{max(v) - min(v):.0f}", xy=(1.0, r), xycoords=("axes fraction", "data"), xytext=(5, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=5.8, color=INK2,
                    annotation_clip=False)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{CODE[da]}  from  {CODE[db]}" for da, db, *_ in rows], fontsize=6.4, color=INK)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.invert_yaxis()
    ax.set_xlim(40, 100)
    ax.tick_params(axis="x", labelsize=6.4)
    ax.set_xlabel("AUROC (%)", fontsize=6.8)
    ax.annotate("spread", xy=(1.0, 1.0), xycoords="axes fraction", xytext=(5, 3), textcoords="offset points",
                ha="left", va="bottom", fontsize=5.8, color=INK2)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--isolates", default=str(HOME / "AMRPredict_isolates_used.xlsx"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="extended_figure3_draft")
    ap.add_argument("--width", type=float, default=7.8)
    a = ap.parse_args()
    wb = load_workbook(a.isolates, read_only=True)
    rng = np.random.default_rng(3)
    data = {k: pair_table(wb, k, rng) for k in ("mtb", "pseaer")}
    n_m, n_p = len(data["mtb"][0]), len(data["pseaer"][0])
    fig = plt.figure(figsize=(a.width, fig_h(6.6)))
    gs = fig.add_gridspec(2, 1, height_ratios=[n_m, n_p + 3], hspace=0.38)
    axes = {}
    for r, k in enumerate(("mtb", "pseaer")):
        ax = fig.add_subplot(gs[r, 0])
        panel(ax, data[k][0], k, "site" if k == "mtb" else "study")
        axes[k] = ax
    top = header(fig, "Extended Data Figure 3 (draft) — how well one antimicrobial predicts another differs between centres",
                 "AUROC for telling isolates resistant to the first antimicrobial from susceptible ones by the log₂ MIC of the\n"
                 "second, between-class pairs with a pooled AUROC of at least 80%. One dot per collection site (a) or study (b)\n"
                 "with at least 20 resistant and 20 susceptible isolates; black bar, all isolates pooled; line, range; the number\n"
                 "at the right is the range in AUROC points. A single site's AUROC has a 95% bootstrap interval of about 5 to 12\n"
                 "points. The studies of P. aeruginosa also differ in how the MICs were measured.", pad_in=0.30)
    fig.subplots_adjust(top=top - 0.07, bottom=0.075, left=0.20, right=0.93)
    for k, lt, name in (("mtb", "a", "M. tuberculosis · collection sites"), ("pseaer", "b", "P. aeruginosa · studies")):
        pos = axes[k].get_position()
        fig.text(0.02, pos.y1 + 0.012, lt, fontsize=10, fontweight="bold", color=INK, ha="left", va="bottom")
        fig.text(pos.x0 + pos.width / 2, pos.y1 + 0.012, name, ha="center", va="bottom", fontsize=9.0,
                 fontproperties=ITALIC_BOLD, color=INK)
    out = pathlib.Path(a.out_dir)
    save(fig, out, a.stem)
    with open(out / "extended_figure3_values.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["organism", "resistance_to", "mic_of", "group", "auroc_group", "resistant", "susceptible", "auroc_pooled"])
        for k in ("mtb", "pseaer"):
            for da, db, pooled, per in data[k][0]:
                for g, (au, npos, nneg) in per.items():
                    w.writerow([ORG[k], da, db, g, round(100 * au, 1), npos, nneg, round(100 * pooled, 1)])
    for k in ("mtb", "pseaer"):
        sp = [100 * (max(a for a, *_ in per.values()) - min(a for a, *_ in per.values())) for *_, per in data[k][0]]
        print(f"[INFO] {ORG[k]}: {len(sp)} pairs; spread median {np.median(sp):.0f} (range {min(sp):.0f}-{max(sp):.0f}); "
              f"single-group 95% interval width median {np.median(data[k][1]):.1f} (IQR "
              f"{np.percentile(data[k][1], 25):.1f}-{np.percentile(data[k][1], 75):.1f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
