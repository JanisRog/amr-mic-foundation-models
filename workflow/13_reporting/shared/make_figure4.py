#!/usr/bin/env python3
"""Figure 4 (draft): what the reported models rely on, from TreeSHAP on the held-out isolates.

Left column, M. tuberculosis. Right column, P. aeruginosa.
  top     beeswarm of the highest-ranked features of one antimicrobial (isoniazid; ciprofloxacin):
          one point per held-out isolate, gene body above and 300 bp upstream window below
  bottom  per antimicrobial, the share of the mean absolute SHAP value carried by each kind of input
          M. tuberculosis    TB-Profiler's call / loci the TB-Profiler database annotates / loci
                             from the association study
          P. aeruginosa      ARDaP's call / chromosomal loci / acquired genes / AMRFinderPlus columns
                             (levofloxacin and aztreonam, which ARDaP does not call)

SHAP values are those of the reported models. Shares are computed here from the SHAP long tables;
they are written to figure4_values.tsv beside the figure.
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
from matplotlib.patches import Patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_panel_figures import (SUBMISSION, fig_h, GREY, INK, INK2, ITALIC, ITALIC_BOLD, MUTED, ORG_COLOUR, ORG_PALE, PURPLE, SURFACE, header,  # noqa: E402
                                load, pick, save, tidy)
from make_figure3_alt import CODE, spread_labels  # noqa: E402
from make_provisional_figures import NAME, p_shap, p_share  # noqa: E402
from plot_metric_bars import PANELS, SHORT  # noqa: E402
from plot_shap_beeswarm import REGION, draw_beeswarm, locus_label, normalise, rank_loci  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
PA_KINDS = (("ARDaP call", "#8a8880"), ("chromosomal loci", "#6b3fa0"), ("acquired genes", "#0079f2"),
            ("AMRFinderPlus columns", "#cf2a5a"))


def pa_beeswarm(ax, shap_dir, drug="ciprofloxacin", top=10):
    long = normalise(pd.read_csv(pathlib.Path(shap_dir) / f"{drug}_evo2_shap_long.tsv", sep="\t")
                     .fillna({"gene": ""}))
    order, meta = rank_loci(long, top=top)
    sub = long[long["locus"].isin(order)]
    draw_beeswarm(ax, sub, order, s=2.4, alpha=0.45)
    ax.set_yticklabels([locus_label(meta, lo).split(" (")[0] for lo in order], fontsize=6.6)
    ax.tick_params(axis="x", labelsize=6.8)
    ax.set_xlabel("SHAP value  (log₂ MIC)", fontsize=7.2)
    ax.set_title(f"{drug.capitalize()}: the {top} highest-ranked features", loc="left", fontsize=8.4, color=INK,
                 pad=6, fontweight="medium")
    present = set(sub["region"])
    label = {"other": "ARDaP call"}
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=4.2, color=REGION[r][0],
                              label=label.get(r, REGION[r][1])) for r in REGION if r in present],
              loc="upper center", bbox_to_anchor=(0.40, -0.13), ncol=2, fontsize=6.2, handletextpad=0.1,
              columnspacing=0.9)


DAGGER = {"levofloxacin", "aztreonam"}          # P. aeruginosa: AMRFinderPlus columns in place of ARDaP's call


def share_panel(ax, mtb_dir, pa_dir, which=("mtb", "pseaer"), legend=True):
    """Per antimicrobial, the share of the mean |SHAP| carried by the catalogue input and by the embeddings.

    Catalogue input: TB-Profiler's call (M. tuberculosis), ARDaP's call (P. aeruginosa), or the
    AMRFinderPlus presence columns where ARDaP gives no call (levofloxacin, aztreonam). Everything
    else is the embeddings of the target loci. Read from the per-target mean |SHAP| tables, where the
    catalogue input is the feature whose target id starts with `rule`.
    """
    tidy(ax, x=False, y=True)
    pos, ticks, labels, groups, rows = 0.0, [], [], [], []
    for k, d_, name in [t for t in (("mtb", mtb_dir, "M. tuberculosis"), ("pseaer", pa_dir, "P. aeruginosa"))
                        if t[0] in which]:
        P, start = PANELS[k], pos
        for d in P["drugs"]:
            ts = pd.read_csv(pathlib.Path(d_) / f"{d}_evo2_target_summary.tsv", sep="\t")
            tot = ts["mean_abs_shap"].sum()
            rule = ts.loc[ts["target_id"].str.startswith("rule"), "mean_abs_shap"].sum() / tot
            al = 1.0           # the attribution rests on every held-out isolate: bedaquiline and delamanid are shown in full
            ax.bar(pos, rule, 0.70, color=GREY, edgecolor=SURFACE, linewidth=0.6, alpha=al, zorder=3)
            ax.bar(pos, 1 - rule, 0.70, bottom=rule, color=PURPLE, edgecolor=SURFACE, linewidth=0.6, alpha=al,
                   zorder=3)
            if rule > 0.07:
                ax.annotate(f"{rule:.0%}", (pos, rule / 2), ha="center", va="center", fontsize=5.4, color=SURFACE,
                            zorder=4)
            if 1 - rule > 0.07:
                ax.annotate(f"{1 - rule:.0%}", (pos, rule + (1 - rule) / 2), ha="center", va="center",
                            fontsize=5.4, color=SURFACE, zorder=4)
            ticks.append(pos)
            labels.append(SHORT.get(d, d.capitalize()) + (" †" if (k == "pseaer" and d in DAGGER) else ""))
            rows.append((name, d, rule))
            pos += 1.0
        groups.append((name, (start + pos - 1) / 2))
        pos += 1.1
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=6.4, color=INK, rotation=38, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.7, pos - 1.1 - 0.3)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"], fontsize=6.8)
    ax.set_ylabel("share of mean |SHAP|", fontsize=7.2)
    ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    for name, x in (groups if len(which) > 1 else []):
        ax.annotate(name, xy=(x, 1.0), xytext=(0, 5), textcoords="offset points", ha="center", va="bottom",
                    fontsize=9.0, fontproperties=ITALIC_BOLD, color=INK, annotation_clip=False)
    if legend:
        ax.legend(handles=[Patch(facecolor=GREY, edgecolor=SURFACE, label="catalogue input (TB-Profiler or ARDaP call; "
                                                                      "† AMRFinderPlus columns)"),
                       Patch(facecolor=PURPLE, edgecolor=SURFACE, label="embeddings of the target loci")],
                  loc="upper center", bbox_to_anchor=(0.5, -0.36), ncol=2, fontsize=6.6, handlelength=1.1,
                  columnspacing=1.6, borderaxespad=0)
    return rows


def quality_panel(ax, a, rows):
    """The call's share of the attribution against the tool's own quality, one point per antimicrobial.

    Quality is the F1 of the tool alone on the development folds, so it does not depend on the held-out
    isolates the SHAP values are computed on. Bedaquiline and delamanid are included: the development
    folds hold 54 and 62 resistant isolates, against 12 and 15 held out. Spearman correlations are
    printed and drawn in the panel.
    """
    from scipy.stats import spearmanr
    books = {"M. tuberculosis": ("mtb", load(pathlib.Path(a.mtb_workbook))[0]),
             "P. aeruginosa": ("pseaer", load(pathlib.Path(a.pa_workbook))[0])}
    pts = []
    for org, d, rule in rows:
        k, book = books[org]
        P = PANELS[k]
        f1 = pick(book, d, P["baseline"], "development_cv", "SR_f1")
        if f1 is None:                       # no tool call for this antimicrobial
            continue
        pts.append((org, d, rule, f1))
    tidy(ax, x=True, y=True)
    labels = []
    for org, d, rule, f1 in pts:
        mk = "o" if org.startswith("M") else "D"          # organism shapes and colours as in Figure 3
        ax.plot(100 * f1, rule, mk, ms=5.4 if mk == "o" else 4.6, mew=0,
                color=ORG_COLOUR["mtb"] if org.startswith("M") else ORG_COLOUR["pseaer"],
                zorder=4, alpha=0.95)
        labels.append((100 * f1, rule, CODE[d], 0.7, 0.7))
    rho = {}
    for name, sel in (("both organisms", pts), ("M. tuberculosis", [t for t in pts if t[0].startswith("M")]),
                      ("P. aeruginosa", [t for t in pts if t[0].startswith("P")])):
        c, pv = spearmanr([t[3] for t in sel], [t[2] for t in sel])
        rho[name] = (c, pv, len(sel))
        print(f"[INFO] call share vs tool F1 (development), {name}: rho = {c:+.2f}, P = {pv:.3f}, n = {len(sel)}")
    ax.set_xlim(22, 100)
    ax.set_ylim(-0.02, 0.70)
    ax.set_yticks([0, 0.2, 0.4, 0.6])
    ax.set_yticklabels(["0", "20%", "40%", "60%"], fontsize=6.8)
    ax.tick_params(axis="x", labelsize=6.8)
    ax.set_xlabel("F1 of the tool alone (%, development folds)", fontsize=7.0)
    ax.set_ylabel("share of attribution on the catalogue call", fontsize=7.0)
    ax.set_title("Share on the call against the tool's quality", loc="left", fontsize=8.0, color=INK, pad=14,
                 fontweight="semibold")
    txt = (f"Spearman ρ = {rho['both organisms'][0]:+.2f}, P = {rho['both organisms'][1]:.3f} "
           f"(n = {rho['both organisms'][2]})")
    # above the plot, on the line of the organism names in d
    ax.annotate(txt, xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 2), textcoords="offset points",
                ha="left", va="bottom", fontsize=6.4, color=INK2, annotation_clip=False)
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=5.2, color=ORG_COLOUR["mtb"], mew=0,
                              label="M. tuberculosis"),
                       Line2D([], [], marker="D", ls="none", ms=4.6, color=ORG_COLOUR["pseaer"], mew=0,
                              label="P. aeruginosa")],
              loc="lower right", fontsize=6.2, handletextpad=0.2, borderaxespad=0.4,
              prop={"size": 6.2})
    for t in ax.get_legend().get_texts():            # species names in italics
        t.set_fontproperties(ITALIC)
        t.set_fontsize(6.2)
    return labels


def override_panel(ax, a):
    """Per antimicrobial: the isolates on which the framework's S/R call differed from the tool's, as a share of
    the held-out isolates, split into overrides that matched the measured phenotype and those that did not.

    The tool's call per isolate is read from the sign of its SHAP value (positive = resistant call); this
    reproduces the tools' held-out error rates exactly for M. tuberculosis and within 2.6 points for
    P. aeruginosa (intermediate ARDaP calls). Levofloxacin and aztreonam (P. aeruginosa) have no call.
    """
    from make_figure2 import MTB_BOUND, PA_BOUND
    tidy(ax, x=False, y=True)
    pos, ticks, labels, groups, out = 0.0, [], [], [], []
    for k, d_, bounds, name in (("mtb", a.mtb_shap, MTB_BOUND, "M. tuberculosis"),
                                ("pseaer", a.pa_shap, PA_BOUND, "P. aeruginosa")):
        start = pos
        for d, b in bounds.items():
            sh = pd.read_csv(pathlib.Path(d_) / f"{d}_evo2_shap_long.tsv", sep="\t")
            call = sh[sh["target_id"].str.startswith("rule")]
            if call.empty or call["target_id"].nunique() > 1:
                continue                                    # no single tool call
            pr = pd.read_csv(pathlib.Path(d_) / f"{d}_evo2_predictions.tsv", sep="\t")
            if "median" in set(pr["aggregation"]):          # the reported aggregation
                pr = pr[pr["aggregation"] == "median"]
            pr = pr.drop_duplicates("sample_id")
            df = call.set_index("sample_id")[["shap"]].join(pr.set_index("sample_id")[["log2_mic", "raw"]], how="inner")
            c_r, t_r, p_r = df["shap"] > 0, df["log2_mic"] > b, df["raw"] > b
            ov = p_r != c_r
            right, wrong = int((ov & (p_r == t_r)).sum()), int((ov & (p_r != t_r)).sum())
            n = len(df)
            ax.bar(pos, 100 * right / n, 0.70, color=ORG_COLOUR[k], edgecolor=SURFACE, linewidth=0.5, zorder=3)
            ax.bar(pos, 100 * wrong / n, 0.70, bottom=100 * right / n, color=ORG_PALE[k], edgecolor=SURFACE,
                   linewidth=0.5, zorder=3)
            ax.annotate(f"{right}/{right + wrong}", (pos, 100 * (right + wrong) / n), xytext=(0, 2),
                        textcoords="offset points", ha="center", va="bottom", fontsize=4.8, color=INK2, rotation=90)
            ticks.append(pos)
            labels.append(CODE[d])
            out.append((name, d, n, right, wrong))
            pos += 1.0
        groups.append((name, (start + pos - 1) / 2))
        pos += 1.0
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=6.0, color=INK, rotation=90)
    ax.set_xlim(-0.7, pos - 1.3)
    ax.set_ylim(0, 44)
    ax.set_yticks([0, 10, 20, 30, 40])
    ax.set_yticklabels(["0", "10%", "20%", "30%", "40%"], fontsize=6.8)
    ax.set_ylabel("isolates where the framework\noverrode the tool's call", fontsize=7.0, linespacing=1.2)
    ax.axhline(0, color=INK2, lw=0.7, zorder=4)
    for name, x in groups:
        ax.annotate(name, xy=(x, 1.0), xycoords=("data", "axes fraction"), xytext=(0, 2), textcoords="offset points",
                    ha="center", va="bottom", fontsize=6.4, fontproperties=ITALIC, color=INK2)
    ax.set_title("Overrides of the tool's call", loc="left", fontsize=8.0, color=INK, pad=14, fontweight="semibold")
    ax.legend(handles=[Patch(facecolor=INK2, label="matched the measured phenotype"),
                       Patch(facecolor="#d6d4cc", label="did not")], loc="upper left", fontsize=5.8,
              frameon=False, handlelength=1.0, borderaxespad=0.2, bbox_to_anchor=(0.0, 0.90))
    for org in ("M. tuberculosis", "P. aeruginosa"):
        rr = [t for t in out if t[0] == org]
        n, ri, wr = sum(t[2] for t in rr), sum(t[3] for t in rr), sum(t[4] for t in rr)
        print(f"[INFO] overrides {org}: {ri + wr} of {n} isolates ({100 * (ri + wr) / n:.1f}%), "
              f"matching the phenotype {ri} ({100 * ri / max(ri + wr, 1):.0f}%)")
    return out


def shares_only(a, scatter=False) -> int:
    """Alternative: two panels, the attribution in M. tuberculosis (a) and in P. aeruginosa (b).

    With scatter=True a third panel (c) shows how far the call is relied on against the tool's quality.
    """
    fig = plt.figure(figsize=(a.width, fig_h(6.9) if scatter else 3.95))
    gs = fig.add_gridspec(2 if scatter else 1, 2, wspace=0.10, hspace=0.95, height_ratios=[1.0, 1.05] if scatter else [1.0])
    ax1 = fig.add_subplot(gs[0, 0])
    rows = share_panel(ax1, a.mtb_shap, a.pa_shap, which=("mtb",), legend=False)
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    rows += share_panel(ax2, a.mtb_shap, a.pa_shap, which=("pseaer",), legend=False)
    ax3, labels3 = None, None
    if scatter:
        # c: how far the call is relied on, against its quality; d: where the framework overrode it
        gs2 = gs[1, :].subgridspec(1, 2, width_ratios=[0.85, 1.0], wspace=0.30)
        ax_c = fig.add_subplot(gs2[0, 1])
        override_panel(ax_c, a)
        ax3 = fig.add_subplot(gs2[0, 0])
        labels3 = quality_panel(ax3, a, rows)
    ax2.set_ylabel("")
    plt.setp(ax2.get_yticklabels(), visible=False)
    for ax, letter, name in ((ax1, "a", "M. tuberculosis"), (ax2, "b", "P. aeruginosa")):
        ax.set_title(name, loc="center", fontsize=9.5, fontproperties=ITALIC_BOLD, color=INK, pad=8)
        ax.annotate(letter, xy=(-0.03 if letter == "b" else -0.13, 1.10), xycoords="axes fraction", fontsize=10,
                    fontweight="bold", color=INK, ha="right", va="bottom", annotation_clip=False)
    top = header(fig, "Figure 4, alternative (draft) — where the attribution sits",
                 "TreeSHAP on the held-out isolates: share of the mean absolute SHAP value carried by the catalogue "
                 "input and by the embeddings of the target loci.\n† levofloxacin and "
                 "aztreonam, which ARDaP does not call, use the AMRFinderPlus presence columns."
                 + ("\nc: share of the attribution on the call against the tool's F1 on the development folds. "
                    "d: held-out isolates on which the framework's S/R call differed\nfrom the tool's (dark: matched the "
                    "measured phenotype; numbers: matched/all overrides)." if scatter else ""), pad_in=0.28)
    fig.subplots_adjust(top=top - 0.045, bottom=0.075 if scatter else 0.265, left=0.105, right=0.975)
    leg_y = ax1.get_position().y0 - 0.125 if scatter else 0.0
    fig.legend(handles=[Patch(facecolor=GREY, edgecolor=SURFACE, label="catalogue input (TB-Profiler or ARDaP call; "
                                                                     "† AMRFinderPlus columns)"),
                        Patch(facecolor=PURPLE, edgecolor=SURFACE, label="embeddings of the target loci")],
               loc="lower center" if not scatter else "upper center", ncol=2, fontsize=6.8, handlelength=1.1,
               columnspacing=1.6, bbox_to_anchor=(0.5, 0.0 if not scatter else leg_y + 0.03))
    if scatter:
        p3, pc = ax3.get_position(), ax_c.get_position()    # c and d: same top and bottom
        y0, y1 = max(p3.y0, pc.y0), min(p3.y1, pc.y1)
        ax3.set_position([p3.x0, y0, p3.width, y1 - y0])
        ax_c.set_position([pc.x0, y0, pc.width, y1 - y0])
        spread_labels(fig, ax3, labels3)
        pos = ax3.get_position()
        for axx, lt in ((ax3, "c"), (ax_c, "d")):
            axx.annotate(lt, xy=(-0.22 if lt == "c" else -0.20, 1.10), xycoords="axes fraction", fontsize=10,
                         fontweight="bold", color=INK, ha="right", va="bottom", annotation_clip=False)
    out = pathlib.Path(a.out_dir)
    save(fig, out, a.stem or ("figure4_alt_quality_draft" if scatter else "figure4_alt_draft"))
    for name, d, rule in rows:
        print(f"[INFO] {name:16s} {d:24s} catalogue input {rule:.0%}  embeddings {1 - rule:.0%}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-shap", default=str(HOME / "mtb_shap_gwfrule1k/shap_gwfrule1k"))
    ap.add_argument("--mtb-earlier", default=str(HOME),
                    help="target summaries of the catalogue-derived model, to tell the two sources apart")
    ap.add_argument("--pa-shap", default=str(HOME / "pa_reported_20261021/predictions_strict"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stem", default="")
    ap.add_argument("--variant", choices=["full", "shares", "shares_quality"], default="full",
                    help="shares: only the attribution panels, one per organism; shares_quality: plus a panel "
                         "relating the call's share to the tool's quality")
    ap.add_argument("--mtb-workbook", default=str(HOME / "AMRPredict_panel_results_wider.xlsx"))
    ap.add_argument("--pa-workbook", default=str(HOME / "AMRPredict_pseaer_reported_results.xlsx"))
    ap.add_argument("--width", type=float, default=7.8)
    a = ap.parse_args()

    if a.variant in ("shares", "shares_quality"):
        return shares_only(a, scatter=a.variant == "shares_quality")
    fig = plt.figure(figsize=(a.width, 9.2))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.62, height_ratios=[0.78, 1.55])
    ax_a = fig.add_subplot(gs[0, :])
    rows = share_panel(ax_a, a.mtb_shap, a.pa_shap)
    ax_b = fig.add_subplot(gs[1, 0]); p_shap(ax_b, a)
    ax_c = fig.add_subplot(gs[1, 1]); pa_beeswarm(ax_c, a.pa_shap)
    for ax, letter, dx, dy in ((ax_a, "a", -0.07, 1.20), (ax_b, "b", -0.36, 1.08), (ax_c, "c", -0.30, 1.08)):
        ax.annotate(letter, xy=(dx, dy), xycoords="axes fraction", fontsize=10, fontweight="bold", color=INK,
                    ha="right", va="bottom", annotation_clip=False)
    top = header(fig, "Figure 4 (draft) — what the reported models rely on",
                 "TreeSHAP on the held-out isolates. a, share of the mean absolute SHAP value carried by the catalogue "
                 "input and by the embeddings of the target loci\n(† levofloxacin and "
                 "aztreonam, which ARDaP does not call, use the AMRFinderPlus presence columns). b, c, one point per "
                 "isolate; in each\nrow the gene body is above and the 300 bp upstream window below.", pad_in=0.30)
    fig.subplots_adjust(top=top - 0.045, bottom=0.085, left=0.20, right=0.965)
    for ax, title in ((ax_b, "M. tuberculosis"), (ax_c, "P. aeruginosa")):
        pos = ax.get_position()
        fig.text((pos.x0 + pos.x1) / 2, pos.y1 + 0.052, title, ha="center", va="bottom", fontsize=9.5,
                 fontproperties=ITALIC_BOLD, color=INK)
    out = pathlib.Path(a.out_dir)
    save(fig, out, a.stem or "figure4_draft")
    with open(out / "figure4_values.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["organism", "antimicrobial", "catalogue_input_share", "embeddings_share"])
        for name, d, rule in rows:
            w.writerow([name, d, round(rule, 4), round(1 - rule, 4)])
    for name, d, rule in rows:
        print(f"[INFO] {name:16s} {d:24s} catalogue input {rule:.0%}  embeddings {1 - rule:.0%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
