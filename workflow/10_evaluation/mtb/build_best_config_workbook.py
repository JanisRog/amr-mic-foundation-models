#!/usr/bin/env python3
"""Every arm's selected configuration, every metric, with the incumbent beside it.

One row per (drug, arm) carrying the configuration that 1-SE selection actually
chose -- not a sweep summary, not a best-of over metrics. The selection ran once
per run on `cost` with `bal_OT` as tie-break, and every number here comes from
that single configuration's fold file, so a drug cannot look good because one
metric found a different winner.

TB-Profiler 6.3.0 sits in the same tables as another row rather than in a
separate sheet, because the comparison is the point. Two things about that row
must not be read as if it were an arm:

  **It has no ordinal metrics.** A catalogue is a rules engine: it returns
  susceptible or resistant and has no position on the MIC ladder at all. Those
  cells are empty because the quantity does not exist for it, not because it
  scored zero. Writing 0 there would manufacture a win.

  **It is pooled, not per-fold.** Its vMER and MER are computed once over the
  whole development set; the arms' are means over 15 fold fits. The isolates are
  the same, so the levels are comparable, but there is no fold-level spread and
  therefore no paired test and no standard error against it. Arm-vs-arm
  comparisons are paired and testable; arm-vs-catalogue comparisons are a
  difference of two point estimates.

Bedaquiline and delamanid are reported and excluded from panel aggregates. The
catalogue misses 52% and 81% of their resistant isolates respectively, giving
costs of 1.05 and 1.61 against 0.08-0.35 everywhere else, and a mean over ten
drugs would be a statement about those two. Panel rows use the median and the
win count, and name the drugs they rest on.

Development folds only; the holdout is scored once and is not what this compares.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics
import sys

# Three naming conventions, because run_inh_redo.sbatch and run_rif_dev.sbatch
# take no target-set variable and so write <model>_<tag>, while
# run_panel_dev.sbatch writes <model>_<targetset>_<tag>. A suffix missing from
# this table raises nothing -- the arm is simply absent and the n shrinks.
# All four arms embed the TB-Profiler catalogue loci. The arms differ in two
# binary choices, not in whether the catalogue is used:
#
#   target set   `both`         = curated mechanism loci | TB-Profiler loci
#                `gwasfiltered` = `both` | unannotated GWAS loci   (a superset)
#   extra column rule bit (TB-Profiler's own call) or per-mutation columns
#
# `gwasfiltered` ADDS to `both`; it never replaces it. The only GWAS bucket
# excluded is `other` -- loci the catalogue annotates as another drug's
# mechanism, which is what imports co-resistance. Unannotated loci are kept
# deliberately, since dropping them would remove the only route to a mechanism
# the catalogue does not already know.
ARMS = {
    "baseline": (["both_final1k", "tbrule1k", "final1k"],
                 "TB-Profiler + curated loci, + rule bit (reported arm)"),
    "arm75":    (["both_arm75", "arm75"],
                 "TB-Profiler + curated loci, + per-mutation columns"),
    "gwfrule":  (["gwasfiltered_gwfrule1k"],
                 "the above plus unannotated GWAS loci, + rule bit"),
    "gwfsnp":   (["gwasfiltered_gwfsnp1k"],
                 "the above plus unannotated GWAS loci, + mutation columns"),
}
MODELS = ("evo2", "carbon", "dnabert")
# Declared primary before any run (§9). Every reported number uses it.
# Switching the reported model was considered and rejected once already,
# at a gap of 0.0013; the best-model columns here are a diagnostic and
# must not become a per-drug selection.
PRIMARY = "evo2"

# A difference smaller than the arm's own fold-to-fold standard error is not a
# verdict. The band is per drug and per arm, not a flat number: the standard
# error of cost runs from 0.0027 (rifampicin) to 0.0133 (amikacin), so a single
# threshold would call isoniazid's 2.9-sigma separation a tie while treating
# amikacin's noise as a result. One SE is the same rule selection already uses.
TIE_SE_MULTIPLE = 1.0

# metric -> lower is better. The SRt_* family is excluded: it comes from a
# fitted S/R threshold that §11 measured and rejected, and putting it beside
# breakpoint-based metrics would mix two decision rules in one table.
METRICS = {"cost": True, "vMER": True, "MER": True,
           "SR_f1": False, "SR_acc": False, "SR_AUROC_raw": False,
           "exact": False, "exact_bal": False,
           "OT": False, "bal_OT": False, "bal_OT_common": False,
           "MAE": True, "RMSE": True, "R2": False, "spearman": False}

# The five the catalogue can be scored on. Everything else is ladder geometry.
BASELINE_METRICS = ("cost", "vMER", "MER", "SR_f1", "SR_acc")

# Fold-directory name -> the summaries directory holding that drug's baseline.
# The codes are the rule-bit file codes, not the first three letters of the drug
# (moxifloxacin is mxf), and isoniazid's directory is cv_final_summaries.
SUMMARIES = {
    "isoniazid_final": "cv_final_summaries", "rifampicin": "cv_rif_summaries",
    "ethambutol": "cv_emb_summaries", "ethionamide": "cv_eth_summaries",
    "kanamycin": "cv_kan_summaries", "amikacin": "cv_ami_summaries",
    "levofloxacin": "cv_lev_summaries", "moxifloxacin": "cv_mxf_summaries",
    "bedaquiline": "cv_bdq_summaries", "delamanid": "cv_dlm_summaries",
}
ORDER = list(SUMMARIES)
NO_VERDICT = {"bedaquiline", "delamanid"}


def pretty(drug: str) -> str:
    return drug[:-len("_final")] if drug.endswith("_final") else drug


def load_folds(run_dir: pathlib.Path, aggregation: str):
    """The selected configuration's per-fold metrics, keyed by fold."""
    f = run_dir / "folds.tsv"
    if not f.exists():
        return None
    out = {}
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["aggregation"] != aggregation:
            continue
        m = {k: float(r[k]) for k in METRICS
             if k in r and r[k] not in ("", "NA", "nan", None)}
        # Recomputed rather than read, so the declared objective is the one
        # being reported even if a run's own cost column used other weights.
        m["cost"] = 2 * float(r["vMER"]) + float(r["MER"])
        out[r["fold"]] = m
    return out or None


def cost_se(folds: dict) -> float | None:
    """Standard error of cost across this configuration's folds."""
    v = [f["cost"] for f in folds.values() if "cost" in f]
    if len(v) < 3:
        return None
    return statistics.stdev(v) / len(v) ** 0.5


def find(root, drug, model, arm):
    for suffix in ARMS[arm][0]:
        d = root / f"cv_{drug}" / "runs" / f"{model}_{suffix}"
        if (d / "folds.tsv").exists():
            return d
    return None


def selected_tag(cfg_dir: pathlib.Path, drug, model, arm):
    """The configuration 1-SE selection chose, for the provenance column."""
    import json
    for suffix in ARMS[arm][0]:
        p = cfg_dir / f"{drug}__{model}__{suffix}.json"
        if p.exists():
            try:
                return json.load(open(p))[0]["tag"]
            except Exception:
                return ""
    return ""


def load_baseline(summaries_root: pathlib.Path, drug: str, split: str):
    """Every baseline label for one drug on one split: TB-Profiler and the
    catalogue AF variants where they were computed."""
    p = summaries_root / SUMMARIES[drug] / "baseline_rows.tsv"
    if not p.exists():
        return {}
    out = {}
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["split"] != split:
            continue
        out[r["label"]] = {k: float(r[k]) for k in BASELINE_METRICS if r.get(k)}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True,
                    help="slim_arms directory (selected configs only)")
    ap.add_argument("--summaries", required=True,
                    help="mtb_tbprofiler directory holding cv_*_summaries/")
    ap.add_argument("--aggregation", default="median")
    ap.add_argument("--split", default="development_cv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font, PatternFill
        from openpyxl.utils import get_column_letter
    except ImportError:
        raise SystemExit("[ERROR] needs openpyxl (not present in the cluster env "
                         "-- run this locally)")

    root = pathlib.Path(args.root)
    summ = pathlib.Path(args.summaries)

    present = {d.name[3:] for d in root.glob("cv_*") if d.is_dir()}
    unknown = present - set(SUMMARIES)
    if unknown:
        raise SystemExit(f"[ERROR] fold directories with no summaries mapping: "
                         f"{sorted(unknown)} -- add them to SUMMARIES")
    drugs = [d for d in ORDER if d in present or (summ / SUMMARIES[d] /
                                                  "baseline_rows.tsv").exists()]

    # (drug, model, arm) -> {metric: mean over folds}, and the fold count.
    cell, nfolds, tags, ses = {}, {}, {}, {}
    for drug in drugs:
        for model in MODELS:
            for arm in ARMS:
                d = find(root, drug, model, arm)
                if d is None:
                    continue
                v = load_folds(d, args.aggregation)
                if not v:
                    continue
                cell[(drug, model, arm)] = {
                    k: statistics.mean(f[k] for f in v.values() if k in f)
                    for k in METRICS
                    if any(k in f for f in v.values())}
                nfolds[(drug, model, arm)] = len(v)
                ses[(drug, model, arm)] = cost_se(v)
                tags[(drug, model, arm)] = selected_tag(
                    root / "configs", drug, model, arm)

    base = {d: load_baseline(summ, d, args.split) for d in drugs}
    # The catalogue's holdout cost, carried as a reference column only. It is a
    # different split and must never be differenced against a development
    # number; it is here because the catalogue itself moves between the two
    # (kanamycin 0.354 -> 0.293, levofloxacin 0.145 -> 0.194), so a development
    # verdict does not imply the holdout verdict.
    base_hold = {d: load_baseline(summ, d, "holdout") for d in drugs}
    missing_base = [d for d in drugs if "TB-Profiler 6.3.0" not in base[d]]
    if missing_base:
        print(f"[WARN] no TB-Profiler row on split '{args.split}' for "
              f"{[pretty(d) for d in missing_base]}")

    def common_models(drug):
        """The encoders present in every arm this drug has.

        Pooling whatever exists makes rows within a drug incomparable:
        moxifloxacin's arm75 lost its carbon and dnabert fits, so a
        three-model mean sat beside a one-model mean and showed the gwfsnp
        control differing by 0.013 when at the per-model level it is exact to
        six decimals. Restricting to the common set costs resolution and keeps
        the comparison a comparison."""
        sets = [{m for m in MODELS if (drug, m, a) in cell} for a in ARMS]
        sets = [x for x in sets if x]
        if not sets:
            return []
        keep = set.intersection(*sets)
        return [m for m in MODELS if m in keep]

    def pooled_se(drug, arm):
        """The primary encoder's own standard error across its folds, matching
        the basis `pooled` reports."""
        return ses.get((drug, PRIMARY, arm))

    def pooled(drug, arm):
        """The declared primary encoder's numbers -- the reported basis.

        This was a mean over the three encoders and that was the wrong default:
        it made isoniazid's baseline cost read 0.0923 where the manuscript
        reports evo2's 0.0866, a quantity that appears in no table anywhere.
        A cross-encoder mean is nobody's reported figure. evo2 is present in
        every arm of every drug, so the primary basis is also automatically
        comparable across arms, which the common-model restriction had to work
        around. The mean is kept, on its own sheet, clearly named."""
        return cell.get((drug, PRIMARY, arm))

    def encoder_mean(drug, arm):
        """Mean over the encoders present in every arm of this drug. Separate
        sheet, because mixing a one-encoder mean with three-encoder means makes
        rows incomparable -- see common_models."""
        vs = [cell[(drug, m, arm)] for m in common_models(drug)
              if (drug, m, arm) in cell]
        if not vs:
            return None
        return {k: statistics.mean(v[k] for v in vs if k in v)
                for k in METRICS if any(k in v for v in vs)}

    arms = list(ARMS)
    wb = Workbook()
    HEAD = Font(bold=True, color="FFFFFF", size=10)
    FILL = PatternFill("solid", fgColor="2A4D69")
    DRUGFILL = PatternFill("solid", fgColor="DCE6F1")
    BASEFILL = PatternFill("solid", fgColor="FDE9D9")
    BOLD = Font(bold=True, size=10)
    MUTED = Font(italic=True, size=9, color="808080")
    NUM = "0.0000"

    def header(ws, cols, widths):
        ws.append(cols)
        for i, c in enumerate(cols, 1):
            cl = ws.cell(row=1, column=i)
            cl.font, cl.fill = HEAD, FILL
            cl.alignment = Alignment(horizontal="center", wrap_text=True)
            ws.column_dimensions[get_column_letter(i)].width = widths[i - 1]
        ws.freeze_panes = ws.cell(row=2, column=1)

    def write_metrics(ws, row, start, vals):
        for j, k in enumerate(METRICS):
            c = ws.cell(row=row, column=start + j)
            if vals and k in vals:
                c.value = round(vals[k], 6)
                c.number_format = NUM
            else:
                c.value = None

    mcols = list(METRICS)

    # ---- by_drug: the sheet to read first -----------------------------------
    ws = wb.active
    ws.title = "by_drug"
    header(ws, ["drug", "arm", "what it is", "basis", "folds"]
           + mcols, [17, 11, 44, 16, 7] + [10] * len(mcols))
    r = 2
    for drug in drugs:
        for arm in arms:
            v = pooled(drug, arm)
            have = sum(1 for m in MODELS if (drug, m, arm) in cell)
            ws.cell(row=r, column=1, value=pretty(drug) if arm == arms[0] else "")
            ws.cell(row=r, column=2, value=arm)
            ws.cell(row=r, column=3, value=ARMS[arm][1])
            ws.cell(row=r, column=4,
                    value=(f"{PRIMARY} (primary)" if have else None))
            ws.cell(row=r, column=5,
                    value=nfolds.get((drug, PRIMARY, arm)))
            if arm == "baseline":
                ws.cell(row=r, column=1).font = BOLD
                for c in range(1, 6 + len(mcols)):
                    ws.cell(row=r, column=c).fill = DRUGFILL
            write_metrics(ws, r, 6, v)
            if v is None:
                ws.cell(row=r, column=3).font = MUTED
                ws.cell(row=r, column=3, value=ARMS[arm][1] + "  — not run yet")
            r += 1
        b = base[drug].get("TB-Profiler 6.3.0")
        ws.cell(row=r, column=2, value="TB-Profiler")
        ws.cell(row=r, column=3,
                value="catalogue 6.3.0 — S/R only, pooled over development")
        for c in range(1, 6 + len(mcols)):
            ws.cell(row=r, column=c).fill = BASEFILL
        write_metrics(ws, r, 6, b)
        if drug in NO_VERDICT:
            ws.cell(row=r, column=3,
                    value="catalogue 6.3.0 — misses most resistant isolates; "
                          "excluded from panel aggregates")
        r += 1

    # ---- vs_tbprofiler: the five metrics the comparison can actually use ----
    ws2 = wb.create_sheet("vs_tbprofiler")
    dcols = [f"Δ{k}" for k in BASELINE_METRICS]
    header(ws2, ["drug", "arm"] + list(BASELINE_METRICS) + ["TBP cost"] + dcols
           + ["SE(cost)", "Δcost / SE", "verdict (Δcost)",
              "TBP cost\n(holdout, ref only)"],
           [17, 11] + [10] * 5 + [10] + [10] * 5 + [10, 11, 16, 14])
    r = 2
    wins = {a: [0, 0, 0] for a in arms}
    for drug in drugs:
        b = base[drug].get("TB-Profiler 6.3.0")
        for arm in arms:
            v = pooled(drug, arm)
            ws2.cell(row=r, column=1, value=pretty(drug) if arm == arms[0] else "")
            ws2.cell(row=r, column=2, value=arm)
            for j, k in enumerate(BASELINE_METRICS):
                c = ws2.cell(row=r, column=3 + j)
                if v and k in v:
                    c.value, c.number_format = round(v[k], 6), NUM
            if b:
                c = ws2.cell(row=r, column=8, value=round(b["cost"], 6))
                c.number_format = NUM
                for j, k in enumerate(BASELINE_METRICS):
                    if v and k in v and k in b:
                        c = ws2.cell(row=r, column=9 + j,
                                     value=round(v[k] - b[k], 6))
                        c.number_format = "+0.0000;-0.0000"
                if v:
                    d_cost = v["cost"] - b["cost"]
                    se = pooled_se(drug, arm)
                    if se:
                        c = ws2.cell(row=r, column=14, value=round(se, 6))
                        c.number_format = NUM
                        c = ws2.cell(row=r, column=15,
                                     value=round(d_cost / se, 2))
                        c.number_format = "+0.00;-0.00"
                    tie = se is not None and abs(d_cost) < TIE_SE_MULTIPLE * se
                    verdict = ("tie (<1 SE)" if tie else
                               "model" if d_cost < 0 else "catalogue")
                    ws2.cell(row=r, column=16, value=verdict)
                    if drug not in NO_VERDICT:
                        wins[arm][0 if d_cost < 0 else 1] += 1
                        if tie:
                            wins[arm][2] += 1
                if drug in base_hold and "TB-Profiler 6.3.0" in base_hold[drug]:
                    c = ws2.cell(row=r, column=17,
                                 value=round(base_hold[drug]
                                             ["TB-Profiler 6.3.0"]["cost"], 6))
                    c.number_format, c.font = NUM, MUTED
            r += 1
        r += 1
    ws2.cell(row=r + 1, column=1,
             value=f"Costs and SE are the declared primary encoder ({PRIMARY}); "
                   f"encoder_mean has the pooled view and the two differ "
                   f"(isoniazid baseline 0.0866 against 0.0923). "
                   "Selection ran on these folds, so development cost flatters "
                   "the arms and not the catalogue, which had no selection step. "
                   "The holdout column is a different split: reference only, "
                   "never to be differenced against the columns left of "
                   "it.").font = MUTED

    # ---- panel: medians and win counts, never a mean over ten drugs --------
    ws3 = wb.create_sheet("panel")
    header(ws3, ["arm", "what it is", "drugs"] + [f"median {k}" for k in mcols],
           [11, 38, 8] + [12] * len(mcols))
    scored = [d for d in drugs if d not in NO_VERDICT
              and all(pooled(d, a) for a in arms)]
    r = 2
    for arm in arms:
        ws3.cell(row=r, column=1, value=arm)
        ws3.cell(row=r, column=2, value=ARMS[arm][1])
        ws3.cell(row=r, column=3, value=len(scored))
        for j, k in enumerate(mcols):
            vals = [pooled(d, arm)[k] for d in scored if k in pooled(d, arm)]
            if vals:
                c = ws3.cell(row=r, column=4 + j,
                             value=round(statistics.median(vals), 6))
                c.number_format = NUM
        r += 1
    ws3.cell(row=r, column=1, value="TB-Profiler")
    ws3.cell(row=r, column=2, value="catalogue 6.3.0 (S/R only)")
    ws3.cell(row=r, column=3, value=len(scored))
    for j, k in enumerate(mcols):
        vals = [base[d]["TB-Profiler 6.3.0"][k] for d in scored
                if k in base[d].get("TB-Profiler 6.3.0", {})]
        if vals:
            c = ws3.cell(row=r, column=4 + j,
                         value=round(statistics.median(vals), 6))
            c.number_format = NUM
    for c in range(1, 4 + len(mcols)):
        ws3.cell(row=r, column=c).fill = BASEFILL
    r += 3
    ws3.cell(row=r, column=1, value="Drugs beating TB-Profiler on cost "
                                    f"(of {len(scored)}, bedaquiline and "
                                    "delamanid excluded)").font = BOLD
    r += 1
    for arm in arms:
        ws3.cell(row=r, column=1, value=arm)
        n = wins[arm][0] + wins[arm][1]
        ws3.cell(row=r, column=2,
                 value=f"{wins[arm][0]} of {n}"
                       + (f"  ({wins[arm][2]} of {n} within 1 SE either way)"
                          if wins[arm][2] else ""))
        r += 1

    ws3.cell(row=r + 1, column=1,
             value="Median, not mean: costs span 0.07 to 1.61 across the panel "
                   "and a mean would report the two drugs the catalogue cannot "
                   "call.").font = MUTED

    # ---- by_drug_model: nothing pooled -------------------------------------
    ws4 = wb.create_sheet("by_drug_model")
    header(ws4, ["drug", "arm", "model", "cost rank", "selected config",
                 "folds"] + mcols, [17, 11, 10, 10, 30, 7] + [10] * len(mcols))
    r = 2
    for drug in drugs:
        for arm in arms:
            have = [m for m in MODELS if (drug, m, arm) in cell]
            if not have:
                continue
            order = sorted(have, key=lambda m: cell[(drug, m, arm)]["cost"])
            for model in have:
                ws4.cell(row=r, column=1, value=pretty(drug))
                ws4.cell(row=r, column=2, value=arm)
                ws4.cell(row=r, column=3,
                         value=model + ("  (primary)" if model == PRIMARY else ""))
                rank = order.index(model) + 1
                c = ws4.cell(row=r, column=4,
                             value=f"{rank} of {len(have)}"
                                   + ("  ← best" if rank == 1 else ""))
                if rank == 1:
                    c.font = BOLD
                if model == PRIMARY:
                    ws4.cell(row=r, column=3).font = BOLD
                ws4.cell(row=r, column=5, value=tags.get((drug, model, arm), ""))
                ws4.cell(row=r, column=6, value=nfolds[(drug, model, arm)])
                write_metrics(ws4, r, 7, cell[(drug, model, arm)])
                r += 1
            r += 1

    # ---- encoder_mean: the pooled view, named so it cannot be mistaken ----
    ws8 = wb.create_sheet("encoder_mean")
    header(ws8, ["drug", "arm", "encoders\n(used/have)"] + mcols,
           [17, 11, 13] + [10] * len(mcols))
    r = 2
    for drug in drugs:
        for arm in arms:
            v = encoder_mean(drug, arm)
            have = sum(1 for m in MODELS if (drug, m, arm) in cell)
            used = sum(1 for m in common_models(drug) if (drug, m, arm) in cell)
            ws8.cell(row=r, column=1, value=pretty(drug) if arm == arms[0] else "")
            ws8.cell(row=r, column=2, value=arm)
            ws8.cell(row=r, column=3, value=(f"{used}/{have}" if have else None))
            write_metrics(ws8, r, 4, v)
            r += 1
        r += 1
    ws8.cell(row=r + 1, column=1,
             value="Mean over the encoders present in every arm of that drug. "
                   "This is NOT the reported basis -- no table anywhere reports "
                   "a cross-encoder mean. by_drug carries the declared primary "
                   f"({PRIMARY}). Isoniazid's baseline reads 0.0923 here and "
                   "0.0866 there for exactly this reason.").font = MUTED

    # ---- best_model: a diagnostic on the declared primary, not a selection --
    ws7 = wb.create_sheet("best_model")
    header(ws7, ["drug", "arm"] + [f"{m} cost" for m in MODELS]
           + ["lowest cost", "spread", "gap: primary − best", "SE(cost)",
              "spread / SE", "separable?"],
           [17, 11] + [11] * 3 + [12, 10, 17, 10, 11, 12])
    r = 2
    tally = {m: 0 for m in MODELS}
    gaps, complete = [], 0
    for drug in drugs:
        for arm in arms:
            have = [m for m in MODELS if (drug, m, arm) in cell]
            if len(have) < 2:
                continue
            v = {m: cell[(drug, m, arm)]["cost"] for m in have}
            best = min(v, key=v.get)
            spread = max(v.values()) - min(v.values())
            se = pooled_se(drug, arm)
            ws7.cell(row=r, column=1, value=pretty(drug) if arm == arms[0] else "")
            ws7.cell(row=r, column=2, value=arm)
            for j, m in enumerate(MODELS):
                if m in v:
                    c = ws7.cell(row=r, column=3 + j, value=round(v[m], 6))
                    c.number_format = NUM
                    if m == best:
                        c.font = BOLD
            ws7.cell(row=r, column=6,
                     value=best + ("  (primary)" if best == PRIMARY else ""))
            c = ws7.cell(row=r, column=7, value=round(spread, 6))
            c.number_format = NUM
            if PRIMARY in v:
                c = ws7.cell(row=r, column=8, value=round(v[PRIMARY] - v[best], 6))
                c.number_format = NUM
            if se:
                c = ws7.cell(row=r, column=9, value=round(se, 6))
                c.number_format = NUM
                c = ws7.cell(row=r, column=10, value=round(spread / se, 2))
                c.number_format = "0.00"
                ws7.cell(row=r, column=11,
                         value="yes" if spread > se else "no (< 1 SE)")
            if len(have) == len(MODELS):
                tally[best] += 1
                complete += 1
                if PRIMARY in v:
                    gaps.append(v[PRIMARY] - v[best])
            r += 1
        r += 1
    r += 1
    ws7.cell(row=r, column=1,
             value=f"Lowest cost over the {complete} (drug, arm) cells with all "
                   f"three encoders").font = BOLD
    r += 1
    for m in MODELS:
        ws7.cell(row=r, column=1,
                 value=m + ("  (declared primary)" if m == PRIMARY else ""))
        ws7.cell(row=r, column=2, value=f"{tally[m]} of {complete}")
        r += 1
    if gaps:
        r += 1
        ws7.cell(row=r, column=1, value="primary − best, median").font = BOLD
        ws7.cell(row=r, column=2, value=round(statistics.median(gaps), 6)
                 ).number_format = NUM
        ws7.cell(row=r + 1, column=1, value="primary − best, worst")
        ws7.cell(row=r + 1, column=2, value=round(max(gaps), 6)
                 ).number_format = NUM
        ws7.cell(row=r + 3, column=1,
                 value="evo2 was declared primary before any run and every "
                       "reported number uses it. These columns are a diagnostic: "
                       "picking the lowest-cost encoder per drug after seeing "
                       "the results would be a selection step, and none of the "
                       "margins here would survive it.").font = MUTED

    # ---- catalogue_variants: only where they were computed -----------------
    ws5 = wb.create_sheet("catalogue_variants")
    header(ws5, ["drug", "label"] + list(BASELINE_METRICS),
           [17, 20] + [10] * 5)
    r = 2
    for drug in drugs:
        for label, v in base[drug].items():
            if label == "TB-Profiler 6.3.0":
                continue
            ws5.cell(row=r, column=1, value=pretty(drug))
            ws5.cell(row=r, column=2, value=label)
            for j, k in enumerate(BASELINE_METRICS):
                if k in v:
                    c = ws5.cell(row=r, column=3 + j, value=round(v[k], 6))
                    c.number_format = NUM
            r += 1
    if r == 2:
        ws5.cell(row=2, column=1,
                 value="No allele-fraction variants were computed for any drug.")

    # ---- notes --------------------------------------------------------------
    ws6 = wb.create_sheet("notes")
    ws6.column_dimensions["A"].width = 24
    ws6.column_dimensions["B"].width = 108
    rows = [
        ("What this is",
         "The configuration 1-SE selection chose for each run -- one per (drug, "
         "model, arm) -- with every metric that configuration's fold file "
         f"carries. Aggregation '{args.aggregation}', split '{args.split}'."),
        ("Selection",
         "Minimise cost = 2*vMER + MER; ties within one standard error broken on "
         "bal_OT. Declared before the runs. Every metric here comes from that one "
         "configuration, so no drug can look good because a different metric "
         "picked a different winner."),
        ("TB-Profiler: no ladder",
         "Blank ordinal cells (exact, OT, bal_OT, MAE, RMSE, R2, spearman) are "
         "blank because a rules engine has no position on the MIC ladder. The "
         "quantity does not exist for it. A zero there would be a manufactured "
         "win and must not be filled in."),
        ("TB-Profiler: pooled",
         "Its vMER and MER are computed once over the whole development set; the "
         "arms' are means over 15 fold fits on the same isolates. Levels are "
         "comparable; there is no fold-level spread, so no paired test and no "
         "standard error against it."),
        ("Development folds flatter the arms",
         "1-SE selection chose each arm's configuration on these same folds; "
         "TB-Profiler had no selection step. The development comparison is "
         "therefore biased toward the arms, and a loss here is a stronger result "
         "than a loss on the holdout would be. The pre-registered comparison is "
         "the holdout one and is not in this workbook."),
        ("Ties",
         "A difference smaller than the arm's own standard error of cost across "
         "its folds is reported as a tie, the same one-SE rule selection uses. "
         "The band is per drug: SE(cost) runs from 0.0027 for rifampicin to "
         "0.0133 for amikacin, so a flat threshold would call isoniazid's "
         "2.9-sigma separation a tie and amikacin's noise a result. The Δcost/SE "
         "column is there to be read directly rather than taken on trust."),
        ("Which comparisons are testable",
         "Arm vs arm shares drugs, folds and isolates -- paired Wilcoxon applies "
         "(see compare_arms.py). Arm vs catalogue is a difference of two point "
         "estimates and carries no p-value."),
        ("Bedaquiline, delamanid",
         "Reported, excluded from panel aggregates. The catalogue misses 51.9% "
         "and 80.6% of resistant isolates, costing 1.05 and 1.61 against "
         "0.08-0.35 elsewhere. Neither drug carries a verdict."),
        ("Arms all use the catalogue",
         "Every arm embeds the TB-Profiler loci. The arms differ in whether the "
         "unannotated GWAS loci are added on top (gwasfiltered is a superset of "
         "both, never a replacement) and in whether the extra column is "
         "TB-Profiler's own call or per-mutation columns. Feature counts show "
         "it: rifampicin 33 -> 257, ethambutol 57 -> 289, levofloxacin 17 -> 49, "
         "moxifloxacin 17 -> 17."),
        ("Which basis each sheet uses",
         f"by_drug, vs_tbprofiler and panel report the declared primary encoder "
         f"({PRIMARY}) -- the reported basis. encoder_mean carries the mean over "
         f"encoders; best_model and by_drug_model carry every encoder "
         f"separately. The same cell therefore has two legitimate values and "
         f"they must not be quoted interchangeably: isoniazid's baseline cost is "
         f"0.0866 on the primary and 0.0923 as a three-encoder mean. Only the "
         f"first appears in any reported table."),
        ("Why not the encoder mean",
         "It was the original default here and it was wrong. A mean across "
         "encoders is a quantity no table in the manuscript reports, and it "
         "also had to be restricted to encoders common to every arm -- "
         "moxifloxacin's arm75 lost its carbon and dnabert fits, so pooling "
         "whatever existed put a one-encoder mean beside three-encoder means "
         "and made its gwfsnp control appear to differ by 0.013 when per "
         "encoder it is exact to six decimals. The primary is present in every "
         "arm of every drug, so that problem does not arise on the reported "
         "basis. encoder_mean keeps the used/have column for the same reason."),
        ("Which model is best",
         "evo2 was declared primary in advance (§9) and holds: it has the "
         "lowest cost in 23 of the 35 (drug, arm) cells with all three "
         "encoders, and the median gap between it and the per-cell winner is "
         "0.0000. The best_model sheet reports this as a diagnostic. It is not "
         "a selection -- choosing the encoder per drug after seeing the results "
         "would be exactly the step the evaluation contract forbids, and §13 "
         "already considered and rejected switching the reported model at a gap "
         "of 0.0013."),
        ("Most encoder differences are not separable",
         "The spread/SE column says when a 'best' model is real. For the drugs "
         "with a verdict the across-encoder spread is 0.0025-0.0119 while "
         "SE(cost) is 0.0026-0.0143, so the three encoders are mostly within "
         "one standard error of each other. Where evo2's margin is largest is "
         "informative in the other direction: it wins 7 of 8 arm75 cells and 8 "
         "of 9 gwfsnp cells but only 3 of 9 gwfrule cells. With the rule bit "
         "present the encoders converge, which is §6's point -- a dominant "
         "shared feature is the condition under which different encoders agree, "
         "not evidence that encoders are interchangeable."),
        ("Controls",
         "Moxifloxacin has no unannotated GWAS loci, so gwasfiltered equals both "
         "for it and two pairs must reproduce exactly: gwfrule against baseline "
         "and gwfsnp against arm75. Both are exact to six decimals on all 15 "
         "folds for every model present. A non-zero difference would mean the "
         "cache or config path differed, not the target set."),
        ("SRt_* excluded",
         "The fitted-threshold family is left out: §11 measured and rejected it, "
         "and showing it beside breakpoint metrics would mix two decision rules."),
        ("What this cannot settle",
         "Development folds only. §10 records that a co-resistance marker helps "
         "on development folds and gives the gain back across lineages, so a "
         "wider target set scoring better here cannot be distinguished from one "
         "that imported more co-resistance. Leave-one-lineage-out is the "
         "discriminating test and has not been run."),
        ("Coverage",
         "Arms present for " + ", ".join(pretty(d) for d in drugs
                                         if all(pooled(d, a) for a in arms))
         + ". Absent: " + (", ".join(pretty(d) for d in drugs
                                     if not all(pooled(d, a) for a in arms))
                           or "none") + "."),
    ]
    for i, (k, v) in enumerate(rows, 1):
        ws6.cell(row=i, column=1, value=k).font = BOLD
        c = ws6.cell(row=i, column=2, value=v)
        c.alignment = Alignment(wrap_text=True, vertical="top")
        ws6.row_dimensions[i].height = 14 * (1 + len(v) // 95)

    wb.save(args.out)
    print(f"[OK] {args.out}")
    print(f"     {len(drugs)} drugs, {len(arms)} arms, {len(METRICS)} metrics; "
          f"{len(cell)} (drug, model, arm) cells")
    print(f"     panel aggregates rest on {len(scored)} drugs: "
          f"{', '.join(pretty(d) for d in scored)}")
    for arm in arms:
        print(f"     {arm:<10} beats TB-Profiler on cost in "
              f"{wins[arm][0]} of {wins[arm][0] + wins[arm][1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
