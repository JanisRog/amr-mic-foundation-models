#!/usr/bin/env python3
"""Holdout: the reported model against its baseline, one row per drug.

Both panels answer to the same declared objective -- cost = 2*vMER + MER -- so
the table leads with it and carries the components beside it, because a cost
alone cannot say whether a loss is missed resistance or over-calling.

Two things this script will not do, both of which would flatter the result:

  * pick the model on the holdout. MTB declared evo2 in advance; PA picks per
    drug on DEVELOPMENT cost. Either way the choice is made before these numbers
    are looked at.
  * quietly choose the baseline convention. ARDaP reports Sensitive /
    Intermediate / Resistant, and folding Intermediate to R or to S moves the
    number; `I=R` is shown because it is the one that flatters the comparator.

Ladder metrics are blank for the baselines: a two-valued S/R call has no position
on the MIC ladder, and printing a zero there would read as a score rather than an
absence.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

MODELS = ("evo2", "carbon", "dnabert")
PANELS = {
    "mtb": ["isoniazid", "rifampicin", "ethambutol", "ethionamide", "kanamycin",
            "amikacin", "levofloxacin", "moxifloxacin", "bedaquiline", "delamanid"],
    "pseaer": ["ciprofloxacin", "levofloxacin", "tobramycin", "imipenem",
               "meropenem", "ceftazidime", "cefepime", "aztreonam", "amikacin",
               "piperacillin-tazobactam"],
}
# MTB's declared comparator. PA's comes from baseline_rows.tsv, not the workbook.
MTB_BASELINE = "TB-Profiler 6.3.0"
# Too few resistant isolates for a verdict; carried in the table with their
# numbers but kept out of the tally, which is what the panel table in the MTB
# protocol does. Counting them would let two underpowered drugs set the headline.
UNDERPOWERED = {"bedaquiline", "delamanid"}
# Below this, a per-drug difference is inside what a resplit moves -- measured on
# the PA comparator in section 12, where ARDaP shifted 0.145 on a cohort that lost
# one isolate. Reported as a second line, never as the headline: the sign count is
# the declared accounting and this is the caveat on it.
NOISE = 0.01
# A difference that rounds to 0.0000 at the precision printed is a tie, not a
# loss. The fluoroquinolones are the case: the model reproduces the catalogue
# exactly, because the whole signal is gyrA and the rule bit already carries it.
TIE = 5e-5


def load_summary(xlsx: pathlib.Path):
    import openpyxl
    wb = openpyxl.load_workbook(xlsx, data_only=True)
    rows = list(wb["summary"].iter_rows(values_only=True))
    head = [str(c) for c in rows[0]]
    return [dict(zip(head, r)) for r in rows[1:] if r[0] is not None]


def get(summary, drug, model, split, field):
    for r in summary:
        if r["drug"] == drug and r["model"] == model and r["split"] == split:
            v = r.get(field)
            return None if v in ("", None) else float(v)
    return None


def cost(vmer, mer):
    return None if vmer is None or mer is None else 2 * vmer + mer


def reported_model(summary, drug, declared):
    """Declared model if there is one; otherwise lowest DEVELOPMENT cost."""
    if declared:
        return declared
    best, bestc = None, float("inf")
    for m in MODELS:
        c = cost(get(summary, drug, m, "development_cv", "vMER"),
                 get(summary, drug, m, "development_cv", "MER"))
        if c is not None and c < bestc:
            best, bestc = m, c
    return best


def ardap_row(base_dir: pathlib.Path, drug: str, label="ARDaP (I=R)"):
    f = base_dir / f"cv_{drug}" / "baseline_rows.tsv"
    if not f.exists():
        return None
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["split"] == "holdout" and r["label"] == label:
            return {k: (float(r[k]) if r.get(k) not in ("", None) else None)
                    for k in ("vMER", "MER", "SR_f1", "cost", "SR_AUROC_raw")}
    return None


def fmt(v, w=7, p=4):
    return f"{'—':>{w}}" if v is None else f"{v:>{w}.{p}f}"


def cell(model_v, base_v, p=3, w=16, pct=False):
    """value (Δ vs baseline), with the bracket omitted when there is no baseline.

    As percentages the difference is in **percentage points**, not a relative
    change: 3.5% against 5.0% is written −1.5, and calling that "−30%" would be a
    different and much larger-sounding claim about the same two numbers.
    """
    if model_v is None:
        return f"{'—':>{w}}"
    k, unit, dp = (100.0, "%", 1) if pct else (1.0, "", p)
    if base_v is None:
        return f"{model_v * k:.{dp}f}{unit}".rjust(w)
    d = (model_v - base_v) * k
    # A signed "-0.0" reads as a difference the precision cannot actually show.
    txt = f"{0:.{dp}f}" if abs(d) < 0.5 * 10 ** -dp else f"{d:+.{dp}f}"
    return f"{model_v * k:.{dp}f}{unit} ({txt})".rjust(w)


def deltas_table(summary, bdir, args) -> int:
    """One-tier accuracy, vMER, MER, F1 — each with its delta to the baseline.

    Sign convention is always model minus baseline, so the *good* direction is
    not the same for every column: negative is better for vMER and MER, positive
    for F1. Stated under the table rather than encoded in the sign, because
    flipping some columns to "improvement" would make two tables in this document
    disagree about what a minus means.

    One-tier accuracy carries no bracket. The baselines are S/R callers and have
    no position on the MIC ladder at all; a delta there would need a number that
    does not exist.
    """
    pct = args.percent
    dlab = "Δ in percentage points" if pct else "Δ"
    print(f"\n{args.title or args.panel}  —  holdout, model value "
          f"({dlab} vs baseline)")
    print(f"{'drug':<25}{'model':<9}{'one-tier':>10}{'vMER':>18}"
          f"{'MER':>18}{'F1':>18}")
    print("-" * 98)

    for drug in PANELS[args.panel]:
        m = reported_model(summary, drug, args.declared_model)
        if m is None:
            continue
        g = lambda f: get(summary, drug, m, "holdout", f)
        if bdir is not None:
            b = ardap_row(bdir, drug, args.baseline_label)
            bv, bm, bf = (b["vMER"], b["MER"], b["SR_f1"]) if b else (None,) * 3
        else:
            bv = get(summary, drug, MTB_BASELINE, "holdout", "vMER")
            bm = get(summary, drug, MTB_BASELINE, "holdout", "MER")
            bf = get(summary, drug, MTB_BASELINE, "holdout", "SR_f1")
        flag = " †" if drug in UNDERPOWERED else "  "
        ot = g("OT")
        ot_txt = (f"{'—':>10}" if ot is None else
                  (f"{ot * 100:>9.1f}%" if pct else f"{ot:>10.3f}"))
        print(f"{drug:<23}{flag}{m:<9}" + ot_txt
              + cell(g("vMER"), bv, 3, 18, pct)
              + cell(g("MER"), bm, 3, 18, pct)
              + cell(g("SR_f1"), bf, 3, 18, pct))

    print("-" * 98)
    print(f"{dlab} is model − baseline: NEGATIVE is better for vMER and MER, "
          f"POSITIVE for F1.")
    print("One-tier accuracy has no Δ — the baselines are S/R callers with no "
          "position on the MIC ladder.")
    print("A blank bracket means the baseline does not cover that drug.")
    if UNDERPOWERED & set(PANELS[args.panel]):
        print(f"† {', '.join(sorted(UNDERPOWERED & set(PANELS[args.panel])))}: "
              f"too few resistant isolates for a verdict.")
    print()
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workbook", required=True)
    ap.add_argument("--panel", choices=sorted(PANELS), required=True)
    ap.add_argument("--declared-model", default="",
                    help="MTB declares evo2 in advance; PA leaves this empty and "
                         "the model is chosen per drug on development cost")
    ap.add_argument("--baseline-dir", default="",
                    help="tree holding cv_<drug>/baseline_rows.tsv (PA/ARDaP)")
    ap.add_argument("--baseline-label", default="ARDaP (I=R)")
    ap.add_argument("--title", default="")
    ap.add_argument("--format", choices=("full", "deltas"), default="full",
                    help="deltas: one-tier accuracy, vMER, MER and F1, each with "
                         "the difference from the baseline in brackets")
    ap.add_argument("--percent", action="store_true",
                    help="report as percentages; the bracketed difference is "
                         "then in percentage points")
    args = ap.parse_args()

    summary = load_summary(pathlib.Path(args.workbook))
    bdir = pathlib.Path(args.baseline_dir) if args.baseline_dir else None

    if args.format == "deltas":
        return deltas_table(summary, bdir, args)

    print(f"\n{args.title or args.panel}  —  holdout, model vs baseline")
    print(f"{'':<25}{'':<9}{'reported model':^39}{'baseline':^34}")
    print(f"{'drug':<25}{'model':<9}{'cost':>8}{'vMER':>8}{'MER':>8}"
          f"{'F1':>8}{'bal_OT':>8}{'OT':>7}   "
          f"{'cost':>8}{'vMER':>8}{'MER':>8}{'F1':>7}   {'Δcost':>8}")
    print("-" * 127)

    wins = losses = exact_ties = covered = big = 0
    for drug in PANELS[args.panel]:
        m = reported_model(summary, drug, args.declared_model)
        if m is None:
            print(f"{drug:<25}{'—':<9}  (no model rows)")
            continue
        g = lambda f: get(summary, drug, m, "holdout", f)
        mc = cost(g("vMER"), g("MER"))

        if bdir is not None:
            b = ardap_row(bdir, drug, args.baseline_label)
            bc, bv, bm, bf = ((b["cost"], b["vMER"], b["MER"], b["SR_f1"])
                              if b else (None, None, None, None))
        else:
            bv = get(summary, drug, MTB_BASELINE, "holdout", "vMER")
            bm = get(summary, drug, MTB_BASELINE, "holdout", "MER")
            bf = get(summary, drug, MTB_BASELINE, "holdout", "SR_f1")
            bc = cost(bv, bm)

        d = None if (mc is None or bc is None) else mc - bc
        flag = " †" if drug in UNDERPOWERED else "  "
        if d is not None and drug not in UNDERPOWERED:
            covered += 1
            if abs(d) < TIE:
                exact_ties += 1
            elif d < 0:
                wins += 1
            else:
                losses += 1
            big += abs(d) > NOISE
        print(f"{drug:<23}{flag}{m:<9}{fmt(mc,8)}{fmt(g('vMER'),8)}"
              f"{fmt(g('MER'),8)}{fmt(g('SR_f1'),8)}{fmt(g('bal_OT'),8)}"
              f"{fmt(g('OT'),7,3)}   {fmt(bc,8)}{fmt(bv,8)}{fmt(bm,8)}"
              f"{fmt(bf,7,3)}   " + ("       —" if d is None else f"{d:>+8.4f}"))

    print("-" * 127)
    print(f"Lower cost on {wins} of {covered} — {losses} to the baseline"
          + (f", {exact_ties} reproducing it exactly" if exact_ties else ""))
    print(f"  of which {big} differ by more than {NOISE:g} — the rest are inside "
          f"what a resplit moves (§12).")
    if any(d in PANELS[args.panel] for d in UNDERPOWERED):
        print(f"  † {', '.join(sorted(UNDERPOWERED & set(PANELS[args.panel])))}"
              f": too few resistant isolates for a verdict; shown, not counted.")
    print("Ladder metrics (bal_OT, OT) are model-only: an S/R call has no "
          "position on the MIC ladder.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
