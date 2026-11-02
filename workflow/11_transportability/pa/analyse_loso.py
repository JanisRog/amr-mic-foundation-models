#!/usr/bin/env python3
"""Leave-one-study-out against matched random cross-validation.

The question is not "how well does the model do on study S" -- a collection that
is 99% resistant is easy for any constant predictor and a collection with two
susceptible isolates is unscoreable. The question is how much performance is lost
**on the same isolates** when the model can no longer see that collection during
training. That difference is the transport cost; the level is study difficulty.

So every row compares two fits over one study's isolates:

  LOSO    the reported config, refit with that whole collection held out
  CV      the same config under ordinary 5-fold x 3-repeat cross-validation

Under CV each isolate has exactly one out-of-fold prediction per repeat, so the
metric is computed within each repeat and averaged across the three. Pooling the
three predictions per isolate first would be an ensemble and would flatter CV,
which is the arm we are trying not to flatter.

Fold -1 is reported separately and never counted. It holds out every collection
too small to score, pooled, which is a real fit but not a per-study answer.
"""

from __future__ import annotations

import argparse
import csv
import pathlib

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.environ.get("MTB_SCRIPTS") or str(pathlib.Path(__file__).resolve().parent.parent / "mtb_tbprofiler"))
from plot_mic_scatter import metrics

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from pick_reported_model import DESIGNATED  # noqa: E402

BOUND = {"ciprofloxacin": -1.0, "levofloxacin": 1.0, "tobramycin": 1.0,
         "imipenem": 2.0, "meropenem": 3.0, "ceftazidime": 3.0, "cefepime": 3.0,
         "aztreonam": 4.0, "amikacin": 4.0, "piperacillin-tazobactam": 4.0}


def score(df, boundary):
    if df.empty:
        return None
    m = metrics(df["log2_mic"].to_numpy(float), df["pred"].to_numpy(float), boundary)
    m["cost"] = 2 * m["vMER"] + m["MER"]
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--loso-dir", required=True)
    ap.add_argument("--metadata-root", default="",
                    help="holds <drug>/metadata.tsv (cluster layout)")
    ap.add_argument("--metadata-flat", default="",
                    help="one TSV with 'Accession number' and 'Study' for every isolate; use instead of "
                         "--metadata-root when the per-drug files are not on this machine")
    ap.add_argument("--aggregation", default="median")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    root = pathlib.Path(args.loso_dir)
    rows, pooled = [], []

    for drug in sorted(BOUND):
        lo = sorted(root.glob(f"{drug}_*_loso_predictions.tsv"))
        cv = sorted(root.glob(f"{drug}_*_cv_predictions.tsv"))
        if not lo or not cv:
            print(f"[skip] {drug}: no leave-one-study-out run "
                  f"(single-collection drug)")
            continue
        # Several models can have a run on disk (carbon for amikacin and meropenem
        # predates 2026-10-13). Read the reported one, never whichever sorts first.
        want = [f for f in lo if f.name == f"{drug}_{DESIGNATED}_loso_predictions.tsv"]
        if want:
            lo = want
            cv = [root / f"{drug}_{DESIGNATED}_cv_predictions.tsv"]
        model = lo[0].name[len(drug) + 1:-len("_loso_predictions.tsv")]
        b = BOUND[drug]

        if args.metadata_flat:
            study = {r["Accession number"]: (r.get("Study") or "unknown") for r in
                     csv.DictReader(open(args.metadata_flat), delimiter="\t")}
        else:
            study = {r["accession"]: (r.get("study") or "unknown") for r in
                     csv.DictReader(open(pathlib.Path(args.metadata_root) / drug /
                                         "metadata.tsv"), delimiter="\t")}
        fold_study = {int(r["fold"]): r["study"] for r in
                      csv.DictReader(open(root / f"{drug}_studies.tsv"),
                                     delimiter="\t") if int(r["fold"]) >= 0}

        L = pd.read_csv(lo[0], sep="\t")
        L = L[L["aggregation"] == args.aggregation].copy()
        L["fold_i"] = L["fold"].str.rsplit("_f", n=1).str[1].astype(int)

        C = pd.read_csv(cv[0], sep="\t")
        C = C[C["aggregation"] == args.aggregation].copy()
        C["repeat"] = C["fold"].str.extract(r"_r(\d+)_")[0].astype(int)
        C["study"] = C["sample_id"].map(study)

        # The development set, and its resistance rate, taken from the CV arm:
        # it covers every development isolate exactly once per repeat.
        one_rep = C[C["repeat"] == C["repeat"].min()]
        dev_n = len(one_rep)
        dev_res = int((one_rep["log2_mic"] > b).sum())

        for fi, st in sorted(fold_study.items()):
            sub_l = L[L["fold_i"] == fi]
            ml = score(sub_l, b)
            # CV: the metric within each repeat, then averaged. Each isolate has
            # one out-of-fold prediction per repeat.
            per_rep = [score(g, b) for _, g in
                       C[C["study"] == st].groupby("repeat")]
            per_rep = [m for m in per_rep if m]
            if ml is None or not per_rep:
                continue
            mc = {k: float(np.mean([m[k] for m in per_rep]))
                  for k in ("cost", "vMER", "MER", "SR_f1", "bal_OT", "OT")}
            # Two covariates that decide how the transport cost is read.
            # frac_removed: how much of the training set this collection was, so
            #   a cost that is really a sample-size effect is visible as one.
            # rate_gap: how far its resistance rate sits from the remainder,
            #   which is what the cohort-prediction hypothesis says should drive
            #   the cost.
            rest_n = dev_n - ml["n"]
            rest_rate = (dev_res - ml["n_resistant"]) / max(1, rest_n)
            rows.append(dict(drug=drug, model=model, study=st, n=ml["n"],
                             n_res=ml["n_resistant"], dev_n=dev_n,
                             frac_removed=ml["n"] / dev_n,
                             res_rate=ml["n_resistant"] / ml["n"],
                             rate_gap=abs(ml["n_resistant"] / ml["n"] - rest_rate),
                             **{f"loso_{k}": ml[k] for k in
                                ("cost", "vMER", "MER", "SR_f1", "bal_OT", "OT")},
                             **{f"cv_{k}": v for k, v in mc.items()}))

        sub = L[L["fold_i"] < 0]
        if not sub.empty:
            m = score(sub, b)
            pooled.append((drug, m["n"], m["cost"], m["SR_f1"]))

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("[ERROR] nothing scored")
    df["d_cost"] = df["loso_cost"] - df["cv_cost"]
    df["d_f1"] = df["loso_SR_f1"] - df["cv_SR_f1"]
    df["d_vMER"] = df["loso_vMER"] - df["cv_vMER"]

    print(f"\nLeave-one-study-out vs matched random CV, same isolates "
          f"({args.aggregation} aggregation)")
    print(f"{'drug':<24}{'study':<26}{'n':>5}{'R%':>5}"
          f"{'cost LOSO':>11}{'CV':>8}{'Δ':>8}{'F1 LOSO':>9}{'CV':>7}{'Δ':>8}")
    print("-" * 119)
    for drug, g in df.groupby("drug"):
        for _, r in g.iterrows():
            print(f"{r.drug:<24}{r.study.split(' et')[0].split(' (')[0][:25]:<26}"
                  f"{int(r.n):>5}{r.n_res / r.n:>5.0%}"
                  f"{r.loso_cost:>11.3f}{r.cv_cost:>8.3f}{r.d_cost:>+8.3f}"
                  f"{r.loso_SR_f1:>9.3f}{r.cv_SR_f1:>7.3f}{r.d_f1:>+8.3f}")
    print("-" * 119)

    worse = int((df.d_cost > 0).sum())
    print(f"\n{len(df)} (drug, study) pairs across {df.drug.nunique()} drugs.")
    print(f"Transport cost: median Δcost {df.d_cost.median():+.3f}, "
          f"mean {df.d_cost.mean():+.3f}, worst {df.d_cost.max():+.3f} "
          f"({worse} of {len(df)} worse under LOSO).")
    print(f"Δ S/R F1: median {df.d_f1.median():+.3f}, worst {df.d_f1.min():+.3f}.")
    print("\nPer drug (median Δcost across its studies):")
    for drug, g in df.groupby("drug"):
        print(f"  {drug:<24}{g.d_cost.median():>+8.3f}   "
              f"({len(g)} studies, worst {g.d_cost.max():+.3f})")

    if pooled:
        print("\nPooled fold (-1): the collections too small to score on their "
              "own, held out together.")
        for drug, n, c, f1 in pooled:
            print(f"  {drug:<24}n={n:<6}cost {c:.3f}   F1 {f1:.3f}")

    if args.out:
        df.to_csv(args.out, sep="\t", index=False)
        print(f"\n[OK] {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
