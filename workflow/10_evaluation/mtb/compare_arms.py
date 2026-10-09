#!/usr/bin/env python3
"""Paired comparison between feature arms, fold by fold.

The arms share drugs, folds and isolates, so the comparison is **paired** and the
test is Wilcoxon signed-rank on the per-fold differences. An unpaired test would
be wrong here and badly so: cost ranges from 0.08 (isoniazid) to 1.47 (delamanid)
across drugs while an arm difference is of order 0.01, so between-drug variance
would swamp the effect entirely.

Two levels, answering different questions:

  per drug    15 paired fold differences (5 folds x 3 repeats). Those folds are
              not independent -- each isolate appears in three validation folds
              and folds within a repeat partition the same data -- so the
              fold-win count is the primary statement and the p-value supports it.
  per panel   one paired difference per drug. Independent, but n = 10 means the
              smallest attainable two-sided p is 2/2^10 = 0.002, close to a sign
              test on the win count.

Reported for each metric in both directions, because §14's panel disagrees with
itself: the models shift toward specificity, which F1 rewards and the 2:1
very-major-error weighting punishes, and §10 records that wider target sets buy
MAE and cost bal_OT_common. An arm that wins one and loses the other is the
expected result, not an anomaly.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import statistics
import sys

# Isoniazid and rifampicin carry different suffixes: run_inh_redo.sbatch and
# run_rif_dev.sbatch both write <model>_<tag> (tbrule1k, final1k) while
# run_panel_dev.sbatch writes <model>_<targetset>_<tag> (both_final1k).
ARMS = {
    "baseline": ["both_final1k", "tbrule1k", "final1k"],
    "arm75":    ["both_arm75", "arm75"],
    "gwfrule":  ["gwasfiltered_gwfrule1k"],
    "gwfsnp":   ["gwasfiltered_gwfsnp1k"],
}
MODELS = ("evo2", "carbon", "dnabert")
METRICS = ("cost", "bal_OT_common", "bal_OT", "MAE", "vMER", "MER", "SR_f1")


def load(run_dir: pathlib.Path, aggregation: str):
    """fold -> metrics, for the one config this run was slimmed to."""
    f = run_dir / "folds.tsv"
    if not f.exists():
        return None
    out = {}
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["aggregation"] != aggregation:
            continue
        m = {k: float(r[k]) for k in r
             if k in METRICS and r[k] not in ("", "NA", None)}
        m["cost"] = 2 * float(r["vMER"]) + float(r["MER"])
        out[r["fold"]] = m
    return out or None


def find(root: pathlib.Path, drug: str, model: str, arm: str):
    for suffix in ARMS[arm]:
        d = root / f"cv_{drug}" / "runs" / f"{model}_{suffix}"
        if (d / "folds.tsv").exists():
            return d
    return None


def wilcoxon(diffs):
    try:
        from scipy.stats import wilcoxon as w
        if len(diffs) >= 6 and any(d != 0 for d in diffs):
            return float(w(diffs).pvalue)
    except Exception:
        pass
    return float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="slim_arms directory")
    ap.add_argument("--a", required=True, choices=sorted(ARMS))
    ap.add_argument("--b", required=True, choices=sorted(ARMS))
    ap.add_argument("--metric", default="cost", choices=METRICS)
    ap.add_argument("--aggregation", default="median")
    ap.add_argument("--model", default="", help="one model, or all three pooled")
    args = ap.parse_args()

    root = pathlib.Path(args.root)
    drugs = sorted(d.name[3:] for d in root.glob("cv_*") if d.is_dir())
    models = (args.model,) if args.model else MODELS
    lower_better = args.metric in ("cost", "MAE", "vMER", "MER")

    print(f"\n{args.a} vs {args.b}   metric {args.metric}   "
          f"({args.aggregation} aggregation; "
          f"{'lower' if lower_better else 'higher'} is better)")
    print(f"{'drug':<18}{'model':<9}{args.a:>10}{args.b:>10}{'Δ':>9}"
          f"{'folds':>8}{'p':>9}")
    print("-" * 73)

    per_drug = collections.defaultdict(list)
    rows = 0
    for drug in drugs:
        for model in models:
            da, db = find(root, drug, model, args.a), find(root, drug, model, args.b)
            if da is None or db is None:
                continue
            A, B = load(da, args.aggregation), load(db, args.aggregation)
            if not A or not B:
                continue
            shared = sorted(set(A) & set(B))
            if not shared:
                continue
            diffs = [A[f][args.metric] - B[f][args.metric] for f in shared]
            ma = statistics.mean(A[f][args.metric] for f in shared)
            mb = statistics.mean(B[f][args.metric] for f in shared)
            wins = sum(1 for d in diffs if (d < 0) == lower_better and d != 0)
            p = wilcoxon(diffs)
            per_drug[drug].append(ma - mb)
            rows += 1
            print(f"{drug:<18}{model:<9}{ma:>10.4f}{mb:>10.4f}{ma - mb:>+9.4f}"
                  f"{wins:>5}/{len(shared):<3}{p:>9.4f}")

    if not rows:
        raise SystemExit("[ERROR] no comparable runs found")
    print("-" * 73)

    # Per panel: one difference per drug, averaged over whichever models exist.
    drug_d = {d: statistics.mean(v) for d, v in per_drug.items() if v}
    vals = list(drug_d.values())
    better = sum(1 for v in vals if (v < 0) == lower_better)
    p = wilcoxon(vals)
    print(f"\nAcross {len(vals)} drugs (one paired difference each, models pooled):")
    print(f"  median Δ {statistics.median(vals):+.4f}   "
          f"{args.a} better on {better}/{len(vals)}   Wilcoxon p = {p:.4f}")
    print(f"  n = {len(vals)}; the smallest attainable two-sided p at this n is "
          f"{2 / 2 ** len(vals):.4f}")
    worst = max(drug_d.items(), key=lambda kv: abs(kv[1]))
    print(f"  largest movement: {worst[0]} {worst[1]:+.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
