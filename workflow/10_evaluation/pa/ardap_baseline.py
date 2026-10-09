#!/usr/bin/env python3
"""ARDaP's baseline rows, both Intermediate conventions and its threshold-free AUROC.

`catalogue_baseline_rows.py` takes a binary call column, and ARDaP does not give
one: it reports Sensitive / Intermediate / Resistant plus a graded `Level`.
Collapsing Intermediate is a modelling choice, not a parsing detail, and it moves
the numbers -- ceftazidime is 25.9% resistant with Intermediate folded to S and
35.0% with it folded to R. So both are computed and reported side by side rather
than one being chosen quietly.

`Level` matters more. A tool that ranks isolates correctly but draws its line in
the wrong place is a different object from one that cannot separate the classes,
and cost alone cannot tell them apart. Section 13 of the MTB protocol turned on
this: rifampicin's model had AUROC 0.984 alongside a mediocre vMER, and the
honest conclusion was that the ranking was nearly right and the threshold was
not. If ARDaP ranks cefepime well while calling almost everything sensitive,
quoting its binary cost of ~1.8 as the bar overstates what a model must beat.

Emits one `baseline_rows.tsv` per drug in the format the workbook already reads,
with three labels per split: the two Intermediate conventions and an AUROC-only
row carrying no error rates, because a ranking has no operating point.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys


def auroc(scores: list[float], labels: list[bool]) -> float:
    """Rank-based AUROC with ties averaged; nan when a class is absent."""
    n1 = sum(labels)
    n0 = len(labels) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    s1 = sum(r for r, l in zip(ranks, labels) if l)
    return (s1 - n1 * (n1 + 1) / 2) / (n1 * n0)


def metrics(pairs: list[tuple[bool, bool]]) -> dict:
    tp = sum(1 for t, c in pairs if t and c)
    fn = sum(1 for t, c in pairs if t and not c)
    fp = sum(1 for t, c in pairs if not t and c)
    tn = sum(1 for t, c in pairs if not t and not c)
    v = fn / (tp + fn) if tp + fn else float("nan")
    m = fp / (fp + tn) if fp + tn else float("nan")
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return {"vMER": v, "MER": m, "SR_f1": f1,
            "SR_acc": (tp + tn) / len(pairs) if pairs else float("nan"),
            "cost": 2 * v + m, "_n": len(pairs), "_R": tp + fn}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds-final", required=True)
    ap.add_argument("--calls", required=True, help="ardap_calls.tsv")
    ap.add_argument("--drug", required=True)
    ap.add_argument("--boundary", type=float, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    calls = {r["sample_id"]: r for r in
             csv.DictReader(open(args.calls), delimiter="\t")}
    if args.drug not in next(iter(calls.values())):
        raise SystemExit(f"[ERROR] no {args.drug!r} column in {args.calls}")

    split, mic = {}, {}
    for r in csv.DictReader(open(args.folds_final), delimiter="\t"):
        split[r["sample_id"]] = ("development_cv" if int(r["fold"]) < 0
                                 else "holdout")
        mic[r["sample_id"]] = float(r["log2_mic"])

    rows = []
    for sp in ("development_cv", "holdout"):
        ids = [s for s in split if split[s] == sp and s in calls
               and calls[s].get(args.drug) in ("S", "I", "R")]
        miss = sum(1 for s in split if split[s] == sp and s not in calls)
        truth = [mic[s] > args.boundary for s in ids]
        for label, resistant_set in (("ARDaP (I=R)", {"I", "R"}),
                                     ("ARDaP (I=S)", {"R"})):
            m = metrics([(t, calls[s][args.drug] in resistant_set)
                         for s, t in zip(ids, truth)])
            print(f"[INFO] {label:<14}{sp:<16}{m['_n']:>6,} isolates, "
                  f"{m['_R']:>5,} resistant   vMER {m['vMER']:.3f}  "
                  f"MER {m['MER']:.3f}  cost {m['cost']:.4f}")
            rows.append({"label": label, "split": sp,
                         **{k: f"{v:.6f}" for k, v in m.items()
                            if not k.startswith("_")}})
        lv = [s for s in ids if (calls[s].get(args.drug + "__level") or "").strip()
              not in ("", "NA")]
        a = auroc([float(calls[s][args.drug + "__level"]) for s in lv],
                  [mic[s] > args.boundary for s in lv])
        print(f"[INFO] {'ARDaP Level':<14}{sp:<16}{len(lv):>6,} isolates, "
              f"AUROC {a:.4f}")
        rows.append({"label": "ARDaP Level (ranking only)", "split": sp,
                     "SR_AUROC_raw": f"{a:.6f}"})
        if miss:
            print(f"[WARN] {miss:,} {sp} isolate(s) absent from the ARDaP table")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["label", "split", "vMER", "MER", "SR_f1", "SR_acc", "cost",
            "SR_AUROC_raw"]
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
