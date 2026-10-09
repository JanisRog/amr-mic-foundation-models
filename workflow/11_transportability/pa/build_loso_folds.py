#!/usr/bin/env python3
"""Leave-one-study-out folds for one P. aeruginosa drug.

Section 1 declared study a post-hoc check rather than a fold constraint. This
builds the check. The panel pools eight collections whose resistance rates differ
enormously for the same drug -- ciprofloxacin runs from 18% resistant in
Cortes-Lara to 98% in Barrio-Tofino -- so under a random split every test isolate
has study-mates in training and a model can score well by recognising which
collection an isolate came from and predicting that collection's base rate. That
is the §1 failure mode in a new costume.

**Development isolates only.** The holdout is not touched: this is a declared
post-hoc check and must not spend the one-shot evaluation.

Fold convention, matching `cv_folds_final.tsv` and what `cv_run_architecture.py`
expects:

  fold = -1   always training, never scored
  fold = k    study k is the validation fold

Studies too small to score are assigned -1 rather than dropped. They are real
training data and removing them would change the model being tested as well as
the test; they simply never take a turn as the held-out collection.

A study is scoreable only with enough isolates AND enough of both classes. A
collection that is 100% resistant yields no major-error rate and an undefined
specificity, and reporting a cost for it would be reporting half a metric.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True,
                    help="<drug>/metadata.tsv, for the study column")
    ap.add_argument("--folds-final", required=True,
                    help="cv_folds_final.tsv; fold < 0 is the development set")
    ap.add_argument("--drug", required=True)
    ap.add_argument("--boundary", type=float, required=True,
                    help="log2 MIC above which an isolate is resistant; needed to "
                         "check both classes are represented, not to label")
    ap.add_argument("--min-isolates", type=int, default=40)
    ap.add_argument("--min-per-class", type=int, default=8)
    ap.add_argument("--min-train", type=int, default=100,
                    help="isolates that must REMAIN after a study is held out. "
                         "Tested because the held-out side being large enough "
                         "says nothing about whether anything is left to train "
                         "on: a drug drawn from one collection passes every "
                         "check on the study itself and then fits a PCA on zero "
                         "rows")
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", default="",
                    help="optional TSV recording which studies were scoreable")
    args = ap.parse_args()

    idc = None
    meta = list(csv.DictReader(open(args.metadata), delimiter="\t"))
    if meta:
        idc = next((c for c in meta[0] if c.lower() in
                    ("accession", "sample_id", "sample", "id")), None)
    if idc is None:
        raise SystemExit(f"[ERROR] no id column in {args.metadata}")
    study = {r[idc]: (r.get("study") or "unknown") for r in meta}

    dev = [(r["sample_id"], float(r["log2_mic"])) for r in
           csv.DictReader(open(args.folds_final), delimiter="\t")
           if int(r["fold"]) < 0]
    if not dev:
        raise SystemExit(f"[ERROR] no development rows in {args.folds_final}")

    missing = [s for s, _ in dev if s not in study]
    if missing:
        raise SystemExit(f"[ERROR] {len(missing)} development isolate(s) absent "
                         f"from {args.metadata}; first: {missing[:5]}")

    per = collections.defaultdict(lambda: [0, 0])     # study -> [resistant, susceptible]
    for s, v in dev:
        per[study[s]][0 if v > args.boundary else 1] += 1

    # Largest first, so fold 0 is the biggest collection and the numbering is
    # stable to read across drugs.
    tot_R = sum(v[0] for v in per.values())
    tot_S = sum(v[1] for v in per.values())

    order = sorted(per, key=lambda k: -sum(per[k]))
    fold_of, rows = {}, []
    k = 0
    for st in order:
        R, S = per[st]
        n = R + S
        # The remainder has to be a usable training set. Checking only the
        # held-out side passes a single-collection drug, which then fits a PCA on
        # zero rows and dies several minutes into the job.
        tr_R, tr_S = tot_R - R, tot_S - S
        tr_n = tr_R + tr_S
        why = ""
        if n < args.min_isolates:
            why = f"n<{args.min_isolates}"
        elif R < args.min_per_class:
            why = f"resistant<{args.min_per_class}"
        elif S < args.min_per_class:
            why = f"susceptible<{args.min_per_class}"
        elif tr_n < args.min_train:
            why = f"only {tr_n} would remain to train on"
        elif tr_R < args.min_per_class or tr_S < args.min_per_class:
            why = f"training remainder {tr_R}R/{tr_S}S lacks a class"
        ok = why == ""
        if ok:
            fold_of[st] = k
            k += 1
        else:
            fold_of[st] = -1
        rows.append((st, n, R, S, fold_of[st], why))

    print(f"[INFO] {args.drug}: {len(dev):,} development isolates, "
          f"{len(order)} studies, {k} scoreable")
    print(f"{'study':<36}{'n':>6}{'R':>6}{'S':>6}{'R%':>6}{'fold':>6}  note")
    for st, n, R, S, f, why in rows:
        print(f"{st[:35]:<36}{n:>6}{R:>6}{S:>6}{R/n:>6.0%}{f:>6}  {why}")

    if k == 0:
        print(f"[WARN] {args.drug}: no study is scoreable — leave-one-study-out "
              f"cannot be run for this drug. This is a property of the cohort "
              f"(a single-collection drug), not a failure.")
    elif k == 1:
        print(f"[WARN] {args.drug}: only one scoreable study, so the check rests "
              f"on a single held-out collection. Report it as such.")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "repeat", "fold", "log2_mic"])
        for s, v in sorted(dev):
            w.writerow([s, 0, fold_of[study[s]], f"{v:g}"])
    print(f"[OK] {out}")

    if args.summary:
        sp = pathlib.Path(args.summary)
        sp.parent.mkdir(parents=True, exist_ok=True)
        with sp.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["drug", "study", "n", "n_resistant", "n_susceptible",
                        "fold", "excluded_because"])
            for st, n, R, S, f, why in rows:
                w.writerow([args.drug, st, n, R, S, f, why])
        print(f"[OK] {sp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
