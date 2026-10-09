#!/usr/bin/env python3
"""Is the leave-one-lineage-out run valid? Deliberately does not say who won.

§19 commits the baseline-vs-gwasfiltered verdict to a single read on the
complete panel, and ethionamide's development arms are not in yet. But the run
can be checked for validity now, and should be: a run that turns out to be
broken after the verdict is read has cost the one read.

So this script computes only quantities that carry no information about
gwfrule - baseline on the drugs under test:

  1. fold shape      lineage folds must be strongly unequal (~50/35/8/5). Four
                     near-equal folds mean the grouping did not take and the
                     run is random folds wearing a lineage label.
  2. size control    the lineage folds and the size-matched random folds must
                     have the SAME fold sizes, hence the same training fraction
                     on every fold. This is the control §19 actually rests on:
                     the two arms then differ only in how isolates are grouped.

                     The size-matched numbers are NOT required to match the
                     ordinary random 5x3 numbers, and that comparison was wrong
                     to enforce. A 4-fold split sized [3688, 2512, 521, 356]
                     trains on 48% of the data when its largest fold is held
                     out, against a flat 80% for 5x3, and the unweighted mean
                     over four folds gives that worst fold a quarter of the
                     weight. Isoniazid moves -0.179 of bal_OT between the two
                     for that reason alone. §2 reported size-matched matching
                     its 15-fold numbers on the isoniazid cohort; it does not
                     reproduce on these cohorts and the 5x3 figures are printed
                     for information only.
  3. robust metrics  §2 found vMER transports across lineages intact. A large
                     move is a sign something other than composition changed.
  4. moxifloxacin    its two arms are identical by construction, so its lineage
                     delta must be 0. §19 pre-registers this as a wiring check;
                     it is excluded from the verdict median, so reporting it
                     leaks nothing about the drugs that are not.

Everything else is withheld. --i-accept-reading-early prints the verdict too,
and requires saying so on the command line, because that is a decision about
the evaluation contract and not a formatting preference.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics
import sys

ARMS = ("baseline", "gwfrule")
DEV_SUFFIX = {"baseline": ("both_final1k", "tbrule1k", "final1k"),
              "gwfrule": ("gwasfiltered_gwfrule1k",)}
CONTROL = "moxifloxacin"
PRIMARY_ORDINAL = "bal_OT_common"


def folds_file(p: pathlib.Path, aggregation="median"):
    f = p / "folds.tsv"
    if not f.exists():
        return None
    rows = [r for r in csv.DictReader(open(f), delimiter="\t")
            if r["aggregation"] == aggregation]
    return rows or None


def fold_sizes(d: pathlib.Path):
    """Isolate count per non-holdout fold, largest first."""
    f = d / "cv_folds.tsv"
    if not f.exists():
        return None
    c = {}
    for r in csv.DictReader(open(f), delimiter="\t"):
        k = int(r["fold"])
        if k >= 0:
            c[k] = c.get(k, 0) + 1
    return sorted(c.values(), reverse=True) or None


def mean_of(rows, key):
    if key == "cost":
        v = [2 * float(r["vMER"]) + float(r["MER"]) for r in rows]
    else:
        v = [float(r[key]) for r in rows if r.get(key) not in ("", "NA", None)]
    return statistics.mean(v) if v else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lineage-root", required=True,
                    help="directory holding cv_<drug>/lineage/...")
    ap.add_argument("--dev-root", required=True,
                    help="slim_arms, for the ordinary random 5x3 numbers")
    ap.add_argument("--model", default="evo2")
    ap.add_argument("--tolerance", type=float, default=0.02,
                    help="size-matched vs random 5x3 agreement on "
                         f"{PRIMARY_ORDINAL}")
    ap.add_argument("--i-accept-reading-early", action="store_true",
                    help="also print the verdict, on an incomplete panel, "
                         "spending §19's single read")
    args = ap.parse_args()

    lin = pathlib.Path(args.lineage_root)
    dev = pathlib.Path(args.dev_root)
    M = args.model
    drugs = sorted(d.name[3:] for d in lin.glob("cv_*") if d.is_dir())
    if not drugs:
        raise SystemExit(f"[ERROR] no cv_* directories under {lin}")

    problems = []

    # ---- 0. same model as development ------------------------------------
    # The first lineage run passed every check below while fitting the wrong
    # models: a hardcoded --top-k 8 cut every arm to 8 targets, so gwasfiltered
    # (35-72 targets in development) and baseline became the same size. Fold
    # shape, the size control and the moxifloxacin wiring check are all blind
    # to that -- moxifloxacin has 4 targets, under the cap. Feature parity with
    # the development run is checked first because nothing after it means
    # anything if it fails.
    def nfeat(run):
        r = folds_file(run)
        return int(float(r[0]["n_features"])) if r and "n_features" in r[0] else None
    print("0. Feature count: lineage and size-matched refits vs development")
    print(f"   {'drug':<18}{'arm':<10}{'dev':>6}{'lineage':>9}{'size-mat':>10}")
    for d in drugs:
        for a in ARMS:
            dv = None
            for suf in DEV_SUFFIX[a]:
                dv = nfeat(dev / f"cv_{d}" / "runs" / f"{M}_{suf}")
                if dv:
                    break
            lo = nfeat(lin / f"cv_{d}" / "lineage" / "runs" / f"{M}_lolo_{a}")
            sm = nfeat(lin / f"cv_{d}" / "lineage" / "sizematched" / "runs"
                       / f"{M}_sm_{a}")
            bad = dv is None or lo != dv or sm != dv
            print(f"   {d[:17]:<18}{a:<10}{str(dv):>6}{str(lo):>9}{str(sm):>10}"
                  + ("   <-- not the development model" if bad else ""))
            if bad:
                problems.append(f"{d}/{a}: n_features dev={dv} lineage={lo} "
                                f"size-matched={sm}")
    print()

    # ---- 1. fold shape ----------------------------------------------------
    print("1. Lineage fold shape  (expect strongly unequal, ~50/35/8/5)")
    print(f"   {'drug':<18}{'n folds':>8}{'sizes':>34}{'max/min':>9}")
    for d in drugs:
        f = lin / f"cv_{d}" / "lineage" / "cv_folds.tsv"
        if not f.exists():
            print(f"   {d[:17]:<18}{'-':>8}  no cv_folds.tsv")
            problems.append(f"{d}: no lineage cv_folds.tsv")
            continue
        c = {}
        for r in csv.DictReader(open(f), delimiter="\t"):
            k = int(r["fold"])
            if k >= 0:
                c[k] = c.get(k, 0) + 1
        sizes = sorted(c.values(), reverse=True)
        ratio = sizes[0] / sizes[-1] if sizes and sizes[-1] else float("nan")
        flag = "" if ratio >= 3 else "   <-- too even, grouping may not have taken"
        print(f"   {d[:17]:<18}{len(sizes):>8}{str(sizes):>34}{ratio:>9.1f}{flag}")
        if ratio < 3:
            problems.append(f"{d}: lineage folds nearly equal (max/min {ratio:.1f})")

    # ---- 2. the control §19 rests on -------------------------------------
    print("\n2. Lineage vs size-matched FOLD SIZES  (must match: same training "
          "fraction,\n   so the two differ only in grouping)")
    print(f"   {'drug':<18}{'lineage':>22}{'size-matched':>22}{'max|Δ|':>8}"
          f"{'train% worst':>13}")
    for d in drugs:
        a = fold_sizes(lin / f"cv_{d}" / "lineage")
        b = fold_sizes(lin / f"cv_{d}" / "lineage" / "sizematched")
        if a is None or b is None:
            which = "lineage" if a is None else "size-matched"
            print(f"   {d[:17]:<18}{'cv_folds.tsv missing: ' + which:>44}")
            problems.append(f"{d}: {which} cv_folds.tsv missing -- cannot "
                            f"verify the size control")
            continue
        if len(a) != len(b):
            print(f"   {d[:17]:<18}{str(a):>22}{str(b):>22}"
                  f"{'different n':>21}")
            problems.append(f"{d}: {len(a)} lineage folds against {len(b)} "
                            f"size-matched -- not a size control")
            continue
        worst = max(abs(x - y) for x, y in zip(a, b))
        tot = sum(a)
        flag = "" if worst <= 2 else "   <-- sizes differ"
        print(f"   {d[:17]:<18}{str(a):>22}{str(b):>22}{worst:>8}"
              f"{(tot - a[0]) / tot:>12.0%}{flag}")
        if worst > 2:
            problems.append(f"{d}: lineage and size-matched fold sizes differ "
                            f"by up to {worst} isolates")

    # Printed for information. A mismatch here is expected: see the module
    # docstring -- the fold count and training fraction differ, so this is not
    # a control and is not allowed to fail the run.
    print(f"\n2b. Size-matched vs ordinary 5x3 on {PRIMARY_ORDINAL} "
          f"(INFORMATION ONLY --\n    different fold count and training "
          f"fraction; not a validity condition)")
    print(f"   {'drug':<18}{'arm':<10}{'5x3':>9}{'size-mat':>10}{'Δ':>9}")
    for d in drugs:
        for a in ARMS:
            devrows = None
            for suf in DEV_SUFFIX[a]:
                devrows = folds_file(dev / f"cv_{d}" / "runs" / f"{M}_{suf}")
                if devrows:
                    break
            smrows = folds_file(lin / f"cv_{d}" / "lineage" / "sizematched"
                                / "runs" / f"{M}_sm_{a}")
            if not devrows or not smrows:
                continue
            v1 = mean_of(devrows, PRIMARY_ORDINAL)
            v2 = mean_of(smrows, PRIMARY_ORDINAL)
            print(f"   {d[:17]:<18}{a:<10}{v1:>9.4f}{v2:>10.4f}{v2 - v1:>+9.4f}")

    # ---- 3. robust metrics transport --------------------------------------
    print("\n3. vMER on lineage folds vs size-matched  (same fold sizes, so a "
          "move here\n   is population shift, not training size). §2 found "
          "vMER transporting intact\n   on the isoniazid cohort: 0.074 -> "
          "0.072.")
    print(f"   {'drug':<18}{'arm':<10}{'size-mat':>10}{'lineage':>9}{'Δ':>9}")
    for d in drugs:
        for a in ARMS:
            sm = folds_file(lin / f"cv_{d}" / "lineage" / "sizematched" / "runs"
                            / f"{M}_sm_{a}")
            lo = folds_file(lin / f"cv_{d}" / "lineage" / "runs"
                            / f"{M}_lolo_{a}")
            if not sm or not lo:
                continue
            v1, v2 = mean_of(sm, "vMER"), mean_of(lo, "vMER")
            print(f"   {d[:17]:<18}{a:<10}{v1:>10.4f}{v2:>9.4f}{v2 - v1:>+9.4f}")

    # ---- 4. the wiring control -------------------------------------------
    print(f"\n4. {CONTROL} wiring control  (arms identical by construction, "
          f"Δ must be 0)")
    ok = None
    for d in drugs:
        if d != CONTROL:
            continue
        rows = {a: folds_file(lin / f"cv_{d}" / "lineage" / "runs"
                              / f"{M}_lolo_{a}") for a in ARMS}
        if not all(rows.values()):
            print(f"   {CONTROL}: run missing, cannot check")
            problems.append(f"{CONTROL}: wiring control run missing")
            break
        per = {}
        for a in ARMS:
            per[a] = {r["fold"]: mean_of([r], PRIMARY_ORDINAL) for r in rows[a]}
        shared = sorted(set(per["baseline"]) & set(per["gwfrule"]))
        worst = max(abs(per["gwfrule"][f] - per["baseline"][f]) for f in shared) \
            if shared else float("nan")
        ok = worst < 1e-9
        print(f"   folds compared: {len(shared)}   max |Δ {PRIMARY_ORDINAL}| = "
              f"{worst:.9f}   {'OK' if ok else 'FAILED'}")
        if not ok:
            problems.append(f"{CONTROL}: identical arms differ by {worst:.2e} "
                            f"-- the run, not the arm, is wrong")

    print("\n" + "=" * 72)
    if problems:
        print("VALIDITY: problems found\n")
        for p in problems:
            print(f"  - {p}")
        print("\nFix these before the verdict is read. A verdict on a broken "
              "run spends §19's single read for nothing.")
    else:
        print("VALIDITY: all checks passed. The run is sound.")
        print("\nThe verdict is NOT printed here. §19 commits it to one read on "
              "the complete panel, and ethionamide's development arms are not "
              "in yet. When they land: re-run slim_mtb_runs.sh, submit the "
              "lineage array with --array=3, then analyse_lineage_arms.py.")

    if args.i_accept_reading_early:
        print("\n" + "=" * 72)
        print("Reading early, as explicitly requested. This spends §19's single "
              "read on an incomplete panel; the n must be stated wherever the "
              "result is quoted, and adding ethionamide later is a second look.")
        print("Run: analyse_lineage_arms.py --root <lineage root> --model "
              + M)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
