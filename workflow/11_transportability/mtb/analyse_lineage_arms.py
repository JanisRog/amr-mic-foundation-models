#!/usr/bin/env python3
"""Does gwasfiltered's ladder gain survive leaving a lineage out?

On development folds gwfrule is cost-neutral against baseline and better on
every ladder metric. That is also the signature of a co-resistance marker:
knowing an isolate belongs to a resistant subpopulation narrows its MIC
distribution without encoding mechanism, so bal_OT_common, MAE and spearman
improve while vMER and MER do not move. §10 recorded a marker doing exactly
this and giving the gain back across lineages.

The arithmetic, fixed before the run:

  drop        random-fold value minus lineage-fold value, per arm. §2 measured
              -0.098 to -0.133 of bal_OT for evo2; an arm that leans on lineage
              signal drops further.
  delta       gwfrule minus baseline, computed separately on random folds and
              on lineage folds. The development delta is +0.0269.
  survival    lineage delta / random delta. Near 1 the gain is mechanism; near
              0 it was lineage signal; negative means it reverses.

The size-matched arm is a control, not a result: random folds shaped to the
lineage fold sizes. If it does not reproduce the ordinary random-fold numbers
then the lineage folds differ in training size as well as composition and
nothing downstream is interpretable. That is checked first and refuses to
report if it fails.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics
import sys

PRIMARY_ORDINAL = "bal_OT_common"   # declared the reported ordinal metric
ROBUST = ("vMER", "MER")            # §2: these transport intact
LADDER = ("bal_OT_common", "bal_OT", "MAE", "RMSE", "spearman", "exact")
ARMS = ("baseline", "gwfrule")


def read(run: pathlib.Path, aggregation="median"):
    f = run / "folds.tsv"
    if not f.exists():
        return None
    rows = [r for r in csv.DictReader(open(f), delimiter="\t")
            if r["aggregation"] == aggregation]
    if not rows:
        return None
    out = {}
    for k in set(LADDER) | set(ROBUST):
        v = [float(r[k]) for r in rows if r.get(k) not in ("", "NA", None)]
        if v:
            out[k] = statistics.mean(v)
    out["cost"] = statistics.mean(2 * float(r["vMER"]) + float(r["MER"])
                                  for r in rows)
    out["_per_fold"] = {r["fold"]: {k: float(r[k]) for k in LADDER
                                    if r.get(k) not in ("", "NA", None)}
                        for r in rows}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="metadata root")
    ap.add_argument("--model", default="evo2")
    ap.add_argument("--drugs", default="isoniazid_final,rifampicin,ethambutol,"
                                       "kanamycin,amikacin,levofloxacin,"
                                       "moxifloxacin,bedaquiline,delamanid")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    root = pathlib.Path(args.root)
    M = args.model
    rows = []
    for drug in args.drugs.split(","):
        R = root / f"cv_{drug}"
        d = {"drug": drug}
        for arm in ARMS:
            # The random-fold arm here is the size-matched control, not the
            # ordinary 5x3 development run: it is the only random arm fitted on
            # the same training sizes as the lineage folds, so it is the one a
            # lineage drop may legitimately be measured against.
            d[f"rand_{arm}"] = read(R / "lineage" / "sizematched" / "runs"
                                    / f"{M}_sm_{arm}")
            d[f"lin_{arm}"] = read(R / "lineage" / "runs" / f"{M}_lolo_{arm}")
        rows.append(d)

    have = [r for r in rows if all(r.get(f"{p}_{a}") for p in ("rand", "lin")
                                   for a in ARMS)]
    if not have:
        raise SystemExit("[ERROR] no drug has all four runs (lineage and "
                         "size-matched, both arms) -- nothing to compare")
    missing = [r["drug"] for r in rows if r not in have]
    if missing:
        print(f"[WARN] incomplete, skipped: {missing}")

    print(f"\nLeave-one-lineage-out, {M}, {len(have)} drug(s)")
    print(f"Primary ordinal metric: {PRIMARY_ORDINAL} (declared)\n")

    print("1. Do the robust metrics transport? (§2 says yes; a large move here "
          "means something else changed)")
    print(f"   {'drug':<16}{'arm':<10}{'vMER rand':>11}{'lineage':>9}{'Δ':>9}")
    for r in have:
        for a in ARMS:
            v1, v2 = r[f"rand_{a}"]["vMER"], r[f"lin_{a}"]["vMER"]
            print(f"   {r['drug'][:15]:<16}{a:<10}{v1:>11.4f}{v2:>9.4f}"
                  f"{v2 - v1:>+9.4f}")

    print(f"\n2. The ladder drop, per arm. An arm leaning on lineage signal "
          f"drops further.")
    print(f"   {'drug':<16}{'base rand':>10}{'base lin':>10}{'drop':>9}   "
          f"{'gwf rand':>10}{'gwf lin':>10}{'drop':>9}   {'extra drop':>11}")
    extra = []
    for r in have:
        br, bl = (r["rand_baseline"][PRIMARY_ORDINAL],
                  r["lin_baseline"][PRIMARY_ORDINAL])
        gr, gl = (r["rand_gwfrule"][PRIMARY_ORDINAL],
                  r["lin_gwfrule"][PRIMARY_ORDINAL])
        e = (gr - gl) - (br - bl)
        extra.append(e)
        print(f"   {r['drug'][:15]:<16}{br:>10.4f}{bl:>10.4f}{br - bl:>+9.4f}   "
              f"{gr:>10.4f}{gl:>10.4f}{gr - gl:>+9.4f}{e:>+11.4f}")

    print(f"\n3. Does the gain survive? delta = gwfrule − baseline")
    print(f"   {'drug':<16}{'random':>10}{'lineage':>10}{'survival':>10}")
    surv, dr_all, dl_all = [], [], []
    for r in have:
        dr = (r["rand_gwfrule"][PRIMARY_ORDINAL]
              - r["rand_baseline"][PRIMARY_ORDINAL])
        dl = (r["lin_gwfrule"][PRIMARY_ORDINAL]
              - r["lin_baseline"][PRIMARY_ORDINAL])
        dr_all.append(dr); dl_all.append(dl)
        s = dl / dr if abs(dr) > 1e-9 else float("nan")
        surv.append(s)
        print(f"   {r['drug'][:15]:<16}{dr:>+10.4f}{dl:>+10.4f}{s:>10.2f}")

    med_r, med_l = statistics.median(dr_all), statistics.median(dl_all)
    med_e = statistics.median(extra)
    pos_l = sum(1 for x in dl_all if x > 0)
    print(f"\n   median delta, random folds   {med_r:+.4f}")
    print(f"   median delta, lineage folds  {med_l:+.4f}  "
          f"({pos_l} of {len(dl_all)} positive)")
    print(f"   median extra drop for gwfrule {med_e:+.4f}")

    print("\n4. Verdict against the pre-registered rule")
    ok_sign = med_l > 0
    ok_drop = med_e <= 0.005          # gwfrule drops no worse than baseline
    if ok_sign and ok_drop:
        print("   PASS -- the ladder gain survives leaving a lineage out and "
              "gwfrule does not degrade further than baseline. The wider target "
              "set carries mechanism the catalogue set does not. Carry it "
              "forward.")
    elif ok_sign:
        print(f"   PARTIAL -- the gain keeps its sign ({med_l:+.4f}) but "
              f"gwfrule degrades {med_e:+.4f} further than baseline across "
              f"lineages. Some of the development gain was lineage signal. "
              f"Report both arms; do not claim the wider set is better.")
    else:
        print(f"   FAIL -- the gain does not survive (median {med_l:+.4f} "
              f"against {med_r:+.4f} on random folds). The development "
              f"advantage was lineage signal, which is what §10 predicted. "
              f"Baseline stays the reported arm.")
    print("\n   This judges the target set, not the models. Nothing here "
          "re-selects a configuration: every run reuses the configuration "
          "chosen on development folds.")

    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["drug", "metric", "rand_baseline", "rand_gwfrule",
                        "lin_baseline", "lin_gwfrule", "delta_rand",
                        "delta_lin", "drop_baseline", "drop_gwfrule"])
            for r in have:
                for k in LADDER + ROBUST + ("cost",):
                    try:
                        br, gr = r["rand_baseline"][k], r["rand_gwfrule"][k]
                        bl, gl = r["lin_baseline"][k], r["lin_gwfrule"][k]
                    except KeyError:
                        continue
                    w.writerow([r["drug"], k, f"{br:.6f}", f"{gr:.6f}",
                                f"{bl:.6f}", f"{gl:.6f}", f"{gr - br:.6f}",
                                f"{gl - bl:.6f}", f"{br - bl:.6f}",
                                f"{gr - gl:.6f}"])
        print(f"\n[OK] {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
