#!/usr/bin/env python3
"""ARDaP's own per-drug S/R calls, as the P. aeruginosa comparator.

MTB's baseline was TB-Profiler's own drug call, and section 11 records what
substituting a reconstruction for it cost: 0.077 of development cost on
ethionamide, unnoticed until someone asked. AMRFinder cannot play that role here
because it reports gene hits with drug classes, not a per-drug S/R decision.
ARDaP can, and has already run on every isolate.

The call lives in `patientDrugSusceptibilityData.csv`, which Nextflow writes
inside its `work/<hash>/` scratch rather than to `Outputs/`. The published
`AbR_output.final.txt` is variant-level detail -- gene, change, effect, level --
from which a drug call could be *derived*, and deriving one is precisely the
substitution this project has already paid for once. So the scratch file is the
source, with two consequences handled explicitly:

  **Multiple copies.** Nextflow leaves one per process instance, typically two
  per isolate. They are compared, and a disagreement is fatal rather than
  resolved by picking the newer -- if the same tool produced two different calls
  for one isolate, neither is trustworthy without knowing why.
  **Scratch is deletable.** Once gathered here the calls are independent of the
  work directory, which is the point of running this before any cleanup.

ARDaP's database covers 8 of the panel's 10 drugs. Aztreonam and levofloxacin
are absent and get no comparator; ciprofloxacin's call is *not* substituted for
levofloxacin, however tempting the shared mechanism makes it.

**ARDaP emits three statuses, not two.** An earlier version of this script kept
only Sensitive and Resistant, which silently discarded every Intermediate call
and made ARDaP look as though it had declined to call 10% of isolates -- an
asymmetry that was an artefact of the parser, not a property of the tool. All
three are kept here, and how Intermediate maps onto a binary phenotype is left
to the consumer rather than decided in the parser.

`Level` is also kept. It is a graded score (-90 to 505) and Status is not a pure
function of it -- the same Level appears as both Intermediate and Resistant for
one drug -- so there is determinant-specific logic on top. But it ranks isolates,
which means ARDaP has an AUROC, and that separates *discrimination* from
*threshold placement*. Section 13 of the MTB protocol turned on exactly that
distinction: a comparator with good ranking and a badly placed cut is a different
object from one that cannot tell the classes apart.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import sys

# ARDaP drug name -> panel drug name. Only the eight it actually reports.
DRUG_MAP = {
    "ciprofloxacin": "ciprofloxacin",
    "meropenem": "meropenem",
    "imipenem": "imipenem",
    "ceftazidime": "ceftazidime",
    "cefepime": "cefepime",
    "piperacillin+tazobactam": "piperacillin-tazobactam",
    "amikacin": "amikacin",
    "tobramycin": "tobramycin",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ardap-root", required=True,
                    help="results/ardap/all_samples")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-samples", type=int, default=1000)
    args = ap.parse_args()

    root = pathlib.Path(args.ardap_root)
    sample_dirs = sorted((root / "samples").iterdir()) if (root / "samples").is_dir() else []
    if not sample_dirs:
        raise SystemExit(f"[ERROR] no samples under {root}")
    print(f"[INFO] {len(sample_dirs):,} sample directories")

    calls: dict[str, dict[str, str]] = {}
    seen_drugs: collections.Counter = collections.Counter()
    seen_status: collections.Counter = collections.Counter()
    no_file, conflict = [], []

    for i, sd in enumerate(sample_dirs, 1):
        if not sd.is_dir():
            continue
        sample = sd.name
        found = sorted(sd.rglob("patientDrugSusceptibilityData.csv"))
        if not found:
            no_file.append(sample)
            continue
        per_copy = []
        for f in found:
            one: dict[str, str] = {}
            with f.open() as fh:
                for r in csv.DictReader(fh):
                    d = (r.get("Drug") or "").strip().lower()
                    st = (r.get("Status") or "").strip().lower()
                    if not d or st not in ("sensitive", "intermediate", "resistant"):
                        continue
                    seen_drugs[d] += 1
                    seen_status[st] += 1
                    if d in DRUG_MAP:
                        lvl = (r.get("Level") or "").strip()
                        one[DRUG_MAP[d]] = st[0].upper()      # S / I / R
                        one[DRUG_MAP[d] + "__level"] = lvl
            per_copy.append(one)
        first = per_copy[0]
        if any(c != first for c in per_copy[1:]):
            conflict.append(sample)
            continue
        if first:
            calls[sample] = first
        if i % 500 == 0:
            print(f"[INFO] {i:,}/{len(sample_dirs):,} read")

    if no_file:
        print(f"[WARN] {len(no_file):,} sample(s) had no drug-susceptibility file, "
              f"e.g. {no_file[:3]}")
    if conflict:
        raise SystemExit(
            f"[ERROR] {len(conflict):,} sample(s) have copies of "
            f"patientDrugSusceptibilityData.csv that disagree, e.g. "
            f"{conflict[:3]}.\n         Two different calls from one tool for one "
            f"isolate; resolve the cause\n         rather than picking a copy.")
    if len(calls) < args.min_samples:
        raise SystemExit(f"[ERROR] only {len(calls):,} samples with calls")

    print(f"\n[INFO] ARDaP reports {len(seen_drugs)} drugs; "
          f"{len(DRUG_MAP)} map onto the panel")
    for d, n in seen_drugs.most_common():
        mark = f"  -> {DRUG_MAP[d]}" if d in DRUG_MAP else "   (not in the panel)"
        print(f"       {d:<28}{n:>7,}{mark}")
    missing = {"aztreonam", "levofloxacin"}
    print(f"[INFO] no ARDaP comparator for: {', '.join(sorted(missing))}")

    print(f"[INFO] statuses seen: " +
          ", ".join(f"{k}={v:,}" for k, v in seen_status.most_common()))

    drugs = sorted(DRUG_MAP.values())
    cols = [c for d in drugs for c in (d, d + "__level")]
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id"] + cols)
        for s in sorted(calls):
            w.writerow([s] + [calls[s].get(c, "") for c in cols])

    print(f"\n{'drug':<26}{'called':>8}{'S':>7}{'I':>7}{'R':>7}"
          f"{'R rate':>8}{'R+I rate':>10}")
    for d in drugs:
        v = [calls[s].get(d, "") for s in calls]
        n = sum(1 for x in v if x in ("S", "I", "R"))
        S, I, R = (sum(1 for x in v if x == k) for k in ("S", "I", "R"))
        print(f"{d:<26}{n:>8,}{S:>7,}{I:>7,}{R:>7,}"
              f"{(R/n if n else 0):>8.1%}{((R+I)/n if n else 0):>10.1%}")
    print("\nIntermediate is reported, not resolved. Mapping it to R or to S is a"
          "\nchoice that changes every downstream number, so it belongs in the"
          "\nanalysis that states it, not in this parser.")
    print(f"\n[OK] {out}  ({len(calls):,} isolates)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
