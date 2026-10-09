#!/usr/bin/env python3
"""Per-drug labels for P. aeruginosa: clamped, operator-aware, Rebelo excluded.

MTB never needed this. CRyPTIC plates give an exact MIC for every isolate, so
`OT`, `exact` and `MAE` were measured against a true value. Here 7-43% of MICs
are censored -- recorded as `>16` or `<=1` -- and a censored value is a *bound*,
not a measurement. Scoring ordinal accuracy against a bound measures agreement
with something that is not the truth, and the MIC ladder is exactly what the MTB
panel established as the contribution no baseline competes with. So it has to be
handled before any label is written, not explained afterwards.

The rule is the one the cipro package already used, and it is better than
dropping every censored row. Clamp the ladder to `[LO, HI]` and keep an isolate
when its value is unambiguous *after clamping*:

  - exact (`=`)                                         -> keep, clipped
  - high-censored (`>`, `>=`) at or above `HI`          -> keep, clips to `HI`
  - low-censored  (`<`, `<=`) at or below `LO`          -> keep, clips to `LO`
  - censored bound strictly inside `[LO, HI]`           -> **drop**

`>16` with a ceiling of 16 lands on 16 either way, so it carries no ambiguity.
`>2` with the same ceiling could be 4, 8, 16 or more, so it is unknown and goes.

`LO`/`HI` default to the range of the *exact* measurements for that drug, which
is what ciprofloxacin used ([-4, 3]) and what the censoring pattern implies:
studies censor at their plate edges, so the exact range is the plate.

Rebelo et al. (2022) is excluded entirely -- 72% censored across 168 rows.

One consequence worth stating rather than discovering: meropenem loses the 184
isolates recorded `>8`. They are unambiguously resistant (the breakpoint is 8)
but their position on the ladder is unknown, and this rule is about the ladder.
Reporting S/R on a larger cohort than the ordinal metrics would be defensible
too; it is not what the previous pipeline did, so it is not what this does.
"""

from __future__ import annotations

import argparse
import collections
import csv
import math
import pathlib
import re
import sys

EXCLUDE_STUDIES = {"Rebelo et al. (2022)"}

# Two isolates, dropped by accession, for the reason the source package gives
# them: `explicit_noninteger_mic`. Their plates read 6 and 3 mg/L -- gradient-strip
# values off the doubling series -- and each would otherwise be the sole occupant
# of its MIC class. `snap()` builds its allowed set from the whole cache, so a
# class with no holdout isolate in it is still a legal prediction target for every
# holdout isolate: 19 of 315 amikacin and 49 of 393 ceftazidime holdout predictions
# landed on these two phantom rungs before they were removed.
#
# Dropping rather than rounding is the source package's own call, and it is the
# better one: rounding 6 up to 8 records a concentration the lab did not test.
EXPLICIT_EXCLUSIONS = {"amikacin": {"ERR2103531"},
                       "ceftazidime": {"SRR8737513"}}

# log2 MIC above which an isolate is resistant. Derived from the R/S labels
# already in the metadata (each drug's max susceptible sits one dilution below
# its min resistant) and cross-checked against the cipro package's hardcoded
# SR_BOUNDARY = -1.0, which this reproduces exactly.
BOUNDARY = {"ciprofloxacin": -1.0, "levofloxacin": 1.0, "tobramycin": 1.0,
            "imipenem": 2.0, "meropenem": 3.0, "ceftazidime": 3.0,
            "cefepime": 3.0, "aztreonam": 4.0, "amikacin": 4.0,
            "piperacillin-tazobactam": 4.0}

HIGH = {">", ">="}
LOW = {"<", "<="}
# Only a *strict* bound hides a dilution. `>=16` asserts that 16 itself is
# consistent with the reading, so 16 is the nearest consistent value and moving
# it to 32 would record a dilution the lab never claimed. `>16` rules 16 out, so
# the nearest consistent value is one dilution up. The operator is in the data;
# reading it costs nothing, and it is the difference between transcribing a
# measurement and inventing one.
STRICT = {">", "<"}

# Two schemas reach this script. The AMRFinder package splits the measurement
# into `log2_mic` + `mic_operator`; metadata/parsed/<drug>_direct.tsv keeps it
# whole in `mic_raw` (">8", "<=0.5", "0.25") with the label in `rs`. Normalising
# here rather than preprocessing means one rule, one place, either input.
RAW = re.compile(r"^\s*(<=|>=|<|>)?\s*([0-9]*\.?[0-9]+)\s*$")


def measurement(r):
    """row -> (operator, log2 MIC), or None when the row carries no usable MIC."""
    if r.get("log2_mic") not in (None, ""):
        try:
            return (r.get("mic_operator") or "=", float(r["log2_mic"]))
        except ValueError:
            return None
    m = RAW.match(r.get("mic_raw") or "")
    if not m:
        return None
    v = float(m.group(2))
    return (m.group(1) or "=", math.log2(v)) if v > 0 else None

# The clamp per drug, in log2 mg/L. These were passed as --clamp on the original
# runs and recorded nowhere, so re-running without them produced different
# cohorts -- ciprofloxacin kept 1,288 instead of 1,681, because a clamp derived
# from the exact range is wider than the one actually used and a wider clamp
# drops more censored bounds as "inside".
#
# Recovered by reading the emitted <drug>_log2_mic.tsv files: every kept value is
# clipped into [LO, HI], so the min and max of the cohort *are* the clamp. All
# ten reproduce the protocol's cohort sizes exactly, and ciprofloxacin's [-4, 3]
# matches the value the protocol records independently.
CLAMP = {"amikacin": (0.0, 5.0), "aztreonam": (0.0, 6.0),
         "cefepime": (-1.0, 5.0), "ceftazidime": (-1.0, 5.0),
         "ciprofloxacin": (-4.0, 3.0), "imipenem": (-1.0, 5.0),
         "levofloxacin": (-3.0, 3.0), "meropenem": (-3.0, 5.0),
         "piperacillin-tazobactam": (1.0, 6.0), "tobramycin": (-2.0, 4.0)}


def _clip(v, lo, hi):
    return min(max(v, lo), hi)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True,
                    help="<drug>/metadata.tsv from the AMRFinder package")
    ap.add_argument("--drug", required=True)
    ap.add_argument("--clamp", default="",
                    help="'LO,HI' to override; defaults to the exact-measurement range")
    ap.add_argument("--boundary", type=float, default=None,
                    help="overrides the table above")
    ap.add_argument("--snap-dilutions", action="store_true",
                    help="round an off-series MIC up to the next doubling "
                         "dilution instead of failing. The two known ones are "
                         "dropped by accession (EXPLICIT_EXCLUSIONS), so anything "
                         "reaching this is new, and a new one deserves a decision "
                         "rather than a default")
    ap.add_argument("--keep-off-series", action="store_true",
                    help="admit off-series MICs as their own class. Each becomes "
                         "a class of one that bal_OT weights like a class of 500, "
                         "and a legal prediction target on the holdout even when "
                         "no holdout isolate occupies it")
    ap.add_argument("--censoring", choices=("drop", "impute"), default="drop",
                    help="drop: the clipped-unambiguous rule -- a censored bound "
                         "strictly inside the clamp is discarded, which is what "
                         "the first panel did. impute: no isolate is discarded; "
                         "a strict high bound ('>') moves one dilution up and a "
                         "strict low bound ('<') one down, inclusive bounds "
                         "('>=', '<=') are taken at face value, then everything "
                         "is clipped into the clamp")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.metadata), delimiter="\t"))
    if not rows:
        raise SystemExit(f"[ERROR] {args.metadata} is empty")
    idc = next((c for c in rows[0] if c.lower() in
                ("accession", "sample_id", "sample", "id")), None)
    if idc is None:
        raise SystemExit(f"[ERROR] no id column: {list(rows[0])[:6]}")
    bound = args.boundary if args.boundary is not None else BOUNDARY.get(args.drug)
    if bound is None:
        raise SystemExit(f"[ERROR] no breakpoint for {args.drug!r}; pass --boundary")

    kept_study, dropped_study = [], collections.Counter()
    dropped_explicit: list[str] = []
    for r in rows:
        if r.get("log2_mic", "") in ("", "NA"):
            continue
        if r.get("study", "") in EXCLUDE_STUDIES:
            dropped_study[r.get("study", "")] += 1
            continue
        if r[idc] in EXPLICIT_EXCLUSIONS.get(args.drug, ()):
            dropped_explicit.append(r[idc])
            continue
        kept_study.append(r)
    for s, n in dropped_study.items():
        print(f"[INFO] excluded {n:,} rows from {s}")
    if dropped_explicit:
        print(f"[INFO] excluded {len(dropped_explicit)} isolate(s) by accession "
              f"(non-integer MIC): {', '.join(dropped_explicit)}")

    meas = {id(r): measurement(r) for r in kept_study}
    nomic = sum(1 for v in meas.values() if v is None)
    if nomic:
        print(f"[INFO] {nomic:,} rows carry no usable MIC")
    exact = [v for op, v in filter(None, meas.values()) if op == "="]
    if not exact:
        raise SystemExit("[ERROR] no exact measurements to derive a clamp from")
    if args.clamp:
        lo, hi = (float(x) for x in args.clamp.split(","))
    elif args.drug in CLAMP:
        lo, hi = CLAMP[args.drug]
    else:
        lo, hi = min(exact), max(exact)
        print(f"[WARN] no recorded clamp for {args.drug}; derived [{lo:g}, {hi:g}] "
              f"from the exact range. That is not what the panel used for the ten "
              f"drugs in CLAMP and will not reproduce their cohorts.")
    print(f"[INFO] {args.drug}: clamp [{lo:g}, {hi:g}] "
          f"(from {len(exact):,} exact measurements), breakpoint {bound:g}")
    if not lo < bound < hi:
        print(f"[WARN] the breakpoint {bound:g} is not strictly inside the clamp; "
              f"every isolate would fall on one side")

    out_rows, drop = [], collections.Counter()
    imputed = collections.Counter()
    snapped = 0
    for r in kept_study:
        if meas[id(r)] is None:
            continue
        op, v = meas[id(r)]
        if abs(v - round(v)) > 1e-6 and not args.keep_off_series:
            if not args.snap_dilutions:
                raise SystemExit(
                    f"[ERROR] {r[idc]} has an off-series MIC ({2**v:g} mg/L, "
                    f"log2 {v:g}) and is not in EXPLICIT_EXCLUSIONS.\n"
                    f"        It would be the only isolate at its rung -- a class "
                    f"of one that bal_OT weights like a class of 500, and a legal\n"
                    f"        prediction target on the holdout with no holdout "
                    f"isolate in it. Decide: add it to EXPLICIT_EXCLUSIONS, or\n"
                    f"        pass --snap-dilutions to round it up, or "
                    f"--keep-off-series to admit it as its own class.")
            # up to the next doubling dilution, the convention for gradient-strip
            # MICs.
            v = math.ceil(v - 1e-9)
            snapped += 1
        if op == "=":
            pass
        elif op in HIGH:
            if args.censoring == "impute":
                # `>` excludes v, so the nearest consistent value is one dilution
                # up; `>=` admits v, so v stays. Clipping afterwards means a strict
                # bound already at the ceiling stays there too -- the clamp's whole
                # meaning is "at or beyond this", so there is nowhere further to
                # move it, and that isolate's value is not an assumption either.
                step = 1.0 if op in STRICT else 0.0
                imputed["high" if step else "high_incl"] += 1
                imputed["moved_high"] += _clip(v + step, lo, hi) != _clip(v, lo, hi)
                v = v + step
            elif v < hi:
                drop["high-censored inside the clamp"] += 1
                continue
        elif op in LOW:
            if args.censoring == "impute":
                step = 1.0 if op in STRICT else 0.0
                imputed["low" if step else "low_incl"] += 1
                imputed["moved_low"] += _clip(v - step, lo, hi) != _clip(v, lo, hi)
                v = v - step
            elif v > lo:
                drop["low-censored inside the clamp"] += 1
                continue
        else:
            drop[f"unknown operator {op!r}"] += 1
            continue
        out_rows.append((r[idc], min(max(v, lo), hi)))

    if snapped:
        print(f"[INFO] snapped {snapped:,} off-series MIC(s) up to the next "
              f"doubling dilution")
    for k, n in drop.most_common():
        print(f"[INFO] dropped {n:,}: {k}")
    if imputed:
        n = sum(imputed[k] for k in ("high", "low", "high_incl", "low_incl"))
        incl = imputed["high_incl"] + imputed["low_incl"]
        moved = imputed["moved_high"] + imputed["moved_low"]
        # Counting censored rows overstates the assumption twice over. An
        # inclusive bound is transcribed, not imputed -- nothing was assumed for
        # it. And a strict bound already at the clamp edge moves one step and is
        # clipped straight back, so its value is unchanged too. Only `moved`
        # rests on the rule.
        print(f"[INFO] {n:,} censored bounds "
              f"({imputed['high']:,} strict-high up, {imputed['low']:,} "
              f"strict-low down, {incl:,} inclusive at face value); "
              f"{moved:,} changed value\n"
              f"       ({moved / max(1, len(out_rows)):.1%} of the cohort has a "
              f"MIC that rests on the imputation rule; the rest are\n"
              f"       inclusive bounds or sat at a clamp edge, where the step "
              f"is absorbed by clipping)")
    if not out_rows:
        raise SystemExit("[ERROR] nothing survived the filter")

    hist = collections.Counter(v for _, v in out_rows)
    R = sum(n for v, n in hist.items() if v > bound)
    n = len(out_rows)
    print(f"\n{'log2 MIC':>9}{'n':>8}{'S/R':>5}{'~holdout':>10}")
    small = []
    for v in sorted(hist):
        h = round(hist[v] * 0.2)
        flag = "  <-- single digits held out" if h < 10 else ""
        if h < 10:
            small.append(v)
        print(f"{v:>9.1f}{hist[v]:>8,}{'R' if v > bound else 'S':>5}{h:>10,}{flag}")
    print(f"\n[INFO] {n:,} isolates kept of {len(kept_study):,} "
          f"({n/len(kept_study):.1%}); {R:,} resistant ({R/n:.1%})")
    if small:
        print(f"[WARN] {len(small)} class(es) under 10 in a 20% holdout: "
              f"{', '.join(f'{v:g}' for v in small)} — report bal_OT_common.")

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{args.drug}_log2_mic.tsv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "log2_mic"])
        for s, v in sorted(out_rows):
            w.writerow([s, f"{v:g}"])
    print(f"[OK] {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
