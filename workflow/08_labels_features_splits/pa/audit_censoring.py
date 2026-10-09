#!/usr/bin/env python3
"""What the censoring rule discards, per drug, and how much of it is decidable.

The protocol records kept-counts and rounded percentages; it does not record what
was dropped or why. This recomputes the rule against the parsed metadata and
reports the whole ledger, including the column the decision turns on: **how many
dropped isolates carry a recorded R/S label.**

A high-censored bound strictly inside the clamp is dropped because its ladder
position is unknown. Its S/R status is not unknown -- `>8` is resistant whatever
the true value -- and the source table records that status directly in `rs`, so
this counts rather than infers it. Those isolates are the ones the protocol's
suggested change ("keep for S/R, drop for ordinal") would recover.

Reads `<drug>_direct.tsv` from metadata/parsed, where `mic_raw` holds the
operator inline (`>8`, `<=0.5`) and `rs` the label. The clamp is derived the way
the builder derives it: from the exact measurements only.

Nothing is written. The operator sets, excluded studies and breakpoints are
imported from the builder rather than restated, so this cannot drift from what
produced the cohorts.
"""

from __future__ import annotations

import argparse
import collections
import csv
import math
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from prepare_pseaer_drug_inputs import (BOUNDARY, CLAMP, EXCLUDE_STUDIES,
                                        HIGH, LOW)

RAW = re.compile(r"^\s*(<=|>=|<|>)?\s*([0-9]*\.?[0-9]+)\s*$")


def parse(mic_raw: str):
    """'>8' -> ('>', 3.0). Returns None when the cell is not a MIC."""
    m = RAW.match(mic_raw or "")
    if not m:
        return None
    op, val = m.group(1) or "=", float(m.group(2))
    if val <= 0:
        return None
    return op, math.log2(val)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parsed-dir", required=True,
                    help="metadata/parsed, holding <drug>_direct.tsv")
    ap.add_argument("--drug", action="append", default=[])
    args = ap.parse_args()

    root = pathlib.Path(args.parsed_dir)
    print(f"{'drug':<24}{'rows':>7}{'no MIC':>8}{'study':>7}{'exact':>7}"
          f"{'kept':>7}{'drop hi':>9}{'drop lo':>9}{'  dropped R':>12}{'S':>6}")
    tot = collections.Counter()
    for drug in (args.drug or sorted(BOUNDARY)):
        # the parsed tables use underscores where the pipeline uses hyphens
        p = root / f"{drug.replace('-', '_')}_direct.tsv"
        if not p.exists():
            print(f"{drug:<24}  no file at {p.name}")
            continue
        rows = list(csv.DictReader(open(p), delimiter="\t"))
        parsed = [(r, parse(r.get("mic_raw", ""))) for r in rows]
        nomic = sum(1 for _, v in parsed if v is None)
        usable = [(r, v) for r, v in parsed if v is not None]
        excl = sum(1 for r, _ in usable if r.get("study", "") in EXCLUDE_STUDIES)
        keep = [(r, v) for r, v in usable if r.get("study", "") not in EXCLUDE_STUDIES]

        exact = [v for _, (op, v) in keep if op == "="]
        if not exact:
            print(f"{drug:<24}  no exact measurements")
            continue
        # the clamp the panel used, not one re-derived here -- see CLAMP
        lo, hi = CLAMP.get(drug, (min(exact), max(exact)))

        kept = dhi = dlo = dR = dS = 0
        for r, (op, v) in keep:
            if op == "=":
                kept += 1
                continue
            drop = (op in HIGH and v < hi) or (op in LOW and v > lo)
            if not drop:
                kept += 1
                continue
            dhi += op in HIGH
            dlo += op in LOW
            lab = (r.get("rs") or "").strip().upper()[:1]
            dR += lab == "R"
            dS += lab == "S"
        print(f"{drug:<24}{len(rows):>7,}{nomic:>8,}{excl:>7,}{len(exact):>7,}"
              f"{kept:>7,}{dhi:>9,}{dlo:>9,}{dR:>12,}{dS:>6,}")
        tot["hi"] += dhi; tot["lo"] += dlo; tot["R"] += dR; tot["S"] += dS

    print(f"\n{tot['hi'] + tot['lo']:,} dropped as censored-inside-clamp "
          f"({tot['hi']:,} high, {tot['lo']:,} low).")
    print(f"{tot['R']:,} of them are labelled R and {tot['S']:,} labelled S in the "
          f"source table —\ntheir S/R status is on record even though their MIC "
          f"is not, which is what the\nprotocol's 'keep for S/R, drop for ordinal' "
          f"change would recover.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
