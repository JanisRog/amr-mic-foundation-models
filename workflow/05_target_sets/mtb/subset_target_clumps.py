#!/usr/bin/env python3
"""Cut the all-gene clumps file down to a fixed target list.

`13l_build_snp_gene_consensus_chain_sequences.py` reads target intervals from a
clumps TSV, and the all-gene one has 8,243 of them. The final model uses a fixed
list of a few dozen, and re-deriving their coordinates independently would be a
second implementation of something already on disk -- and a chance for the two to
disagree about where a gene starts. So this filters rather than rebuilds: every
retained row is byte-identical to the row the 8,243-target extraction used, which
is what makes the new sequences comparable to the old ones.

Fails loudly if any requested target is absent, because a silently shorter clumps
file produces a silently smaller embedding block, which would then be diagnosed as
a modelling result.
"""

from __future__ import annotations

import argparse
import pathlib
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clumps", required=True, help="all_gene_target_clumps.tsv")
    ap.add_argument("--targets", required=True, help="one target_id per line")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    want, seen = [], set()
    for line in pathlib.Path(args.targets).read_text().splitlines():
        t = line.split("\t")[0].strip()
        if t and not t.startswith("#") and t != "target_id" and t not in seen:
            seen.add(t)
            want.append(t)
    if not want:
        raise SystemExit(f"[ERROR] {args.targets} is empty")

    src = pathlib.Path(args.clumps).read_text().splitlines()
    if not src:
        raise SystemExit(f"[ERROR] {args.clumps} is empty")
    header = src[0]
    if header.split("\t")[0] != "target_id":
        raise SystemExit(f"[ERROR] first column of {args.clumps} is "
                         f"{header.split(chr(9))[0]!r}, expected 'target_id'")

    by_id = {}
    for line in src[1:]:
        if line.strip():
            by_id.setdefault(line.split("\t")[0], line)

    missing = [t for t in want if t not in by_id]
    if missing:
        raise SystemExit(f"[ERROR] {len(missing)} of {len(want)} targets are not in "
                         f"the clumps file: {missing[:8]}")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(header + "\n")
        for t in want:                      # caller's order, not file order
            fh.write(by_id[t] + "\n")
    print(f"[OK] {out}  {len(want)} of {len(by_id):,} targets kept")
    loci = {t.split(":")[-2] if len(t.split(":")) >= 2 else t for t in want}
    print(f"     {len(loci)} loci: {' '.join(sorted(loci))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
