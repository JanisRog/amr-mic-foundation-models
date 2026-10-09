#!/usr/bin/env python3
"""Drop acquired-sequence targets whose gene is already a curated target.

Arm C adds an embedding of each AMRFinder hit's sequence. For a point mutation the
hit spans the whole gene, so `amrseq_gyrA_T83I` is an embedding of *gyrA* -- a gene
already embedded as a curated target, and embedded again as `amrseq_gyrA_D87N` for
any isolate carrying both mutations. Ciprofloxacin would therefore carry three
correlated 4-PC blocks of gyrA. That is not merely redundant: it dilutes
`max_features` sampling and makes any gain unattributable between the two
representations.

So an `amrseq_` target survives only if its gene is *not* in the drug's curated
locus set. What remains is the genuinely new material -- acquired genes absent
from PAO1 -- which is the whole reason arm C exists.

The gene behind a family:

  POINT families      the part before the first underscore: gyrA_T83I -> gyrA,
                      oprD_LOF -> oprD, mexZ_LOF -> mexZ
  AMR families        the family name itself: crpP, blaVIM, sul1, aph(3')-IIb

One alias is needed. Bakta annotates chromosomal *ampC* as `blaPDC-1`, so the
curated symbol for that locus is `blaPDC-1` while AMRFinder reports `blaPDC-3`,
`blaPDC-5`, `blaPDC-35` and so on. Those are alleles of the same chromosomal gene
we already embed, so they are matched on the `blaPDC` prefix rather than exactly.
Missing that would have left five copies of ampC in the cephalosporin arms.

Note what this does NOT drop: `blaOXA-50` and relatives. They are intrinsic to
P. aeruginosa but are not in the curated target set, so their sequence is new
information here even though their *presence* is near-universal.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

PREFIX_ALIASES = ("blaPDC",)


def gene_of(family: str) -> str:
    """family name (without the amrseq_ prefix) -> gene symbol"""
    for p in PREFIX_ALIASES:
        if family == p or family.startswith(p + "-"):
            return p
    gene = family.split("_")[0]
    for p in PREFIX_ALIASES:
        if gene == p or gene.startswith(p + "-"):
            return p
    return gene


def norm(sym: str) -> str:
    sym = (sym or "").strip()
    for p in PREFIX_ALIASES:
        if sym == p or sym.startswith(p + "-"):
            return p
    return sym.split("_")[0]          # mexZ_G195E -> mexZ


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clumps", required=True)
    ap.add_argument("--curated-targets", required=True,
                    help="the drug's curated target list")
    ap.add_argument("--amrseq-targets", required=True,
                    help="candidate amrseq_* targets for this drug")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sym_of_locus = {}
    with open(args.clumps) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sym_of_locus[r["reference_gene_id"]] = r["reference_gene_symbol"]

    curated = set()
    for l in open(args.curated_targets):
        t = l.strip()
        if not t or t.startswith("#"):
            continue
        lt = t.split(":")[0]
        s = norm(sym_of_locus.get(lt, ""))
        if s:
            curated.add(s)

    kept, dropped = [], []
    for l in open(args.amrseq_targets):
        t = l.strip()
        if not t or t.startswith("#"):
            continue
        fam = t[len("amrseq_"):] if t.startswith("amrseq_") else t
        (dropped if gene_of(fam) in curated else kept).append((t, gene_of(fam)))

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(f"{t}\n" for t, _ in kept))

    print(f"[INFO] curated genes: {len(curated)}  ({', '.join(sorted(curated))})")
    print(f"[OK] {out}: {len(kept)} kept, {len(dropped)} dropped as already curated")
    if dropped:
        print("     dropped: " + ", ".join(f"{t}({g})" for t, g in dropped[:12])
              + (" ..." if len(dropped) > 12 else ""))
    if kept:
        print("     kept   : " + ", ".join(t for t, _ in kept[:12])
              + (" ..." if len(kept) > 12 else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
