#!/usr/bin/env python3
"""Rebuild a Snippy VCF at a chosen allele fraction, for `bcftools consensus`.

Today's consensus sequences come from Snippy's `snps.vcf`, where freebayes' ploidy-2
genotyping has already reduced calls to `1/1`. That behaves as a step at allele
fraction ~0.9, and section 7 measured what it costs: at 0.2 every drug's vMER falls
and none loses more than 0.007 of MER -- isoniazid gains 0.011, bedaquiline 0.242.
The threshold was never acted on because it meant re-embedding 8,243 targets. On a
fixed set of ~34 it is cheap.

Nothing is re-called. `snps.raw.vcf` is freebayes' pre-filter output and carries AO
and DP at every site, so this is a parse and a rewrite, using **exactly** the rules
`tools/sweep_consensus_threshold.py` used to justify 0.2 -- QUAL >= 100, DP >= 10,
fraction = AO/DP per alt allele. Change those and it is no longer the same
experiment.

Two things need deciding here that did not exist at 0.9, and both are recorded in
the output header rather than left implicit:

**Competing alleles at one position.** At 0.9 at most one allele could pass. At 0.2
an isolate can pass two -- say 0.45 and 0.35 -- and once both are forced to `1/1`,
`bcftools consensus` applies whichever comes first and warns. The highest-AO allele
wins here, which matches how the sweep resolved competing fractions and keeps the
consensus a single well-defined sequence.

**Overlapping records.** Freebayes emits MNPs and haplotype calls, so a 4 bp record
at 100 and a SNP at 102 can both pass and cannot both be applied. Records are
accepted greedily in descending AO and any that overlaps an accepted one is
dropped, so the result does not depend on file order.

A consensus built this way is a *majority* sequence, not a diploid one: an isolate
with a 30% resistant subpopulation is rendered as fully resistant. That is the
intended reading -- section 7 found 14 of 14 variants where minority carriers
behave phenotypically like full carriers -- but it is a modelling choice about what
a mixed infection *is*, and it belongs in the methods.
"""

from __future__ import annotations

import argparse
import gzip
import pathlib
import sys


def open_text(p: pathlib.Path):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "r")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-vcf", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--threshold", type=float, required=True,
                    help="minimum AO/DP for an allele to enter the consensus")
    ap.add_argument("--min-qual", type=float, default=100.0,
                    help="Snippy's --minqual; the sweep used 100")
    ap.add_argument("--min-dp", type=int, default=10,
                    help="Snippy's --mincov; the sweep used 10")
    ap.add_argument("--stats", default="", help="optional one-line TSV of counts")
    args = ap.parse_args()

    src = pathlib.Path(args.in_vcf)
    if not src.exists():
        raise SystemExit(f"[ERROR] missing {src}")

    header: list[str] = []
    cand: list[tuple[int, int, int, str]] = []   # (ao, pos, end, line)
    n_records = n_qual = n_dp = n_frac = 0

    with open_text(src) as fh:
        for raw in fh:
            if raw.startswith("##"):
                header.append(raw.rstrip("\n"))
                continue
            if raw.startswith("#CHROM"):
                header.append(raw.rstrip("\n"))
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            n_records += 1
            try:
                qual = float(f[5])
            except ValueError:
                continue
            if qual < args.min_qual:
                continue
            n_qual += 1
            fmt, vals = f[8].split(":"), f[9].split(":")
            d = dict(zip(fmt, vals))
            try:
                dp = int(d.get("DP", 0))
            except ValueError:
                continue
            if dp < args.min_dp:
                continue
            n_dp += 1
            alts = f[4].split(",")
            aos = [x for x in d.get("AO", "").split(",") if x not in ("", ".")]
            if len(aos) != len(alts):
                continue
            best = None
            for alt, ao in zip(alts, aos):
                try:
                    ao_i = int(ao)
                except ValueError:
                    continue
                if dp <= 0 or ao_i / dp < args.threshold:
                    continue
                if best is None or ao_i > best[0]:
                    best = (ao_i, alt)
            if best is None:
                continue
            n_frac += 1
            ao_i, alt = best
            pos, ref = int(f[1]), f[3]
            line = "\t".join([f[0], f[1], f[2], ref, alt, f[5], "PASS",
                              f[7], "GT:DP:AO", f"1/1:{dp}:{ao_i}"])
            cand.append((ao_i, pos, pos + len(ref) - 1, line))

    # greedy by descending AO so the outcome does not depend on file order
    kept: list[tuple[int, int, str]] = []
    for ao_i, pos, end, line in sorted(cand, key=lambda c: (-c[0], c[1])):
        if any(pos <= e and s <= end for s, e, _ in kept):
            continue
        kept.append((pos, end, line))
    kept.sort(key=lambda k: k[0])

    out = pathlib.Path(args.out_vcf)
    out.parent.mkdir(parents=True, exist_ok=True)
    note = [
        f"##afFilter_threshold={args.threshold}",
        f"##afFilter_minQual={args.min_qual}",
        f"##afFilter_minDP={args.min_dp}",
        "##afFilter_rule=AO/DP per alt; highest-AO allele per position; "
        "overlapping records resolved greedily by AO",
        "##afFilter_source=" + str(src),
    ]
    with out.open("w") as fh:
        for h in header:
            if h.startswith("#CHROM"):
                for n in note:
                    fh.write(n + "\n")
            fh.write(h + "\n")
        for _, _, line in kept:
            fh.write(line + "\n")

    msg = (f"{src.parent.name}\t{n_records}\t{n_qual}\t{n_dp}\t{n_frac}\t"
           f"{len(kept)}\t{n_frac - len(kept)}")
    print(f"[OK] {out}  records {n_records:,} -> qual {n_qual:,} -> dp {n_dp:,} "
          f"-> frac {n_frac:,} -> kept {len(kept):,} "
          f"({n_frac - len(kept):,} dropped as overlapping)")
    if args.stats:
        p = pathlib.Path(args.stats)
        p.parent.mkdir(parents=True, exist_ok=True)
        new = not p.exists()
        with p.open("a") as fh:
            if new:
                fh.write("sample\tn_records\tpass_qual\tpass_dp\tpass_frac\t"
                         "kept\tdropped_overlap\n")
            fh.write(msg + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
