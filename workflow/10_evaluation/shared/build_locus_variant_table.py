#!/usr/bin/env python3
"""Every variant inside a model's target loci, one binary column each.

The comparator for "embeddings only": the same loci, the same variant calls (the
AF 0.1 Snippy VCFs the consensus sequences were built from), but each variant as a
column instead of the locus as an embedding. Nothing is curated -- a variant is a
column if enough development isolates carry it -- so the two models differ in the
representation and nothing else.

Generalises the earlier P. aeruginosa variant-binary builder (working repository,
pseaer/) in three ways, so it serves both organisms
and a held-out score:

  rows and floor  rows are every isolate of the final fold file (development and
                  held-out); the carrier floor and ceiling are counted on the
                  development rows only (fold == --floor-fold), so the column set
                  is fixed before any held-out isolate is seen
  intervals       merged before lookup. A gene's upstream window can overlap its
                  neighbour's gene body; looking up only the last-starting interval
                  silently drops variants inside an earlier, longer one
  extra block     --extra-tsv prepends per-isolate columns (P. aeruginosa: acquired-
                  gene presence). --drop-point-from removes AMRFinder's curated
                  point mutations from it, which are variants already counted here

A variant is keyed by position, reference and first alternate allele (the AF 0.1
VCFs are forced 1/1 on the highest-AO allele). SNPs, MNPs and indels are all kept.
Isolates without a VCF stop the build unless --allow-missing: a missing file would
otherwise become a row of zeros, indistinguishable from a wild-type isolate.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import csv
import gzip
import json
import pathlib
import sys


def open_text(p: pathlib.Path):
    with open(p, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(p, "rt") if magic == b"\x1f\x8b" else open(p, "r")


def merged_intervals(clumps: pathlib.Path, wanted: set[str]):
    iv, found = [], set()
    with clumps.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["target_id"] in wanted:
                iv.append((int(r["target_start"]), int(r["target_end"])))
                found.add(r["target_id"])
    iv.sort()
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out, found


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snippy-runs", required=True)
    ap.add_argument("--vcf-name", default="snps.af010.vcf.gz")
    ap.add_argument("--clumps", required=True, help="target_id, target_start, target_end")
    ap.add_argument("--targets", required=True, help="the model's target list (first column)")
    ap.add_argument("--rows", required=True,
                    help="cv_folds_final.tsv: sample_id and fold; every isolate becomes a row")
    ap.add_argument("--floor-fold", default="-1",
                    help="fold value of the development rows, on which carriers are counted")
    ap.add_argument("--min-carriers", type=int, default=10)
    ap.add_argument("--max-carrier-frac", type=float, default=0.95)
    ap.add_argument("--extra-tsv", default="", help="per-isolate columns to prepend")
    ap.add_argument("--drop-point-from", default="",
                    help="amrfinder_columns.tsv; drops amr_<family> columns with subtype POINT")
    ap.add_argument("--allow-missing", action="store_true")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    want = {l.split("\t")[0].strip() for l in open(args.targets)
            if l.strip() and not l.startswith("#")}
    want.discard("target_id")
    iv, found = merged_intervals(pathlib.Path(args.clumps), want)
    if not iv:
        raise SystemExit(f"[ERROR] none of the {len(want)} targets are in {args.clumps}")
    starts = [s for s, _ in iv]
    print(f"[INFO] {len(found)} of {len(want)} targets have coordinates; "
          f"{len(iv)} merged intervals, {sum(e - s + 1 for s, e in iv):,} bp")
    if len(found) < len(want):
        print(f"[INFO] without coordinates (not in this reference, e.g. acquired genes): "
              f"{', '.join(sorted(want - found)[:8])}{' ...' if len(want - found) > 8 else ''}")

    ids, dev, seen = [], set(), set()
    for r in csv.DictReader(open(args.rows), delimiter="\t"):
        s = r["sample_id"]
        if s not in seen:
            seen.add(s)
            ids.append(s)
        if r["fold"] == args.floor_fold:
            dev.add(s)
    if not dev:
        raise SystemExit(f"[ERROR] no rows with fold == {args.floor_fold} in {args.rows}")
    print(f"[INFO] {len(ids):,} isolates, {len(dev):,} development (carriers counted on these)")

    carriers: dict[str, set[str]] = collections.defaultdict(set)
    root = pathlib.Path(args.snippy_runs)
    missing = []
    for i, s in enumerate(ids, 1):
        p = root / s / args.vcf_name
        if not p.exists():
            missing.append(s)
            continue
        with open_text(p) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 5)
                if len(f) < 5 or f[4] in (".", ""):
                    continue
                pos = int(f[1])
                j = bisect.bisect_right(starts, pos) - 1
                if j < 0 or pos > iv[j][1]:
                    continue
                carriers[f"var_{pos}_{f[3]}>{f[4].split(',')[0]}"].add(s)
        if i % 1000 == 0:
            print(f"[INFO] {i:,}/{len(ids):,} VCFs read")
    if missing:
        msg = f"{len(missing):,} isolate(s) have no {args.vcf_name}, e.g. {missing[0]}"
        if not args.allow_missing:
            raise SystemExit(f"[ERROR] {msg}")
        print(f"[WARN] {msg}; written as rows of zeros")

    n_dev = len(dev)
    dev_count = {v: len(c & dev) for v, c in carriers.items()}
    keep = sorted((v for v, k in dev_count.items()
                   if args.min_carriers <= k <= args.max_carrier_frac * n_dev),
                  key=lambda v: (int(v.split("_")[1]), v))
    too_few = sum(1 for k in dev_count.values() if k < args.min_carriers)
    too_many = sum(1 for k in dev_count.values() if k > args.max_carrier_frac * n_dev)
    print(f"[INFO] {len(carriers):,} distinct variants -> {len(keep):,} columns "
          f"({too_few:,} under {args.min_carriers} development carriers, "
          f"{too_many:,} over {args.max_carrier_frac:.0%})")
    if not keep:
        raise SystemExit("[ERROR] no variant passes the filter")

    pre_cols, pre = [], {}
    if args.extra_tsv:
        rows = list(csv.DictReader(open(args.extra_tsv), delimiter="\t"))
        idc = list(rows[0])[0]
        pre_cols = [c for c in rows[0] if c != idc]
        if args.drop_point_from:
            meta = list(csv.DictReader(open(args.drop_point_from), delimiter="\t"))
            point = {f"amr_{r['family']}" for r in meta if r["subtype"] == "POINT"}
            # A column the metadata does not describe cannot be classified; keeping
            # it could keep a curated point mutation. Stop instead.
            unknown = [c for c in pre_cols if c not in {f"amr_{r['family']}" for r in meta}]
            if unknown:
                raise SystemExit(f"[ERROR] {len(unknown)} column(s) of {args.extra_tsv} are not in "
                                 f"{args.drop_point_from}, e.g. {unknown[0]}; wrong metadata file?")
            dropped = [c for c in pre_cols if c in point]
            pre_cols = [c for c in pre_cols if c not in point]
            print(f"[INFO] extra block: dropped {len(dropped)} curated point-mutation column(s)"
                  + (f" ({', '.join(dropped[:6])}{' ...' if len(dropped) > 6 else ''})" if dropped else ""))
        pre = {r[idc]: r for r in rows}
        absent = [s for s in ids if s not in pre]
        if absent:
            raise SystemExit(f"[ERROR] {len(absent):,} isolate(s) missing from {args.extra_tsv}, "
                             f"e.g. {absent[0]}")
        print(f"[INFO] prepending {len(pre_cols)} column(s) from {pathlib.Path(args.extra_tsv).name}")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id"] + pre_cols + keep)
        for s in ids:
            w.writerow([s] + [pre[s][c] for c in pre_cols]
                       + [1 if s in carriers[v] else 0 for v in keep])
    manifest = {
        "targets": args.targets, "targets_requested": len(want), "targets_with_coordinates": len(found),
        "merged_intervals": len(iv), "bp": sum(e - s + 1 for s, e in iv),
        "vcf_name": args.vcf_name, "rows": len(ids), "development_rows": n_dev,
        "min_carriers": args.min_carriers, "max_carrier_frac": args.max_carrier_frac,
        "distinct_variants": len(carriers), "variant_columns": len(keep),
        "extra_tsv": args.extra_tsv, "extra_columns": len(pre_cols), "missing_vcfs": len(missing),
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    cc = sorted(dev_count[v] for v in keep)
    print(f"[OK] {out}: {len(ids):,} isolates x {len(pre_cols) + len(keep)} columns "
          f"({len(pre_cols)} extra + {len(keep)} variants)")
    print(f"     development carriers per variant: min {cc[0]}, median {cc[len(cc) // 2]}, max {cc[-1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
