#!/usr/bin/env python3
"""Per-drug slices of the AMRFinder presence matrix, for --extra-feature-tsv.

`cv_run_architecture.py` appends every column of the TSV it is given, so the
matrix has to be sliced per drug before it is handed over -- otherwise
ciprofloxacin would be trained on carbapenemase columns and the comparison
between drugs would stop meaning anything. The per-drug column lists come from
AMRFinder's own Class/Subclass, so the slice is the tool's drug assignment, not
one invented here.

Isolates absent from the matrix are written as all-zero rather than dropped: a
sample AMRFinder never reported on carries no known determinant, which is what
zero means. Dropping it instead would silently shrink the cohort and make this
arm incomparable to the embeddings-only arm it is meant to be measured against.

**The carrier floor has to be applied here, not upstream.** The presence matrix
was filtered at >=10 carriers across all 2,126 isolates, but each drug's
development cohort is a subset -- 480 for cefepime against 1,456 for
ceftazidime -- so a gene with 10 carriers overall can have two, one or none in a
given drug's folds. The first run of this script produced cefepime columns with a
median of 8 carriers and a minimum of zero: constant columns that add width and
dilute every split. The floor is therefore recomputed against the isolates this
drug is actually trained on.
"""
from __future__ import annotations
import argparse, csv, pathlib, sys

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--presence", required=True)
    ap.add_argument("--columns", required=True, help="amr_columns_by_drug/<drug>.txt")
    ap.add_argument("--samples", required=True,
                    help="development fold TSV; decides which columns pass the "
                         "carrier floor")
    ap.add_argument("--emit-samples", default="",
                    help="fold TSV whose isolates get rows; defaults to "
                         "--samples. Pass cv_folds_final.tsv so the holdout is "
                         "covered without its isolates entering the floor")
    ap.add_argument("--min-carriers", type=int, default=10,
                    help="counted within this drug's cohort, not across the panel")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    want = [l.strip() for l in open(args.columns) if l.strip()]
    rows = {r["sample_id"]: r for r in
            csv.DictReader(open(args.presence), delimiter="\t")}
    if not rows:
        raise SystemExit(f"[ERROR] {args.presence} is empty")
    missing = [c for c in want if c not in next(iter(rows.values()))]
    if missing:
        raise SystemExit(f"[ERROR] {len(missing)} column(s) not in the matrix, "
                         f"e.g. {missing[:3]}")

    def cohort(path):
        ids, seen = [], set()
        for r in csv.DictReader(open(path), delimiter="\t"):
            s = r["sample_id"]
            if s not in seen:
                seen.add(s); ids.append(s)
        return ids

    # Columns are chosen on development; rows are written for whichever cohort
    # the consumer needs. The first version conflated the two and wrote rows for
    # the development isolates only, so the holdout scoring aborted -- correctly
    # -- on 339 isolates absent from the feature table. Widening --samples
    # instead would have let holdout isolates into the carrier floor, which is
    # feature selection seeing the test set.
    ids = cohort(args.samples)
    emit = cohort(args.emit_samples) if args.emit_samples else ids
    if args.emit_samples:
        extra = len(set(emit) - set(ids))
        print(f"[INFO] floor on {len(ids):,} development isolates; "
              f"writing {len(emit):,} rows ({extra:,} held out)")
    absent = [s for s in emit if s not in rows]

    carr_all = {c: sum(1 for s in ids if rows.get(s) and rows[s][c] == "1")
                for c in want}
    kept = [c for c in want if carr_all[c] >= args.min_carriers]
    dropped = [c for c in want if carr_all[c] < args.min_carriers]
    if not kept:
        raise SystemExit(f"[ERROR] no column reaches {args.min_carriers} carriers "
                         f"in this drug's {len(ids):,} isolates")
    want = kept

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id"] + want)
        for s in emit:
            r = rows.get(s)
            w.writerow([s] + [(r[c] if r else "0") for c in want])
    carr = {c: sum(1 for s in ids if rows.get(s) and rows[s][c] == "1")
            for c in want}
    # Report what was written, not what the floor was computed on -- the two
    # differ now and a summary line naming the smaller number for the larger
    # file is how a truncated feature table goes unnoticed.
    print(f"[OK] {out}: {len(emit):,} isolates x {len(want)} columns"
          + (f" (floor on {len(ids):,})" if len(emit) != len(ids) else "")
          + (f"  ({len(absent):,} absent from AMRFinder -> all zero)" if absent else ""))
    print(f"     carriers: min {min(carr.values())}, median "
          f"{sorted(carr.values())[len(carr)//2]}, max {max(carr.values())}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
