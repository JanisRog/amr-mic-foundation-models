#!/usr/bin/env python3
"""The two per-drug tables every downstream step needs, built the same way twice.

Isoniazid's and ethionamide's were made by hand, and one of them went wrong in a
way that survived into scored results: three files differing only by suffix, two
of them reconstructions of TB-Profiler's call rather than the call itself, and a
glob picked the wrong one. Ten drugs remain. This makes both tables from the two
authoritative sources, names them predictably, and records in the header where
each column came from.

  **`<drug>_log2_mic.tsv`** — `sample_id`, `log2_mic`. `cv_build_folds.py` wants
  a `log2_mic` column and the panel MIC table stores one column per drug
  (`RIF_int`, `INH_int`, ...), so the extraction has to happen somewhere; doing it
  here means the fold file and every later scorer read the same numbers.
  **`rule_bit_<code>_tbprofiler.tsv`** — `sample_id`, `rule_binary_<drug>`, taken
  from TB-Profiler's own multi-drug call matrix. **Not** reconstructed from
  "carries a catalogued variant at allele fraction >= t"; that substitution cost
  ethionamide 0.077 of development cost (section 11).

The MIC histogram is printed with the S/R split the given boundary implies,
because the boundary is the one parameter that fails silently: `SR_BOUNDARY`
defaults to isoniazid's -3.0, and a drug whose classes straddle -3 scores against
the wrong threshold without erroring. Classes small enough to be single digits in
a 20% holdout are flagged for the same reason -- bal_OT averages classes equally,
so a class of 40 decides a seventh of the headline on eight isolates.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import sys


def read_id_column(rows: list[dict], path: str) -> str:
    for c in rows[0]:
        if c.lower() in ("sample_id", "sample", "id", "sample_ids"):
            return c
    raise SystemExit(f"[ERROR] no sample id column in {path}: {list(rows[0])[:5]}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mic-tsv", required=True, help="mtb_all_drug_log2_mic.tsv")
    ap.add_argument("--mic-column", default="",
                    help="defaults to <CODE>_int, e.g. RIF_int")
    ap.add_argument("--tbprofiler-binary", required=True,
                    help="the MULTI-drug call matrix, not a single-drug rule bit")
    ap.add_argument("--drug", required=True, help="column in the matrix, e.g. rifampicin")
    ap.add_argument("--code", required=True, help="short code, e.g. RIF")
    ap.add_argument("--boundary", type=float, required=True,
                    help="log2 MIC above which an isolate is resistant")
    ap.add_argument("--holdout-frac", type=float, default=0.2)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    col = args.mic_column or f"{args.code}_int"
    mrows = list(csv.DictReader(open(args.mic_tsv), delimiter="\t"))
    mid = read_id_column(mrows, args.mic_tsv)
    if col not in mrows[0]:
        raise SystemExit(f"[ERROR] no {col!r} in {args.mic_tsv}; have "
                         f"{[c for c in mrows[0] if c.endswith('_int')]}")
    mic = {r[mid]: float(r[col]) for r in mrows
           if r.get(col, "") not in ("", "NA", "nan")}
    print(f"[INFO] {len(mic):,} isolates with a {args.code} MIC ({col})")

    trows = list(csv.DictReader(open(args.tbprofiler_binary), delimiter="\t"))
    tid = read_id_column(trows, args.tbprofiler_binary)
    if args.drug not in trows[0]:
        raise SystemExit(f"[ERROR] no {args.drug!r} column in "
                         f"{args.tbprofiler_binary}. This must be the multi-drug "
                         f"matrix; a single-drug rule bit has only its own column: "
                         f"{list(trows[0])[:6]}")
    call = {r[tid]: int(float(r[args.drug])) for r in trows
            if r.get(args.drug, "") not in ("", "NA", "nan")}
    print(f"[INFO] {len(call):,} isolates with a TB-Profiler {args.drug} call")

    both = sorted(set(mic) & set(call))
    print(f"[INFO] {len(both):,} with both "
          f"({len(set(mic) - set(call)):,} MIC without a call, "
          f"{len(set(call) - set(mic)):,} call without a MIC)")
    if not both:
        raise SystemExit("[ERROR] the two tables share no isolates — check the id columns")

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    mic_path = out / f"{args.drug}_log2_mic.tsv"
    with mic_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "log2_mic"])
        for s in both:
            w.writerow([s, f"{mic[s]:g}"])
    rule_path = out / f"rule_bit_{args.code.lower()}_tbprofiler.tsv"
    with rule_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", f"rule_binary_{args.drug}"])
        for s in both:
            w.writerow([s, call[s]])

    hist = collections.Counter(mic[s] for s in both)
    R = sum(n for m, n in hist.items() if m > args.boundary)
    print(f"\n{'log2 MIC':>9}{'n':>8}{'S/R':>5}{'est. holdout n':>16}")
    small = []
    for m in sorted(hist):
        h = round(hist[m] * args.holdout_frac)
        flag = "  <-- single digits in the holdout" if h < 10 else ""
        if h < 10:
            small.append(m)
        print(f"{m:>9.1f}{hist[m]:>8,}{'R' if m > args.boundary else 'S':>5}"
              f"{h:>16,}{flag}")
    print(f"\n[INFO] boundary {args.boundary:g}: {R:,} resistant / {len(both) - R:,} "
          f"susceptible ({R / len(both):.1%} resistant)")
    tb_R = sum(1 for s in both if call[s] == 1)
    print(f"[INFO] TB-Profiler calls {tb_R:,} resistant ({tb_R / len(both):.1%})")
    agree = sum(1 for s in both if (mic[s] > args.boundary) == (call[s] == 1))
    print(f"[INFO] phenotype and catalogue agree on {agree:,}/{len(both):,} "
          f"({agree / len(both):.1%}) at this boundary")

    if small:
        print(f"\n[WARN] {len(small)} MIC class(es) would be under 10 isolates in a "
              f"{args.holdout_frac:.0%} holdout: "
              f"{', '.join(f'{m:g}' for m in small)}.\n"
              f"       bal_OT weights classes equally, so each of these decides "
              f"1/{len(hist)} of the\n       headline on single digits. Report "
              f"bal_OT_common alongside it.")
    print(f"\n[OK] {mic_path}\n[OK] {rule_path}")
    print("\nThe boundary is not checked against anything here -- it is the one "
          "parameter that\nfails silently. Confirm it against antibiotics.json "
          "before building folds.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
