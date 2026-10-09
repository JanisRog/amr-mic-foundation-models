#!/usr/bin/env python3
"""Per-drug log2 MIC tables from the CRyPTIC supplementary data.

Source: Supplementary Table 3 of the CRyPTIC consortium GWAS study (PLoS Biol
20:e3001755), which carries one row per isolate, one column per antimicrobial as
a three-letter code, and the ENA sample accession.

Written because three antimicrobials in that table were never brought into this
project -- linezolid, clofazimine and rifabutin -- and adding a drug should not
mean reconstructing where its phenotypes came from. The same script regenerates
any of the thirteen.

Emits `<drug>_log2_mic.tsv` with `sample_id` and `log2_mic`, the format
`cv_build_folds.py` expects, plus the resistant count at a supplied breakpoint so
that a drug too sparse to evaluate is visible before a sweep is spent on it.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

CODES = {"amikacin": "AMI", "bedaquiline": "BDQ", "clofazimine": "CFZ",
         "delamanid": "DLM", "ethambutol": "EMB", "ethionamide": "ETH",
         "isoniazid": "INH", "kanamycin": "KAN", "levofloxacin": "LEV",
         "linezolid": "LZD", "moxifloxacin": "MXF", "rifabutin": "RFB",
         "rifampicin": "RIF"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True,
                    help="journal.pbio.3001755.s002.xlsx")
    ap.add_argument("--drug", required=True, choices=sorted(CODES))
    ap.add_argument("--boundary", type=float, default=None,
                    help="log2 MIC above which an isolate is resistant. Only used "
                         "to report the resistant count; no default, because a "
                         "wrong breakpoint produces plausible-looking error rates")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    try:
        from openpyxl import load_workbook
    except ImportError:
        raise SystemExit("[ERROR] needs openpyxl (pip install openpyxl)")

    wb = load_workbook(args.xlsx, read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))

    # The sheet opens with four lines of prose before the header, so the header is
    # located by content rather than by a fixed index -- a row number would break
    # silently if the file were re-exported.
    hi = next((i for i, r in enumerate(rows)
               if r and str(r[0]).strip() == "Sample_IDs"), None)
    if hi is None:
        raise SystemExit(f"[ERROR] no 'Sample_IDs' header row in {args.xlsx}")
    hdr = [str(c).strip() if c is not None else "" for c in rows[hi]]
    code = CODES[args.drug]
    if code not in hdr:
        raise SystemExit(f"[ERROR] no column {code!r}; columns are {hdr}")
    ic, jc = hdr.index("Sample_IDs"), hdr.index(code)
    acc = hdr.index("ENA_SAMPLE") if "ENA_SAMPLE" in hdr else None

    out_rows, blank = [], 0
    for r in rows[hi + 1:]:
        if not r or r[ic] is None:
            continue
        v = r[jc] if jc < len(r) else None
        if v in (None, "", "NA"):
            blank += 1
            continue
        out_rows.append((str(r[ic]).strip(), float(v),
                         str(r[acc]).strip() if acc is not None else ""))
    if not out_rows:
        raise SystemExit(f"[ERROR] no {args.drug} measurements found")

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{args.drug}_log2_mic.tsv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "log2_mic"])
        for s, v, _ in sorted(out_rows):
            w.writerow([s, f"{v:g}"])

    print(f"[INFO] {args.drug} ({code}): {len(out_rows):,} isolates with an MIC, "
          f"{blank:,} without")
    if args.boundary is not None:
        R = sum(1 for _, v, _ in out_rows if v > args.boundary)
        hold = round(R * 0.2)
        print(f"[INFO] breakpoint {args.boundary:g}: {R:,} resistant "
              f"({R / len(out_rows):.1%}), approximately {hold} in a 20% holdout")
        if hold < 10:
            print(f"[WARN] fewer than ten resistant isolates would be held out. "
                  f"Error rates computed on that many are counts rather than "
                  f"estimates; decide whether this drug is scoreable before "
                  f"spending a sweep on it.")
    else:
        print("[INFO] no --boundary given, so no resistant count is reported")
    if acc is not None:
        print(f"[INFO] ENA sample accessions are in the source table "
              f"(e.g. {out_rows[0][2]}) for the data-availability statement")
    print(f"[OK] {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
