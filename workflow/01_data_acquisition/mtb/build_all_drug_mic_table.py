#!/usr/bin/env python3
"""One log2 MIC table for all 13 drugs, from the CRYPTIC metadata workbook.

Only isoniazid MICs were ever moved to the cluster, which is why every drug-level
question so far has had to wait. The workbook holds all of them; this converts it
once so the comparison scripts take a plain TSV.

Columns are the drug abbreviations (AMI BDQ CFZ DLM EMB ETH INH KAN LEV LZD MXF
RFB RIF); values are the raw log2 MIC, plus an integer-rounded copy suffixed
_int. Both are kept because the established runs use the rounded series -- the
isoniazid breakpoint is -3.0 there against -3.32193 raw -- and mixing them
silently shifts every threshold by a third of a dilution.

Verify with --check against an existing sample_splits_*.tsv, which carries both
log2_mic and log2_mic_raw: if this reproduces them for isoniazid it reproduces
them for the rest, since one workbook column feeds all of it.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys
import xml.etree.ElementTree as ET
import zipfile

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
DRUGS = ("AMI", "BDQ", "CFZ", "DLM", "EMB", "ETH", "INH",
         "KAN", "LEV", "LZD", "MXF", "RFB", "RIF")


def read_sheet(path: pathlib.Path) -> list[list[str]]:
    z = zipfile.ZipFile(path)
    shared: list[str] = []
    if "xl/sharedStrings.xml" in z.namelist():
        ss = ET.fromstring(z.read("xl/sharedStrings.xml"))
        shared = ["".join(t.text or "" for t in si.iter(NS + "t"))
                  for si in ss.iter(NS + "si")]
    names = [n for n in z.namelist() if n.startswith("xl/worksheets/sheet")]
    if len(names) != 1:
        raise SystemExit(f"[ERROR] expected one worksheet, found {names}")
    sheet = ET.fromstring(z.read(names[0]))

    rows: list[list[str]] = []
    for row in sheet.iter(NS + "row"):
        cells: list[str] = []
        for c in row.iter(NS + "c"):
            # cell references skip empty cells, so place by column letter
            ref = c.get("r") or ""
            col = 0
            for ch in ref:
                if ch.isalpha():
                    col = col * 26 + (ord(ch.upper()) - 64)
                else:
                    break
            while len(cells) < col - 1:
                cells.append("")
            v = c.find(NS + "v")
            if c.get("t") == "inlineStr":
                text = "".join(t.text or "" for t in c.iter(NS + "t"))
            elif v is None:
                text = ""
            elif c.get("t") == "s":
                text = shared[int(v.text)]
            else:
                text = v.text or ""
            cells.append(text)
        rows.append(cells)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata-xlsx", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--check", default="",
                    help="a sample_splits_*.tsv with log2_mic and log2_mic_raw")
    ap.add_argument("--check-drug", default="INH")
    args = ap.parse_args()

    rows = read_sheet(pathlib.Path(args.metadata_xlsx))
    header = rows[0]
    idx = {name: i for i, name in enumerate(header)}
    for want in ("Sample_IDs", *DRUGS):
        if want not in idx:
            raise SystemExit(f"[ERROR] column {want!r} missing; header is {header}")
    extra = [c for c in ("ENA_SAMPLE", "subject") if c in idx]

    out_rows = []
    counts = {d: 0 for d in DRUGS}
    for r in rows[1:]:
        if len(r) <= idx["Sample_IDs"] or not r[idx["Sample_IDs"]].strip():
            continue
        rec = {"sample_id": r[idx["Sample_IDs"]].strip()}
        for c in extra:
            rec[c] = r[idx[c]] if len(r) > idx[c] else ""
        for d in DRUGS:
            raw = r[idx[d]].strip() if len(r) > idx[d] else ""
            try:
                v = float(raw)
            except ValueError:
                rec[d] = ""
                rec[d + "_int"] = ""
                continue
            counts[d] += 1
            rec[d] = f"{v:g}"
            # the established runs round to the integer ladder; float round() is
            # banker's rounding, which would send -2.5 to -2 and 2.5 to 2
            rec[d + "_int"] = f"{int(v + 0.5) if v >= 0 else -int(-v + 0.5):g}"
        out_rows.append(rec)

    fields = ["sample_id", *extra]
    for d in DRUGS:
        fields += [d, d + "_int"]
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    print(f"[OK] {len(out_rows):,} samples -> {out}")
    print(f"\n{'drug':<6}{'with MIC':>10}")
    for d in DRUGS:
        print(f"{d:<6}{counts[d]:>10,}")

    if args.check:
        by_id = {r["sample_id"]: r for r in out_rows}
        n = ok_raw = ok_int = miss = 0
        bad = []
        with open(args.check) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if "log2_mic_raw" not in r:
                    raise SystemExit("[ERROR] --check file has no log2_mic_raw column")
                s = r["sample_id"]
                if s not in by_id or not by_id[s][args.check_drug]:
                    miss += 1
                    continue
                n += 1
                a = abs(float(by_id[s][args.check_drug]) - float(r["log2_mic_raw"])) < 1e-6
                b = abs(float(by_id[s][args.check_drug + "_int"]) - float(r["log2_mic"])) < 1e-6
                ok_raw += a
                ok_int += b
                if not (a and b) and len(bad) < 5:
                    bad.append((s, by_id[s][args.check_drug],
                                r["log2_mic_raw"], by_id[s][args.check_drug + "_int"],
                                r["log2_mic"]))
        print(f"\n[CHECK] {args.check_drug} against {pathlib.Path(args.check).name}: "
              f"{n:,} compared, raw match {ok_raw:,}, integer match {ok_int:,}, "
              f"{miss:,} not in the workbook")
        for row in bad:
            print(f"  mismatch {row}")
        if n and (ok_raw != n or ok_int != n):
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
