#!/usr/bin/env python3
"""The curated gene / intergenic-region list for one drug, without openpyxl.

The cluster python has no openpyxl and installing into a shared analysis env to
read one table is a poor trade, so this parses the xlsx directly -- it is a zip of
XML. `tools/build_all_drug_mic_table.py` already does this for the MIC workbook,
but that one has a single worksheet and this workbook does not, so the sheet name
has to be resolved through xl/workbook.xml and its relationships rather than by
taking the only file in xl/worksheets/.

Why this list matters: TB-Profiler's tbdb.bed is keyed on genes, so a target set
built from it has no fabG1, no fabD and no oxyR' -- and embeddings restricted to
those genes miss 26% of resistant isolates, because the inhA promoter mutations
live in an interval the gene list never names. A curated table with a column for
*intergenic regions* names them, and it comes from a published paper rather than
from anyone's judgement about this cohort, which is what makes it usable as a
fixed target list without leaking anything.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
PKG = "{http://schemas.openxmlformats.org/package/2006/relationships}"


def sheet_path(z: zipfile.ZipFile, want: str) -> str:
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target_of = {r.get("Id"): r.get("Target")
                 for r in rels.iter(PKG + "Relationship")}
    names = []
    for sh in wb.iter(NS + "sheet"):
        name = sh.get("name") or ""
        names.append(name)
        if name.strip().casefold() == want.strip().casefold():
            # targets come either relative to xl/ ("worksheets/sheet2.xml") or
            # absolute from the package root ("/xl/worksheets/sheet2.xml")
            t = (target_of.get(sh.get(RNS + "id"), "") or "").lstrip("/")
            return t if t.startswith("xl/") else "xl/" + t
    raise SystemExit(f"[ERROR] no sheet named {want!r}; sheets are {names}")


def read_sheet(path: pathlib.Path, sheet: str) -> list[list[str]]:
    z = zipfile.ZipFile(path)
    shared: list[str] = []
    if "xl/sharedStrings.xml" in z.namelist():
        ss = ET.fromstring(z.read("xl/sharedStrings.xml"))
        shared = ["".join(t.text or "" for t in si.iter(NS + "t"))
                  for si in ss.iter(NS + "si")]
    ws = ET.fromstring(z.read(sheet_path(z, sheet)))
    rows: list[list[str]] = []
    for row in ws.iter(NS + "row"):
        cells: list[str] = []
        for c in row.iter(NS + "c"):
            ref = c.get("r") or ""            # empty cells are simply absent
            col = 0
            for ch in ref:
                if not ch.isalpha():
                    break
                col = col * 26 + (ord(ch.upper()) - 64)
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
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--sheet", default="Long format")
    ap.add_argument("--drug", default="Isoniazid")
    ap.add_argument("--drug-column", default="Drug")
    ap.add_argument("--entry-column", default="Gene or intergenic region")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    rows = read_sheet(pathlib.Path(args.xlsx), args.sheet)
    if not rows:
        raise SystemExit("[ERROR] empty sheet")
    header = [h.strip() for h in rows[0]]
    idx = {h: i for i, h in enumerate(header)}
    for want in (args.drug_column, args.entry_column):
        if want not in idx:
            raise SystemExit(f"[ERROR] column {want!r} missing; header is {header}")
    di, ei = idx[args.drug_column], idx[args.entry_column]

    drugs = sorted({r[di].strip() for r in rows[1:]
                    if len(r) > di and r[di].strip()})
    print(f"sheet {args.sheet!r}: {len(rows)-1} rows, drugs = {drugs}")

    entries, seen = [], set()
    for r in rows[1:]:
        if len(r) <= max(di, ei):
            continue
        if r[di].strip().casefold() != args.drug.casefold():
            continue
        e = r[ei].strip()
        if e and e.lower() != "nan" and e not in seen:
            seen.add(e)
            entries.append(e)
    if not entries:
        raise SystemExit(f"[ERROR] no entries for drug={args.drug!r}")

    inter = [e for e in entries if re.search(r"[:/]|inter", e, re.I)]
    print(f"\n{args.drug}: {len(entries)} entries, "
          f"{len(inter)} of them intergenic")
    for e in entries:
        print(f"    {e}{'   <- intergenic' if e in inter else ''}")

    if args.out:
        p = pathlib.Path(args.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["entry", "is_intergenic"])
            for e in entries:
                w.writerow([e, int(e in inter)])
        print(f"\n[OK] {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
