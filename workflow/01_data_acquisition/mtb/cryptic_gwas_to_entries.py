#!/usr/bin/env python3
"""Per-drug entry lists from the published CRyPTIC GWAS hits.

Source: S1 Data of the CRyPTIC GWAS study (PLoS Biol 20:e3001755) -- the manually
curated interpretation of its significant oligopeptide and oligonucleotide
signals, from GEMMA linear mixed models at Bonferroni thresholds. 771 rows across
thirteen antimicrobials.

This file is where the project's curated drug-gene table came from; §10 recorded
that provenance as unverified until 2026-10-08. Reading it directly, rather than a
transcription of it, is the point of this script: the table was hand-made for
isoniazid and never regenerated for the rest of the panel, which is why target
selection is inconsistent across drugs.

Emits `read_curated_drug_genes.py`'s format -- `entry`, `is_intergenic` -- so the
output feeds `build_panel_target_sets.py` unchanged and inherits its locus
resolution, including the three annotation traps and the adjacency promotion that
is the only reason fabG1 enters isoniazid's set at all.

**It does not filter.** The published hits include co-resistance markers: rpoB
appears under isoniazid with six significant signals, because 94.3% of
rifampicin-resistant isolates are also isoniazid-resistant. Filtering happens
downstream, where the catalogue annotation says which drug a gene is a mechanism
for; doing it here would duplicate that logic in a second place.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import sys

CODES = {"AMI": "amikacin", "BDQ": "bedaquiline", "CFZ": "clofazimine",
         "DLM": "delamanid", "EMB": "ethambutol", "ETH": "ethionamide",
         "INH": "isoniazid", "KAN": "kanamycin", "LEV": "levofloxacin",
         "LZD": "linezolid", "MXF": "moxifloxacin", "RFB": "rifabutin",
         "RIF": "rifampicin"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True, help="journal.pbio.3001755.s001.xlsx")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-neglog10p", type=float, default=0.0,
                    help="additional threshold on top of the paper's own "
                         "Bonferroni cut. 0 keeps every published hit, which is "
                         "the default: the rows are already the significant set")
    args = ap.parse_args()

    try:
        from openpyxl import load_workbook
    except ImportError:
        raise SystemExit("[ERROR] needs openpyxl (pip install openpyxl)")

    wb = load_workbook(args.xlsx, read_only=True, data_only=True)
    rows = list(wb[wb.sheetnames[0]].iter_rows(values_only=True))
    hi = next((i for i, r in enumerate(rows)
               if r and str(r[0]).strip().lower() == "drug"), None)
    if hi is None:
        raise SystemExit(f"[ERROR] no 'Drug' header row in {args.xlsx}")
    hdr = [str(c).strip() if c is not None else "" for c in rows[hi]]
    gi = next(i for i, h in enumerate(hdr) if h.lower().startswith("gene"))
    pi = next((i for i, h in enumerate(hdr) if "neglog10" in h.lower()), None)

    per: dict[str, set[str]] = collections.defaultdict(set)
    dropped = 0
    unknown = set()
    for r in rows[hi + 1:]:
        if not r or r[0] is None:
            continue
        code = str(r[0]).strip().upper()
        if code not in CODES:
            unknown.add(code)
            continue
        if pi is not None and args.min_neglog10p > 0:
            try:
                if float(r[pi]) < args.min_neglog10p:
                    dropped += 1
                    continue
            except (TypeError, ValueError):
                pass
        entry = str(r[gi]).strip()
        if entry and entry.lower() != "none":
            per[CODES[code]].add(entry)

    if unknown:
        print(f"[WARN] unrecognised drug codes ignored: {sorted(unknown)}")
    if dropped:
        print(f"[INFO] {dropped} hit(s) below --min-neglog10p {args.min_neglog10p}")

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    print(f"{'drug':<16}{'entries':>8}  intergenic  split")
    for drug in sorted(per):
        raw = sorted(per[drug])
        # ':' and '/' are NOT the same notation, and treating them alike is a
        # genome-sweeping bug rather than a cosmetic one.
        #
        #   'proA:ahpC'    an intergenic REGION between two flanking genes. The
        #                  loci lying between them belong to it -- that is how
        #                  oxyR' reaches the ahpC promoter set.
        #   'cysA2/cysA3'  an AMBIGUOUS assignment between paralogues: the signal
        #                  maps to one of them. They are not flanking anything.
        #
        # Marking 'cysA2/cysA3' intergenic made the builder span cysA2 (~Rv0815c)
        # to cysA3 (Rv3117) and admit every locus between -- 4,776 rifabutin
        # targets from twenty published hits. A '/' entry is therefore split into
        # its alternatives, each an ordinary gene entry.
        entries, split = [], 0
        for e in raw:
            if "/" in e and ":" not in e:
                parts = [x.strip() for x in e.split("/") if x.strip()]
                entries.extend(parts)
                split += 1
            else:
                entries.append(e)
        entries = sorted(set(entries))
        inter = {e for e in entries if ":" in e}
        p = out / f"{drug}_entries.tsv"
        with p.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["entry", "is_intergenic"])
            for e in entries:
                w.writerow([e, int(e in inter)])
        print(f"{drug:<16}{len(entries):>8}  {len(inter):<11} {split}")

    print(f"\n[OK] {len(per)} drug(s) under {out}")
    print("Pass these to build_panel_target_sets.py as repeated --curated-entries "
          "DRUG=FILE; the co-resistance filter is its annotation step, not this one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
