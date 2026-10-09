#!/usr/bin/env python3
"""Extract per-sample lineage calls from TB-Profiler result JSONs.

The gather step keeps only AMR gene hits and resistance mutations, so lineage --
the covariate needed to separate resistance signal from population structure --
is discarded. This recovers it.

Tolerant of TB-Profiler schema drift: main lineage has been called "main_lin"
and "main_lineage", sublineage "sublin" and "sub_lineage", and the per-call list
uses "lin"/"lineage" with "frac"/"fraction" depending on version.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Optional

MAIN_KEYS = ("main_lin", "main_lineage", "mainLineage")
SUB_KEYS = ("sublin", "sub_lineage", "subLineage")
LIST_KEYS = ("lineage", "lineages")
LIN_NAME_KEYS = ("lin", "lineage", "name")
FRAC_KEYS = ("frac", "fraction")


def first_str(obj: dict, keys) -> str:
    for k in keys:
        v = obj.get(k)
        if isinstance(v, str) and v.strip():
            return v.strip()
    return ""


def from_list(obj: dict) -> tuple[str, str]:
    """Fall back to the per-call list: highest-fraction call wins."""
    for lk in LIST_KEYS:
        entries = obj.get(lk)
        if not isinstance(entries, list) or not entries:
            continue
        best, best_frac = "", -1.0
        names = []
        for e in entries:
            if not isinstance(e, dict):
                if isinstance(e, str):
                    names.append(e)
                continue
            name = first_str(e, LIN_NAME_KEYS)
            if not name:
                continue
            names.append(name)
            frac = 0.0
            for fk in FRAC_KEYS:
                try:
                    frac = float(e.get(fk, 0.0))
                    break
                except (TypeError, ValueError):
                    pass
            if frac > best_frac:
                best, best_frac = name, frac
        if names:
            # main = shortest name (lineage4 is the parent of lineage4.3.4.2)
            main = min(names, key=len)
            return main, best or max(names, key=len)
    return "", ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True,
                    help="directory holding TB-Profiler *.results.json files")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pattern", default="*.results.json")
    args = ap.parse_args()

    files = sorted(Path(args.results_dir).rglob(args.pattern))
    if not files:
        raise SystemExit(f"[ERROR] no files matching {args.pattern} under {args.results_dir}")
    print(f"[INFO] {len(files):,} result files")

    rows, missing, key_use = [], 0, Counter()
    for f in files:
        try:
            obj = json.loads(f.read_text())
        except Exception as exc:  # noqa: BLE001
            print(f"[WARN] unreadable {f.name}: {exc!r}")
            continue
        sample = first_str(obj, ("id", "sample", "sample_id")) or f.name.split(".results.json")[0]
        main = first_str(obj, MAIN_KEYS)
        sub = first_str(obj, SUB_KEYS)
        if main:
            key_use["explicit"] += 1
        else:
            main, sub2 = from_list(obj)
            sub = sub or sub2
            if main:
                key_use["from_list"] += 1
        if not main:
            missing += 1
            main = "unknown"
        # mixed calls (lineage2;lineage4) are mixed infection or contamination;
        # they would otherwise form many 2-sample categories
        if ";" in main:
            coarse = "mixed"
        else:
            # coarse grouping: lineage4.3.4.2 -> lineage4
            coarse = main.split(".")[0] if main.startswith("lineage") else main
        rows.append({"sample_id": sample, "main_lineage": coarse,
                     "lineage": main, "sub_lineage": sub or main})

    # the same sample can appear under several result directories; keep one row
    by_sample: dict[str, dict] = {}
    conflicts = 0
    for r in rows:
        prev = by_sample.get(r["sample_id"])
        if prev is None:
            by_sample[r["sample_id"]] = r
        elif prev["lineage"] != r["lineage"]:
            conflicts += 1
    n_dupes = len(rows) - len(by_sample)
    rows = [by_sample[k] for k in sorted(by_sample)]
    if n_dupes:
        print(f"[INFO] collapsed {n_dupes:,} duplicate result files -> {len(rows):,} unique samples")
    if conflicts:
        print(f"[WARN] {conflicts:,} samples had conflicting lineage calls across files; kept the first")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample_id", "main_lineage", "lineage", "sub_lineage"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"[INFO] lineage source: {dict(key_use)}")
    if missing:
        print(f"[WARN] {missing:,} samples had no lineage call -> 'unknown'")
    print("\nmain_lineage distribution:")
    for lin, n in Counter(r["main_lineage"] for r in rows).most_common():
        print(f"   {lin:<16} {n:>7,}  ({100 * n / len(rows):.1f}%)")
    print(f"\n[OK] {len(rows):,} samples -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
