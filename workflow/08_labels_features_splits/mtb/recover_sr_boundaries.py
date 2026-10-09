#!/usr/bin/env python3
"""Recover each drug's S/R boundary by solving for it, not by remembering it.

`SR_BOUNDARY` is the one parameter in this pipeline that fails silently: a drug
scored at another drug's threshold produces a plausible number and no error. It
is recorded in the protocol for two drugs out of ten; the rest exist only in
submitted job scripts and logs that have since been cleaned.

It does not have to be remembered, because it is determined by data already on
disk. TB-Profiler's baseline row for a drug -- vMER and MER on a stated split --
is a function of exactly three things: the MIC per isolate, the catalogue's call
per isolate, and the boundary. Two of those are files. So every candidate
boundary between the drug's MIC classes is tried, and the one that reproduces the
accepted vMER and MER is the one that was used.

This is evidence rather than recall, and it refuses rather than guesses:

  **no candidate matches**      the accepted row was not produced from these two
                                files, or the split convention differs. Reported,
                                not papered over with a nearest fit.
  **more than one matches**     the boundary is not identified by this drug's data
                                -- possible when no isolate sits between two
                                candidates. All matches are listed; the protocol
                                should then record the interval, since any value
                                in it is equally consistent with what was run.

The split convention is `fold < 0` development, `fold >= 0` holdout, taken from
`catalogue_baseline_rows.py` rather than restated -- getting it backwards swaps
the two rows and both would still look reasonable.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

CODES = {"isoniazid": "INH", "rifampicin": "RIF", "ethambutol": "EMB",
         "ethionamide": "ETH", "amikacin": "AMI", "kanamycin": "KAN",
         "levofloxacin": "LEV", "moxifloxacin": "MXF", "bedaquiline": "BDQ",
         "delamanid": "DLM"}
TOL = 5e-4


def read_tsv(p: pathlib.Path) -> list[dict]:
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def id_col(rows: list[dict]) -> str:
    for c in rows[0]:
        if c.lower() in ("sample_id", "sample", "id"):
            return c
    raise SystemExit(f"[ERROR] no sample id column in {list(rows[0])[:5]}")


def rates(mic, call, ids, boundary):
    tr = [s for s in ids if mic[s] > boundary]
    ts = [s for s in ids if mic[s] <= boundary]
    if not tr or not ts:
        return None
    vmer = sum(1 for s in tr if call[s] == 0) / len(tr)
    mer = sum(1 for s in ts if call[s] == 1) / len(ts)
    return vmer, mer


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mic-tsv", required=True)
    ap.add_argument("--tbprofiler-binary", required=True)
    ap.add_argument("--metadata-root", required=True,
                    help="holds cv_<drug>/ with cv_folds_final.tsv and "
                         "baseline_rows.tsv")
    ap.add_argument("--fold-dir", action="append", default=[],
                    metavar="DRUG=PATH",
                    help="override for a drug whose folds are not in "
                         "cv_<drug>, e.g. isoniazid=.../cv_isoniazid_final")
    ap.add_argument("--drug", action="append", default=[],
                    help="default: all ten")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    override = dict(s.split("=", 1) for s in args.fold_dir)
    root = pathlib.Path(args.metadata_root)
    mrows = read_tsv(pathlib.Path(args.mic_tsv))
    mid = id_col(mrows)
    crows = read_tsv(pathlib.Path(args.tbprofiler_binary))
    cid = id_col(crows)

    out_rows = []
    for drug in (args.drug or sorted(CODES)):
        d = pathlib.Path(override.get(drug, str(root / f"cv_{drug}")))
        folds_p, base_p = d / "cv_folds_final.tsv", d / "baseline_rows.tsv"
        if not folds_p.exists() or not base_p.exists():
            print(f"[WARN] {drug}: missing {folds_p.name} or {base_p.name} in {d}")
            continue
        col = f"{CODES[drug]}_int"
        mic = {r[mid]: float(r[col]) for r in mrows
               if r.get(col, "") not in ("", "NA", "nan")}
        call = {r[cid]: int(float(r[drug])) for r in crows
                if r.get(drug, "") not in ("", "NA", "nan")}

        split = {}
        for r in read_tsv(folds_p):
            s = r["sample_id"]
            if s in mic and s in call:
                split[s] = "development_cv" if int(r["fold"]) < 0 else "holdout"
        want = {r["split"]: (float(r["vMER"]), float(r["MER"]))
                for r in read_tsv(base_p) if r["label"].startswith("TB-Profiler")}
        if not want:
            print(f"[WARN] {drug}: no TB-Profiler rows in {base_p.name}")
            continue

        classes = sorted({mic[s] for s in split})
        cands = [(a + b) / 2 for a, b in zip(classes, classes[1:])] + classes
        hits = []
        for b in sorted(set(cands)):
            ok = True
            for sp, (wv, wm) in want.items():
                ids = [s for s, v in split.items() if v == sp]
                got = rates(mic, call, ids, b)
                if got is None or abs(got[0] - wv) > TOL or abs(got[1] - wm) > TOL:
                    ok = False
                    break
            if ok:
                hits.append(b)

        n = {sp: sum(1 for v in split.values() if v == sp) for sp in want}
        if not hits:
            print(f"[FAIL] {drug:<14} no boundary reproduces the accepted row "
                  f"({', '.join(f'{k} n={v:,}' for k, v in n.items())})")
        elif len(hits) > 1:
            lo, hi = min(hits), max(hits)
            print(f"[AMBIG] {drug:<13} any boundary in [{lo:g}, {hi:g}] fits — "
                  f"no isolate separates them; record the interval")
            out_rows.append((drug, f"[{lo:g}, {hi:g}]", n))
        else:
            print(f"[OK]   {drug:<14} boundary {hits[0]:>6g}   "
                  f"({', '.join(f'{k} n={v:,}' for k, v in n.items())})")
            out_rows.append((drug, f"{hits[0]:g}", n))

    if args.out and out_rows:
        p = pathlib.Path(args.out)
        with p.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["drug", "sr_boundary", "n_development", "n_holdout"])
            for drug, b, n in out_rows:
                w.writerow([drug, b, n.get("development_cv", ""),
                            n.get("holdout", "")])
        print(f"\n[OK] {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
