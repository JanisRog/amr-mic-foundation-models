#!/usr/bin/env python3
"""The incumbent's score on both splits, from the fold file and its own calls.

Every drug needs one row per split saying what TB-Profiler achieves on exactly
the isolates the model was scored on. Isoniazid's and ethionamide's were produced
interactively and the invocation was not recorded; reconstructing it later cost a
round of archaeology and a failed run against a `--splits` argument whose format
no surviving file matches. Ten drugs remain, so the computation is written down
here instead, taking only files the pipeline already produces.

The split definition is read from the final fold file rather than restated:
`fold < 0` is development and `fold >= 0` is the holdout, which is the convention
`--only-fold 0:0` relies on. Getting that backwards would silently swap the two
numbers, so the isolate counts are printed and should be checked against the
fold builder's own output.

Error definitions, in the direction clinical microbiology uses them:

  **vMER** — very major error: truly resistant, called susceptible. 1 − sensitivity.
  **MER**  — major error: truly susceptible, called resistant. 1 − specificity.
  **cost** — 2·vMER + MER, the pre-registered objective. The 2:1 weighting is a
  choice, not a fact, and it decides verdicts where cost and F1 disagree.

`--expect` exists because a metric script that is never checked against a known
answer is a metric script nobody should believe. Pass a previously accepted
baseline table and it verifies agreement to four decimals before writing anything.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys


def read_tsv(path: str) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def id_col(rows: list[dict], path: str) -> str:
    for c in rows[0]:
        if c.lower() in ("sample_id", "sample", "id", "sample_ids"):
            return c
    raise SystemExit(f"[ERROR] no sample id column in {path}: {list(rows[0])[:5]}")


def metrics(pairs: list[tuple[bool, bool]]) -> dict:
    """pairs of (truly_resistant, called_resistant)."""
    tp = sum(1 for t, c in pairs if t and c)
    fn = sum(1 for t, c in pairs if t and not c)
    fp = sum(1 for t, c in pairs if not t and c)
    tn = sum(1 for t, c in pairs if not t and not c)
    vmer = fn / (tp + fn) if tp + fn else 0.0
    mer = fp / (fp + tn) if fp + tn else 0.0
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return {"vMER": vmer, "MER": mer, "SR_f1": f1,
            "SR_acc": (tp + tn) / len(pairs) if pairs else 0.0,
            "cost": 2 * vmer + mer,
            "_n": len(pairs), "_R": tp + fn}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds-final", required=True,
                    help="cv_folds_final.tsv: fold < 0 is development, "
                         "fold >= 0 is the holdout")
    ap.add_argument("--calls", action="append", required=True,
                    metavar="LABEL=PATH:COLUMN",
                    help="repeatable, e.g. "
                         "'TB-Profiler 6.3.0=/path/rule_bit_rif_tbprofiler.tsv:"
                         "rule_binary_rifampicin'")
    ap.add_argument("--boundary", type=float, required=True)
    ap.add_argument("--expect", default="",
                    help="a baseline_rows.tsv whose numbers this must reproduce")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    frows = read_tsv(args.folds_final)
    fid = id_col(frows, args.folds_final)
    if "log2_mic" not in frows[0] or "fold" not in frows[0]:
        raise SystemExit(f"[ERROR] {args.folds_final} needs log2_mic and fold: "
                         f"{list(frows[0])[:6]}")
    split = {}
    mic = {}
    for r in frows:
        s = r[fid]
        split[s] = "development_cv" if int(r["fold"]) < 0 else "holdout"
        mic[s] = float(r["log2_mic"])
    n_dev = sum(1 for v in split.values() if v == "development_cv")
    print(f"[INFO] {len(split):,} isolates: {n_dev:,} development, "
          f"{len(split) - n_dev:,} holdout")
    if n_dev < len(split) - n_dev:
        print("[WARN] the development split is smaller than the holdout — "
              "fold<0 may not mean\n       development in this file. Check "
              "against the fold builder's output.")

    out_rows = []
    for spec in args.calls:
        label, _, rest = spec.partition("=")
        path, _, column = rest.rpartition(":")
        crows = read_tsv(path)
        cid = id_col(crows, path)
        if column not in crows[0]:
            raise SystemExit(f"[ERROR] no {column!r} in {path}: "
                             f"{list(crows[0])[:6]}")
        call = {r[cid]: int(float(r[column])) for r in crows
                if r.get(column, "") not in ("", "NA", "nan")}
        for sp in ("development_cv", "holdout"):
            pairs = [(mic[s] > args.boundary, call[s] == 1)
                     for s in split if split[s] == sp and s in call]
            m = metrics(pairs)
            print(f"[INFO] {label:<20}{sp:<16}{m['_n']:>6,} isolates, "
                  f"{m['_R']:>5,} resistant   cost {m['cost']:.6f}")
            out_rows.append({"label": label, "split": sp,
                             **{k: f"{v:.6f}" for k, v in m.items()
                                if not k.startswith("_")}})

    if args.expect:
        want = {(r["label"], r["split"]): r for r in read_tsv(args.expect)}
        bad = []
        for r in out_rows:
            w = want.get((r["label"], r["split"]))
            if not w:
                continue
            for k in ("vMER", "MER", "SR_f1", "SR_acc", "cost"):
                if k in w and abs(float(w[k]) - float(r[k])) > 5e-5:
                    bad.append(f"{r['label']}/{r['split']}/{k}: "
                               f"expected {float(w[k]):.6f}, got {float(r[k]):.6f}")
        checked = sum(1 for r in out_rows if (r["label"], r["split"]) in want)
        if bad:
            raise SystemExit("[ERROR] does not reproduce the accepted baseline:\n  "
                             + "\n  ".join(bad))
        print(f"[OK] reproduces {checked} accepted row(s) in {args.expect}")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["label", "split", "vMER", "MER",
                                       "SR_f1", "SR_acc", "cost"])
        w.writeheader()
        w.writerows(out_rows)
    print(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
