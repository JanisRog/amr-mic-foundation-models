#!/usr/bin/env python3
"""ARDaP's own per-drug call as one feature column, optionally after another block.

The P. aeruginosa counterpart of the TB-Profiler call in the M. tuberculosis model
(exploratory arms CR and DR). The value is ARDaP's Status as written by
gather_ardap_calls.py, encoded S = 0, I = 1, R = 2: the tool's own output with
nothing derived and nothing lost, so the forest can treat Intermediate either way.

Every isolate of the cohort (cv_folds_final.tsv) keeps its row. An isolate ARDaP
gave no call (absent from the calls file, or NA) is coded -1, "no call": filling it
with S would invent a call ARDaP did not make, and dropping it would change the
cohort the arms are compared on. Allowed only up to --max-no-call (default 1%) of
the cohort; more stops the build, since it would point at a failed ARDaP run.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

CODE = {"S": 0, "I": 1, "R": 2}
NO_CALL = -1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="ardap_calls.tsv (gather_ardap_calls.py)")
    ap.add_argument("--drug", required=True)
    ap.add_argument("--rows", required=True, help="cv_folds_final.tsv: every isolate needs a call")
    ap.add_argument("--merge-with", default="", help="feature TSV to keep in front (arm DR: presence)")
    ap.add_argument("--max-no-call", type=float, default=0.01,
                    help="largest fraction of the cohort allowed without a call")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    calls = {r["sample_id"]: r for r in csv.DictReader(open(a.calls), delimiter="\t")}
    if a.drug not in next(iter(calls.values())):
        raise SystemExit(f"[ERROR] no {a.drug!r} column in {a.calls}: ARDaP does not cover it")
    ids = list(dict.fromkeys(r["sample_id"] for r in csv.DictReader(open(a.rows), delimiter="\t")))
    value = {s: CODE.get((calls.get(s) or {}).get(a.drug, ""), NO_CALL) for s in ids}
    nocall = [s for s in ids if value[s] == NO_CALL]
    if len(nocall) > a.max_no_call * len(ids):
        raise SystemExit(f"[ERROR] {a.drug}: {len(nocall)} of {len(ids)} isolates have no ARDaP call "
                         f"(limit {a.max_no_call:.0%}), e.g. {nocall[0]}")

    pre_cols, pre = [], {}
    if a.merge_with:
        rows = list(csv.DictReader(open(a.merge_with), delimiter="\t"))
        idc = list(rows[0])[0]
        pre_cols = [c for c in rows[0] if c != idc]
        pre = {r[idc]: r for r in rows}
        miss = [s for s in ids if s not in pre]
        if miss:
            raise SystemExit(f"[ERROR] {len(miss)} isolate(s) missing from {a.merge_with}, e.g. {miss[0]}")

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id"] + pre_cols + [f"ardap_call_{a.drug}"])
        for s in ids:
            w.writerow([s] + [pre[s][c] for c in pre_cols] + [value[s]])
    n = {k: sum(1 for s in ids if value[s] == v) for k, v in CODE.items()}
    print(f"[OK] {out}: {len(ids)} isolates, S {n['S']} / I {n['I']} / R {n['R']} / no call {len(nocall)}"
          + (f" ({', '.join(nocall[:3])}{' ...' if len(nocall) > 3 else ''})" if nocall else "")
          + (f", after {len(pre_cols)} column(s) of {pathlib.Path(a.merge_with).name}" if pre_cols else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
