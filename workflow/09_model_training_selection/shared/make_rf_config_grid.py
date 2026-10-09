#!/usr/bin/env python3
"""A larger random-forest search space, sampled rather than crossed.

Section 9 settled the learner: random forest, for consistency rather than mean --
its spread across configs is the tightest of the four families (sd 0.008-0.016
against extra-trees' 0.064-0.084). That same tightness caps what a bigger search
can buy, so this exists to be thorough, not because a large gain is expected.

Sampled, not a full grid. A crossed grid over these axes spends most of its budget
on the axis with the most levels and leaves the others coarsely covered; random
search over the same budget covers each axis at the full budget's resolution,
which matters here because `max_features` is the axis that actually moves a forest
and `n_estimators` is the one that barely does.

`max_features` is log-uniform because with 49 features the interesting range is
0.03-0.6, i.e. 1 to 30 features, and uniform sampling would put half the budget
above 0.3 where the curve is flat.

Tags follow the existing convention exactly -- `rf_ne800_dNone_leaf5_sp2_mf0.178_
bs1ms0.8` -- so runs from this file and the old one are comparable by name, and a
config selected here can be found in either.

**The reason to be careful with the output.** Selecting the best of N on the
development folds sits above the truth by roughly the fold-mean's standard error
times the expected maximum of N draws: about 0.017 at N=55 and 0.025 at N=1000.
The holdout stays unbiased either way, but the config carried into it is more
likely to be a lottery winner. Use select_final_config.py's 1-SE rule rather than
the raw minimum.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import random
import sys


def tag_of(c: dict) -> str:
    d = "None" if c["max_depth"] is None else str(c["max_depth"])
    t = (f"rf_ne{c['n_estimators']}_d{d}_leaf{c['min_samples_leaf']}"
         f"_sp{c['min_samples_split']}_mf{c['max_features']}"
         f"_bs{int(c['bootstrap'])}")
    if c["bootstrap"] and c.get("max_samples"):
        t += f"ms{c['max_samples']}"
    return t


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--n-features", type=int, default=49,
                    help="only used to report how many features max_features means")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    NE = [200, 300, 400, 500, 600, 700, 800, 1000]
    DEPTH = [None, 8, 10, 12, 15, 20, 25, 30]
    LEAF = [1, 2, 3, 5, 8, 12]
    SPLIT = [2, 5, 10, 15, 20, 30]
    MS = [0.5, 0.6, 0.7, 0.8, 0.9]

    seen: set[str] = set()
    out: list[dict] = []
    guard = 0
    while len(out) < args.n and guard < args.n * 200:
        guard += 1
        boot = rng.random() < 0.5
        c = {
            "family": "rf",
            "n_estimators": rng.choice(NE),
            "max_depth": rng.choice(DEPTH),
            "min_samples_leaf": rng.choice(LEAF),
            "min_samples_split": rng.choice(SPLIT),
            # log-uniform on [0.03, 0.6]
            "max_features": round(math.exp(rng.uniform(math.log(0.03),
                                                       math.log(0.6))), 3),
            "bootstrap": boot,
        }
        if boot:
            c["max_samples"] = rng.choice(MS)
        # min_samples_split below 2*min_samples_leaf cannot bind; drop the
        # duplicate region rather than spending budget re-testing one setting
        if c["min_samples_split"] < 2 * c["min_samples_leaf"]:
            c["min_samples_split"] = 2 * c["min_samples_leaf"]
        t = tag_of(c)
        if t in seen:
            continue
        seen.add(t)
        out.append({"tag": t, **c})

    if len(out) < args.n:
        raise SystemExit(f"[ERROR] only {len(out)} unique configs from "
                         f"{guard} draws; widen the axes")

    p = pathlib.Path(args.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2) + "\n")
    mf = sorted(c["max_features"] for c in out)
    print(f"[OK] {p}  {len(out)} unique rf configs")
    print(f"     max_features {mf[0]}-{mf[-1]}  "
          f"({max(1, round(mf[0]*args.n_features))}-"
          f"{round(mf[-1]*args.n_features)} of {args.n_features} features)")
    print(f"     bootstrap: {sum(c['bootstrap'] for c in out)} of {len(out)}")
    print(f"     unlimited depth: {sum(c['max_depth'] is None for c in out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
