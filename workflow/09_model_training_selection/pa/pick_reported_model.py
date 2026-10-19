#!/usr/bin/env python3
"""Which model a drug reports: the designated primary model, Evo 2.

Until 2026-10-13 this returned the model with the lowest mean DEVELOPMENT cost,
which chose Evo 2 for eight antimicrobials and Carbon for amikacin and meropenem.
The paper now reports one model throughout -- the primary model designated for
M. tuberculosis -- so the figures, the workbook, the SHAP jobs and the
study-held-out jobs, which all ask this script, answer Evo 2. The development
costs are still computed and printed with --verbose, and the earlier rule is
kept behind --by-development-cost.

The earlier rule, for the record:

Kept as a file rather than inlined in the figure driver so the figures, the
workbook and the protocol tables all answer the question the same way. Cost is
2*vMER + MER, averaged over folds of the selected config, median aggregation --
the panel's declared objective.

Never the holdout. Picking the best of three models by their held-out numbers is
the max-over-configs bias arriving one level above the configs, and it would be
invisible in the figure it produced.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

MODELS = ("evo2", "carbon", "dnabert")
DESIGNATED = "evo2"


def dev_cost(folds: pathlib.Path, config: str | None = None) -> float | None:
    rows = [r for r in csv.DictReader(open(folds), delimiter="\t")
            if r["aggregation"] == "median"
            and (config is None or r["config"] == config)]
    if not rows:
        return None
    return sum(2 * float(r["vMER"]) + float(r["MER"]) for r in rows) / len(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, help="slim/cv_<drug>")
    ap.add_argument("--tag", default="armD1k")
    ap.add_argument("--config-dir", default="",
                    help="final_configs_*/ . Without it the fold rows are used "
                         "as they stand, which is right for a slimmed tree "
                         "(already one config) and WRONG for a full run "
                         "directory, where folds.tsv holds all 1,000 configs and "
                         "the mean would be over the grid rather than over the "
                         "config that was actually selected")
    ap.add_argument("--drug", default="",
                    help="required with --config-dir, to name the config file")
    ap.add_argument("--by-development-cost", action="store_true",
                    help="the earlier rule: lowest development cost of the three")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    d = pathlib.Path(args.run_dir)
    if args.config_dir and not args.drug:
        raise SystemExit("[ERROR] --config-dir needs --drug")
    scored = {}
    for m in MODELS:
        f = d / f"runs/{m}_{args.tag}/folds.tsv"
        if not f.exists():
            continue
        cfg = None
        if args.config_dir:
            import json
            import re
            found = re.match(r"(arm[A-Z]+|reported)", args.tag)   # armCR1k -> armCR; reported1k -> reported
            arm = found.group(0) if found else "armD"
            j = pathlib.Path(args.config_dir) / \
                f"final_config_{args.drug}_{m}_{arm}.json"
            if not j.exists():
                raise SystemExit(f"[ERROR] missing {j}")
            cfg = json.load(open(j))[0]["tag"]
        c = dev_cost(f, cfg)
        if c is not None:
            scored[m] = c
    if not scored:
        raise SystemExit(f"[ERROR] no development folds under {d}")
    if args.by_development_cost:
        best = min(scored, key=scored.get)
    elif DESIGNATED in scored:
        best = DESIGNATED
    else:
        raise SystemExit(f"[ERROR] no development folds for {DESIGNATED} under {d}")
    if args.verbose:
        for m, c in sorted(scored.items(), key=lambda kv: kv[1]):
            print(f"  {m:<8} {c:.4f}{'  <-- reported' if m == best else ''}",
                  file=sys.stderr)
    print(best)
    return 0


if __name__ == "__main__":
    sys.exit(main())
