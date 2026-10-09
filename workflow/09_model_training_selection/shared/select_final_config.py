#!/usr/bin/env python3
"""Pick the one config that goes to the holdout, by the 1-SE rule.

Taking the lowest development cost is the obvious move and the wrong one at this
search size. The quantity being minimised is a mean over 15 folds, so it carries a
standard error; the minimum of N such means sits below the truth by roughly that
error times the expected minimum of N draws -- about 0.017 at N=55 and 0.025 at
N=1000. The holdout number is unbiased regardless, since it is scored once, but
the *config* carried into it is increasingly likely to be one that won a lottery
on these particular folds.

So the standard error of the best config's per-fold cost defines a band, and every
config inside it is a candidate. For isoniazid that band held **496 of 1000**
configs, which is the real finding: cost cannot distinguish half the search space.

Something has to break that tie, and the choice matters more than it sounds --
across the tied set bal_OT ranged 0.681 to 0.773. Cost is computed at one S/R
threshold and is blind to where on the ladder a prediction lands, so leaving the
tie to "whichever happened to be lowest on cost" decides a 0.09-wide spread on a
metric that cannot see it.

Two tie-breaks are implemented and the first is the default.

  **bal_OT** — the best ordinal accuracy among the tied configs. Section 9 argued
  nothing but cost belongs in the objective because bal_OT correlates -0.85 with
  vMER across all configs; conditional on equal cost that collinearity is gone, so
  within the tie it carries independent information. For isoniazid it recovered
  MIC class 0.0 from OT 0.363 to 0.646 -- 36 isolates the minimum-cost config
  essentially could not place -- for 0.40 SE of cost. Note that bal_OT is then
  max-selected over the tied set, so its *development* value is optimistic; the
  holdout is unaffected because it is scored once.

  **regularised** — larger minimum leaf, then shallower, then fewer features per
  split. Kept for comparison, not recommended: with hundreds of configs tied it
  goes to the extreme corner of the grid, and heavily regularised forests
  under-fit an ordinal ladder while keeping S/R errors low. On isoniazid it lost
  0.023 of bal_OT and drove MIC class 5.0 to OT exactly zero.

Writes a one-config JSON, because scoring more than one config on the holdout is
the max-over-configs bias (section 3) reappearing at the last possible step.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import statistics
import sys


def cost(r, w_vmer: float, w_mer: float, pre: str = "") -> float:
    return w_vmer * float(r[pre + "vMER"]) + w_mer * float(r[pre + "MER"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True,
                    help="a cv run directory with summary.tsv and folds.tsv")
    ap.add_argument("--configs", required=True,
                    help="the JSON the run was scored with")
    ap.add_argument("--aggregation", default="median")
    ap.add_argument("--family", default="rf")
    ap.add_argument("--w-vmer", type=float, default=2.0)
    ap.add_argument("--w-mer", type=float, default=1.0)
    ap.add_argument("--se-multiple", type=float, default=1.0,
                    help="0 reproduces the plain minimum, for comparison")
    ap.add_argument("--cost-prefix", default="",
                    help="DIAGNOSTIC. 'SRt_' selects on a fitted S/R cut, which "
                         "should not be used for a reported model: snapping "
                         "already places the decision on the clinical breakpoint "
                         "deliberately, and a fitted cut would call an isolate "
                         "resistant whose predicted MIC is susceptible. Leave "
                         "empty unless you are quantifying how much of a poor cost "
                         "is miscalibration rather than discrimination")
    ap.add_argument("--objective", choices=("cost", "bal_OT"), default="cost",
                    help="what the 1-SE band is built around. 'cost' is the "
                         "panel's declared primary: minimise 2*vMER+MER. "
                         "'bal_OT' selects for the MIC ladder instead -- maximise "
                         "balanced one-tier accuracy, ties inside 1 SE broken on "
                         "lower cost. That is a SEPARATE declared analysis "
                         "answering a different question, not an alternative "
                         "headline for the same one. Choosing between the two "
                         "objectives after seeing their holdouts would be the "
                         "max-over-configs bias at the last possible step")
    ap.add_argument("--tie-break", choices=("bal_OT", "regularised", "none"),
                    default="bal_OT",
                    help="how to choose among configs the cost cannot separate")
    ap.add_argument("--tie-metric", default="bal_OT",
                    help="the column --tie-break bal_OT maximises")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    d = pathlib.Path(args.run_dir)
    S = [r for r in csv.DictReader(open(d / "summary.tsv"), delimiter="\t")
         if r["aggregation"] == args.aggregation
         and r["config"].startswith(args.family + "_")]
    if not S:
        raise SystemExit(f"[ERROR] no {args.family}/{args.aggregation} rows in {d}")
    if S and args.cost_prefix and (args.cost_prefix + "vMER") not in S[0]:
        raise SystemExit(f"[ERROR] no {args.cost_prefix}vMER column in "
                         f"{d}/summary.tsv — rerun with the current scorer")
    mean_cost = {r["config"].removesuffix("_" + args.aggregation):
                 cost(r, args.w_vmer, args.w_mer, args.cost_prefix)
                 for r in S}

    per_fold: dict[str, list[float]] = {}
    for r in csv.DictReader(open(d / "folds.tsv"), delimiter="\t"):
        if r["aggregation"] == args.aggregation:
            per_fold.setdefault(r["config"], []).append(
                cost(r, args.w_vmer, args.w_mer, args.cost_prefix))

    if args.objective == "bal_OT":
        # same 1-SE machinery, maximised. Per-fold values come from the same
        # folds.tsv, so the band means what it means for cost.
        per_fold = {}
        for r in csv.DictReader(open(d / "folds.tsv"), delimiter="\t"):
            if r["aggregation"] == args.aggregation:
                per_fold.setdefault(r["config"], []).append(
                    float(r[args.tie_metric]))
        score = {c: statistics.mean(v) for c, v in per_fold.items()}
    else:
        score = mean_cost

    lower_is_better = args.objective == "cost"
    best = (min if lower_is_better else max)(score, key=lambda c: score[c])
    folds = per_fold.get(best, [])
    if len(folds) < 3:
        raise SystemExit(f"[ERROR] only {len(folds)} folds for {best}")
    se = statistics.stdev(folds) / (len(folds) ** 0.5)
    if lower_is_better:
        ceiling = score[best] + args.se_multiple * se
        cand = [c for c, v in score.items() if v <= ceiling]
    else:
        ceiling = score[best] - args.se_multiple * se
        cand = [c for c, v in score.items() if v >= ceiling]
    print(f"[INFO] {len(score)} configs scored on {len(folds)} folds")
    print(f"[INFO] objective {args.objective}: best {score[best]:.4f}  "
          f"SE {se:.4f}  ceiling {ceiling:.4f}")
    print(f"[INFO] {len(cand)} configs within {args.se_multiple:g} SE")

    by_tag = {c["tag"]: c for c in json.loads(pathlib.Path(args.configs).read_text())}
    missing = [c for c in cand if c not in by_tag]
    if missing:
        raise SystemExit(f"[ERROR] {len(missing)} selected tags are not in "
                         f"{args.configs}, e.g. {missing[:3]}")

    def regularisation(tag: str):
        c = by_tag[tag]
        depth = c["max_depth"] if c["max_depth"] is not None else 10**6
        return (-c["min_samples_leaf"], depth, c["max_features"],
                c["n_estimators"])

    tie_val = {r["config"].removesuffix("_" + args.aggregation):
               float(r[args.tie_metric]) for r in S
               if args.tie_metric in r}
    if args.objective == "bal_OT":
        # the ladder is the objective, so the tie-break is the other metric
        chosen = min(cand, key=lambda c: mean_cost[c])
        print("[INFO] tie broken on lowest cost among the tied ladders")
    elif args.tie_break == "none":
        chosen = best
    elif args.tie_break == "regularised":
        chosen = sorted(cand, key=regularisation)[0]
    else:
        if not tie_val:
            raise SystemExit(f"[ERROR] no {args.tie_metric} column in summary.tsv")
        chosen = max(cand, key=lambda c: (tie_val.get(c, float('-inf')),
                                          -mean_cost[c]))

    span = sorted(tie_val[c] for c in cand if c in tie_val)
    print(f"[INFO] {args.tie_metric} across the tied set: "
          f"{span[0]:.3f} to {span[-1]:.3f}  "
          f"— what the cost objective leaves undetermined")
    print(f"\n[INFO] plain minimum cost : {best}  "
          f"(cost {mean_cost[best]:.4f}, {args.tie_metric} "
          f"{tie_val.get(best, float('nan')):.4f})")
    print(f"[INFO] chosen ({args.tie_break}){'':>4}: {chosen}  "
          f"(cost {mean_cost[chosen]:.4f}, {args.tie_metric} "
          f"{tie_val.get(chosen, float('nan')):.4f})")
    if chosen != best:
        dc = mean_cost[chosen] - mean_cost[best]
        db = tie_val.get(chosen, 0) - tie_val.get(best, 0)
        print(f"       {db:+.4f} {args.tie_metric} for {dc:+.4f} cost "
              f"({dc/se:.2f} SE)")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps([by_tag[chosen]], indent=2) + "\n")
    print(f"[OK] {out}\n{json.dumps(by_tag[chosen], indent=2)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
