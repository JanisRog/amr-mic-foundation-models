#!/usr/bin/env python3
"""Per-isolate SHAP values for the reported model, summed to one per target.

The model's features are four principal components per target plus the rule bit,
so a feature-level attribution says "PC2 of Rv1908c:gene mattered", which is not
a statement about biology. Summing a target's components gives the quantity that
is interpretable -- SHAP is additive in the features, so the sum over a target's
columns is exactly that target's contribution to the prediction, with no
approximation introduced by the grouping.

**This refits the reported model rather than loading it.** The fit is
deterministic given the cache, the config and the seed, and the same refit has
already been shown to reproduce the reported holdout row to five decimals. The
predictions are written alongside the SHAP values so the plotting step can check
that again before drawing -- an attribution for a model that is not the reported
one would be worse than no attribution.

Targets keep their `<locus>:<region>` ids, so gene body and upstream window stay
separate. They are separate features; collapsing them here would hide the one
distinction this panel exists to make -- for isoniazid the promoter windows carry
much of the catalogue's signal.

TreeSHAP on a deep forest can be slow. The job times a small batch first and
prints the projected total before committing to it, so an overnight surprise
becomes a decision made in the first minute.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from cv_run_architecture import (fixed_sr_thresholds, make_model, slice_columns,
                                 snap)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", required=True, help="the HOLDOUT cache")
    ap.add_argument("--configs", required=True, help="one-config JSON")
    ap.add_argument("--extra-feature-tsv", default="")
    ap.add_argument("--mode", default="shared4")
    ap.add_argument("--top-k", type=int, required=True)
    ap.add_argument("--sr-boundary", type=float, required=True)
    ap.add_argument("--clumps", default="",
                    help="optional TSV with reference_gene_id/_symbol, to label "
                         "loci with gene names")
    ap.add_argument("--max-explain", type=int, default=0,
                    help="subsample the explained isolates; 0 explains all")
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    import shap

    if args.clumps and not pathlib.Path(args.clumps).exists():
        print(f"[WARN] no clumps file at {args.clumps} — targets will be "
              f"labelled by locus tag only")
        args.clumps = ""

    cache = pathlib.Path(args.cache_dir)
    folds = sorted(cache.glob("fold_r*_f*.npz"))
    if len(folds) != 1:
        raise SystemExit(f"[ERROR] expected one fold file in {cache}, found "
                         f"{len(folds)}; this is meant for a holdout cache")
    d = np.load(folds[0], allow_pickle=True)
    cols = slice_columns(d["n_pcs"], d["pc_mask"], args.top_k, args.mode)
    x = d["x"][:, cols]
    names = [str(n) for n in np.asarray(d["feature_names"])[cols]]
    y, is_train = d["y"], d["is_train"]
    samples = np.array([str(s) for s in d["samples"]])

    if args.extra_feature_tsv:
        e = pd.read_csv(args.extra_feature_tsv, sep="\t")
        e[e.columns[0]] = e[e.columns[0]].astype(str)
        e = e.set_index(e.columns[0])
        miss = [s for s in samples if s not in e.index]
        if miss:
            raise SystemExit(f"[ERROR] {len(miss):,} cache samples absent from "
                             f"the extra-feature TSV, e.g. {miss[0]!r}")
        x = np.hstack([x, e.loc[list(samples)].to_numpy(np.float32)])
        names += [f"rule::{c}" for c in e.columns]
    print(f"[INFO] {x.shape[0]:,} isolates x {x.shape[1]} features "
          f"({int(is_train.sum()):,} train / {int((~is_train).sum()):,} holdout)")

    cfg = json.load(open(args.configs))
    if isinstance(cfg, list):
        if len(cfg) != 1:
            raise SystemExit(f"[ERROR] {len(cfg)} configs in {args.configs}")
        cfg = cfg[0]
    print(f"[INFO] config {cfg['tag']}")
    model = make_model(cfg, args.seed, args.n_jobs)
    t0 = time.time()
    model.fit(x[is_train], y[is_train])
    print(f"[INFO] fitted in {time.time() - t0:.0f}s")

    xte, yte = x[~is_train], y[~is_train]
    ste = samples[~is_train]
    if args.max_explain and args.max_explain < len(xte):
        rng = np.random.default_rng(args.seed)
        keep = np.sort(rng.choice(len(xte), args.max_explain, replace=False))
        xte, yte, ste = xte[keep], yte[keep], ste[keep]
        print(f"[WARN] explaining a {len(xte):,}-isolate subsample, not the "
              f"whole holdout — the beeswarm is then a sample of the holdout")

    # Predictions first: they are what lets the plot verify this refit is the
    # reported model before it draws an attribution for it.
    allowed = sorted(np.unique(y).tolist())
    thresholds = fixed_sr_thresholds(allowed, args.sr_boundary)
    per_tree = np.vstack([t.predict(xte) for t in model.estimators_])
    out = pathlib.Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for agg, raw in (("mean", per_tree.mean(axis=0)),
                     ("median", np.median(per_tree, axis=0))):
        pred = snap(raw, allowed, thresholds)
        rows += [{"sample_id": s, "fold": folds[0].stem, "config": cfg["tag"],
                  "aggregation": agg, "log2_mic": float(t), "raw": float(r),
                  "pred": float(p)}
                 for s, t, r, p in zip(ste, yte, raw, pred)]
    pd.DataFrame(rows).to_csv(f"{out}_predictions.tsv", sep="\t", index=False)
    print(f"[OK] {out}_predictions.tsv")

    expl = shap.TreeExplainer(model)
    probe = min(40, len(xte))
    t0 = time.time()
    expl.shap_values(xte[:probe], check_additivity=False)
    rate = (time.time() - t0) / probe
    print(f"[INFO] TreeSHAP {rate * 1000:.0f} ms/isolate — "
          f"{rate * len(xte) / 60:.1f} min projected for {len(xte):,}")
    t0 = time.time()
    sv = expl.shap_values(xte, check_additivity=False)
    print(f"[INFO] done in {(time.time() - t0) / 60:.1f} min")
    sv = np.asarray(sv)
    if sv.ndim != 2 or sv.shape != xte.shape:
        raise SystemExit(f"[ERROR] unexpected SHAP shape {sv.shape} for "
                         f"x {xte.shape}")

    # Feature -> target. SHAP is additive, so a target's contribution is exactly
    # the sum of its columns; no weighting or normalisation is involved.
    groups: dict[str, list[int]] = {}
    for j, n in enumerate(names):
        t = n.rsplit("_pc", 1)[0] if "_pc" in n else n
        groups.setdefault(t, []).append(j)
    print(f"[INFO] {len(groups)} targets from {len(names)} features")

    symbol = {}
    if args.clumps:
        try:
            with open(args.clumps) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    gid = r.get("reference_gene_id") or r.get("locus_tag")
                    if gid:
                        symbol[gid] = (r.get("reference_gene_symbol")
                                       or r.get("gene") or "")
            print(f"[INFO] {len(symbol):,} gene symbols from "
                  f"{pathlib.Path(args.clumps).name}")
        except (OSError, csv.Error, KeyError) as e:
            print(f"[WARN] could not read symbols from {args.clumps}: {e}")

    long_rows = []
    for t, js in groups.items():
        locus = t.split(":")[0]
        # PA carries two kinds MTB does not: acquired-gene sequence targets
        # (amrseq_<family>) and a presence block of many columns rather than one
        # rule bit. Both get their own category so the beeswarm can say which
        # *kind* of feature is doing the work, which is the arm question.
        region = ("acquired" if t.startswith("amrseq_")
                  else "other" if t.startswith("rule::")
                  else "upstream" if "upstream" in t
                  else "gene")
        g = symbol.get(locus, "")
        for s, v in zip(ste, sv[:, js].sum(axis=1)):
            long_rows.append({"sample_id": s, "target_id": t, "locus": locus,
                              "gene": g, "region": region, "shap": float(v)})
    long = pd.DataFrame(long_rows)
    long.to_csv(f"{out}_shap_long.tsv", sep="\t", index=False)

    summ = (long.groupby(["target_id", "locus", "gene", "region"], as_index=False)
                .agg(mean_abs_shap=("shap", lambda v: float(np.abs(v).mean())),
                     mean_shap=("shap", "mean"))
                .sort_values("mean_abs_shap", ascending=False))
    summ.to_csv(f"{out}_target_summary.tsv", sep="\t", index=False)
    print(f"[OK] {out}_shap_long.tsv   ({len(long):,} rows)")
    print(f"[OK] {out}_target_summary.tsv\n")
    print(summ.head(15).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
