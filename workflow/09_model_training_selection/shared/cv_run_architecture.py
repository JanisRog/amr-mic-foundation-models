#!/usr/bin/env python3
"""One architecture, all 15 partitions, one log.

Reads the cached folds from cv_precompute_fold_features.py, slices them to this
architecture's target count and PC selection, and scores every RF config across
the repeated CV. Metrics are averaged over the 15 partitions, so `summary.tsv`
carries one row per (config, aggregation) rather than per fold.

Mean and median tree aggregation come from the *same fitted forest* -- median is a
different way of reading the trees, not a different model -- so 50 configs give
100 scored rows for the cost of 50 fits. Doing it as two runs would double the
compute for nothing.

The three PCA modes differ only in what is sliced out of the cache:

  shared4   the shared-PCA cache, 4 PCs per target
  ptvar     the per-target-variance cache, every PC it kept (99% var, cap 25)
  ptmodel   the same cache with the per-PC importance mask applied, so each
            target contributes only the PCs a model actually used

Per-class accuracy is written alongside, because bal_OT is effectively decided by
the rare classes (section 3) and a headline number hides which class moved.
"""

from __future__ import annotations

import argparse
import inspect
import json
import pathlib
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import (ExtraTreesRegressor, HistGradientBoostingRegressor,
                              RandomForestRegressor)
from sklearn.linear_model import Ridge
from sklearn.metrics import (f1_score, mean_absolute_error, mean_squared_error,
                             r2_score, roc_auc_score)


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def balanced_recall(y_true: np.ndarray, ok_mask: np.ndarray,
                    min_n: int = 0) -> float:
    """Mean per-class recall, optionally over classes with enough support.

    With 11 MIC classes each carries ~9% of bal_OT regardless of how many
    isolates it holds, so a class with one isolate per fold injects a Bernoulli
    draw worth ~0.09 of the metric. min_n drops those, giving a companion number
    that moves only when something real moves.
    """
    values = [float(np.mean(ok_mask[y_true == c])) for c in np.unique(y_true)
              if int((y_true == c).sum()) >= min_n]
    return float(np.mean(values)) if values else float("nan")


def fixed_sr_thresholds(allowed: list[float], boundary: float) -> np.ndarray:
    values = np.asarray(sorted(allowed), dtype=float)
    thresholds = (values[:-1] + values[1:]) / 2.0
    candidate = np.flatnonzero((values[:-1] <= boundary) & (values[1:] > boundary))
    if len(candidate) != 1:
        raise ValueError(f"no single S/R interval for boundary {boundary:g} "
                         f"in {allowed}")
    thresholds[candidate[0]] = boundary
    return thresholds


def snap(raw: np.ndarray, allowed: list[float], thresholds: np.ndarray) -> np.ndarray:
    values = np.asarray(sorted(allowed), dtype=float)
    return values[np.searchsorted(thresholds, raw, side="left")]


def tuned_cut(raw_tr: np.ndarray, y_tr: np.ndarray, boundary: float,
              w_vmer: float = 2.0, w_mer: float = 1.0) -> float:
    """DIAGNOSTIC ONLY. The cut that would minimise the objective, from training.

    Do not select or report on this. `fixed_sr_thresholds` already places the S/R
    decision deliberately: snapping uses midpoints between MIC classes everywhere
    except the interval containing the breakpoint, where the threshold is the
    breakpoint itself, so a raw prediction of 2.01 snaps up to 3.0 and is called
    resistant rather than rounding down to susceptible. The decision already sits
    where the clinic puts it and already leans toward sensitivity there.

    Moving it further breaks the model's own coherence. At a fitted cut of 1.28 an
    isolate whose predicted MIC snaps to 1.0 would be called resistant -- a
    susceptible MIC and a resistant call for the same isolate, which is
    indefensible for a model whose entire claim is that it predicts an MIC.

    What it is good for is separating two explanations of a poor cost. On
    ethionamide the fixed cut gives vMER 0.186 against MER 0.047 while TB-Profiler
    sits at 0.122 / 0.110, and cost 0.419 against 0.353. The fitted cut reaches
    0.315 on the same predictions, with AUROC 0.934 -- so the ranking is sound and
    the loss is the regression shrinking toward the mean when the resistant class
    is 13% of the cohort, under-predicting the tail. That is a calibration failure
    in the predicted MICs, which also drags MAE and bal_OT and which relocating a
    threshold cannot fix.
    """
    r = y_tr > boundary
    if r.all() or not r.any():
        return boundary
    n_r, n_s = int(r.sum()), int((~r).sum())
    cand = np.unique(np.concatenate([raw_tr, [boundary]]))
    best, best_cost = boundary, float("inf")
    for c in cand:
        pred = raw_tr > c
        cost = (w_vmer * (r & ~pred).sum() / n_r
                + w_mer * (~r & pred).sum() / n_s)
        if cost < best_cost:
            best, best_cost = float(c), cost
    return best


def metrics_for(y_true: np.ndarray, y_pred: np.ndarray, raw: np.ndarray,
                boundary: float, min_class_n: int = 0,
                cut: float | None = None) -> dict:
    err = np.abs(y_true - y_pred)
    true_r, pred_r = y_true > boundary, y_pred > boundary
    n_r, n_s = int(true_r.sum()), int((~true_r).sum())
    false_s = int((true_r & ~pred_r).sum())
    false_r = int((~true_r & pred_r).sum())
    auroc = float("nan")
    if len(np.unique(true_r)) == 2:
        auroc = float(roc_auc_score(true_r.astype(int), raw))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        rho = spearmanr(y_true, y_pred).correlation
    return {"exact": float(np.mean(err == 0)),
            "exact_bal": balanced_recall(y_true, err == 0),
            "OT": float(np.mean(err <= 1)),
            "bal_OT": balanced_recall(y_true, err <= 1),
            "bal_OT_common": balanced_recall(y_true, err <= 1, min_class_n),
            "MAE": float(mean_absolute_error(y_true, y_pred)),
            "RMSE": float(np.sqrt(mean_squared_error(y_true, y_pred))),
            "R2": float(r2_score(y_true, y_pred)),
            "spearman": float(rho) if rho == rho else float("nan"),
            "vMER": float(false_s / n_r) if n_r else float("nan"),
            "MER": float(false_r / n_s) if n_s else float("nan"),
            "SR_acc": float(np.mean(true_r == pred_r)),
            "SR_f1": float(f1_score(true_r.astype(int), pred_r.astype(int),
                                    zero_division=0)),
            "SR_AUROC_raw": auroc,
            **_tuned_metrics(y_true, raw, boundary, cut)}


def _tuned_metrics(y_true: np.ndarray, raw: np.ndarray, boundary: float,
                   cut: float | None) -> dict:
    """The same S/R numbers at a cut fitted on training rows, reported beside the
    fixed-boundary ones so the two are visible in one table."""
    if cut is None:
        return {}
    true_r, pred_r = y_true > boundary, raw > cut
    n_r, n_s = int(true_r.sum()), int((~true_r).sum())
    return {
        "SRt_cut": float(cut),
        "SRt_vMER": float((true_r & ~pred_r).sum() / n_r) if n_r else float("nan"),
        "SRt_MER": float((~true_r & pred_r).sum() / n_s) if n_s else float("nan"),
        "SRt_acc": float(np.mean(true_r == pred_r)),
        "SRt_f1": float(f1_score(true_r.astype(int), pred_r.astype(int),
                                 zero_division=0)),
    }


_HGB_PARAMS = set(inspect.signature(HistGradientBoostingRegressor).parameters)


def make_model(cfg: dict, seed: int, n_jobs: int):
    """One config, any family.

    The grid was built to compare representations, which silently assumed random
    forests are the right learner for PCA'd embedding blocks. They may not be:
    extremely randomised trees often do better on noisy high-dimensional inputs,
    boosting handles many weak features differently, and a ridge tells us whether
    any of the ensembling is earning its keep at all. Family is a config field so
    all of them are scored on identical folds.
    """
    fam = cfg.get("family", "rf")
    if fam in ("rf", "et"):
        kw = dict(n_estimators=int(cfg["n_estimators"]),
                  max_depth=cfg.get("max_depth"),
                  min_samples_leaf=int(cfg.get("min_samples_leaf", 1)),
                  min_samples_split=int(cfg.get("min_samples_split", 2)),
                  max_features=cfg.get("max_features"),
                  bootstrap=bool(cfg.get("bootstrap", False)),
                  random_state=seed, n_jobs=n_jobs)
        if kw["bootstrap"] and cfg.get("max_samples"):
            kw["max_samples"] = float(cfg["max_samples"])
        return (RandomForestRegressor if fam == "rf" else ExtraTreesRegressor)(**kw)
    if fam == "hgb":
        kw = dict(max_iter=int(cfg.get("max_iter", 300)),
                  learning_rate=float(cfg.get("learning_rate", 0.1)),
                  max_depth=cfg.get("max_depth"),
                  min_samples_leaf=int(cfg.get("min_samples_leaf", 20)),
                  l2_regularization=float(cfg.get("l2", 0.0)),
                  early_stopping=False, random_state=seed)
        # max_features arrived in scikit-learn 1.5; passing it to an older build
        # is a TypeError twelve hours into an overnight grid
        if "max_features" in _HGB_PARAMS:
            kw["max_features"] = float(cfg.get("max_features", 1.0))
        return HistGradientBoostingRegressor(**kw)
    if fam == "ridge":
        # Ridge takes random_state only for the sag/saga solvers, and the default
        # solver ignores it in some versions and rejects it in others
        return Ridge(alpha=float(cfg.get("alpha", 1.0)))
    raise SystemExit(f"[ERROR] unknown family {fam!r} in config {cfg.get('tag')}")


def collect_predictions(rows, ids, y_true, raw, pred, config, agg, fold):
    """One row per held-out isolate: truth, raw prediction, snapped class."""
    if ids is None:
        raise SystemExit("[ERROR] --save-predictions needs the cache's sample "
                         "ids; this code path did not pass them")
    if not (len(ids) == len(y_true) == len(raw) == len(pred)):
        raise SystemExit(f"[ERROR] prediction arrays disagree in length: "
                         f"{len(ids)} ids, {len(y_true)} truths, {len(raw)} raw, "
                         f"{len(pred)} snapped")
    rows.extend({"sample_id": str(i), "fold": fold, "config": config,
                 "aggregation": agg, "log2_mic": float(t), "raw": float(r),
                 "pred": float(p)}
                for i, t, r, p in zip(ids, y_true, raw, pred))


def write_predictions(rows, path: str) -> None:
    if not path:
        return
    out = pathlib.Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    log(f"[OK] {out}: {len(df):,} predictions over "
        f"{df['sample_id'].nunique():,} isolates, "
        f"{df['config'].nunique()} config(s) x {df['aggregation'].nunique()} "
        f"aggregation(s)")


def slice_columns(n_pcs: np.ndarray, pc_mask: np.ndarray, top_k: int,
                  mode: str) -> np.ndarray:
    """Column indices for the first top_k targets under this PCA mode."""
    ends = np.cumsum(n_pcs)
    starts = ends - n_pcs
    keep: list[np.ndarray] = []
    for t in range(min(top_k, len(n_pcs))):
        cols = np.arange(starts[t], ends[t])
        if mode == "ptmodel":
            cols = cols[pc_mask[starts[t]:ends[t]]]
        keep.append(cols)
    return np.concatenate(keep) if keep else np.array([], dtype=int)


def score_fold(xtr, ytr, xte, yte, configs, seed, n_jobs, boundary,
               min_class_n, label, fold_rows, class_rows, fold_name, n_feat,
               pred_rows=None, val_ids=None):
    """Fit every config on one fold and record both aggregations.

    `pred_rows` collects one row per held-out isolate per config: the true MIC,
    the continuous prediction and the snapped class. The summary metrics cannot
    be inverted back to that -- a confusion plot, a per-isolate error listing or
    any question of the form "which isolates does it miss" needs the predictions
    themselves, and re-fitting to recover them invites a silent drift from the
    numbers already reported. Off unless asked for: it is one row per isolate per
    config, which for a 1,000-config sweep is millions of rows.
    """
    allowed = sorted(np.unique(np.concatenate([ytr, yte])).tolist())
    thresholds = fixed_sr_thresholds(allowed, boundary)
    for cfg in configs:
        model = make_model(cfg, seed, n_jobs)
        model.fit(xtr, ytr)
        if cfg.get("family", "rf") in ("rf", "et"):
            per_tree = np.vstack([t.predict(xte) for t in model.estimators_])
            tr_tree = np.vstack([t.predict(xtr) for t in model.estimators_])
            preds = (("mean", per_tree.mean(axis=0), tr_tree.mean(axis=0)),
                     ("median", np.median(per_tree, axis=0),
                      np.median(tr_tree, axis=0)))
        else:
            preds = (("mean", model.predict(xte), model.predict(xtr)),)
        for agg, raw, raw_tr in preds:
            pred = snap(raw, allowed, thresholds)
            m = metrics_for(yte, pred, raw, boundary, min_class_n,
                            cut=tuned_cut(raw_tr, ytr, boundary))
            m.update(config=cfg["tag"], aggregation=agg, fold=fold_name,
                     n_features=n_feat)
            fold_rows.append(m)
            if pred_rows is not None:
                collect_predictions(pred_rows, val_ids, yte, raw, pred,
                                    cfg["tag"], agg, fold_name)
            err = np.abs(yte - pred)
            for c in allowed:
                mask = yte == c
                if mask.any():
                    class_rows.append({"config": cfg["tag"], "aggregation": agg,
                                       "fold": fold_name, "log2_mic": c,
                                       "n": int(mask.sum()),
                                       "exact": float(np.mean(err[mask] == 0)),
                                       "OT": float(np.mean(err[mask] <= 1))})


def run_feature_baseline(args) -> int:
    """A baseline on the same folds: features come from a matrix, not a cache.

    The point of the split is which metric each baseline answers for. TB-Profiler's
    drug-level call is the incumbent for vMER and MER and needs no model at all. For
    the MIC-ladder metrics the fair comparator is a model over the catalogue's
    mutation columns, which can express that katG S315T is high-level and an inhA
    promoter mutation is not. Catalogue-as-MIC gives every resistant isolate one
    value, so it cannot place anything on the ladder and beating it proves nothing.
    """
    if not args.folds_tsv:
        raise SystemExit("[ERROR] --feature-tsv needs --folds-tsv")
    folds = pd.read_csv(args.folds_tsv, sep="\t")
    folds["sample_id"] = folds["sample_id"].astype(str)
    feat = pd.read_csv(args.feature_tsv, sep="\t")
    feat[feat.columns[0]] = feat[feat.columns[0]].astype(str)
    feat = feat.set_index(feat.columns[0])
    samples = [s for s in sorted(folds["sample_id"].unique()) if s in feat.index]
    missing = folds["sample_id"].nunique() - len(samples)
    log(f"baseline {args.label}: {feat.shape[1]:,} columns, {len(samples):,} of "
        f"{folds['sample_id'].nunique():,} fold samples present"
        + (f" ({missing:,} absent from the matrix)" if missing else ""))
    if missing:
        log("  absent samples are dropped, so this baseline is scored on a subset "
            "— report the n alongside it")

    x_all = feat.loc[samples].to_numpy(np.float32)
    y_map = folds.drop_duplicates("sample_id").set_index("sample_id")["log2_mic"]
    y_all = y_map.reindex(samples).to_numpy(float)
    configs = json.loads(pathlib.Path(args.configs).read_text())
    if args.limit_configs:
        configs = configs[:args.limit_configs]

    out = pathlib.Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    fold_rows, class_rows = [], []
    pred_rows = [] if args.save_predictions else None
    pairs = sorted({(int(r), int(f)) for r, f in zip(folds["repeat"], folds["fold"])
                    if int(f) >= 0})
    for i, (rep, fold) in enumerate(pairs, 1):
        val = set(folds.loc[(folds["repeat"] == rep) & (folds["fold"] == fold),
                            "sample_id"])
        is_train = np.array([s not in val for s in samples])
        val_ids = [s for s, tr in zip(samples, is_train) if not tr]
        xtr, xte = x_all[is_train], x_all[~is_train]
        # a column constant in training carries nothing and can only add noise
        keep = xtr.std(axis=0) > 0
        log(f"[{i}/{len(pairs)}] r{rep}_f{fold}  train {int(is_train.sum()):,}  "
            f"val {int((~is_train).sum()):,}  features {int(keep.sum()):,}")
        score_fold(xtr[:, keep], y_all[is_train], xte[:, keep], y_all[~is_train],
                   configs, args.seed + i, args.n_jobs, args.sr_boundary,
                   args.min_class_n, args.label, fold_rows, class_rows,
                   f"fold_r{rep}_f{fold}", int(keep.sum()),
                   pred_rows=pred_rows, val_ids=val_ids)

    write_predictions(pred_rows, args.save_predictions)
    df = pd.DataFrame(fold_rows)
    df.to_csv(out / "folds.tsv", sep="\t", index=False)
    cols = [c for c in df.columns if c not in ("config", "aggregation", "fold")]
    summary = df.groupby(["config", "aggregation"], as_index=False)[cols].mean()
    summary["config"] = summary["config"] + "_" + summary["aggregation"]
    summary.sort_values("bal_OT", ascending=False).to_csv(
        out / "summary.tsv", sep="\t", index=False)
    pd.DataFrame(class_rows).groupby(
        ["config", "aggregation", "log2_mic"], as_index=False)[
        ["n", "exact", "OT"]].mean().to_csv(out / "per_class.tsv", sep="\t",
                                            index=False)
    best = summary.sort_values("bal_OT", ascending=False).iloc[0]
    log(f"\n{args.label}: bal_OT {best['bal_OT']:.3f}  OT {best['OT']:.3f}  "
        f"MAE {best['MAE']:.3f}  vMER {best['vMER']:.3f}  MER {best['MER']:.3f}")
    log(f"[OK] {out}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default="")
    ap.add_argument("--feature-tsv", default="",
                    help="score a baseline instead of an embedding cache: a "
                         "sample_id + feature-columns matrix, e.g. the "
                         "mutation-binary or allele-fraction catalogue features. "
                         "Needs --folds-tsv. Nothing is fitted per fold except the "
                         "drop of columns constant in that fold's training rows")
    ap.add_argument("--folds-tsv", default="",
                    help="required with --feature-tsv; the same fold file the "
                         "embedding arms used, so the comparison is paired")
    ap.add_argument("--save-predictions", default="",
                    help="also write per-isolate predictions to this TSV: "
                         "sample_id, fold, config, aggregation, log2_mic, raw, "
                         "pred. One row per held-out isolate per config per "
                         "aggregation, so pass a short --configs when using it")
    ap.add_argument("--extra-feature-tsv", default="",
                    help="append these columns to an embedding cache, e.g. the "
                         "mutation-binary matrix. The question it answers is "
                         "whether embeddings *add* to a catalogue, which is "
                         "different from whether they beat one — and section 2 "
                         "says they do not beat one for isoniazid")
    ap.add_argument("--mode", choices=("shared4", "ptvar", "ptmodel"),
                    default="shared4",
                    help="how to slice an embedding cache; ignored with "
                         "--feature-tsv, whose columns are already final")
    ap.add_argument("--top-k", type=int, default=64,
                    help="targets to keep from the cache; ignored with "
                         "--feature-tsv")
    ap.add_argument("--configs", required=True)
    ap.add_argument("--sr-boundary", type=float, required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--min-class-n", type=int, default=10,
                    help="classes with fewer validation isolates than this are "
                         "excluded from bal_OT_common; bal_OT itself keeps them")
    ap.add_argument("--limit-configs", type=int, default=0)
    ap.add_argument("--limit-folds", type=int, default=0)
    args = ap.parse_args()

    if args.feature_tsv:
        return run_feature_baseline(args)
    if not args.cache_dir:
        raise SystemExit("[ERROR] give --cache-dir or --feature-tsv")
    if "--mode" not in sys.argv or "--top-k" not in sys.argv:
        log("[WARN] --mode/--top-k not given; defaulting to "
            f"{args.mode} / {args.top_k} for this cache")
    cache = pathlib.Path(args.cache_dir)
    manifest = json.loads((cache / "precompute_manifest.json").read_text())
    if args.mode == "shared4" and manifest["pca_mode"] != "shared4":
        raise SystemExit(f"[ERROR] mode shared4 needs a shared4 cache, "
                         f"{cache} is {manifest['pca_mode']}")
    if args.mode in ("ptvar", "ptmodel") and manifest["pca_mode"] != "ptvar":
        raise SystemExit(f"[ERROR] mode {args.mode} needs a ptvar cache, "
                         f"{cache} is {manifest['pca_mode']}")
    if args.top_k > manifest["top_targets"]:
        raise SystemExit(f"[ERROR] cache holds {manifest['top_targets']} targets, "
                         f"--top-k {args.top_k} asks for more")

    configs = json.loads(pathlib.Path(args.configs).read_text())
    if args.limit_configs:
        configs = configs[:args.limit_configs]
    fold_files = sorted(cache.glob("fold_r*_f*.npz"))
    if args.limit_folds:
        fold_files = fold_files[:args.limit_folds]
    if not fold_files:
        raise SystemExit(f"[ERROR] no fold caches in {cache}")

    log(f"architecture {args.label}")
    log(f"  cache {cache}  mode {args.mode}  top_k {args.top_k}")
    log(f"  {len(configs)} configs x {len(fold_files)} partitions x 2 aggregations")

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fold_rows, class_rows = [], []
    pred_rows = [] if args.save_predictions else None

    extra = None
    if args.extra_feature_tsv:
        e = pd.read_csv(args.extra_feature_tsv, sep="\t")
        e[e.columns[0]] = e[e.columns[0]].astype(str)
        extra = e.set_index(e.columns[0])
        log(f"appending {extra.shape[1]:,} columns from "
            f"{pathlib.Path(args.extra_feature_tsv).name}")

    for fi, fpath in enumerate(fold_files):
        t0 = time.time()
        d = np.load(fpath, allow_pickle=True)
        sample_ids = np.array([str(v) for v in d["samples"]])
        cols = slice_columns(d["n_pcs"], d["pc_mask"], args.top_k, args.mode)
        x = d["x"][:, cols]
        y, is_train = d["y"], d["is_train"]
        if extra is not None:
            ids = list(sample_ids)
            miss = [s for s in ids if s not in extra.index]
            if miss:
                raise SystemExit(f"[ERROR] {len(miss):,} cache samples are absent "
                                 f"from --extra-feature-tsv, e.g. {miss[0]!r}; the "
                                 f"combined arm must cover the same isolates as "
                                 f"the embedding arm it is compared against")
            x = np.hstack([x, extra.loc[ids].to_numpy(np.float32)])
        allowed = sorted(np.unique(y).tolist())
        thresholds = fixed_sr_thresholds(allowed, args.sr_boundary)
        xtr, ytr = x[is_train], y[is_train]
        xte, yte = x[~is_train], y[~is_train]
        counts = {c: int((yte == c).sum()) for c in allowed}
        thin = [c for c, n in counts.items() if n < args.min_class_n]
        log(f"[{fi+1}/{len(fold_files)}] {fpath.name}  features {x.shape[1]:,}  "
            f"train {xtr.shape[0]:,}  val {xte.shape[0]:,}")
        log(f"    val class counts {counts}")
        if thin:
            log(f"    classes below --min-class-n {args.min_class_n}: {thin} "
                f"— kept in bal_OT, excluded from bal_OT_common")

        for ci, cfg in enumerate(configs):
            model = make_model(cfg, args.seed + fi, args.n_jobs)
            model.fit(xtr, ytr)
            # both aggregations off one fit: median is a different read of the
            # same trees, and the per-tree matrix is computed once for both.
            # Only bagged ensembles have trees to read that way -- boosting and
            # ridge produce one prediction, so they are scored once as "mean".
            if cfg.get("family", "rf") in ("rf", "et"):
                per_tree = np.vstack([t.predict(xte) for t in model.estimators_])
                tr_tree = np.vstack([t.predict(xtr) for t in model.estimators_])
                preds = (("mean", per_tree.mean(axis=0), tr_tree.mean(axis=0)),
                         ("median", np.median(per_tree, axis=0),
                          np.median(tr_tree, axis=0)))
            else:
                preds = (("mean", model.predict(xte), model.predict(xtr)),)
            for agg, raw, raw_tr in preds:
                pred = snap(raw, allowed, thresholds)
                m = metrics_for(yte, pred, raw, args.sr_boundary,
                                args.min_class_n,
                                cut=tuned_cut(raw_tr, ytr, args.sr_boundary))
                m.update(config=cfg["tag"], aggregation=agg,
                         fold=fpath.stem, n_features=int(x.shape[1]))
                fold_rows.append(m)
                if pred_rows is not None:
                    collect_predictions(pred_rows, sample_ids[~is_train], yte,
                                        raw, pred, cfg["tag"], agg, fpath.stem)
                err = np.abs(yte - pred)
                for c in allowed:
                    mask = yte == c
                    if mask.any():
                        class_rows.append({"config": cfg["tag"], "aggregation": agg,
                                           "fold": fpath.stem, "log2_mic": c,
                                           "n": int(mask.sum()),
                                           "exact": float(np.mean(err[mask] == 0)),
                                           "OT": float(np.mean(err[mask] <= 1))})
            if (ci + 1) % 10 == 0:
                log(f"    {ci+1}/{len(configs)} configs, {time.time()-t0:.0f}s")
        log(f"    fold done in {time.time()-t0:.0f}s")

    write_predictions(pred_rows, args.save_predictions)
    folds_df = pd.DataFrame(fold_rows)
    folds_df.to_csv(out / "folds.tsv", sep="\t", index=False)

    metric_cols = [c for c in folds_df.columns
                   if c not in ("config", "aggregation", "fold")]
    summary = (folds_df.groupby(["config", "aggregation"], as_index=False)[metric_cols]
               .mean())
    summary["config"] = summary["config"] + "_" + summary["aggregation"]
    summary = summary.sort_values("bal_OT", ascending=False)
    summary.to_csv(out / "summary.tsv", sep="\t", index=False)

    cls = pd.DataFrame(class_rows)
    (cls.groupby(["config", "aggregation", "log2_mic"], as_index=False)
        [["n", "exact", "OT"]].mean()).to_csv(out / "per_class.tsv",
                                              sep="\t", index=False)

    log(f"\n{args.label}: {len(summary)} scored configs")
    show = ["config", "bal_OT", "bal_OT_common", "OT", "MAE", "vMER", "MER",
            "SR_f1"]
    print(summary[show].head(10).to_string(index=False))
    print(f"\nbal_OT   best {summary['bal_OT'].max():.3f}   "
          f"median top10 {summary['bal_OT'].nlargest(10).median():.3f}   "
          f"median all {summary['bal_OT'].median():.3f}")
    by_agg = summary.groupby("aggregation")["bal_OT"].median()
    print("median bal_OT by aggregation: "
          + "  ".join(f"{k} {v:.3f}" for k, v in by_agg.items()))
    fam_of = {c["tag"]: c.get("family", "rf") for c in configs}
    summary["family"] = summary["config"].str.rsplit("_", n=1).str[0].map(fam_of)
    if summary["family"].notna().any():
        print("\nby family (median / best bal_OT):")
        for fam, g in summary.groupby("family"):
            print(f"  {fam:<6}{len(g):>4} scored   median {g['bal_OT'].median():.3f}"
                  f"   best {g['bal_OT'].max():.3f}")
    log(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
