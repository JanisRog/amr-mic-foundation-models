#!/usr/bin/env python3
"""Fit each fold's projection and selection once, cache it, reuse it everywhere.

Every architecture in the grid has to see the same folds and the same PCA fits, or
differences between them are partly differences between their splits. And every
fit has to happen inside the fold's training rows only: the previous packages
built one feature matrix on all training data and then cross-validated on it,
which lets the selection see the labels of the samples it is later scored on.
That does not merely inflate the numbers, it biases the comparison along the axis
being compared -- mRMR 64 gains more from having seen the labels than mRMR 256
does, so it wins for the wrong reason.

So: for each of the 15 partitions, fit the projection on that partition's training
occurrences, rank targets by mRMR on those same training rows, and cache the
result. An architecture run then loads a fold and slices it.

Two economies make this affordable:

  * only the top --top-targets targets are cached. mRMR 64 / 128 / 256 are
    prefixes of one ranking, so three architectures share one cache.
  * the per-target-variance cache also carries a per-PC importance mask, so the
    model-selected-PC mode is the same cache with the mask applied rather than a
    third precompute.

Per-target PCA is fitted on the target's *unique* sequences weighted by how often
they occur in training, which is the same estimator as duplicating rows but over
~65 sequences instead of ~2,800.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import sys
import time

import numpy as np
import pandas as pd


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def weighted_pca(x: np.ndarray, w: np.ndarray, max_components: int,
                 var_threshold: float) -> tuple[np.ndarray, np.ndarray, int]:
    """Occurrence-weighted PCA of unique rows. Returns (mean, components, n_kept).

    Weighting by occurrence count is what makes this identical to running PCA on
    the duplicated sample-level matrix; without it a sequence carried by one
    isolate would pull as hard as the wild type carried by two thousand.

    Two regimes. Per-target, the matrix is a target's ~65 unique alleles by 4096
    dimensions, and a full economy SVD is instant. Shared across all targets it is
    ~282,000 x 4096, where a full SVD would solve for all 4,096 singular vectors
    and allocate a 4.6 GB U it never uses, to keep four of them -- so that case
    goes through a randomized SVD for exactly the components requested. Total
    variance for the threshold comes from the Frobenius norm, which is available
    either way.
    """
    total = float(w.sum())
    mean = (w[:, None] * x).sum(axis=0) / total
    xc = (x - mean) * np.sqrt(w)[:, None]
    total_var = float((xc.astype(np.float64) ** 2).sum())
    if total_var <= 0:
        return mean, np.zeros((0, x.shape[1]), dtype=np.float32), 0

    k = min(max_components, *xc.shape)
    if min(xc.shape) > 512 and k <= min(xc.shape) // 4:
        from sklearn.utils.extmath import randomized_svd
        _, s, vt = randomized_svd(xc, n_components=k, n_iter=7, random_state=0)
    else:
        try:
            _, s, vt = np.linalg.svd(xc, full_matrices=False)
        except np.linalg.LinAlgError:
            return mean, np.zeros((0, x.shape[1]), dtype=np.float32), 0
    var = s.astype(np.float64) ** 2
    ratio = np.cumsum(var) / total_var
    n_keep = int(np.searchsorted(ratio, var_threshold) + 1)
    n_keep = max(1, min(n_keep, max_components, vt.shape[0]))
    return mean, vt[:n_keep].astype(np.float32), n_keep


# same column conventions the ML-root validator uses, so one file format is
# understood in one way across the project
ORIG_NAMES = ("original_sequence_id", "orig_seq_id", "source_seq_id", "seq_id")
PART_NAMES = ("split_sequence_id", "part_seq_id", "chunk_seq_id", "split_id")


def pick_col(header: list[str], names, exclude=()) -> str | None:
    for n in names:
        if n in header and n not in exclude:
            return n
    return None


def load_split_maps(paths, seq2row: dict[str, int]) -> dict[str, list[int]]:
    """Sequences too long to embed whole were cut in half before embedding.

    The index holds the halves, the sample-target map holds the originals, so an
    original resolves to the mean of its parts -- the same convention
    build_multimodel_embedding_block.py uses. Without this every long gene reads
    as a missing embedding, which is exactly how the first real run failed.
    """
    out: dict[str, list[int]] = {}
    for path in paths:
        with open(path) as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            header = rdr.fieldnames or []
            part_col = pick_col(header, PART_NAMES)
            orig_col = pick_col(header, ORIG_NAMES, exclude=(part_col,))
            if not part_col or not orig_col:
                raise SystemExit(f"[ERROR] {path}: cannot identify the original "
                                 f"and part columns in {header}")
            groups: dict[str, list[str]] = {}
            for r in rdr:
                groups.setdefault(str(r[orig_col]), []).append(str(r[part_col]))
        usable = incomplete = 0
        for orig, parts in groups.items():
            if not orig or orig in seq2row:
                continue
            if all(x in seq2row for x in parts):
                out[orig] = [seq2row[x] for x in parts]
                usable += 1
            else:
                incomplete += 1
        print(f"[INFO] split map {pathlib.Path(path).name}: {len(groups):,} split "
              f"sequences ({orig_col} -> {part_col}), {usable:,} resolvable"
              + (f", {incomplete:,} with missing parts" if incomplete else ""))
    return out


def project(x: np.ndarray, mean: np.ndarray, comp: np.ndarray,
            batch: int = 200_000) -> np.ndarray:
    """Centre and project, in batches, with the result checked rather than trusted.

    Batched because `x - mean` on the full gathered embedding matrix would
    allocate a second copy of it -- 8 GB for evo2 at this cohort size, which is
    what turned a comfortable job into an OOM.

    numpy 2 also emits a spurious "divide by zero encountered in matmul" on some
    of these small BLAS calls even when every input and output is finite.
    Silencing it blind would silence a real one too, so the warning is suppressed
    and the output asserted finite instead.
    """
    if comp.shape[0] == 0:
        return np.zeros((x.shape[0], 0), dtype=np.float32)
    out = np.empty((x.shape[0], comp.shape[0]), dtype=np.float32)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        for i in range(0, x.shape[0], batch):
            out[i:i + batch] = (x[i:i + batch] - mean) @ comp.T
    if not np.isfinite(out).all():
        raise SystemExit("[ERROR] non-finite values in a projected block")
    return out


def standardise(m: np.ndarray) -> np.ndarray:
    """Column-wise z-score; constant columns become zeros so they never rank."""
    mu = m.mean(axis=0)
    sd = m.std(axis=0)
    out = np.zeros_like(m, dtype=np.float32)
    ok = sd > 1e-12
    out[:, ok] = ((m[:, ok] - mu[ok]) / sd[ok]).astype(np.float32)
    return out


def best_pc(block_train: np.ndarray, y_train: np.ndarray) -> tuple[int, float]:
    """The PC of one target most correlated with the label, and that |corr|.

    Done through standardise() rather than np.corrcoef because a target that is
    monomorphic in this fold's training isolates has constant PC columns, and
    corrcoef divides by their zero standard deviation -- a NaN that is correct in
    meaning but arrives as a warning on every fold. Here a constant column
    becomes zeros, so its correlation is 0 and mRMR drops it, which is the same
    outcome stated rather than repaired.
    """
    z = standardise(block_train)
    yz = standardise(y_train.reshape(-1, 1))[:, 0]
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        c = np.abs(z.T @ yz) / len(yz)
    if not c.size:
        return 0, 0.0
    return int(np.argmax(c)), float(c.max())


def mrmr_order(rep_train: np.ndarray, y_train: np.ndarray, k: int,
               diag: dict | None = None) -> np.ndarray:
    """Greedy mRMR in difference form over one representative vector per target.

    Relevance is |corr| with the label, redundancy the mean |corr| with what is
    already chosen. Both are dot products once the columns are standardised, so
    selecting 256 of 8,243 costs a handful of matrix-vector products.
    """
    z = standardise(rep_train)
    n = z.shape[0]
    yz = standardise(y_train.reshape(-1, 1))[:, 0]
    # same spurious numpy 2 matmul warning as in project(); the inputs are
    # standardised and finite, and the result is checked below
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        relevance = np.abs(z.T @ yz) / n
    if not np.isfinite(relevance).all():
        raise SystemExit("[ERROR] non-finite relevance in mRMR")
    alive = relevance > 0
    if not alive.any():
        raise SystemExit("[ERROR] no target correlates with the label at all")

    if diag is not None:
        diag["relevance"] = relevance.copy()
    chosen: list[int] = []
    red_sum = np.zeros(z.shape[1], dtype=np.float64)
    score = relevance.copy()
    score[~alive] = -np.inf
    for step in range(min(k, int(alive.sum()))):
        pick = int(np.argmax(score))
        chosen.append(pick)
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            red_sum += np.abs(z.T @ z[:, pick]) / n
        score = relevance - red_sum / (step + 1)
        score[~alive] = -np.inf
        score[chosen] = -np.inf
    if diag is not None:
        diag["redundancy"] = red_sum / max(len(chosen), 1)
    return np.array(chosen, dtype=np.int64)


def pc_importance_mask(x: np.ndarray, y: np.ndarray, seed: int) -> np.ndarray:
    """Which PCs of one target a model actually uses.

    A small forest on that target's own PCs, keeping any PC whose importance
    beats the uniform share. Applied after mRMR has narrowed the target list, so
    it is a few hundred tiny fits per fold rather than 8,243.
    """
    from sklearn.ensemble import ExtraTreesRegressor
    n_pc = x.shape[1]
    if n_pc <= 1:
        return np.ones(n_pc, dtype=bool)
    est = ExtraTreesRegressor(n_estimators=60, max_depth=8, random_state=seed,
                              n_jobs=1)
    est.fit(x, y)
    imp = est.feature_importances_
    mask = imp > (1.0 / n_pc)
    if not mask.any():                      # never return an empty target
        mask[int(np.argmax(imp))] = True
    return mask


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="per-model ML root")
    ap.add_argument("--embedding-subdir", default="ml_inputs/embeddings/tb_profiler")
    ap.add_argument("--folds-tsv", required=True)
    ap.add_argument("--pca-mode", choices=("shared4", "ptvar"), required=True)
    ap.add_argument("--pca-dim", type=int, default=4)
    ap.add_argument("--var-threshold", type=float, default=0.99)
    ap.add_argument("--max-components", type=int, default=25)
    ap.add_argument("--top-targets", type=int, default=256)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--save-ranking", action="store_true",
                    help="write every target's mRMR relevance, mean redundancy and "
                         "final rank. The two columns answer separately whether a "
                         "target was judged uninformative or merely crowded out")
    ap.add_argument("--save-projection", action="store_true",
                    help="write the fitted per-target means and components beside "
                         "the fold cache. Needed when the feature set has to be "
                         "applied to isolates that are not in this run at all -- a "
                         "deployment, or a second cohort -- since the cache alone "
                         "only holds the samples it was built from")
    ap.add_argument("--only-fold", default="",
                    help="'repeat:fold' to build a single partition, for timing")
    ap.add_argument("--limit-targets", type=int, default=0,
                    help="first N targets only, for a smoke test")
    ap.add_argument("--chunk-rows", type=int, default=5_000_000,
                    help="rows per chunk when streaming the sample-target map")
    ap.add_argument("--exclude-loci", default="",
                    help="locus tags whose targets are removed before ranking, e.g. "
                         "Rv0682,Rv3919c,Rv2043c. Separate with commas, spaces or "
                         "semicolons, or pass @path to read one per line -- sbatch "
                         "--export treats commas as its own separator, so a "
                         "comma-joined list silently truncates at the first one. "
                         "Both the :gene and :upstream_300bp target of each locus "
                         "is dropped. Section 1: co-resistance markers predict well "
                         "and do not transport, so a model meant to carry across "
                         "settings must not be allowed to read them")
    ap.add_argument("--fixed-targets", default="",
                    help="path to a target list -- one target_id per line, or a "
                         "TSV whose first column is target_id -- used verbatim in "
                         "every fold instead of running mRMR. The final evaluation "
                         "needs one feature set, not fifteen: selection refitted "
                         "per fold is what makes the CV honest, but it leaves no "
                         "single model to report. Fixing the list moves the "
                         "selection outside the run, so it must come from data "
                         "that never reaches the holdout -- see section 4. PCA is "
                         "still fitted inside each fold's training rows, which "
                         "stays honest because it never sees a label")
    ap.add_argument("--zero-fill-targets", default="",
                    help="targets where a missing sample-target cell is a zero "
                         "vector instead of a fatal error. For an acquired "
                         "resistance gene absence is the normal state, not a "
                         "failure -- but the check this exempts is the one that "
                         "caught 27 isolates with no embeddings for rifampicin, "
                         "so it is exempted per target rather than switched off. "
                         "Every target outside this list must still be present "
                         "for every isolate")
    ap.add_argument("--mrmr-restrict-samples", default="",
                    help="compute the mRMR ranking on the training isolates that "
                         "are also in this cohort file, while the PCA and the "
                         "cached features still cover everyone. For isoniazid the "
                         "cohort is the rifampicin-susceptible set: selecting "
                         "targets there keeps rpoB out of the ranking (section 1), "
                         "and training on everyone keeps the sample size. The "
                         "restriction is applied *after* the fold's training rows, "
                         "so nothing from a validation fold enters the selection")
    ap.add_argument("--split-map", action="append", default=[],
                    help="split-long-sequence map; auto-discovered in the "
                         "embedding directory when not given")
    args = ap.parse_args()

    root = pathlib.Path(args.root)
    emb_dir = root / args.embedding_subdir
    index_path = emb_dir / "embedding_index.tsv"
    map_path = emb_dir / "sample_target_seqid_map.tsv"
    emb_path = emb_dir / "embeddings.npy"
    for p in (index_path, map_path, emb_path):
        if not p.exists():
            raise SystemExit(f"[ERROR] missing {p}")

    folds = pd.read_csv(args.folds_tsv, sep="\t")
    folds["sample_id"] = folds["sample_id"].astype(str)
    samples = sorted(folds["sample_id"].unique())
    sample_pos = {s: i for i, s in enumerate(samples)}
    log(f"{len(samples):,} development samples, "
        f"{folds['repeat'].nunique()} repeats x {folds['fold'].nunique()} folds")

    log("loading embedding index")
    idx = pd.read_csv(index_path, sep="\t", usecols=["row", "seq_id"])
    seq2row = dict(zip(idx["seq_id"].astype(str), idx["row"].astype(int)))
    n_real_rows = max(seq2row.values()) + 1
    log(f"{len(seq2row):,} embedded sequences")

    split_paths = args.split_map or sorted(
        str(p) for p in emb_dir.glob("unique_window_sequences.split_long_*.map.tsv"))
    if not split_paths:
        log("[WARN] no split map found — long sequences will not resolve")
    splits = load_split_maps(split_paths, seq2row)
    virtual_parts: list[list[int]] = []
    for orig, part_rows in splits.items():
        seq2row[orig] = n_real_rows + len(virtual_parts)
        virtual_parts.append(part_rows)
    if virtual_parts:
        log(f"{len(virtual_parts):,} sequences resolve as the mean of their parts")

    # The map is ~84M rows across the full cohort. Reading it whole and then
    # filtering costs tens of gigabytes of Python strings before a single row of
    # it is used, so it is streamed and reduced to integers chunk by chunk.
    log(f"streaming {map_path.name} in chunks")
    target_pos: dict[str, int] = {}
    cells: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    probe: list[tuple[str, str, str]] = []
    kept = seen = 0
    for chunk in pd.read_csv(map_path, sep="\t",
                             usecols=["sample_id", "target_id", "seq_id"],
                             dtype=str, chunksize=args.chunk_rows):
        seen += len(chunk)
        chunk = chunk[chunk["sample_id"].isin(sample_pos)]
        if chunk.empty:
            continue
        for t in chunk["target_id"].unique():
            if t not in target_pos:
                target_pos[t] = len(target_pos)
        rows = chunk["seq_id"].map(seq2row)
        if rows.isna().any():
            missing = chunk.loc[rows.isna(), "seq_id"].iloc[0]
            raise SystemExit(f"[ERROR] sequence {missing!r} has no embedding row")
        cells.append((chunk["sample_id"].map(sample_pos).to_numpy(np.int32),
                      chunk["target_id"].map(target_pos).to_numpy(np.int32),
                      rows.to_numpy(np.int32)))
        # keep a few (sample, target) -> seq_id pairs to verify the assembled
        # matrix against the file afterwards
        if len(probe) < 200:
            take = chunk.head(3)
            for sid, tid, seq in zip(take["sample_id"], take["target_id"],
                                     take["seq_id"]):
                probe.append((sid, tid, seq))
        kept += len(chunk)
        if seen % (20 * args.chunk_rows) == 0:
            log(f"    {seen:,} rows read, {kept:,} kept")
    log(f"{seen:,} map rows read, {kept:,} for this cohort, "
        f"{len(target_pos):,} targets")

    # The codes handed out during streaming are in the order targets were first
    # encountered; `targets` is sorted. Those two orders are not the same, so the
    # codes MUST be remapped -- doing it only under --limit-targets, as this once
    # did, silently attaches every target's name to a different target's column.
    insertion = target_pos
    targets = sorted(insertion)
    if args.exclude_loci:
        spec = args.exclude_loci
        if spec.startswith("@"):
            spec = pathlib.Path(spec[1:]).read_text()
        drop = {x for x in re.split(r"[,;\s]+", spec) if x}
        before = len(targets)
        targets = [t for t in targets if t.split(":")[0] not in drop]
        gone = before - len(targets)
        log(f"excluding {gone} targets from {len(drop)} loci: {sorted(drop)}")
        if not gone:
            raise SystemExit(f"[ERROR] --exclude-loci matched nothing; locus tags "
                             f"are matched against the part before ':' , e.g. "
                             f"Rv0682 for Rv0682:gene")
    if args.limit_targets:
        targets = targets[:args.limit_targets]
    target_pos = {t: i for i, t in enumerate(targets)}

    # A fixed list is resolved here, against the same `targets` the codes are
    # remapped onto, so a name that survived --exclude-loci is the only kind that
    # can be asked for. Order is the caller's: it is a ranking, and truncating it
    # later must mean the same thing as asking for fewer targets now.
    fixed_idx: np.ndarray | None = None
    if args.fixed_targets:
        raw = pathlib.Path(args.fixed_targets.lstrip("@")).read_text().splitlines()
        want = []
        for line in raw:
            tok = line.split("\t")[0].strip()
            if not tok or tok.startswith("#") or tok == "target_id":
                continue
            if tok not in want:
                want.append(tok)
        if not want:
            raise SystemExit(f"[ERROR] --fixed-targets {args.fixed_targets} is empty")
        missing = [t for t in want if t not in target_pos]
        if missing:
            raise SystemExit(
                f"[ERROR] {len(missing)} of {len(want)} fixed targets are not in "
                f"this run's target pool, e.g. {missing[:5]}. A target excluded by "
                f"--exclude-loci cannot also be requested by --fixed-targets")
        fixed_idx = np.array([target_pos[t] for t in want], dtype=np.int64)
        log(f"fixed target list: {len(fixed_idx)} targets from "
            f"{args.fixed_targets}; mRMR will not run")
    new_of = np.full(len(insertion), -1, dtype=np.int32)
    for t, i in target_pos.items():
        new_of[insertion[t]] = i

    rowmat = np.full((len(samples), len(targets)), -1, dtype=np.int32)
    for s_code, t_code, row in cells:
        mapped = new_of[t_code]
        keep = mapped >= 0
        rowmat[s_code[keep], mapped[keep]] = row[keep]
    del cells

    # verify the assembled matrix against pairs read straight from the file
    bad = 0
    for sid, tid, seq in probe:
        if tid not in target_pos:
            continue
        if rowmat[sample_pos[sid], target_pos[tid]] != seq2row[seq]:
            bad += 1
    if bad:
        raise SystemExit(f"[ERROR] {bad} of {len(probe)} verified cells disagree "
                         f"with the map — target columns are misaligned")
    log(f"verified {len(probe)} sample-target cells against the map")
    zero_ok = np.zeros(len(targets), dtype=bool)
    if args.zero_fill_targets:
        want_z = {l.split("\t")[0].strip()
                  for l in open(args.zero_fill_targets.lstrip("@"))
                  if l.strip() and not l.startswith("#")}
        want_z.discard("target_id")
        for tname, i in target_pos.items():
            if tname in want_z:
                zero_ok[i] = True
        log(f"{int(zero_ok.sum())} of {len(targets)} targets may be absent "
            f"(zero-filled)")
    missing = rowmat < 0
    hard = missing & ~zero_ok[None, :]
    if hard.any():
        n = int(hard.sum())
        cols = np.flatnonzero(hard.any(axis=0))[:3]
        raise SystemExit(f"[ERROR] {n:,} sample-target cells have no sequence; "
                         f"this builder expects exactly one per pair. "
                         f"e.g. {[targets[c] for c in cols]}. Add them to "
                         f"--zero-fill-targets only if absence is meaningful "
                         f"there, not to silence a coverage gap")
    n_zero = int(missing.sum())
    if n_zero:
        log(f"{n_zero:,} cells zero-filled "
            f"({n_zero / missing.size:.1%} of the matrix)")

    # one gather of every row this cohort touches, then all 15 folds are compute
    used = np.unique(rowmat)
    store = np.load(emb_path, mmap_mode="r")
    gb = used.size * store.shape[1] * 4 / 1e9
    log(f"gathering {used.size:,} unique sequence embeddings ({gb:.1f} GB, "
        f"dim {store.shape[1]})")
    emb = np.empty((used.size, store.shape[1]), dtype=np.float32)
    real = used < n_real_rows
    step = 50_000
    real_idx = np.flatnonzero(real)
    for i in range(0, real_idx.size, step):
        sel = real_idx[i:i + step]
        emb[sel] = store[used[sel]]
    for j in np.flatnonzero(~real):
        parts = virtual_parts[int(used[j]) - n_real_rows]
        emb[j] = np.asarray(store[parts], dtype=np.float32).mean(axis=0)
    if (~real).any():
        log(f"    {int((~real).sum()):,} of them averaged from split parts")
    remap = np.full(int(used.max()) + 1, -1, dtype=np.int32)
    remap[used] = np.arange(used.size, dtype=np.int32)
    # -1 would index the end of `remap` under numpy's negative indexing and give a
    # silently wrong row, so missing cells are clamped and then redirected to an
    # explicit zero row appended to `emb`.
    local = remap[np.maximum(rowmat, 0)]       # samples x targets -> row in emb
    if n_zero:
        emb = np.vstack([emb, np.zeros((1, emb.shape[1]), dtype=np.float32)])
        local[missing] = emb.shape[0] - 1
    del store, rowmat

    mrmr_cohort = None
    if args.mrmr_restrict_samples:
        mrmr_cohort = {l.strip() for l in open(args.mrmr_restrict_samples)
                       if l.strip()}
        n = sum(1 for s in samples if s in mrmr_cohort)
        log(f"mRMR cohort: {n:,} of {len(samples):,} development isolates")
        if not n:
            raise SystemExit("[ERROR] no development isolate is in the mRMR cohort")

    y_all = (folds.drop_duplicates("sample_id")
             .set_index("sample_id")["log2_mic"].reindex(samples).to_numpy(float))

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "targets.txt").write_text("\n".join(targets) + "\n")

    pairs = sorted({(int(r), int(f)) for r, f in
                    zip(folds["repeat"], folds["fold"])})
    if args.only_fold:
        r, _, f = args.only_fold.partition(":")
        pairs = [(int(r), int(f))]

    for rep, fold in pairs:
        t0 = time.time()
        val_ids = set(folds.loc[(folds["repeat"] == rep) & (folds["fold"] == fold),
                                "sample_id"])
        is_train = np.array([s not in val_ids for s in samples])
        y_train = y_all[is_train]
        log(f"--- repeat {rep} fold {fold}: {is_train.sum():,} train, "
            f"{(~is_train).sum():,} val")

        # ---- projection, fitted on training rows only --------------------
        rep_vec = np.zeros((len(samples), len(targets)), dtype=np.float32)
        shared_mean = shared_comp = None
        if args.pca_mode == "shared4":
            train_rows = local[is_train].ravel()
            counts = np.bincount(train_rows, minlength=emb.shape[0]).astype(np.float64)
            nz = counts > 0
            shared_mean, shared_comp, _ = weighted_pca(
                emb[nz], counts[nz], args.pca_dim, 1.0)
            proj_all = project(emb, shared_mean, shared_comp)    # unique rows x dim
            log(f"    shared PCA{shared_comp.shape[0]} on "
                f"{int(counts.sum()):,} training occurrences")
            flat = 0
            for ti in range(len(targets)):
                block = proj_all[local[:, ti]]
                pc, corr = best_pc(block[is_train], y_train)
                flat += corr == 0.0
                rep_vec[:, ti] = block[:, pc]
            log(f"    {flat:,} of {len(targets):,} targets are monomorphic in this "
                f"fold's training isolates and cannot be ranked")
        else:
            flat = 0
            n_pc_all = np.zeros(len(targets), dtype=np.int32)
            for ti in range(len(targets)):
                col = local[:, ti]
                uniq, inv = np.unique(col[is_train], return_inverse=True)
                w = np.bincount(inv).astype(np.float64)
                mean, comp, n_keep = weighted_pca(emb[uniq], w,
                                                  args.max_components,
                                                  args.var_threshold)
                n_pc_all[ti] = n_keep
                if n_keep == 0:
                    flat += 1
                    continue
                block = project(emb[col], mean, comp)
                pc, corr = best_pc(block[is_train], y_train)
                flat += corr == 0.0
                rep_vec[:, ti] = block[:, pc]
            log(f"    per-target PCA: median {int(np.median(n_pc_all))} PCs, "
                f"max {int(n_pc_all.max())}")
            log(f"    {flat:,} of {len(targets):,} targets are monomorphic in this "
                f"fold's training isolates and cannot be ranked")

        # ---- target selection: a fixed list, or mRMR on training rows ----
        order = fixed_idx
        if order is not None:
            log(f"    using the fixed list of {len(order)} targets")
        sel_mask = is_train
        if order is None and mrmr_cohort is not None:
            sel_mask = np.array([t and s in mrmr_cohort
                                 for t, s in zip(is_train, samples)])
            log(f"    mRMR selection restricted to {int(sel_mask.sum()):,} of "
                f"{int(is_train.sum()):,} training isolates in the cohort")
            if sel_mask.sum() < 50:
                raise SystemExit("[ERROR] fewer than 50 training isolates in the "
                                 "mRMR cohort — the ranking would be noise")
        diag = {} if (args.save_ranking and order is None) else None
        if order is None:
            order = mrmr_order(rep_vec[sel_mask], y_all[sel_mask],
                               args.top_targets, diag)
            log(f"    mRMR ranked {len(order)} targets in {time.time()-t0:.0f}s")
        if diag is not None:
            rank_of = {int(t): i + 1 for i, t in enumerate(order)}
            rel = diag["relevance"]
            by_rel = np.argsort(-rel)
            rel_rank = {int(t): i + 1 for i, t in enumerate(by_rel)}
            with (out / f"ranking_r{rep}_f{fold}.tsv").open("w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["target_id", "relevance", "relevance_rank",
                            "mean_redundancy", "mrmr_rank"])
                for ti in by_rel:
                    w.writerow([targets[ti], f"{rel[ti]:.6f}", rel_rank[int(ti)],
                                f"{diag['redundancy'][ti]:.6f}",
                                rank_of.get(int(ti), "")])
            log(f"    wrote ranking_r{rep}_f{fold}.tsv "
                f"(relevance for all {len(rel):,} targets)")

        # ---- features for the selected targets ---------------------------
        blocks, names, n_pcs, masks = [], [], [], []
        zero_width: list[str] = []
        for rank, ti in enumerate(order):
            col = local[:, ti]
            if args.pca_mode == "shared4":
                block = project(emb[col], shared_mean, shared_comp)
            else:
                uniq, inv = np.unique(col[is_train], return_inverse=True)
                w = np.bincount(inv).astype(np.float64)
                mean, comp, n_keep = weighted_pca(emb[uniq], w,
                                                  args.max_components,
                                                  args.var_threshold)
                block = project(emb[col], mean, comp)
            block = block.astype(np.float32)
            # Only reachable with --fixed-targets: mRMR requires relevance > 0 so
            # it never picks a target that is monomorphic in the training rows,
            # but a fixed list may name one in some folds and not others. It
            # contributes no columns; dropping it silently would make the feature
            # count differ per fold without saying so, hence the count below.
            if block.shape[1] == 0:
                zero_width.append(targets[ti])
            blocks.append(block)
            n_pcs.append(block.shape[1])
            names.extend(f"{targets[ti]}_pc{j}" for j in range(block.shape[1]))
            masks.append(pc_importance_mask(block[is_train], y_train,
                                            args.seed + rank)
                         if args.pca_mode == "ptvar" and block.shape[1]
                         else np.ones(block.shape[1], dtype=bool))
        if zero_width:
            log(f"    {len(zero_width)} fixed targets are monomorphic in this "
                f"fold's training isolates and contribute no columns: "
                f"{zero_width[:5]}")

        x = np.hstack(blocks)
        if args.save_projection:
            proj = {"target_ids": np.array([targets[t] for t in order], dtype=object),
                    "pca_mode": args.pca_mode}
            if args.pca_mode == "shared4":
                proj["shared_mean"] = shared_mean
                proj["shared_components"] = shared_comp
            else:
                for rank, ti in enumerate(order):
                    col = local[:, ti]
                    uniq, inv = np.unique(col[is_train], return_inverse=True)
                    w_ = np.bincount(inv).astype(np.float64)
                    mean, comp, _ = weighted_pca(emb[uniq], w_,
                                                 args.max_components,
                                                 args.var_threshold)
                    proj[f"mean_{rank}"] = mean
                    proj[f"comp_{rank}"] = comp
            np.savez_compressed(out / f"projection_r{rep}_f{fold}.npz", **proj)
            log(f"    wrote projection_r{rep}_f{fold}.npz "
                f"({len(order)} targets, mode {args.pca_mode})")
        np.savez_compressed(
            out / f"fold_r{rep}_f{fold}.npz",
            x=x, y=y_all, is_train=is_train,
            samples=np.array(samples, dtype=object),
            target_ids=np.array([targets[t] for t in order], dtype=object),
            n_pcs=np.array(n_pcs, dtype=np.int32),
            pc_mask=np.concatenate(masks),
            feature_names=np.array(names, dtype=object))
        log(f"    wrote fold_r{rep}_f{fold}.npz  x={x.shape}  "
            f"kept PCs {int(np.concatenate(masks).sum())}/{x.shape[1]}  "
            f"({time.time()-t0:.0f}s)")

    (out / "precompute_manifest.json").write_text(json.dumps({
        "root": str(root), "pca_mode": args.pca_mode, "pca_dim": args.pca_dim,
        "var_threshold": args.var_threshold, "max_components": args.max_components,
        "top_targets": args.top_targets, "n_samples": len(samples),
        "n_targets_total": len(targets), "folds_tsv": str(args.folds_tsv),
        "fixed_targets": args.fixed_targets or None,
        "n_fixed_targets": None if fixed_idx is None else int(len(fixed_idx)),
        "seed": args.seed}, indent=2) + "\n")
    log(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
