#!/usr/bin/env python3
"""Concatenate several models' embeddings for the same sequences into one block.

The sweep has two block slots (tb_profiler, gwas), so a 3-model combination
cannot use one slot per model. Instead this merges along the *feature* axis:
row i becomes [model1_i | model2_i | ...], keyed on the original (pre-split)
seq_id, so the existing sample_target_seqid_map resolves it unchanged and no
split map is needed downstream.

Two things this has to get right:

1. **Split representation differs per model.** evo2 and DNABERT-2 embedded the
   halved FASTA, Carbon embedded whole sequences. For a sequence that was split,
   the halves are averaged back into one vector so every model contributes a
   row for the same original id.

2. **Scale differs wildly per model** -- median row norms are ~22 (evo2), ~71
   (Carbon), ~4.1 (DNABERT-2). Concatenating raw would let one model dominate the
   shared PCA on magnitude alone, so each block is rescaled by a single constant
   (which cannot change within-model structure).

   The first version matched median row norms, and that was the wrong target.
   The sweep's PCA centers before decomposing, so blocks compete on variance
   AFTER the shared mean is removed -- and evo2's embeddings are collapsed (mean
   pairwise cosine 0.998), meaning most of their norm is a common mean vector
   carrying no information. Norm-matched, evo2 lost every component to Carbon:
   katG relevance fell from 0.809 (evo2 alone) to 0.291 (carbon+evo2), and the
   merged selection shared 209/256 targets with carbon-alone against 108 with
   evo2-alone. Default is now centered_rms.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_index(d: Path) -> dict[str, int]:
    idx = pd.read_csv(d / "embedding_index.tsv", sep="\t", usecols=["row", "seq_id"])
    return dict(zip(idx["seq_id"].astype(str), idx["row"].astype(int)))


def load_split_map(path: Path | None) -> dict[str, list[str]]:
    """original_sequence_id -> [split ids], tolerating the column-name variants."""
    if path is None:
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    orig_names = ("original_sequence_id", "orig_seq_id", "source_seq_id", "seq_id")
    part_names = ("split_sequence_id", "part_seq_id", "chunk_seq_id", "split_id")
    part = next((c for c in part_names if c in df.columns), None)
    orig = next((c for c in orig_names if c in df.columns and c != part), None)
    if part is None or orig is None:
        raise SystemExit(f"[ERROR] cannot identify columns in {path}: {list(df.columns)}")
    out: dict[str, list[str]] = {}
    for o, p in zip(df[orig], df[part]):
        out.setdefault(str(o), []).append(str(p))
    return out


def resolvable_ids(direct: dict[str, int], splits: dict[str, list[str]]) -> set[str]:
    ids = set(direct)
    for orig, parts in splits.items():
        if all(p in direct for p in parts):
            ids.add(orig)
    return ids


def block_scale(emb, direct: dict[str, int], ids: list[str], sample: int, seed: int,
                mode: str) -> float:
    """Divisor that puts this model's block on a comparable footing.

    median_norm equalises magnitude. That is the wrong target: the sweep fits an
    occurrence-weighted PCA, which centers first, so what decides how many
    components a block wins is its variance AFTER the shared mean is removed.
    evo2's embeddings are collapsed (mean pairwise cosine 0.998) -- nearly all of
    their norm is a common mean vector carrying no information -- so norm-matching
    left evo2 with far less centered variance than Carbon and it lost every
    component. centered_rms equalises the quantity that actually competes.
    """
    rng = np.random.default_rng(seed)
    pick = ids if len(ids) <= sample else [ids[i] for i in rng.choice(len(ids), sample, replace=False)]
    rows = sorted({direct[i] for i in pick if i in direct})
    if not rows:
        return 1.0
    x = np.asarray(emb[rows], dtype=np.float32)
    if mode == "median_norm":
        v = float(np.median(np.linalg.norm(x, axis=1)))
    else:  # centered_rms
        x = x - x.mean(axis=0, keepdims=True)
        v = float(np.sqrt((x.astype(np.float64) ** 2).sum(axis=1).mean()))
    return v if v > 0 else 1.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", action="append", required=True, metavar="NAME=DIR",
                    help="repeatable, in the order the dimensions should be concatenated")
    ap.add_argument("--split-map", action="append", default=[], metavar="NAME=PATH",
                    help="repeatable; only for models that embedded a split fasta")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--normalize", choices=["centered_rms", "median_norm", "none"],
                    default="centered_rms",
                    help="centered_rms equalises post-centering variance, which is what the "
                         "downstream PCA competes on; median_norm equalises magnitude only")
    ap.add_argument("--norm-sample", type=int, default=5000)
    ap.add_argument("--chunk", type=int, default=4096)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    models: list[tuple[str, Path]] = []
    for spec in args.model:
        name, _, d = spec.partition("=")
        if not name or not d:
            raise SystemExit(f"[ERROR] --model expects NAME=DIR, got {spec!r}")
        models.append((name, Path(d)))
    smaps = {}
    for spec in args.split_map:
        name, _, p = spec.partition("=")
        smaps[name] = Path(p)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    embs, directs, splits, dims = {}, {}, {}, {}
    shared: set[str] | None = None
    for name, d in models:
        e = np.load(d / "embeddings.npy", mmap_mode="r")
        direct = load_index(d)
        if len(direct) != e.shape[0]:
            raise SystemExit(f"[ERROR] {name}: index rows {len(direct)} != embedding rows {e.shape[0]}")
        sp = load_split_map(smaps.get(name))
        ids = resolvable_ids(direct, sp)
        embs[name], directs[name], splits[name], dims[name] = e, direct, sp, int(e.shape[1])
        print(f"[INFO] {name}: rows={e.shape[0]:,} dim={e.shape[1]} "
              f"split_originals={len(sp):,} resolvable={len(ids):,}")
        shared = ids if shared is None else (shared & ids)

    assert shared is not None
    ordered = sorted(shared)
    total_dim = sum(dims[n] for n, _ in models)
    print(f"[INFO] shared sequences={len(ordered):,}  merged dim={total_dim}")
    if not ordered:
        raise SystemExit("[ERROR] no sequences shared across all models")

    scale = {}
    for name, _ in models:
        if args.normalize == "none":
            scale[name] = 1.0
        else:
            m = block_scale(embs[name], directs[name], ordered, args.norm_sample,
                            args.seed, args.normalize)
            scale[name] = 1.0 / m
            print(f"[INFO] {name}: {args.normalize} {m:.4f} -> scale {scale[name]:.6f}")

    merged = np.lib.format.open_memmap(
        out / "embeddings.npy", mode="w+", dtype=np.float32, shape=(len(ordered), total_dim))

    n_averaged = {n: 0 for n, _ in models}
    for start in range(0, len(ordered), args.chunk):
        block = ordered[start:start + args.chunk]
        col = 0
        for name, _ in models:
            e, direct, sp, dim = embs[name], directs[name], splits[name], dims[name]
            buf = np.empty((len(block), dim), dtype=np.float32)
            for k, sid in enumerate(block):
                r = direct.get(sid)
                if r is not None:
                    buf[k] = e[r]
                else:
                    rows = [direct[p] for p in sp[sid]]
                    buf[k] = np.asarray(e[sorted(rows)], dtype=np.float32).mean(axis=0)
                    n_averaged[name] += 1
            merged[start:start + len(block), col:col + dim] = buf * scale[name]
            col += dim
        if start % (args.chunk * 20) == 0:
            print(f"  {start + len(block):,}/{len(ordered):,}")
    merged.flush()
    del merged

    pd.DataFrame({"row": np.arange(len(ordered)), "seq_id": ordered}).to_csv(
        out / "embedding_index.tsv", sep="\t", index=False)

    manifest = {
        "models": [{"name": n, "dir": str(d), "dim": dims[n], "scale": scale[n],
                    "halves_averaged": n_averaged[n]} for n, d in models],
        "normalize": args.normalize,
        "n_sequences": len(ordered),
        "merged_dim": total_dim,
        "note": "rows keyed on original (pre-split) seq_id; no split map needed downstream",
    }
    (out / "merge_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"[OK] {out}/embeddings.npy  shape=({len(ordered)}, {total_dim})")
    for n, _ in models:
        print(f"[OK] {n}: {n_averaged[n]:,} sequences reconstructed by averaging halves")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
