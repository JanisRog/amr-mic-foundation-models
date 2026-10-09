#!/usr/bin/env python3
"""Fuse two embedding blocks into the single layout the precompute reads.

`cv_precompute_fold_features.py` reads one `embeddings.npy`, one
`embedding_index.tsv` and one `sample_target_seqid_map.tsv` from a single
directory. The chromosomal targets and the acquired-gene targets were extracted
and embedded separately -- different sources (PAO1 windows versus AMRFinder hits
on each isolate's own assembly) and different runs -- so they arrive as two such
triples. This concatenates them rather than re-embedding 15,391 sequences to get
one index.

Row numbers in the second block are shifted by the size of the first, which is
the whole trick. Everything else is bookkeeping with two checks worth having:

  **A seq_id appearing in both blocks must carry the same embedding.** The two
  runs used the same model and the same settings, so an id collision means the
  identical sequence was embedded twice and both rows should agree. They are
  compared rather than assumed, because a mismatch would mean the blocks came
  from different models or different hooks, and the merged file would silently
  mix two embedding spaces.

  **Target ids must not collide.** Chromosomal targets are `<locus_tag>:gene`
  and acquired ones are `amrseq_<family>`, so a clash would mean one of the
  extractors changed its naming and the two blocks are no longer distinguishable
  downstream.

The map is deliberately left sparse for the acquired targets: an isolate without
the gene has no row, and the precompute's `--zero-fill-targets` decides what that
means. Writing zeros here instead would bake a modelling choice into a file
format.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

import numpy as np


def read_index(p: pathlib.Path) -> list[str]:
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    if not rows or "seq_id" not in rows[0]:
        raise SystemExit(f"[ERROR] {p} has no seq_id column")
    out = [""] * len(rows)
    for r in rows:
        out[int(r["row"])] = r["seq_id"]
    if any(s == "" for s in out):
        raise SystemExit(f"[ERROR] {p} has gaps in its row numbering")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", action="append", required=True, metavar="DIR",
                    help="directory with embeddings.npy, embedding_index.tsv "
                         "and sample_target_seqid_map.tsv; repeatable, order kept")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    blocks = [pathlib.Path(b) for b in args.block]
    if len(blocks) < 2:
        raise SystemExit("[ERROR] give at least two --block directories")

    arrays, index, offset = [], [], 0
    seen: dict[str, tuple[int, int]] = {}
    rows_out: list[tuple[str, str, str]] = []
    targets_of: dict[str, set[str]] = {}

    for b in blocks:
        emb_p, idx_p = b / "embeddings.npy", b / "embedding_index.tsv"
        map_p = b / "sample_target_seqid_map.tsv"
        for f in (emb_p, idx_p, map_p):
            if not f.exists():
                raise SystemExit(f"[ERROR] missing {f}")
        a = np.load(emb_p, mmap_mode="r")
        ids = read_index(idx_p)
        if len(ids) != a.shape[0]:
            raise SystemExit(f"[ERROR] {b.name}: index has {len(ids):,} rows, "
                             f"array has {a.shape[0]:,}")
        if arrays and a.shape[1] != arrays[0].shape[1]:
            raise SystemExit(f"[ERROR] {b.name} has dim {a.shape[1]}, the first "
                             f"block has {arrays[0].shape[1]} — different models "
                             f"or different hooks cannot be concatenated")
        cur = np.asarray(a, dtype=np.float32)
        dup = 0
        for i, s in enumerate(ids):
            if s in seen:
                bi, bj = seen[s]
                if not np.allclose(arrays[bi][bj], cur[i], atol=1e-5):
                    raise SystemExit(
                        f"[ERROR] seq_id {s} appears in two blocks with "
                        f"different embeddings — the blocks are not from the "
                        f"same model/hook")
                dup += 1
            else:
                seen[s] = (len(arrays), i)
        arrays.append(cur)
        index.extend(ids)

        tset = set()
        with map_p.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows_out.append((r["sample_id"], r["target_id"], r["seq_id"]))
                tset.add(r["target_id"])
        for other, ots in targets_of.items():
            clash = tset & ots
            if clash:
                raise SystemExit(f"[ERROR] target ids appear in both {other} and "
                                 f"{b.name}: {sorted(clash)[:3]}")
        targets_of[b.name] = tset
        print(f"[INFO] {b.name}: {a.shape[0]:,} sequences (dim {a.shape[1]}), "
              f"{len(tset):,} targets, {dup:,} seq_ids shared with an earlier block")
        offset += a.shape[0]

    emb = np.vstack(arrays)
    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "embeddings.npy", emb)
    with (out / "embedding_index.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["row", "seq_id"])
        for i, s in enumerate(index):
            w.writerow([i, s])
    with (out / "sample_target_seqid_map.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "target_id", "seq_id"])
        w.writerows(rows_out)

    n_t = sum(len(v) for v in targets_of.values())
    n_s = len({s for s, _, _ in rows_out})
    print(f"\n[OK] {out}")
    print(f"     {emb.shape[0]:,} sequences (dim {emb.shape[1]}), {n_t:,} targets, "
          f"{len(rows_out):,} cells over {n_s:,} isolates")
    print(f"     density {len(rows_out) / (n_s * n_t):.1%} — the acquired targets "
          f"are sparse by design;\n     pass their ids to --zero-fill-targets.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
