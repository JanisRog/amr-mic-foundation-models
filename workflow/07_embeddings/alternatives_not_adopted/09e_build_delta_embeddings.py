#!/usr/bin/env python
"""Build delta embeddings: each sequence minus its target's anchor allele.

Whole-gene mean-pooled embeddings are dominated by gene identity, so a single
causal SNP (one 6-mer token of ~217) is buried. Subtracting a per-target anchor
cancels the shared gene-identity component and the anisotropic common direction,
leaving what the model thinks the variant did.

Every seq_id belongs to exactly one target, so the output keeps the same row
order and index as the input and is a drop-in replacement for embeddings.npy.

Anchors:
  modal     most frequent allele among TRAINING samples (leak-free, default)
  reference explicit target -> seq_id table (e.g. H37Rv wild-type windows)
  mean      per-target mean over training occurrences
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import shutil
from pathlib import Path
from typing import Dict, List

import numpy as np


def read_index(path: Path) -> Dict[str, int]:
    seq2row: Dict[str, int] = {}
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        ci = header.index("seq_id")
        cr = header.index("row")
        for line in handle:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(ci, cr):
                seq2row[p[ci]] = int(p[cr])
    return seq2row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb-dir", required=True, help="dir with embeddings.npy + embedding_index.tsv")
    ap.add_argument("--map-path", default="", help="sample_target_seqid_map.tsv (default: inside --emb-dir)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--split-file", default="", help="restrict anchor estimation to split==train")
    ap.add_argument("--anchor", default="modal", choices=["modal", "reference", "mean"])
    ap.add_argument("--reference-map", default="", help="TSV with target_id and seq_id, for --anchor reference")
    ap.add_argument("--batch-size", type=int, default=8192)
    args = ap.parse_args()

    emb_dir = Path(args.emb_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    idx_path = emb_dir / "embedding_index.tsv"
    map_path = Path(args.map_path) if args.map_path else emb_dir / "sample_target_seqid_map.tsv"

    seq2row = read_index(idx_path)
    X = np.load(emb_dir / "embeddings.npy", mmap_mode="r")
    n, d = X.shape
    print(f"[INFO] embeddings {X.shape}  index {len(seq2row)}")

    train_ids = None
    if args.split_file:
        train_ids = set()
        with open(args.split_file) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if str(row.get("split", "")).strip().lower() == "train":
                    train_ids.add(str(row["sample_id"]))
        print(f"[INFO] training samples for anchor estimation: {len(train_ids):,}")

    seq_target: Dict[str, str] = {}
    counts: Dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    with open(map_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            s, t = row["seq_id"], row["target_id"]
            seq_target[s] = t
            if train_ids is None or str(row["sample_id"]) in train_ids:
                counts[t][s] += 1
    targets = sorted(counts)
    print(f"[INFO] targets {len(targets)}  sequences mapped {len(seq_target):,}")

    anchor_vec: Dict[str, np.ndarray] = {}
    anchor_seq: Dict[str, str] = {}

    if args.anchor == "reference":
        if not args.reference_map:
            raise SystemExit("[ERROR] --anchor reference requires --reference-map")
        with open(args.reference_map) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                anchor_seq[row["target_id"]] = row["seq_id"]
        missing = [t for t in targets if t not in anchor_seq or anchor_seq[t] not in seq2row]
        if missing:
            raise SystemExit(f"[ERROR] {len(missing)} targets lack a usable reference seq_id, e.g. {missing[:3]}")
        for t in targets:
            anchor_vec[t] = np.asarray(X[seq2row[anchor_seq[t]]], dtype=np.float64)
    elif args.anchor == "modal":
        for t in targets:
            s = counts[t].most_common(1)[0][0]
            anchor_seq[t] = s
            anchor_vec[t] = np.asarray(X[seq2row[s]], dtype=np.float64)
    else:  # mean
        for t in targets:
            tot = 0
            acc = np.zeros(d, dtype=np.float64)
            for s, c in counts[t].items():
                r = seq2row.get(s)
                if r is None:
                    continue
                acc += np.asarray(X[r], dtype=np.float64) * c
                tot += c
            anchor_seq[t] = f"<mean of {len(counts[t])} alleles>"
            anchor_vec[t] = acc / max(tot, 1)

    rows_target: List[str] = [""] * n
    for s, r in seq2row.items():
        t = seq_target.get(s)
        if t is not None:
            rows_target[r] = t
    unmapped = sum(1 for t in rows_target if not t)
    if unmapped:
        print(f"[WARN] {unmapped} embedding rows have no target in the map; left unchanged")

    A = np.zeros((n, d), dtype=np.float32)
    for t in targets:
        A[[r for r, tt in enumerate(rows_target) if tt == t]] = anchor_vec[t].astype(np.float32)

    out_npy = out_dir / "embeddings.npy"
    D = np.lib.format.open_memmap(out_npy, mode="w+", dtype=np.float32, shape=(n, d))
    zero_rows = 0
    for start in range(0, n, args.batch_size):
        end = min(start + args.batch_size, n)
        block = np.asarray(X[start:end], dtype=np.float32) - A[start:end]
        D[start:end] = block
        zero_rows += int((np.abs(block).max(axis=1) == 0).sum())
    D.flush()

    shutil.copyfile(idx_path, out_dir / "embedding_index.tsv")
    if map_path.exists() and not (out_dir / "sample_target_seqid_map.tsv").exists():
        try:
            (out_dir / "sample_target_seqid_map.tsv").symlink_to(map_path.resolve())
        except OSError:
            pass

    with (out_dir / "anchor_report.tsv").open("w") as handle:
        handle.write("target_id\tanchor_seq_id\tn_alleles_train\tanchor_share\n")
        for t in targets:
            tot = sum(counts[t].values())
            top = counts[t].most_common(1)[0][1]
            handle.write(f"{t}\t{anchor_seq[t]}\t{len(counts[t])}\t{top / max(tot, 1):.4f}\n")

    info = {
        "anchor": args.anchor,
        "n_sequences": int(n),
        "embedding_dim": int(d),
        "n_targets": len(targets),
        "n_rows_exactly_zero": int(zero_rows),
        "split_file": args.split_file,
        "source_emb_dir": str(emb_dir),
    }
    (out_dir / "delta_manifest.json").write_text(json.dumps(info, indent=2, sort_keys=True))
    print(json.dumps(info, indent=2, sort_keys=True))
    print(f"\n[OK] wrote {out_npy}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
