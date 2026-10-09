#!/usr/bin/env python3
"""Collapse per-sample target sequences into the layout the embedding step wants.

Each isolate wrote its own `<sample>.targets.tsv` and `<sample>.sequences.fasta`.
Embedding runs on *distinct sequences*, not on isolates: at 62 targets across
~2,100 isolates that is 130,000 extractions, but conserved targets collapse hard
-- for MTB the same step reduced 8,243 targets x 10,170 isolates to 537,768
unique sequences, and on a curated set the ratio is far better still. Embedding
the duplicates would be the single most wasteful thing in the pipeline.

Emits the three files every downstream script already expects, byte-compatible
with the MTB layout:

    unique_window_sequences.fasta   one record per distinct sequence
    sample_target_seqid_map.tsv     sample_id, target_id, seq_id
    embedding_index.tsv             row, seq_id  -- written after embedding, but
                                    the row order is fixed here so the two cannot
                                    disagree

Two hard failures, both learned from MTB. An isolate missing any target is fatal:
`cv_precompute_fold_features.py` refuses to build a cache when a sample-target
cell has no sequence, and finding out here is cheaper than finding out from a
cell count that describes the map rather than the drug. And a seq_id appearing
with two different sequences is fatal: that would mean a sha1 collision or a
corrupted file, and silently keeping one would put the wrong sequence under some
isolate's target.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True,
                    help="directory of <sample>.targets.tsv / .sequences.fasta")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--expect-targets", type=int, default=0,
                    help="every isolate must have exactly this many targets")
    ap.add_argument("--min-samples", type=int, default=1)
    args = ap.parse_args()

    src = pathlib.Path(args.in_dir)
    tsvs = sorted(src.glob("*.targets.tsv"))
    if len(tsvs) < args.min_samples:
        raise SystemExit(f"[ERROR] only {len(tsvs)} sample files in {src}")
    print(f"[INFO] {len(tsvs):,} isolates")

    seqs: dict[str, str] = {}
    rows: list[tuple[str, str, str]] = []
    per_sample: dict[str, int] = {}
    lengths: dict[str, list[int]] = collections.defaultdict(list)

    for i, tp in enumerate(tsvs, 1):
        sample = tp.name[: -len(".targets.tsv")]
        fp = src / f"{sample}.sequences.fasta"
        if not fp.exists():
            raise SystemExit(f"[ERROR] {sample}: no sequences.fasta beside its targets.tsv")
        name, buf = "", []
        local: dict[str, str] = {}
        with fp.open() as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name:
                        local[name] = "".join(buf)
                    name, buf = line[1:].split()[0], []
                else:
                    buf.append(line)
        if name:
            local[name] = "".join(buf)
        for sid, s in local.items():
            prev = seqs.get(sid)
            if prev is not None and prev != s:
                raise SystemExit(
                    f"[ERROR] {sid} maps to two different sequences "
                    f"(lengths {len(prev)} and {len(s)}), first seen before "
                    f"{sample}. A sha1 collision or a corrupted file; either way "
                    f"some isolate would get the wrong sequence.")
            seqs[sid] = s
        n = 0
        with tp.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows.append((r["sample_id"], r["target_id"], r["seq_id"]))
                lengths[r["target_id"]].append(int(r.get("length", 0) or 0))
                n += 1
        per_sample[sample] = n
        if i % 500 == 0:
            print(f"[INFO] {i:,}/{len(tsvs):,} read")

    counts = collections.Counter(per_sample.values())
    expect = args.expect_targets or counts.most_common(1)[0][0]
    short = {s: n for s, n in per_sample.items() if n != expect}
    if short:
        raise SystemExit(
            f"[ERROR] {len(short):,} isolate(s) do not have exactly {expect} "
            f"targets, e.g. " +
            ", ".join(f"{s}={n}" for s, n in list(short.items())[:5]) +
            "\n         The precompute refuses a cache with a missing "
            "sample-target cell.\n         Re-run those array tasks rather than "
            "dropping the isolates silently.")

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    order = sorted(seqs)
    with (out / "unique_window_sequences.fasta").open("w") as fh:
        for sid in order:
            fh.write(f">{sid}\n{seqs[sid]}\n")
    with (out / "sample_target_seqid_map.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "target_id", "seq_id"])
        w.writerows(rows)
    with (out / "embedding_row_order.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["row", "seq_id"])
        for i, sid in enumerate(order):
            w.writerow([i, sid])

    total = len(rows)
    print(f"\n[INFO] {total:,} sample-target cells -> {len(seqs):,} distinct "
          f"sequences ({len(seqs)/total:.1%} of the naive count)")
    print(f"{'target':<28}{'distinct':>9}{'isolates':>10}{'bp':>8}")
    by_t = collections.defaultdict(set)
    for s, t, sid in rows:
        by_t[t].add(sid)
    for t in sorted(by_t, key=lambda t: -len(by_t[t]))[:12]:
        L = lengths[t]
        print(f"{t:<28}{len(by_t[t]):>9}{len(tsvs):>10}"
              f"{(sum(L)//len(L)) if L else 0:>8}")
    print(f"\n[OK] {out}")
    print("Embed unique_window_sequences.fasta in the row order of "
          "embedding_row_order.tsv,\nthen write embeddings.npy and "
          "embedding_index.tsv beside the map.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
