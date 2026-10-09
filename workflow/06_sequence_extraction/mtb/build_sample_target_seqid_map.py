#!/usr/bin/env python3
"""Join 13l's per-sample sequences to their sequence ids, for the ML root.

`13l_build_snp_gene_consensus_chain_sequences.py` writes `sample_target_sequences
.tsv` with the sequence spelled out, and separately a table of unique sequences
with their ids. The precompute reads neither: it wants `sample_target_seqid_map
.tsv`, one row per sample x target carrying the *id* of the sequence, so the
embedding matrix can be indexed without ever materialising a sequence again.

An equivalent map exists for the all-gene block, built on 2026-08-13 by something
that left a manifest but no script. This reproduces it column for column and,
more importantly, reproduces its checks -- `rows_missing_seq_id` and
`seq_ids_missing_from_embedding_index` both had to be 0 there, and they are the
only thing standing between a silently truncated block and a modelling result
that looks like a finding.

Streams the sequence table, which for the all-gene case is 75 GB. Only the unique
sequences are held in memory, and they are the small side by construction: 532,255
against 84M rows there, 4,264 against 346k here.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import pathlib
import sys

OUT_COLUMNS = ["sample_id", "target_id", "seq_id", "target_type",
               "projection_status", "consensus_status", "sequence_length"]


def load_unique(tsv: pathlib.Path | None,
                fasta: pathlib.Path | None) -> dict[str, str]:
    """sequence -> seq_id, from whichever of the two is given."""
    out: dict[str, str] = {}
    if tsv and tsv.exists():
        with tsv.open() as fh:
            r = csv.DictReader(fh, delimiter="\t")
            if "seq_id" not in r.fieldnames or "sequence" not in r.fieldnames:
                raise SystemExit(f"[ERROR] {tsv} needs seq_id and sequence "
                                 f"columns; has {r.fieldnames}")
            for row in r:
                out[row["sequence"]] = row["seq_id"]
        return out
    if not (fasta and fasta.exists()):
        raise SystemExit("[ERROR] give --unique-tsv or --unique-fasta")
    sid, parts = None, []
    with fasta.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if sid is not None:
                    out["".join(parts)] = sid
                sid, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if sid is not None:
        out["".join(parts)] = sid
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-target-sequences", required=True)
    ap.add_argument("--unique-tsv", default="")
    ap.add_argument("--unique-fasta", default="")
    ap.add_argument("--embedding-index", default="",
                    help="embedding_index.tsv from one model's run; every seq_id "
                         "in the map must appear in it, or the block is short and "
                         "some isolates would silently lose a target")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    src = pathlib.Path(args.sample_target_sequences)
    uniq = load_unique(pathlib.Path(args.unique_tsv) if args.unique_tsv else None,
                       pathlib.Path(args.unique_fasta) if args.unique_fasta
                       else None)
    print(f"[INFO] {len(uniq):,} unique sequences loaded")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = missing = 0
    seen: set[str] = set()
    with src.open() as fh, out.open("w", newline="") as oh:
        r = csv.DictReader(fh, delimiter="\t")
        need = {"sample_id", "target_id", "sequence"}
        absent = need - set(r.fieldnames or ())
        if absent:
            raise SystemExit(f"[ERROR] {src} is missing {sorted(absent)}")
        w = csv.DictWriter(oh, delimiter="\t", lineterminator="\n",
                           fieldnames=OUT_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for row in r:
            rows += 1
            sid = uniq.get(row["sequence"])
            if sid is None:
                missing += 1
                continue
            seen.add(sid)
            w.writerow({
                "sample_id": row["sample_id"],
                "target_id": row["target_id"],
                "seq_id": sid,
                "target_type": row.get("target_type", ""),
                "projection_status": row.get("projection_status", ""),
                "consensus_status": row.get("consensus_status", ""),
                "sequence_length": row.get("sequence_length", ""),
            })
            if rows % 5_000_000 == 0:
                print(f"[INFO] {rows:,} rows")

    idx_missing = 0
    if args.embedding_index:
        p = pathlib.Path(args.embedding_index)
        with p.open() as fh:
            r = csv.DictReader(fh, delimiter="\t")
            col = next((c for c in (r.fieldnames or [])
                        if c in ("seq_id", "sequence_id", "id")), None)
            if col is None:
                raise SystemExit(f"[ERROR] no seq_id column in {p}; "
                                 f"has {r.fieldnames}")
            have = {row[col] for row in r}
        idx_missing = len(seen - have)

    man = {
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "sample_target_sequences": str(src),
        "unique_window_sequences_tsv": args.unique_tsv,
        "unique_window_sequences_fasta": args.unique_fasta,
        "embedding_index": args.embedding_index,
        "out": str(out),
        "counts": {
            "rows": rows,
            "unique_sequences_loaded": len(uniq),
            "rows_missing_seq_id": missing,
            "unique_seq_ids_in_output": len(seen),
            "seq_ids_missing_from_embedding_index": idx_missing,
        },
    }
    pathlib.Path(str(out) + ".manifest.json").write_text(
        json.dumps(man, indent=2) + "\n")
    print(f"[OK] {out}\n     rows={rows:,} unique_in_output={len(seen):,} "
          f"missing_seq_id={missing:,} missing_from_index={idx_missing:,}")
    if missing or idx_missing:
        raise SystemExit("[ERROR] the map and the embeddings do not agree; "
                         "both counts must be 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
