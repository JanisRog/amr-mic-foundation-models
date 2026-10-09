#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path


TABLES = (
    "sample_target_sequences.tsv",
    "sample_target_seqid_map.tsv",
    "drug_resistance_binary.tsv",
    "mutation_features_long.tsv",
)


def read_sample_ids(path: Path) -> list[str]:
    sample_ids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(sample_ids) != len(set(sample_ids)):
        raise RuntimeError(f"duplicate sample IDs in {path}")
    return sample_ids


def append_table(source: Path, destination, expected_header: str | None) -> tuple[str, int]:
    with source.open() as handle:
        header = handle.readline()
        if not header:
            raise RuntimeError(f"empty table: {source}")
        if expected_header is not None and header != expected_header:
            raise RuntimeError(f"header mismatch: {source}")
        if expected_header is None:
            destination.write(header)
        rows = 0
        for line in handle:
            destination.write(line)
            rows += 1
    return header, rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-root", required=True, type=Path)
    parser.add_argument("--sample-list", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--expected-rows-per-sample", type=int, default=145)
    args = parser.parse_args()

    sample_ids = read_sample_ids(args.sample_list)
    args.outdir.mkdir(parents=True, exist_ok=True)
    tmpdir = args.outdir / f".gather_tmp_{os.getpid()}"
    tmpdir.mkdir()

    outputs = {name: (tmpdir / name).open("w") for name in TABLES}
    headers: dict[str, str | None] = {name: None for name in TABLES}
    row_counts = {name: 0 for name in TABLES}
    unique_ids: set[str] = set()
    fasta_path = tmpdir / "unique_window_sequences.fasta"
    completed = 0

    try:
        with fasta_path.open("w") as fasta:
            for index, sample_id in enumerate(sample_ids, 1):
                sample_dir = args.samples_root / sample_id
                for name in TABLES:
                    source = sample_dir / name
                    if not source.is_file() or source.stat().st_size == 0:
                        raise RuntimeError(f"missing/empty sample output: {source}")
                    header, rows = append_table(source, outputs[name], headers[name])
                    if headers[name] is None:
                        headers[name] = header
                    row_counts[name] += rows

                with (sample_dir / "sample_target_sequences.tsv").open() as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    for row in reader:
                        seq_id = row["unique_sequence_id"]
                        if seq_id not in unique_ids:
                            unique_ids.add(seq_id)
                            fasta.write(f">{seq_id}\n{row['sequence']}\n")

                completed += 1
                if index % 250 == 0 or index == len(sample_ids):
                    print(
                        f"[INFO] gathered_samples={index}/{len(sample_ids)} "
                        f"unique_sequences={len(unique_ids)}",
                        flush=True,
                    )
    finally:
        for handle in outputs.values():
            handle.close()

    expected_sequence_rows = len(sample_ids) * args.expected_rows_per_sample
    for name in ("sample_target_sequences.tsv", "sample_target_seqid_map.tsv"):
        if row_counts[name] != expected_sequence_rows:
            raise RuntimeError(
                f"{name} rows={row_counts[name]} expected={expected_sequence_rows}"
            )

    for name in TABLES:
        os.replace(tmpdir / name, args.outdir / name)
    os.replace(fasta_path, args.outdir / fasta_path.name)
    tmpdir.rmdir()

    manifest = {
        "samples_expected": len(sample_ids),
        "samples_gathered": completed,
        "expected_rows_per_sample": args.expected_rows_per_sample,
        "unique_sequences": len(unique_ids),
        "row_counts": row_counts,
    }
    (args.outdir / "gather_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    print(f"[OK] gathered={args.outdir}")


if __name__ == "__main__":
    main()
