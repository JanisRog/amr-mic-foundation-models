#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def read_fasta(path: Path):
    name = None
    seq = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def split_id(original_id: str, start: int, end: int, sequence: str) -> str:
    digest = hashlib.sha256(f"{original_id}:{start}:{end}:{sequence}".encode()).hexdigest()[:20]
    return f"{original_id}__part_{start}_{end}__{digest}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--map-tsv", required=True, type=Path)
    parser.add_argument("--threshold", type=int, default=8192)
    parser.add_argument("--parts", type=int, default=2)
    args = parser.parse_args()

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.map_tsv.parent.mkdir(parents=True, exist_ok=True)
    n_input = n_output = n_split = 0

    with args.output_fasta.open("w") as fasta_out, args.map_tsv.open("w", newline="") as map_out:
        writer = csv.DictWriter(
            map_out,
            delimiter="\t",
            fieldnames=[
                "split_sequence_id", "original_sequence_id", "part_index", "n_parts",
                "start_1based", "end_1based", "original_length", "split_length",
            ],
            lineterminator="\n",
        )
        writer.writeheader()

        for original_id, sequence in read_fasta(args.input_fasta):
            n_input += 1
            if len(sequence) <= args.threshold:
                fasta_out.write(f">{original_id}\n{sequence}\n")
                writer.writerow({
                    "split_sequence_id": original_id,
                    "original_sequence_id": original_id,
                    "part_index": 1,
                    "n_parts": 1,
                    "start_1based": 1,
                    "end_1based": len(sequence),
                    "original_length": len(sequence),
                    "split_length": len(sequence),
                })
                n_output += 1
                continue

            n_split += 1
            part_size = (len(sequence) + args.parts - 1) // args.parts
            for part_index in range(args.parts):
                start0 = part_index * part_size
                end0 = min(len(sequence), start0 + part_size)
                if start0 >= end0:
                    continue
                subseq = sequence[start0:end0]
                sid = split_id(original_id, start0 + 1, end0, subseq)
                fasta_out.write(f">{sid}\n{subseq}\n")
                writer.writerow({
                    "split_sequence_id": sid,
                    "original_sequence_id": original_id,
                    "part_index": part_index + 1,
                    "n_parts": args.parts,
                    "start_1based": start0 + 1,
                    "end_1based": end0,
                    "original_length": len(sequence),
                    "split_length": len(subseq),
                })
                n_output += 1

    print({
        "input_records": n_input,
        "output_records": n_output,
        "split_original_records": n_split,
        "output_fasta": str(args.output_fasta),
        "map_tsv": str(args.map_tsv),
    })


if __name__ == "__main__":
    main()
