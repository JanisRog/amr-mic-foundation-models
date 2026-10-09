#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} MANIFEST_DIR OUTPUT_LIST")

    manifest_dir = Path(sys.argv[1])
    output = Path(sys.argv[2])
    paired: list[str] = []
    invalid: list[str] = []

    for manifest in sorted(manifest_dir.glob("*.fastq_manifest.tsv")):
        sample = manifest.name.removesuffix(".fastq_manifest.tsv")
        indexes: dict[str, int] = {"1": 0, "2": 0}
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                index = row.get("read_index", "").strip()
                fastq = Path(row.get("fastq", "").strip())
                if index in indexes and fastq.is_file() and fastq.stat().st_size > 0:
                    indexes[index] += 1
        if indexes["1"] > 0 and indexes["1"] == indexes["2"]:
            paired.append(sample)
        else:
            invalid.append(sample)

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(f"{sample}\n" for sample in paired))
    invalid_path = output.with_suffix(".invalid.txt")
    invalid_path.write_text("".join(f"{sample}\n" for sample in invalid))

    print(f"[OK] paired={len(paired)} invalid={len(invalid)}")
    print(f"[OK] paired_list={output}")
    print(f"[OK] invalid_list={invalid_path}")


if __name__ == "__main__":
    main()
