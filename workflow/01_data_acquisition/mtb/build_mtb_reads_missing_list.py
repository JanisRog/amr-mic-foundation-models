#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


def safe_name(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip()).strip("._")
    return text or "sample"


def sample_complete(outroot: Path, sample_id: str) -> bool:
    sample_safe = safe_name(sample_id)
    status_path = outroot / "status" / f"{sample_safe}.tsv"
    manifest_path = outroot / "manifests" / f"{sample_safe}.fastq_manifest.tsv"
    if not status_path.is_file() or not manifest_path.is_file():
        return False

    with status_path.open(newline="") as handle:
        status = next(csv.DictReader(handle, delimiter="\t"), {})
    if status.get("status", "").strip() != "complete":
        return False

    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return bool(rows) and all(
        (path := Path(row.get("fastq", ""))).is_file() and path.stat().st_size > 0
        for row in rows
    )


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(f"usage: {sys.argv[0]} METADATA_TSV OUTROOT TODO_TSV")

    metadata_path = Path(sys.argv[1])
    outroot = Path(sys.argv[2])
    todo_path = Path(sys.argv[3])

    with metadata_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []

    sample_col = "sample_id" if "sample_id" in fields else "Sample_IDs"
    todo = [row for row in rows if not sample_complete(outroot, row.get(sample_col, ""))]

    todo_path.parent.mkdir(parents=True, exist_ok=True)
    with todo_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(todo)

    print(f"[OK] expected={len(rows)} complete={len(rows) - len(todo)} todo={len(todo)}")
    print(f"[OK] todo_tsv={todo_path}")


if __name__ == "__main__":
    main()
