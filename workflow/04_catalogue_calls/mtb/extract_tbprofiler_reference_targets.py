#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = ""
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            else:
                records[name].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def parse_attrs(text: str) -> dict[str, str]:
    attrs = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def read_gene_strands(path: Path) -> dict[str, str]:
    strands: dict[str, str] = {}
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_attrs(fields[8])
            gene_id = attrs.get("gene_id") or attrs.get("ID", "").removeprefix("gene:")
            if gene_id:
                strands[gene_id] = fields[6]
    return strands


def read_targets(path: Path, strands: dict[str, str]) -> list[dict[str, object]]:
    targets: list[dict[str, object]] = []
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            gene_id = fields[3]
            targets.append({
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "gene_id": gene_id,
                "gene_name": fields[4],
                "drugs": fields[5] if len(fields) > 5 else "",
                "strand": strands.get(gene_id, "+"),
            })
    return targets


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def seq_id(sequence: str) -> str:
    return "tbref_" + hashlib.sha256(sequence.encode()).hexdigest()[:20]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sample-id", default="H37Rv_tbprofiler_reference")
    parser.add_argument("--upstream-bp", type=int, default=300)
    args = parser.parse_args()

    reference_path = args.db_dir / "tbdb.fasta"
    gff_path = args.db_dir / "tbdb.gff"
    bed_path = args.db_dir / "tbdb.bed"

    references = read_fasta(reference_path)
    strands = read_gene_strands(gff_path)
    targets = read_targets(bed_path, strands)
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    for target in targets:
        chrom = str(target["chrom"])
        if chrom not in references:
            raise RuntimeError(f"target chromosome missing from reference: {chrom}")
        ref = references[chrom]
        gene_start = int(target["start"])
        gene_end = int(target["end"])
        strand = str(target["strand"])
        intervals = [("gene", gene_start, gene_end)]
        if strand == "-":
            up_start, up_end = gene_end + 1, min(len(ref), gene_end + args.upstream_bp)
        else:
            up_start, up_end = max(1, gene_start - args.upstream_bp), gene_start - 1
        if up_start <= up_end:
            intervals.append((f"upstream_{args.upstream_bp}bp", up_start, up_end))

        for target_type, start, end in intervals:
            sequence = ref[start - 1:end]
            if strand == "-":
                sequence = revcomp(sequence)
            unique_id = seq_id(sequence)
            target_id = f'{target["gene_id"]}:{target_type}'
            row = {
                "sample_id": args.sample_id,
                "target_id": target_id,
                "target_type": target_type,
                "gene_id": target["gene_id"],
                "gene_name": target["gene_name"],
                "drugs": target["drugs"],
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "sequence": sequence,
                "sequence_length": len(sequence),
                "n_applied_variants": 0,
                "n_conflicts": 0,
                "unique_sequence_id": unique_id,
            }
            rows.append(row)
            manifest_rows.append({k: row[k] for k in ("target_id", "target_type", "gene_id", "gene_name", "drugs", "chrom", "start", "end", "strand")})

    if not rows:
        raise RuntimeError("No TBProfiler reference targets extracted")

    write_tsv(args.outdir / "sample_target_sequences.tsv", rows, list(rows[0]))
    write_tsv(
        args.outdir / "target_manifest.tsv",
        list({str(row["target_id"]): row for row in manifest_rows}.values()),
        ["target_id", "target_type", "gene_id", "gene_name", "drugs", "chrom", "start", "end", "strand"],
    )
    write_tsv(
        args.outdir / "sample_target_seqid_map.tsv",
        rows,
        ["sample_id", "target_id", "target_type", "gene_id", "gene_name", "drugs", "unique_sequence_id"],
    )

    unique: dict[str, str] = {}
    for row in rows:
        unique[str(row["unique_sequence_id"])] = str(row["sequence"])
    with (args.outdir / "unique_window_sequences.fasta").open("w") as handle:
        for unique_id, sequence in sorted(unique.items()):
            handle.write(f">{unique_id}\n{sequence}\n")

    summary = {
        "sample_id": args.sample_id,
        "targets": len(targets),
        "sequence_rows": len(rows),
        "unique_sequences": len(unique),
        "target_type_counts": {
            "gene": sum(row["target_type"] == "gene" for row in rows),
            f"upstream_{args.upstream_bp}bp": sum(row["target_type"] == f"upstream_{args.upstream_bp}bp" for row in rows),
        },
    }
    (args.outdir / "extract_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"[OK] outdir={args.outdir}")
    print(f"[OK] unique_fasta={args.outdir / 'unique_window_sequences.fasta'}")


if __name__ == "__main__":
    main()
