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
    targets = []
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


def choose_alt(fmt: str, sample: str, alts: list[str]) -> str | None:
    if not alts:
        return None
    if not fmt or not sample:
        return alts[0]
    keys = fmt.split(":")
    values = sample.split(":")
    data = dict(zip(keys, values))
    gt = data.get("GT", "")
    alleles = [int(x) for x in gt.replace("|", "/").split("/") if x.isdigit() and int(x) > 0]
    if alleles and alleles[0] <= len(alts):
        return alts[alleles[0] - 1]
    return alts[0]


def read_variants(path: Path) -> list[dict[str, object]]:
    variants = []
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            alts = fields[4].split(",") if fields[4] not in {"", "."} else []
            alt = choose_alt(fields[8] if len(fields) > 8 else "", fields[9] if len(fields) > 9 else "", alts)
            if not alt or alt.startswith("<") or "[" in alt or "]" in alt or "*" in alt:
                continue
            variants.append({
                "chrom": fields[0],
                "pos": int(fields[1]),
                "ref": fields[3].upper(),
                "alt": alt.upper(),
                "filter": fields[6] if len(fields) > 6 else "",
            })
    return variants


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def apply_variants(
    reference: str,
    chrom: str,
    start: int,
    end: int,
    variants: list[dict[str, object]],
) -> tuple[str, int, int]:
    sequence = reference[start - 1:end]
    selected = [
        v for v in variants
        if v["chrom"] == chrom and start <= int(v["pos"]) and int(v["pos"]) + len(str(v["ref"])) - 1 <= end
    ]
    selected.sort(key=lambda v: int(v["pos"]))
    output = []
    cursor = start
    applied = 0
    conflicts = 0
    for variant in selected:
        pos = int(variant["pos"])
        ref = str(variant["ref"])
        alt = str(variant["alt"])
        if pos < cursor:
            conflicts += 1
            continue
        observed = reference[pos - 1:pos - 1 + len(ref)]
        if observed != ref:
            conflicts += 1
            continue
        output.append(reference[cursor - 1:pos - 1])
        output.append(alt)
        cursor = pos + len(ref)
        applied += 1
    output.append(reference[cursor - 1:end])
    return "".join(output), applied, conflicts


def seq_id(sequence: str) -> str:
    return "tbseq_" + hashlib.sha256(sequence.encode()).hexdigest()[:20]


def collect_drug_names(value: object, under_drug_key: bool = False) -> set[str]:
    names: set[str] = set()
    if isinstance(value, dict):
        for key, item in value.items():
            key_lower = str(key).lower()
            if key_lower in {"drug", "drugs"}:
                names.update(collect_drug_names(item, True))
            elif under_drug_key or key_lower in {"name", "drug_name"}:
                names.update(collect_drug_names(item, under_drug_key))
    elif isinstance(value, list):
        for item in value:
            names.update(collect_drug_names(item, under_drug_key))
    elif under_drug_key and isinstance(value, str) and value.strip():
        names.add(value.strip())
    return names


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def find_single(root: Path, pattern: str) -> Path:
    matches = sorted(root.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"expected one {pattern} under {root}, found {len(matches)}")
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--db-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--upstream-bp", type=int, default=300)
    args = parser.parse_args()

    sample_id = args.sample_dir.name
    reference_path = args.db_dir / "tbdb.fasta"
    gff_path = args.db_dir / "tbdb.gff"
    bed_path = args.db_dir / "tbdb.bed"
    vcf_path = find_single(args.sample_dir / "vcf", "*.targets.vcf.gz")
    json_path = find_single(args.sample_dir / "results", "*.results.json")

    references = read_fasta(reference_path)
    strands = read_gene_strands(gff_path)
    targets = read_targets(bed_path, strands)
    variants = read_variants(vcf_path)
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
            sequence, applied, conflicts = apply_variants(ref, chrom, start, end, variants)
            if strand == "-":
                sequence = revcomp(sequence)
            unique_id = seq_id(sequence)
            target_id = f'{target["gene_id"]}:{target_type}'
            row = {
                "sample_id": sample_id,
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
                "n_applied_variants": applied,
                "n_conflicts": conflicts,
                "unique_sequence_id": unique_id,
            }
            rows.append(row)
            manifest_rows.append({k: row[k] for k in ("target_id", "target_type", "gene_id", "gene_name", "drugs", "chrom", "start", "end", "strand")})

    sequence_fields = list(rows[0])
    write_tsv(args.outdir / "sample_target_sequences.tsv", rows, sequence_fields)
    unique_manifest = {str(row["target_id"]): row for row in manifest_rows}
    write_tsv(
        args.outdir / "target_manifest.tsv",
        list(unique_manifest.values()),
        ["target_id", "target_type", "gene_id", "gene_name", "drugs", "chrom", "start", "end", "strand"],
    )
    write_tsv(
        args.outdir / "sample_target_seqid_map.tsv",
        rows,
        ["sample_id", "target_id", "target_type", "gene_id", "gene_name", "drugs", "unique_sequence_id"],
    )
    unique = {}
    for row in rows:
        unique[str(row["unique_sequence_id"])] = str(row["sequence"])
    with (args.outdir / "unique_window_sequences.fasta").open("w") as handle:
        for unique_id, sequence in sorted(unique.items()):
            handle.write(f">{unique_id}\n{sequence}\n")

    variant_rows = []
    for variant in variants:
        feature = f'{variant["chrom"]}:{variant["pos"]}:{variant["ref"]}>{variant["alt"]}'
        variant_rows.append({"sample_id": sample_id, "feature": feature, "value": 1})
    write_tsv(args.outdir / "mutation_features_long.tsv", variant_rows, ["sample_id", "feature", "value"])

    with json_path.open() as handle:
        result_json = json.load(handle)
    all_drugs = {
        drug.strip()
        for target in targets
        for drug in str(target["drugs"]).split(",")
        if drug.strip()
    }
    resistant_drugs = collect_drug_names(result_json.get("dr_variants", []))
    drug_rows = [{
        "sample_id": sample_id,
        "drug": drug,
        "prediction": "R" if drug in resistant_drugs else "S",
        "resistant_binary": int(drug in resistant_drugs),
    } for drug in sorted(all_drugs | resistant_drugs)]
    write_tsv(
        args.outdir / "drug_resistance_binary.tsv",
        drug_rows,
        ["sample_id", "drug", "prediction", "resistant_binary"],
    )
    summary = {
        "sample_id": sample_id,
        "targets": len(targets),
        "sequence_rows": len(rows),
        "unique_sequences": len(unique),
        "vcf_variants": len(variants),
        "drug_binary_rows": len(drug_rows),
        "rows_with_applied_variants": sum(int(row["n_applied_variants"]) > 0 for row in rows),
        "total_conflicts": sum(int(row["n_conflicts"]) for row in rows),
        "tbprofiler_json_top_level_keys": sorted(result_json),
    }
    (args.outdir / "extract_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"[OK] outdir={args.outdir}")


if __name__ == "__main__":
    main()
