#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read_samples(path: Path) -> list[str]:
    samples = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(samples) != len(set(samples)):
        raise RuntimeError(f"duplicate sample IDs in {path}")
    return samples


def normalized_drugs(value: str) -> tuple[str, ...]:
    return tuple(sorted({drug.strip() for drug in value.split(";") if drug.strip() and drug.strip() != "-"}))


def safe_token(value: str) -> str:
    return value.strip().replace("\t", " ").replace("\n", " ")


def write_matrix(
    path: Path,
    samples: list[str],
    features: list[str],
    sample_features: dict[str, set[str]],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", *features])
        for sample in samples:
            present = sample_features.get(sample, set())
            writer.writerow([sample, *(int(feature in present) for feature in features)])


def write_long(path: Path, rows: set[tuple[str, str, str, str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "feature", "feature_class", "gene_name", "change", "drugs"])
        for row in sorted(rows):
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants-csv", required=True, type=Path)
    parser.add_argument("--sample-list", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    samples = read_samples(args.sample_list)
    sample_set = set(samples)
    args.outdir.mkdir(parents=True, exist_ok=True)

    mutation_by_sample: dict[str, set[str]] = defaultdict(set)
    gene_by_sample: dict[str, set[str]] = defaultdict(set)
    combined_by_sample: dict[str, set[str]] = defaultdict(set)
    long_rows: set[tuple[str, str, str, str, str, str]] = set()
    unknown_samples: set[str] = set()
    input_rows = 0
    drug_annotated_rows = 0

    with args.variants_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"sample", "gene_name", "change", "drugs"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise RuntimeError(f"missing columns: {sorted(missing)}")

        for row in reader:
            input_rows += 1
            sample = safe_token(row["sample"])
            gene = safe_token(row["gene_name"])
            change = safe_token(row["change"])
            drugs = normalized_drugs(row["drugs"])
            if sample not in sample_set:
                unknown_samples.add(sample)
                continue
            if not gene:
                continue

            if drugs:
                feature = f"gene::{gene}"
                gene_by_sample[sample].add(feature)
                combined_by_sample[sample].add(feature)
                long_rows.add((sample, feature, "amr_gene_hit", gene, "", ";".join(drugs)))

            if drugs and change:
                drug_annotated_rows += 1
                feature = f"mutation::{gene}::{change}"
                mutation_by_sample[sample].add(feature)
                combined_by_sample[sample].add(feature)
                long_rows.add((sample, feature, "resistance_mutation", gene, change, ";".join(drugs)))

    mutation_features = sorted({feature for values in mutation_by_sample.values() for feature in values})
    gene_features = sorted({feature for values in gene_by_sample.values() for feature in values})
    combined_features = sorted({feature for values in combined_by_sample.values() for feature in values})

    write_matrix(args.outdir / "tbprofiler_resistance_mutation_binary.tsv", samples, mutation_features, mutation_by_sample)
    write_matrix(args.outdir / "tbprofiler_amr_gene_hit_binary.tsv", samples, gene_features, gene_by_sample)
    write_matrix(args.outdir / "tbprofiler_amrfinder_style_binary.tsv", samples, combined_features, combined_by_sample)
    write_long(args.outdir / "tbprofiler_amrfinder_style_long.tsv", long_rows)

    manifest = {
        "samples": len(samples),
        "input_rows": input_rows,
        "drug_annotated_input_rows": drug_annotated_rows,
        "resistance_mutation_features": len(mutation_features),
        "amr_gene_hit_features": len(gene_features),
        "combined_features": len(combined_features),
        "samples_with_resistance_mutation": sum(bool(mutation_by_sample.get(sample)) for sample in samples),
        "samples_with_amr_gene_hit": sum(bool(gene_by_sample.get(sample)) for sample in samples),
        "unknown_samples": sorted(unknown_samples),
    }
    (args.outdir / "matrix_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    print(f"[OK] outdir={args.outdir}")


if __name__ == "__main__":
    main()
