#!/usr/bin/env python3
"""Build target clumps for every gene + upstream flank in an MTB reference GFF.

Emits the schema consumed by 13l_build_snp_gene_consensus_chain_sequences.py, so
the consensus-chain extraction used for P. aeruginosa applies unchanged to a
whole-genome MTB target set.

Window convention matches extract_tbprofiler_reference_targets.py:
  gene            (gene_start, gene_end)
  upstream_Nbp    strand '-' -> (gene_end+1,   gene_end+N)
                  strand '+' -> (gene_start-N, gene_start-1)
clipped to the contig, and target_id formatted as "<gene>:gene" /
"<gene>:upstream_300bp" to match the existing sample_target_seqid_map.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path
from typing import Dict, List, Tuple

CLUMP_COLUMNS = [
    "target_id",
    "target_type",
    "reference_contig",
    "target_start",
    "target_end",
    "reference_strand",
    "reference_gene_id",
    "reference_gene_symbol",
    "reference_product_name",
    "fallback_anchor_pos",
    "anchors_all",
    "anchors_genic",
    "anchors_upstream200",
]


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def read_contig_lengths(fasta: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    name = None
    n = 0
    with open_maybe_gz(fasta) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def parse_attributes(raw: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in raw.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
        elif " " in item:  # GTF style
            k, v = item.split(" ", 1)
            v = v.strip('"')
        else:
            continue
        out[k.strip().lower()] = v.strip()
    return out


def pick(attrs: Dict[str, str], keys: Tuple[str, ...], default: str = "") -> str:
    for k in keys:
        v = attrs.get(k, "")
        if v:
            return re.sub(r"^(gene|cds|id):", "", v, flags=re.I)
    return default


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--reference-fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--upstream-bp", type=int, default=300)
    ap.add_argument(
        "--feature-types",
        default="gene",
        help="comma-separated GFF column-3 types to keep (e.g. gene or CDS,rRNA,tRNA)",
    )
    ap.add_argument("--min-length", type=int, default=1, help="skip genes shorter than this")
    args = ap.parse_args()

    keep = {t.strip().lower() for t in args.feature_types.split(",") if t.strip()}
    lengths = read_contig_lengths(Path(args.reference_fasta))
    if not lengths:
        raise SystemExit(f"[ERROR] no contigs parsed from {args.reference_fasta}")
    print(f"[INFO] reference contigs: {len(lengths)}  total_bp={sum(lengths.values()):,}")

    genes: List[dict] = []
    seen_ids: Dict[str, int] = {}
    type_counts: Dict[str, int] = {}
    skipped_contig = 0

    with open_maybe_gz(Path(args.gff)) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2].strip().lower()
            type_counts[ftype] = type_counts.get(ftype, 0) + 1
            if ftype not in keep:
                continue
            contig = f[0].strip()
            if contig not in lengths:
                skipped_contig += 1
                continue
            start, end = int(f[3]), int(f[4])
            if end - start + 1 < args.min_length:
                continue
            attrs = parse_attributes(f[8])
            gid = pick(attrs, ("locus_tag", "id", "gene_id", "name", "gene"))
            if not gid:
                gid = f"{contig}_{start}_{end}"
            if gid in seen_ids:  # keep first, disambiguate duplicates
                seen_ids[gid] += 1
                gid = f"{gid}.{seen_ids[gid]}"
            else:
                seen_ids[gid] = 0
            genes.append(
                {
                    "gene_id": gid,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": f[6].strip() or "+",
                    "symbol": pick(attrs, ("gene", "name", "gene_name")),
                    "product": pick(attrs, ("product", "description", "note")),
                }
            )

    print(f"[INFO] GFF feature types seen: {dict(sorted(type_counts.items(), key=lambda x: -x[1])[:8])}")
    print(f"[INFO] genes kept ({'/'.join(sorted(keep))}): {len(genes):,}")
    if skipped_contig:
        print(f"[WARN] {skipped_contig} features on contigs absent from the reference FASTA")
    if not genes:
        raise SystemExit("[ERROR] no features matched --feature-types; check the GFF")

    rows: List[dict] = []
    n_clipped = 0
    n_no_upstream = 0
    up = args.upstream_bp
    for g in genes:
        clen = lengths[g["contig"]]
        base = {
            "reference_contig": g["contig"],
            "reference_strand": g["strand"],
            "reference_gene_id": g["gene_id"],
            "reference_gene_symbol": g["symbol"],
            "reference_product_name": g["product"],
            "anchors_all": "",
            "anchors_genic": "",
            "anchors_upstream200": "",
        }
        gs, ge = max(1, g["start"]), min(clen, g["end"])
        rows.append(
            {
                **base,
                "target_id": f"{g['gene_id']}:gene",
                "target_type": "gene",
                "target_start": gs,
                "target_end": ge,
                "fallback_anchor_pos": gs,
            }
        )
        if g["strand"] == "-":
            us, ue = g["end"] + 1, g["end"] + up
        else:
            us, ue = g["start"] - up, g["start"] - 1
        if us < 1 or ue > clen:
            n_clipped += 1
        us, ue = max(1, us), min(clen, ue)
        if us <= ue:
            rows.append(
                {
                    **base,
                    "target_id": f"{g['gene_id']}:upstream_{up}bp",
                    "target_type": f"upstream_{up}bp",
                    "target_start": us,
                    "target_end": ue,
                    "fallback_anchor_pos": us,
                }
            )
        else:
            n_no_upstream += 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=CLUMP_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    n_gene = sum(1 for r in rows if r["target_type"] == "gene")
    n_up = len(rows) - n_gene
    span = sum(r["target_end"] - r["target_start"] + 1 for r in rows)
    print(f"[OK] targets written : {len(rows):,}  ({n_gene:,} gene + {n_up:,} upstream_{up}bp)")
    if n_clipped:
        print(f"[INFO] upstream windows clipped at a contig edge: {n_clipped}")
    if n_no_upstream:
        print(f"[INFO] genes with no room for an upstream window: {n_no_upstream}")
    print(f"[OK] total target span : {span:,} bp ({100 * span / sum(lengths.values()):.1f}% of reference, overlaps counted)")
    print(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
