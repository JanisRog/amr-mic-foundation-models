#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd


@dataclass(frozen=True)
class TargetClump:
    target_id: str
    target_type: str
    reference_contig: str
    target_start: int
    target_end: int
    reference_strand: str
    reference_gene_id: str
    reference_gene_symbol: str
    reference_product_name: str
    fallback_anchor_pos: Optional[int]
    anchors_all: Tuple[str, ...]
    anchors_genic: Tuple[str, ...]
    anchors_upstream200: Tuple[str, ...]

    @property
    def target_length(self) -> int:
        return self.target_end - self.target_start + 1


@dataclass
class ChainBlock:
    t_start: int
    size: int
    q_start_raw: int


@dataclass
class ChainRecord:
    t_name: str
    t_size: int
    t_strand: str
    t_start: int
    t_end: int
    q_name: str
    q_size: int
    q_strand: str
    q_start: int
    q_end: int
    chain_id: str
    blocks: List[ChainBlock]


def open_text_auto(path: Path):
    return gzip.open(path, "rt") if str(path).lower().endswith(".gz") else path.open("r")


def read_ids(path: Path) -> List[str]:
    out: List[str] = []
    seen = set()
    with path.open() as fh:
        for line in fh:
            sid = line.strip().split("\t")[0]
            if sid and sid not in seen:
                out.append(sid)
                seen.add(sid)
    return out


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    buf: List[str] = []
    with open_text_auto(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2))


def write_fasta(path: Path, rows: Sequence[Tuple[str, str]]) -> None:
    with path.open("w") as fh:
        for seq_id, seq in rows:
            fh.write(f">{seq_id}\n")
            for i in range(0, len(seq), 120):
                fh.write(seq[i : i + 120] + "\n")


def normalize_snp_id(x: str) -> str:
    s = str(x).strip()
    if s.startswith("snp|"):
        parts = s.split("|", 1)
        if len(parts) == 2:
            return parts[1]
    return s


def parse_anchor_list(raw: str) -> Tuple[str, ...]:
    text = str(raw or "").strip()
    if not text:
        return tuple()
    if text.startswith("["):
        try:
            parsed = json.loads(text)
            vals = [normalize_snp_id(x) for x in parsed if str(x).strip()] if isinstance(parsed, list) else [normalize_snp_id(text)]
        except Exception:
            vals = [normalize_snp_id(x) for x in text.split(";") if str(x).strip()]
    else:
        vals = [normalize_snp_id(x) for x in text.split(";") if str(x).strip()]
    out: List[str] = []
    seen = set()
    for val in vals:
        if val and val not in seen:
            out.append(val)
            seen.add(val)
    return tuple(out)


def contig_aliases(name: str) -> List[str]:
    n = str(name).strip()
    out: List[str] = []
    seen = set()

    def add(v: str) -> None:
        v = str(v).strip()
        if v and v not in seen:
            out.append(v)
            seen.add(v)

    add(n)
    add(n[3:] if n.startswith("chr") else "chr" + n)
    add(n.split(".")[0])
    return out


def resolve_contig(name: str, seqs: Dict[str, str]) -> Optional[str]:
    for alias in contig_aliases(name):
        if alias in seqs:
            return alias
    if len(seqs) == 1:
        return next(iter(seqs))
    return None


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def seq_id_for(seq: str) -> str:
    return hashlib.sha1(seq.encode()).hexdigest()[:16]


def load_target_clumps(path: Path, ref_seqs: Dict[str, str]) -> List[TargetClump]:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    rows: List[TargetClump] = []
    for _, row in df.iterrows():
        contig_raw = str(row["reference_contig"]).strip()
        contig = resolve_contig(contig_raw, ref_seqs) or contig_raw
        fallback_anchor_raw = str(row.get("fallback_anchor_pos", "")).strip()
        fallback_anchor_pos = int(fallback_anchor_raw) if fallback_anchor_raw and fallback_anchor_raw.lower() != "nan" else None
        rows.append(
            TargetClump(
                target_id=str(row["target_id"]).strip(),
                target_type=str(row["target_type"]).strip(),
                reference_contig=contig,
                target_start=int(str(row["target_start"]).strip()),
                target_end=int(str(row["target_end"]).strip()),
                reference_strand=str(row.get("reference_strand", ".")).strip() or ".",
                reference_gene_id=str(row.get("reference_gene_id", "")).strip(),
                reference_gene_symbol=str(row.get("reference_gene_symbol", "")).strip(),
                reference_product_name=str(row.get("reference_product_name", "")).strip(),
                fallback_anchor_pos=fallback_anchor_pos,
                anchors_all=parse_anchor_list(row.get("anchors_all", "")),
                anchors_genic=parse_anchor_list(row.get("anchors_genic", "")),
                anchors_upstream200=parse_anchor_list(row.get("anchors_upstream200", "")),
            )
        )
    return rows


def parse_chain_file(path: Path) -> Dict[str, List[ChainRecord]]:
    chains_by_t_name: Dict[str, List[ChainRecord]] = defaultdict(list)
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        if not line.startswith("chain "):
            raise RuntimeError(f"Unexpected chain line in {path}: {line}")
        parts = line.split()
        if len(parts) != 13:
            raise RuntimeError(f"Malformed chain header in {path}: {line}")

        _kw, _score, t_name, t_size, t_strand, t_start, t_end, q_name, q_size, q_strand, q_start, q_end, chain_id = parts
        t_size_i = int(t_size)
        t_start_i = int(t_start)
        t_end_i = int(t_end)
        q_size_i = int(q_size)
        q_start_i = int(q_start)
        q_end_i = int(q_end)

        t_pos = t_start_i
        q_pos_raw = q_start_i if q_strand == "+" else (q_size_i - q_end_i)
        blocks: List[ChainBlock] = []
        while i < len(lines):
            row = lines[i].strip()
            i += 1
            if not row:
                break
            nums = row.split()
            if len(nums) not in {1, 3}:
                raise RuntimeError(f"Malformed chain block line in {path}: {row}")
            size = int(nums[0])
            blocks.append(ChainBlock(t_start=t_pos, size=size, q_start_raw=q_pos_raw))
            if len(nums) == 1:
                break
            t_pos += size + int(nums[1])
            q_pos_raw += size + int(nums[2])

        chains_by_t_name[t_name].append(
            ChainRecord(
                t_name=t_name,
                t_size=t_size_i,
                t_strand=t_strand,
                t_start=t_start_i,
                t_end=t_end_i,
                q_name=q_name,
                q_size=q_size_i,
                q_strand=q_strand,
                q_start=q_start_i,
                q_end=q_end_i,
                chain_id=chain_id,
                blocks=blocks,
            )
        )

    for key in list(chains_by_t_name.keys()):
        chains_by_t_name[key].sort(key=lambda rec: (rec.t_start, rec.t_end, rec.chain_id))
    return chains_by_t_name


def choose_chain(chains_by_t_name: Dict[str, List[ChainRecord]], ref_contig: str, start_1based: int, end_1based: int) -> Optional[ChainRecord]:
    t0 = int(start_1based) - 1
    t1 = int(end_1based)
    candidates: List[ChainRecord] = []
    for alias in contig_aliases(ref_contig):
        candidates.extend(chains_by_t_name.get(alias, []))
    best: Optional[ChainRecord] = None
    best_overlap = 0
    for rec in candidates:
        if rec.t_strand != "+":
            continue
        overlap = max(0, min(t1, rec.t_end) - max(t0, rec.t_start))
        if overlap > best_overlap:
            best = rec
            best_overlap = overlap
    return best


def map_ref_pos_to_consensus(chains_by_t_name: Dict[str, List[ChainRecord]], ref_contig: str, pos_1based: int) -> Optional[Tuple[str, int]]:
    pos0 = int(pos_1based) - 1
    for alias in contig_aliases(ref_contig):
        for rec in chains_by_t_name.get(alias, []):
            if rec.t_strand != "+" or not (rec.t_start <= pos0 < rec.t_end):
                continue
            for block in rec.blocks:
                if block.t_start <= pos0 < block.t_start + block.size:
                    offset = pos0 - block.t_start
                    if rec.q_strand == "+":
                        q_pos0 = block.q_start_raw + offset
                    else:
                        q_pos0 = rec.q_size - 1 - (block.q_start_raw + offset)
                    return rec.q_name, q_pos0 + 1
    return None


def extract_fixed_window(consensus_seqs: Dict[str, str], consensus_contig: str, anchor_pos_1based: int, window_bp: int) -> Tuple[str, int, int]:
    target_len = 2 * window_bp + 1
    contig = resolve_contig(consensus_contig, consensus_seqs)
    if contig is None:
        return "N" * target_len, target_len, target_len
    seq = consensus_seqs[contig]
    anchor_pos0 = int(anchor_pos_1based) - 1
    left = anchor_pos0 - window_bp
    right = anchor_pos0 + window_bp + 1
    pad_left = max(0, -left)
    pad_right = max(0, right - len(seq))
    left = max(0, left)
    right = min(len(seq), right)
    window = ("N" * pad_left) + seq[left:right] + ("N" * pad_right)
    if len(window) < target_len:
        window += "N" * (target_len - len(window))
    elif len(window) > target_len:
        window = window[:target_len]
    return window, pad_left, pad_right


def fallback_anchor_pos(target: TargetClump) -> int:
    if target.fallback_anchor_pos is not None:
        return int(target.fallback_anchor_pos)
    return (target.target_start + target.target_end) // 2


def fallback_window_bp(target: TargetClump) -> int:
    if target.target_length % 2 == 1:
        return target.target_length // 2
    return 25


def project_interval_to_consensus(
    *,
    chains_by_t_name: Dict[str, List[ChainRecord]],
    consensus_seqs: Dict[str, str],
    target: TargetClump,
) -> dict:
    rec = choose_chain(chains_by_t_name, target.reference_contig, target.target_start, target.target_end)
    if rec is None:
        return {
            "sequence": "N" * target.target_length,
            "sequence_length": target.target_length,
            "projection_status": "mapping_failed",
            "consensus_contig": "",
            "projected_start": 0,
            "projected_end": 0,
            "projected_inserted_bp": 0,
            "projected_deleted_ref_bp": target.target_length,
            "n_chain_blocks": 0,
        }
    if rec.q_strand != "+":
        return {
            "sequence": "N" * target.target_length,
            "sequence_length": target.target_length,
            "projection_status": "unsupported_query_strand",
            "consensus_contig": rec.q_name,
            "projected_start": 0,
            "projected_end": 0,
            "projected_inserted_bp": 0,
            "projected_deleted_ref_bp": target.target_length,
            "n_chain_blocks": 0,
        }

    consensus_contig = resolve_contig(rec.q_name, consensus_seqs)
    if consensus_contig is None:
        return {
            "sequence": "N" * target.target_length,
            "sequence_length": target.target_length,
            "projection_status": "missing_consensus_contig",
            "consensus_contig": rec.q_name,
            "projected_start": 0,
            "projected_end": 0,
            "projected_inserted_bp": 0,
            "projected_deleted_ref_bp": target.target_length,
            "n_chain_blocks": 0,
        }

    t0 = int(target.target_start) - 1
    t1 = int(target.target_end)
    consensus_seq = consensus_seqs[consensus_contig]
    parts: List[str] = []
    q_coords: List[int] = []
    covered_ref_bp = 0
    inserted_bp = 0
    n_blocks = 0
    prev_block: Optional[ChainBlock] = None

    for block in rec.blocks:
        b0 = block.t_start
        b1 = block.t_start + block.size
        if prev_block is not None:
            prev_t_end = prev_block.t_start + prev_block.size
            prev_q_end = prev_block.q_start_raw + prev_block.size
            q_gap_start = prev_q_end
            q_gap_end = block.q_start_raw
            if q_gap_end > q_gap_start and t0 <= prev_t_end <= t1:
                ins = consensus_seq[q_gap_start:q_gap_end]
                parts.append(ins)
                q_coords.extend([q_gap_start, q_gap_end])
                inserted_bp += len(ins)

        if b1 <= t0:
            prev_block = block
            continue
        if b0 >= t1:
            break

        o0 = max(t0, b0)
        o1 = min(t1, b1)
        if o1 > o0:
            q0 = block.q_start_raw + (o0 - b0)
            q1 = block.q_start_raw + (o1 - b0)
            parts.append(consensus_seq[q0:q1])
            q_coords.extend([q0, q1])
            covered_ref_bp += o1 - o0
            n_blocks += 1
        prev_block = block

    sequence = "".join(parts).upper()
    if not sequence:
        return {
            "sequence": "N" * target.target_length,
            "sequence_length": target.target_length,
            "projection_status": "mapping_failed",
            "consensus_contig": consensus_contig,
            "projected_start": 0,
            "projected_end": 0,
            "projected_inserted_bp": 0,
            "projected_deleted_ref_bp": target.target_length,
            "n_chain_blocks": 0,
        }

    deleted_bp = max(0, target.target_length - covered_ref_bp)
    return {
        "sequence": sequence,
        "sequence_length": len(sequence),
        "projection_status": "mapped" if deleted_bp == 0 else "mapped_with_ref_deletion",
        "consensus_contig": consensus_contig,
        "projected_start": min(q_coords) + 1 if q_coords else 0,
        "projected_end": max(q_coords) if q_coords else 0,
        "projected_inserted_bp": int(inserted_bp),
        "projected_deleted_ref_bp": int(deleted_bp),
        "n_chain_blocks": int(n_blocks),
    }


def choose_vcf_path(sample_dir: Path, vcf_gz_name: str, vcf_name: str) -> Optional[Path]:
    gz = sample_dir / vcf_gz_name
    plain = sample_dir / vcf_name
    if gz.is_file():
        return gz
    if plain.is_file():
        return plain
    return None


def generate_consensus_and_chain(
    *,
    sample_dir: Path,
    reference_fasta: Path,
    bcftools_bin: Optional[Path],
    consensus_name: str,
    chain_name: str,
    vcf_gz_name: str,
    vcf_name: str,
    regenerate: bool,
) -> Dict[str, object]:
    consensus_path = sample_dir / consensus_name
    chain_path = sample_dir / chain_name
    if consensus_path.is_file() and chain_path.is_file() and not regenerate:
        return {"status": "ready_existing", "consensus_path": consensus_path, "chain_path": chain_path}
    vcf_path = choose_vcf_path(sample_dir, vcf_gz_name, vcf_name)
    if vcf_path is None:
        return {"status": "missing_vcf", "consensus_path": consensus_path, "chain_path": chain_path}
    if bcftools_bin is None:
        return {"status": "missing_bcftools_bin", "consensus_path": consensus_path, "chain_path": chain_path}

    stdout_log = sample_dir / "snps.consensus.chain.stdout.log"
    stderr_log = sample_dir / "snps.consensus.chain.stderr.log"
    cmd = [
        str(bcftools_bin),
        "consensus",
        "-f",
        str(reference_fasta),
        "-c",
        str(chain_path),
        "-o",
        str(consensus_path),
        str(vcf_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    stdout_log.write_text(proc.stdout or "")
    stderr_log.write_text(proc.stderr or "")
    if proc.returncode != 0:
        return {"status": "consensus_failed", "returncode": int(proc.returncode), "consensus_path": consensus_path, "chain_path": chain_path}
    if not consensus_path.is_file() or not chain_path.is_file():
        return {"status": "consensus_missing_outputs", "consensus_path": consensus_path, "chain_path": chain_path}
    return {"status": "ready_regenerated", "consensus_path": consensus_path, "chain_path": chain_path}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-clumps", required=True)
    ap.add_argument("--sample-root", required=True)
    ap.add_argument("--reference-fasta", required=True)
    ap.add_argument("--sample-ids", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--consensus-name", default="snps.consensus.fa")
    ap.add_argument("--chain-name", default="snps.consensus.chain")
    ap.add_argument("--sample-vcf-gz-name", default="snps.vcf.gz")
    ap.add_argument("--sample-vcf-name", default="snps.vcf")
    ap.add_argument("--bcftools-bin", default="")
    ap.add_argument("--regenerate-consensus", action="store_true")
    args = ap.parse_args()

    target_clumps_path = Path(args.target_clumps)
    sample_root = Path(args.sample_root)
    reference_fasta = Path(args.reference_fasta)
    sample_ids_path = Path(args.sample_ids)
    outdir = Path(args.outdir)
    bcftools_bin = Path(args.bcftools_bin) if str(args.bcftools_bin).strip() else None
    outdir.mkdir(parents=True, exist_ok=True)

    sample_tsv_path = outdir / "sample_target_sequences.tsv"
    unique_fasta_path = outdir / "unique_window_sequences.fasta"
    unique_tsv_path = outdir / "unique_window_sequences.tsv"
    per_sample_stats_path = outdir / "per_sample_stats.tsv"
    manifest_path = outdir / "13l_build_snp_gene_consensus_chain_sequences.manifest.json"

    sample_ids = read_ids(sample_ids_path)
    ref_seqs = read_fasta(reference_fasta)
    targets = load_target_clumps(target_clumps_path, ref_seqs)
    if not sample_ids:
        raise RuntimeError("No sample IDs loaded.")
    if not ref_seqs:
        raise RuntimeError("Reference FASTA is empty or unreadable.")
    if not targets:
        raise RuntimeError("No target clumps loaded.")

    seq_counts: Dict[str, int] = {}
    seq_by_id: Dict[str, str] = {}
    per_sample_rows: List[dict] = []
    row_type_counts: Counter[str] = Counter()
    projection_status_counts: Counter[str] = Counter()
    consensus_status_counts: Counter[str] = Counter()
    n_missing_sample_dirs = 0
    n_missing_vcfs = 0
    n_consensus_failures = 0
    n_projected_inserted_bp = 0
    n_projected_deleted_ref_bp = 0

    with sample_tsv_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "target_id",
                "target_type",
                "reference_contig",
                "target_start",
                "target_end",
                "reference_strand",
                "reference_gene_id",
                "reference_gene_symbol",
                "reference_product_name",
                "anchors_all",
                "anchors_genic",
                "anchors_upstream200",
                "sequence",
                "sequence_length",
                "n_applied_variants",
                "n_applied_indels",
                "n_conflicts",
                "projection_status",
                "consensus_status",
                "consensus_contig",
                "projected_start",
                "projected_end",
                "projected_inserted_bp",
                "projected_deleted_ref_bp",
                "n_chain_blocks",
                "source",
            ]
        )

        for sample_id in sample_ids:
            sample_dir = sample_root / sample_id
            if not sample_dir.is_dir():
                n_missing_sample_dirs += 1
                consensus_status = "missing_sample_dir"
                consensus_seqs: Dict[str, str] = {}
                chains_by_t_name: Dict[str, List[ChainRecord]] = {}
            else:
                consensus_result = generate_consensus_and_chain(
                    sample_dir=sample_dir,
                    reference_fasta=reference_fasta,
                    bcftools_bin=bcftools_bin,
                    consensus_name=args.consensus_name,
                    chain_name=args.chain_name,
                    vcf_gz_name=args.sample_vcf_gz_name,
                    vcf_name=args.sample_vcf_name,
                    regenerate=bool(args.regenerate_consensus),
                )
                consensus_status = str(consensus_result["status"])
                if consensus_status == "missing_vcf":
                    n_missing_vcfs += 1
                    consensus_seqs = {}
                    chains_by_t_name = {}
                elif consensus_status.startswith("ready_"):
                    try:
                        consensus_seqs = read_fasta(Path(consensus_result["consensus_path"]))
                        chains_by_t_name = parse_chain_file(Path(consensus_result["chain_path"]))
                    except Exception:
                        n_consensus_failures += 1
                        consensus_status = "consensus_parse_failed"
                        consensus_seqs = {}
                        chains_by_t_name = {}
                else:
                    n_consensus_failures += 1
                    consensus_seqs = {}
                    chains_by_t_name = {}

            consensus_status_counts[consensus_status] += 1
            sample_counts: Counter[str] = Counter()
            sample_status_counts: Counter[str] = Counter()
            sample_inserted_bp = 0
            sample_deleted_ref_bp = 0

            for target in targets:
                if consensus_seqs and chains_by_t_name:
                    if target.target_type.startswith("fallback_w"):
                        anchor_pos = fallback_anchor_pos(target)
                        mapped = map_ref_pos_to_consensus(chains_by_t_name, target.reference_contig, anchor_pos)
                        if mapped is None:
                            result = {
                                "sequence": "N" * target.target_length,
                                "sequence_length": target.target_length,
                                "projection_status": "mapping_failed",
                                "consensus_contig": "",
                                "projected_start": 0,
                                "projected_end": 0,
                                "projected_inserted_bp": 0,
                                "projected_deleted_ref_bp": target.target_length,
                                "n_chain_blocks": 0,
                            }
                        else:
                            q_name, q_pos = mapped
                            seq, pad_left, pad_right = extract_fixed_window(
                                consensus_seqs=consensus_seqs,
                                consensus_contig=q_name,
                                anchor_pos_1based=q_pos,
                                window_bp=fallback_window_bp(target),
                            )
                            result = {
                                "sequence": seq,
                                "sequence_length": len(seq),
                                "projection_status": "mapped" if not (pad_left or pad_right) else "mapped_edge_padded",
                                "consensus_contig": q_name,
                                "projected_start": max(1, q_pos - fallback_window_bp(target)),
                                "projected_end": q_pos + fallback_window_bp(target),
                                "projected_inserted_bp": 0,
                                "projected_deleted_ref_bp": 0,
                                "n_chain_blocks": 1,
                            }
                    else:
                        result = project_interval_to_consensus(
                            chains_by_t_name=chains_by_t_name,
                            consensus_seqs=consensus_seqs,
                            target=target,
                        )
                else:
                    result = {
                        "sequence": "N" * target.target_length,
                        "sequence_length": target.target_length,
                        "projection_status": "missing_consensus",
                        "consensus_contig": "",
                        "projected_start": 0,
                        "projected_end": 0,
                        "projected_inserted_bp": 0,
                        "projected_deleted_ref_bp": target.target_length,
                        "n_chain_blocks": 0,
                    }

                seq = str(result["sequence"]).upper()
                if target.reference_strand == "-":
                    seq = reverse_complement(seq)

                inserted_bp = int(result["projected_inserted_bp"])
                deleted_ref_bp = int(result["projected_deleted_ref_bp"])
                n_projected_inserted_bp += inserted_bp
                n_projected_deleted_ref_bp += deleted_ref_bp
                projection_status = str(result["projection_status"])
                projection_status_counts[projection_status] += 1
                sample_status_counts[projection_status] += 1

                writer.writerow(
                    [
                        sample_id,
                        target.target_id,
                        target.target_type,
                        target.reference_contig,
                        target.target_start,
                        target.target_end,
                        target.reference_strand,
                        target.reference_gene_id,
                        target.reference_gene_symbol,
                        target.reference_product_name,
                        json.dumps(list(target.anchors_all)),
                        json.dumps(list(target.anchors_genic)),
                        json.dumps(list(target.anchors_upstream200)),
                        seq,
                        len(seq),
                        0,
                        0,
                        0,
                        projection_status,
                        consensus_status,
                        str(result["consensus_contig"]),
                        int(result["projected_start"]),
                        int(result["projected_end"]),
                        inserted_bp,
                        deleted_ref_bp,
                        int(result["n_chain_blocks"]),
                        "consensus_chain_13l",
                    ]
                )

                sample_counts[target.target_type] += 1
                row_type_counts[target.target_type] += 1
                sample_inserted_bp += inserted_bp
                sample_deleted_ref_bp += deleted_ref_bp
                seq_id = seq_id_for(seq)
                seq_counts[seq_id] = seq_counts.get(seq_id, 0) + 1
                seq_by_id.setdefault(seq_id, seq)

            sample_row = {
                "sample_id": sample_id,
                "consensus_status": consensus_status,
                "total_rows": int(sum(sample_counts.values())),
                "projected_inserted_bp": int(sample_inserted_bp),
                "projected_deleted_ref_bp": int(sample_deleted_ref_bp),
            }
            for target_type in sorted({t.target_type for t in targets}):
                sample_row[target_type] = int(sample_counts[target_type])
            for status, count in sorted(sample_status_counts.items()):
                sample_row[f"projection_status_{status}"] = int(count)
            per_sample_rows.append(sample_row)

    unique_rows = [(seq_id, seq_by_id[seq_id]) for seq_id in seq_counts.keys()]
    write_fasta(unique_fasta_path, unique_rows)
    pd.DataFrame(
        [{"seq_id": sid, "sequence": seq_by_id[sid], "count": int(seq_counts[sid])} for sid in seq_counts.keys()]
    ).to_csv(unique_tsv_path, sep="\t", index=False)
    pd.DataFrame(per_sample_rows).to_csv(per_sample_stats_path, sep="\t", index=False)

    manifest = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": "ok",
        "args": {
            "target_clumps": str(target_clumps_path),
            "sample_root": str(sample_root),
            "reference_fasta": str(reference_fasta),
            "sample_ids": str(sample_ids_path),
            "outdir": str(outdir),
            "consensus_name": str(args.consensus_name),
            "chain_name": str(args.chain_name),
            "sample_vcf_gz_name": str(args.sample_vcf_gz_name),
            "sample_vcf_name": str(args.sample_vcf_name),
            "regenerate_consensus": bool(args.regenerate_consensus),
        },
        "n_targets": len(targets),
        "target_type_counts": dict(Counter(t.target_type for t in targets)),
        "n_samples": len(sample_ids),
        "n_sample_target_sequences": int(len(sample_ids) * len(targets)),
        "row_type_counts": {k: int(v) for k, v in sorted(row_type_counts.items())},
        "projection_status_counts": {k: int(v) for k, v in sorted(projection_status_counts.items())},
        "consensus_status_counts": {k: int(v) for k, v in sorted(consensus_status_counts.items())},
        "n_unique_sequences": len(seq_counts),
        "n_missing_sample_dirs": int(n_missing_sample_dirs),
        "n_missing_sample_vcfs": int(n_missing_vcfs),
        "n_consensus_generation_failures": int(n_consensus_failures),
        "n_projected_inserted_bp": int(n_projected_inserted_bp),
        "n_projected_deleted_ref_bp": int(n_projected_deleted_ref_bp),
        "output_files": {
            "sample_target_sequences_tsv": str(sample_tsv_path),
            "unique_window_sequences_tsv": str(unique_tsv_path),
            "unique_window_sequences_fasta": str(unique_fasta_path),
            "per_sample_stats_tsv": str(per_sample_stats_path),
            "manifest_json": str(manifest_path),
        },
    }
    write_json(manifest_path, manifest)

    print("[OK] built consensus-chain target sequences")
    print(f"[OK] targets={len(targets)} samples={len(sample_ids)} total_rows={int(len(sample_ids) * len(targets))} unique_sequences={len(seq_counts)}")
    print(f"[OK] projection_status_counts={dict(projection_status_counts)}")
    print(f"[OK] outdir={outdir}")


if __name__ == "__main__":
    main()
