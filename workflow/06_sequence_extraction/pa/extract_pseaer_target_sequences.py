#!/usr/bin/env python3
"""Per-sample target sequences from the PAO1 reference and an allele-filtered VCF.

The MTB extractor reads TB-Profiler's own per-sample VCF and `tbdb.fasta`, so it
does not port. This does the same job from the two things P. aeruginosa actually
has: the PAO1 reference and a Snippy VCF rebuilt at a chosen allele fraction.

Variants are applied **inside each target window only**, never genome-wide. That
avoids needing a coordinate chain: an indel inside a window shifts the rest of
that window, which is the correct sequence for that region, while an indel
elsewhere cannot silently move a window's boundaries. A variant straddling a
boundary is clipped to the window and counted, because dropping it would quietly
return reference sequence for a site that is not reference.

Emits, per sample, one row per target with the sha1 of its sequence, plus the
sequences themselves. A gather step deduplicates across samples so each distinct
sequence is embedded once -- the MTB layout, which every downstream script
already expects:

    unique_window_sequences.fasta      seq_id -> sequence
    sample_target_seqid_map.tsv        sample_id, target_id, seq_id

Sites failing the filter are left as reference, which is what a Snippy consensus
does. The count of applied and skipped variants is written per sample so a
systematically empty sample is visible rather than silently identical to PAO1.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import pathlib
import sys


def open_text(p: pathlib.Path):
    """Decide by magic bytes. `filter_vcf_by_allele_fraction.py` writes plain
    text whatever the filename says, so trusting a `.gz` suffix hands gzip a
    text file and fails with BadGzipFile two steps from the actual cause."""
    with open(p, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(p, "rt") if magic == b"\x1f\x8b" else open(p, "r")


def read_fasta(p: pathlib.Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = ""
    with open_text(p) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def load_targets(clumps: pathlib.Path, wanted: set[str] | None):
    out = []
    with clumps.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if wanted is not None and r["target_id"] not in wanted:
                continue
            out.append((r["target_id"], r["reference_contig"],
                        int(r["target_start"]), int(r["target_end"]),
                        r.get("reference_strand", "+")))
    return out


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, help="PAO1 .fna")
    ap.add_argument("--clumps", required=True, help="all_gene_target_clumps.tsv")
    ap.add_argument("--target-list", default="",
                    help="restrict to these target_ids, one per line")
    ap.add_argument("--vcf", required=True, help="allele-filtered VCF for this sample")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--revcomp-minus", action="store_true",
                    help="return minus-strand targets in gene orientation. Off by "
                         "default: the reference strand is a fixed property of the "
                         "target, so every isolate is consistent either way, and "
                         "flipping only matters if the embedding model is read as "
                         "seeing a coding sequence")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    ref = read_fasta(pathlib.Path(args.reference))
    wanted = None
    if args.target_list:
        wanted = {l.split("\t")[0].strip() for l in open(args.target_list)
                  if l.strip() and not l.startswith("#")}
        wanted.discard("target_id")
    targets = load_targets(pathlib.Path(args.clumps), wanted)
    if wanted is not None and len(targets) != len(wanted):
        got = {t[0] for t in targets}
        raise SystemExit(f"[ERROR] {len(wanted) - len(targets)} requested target(s) "
                         f"absent from the clumps file, e.g. "
                         f"{sorted(wanted - got)[:3]}")
    if not targets:
        raise SystemExit("[ERROR] no targets selected")

    # variants per contig, sorted; only substitutions/indels with a called ALT
    #
    # Snippy names the contig after the reference FASTA header; Bakta renames it
    # in its own annotation, and the clumps file inherits Bakta's name. The two
    # are the same sequence under different labels, so a mismatch produces a run
    # that applies zero variants and returns 62 copies of PAO1 -- a plausible
    # result that is entirely reference. Rather than trust a rename, the VCF's
    # declared contig length is matched against the reference's actual length.
    var: dict[str, list[tuple[int, str, str]]] = collections.defaultdict(list)
    vcf_len: dict[str, int] = {}
    n_var = 0
    with open_text(pathlib.Path(args.vcf)) as fh:
        for line in fh:
            if line.startswith("##contig="):
                body = line.strip()[len("##contig=<"):].rstrip(">")
                kv = dict(x.split("=", 1) for x in body.split(",") if "=" in x)
                if "ID" in kv and "length" in kv:
                    vcf_len[kv["ID"]] = int(kv["length"])
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[4] in (".", ""):
                continue
            alt = f[4].split(",")[0]
            var[f[0]].append((int(f[1]), f[3].upper(), alt.upper()))
            n_var += 1
    for c in var:
        var[c].sort()

    # remap VCF contig names onto reference names, by length
    unknown = [c for c in var if c not in ref]
    if unknown:
        by_len = collections.defaultdict(list)
        for name, seq in ref.items():
            by_len[len(seq)].append(name)
        for c in unknown:
            n = vcf_len.get(c)
            cand = by_len.get(n, []) if n else []
            if len(cand) != 1:
                raise SystemExit(
                    f"[ERROR] VCF contig {c!r} (length {n}) does not match any "
                    f"single reference contig.\n"
                    f"         reference: "
                    + ", ".join(f"{k} ({len(v):,})" for k, v in ref.items())
                    + "\n         Pass a reference whose sequence matches the VCF.")
            print(f"[INFO] contig {c!r} -> {cand[0]!r} (both {n:,} bp)")
            var[cand[0]] = var.pop(c)

    applied = clipped = 0
    rows, seqs = [], {}
    for tid, contig, start, end, strand in targets:
        if contig not in ref:
            raise SystemExit(f"[ERROR] contig {contig!r} not in the reference; "
                             f"have {list(ref)[:3]}")
        chrom = ref[contig]
        start = max(1, start)
        end = min(len(chrom), end)
        pieces, cur = [], start
        for pos, r, a in var.get(contig, ()):
            if pos + len(r) - 1 < start or pos > end:
                continue
            if pos < cur:            # overlapping the previous edit
                continue
            pieces.append(chrom[cur - 1:pos - 1])
            ref_end = pos + len(r) - 1
            if ref_end > end:        # straddles the window edge
                keep = end - pos + 1
                pieces.append(a[:max(1, keep)])
                clipped += 1
                cur = end + 1
                break
            pieces.append(a)
            applied += 1
            cur = ref_end + 1
        pieces.append(chrom[cur - 1:end])
        s = "".join(pieces)
        if args.revcomp_minus and strand == "-":
            s = revcomp(s)
        sid = "sha1:" + hashlib.sha1(s.encode()).hexdigest()
        seqs[sid] = s
        rows.append((args.sample_id, tid, sid, len(s)))

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with (out / f"{args.sample_id}.targets.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "target_id", "seq_id", "length"])
        w.writerows(rows)
    with (out / f"{args.sample_id}.sequences.fasta").open("w") as fh:
        for sid, s in seqs.items():
            fh.write(f">{sid}\n{s}\n")

    print(f"[OK] {args.sample_id}: {len(rows)} targets, {len(seqs)} distinct "
          f"sequences, {applied:,} variants applied, {clipped} clipped at a "
          f"window edge, {n_var:,} in the VCF")
    if applied == 0 and n_var > 0:
        print("[WARN] no variant fell inside any target window — check that the "
              "VCF contig\n       name matches the reference")
    return 0


if __name__ == "__main__":
    sys.exit(main())
