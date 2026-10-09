#!/usr/bin/env python3
"""Acquired-gene sequences from each isolate's own assembly, via AMRFinder hits.

The chromosomal extractor works from PAO1 coordinates, which acquired genes do
not have -- they are absent from the reference entirely. AMRFinder already
located them on each isolate's assembly and reports `Contig id`, `Start`, `Stop`
and `Strand`, so this is a slice rather than a search.

Three choices worth stating, because each could have gone another way:

  **No upstream window.** The chromosomal targets carry gene + 300 bp upstream so
  promoter variation is visible. Acquired genes arrive on mobile elements with
  their own promoters, and in a draft assembly a hit often sits near a contig
  edge -- so a 300 bp request would be silently truncated for some isolates and
  not others, encoding assembly quality as if it were biology. The CDS as
  AMRFinder delimits it is well defined for every hit.

  **One representative per (isolate, family), not a mean of copies.** An isolate
  can carry a family more than once -- integron cassettes, plasmid plus
  chromosome. The existing P. aeruginosa pipeline mean-pooled copies; this takes
  the highest-identity hit instead and keeps **copy number as its own column**.
  Copies within a family are near-identical, so the mean buys little, and one row
  per cell keeps the map format identical to every other target. Where the family
  collapse means an isolate holds two genuinely different alleles (blaVIM-1 and
  blaVIM-2), the copy-number column marks it rather than the mean hiding it.

  **Families, not alleles.** The same collapse the presence matrix uses, so the
  two blocks line up: `amr_blaVIM` as a binary, `amrseq_blaVIM` as a sequence.

Absence is the normal state here -- most isolates lack most families -- so the
map is deliberately sparse and the precompute needs a fill policy for these
targets. That is the consumer's decision, not this script's; it emits only what
exists.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import pathlib
import sys

COLLAPSE_FAMILIES = ("blaVIM", "blaIMP", "blaNDM", "blaGES", "blaKPC", "blaSPM",
                     "blaSIM", "blaAIM", "blaDIM", "blaGIM", "blaCTX-M",
                     "blaPER", "blaVEB")
LOF_TOKENS = ("STOP", "FS", "DEL", "INS", "DUP")


def open_text(p: pathlib.Path):
    with open(p, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(p, "rt") if magic == b"\x1f\x8b" else open(p, "r")


def read_fasta(p: pathlib.Path) -> dict[str, str]:
    seqs, name, buf = {}, "", []
    with open_text(p) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return {k: v.upper() for k, v in seqs.items()}


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNRYKMSWBVDH", "TGCANYRMKSWVBHD"))[::-1]


def family_of(sym: str, subtype: str, hier: str) -> str:
    if subtype == "POINT":
        gene, _, change = sym.partition("_")
        return (f"{gene}_LOF" if gene and any(k in change.upper() for k in LOF_TOKENS)
                else sym)
    fam = hier if hier and hier != "NA" else sym
    for base in COLLAPSE_FAMILIES:
        if fam == base or fam.startswith(base + "-"):
            return base
    return fam


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--amrfinder-dir", required=True)
    ap.add_argument("--assembly-map", required=True,
                    help="TSV of sample_id<TAB>assembly_path")
    ap.add_argument("--families", required=True,
                    help="one family per line, with or without an 'amr_' prefix")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-len", type=int, default=90,
                    help="a slice shorter than this is a coordinate error, not a gene")
    args = ap.parse_args()

    want = {l.strip()[4:] if l.strip().startswith("amr_") else l.strip()
            for l in open(args.families) if l.strip()}
    print(f"[INFO] {len(want)} families wanted")

    asm = {}
    for line in open(args.assembly_map):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] != "sample_id":
            asm[f[0]] = f[1]
    print(f"[INFO] {len(asm):,} assembly paths")

    files = sorted(pathlib.Path(args.amrfinder_dir).rglob("*.amrfinderplus.tsv"))
    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows, copies, seqs = [], [], {}
    no_asm, no_contig, short = [], 0, 0
    for i, f in enumerate(files, 1):
        sample = f.name[: -len(".amrfinderplus.tsv")]
        best: dict[str, tuple] = {}
        n_copy: collections.Counter = collections.Counter()
        with f.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                fam = family_of((r.get("Element symbol") or "").strip(),
                                (r.get("Subtype") or "").strip().upper(),
                                (r.get("Hierarchy node") or "").strip())
                if fam not in want:
                    continue
                try:
                    ident = float(r.get("% Identity to reference") or 0)
                    alen = int(float(r.get("Alignment length") or 0))
                    start, stop = int(r["Start"]), int(r["Stop"])
                except (TypeError, ValueError, KeyError):
                    continue
                n_copy[fam] += 1
                key = (ident, alen)
                if fam not in best or key > best[fam][0]:
                    best[fam] = (key, r["Contig id"], start, stop,
                                 (r.get("Strand") or "+").strip())
        if not best:
            continue
        path = asm.get(sample)
        if not path or not pathlib.Path(path).exists():
            no_asm.append(sample)
            continue
        contigs = read_fasta(pathlib.Path(path))
        for fam, (_, contig, start, stop, strand) in best.items():
            if contig not in contigs:
                no_contig += 1
                continue
            s = contigs[contig][min(start, stop) - 1: max(start, stop)]
            if len(s) < args.min_len:
                short += 1
                continue
            if strand == "-":
                s = revcomp(s)
            sid = "sha1:" + hashlib.sha1(s.encode()).hexdigest()
            seqs[sid] = s
            rows.append((sample, f"amrseq_{fam}", sid))
            copies.append((sample, fam, n_copy[fam]))
        if i % 500 == 0:
            print(f"[INFO] {i:,}/{len(files):,} isolates")

    if no_asm:
        print(f"[WARN] {len(no_asm):,} isolate(s) had hits but no assembly, "
              f"e.g. {no_asm[:3]}")
    if no_contig:
        print(f"[WARN] {no_contig:,} hit(s) named a contig absent from the assembly")
    if short:
        print(f"[WARN] {short:,} slice(s) under {args.min_len} bp, dropped")

    with (out / "sample_target_seqid_map.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "target_id", "seq_id"])
        w.writerows(rows)
    with (out / "unique_window_sequences.fasta").open("w") as fh:
        for sid in sorted(seqs):
            fh.write(f">{sid}\n{seqs[sid]}\n")
    with (out / "copy_number.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "family", "copies"])
        w.writerows(copies)

    per_fam = collections.Counter(t for _, t, _ in rows)
    multi = sum(1 for _, _, n in copies if n > 1)
    print(f"\n[INFO] {len(rows):,} cells, {len(seqs):,} distinct sequences "
          f"({len(seqs)/max(1,len(rows)):.1%} of the naive count)")
    print(f"[INFO] {multi:,} cells where the isolate carries >1 copy")
    print(f"\n{'target':<28}{'isolates':>9}{'distinct':>10}")
    byt = collections.defaultdict(set)
    for _, t, sid in rows:
        byt[t].add(sid)
    for t, n in per_fam.most_common(15):
        print(f"{t:<28}{n:>9}{len(byt[t]):>10}")
    print(f"\n[OK] {out}")
    print("The map is sparse on purpose: an isolate without the gene has no row. "
          "The\nprecompute needs a fill policy for these targets -- that is a "
          "modelling choice,\nnot an extraction one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
