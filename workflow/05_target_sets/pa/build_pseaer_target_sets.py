#!/usr/bin/env python3
"""Per-drug target sets for P. aeruginosa, resolved against the PAO1 annotation.

The MTB panel's single most expensive mistake was treating a catalogue's gene
list as a target set: `tbdb.bed` contains no *fabG1*, so isoniazid models built
from the catalogue's genes missed 23% of resistant isolates because the *inhA*
promoter lives in an interval the gene list never names. The same failure is
waiting here in a different form, and three instances of it are already
confirmed in this annotation:

  ***ampC*** carries `gene=blaPDC-1`. Bakta names it by the PDC allele, so an
  exact match on "ampC" returns nothing and four beta-lactam drugs silently lose
  their principal locus.
  ***mexZ*** carries `gene=mexZ_G195E` -- a *variant allele suffix*. Exact
  matching finds nothing and both aminoglycosides lose their efflux regulator.
  ***oprD*** matches four genes. Three are OprD-family paralogs; only
  NFOBNJ_00992 is the porin whose loss drives carbapenem resistance.

And only 2,701 of 5,682 CDS carry a `gene=` attribute at all, so name-keying
alone cannot work for this annotation.

So every locus is declared with an explicit resolution rule -- a gene name, a
gene-name prefix (for the allele-suffixed ones), a product regex, or a literal
locus tag -- and **a rule matching zero or more than one CDS is a hard error**.
Silence is the failure mode being defended against; a target set that quietly
comes back one locus short still trains, still scores, and still looks fine.

Emits two things:

  **all-gene clumps** in the MTB `CLUMP_COLUMNS` format, so every downstream
  script (precompute, enrichment, coverage audit) works unchanged.
  **per-drug target lists** of `<locus_tag>:gene` and `<locus_tag>:upstream_300bp`,
  the same shape `--fixed-targets` consumes.

The curated mechanism table below is chromosomal only. Acquired genes
(carbapenemases, ESBLs, aminoglycoside-modifying enzymes) are not in PAO1 and
have no coordinates here; they need the AMRFinder route and are deliberately out
of scope for this script.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys

UPSTREAM = 300

CLUMP_COLUMNS = ["target_id", "target_type", "reference_contig", "target_start",
                 "target_end", "reference_strand", "reference_gene_id",
                 "reference_gene_symbol", "reference_product_name",
                 "fallback_anchor_pos", "anchors_all", "anchors_genic",
                 "anchors_upstream200"]

# symbol -> how to find it. ("gene", exact) | ("gene_prefix", p) | ("product", regex)
# | ("locus", tag).  Anything that does not resolve to exactly one CDS is fatal.
RESOLVE: dict[str, tuple[str, str]] = {
    # fluoroquinolone targets
    "gyrA": ("gene", "gyrA"), "gyrB": ("gene", "gyrB"),
    "parC": ("gene", "parC"),
    # the other gene=parE is a 281 bp toxin-antitoxin ParE, not the topoisomerase
    "parE": ("product", r"DNA topoisomerase IV subunit B"),
    # carbapenem porin -- three paralogs share the name
    "oprD": ("product", r"^outer membrane porin OprD$"),
    # AmpC cephalosporinase and its regulators
    "ampC": ("product", r"class C beta-lactamase"),
    "ampR": ("gene", "ampR"),
    "dacB": ("gene", "dacB"),
    "mpl": ("gene", "mpl"),
    "ampD": ("gene", "ampD"),
    # MexAB-OprM and regulators
    "mexA": ("gene", "mexA"), "mexB": ("gene", "mexB"),
    # two genes carry gene=oprM; the operon one sits immediately 3' of mexB
    "oprM": ("product", r"outer membrane channel subunit OprM"),
    "mexR": ("gene", "mexR"), "nalC": ("gene", "nalC"), "nalD": ("gene", "nalD"),
    # MexXY and its repressor (allele-suffixed in this annotation)
    "mexX": ("gene", "mexX"), "mexY": ("gene", "mexY"),
    "mexZ": ("gene_prefix", "mexZ"),
    # MexCD-OprJ / MexEF-OprN and regulators
    "mexC": ("gene", "mexC"), "mexD": ("gene", "mexD"), "oprJ": ("gene", "oprJ"),
    "nfxB": ("gene", "nfxB"),
    "mexE": ("gene", "mexE"), "mexF": ("gene", "mexF"), "oprN": ("gene", "oprN"),
    "mexT": ("gene", "mexT"), "mexS": ("gene", "mexS"),
    # aminoglycoside ribosome / two-component systems
    "fusA1": ("gene", "fusA1"),
    # parRS induces mexXY and is a genuine aminoglycoside mechanism. pmrAB and
    # phoPQ are LPS-modification systems for polymyxins, which are not in this
    # panel -- they were in an earlier draft of this list by reflex and are out.
    "parR": ("gene", "parR"), "parS": ("gene", "parS"),
}

FQ = ["gyrA", "gyrB", "parC", "parE",
      "mexA", "mexB", "oprM", "mexR", "nalC", "nalD",
      "mexC", "mexD", "oprJ", "nfxB", "mexE", "mexF", "oprN", "mexT", "mexS"]
AMPC = ["ampC", "ampR", "dacB", "mpl", "ampD"]
MEXAB = ["mexA", "mexB", "oprM", "mexR", "nalC", "nalD"]
MEXXY = ["mexX", "mexY", "mexZ"]
AG = MEXXY + ["fusA1", "parR", "parS"]

DRUGS: dict[str, list[str]] = {
    "ciprofloxacin": FQ,
    "levofloxacin": FQ,
    "imipenem": ["oprD"] + AMPC + MEXAB,
    "meropenem": ["oprD"] + AMPC + MEXAB + ["mexC", "mexD", "oprJ", "nfxB"],
    "ceftazidime": AMPC + MEXAB,
    "cefepime": AMPC + MEXAB + MEXXY,
    "piperacillin-tazobactam": AMPC + MEXAB,
    "aztreonam": AMPC + MEXAB,
    "amikacin": AG,
    "tobramycin": AG,
}


def load_cds(gff: pathlib.Path) -> list[dict]:
    out = []
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            out.append({"contig": f[0], "start": int(f[3]), "end": int(f[4]),
                        "strand": f[6], "lt": a.get("locus_tag", ""),
                        "gene": a.get("gene", ""),
                        "product": a.get("product", "").replace("%2C", ",")})
    return out


def resolve(symbol: str, cds: list[dict]) -> dict:
    how, key = RESOLVE[symbol]
    if how == "gene":
        hits = [c for c in cds if c["gene"] == key]
    elif how == "gene_prefix":
        hits = [c for c in cds if c["gene"] == key or c["gene"].startswith(key + "_")]
    elif how == "product":
        rx = re.compile(key, re.I)
        hits = [c for c in cds if rx.search(c["product"])]
    else:
        hits = [c for c in cds if c["lt"] == key]
    if len(hits) != 1:
        raise SystemExit(
            f"[ERROR] {symbol!r} ({how}={key!r}) resolved to {len(hits)} CDS, "
            f"need exactly 1.\n" +
            "".join(f"         {h['lt']} gene={h['gene']!r} {h['product'][:70]}\n"
                    for h in hits[:6]) +
            "         A locus that resolves to zero is a locus the model never "
            "sees; one that\n         resolves to many is the wrong gene. Fix the "
            "rule in RESOLVE.")
    return hits[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="PAO1 Bakta GFF3")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--clumps-out", default="",
                    help="all-gene clumps TSV; defaults to <out-dir>/all_gene_target_clumps.tsv")
    args = ap.parse_args()

    cds = load_cds(pathlib.Path(args.gff))
    named = sum(1 for c in cds if c["gene"])
    print(f"[INFO] {len(cds):,} CDS, {named:,} with a gene= attribute "
          f"({named/len(cds):.0%})")

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ---- all-gene clumps, MTB format -------------------------------------
    clumps = []
    for c in cds:
        up_s, up_e = ((c["end"] + 1, c["end"] + UPSTREAM) if c["strand"] == "-"
                      else (max(1, c["start"] - UPSTREAM), c["start"] - 1))
        for kind, s, e in (("gene", c["start"], c["end"]),
                           ("upstream_300bp", up_s, up_e)):
            clumps.append({"target_id": f"{c['lt']}:{kind}", "target_type": kind,
                           "reference_contig": c["contig"], "target_start": s,
                           "target_end": e, "reference_strand": c["strand"],
                           "reference_gene_id": c["lt"],
                           "reference_gene_symbol": c["gene"],
                           "reference_product_name": c["product"],
                           "fallback_anchor_pos": s, "anchors_all": "",
                           "anchors_genic": "", "anchors_upstream200": ""})
    cpath = pathlib.Path(args.clumps_out or out / "all_gene_target_clumps.tsv")
    with cpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CLUMP_COLUMNS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(clumps)
    print(f"[OK] {cpath}  ({len(clumps):,} targets)")

    # ---- resolve every curated symbol before writing anything ------------
    loci = {s: resolve(s, cds) for s in sorted(RESOLVE)}
    print(f"[INFO] all {len(loci)} curated symbols resolved to exactly one CDS")
    for s in ("ampC", "mexZ", "oprD"):
        c = loci[s]
        print(f"       {s:<6} -> {c['lt']}  gene={c['gene'] or '(none)'!r:<16}"
              f"{c['start']}-{c['end']}{c['strand']}")

    # ---- per-drug target lists -------------------------------------------
    print(f"\n{'drug':<26}{'loci':>6}{'targets':>9}")
    tdir = out / "panel_targets"
    tdir.mkdir(exist_ok=True)
    for drug, syms in DRUGS.items():
        seen, lines = set(), []
        for s in syms:
            if s in seen:
                continue
            seen.add(s)
            lt = loci[s]["lt"]
            lines += [f"{lt}:gene", f"{lt}:upstream_300bp"]
        p = tdir / f"targets_{drug}_curated.txt"
        p.write_text("\n".join(lines) + "\n")
        print(f"{drug:<26}{len(seen):>6}{len(lines):>9}")
    print(f"\n[OK] {tdir}")
    print("\nChromosomal loci only. Acquired carbapenemases, ESBLs and AMEs are "
          "absent from\nPAO1 and need the AMRFinder route -- a target set built "
          "from this alone will\nmiss them, which for the beta-lactams is most of "
          "the acquired mechanism.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
