#!/usr/bin/env python3
"""Per-drug target lists for the whole panel, and the union to extract once.

The isoniazid target set took several rounds to settle, and every step of it
generalises: read the curated table's rows for the drug, resolve gene symbols and
intergenic labels to locus tags, drop anything the catalogue annotates for another
drug, and keep the unannotated regions that flank a kept gene because that is
where promoters live. Doing that by hand twelve more times would reproduce the
reasoning and, more likely, reproduce the mistakes.

It emits **candidate sets to test, not one answer**, because no single rule is
right for every drug. Isoniazid's winning set came from the curated table and
adding the catalogue's ten isoniazid genes measurably hurt; amikacin's curated row
resolves to one locus and its real set -- rrs, eis, whiB7 -- is in the catalogue.
Seven candidate sets decided isoniazid; the per-drug development CV decides the
rest.

  **curatedmech**  the curated table's entries that the catalogue annotates for
                   *this* drug, plus unannotated loci rescued two ways: lying
                   between the flanking genes of an intergenic entry (oxyR' is the
                   ahpC promoter and no entry names it), or within 500 bp of a
                   kept locus (fabG1 is 18 bp from inhA, is its own entry, and
                   tbdb.bed does not contain it under any name -- a target set
                   without its upstream window misses 23% of resistant isolates).
  **catgenes**     the catalogue's own genes for the drug, annotated for it and
                   for nothing outside its mechanistic relatives.
  **both**         the union.
  **novel**        unannotated, unrescued. Reported, not embedded into any model
                   by default: isoniazid's six earned nothing (bal_OT_common
                   -0.003 to +0.002, every comparison non-significant).
  **other**        annotated for a different drug. Dropped -- a GWAS for one drug
                   on an MDR-heavy cohort finds the others' markers (section 1),
                   which is why the curated isoniazid rows contain rpoB, rpsL, gid
                   and embB.

Membership tests on `dn in ann`, never `ann & relatives`. Relatives exist so inhA
is not dropped for carrying an ethionamide annotation beside isoniazid; a locus
annotated *only* for the relative -- ethA, ethionamide's prodrug activator -- is
not an isoniazid gene.

Relatives are shared *mechanism*, not co-occurrence: isoniazid and ethionamide
both act through InhA, amikacin and kanamycin both through rrs, the two
fluoroquinolones through gyrA/gyrB, bedaquiline and clofazimine through mmpR5.
Dropping a locus because it is annotated for the drug that shares its target
protein would be wrong.

Symbol resolution uses the all-gene clumps file rather than a hand-written map,
so every symbol the extraction knows about resolves, and anything that does not
is reported rather than silently dropped.

The union of every drug's mech and novel targets is what gets extracted and
embedded -- once, since 13l's cost is per-sample consensus regeneration and does
not scale with target count. Each drug then subsets it with --fixed-targets.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys

# Shared drug *target*, not epidemiological co-occurrence. Being generous here is
# safe: membership requires `dn in ann` first, so a relative can only rescue a
# locus already annotated for the drug itself. Adding capreomycin to linezolid
# lets rrl (capreomycin,linezolid) through and still cannot let tlyA
# (capreomycin-only) in.
_AMINOGLYCOSIDE = {"amikacin", "kanamycin", "capreomycin", "streptomycin"}
_QUINOLONE = {"levofloxacin", "moxifloxacin", "ciprofloxacin", "ofloxacin"}
RELATIVES = {
    # both are prodrugs inhibiting InhA; only the activator differs (katG / ethA)
    "isoniazid": {"ethionamide"},
    "ethionamide": {"isoniazid"},
    # all four act through the 30S rRNA, rrs
    **{d: _AMINOGLYCOSIDE - {d} for d in _AMINOGLYCOSIDE},
    # gyrA / gyrB
    **{d: _QUINOLONE - {d} for d in _QUINOLONE},
    # mmpR5 / mmpL5
    "bedaquiline": {"clofazimine"},
    "clofazimine": {"bedaquiline"},
    # rpoB
    "rifampicin": {"rifabutin"},
    "rifabutin": {"rifampicin"},
    # the 50S rRNA, rrl
    "linezolid": {"capreomycin"},
    # Both nitroimidazoles are activated by Ddn through F420, so ddn and fgd1 /
    # fbiA-D are delamanid mechanism. Asymmetric on purpose: tbdb also tags those
    # genes for clofazimine, where the link is indirect -- F420 affects the redox
    # cycling clofazimine depends on -- so for clofazimine they stay in the
    # multidrug bucket to be tested rather than assumed. The table answers "which
    # other drugs' annotations should not disqualify a locus *for this drug*",
    # and that question genuinely has different answers in the two directions.
    "delamanid": {"pretomanid", "clofazimine"},
    "pretomanid": {"delamanid"},
}
RELATIVES["capreomycin"] = RELATIVES["capreomycin"] | {"linezolid"}


def norm_drug(s: str) -> str:
    return re.sub(r"[^a-z]", "", s.lower())


def load_symbol_map(clumps: pathlib.Path) -> tuple[dict[str, str], set[str]]:
    """symbol -> locus tag, and the set of locus tags that have targets."""
    sym, tags = {}, set()
    with clumps.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = (r.get("reference_gene_id") or "").strip()
            s = (r.get("reference_gene_symbol") or "").strip()
            if gid:
                tags.add(gid)
                if s and s.lower() not in sym:
                    sym[s.lower()] = gid
    return sym, tags


def load_bed_drugs(bed: pathlib.Path) -> dict[str, set[str]]:
    """locus tag -> the drugs tbdb annotates it for."""
    out: dict[str, set[str]] = {}
    for line in bed.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        drugs = {norm_drug(d) for d in re.split(r"[,;]", f[5])
                 if d.strip()} if len(f) > 5 else set()
        out.setdefault(f[3], set()).update(drugs)
        if f[4]:
            out.setdefault(f[4], set()).update(drugs)
    return out


LABEL_CORRECTIONS = {"Rv3430a": "Rv3430c"}   # typo in the curated table


def resolve(token: str, sym: dict[str, str], tags: set[str]) -> str | None:
    t = LABEL_CORRECTIONS.get(token.strip(), token.strip())
    if not t:
        return None
    if t in tags:
        return t
    return sym.get(t.lower())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--curated-entries", required=True, action="append",
                    metavar="DRUG=FILE",
                    help="repeatable; FILE is read_curated_drug_genes.py's --out "
                         "for that drug (entry, is_intergenic)")
    ap.add_argument("--clumps", required=True, help="all_gene_target_clumps.tsv")
    ap.add_argument("--bed", required=True, help="tbdb.bed")
    ap.add_argument("--no-catalogue-genes", action="store_true",
                    help="use only the curated table. Off by default: that table "
                         "is a GWAS table, so for most drugs it lists a handful of "
                         "known genes and a long tail of hits -- amikacin gets one "
                         "mechanism locus from it, and its real set is rrs, eis "
                         "and whiB7, which live in tbdb.bed. Isoniazid is the "
                         "exception where the curated table was the fuller source")
    ap.add_argument("--max-intergenic-span", type=int, default=20,
                    help="loci an intergenic entry may span before it is treated "
                         "as an error. Adjacent flanking genes enclose a handful; "
                         "a span of hundreds means the entry is not intergenic")
    ap.add_argument("--promote-within", type=int, default=500,
                    help="an unannotated locus whose interval lies within this "
                         "many bp of a locus kept for the drug is kept too. That "
                         "is where promoters and operon partners sit: fabG1 is "
                         "18 bp from inhA and carries the -15C>T that half of "
                         "isoniazid resistance runs through, and tbdb.bed does "
                         "not contain it under any name. Sharing a curated "
                         "intergenic entry is not enough -- fabG1 is its own "
                         "entry in the table")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    sym, tags = load_symbol_map(pathlib.Path(args.clumps))
    beddrugs = load_bed_drugs(pathlib.Path(args.bed))
    print(f"[INFO] {len(sym):,} gene symbols, {len(tags):,} locus tags with targets")

    # every target id per locus, in clumps order, plus the locus extent
    per_locus: dict[str, list[str]] = {}
    span: dict[str, tuple[int, int]] = {}
    with open(args.clumps) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = (r.get("reference_gene_id") or "").strip()
            if not gid:
                continue
            per_locus.setdefault(gid, []).append(r["target_id"])
            try:
                s, e = int(r["target_start"]), int(r["target_end"])
            except (KeyError, TypeError, ValueError):
                continue
            lo, hi = span.get(gid, (s, e))
            span[gid] = (min(lo, s), max(hi, e))

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    union: list[str] = []
    seen: set[str] = set()
    summary = []

    for spec in args.curated_entries:
        drug, _, path = spec.partition("=")
        dn = norm_drug(drug)
        rel = {norm_drug(x) for x in RELATIVES.get(dn, set())} | {dn}
        rows = list(csv.DictReader(open(path), delimiter="\t"))

        loci: dict[str, set[str]] = {}     # tag -> the entries that named it
        unresolved = []
        for r in rows:
            e = r["entry"].strip()
            toks = [x for x in re.split(r"[:/]", e) if x.strip()]
            got = []
            for t in toks:
                tag = resolve(t, sym, tags)
                if tag is None:
                    unresolved.append((e, t))
                else:
                    got.append(tag)
                    loci.setdefault(tag, set()).add(e)
            # An intergenic entry names its two flanking genes, so anything lying
            # *between* them is part of the region it refers to and is usually
            # where the mutation is. proA:ahpC names proA and ahpC; oxyR' sits
            # between them, is named by nothing, and is the ahpC promoter.
            if len(got) >= 2 and all(g in span for g in got):
                lo = min(span[g][0] for g in got)
                hi = max(span[g][1] for g in got)
                between = [t for t, (s_, en) in span.items()
                           if t not in got and s_ >= lo and en <= hi]
                # An intergenic entry names two ADJACENT genes, so the region
                # between them holds a handful of loci at most. A span of
                # hundreds means the entry is not intergenic at all -- paralogues
                # written 'cysA2/cysA3' once swept 2,300 loci into rifabutin and
                # were only noticed because the count was absurd. Fail instead.
                if len(between) > args.max_intergenic_span:
                    raise SystemExit(
                        f"[ERROR] {dn}: entry {e!r} spans {len(between):,} loci "
                        f"between {sorted(got)}.\n"
                        f"        An intergenic entry names two adjacent genes; "
                        f"this does not. Check whether the two names are "
                        f"paralogues or alternatives rather than flanks, or raise "
                        f"--max-intergenic-span if the span is genuinely that wide.")
                for t in between:
                    loci.setdefault(t, set()).add(e)
            if not got:
                continue

        # The catalogue's own genes for this drug are a *separate* candidate set,
        # not extra members of the curated one. Isoniazid's winning set came from
        # the curated table and adding the catalogue's ten genes measurably hurt
        # (bal_OT_common 0.693 vs 0.900 without the promoter, and inhcat lost to
        # inhgenes even with it); amikacin's curated row yields one locus and its
        # real set is in the catalogue. No single rule is right for both, so both
        # are emitted and the per-drug CV decides.
        #
        # `dn in ann`, not `ann & rel`: relatives exist so that inhA is not dropped
        # for carrying an ethionamide annotation alongside isoniazid. A locus
        # annotated *only* for the relative -- ethA, ethionamide's prodrug
        # activator -- is not an isoniazid gene and must not be admitted.
        cat_genes = [t for t, ann in beddrugs.items()
                     if t in tags and dn in ann and not (ann - rel)]
        # Annotated for this drug AND for something outside its relatives. Not a
        # decidable class: ndh is isoniazid/delamanid/ethionamide and is a real
        # isoniazid mechanism, while glpK carries six drugs and is a phase-variation
        # locus that turns up for everything. Three annotations against six is not a
        # rule. Reported for a per-drug ruling instead of dropped without trace --
        # it matters most for the drugs where the catalogue, not the curated table,
        # is the main source.
        multidrug = [t for t, ann in beddrugs.items()
                     if t in tags and dn in ann and (ann - rel)]

        mech, other, novel = [], [], []
        for tag in loci:
            ann = beddrugs.get(tag, set())
            if dn in ann:
                mech.append(tag)
            elif ann:
                other.append(tag)
            else:
                novel.append(tag)
        # Two ways an unannotated locus is really part of a kept gene's region:
        # it shares a curated intergenic entry with one (proA/oxyR' for ahpC), or
        # it simply sits next to one (fabG1 is 18 bp from inhA and is its own
        # entry in the table, so the first rule alone would miss it).
        kept_entries = {e for t in mech for e in loci[t]}
        kept_spans = [span[t] for t in mech if t in span]
        def adjacent(t: str) -> bool:
            if t not in span:
                return False
            s, e = span[t]
            return any(s - args.promote_within <= b and a - args.promote_within <= e
                       for a, b in kept_spans)
        promoted = [t for t in novel
                    if (loci[t] & kept_entries) or adjacent(t)]
        novel = [t for t in novel if t not in promoted]
        mech += promoted

        def targets(tl):
            o = []
            for t in sorted(tl):
                o.extend(per_locus.get(t, []))
            bare = set(o)
            return [x for x in o if not (x.startswith("transcript:")
                                         and x.removeprefix("transcript:") in bare)]

        # Candidate sets to test, not one answer. Isoniazid's winner was
        # `curatedmech`; amikacin's curated row is one locus and its answer will
        # be `catgenes` or `both`. The per-drug development CV decides, exactly as
        # seven candidate sets decided isoniazid.
        cands = {
            "curatedmech": mech,
            "catgenes": cat_genes,
            "both": sorted(set(mech) | set(cat_genes)),
            # Everything the source nominates for this drug EXCEPT loci the
            # catalogue annotates as another drug's mechanism. That exclusion is
            # the whole rule: the published GWAS reports rpoB under isoniazid with
            # six signals, because 94.3% of rifampicin-resistant isolates are also
            # isoniazid-resistant, and taking the hits unfiltered imports that.
            #
            # Unannotated loci are kept. The filter cannot tell a genuine novel
            # mechanism from an uncharacterised co-resistance marker, and keeping
            # them is the deliberate choice -- excluding them would discard the
            # only route by which this design can find something the catalogue
            # does not already know.
            #
            # Differs from `both` only by `novel`, and from the union of
            # everything only by `other`.
            "gwasfiltered": sorted(set(mech) | set(cat_genes) | set(novel)),
            "multidrug": multidrug,
            "novel": novel,
            "other": other,
        }
        for name, tl in cands.items():
            p = out / f"targets_{dn}_{name}.txt"
            tt = targets(tl)
            p.write_text("\n".join(tt) + "\n")
            if name != "other":       # everything we might model on gets embedded
                for x in tt:
                    if x not in seen:
                        seen.add(x)
                        union.append(x)
        summary.append({
            "drug": dn,
            "curatedmech": len(mech), "promoted": len(promoted),
            "catgenes": len(cat_genes), "multidrug": len(multidrug),
            "both": len(cands["both"]), "novel": len(novel), "other": len(other),
            "other_loci": sorted(other), "novel_loci": sorted(novel),
            "multidrug_loci": sorted(multidrug), "unresolved": unresolved,
        })

    (out / "targets_panel_union.txt").write_text("\n".join(union) + "\n")

    cols = ["curatedmech", "promoted", "catgenes", "multidrug",
            "both", "novel", "other"]
    print(f"\n{'drug':<16}" + "".join(f"{c:>12}" for c in cols))
    for r in summary:
        print(f"{r['drug']:<16}" + "".join(f"{r[c]:>12}" for c in cols))
    print(f"\nunion for extraction: {len(union)} targets, "
          f"{len({t.split(':')[-2] if len(t.split(':')) >= 2 else t for t in union})} loci")
    print(f"[OK] {out}/targets_panel_union.txt")

    for key, title in (("other_loci", "dropped as another drug's gene"),
                       ("multidrug_loci", "annotated for this drug AND others "
                                          "— needs a per-drug ruling (ndh vs glpK)")):
        print(f"\n{title}:")
        for r in summary:
            if r[key]:
                print(f"  {r['drug']:<14}{' '.join(r[key])}")

    bad = [(r["drug"], u) for r in summary for u in r["unresolved"]]
    if bad:
        print(f"\nUNRESOLVED symbols ({len(bad)}) — silently missing targets:")
        for dn, (entry, tok) in bad[:40]:
            print(f"  {dn:<14}{entry:<30}{tok}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
