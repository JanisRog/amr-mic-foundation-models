#!/usr/bin/env python3
"""Acquired-gene presence from AMRFinderPlus, per isolate and per gene family.

The curated PAO1 target set covers chromosomal mechanisms and nothing else.
Carbapenemases, ESBLs and aminoglycoside-modifying enzymes are not in PAO1, have
no coordinates to extract from, and for several beta-lactams are most of the real
resistance -- which is why ARDaP is at chance on cefepime (AUROC 0.552) and why a
model built only on chromosomal loci would be expected to fail there too. This
builds the missing half as presence/absence, which is cheap, and answers whether
those genes carry signal at all before anything is embedded.

Two rules this script deliberately does **not** apply, both raised and both wrong:

  **No quality filter.** AMRFinder applies curated per-element thresholds -- HMM
  cutoffs and BLAST criteria tuned per gene -- and reports only hits that pass.
  Imposing a further `% Coverage` / `% Identity` cut would substitute our
  judgement for the tool's, which is the mistake that cost ethionamide 0.077 of
  development cost (MTB protocol section 11). If AMRFinder reports it, it is
  present. The `Method` column (EXACTX / ALLELEX / BLASTX / PARTIALX ...) is
  carried through so partial hits stay *visible* rather than being filtered:
  a truncation at a contig break is an artefact, a truncated porin is a
  mechanism, and a coverage threshold cannot tell them apart while a model might.

  **A carrier floor, but derived rather than inherited.** The first version of
  this script had none, on the argument that a random forest sampling
  `max_features` ignores a rare column cheaply. The carriage distribution refutes
  that at this scale: 296 of 383 families are carried by <=1% of isolates, so the
  rare tail is 77% of the feature width, not a handful, and it dilutes every
  split. The floor is set from the fold structure instead of from precedent --
  with 5-fold CV a gene with N carriers appears in a validation fold about 0.2N
  times, so N >= 10 is where an effect can be validated at all. MTB's >=20 on
  10,170 isolates was 0.2% of that cohort; copying the absolute number here would
  be 1%, five times stricter on the smaller dataset.

Families come from `Hierarchy node` for acquired genes (`Subtype=AMR`), so
blaVIM-1 and blaVIM-2 share a column and the allele distinction is left to the
sequence embedding. **Point mutations are kept at allele level** (`Element
symbol`, `Subtype=POINT`): collapsing gyrA_D87G with gyrA_D72G, or oprD_E140STOP
with oprD_E176STOP, would destroy exactly the curated knowledge the embeddings do
not have. Those columns also overlap the chromosomal targets we embed, which is
double-counting on purpose -- AMRFinder's call means "this matches a known
resistance allele", which is a different statement from the sequence itself, and
whether it adds anything is testable by ablation rather than assertion.

Drugs come from AMRFinder's own `Class`/`Subclass`, not from a mapping invented
here. Non-antibiotic resistance is excluded by `Subtype` -- METAL, BIOCIDE and
HEAT -- which is cleaner than blacklisting classes by hand: mercury and chromate
determinants travel on the same plasmids as real resistance genes, so including
them would be pure confounding with plasmid carriage.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import sys

# AMRFinder Class -> panel drugs it plausibly affects. Within-class breadth is
# mechanistically defensible (a beta-lactamase acts on beta-lactams); across
# class it would be the co-resistance route MTB section 1 warns about.
CLASS_DRUGS = {
    "BETA-LACTAM": ["imipenem", "meropenem", "ceftazidime", "cefepime",
                    "piperacillin-tazobactam", "aztreonam"],
    "CARBAPENEM": ["imipenem", "meropenem"],
    "CEPHALOSPORIN": ["ceftazidime", "cefepime"],
    "MONOBACTAM": ["aztreonam"],
    "AMINOGLYCOSIDE": ["amikacin", "tobramycin"],
    "AMIKACIN": ["amikacin"],
    "TOBRAMYCIN": ["tobramycin"],
    "GENTAMICIN": ["amikacin", "tobramycin"],
    "KANAMYCIN": ["amikacin", "tobramycin"],
    "QUINOLONE": ["ciprofloxacin", "levofloxacin"],
    "FLUOROQUINOLONE": ["ciprofloxacin", "levofloxacin"],
    "EFFLUX": ["ciprofloxacin", "levofloxacin", "imipenem", "meropenem",
               "ceftazidime", "cefepime", "piperacillin-tazobactam",
               "aztreonam", "amikacin", "tobramycin"],
}
# Subtype is a cleaner filter than a class blacklist: METAL, BIOCIDE and HEAT are
# whole categories of non-antibiotic resistance, and enumerating their classes by
# hand was catching them only incidentally.
KEEP_SUBTYPES = {"AMR", "POINT"}

# A truncating change is collapsed to one <gene>_LOF column. Allele level is right
# for a substitution -- gyrA_T83I (997 carriers) is the canonical fluoroquinolone
# mutation and must not be blurred with rarer gyrA variants -- but wrong for a
# knockout, where every premature stop has the identical consequence. oprD is the
# case that forced this: its eight surviving truncation columns carry 10-83
# isolates each and are individually unlearnable, while together they are ~488,
# and the one oprD column that did clear the floor was V359L (860 carriers), which
# ARDaP annotates "Natural variation, does not confer AMR" -- a lineage proxy.
# Left alone the feature set carried the oprD confound and dropped the oprD
# mechanism. The same logic applies to every LOF repressor: mexZ, nfxB, mexR.
LOF_TOKENS = ("STOP", "FS", "DEL", "INS", "DUP")

# `Hierarchy node` turns out to be the *allele* (blaVIM-1, blaVIM-2, blaVIM-20),
# not the family, so acquired genes arrive split. For carbapenemases that is the
# wrong granularity for the same reason oprD truncations were: every VIM
# hydrolyses carbapenems, and one column of 88 carriers is learnable where three
# of 13-49 are marginal. Collapsed by explicit list rather than by stripping the
# trailing number, because OXA numbering spans both the intrinsic P. aeruginosa
# oxacillinase (blaOXA-50 and relatives, carried by nearly every isolate) and
# genuinely acquired enzymes (blaOXA-2, blaOXA-10) -- merging those would fuse a
# strain marker with a resistance gene.
COLLAPSE_FAMILIES = ("blaVIM", "blaIMP", "blaNDM", "blaGES", "blaKPC", "blaSPM",
                     "blaSIM", "blaAIM", "blaDIM", "blaGIM", "blaCTX-M",
                     "blaPER", "blaVEB")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--amrfinder-dir", required=True,
                    help="directory of <sample>/<sample>.amrfinderplus.tsv")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-carriers", type=int, default=10,
                    help="a column needs enough carriers to appear in a validation "
                         "fold: at 5-fold CV that is about 0.2N, so 10 is the point "
                         "an effect can be validated. Set 0 to keep everything")
    ap.add_argument("--max-carrier-frac", type=float, default=0.99,
                    help="drop near-universal columns, which cannot discriminate. "
                         "This removes no information -- a gene every isolate has "
                         "splits nothing")
    args = ap.parse_args()

    files = sorted(pathlib.Path(args.amrfinder_dir).rglob("*.amrfinderplus.tsv"))
    if not files:
        raise SystemExit(f"[ERROR] no AMRFinder TSVs under {args.amrfinder_dir}")
    print(f"[INFO] {len(files):,} isolates")

    present: dict[str, set[str]] = collections.defaultdict(set)
    meta: dict[str, dict] = {}
    methods: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    dropped = collections.Counter()
    samples = []

    for i, f in enumerate(files, 1):
        sample = f.name[: -len(".amrfinderplus.tsv")]
        samples.append(sample)
        with f.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sub_t = (r.get("Subtype") or "").strip().upper()
                if sub_t not in KEEP_SUBTYPES:
                    dropped[f"subtype {sub_t or '(blank)'}"] += 1
                    continue
                cls = (r.get("Class") or "").strip().upper()
                if not cls or cls == "NA":
                    dropped["no class"] += 1
                    continue
                sym = (r.get("Element symbol") or "").strip()
                if sub_t == "POINT":
                    gene, _, change = sym.partition("_")
                    up = change.upper()
                    fam = (f"{gene}_LOF"
                           if gene and any(k in up for k in LOF_TOKENS)
                           else sym)
                else:
                    fam = (r.get("Hierarchy node") or "").strip()
                    if not fam or fam == "NA":
                        fam = sym
                    for base in COLLAPSE_FAMILIES:
                        if fam.startswith(base + "-") or fam == base:
                            fam = base
                            break
                if not fam:
                    dropped["no family"] += 1
                    continue
                present[fam].add(sample)
                methods[fam][(r.get("Method") or "").strip()] += 1
                m = meta.setdefault(fam, {"class": cls,
                                          "subclass": (r.get("Subclass") or "").strip().upper(),
                                          "subtype": sub_t, "symbols": set()})
                m["symbols"].add(sym)
        if i % 500 == 0:
            print(f"[INFO] {i:,}/{len(files):,} read")

    for k, n in dropped.most_common():
        print(f"[INFO] skipped {n:,} hits: {k}")

    n_all = len(present)
    too_few = sum(1 for f in present if len(present[f]) < args.min_carriers)
    too_many = sum(1 for f in present
                   if len(present[f]) / len(samples) > args.max_carrier_frac)
    fams = sorted(f for f in present
                  if len(present[f]) >= args.min_carriers
                  and len(present[f]) / len(samples) <= args.max_carrier_frac)
    print(f"[INFO] {n_all:,} families -> {len(fams):,} kept "
          f"({too_few:,} under {args.min_carriers} carriers, "
          f"{too_many:,} over {args.max_carrier_frac:.0%} carriage)")
    n_point = sum(1 for f in fams if meta[f].get("subtype") == "POINT")
    n_lof = sum(1 for f in fams if f.endswith("_LOF"))
    coll = [f for f in fams if f in COLLAPSE_FAMILIES]
    if coll:
        print("[INFO] collapsed carbapenemase/ESBL families: " +
              ", ".join(f"{f}={len(present[f])}" for f in sorted(coll)))
    print(f"[INFO] {n_point} point-mutation columns "
          f"({n_lof} of them collapsed loss-of-function), "
          f"{len(fams) - n_point} acquired-gene columns (family level)")

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    with (out / "amrfinder_presence.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id"] + [f"amr_{f}" for f in fams])
        for s in samples:
            w.writerow([s] + [1 if s in present[f] else 0 for f in fams])

    with (out / "amrfinder_columns.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["family", "subtype", "class", "subclass",
                                       "carriers", "carrier_frac", "drugs",
                                       "methods", "symbols"])
        w.writeheader()
        for f in fams:
            m = meta[f]
            drugs = CLASS_DRUGS.get(m["subclass"]) or CLASS_DRUGS.get(m["class"]) or []
            w.writerow({"family": f, "subtype": m.get("subtype", ""),
                        "class": m["class"], "subclass": m["subclass"],
                        "carriers": len(present[f]),
                        "carrier_frac": round(len(present[f]) / len(samples), 4),
                        "drugs": ",".join(drugs),
                        "methods": ",".join(f"{k}:{v}" for k, v in
                                            methods[f].most_common(3)),
                        "symbols": ",".join(sorted(m["symbols"])[:6])})

    # per-drug column lists, ready for --extra-feature-tsv style use
    bydrug: dict[str, list[str]] = collections.defaultdict(list)
    for f in fams:
        m = meta[f]
        for d in (CLASS_DRUGS.get(m["subclass"]) or CLASS_DRUGS.get(m["class"]) or []):
            bydrug[d].append(f"amr_{f}")
    cdir = out / "amr_columns_by_drug"
    cdir.mkdir(exist_ok=True)
    print(f"\n{'drug':<26}{'columns':>8}{'>=20 carriers':>15}")
    for d in sorted(bydrug):
        (cdir / f"{d}.txt").write_text("\n".join(sorted(bydrug[d])) + "\n")
        big = sum(1 for c in bydrug[d] if len(present[c[4:]]) >= 20)
        print(f"{d:<26}{len(bydrug[d]):>8}{big:>15}")

    print(f"\n{'family':<22}{'class':<16}{'carriers':>9}  symbols")
    for f in sorted(fams, key=lambda f: -len(present[f]))[:15]:
        print(f"{f:<22}{meta[f]['class'][:15]:<16}{len(present[f]):>9}  "
              f"{','.join(sorted(meta[f]['symbols'])[:3])}")
    print(f"\n[OK] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
