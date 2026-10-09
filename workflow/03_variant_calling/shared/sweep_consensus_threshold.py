#!/usr/bin/env python3
"""Which consensus threshold represents the measured MIC best?

The pipeline keeps Snippy. The question is where to cut. Today the cut is
freebayes' ploidy-2 genotype call, which behaves as a step at allele fraction 0.9
(section 7) and discards 953 catalogued calls whose phenotype says they were real.
Lowering it recovers those; lowering it too far writes mapping artefacts into the
sequence every model embeds. The measured MIC is the arbiter.

Nothing is re-called. `snps.raw.vcf` is freebayes' pre-filter output and sits in
every run directory with AO and DP at every site, so a threshold sweep is a parse.
Whichever threshold wins can then be applied by re-filtering that same file and
rebuilding the consensus with bcftools -- still no realignment.

For each candidate threshold, an isolate is called resistant to a drug if it
carries any catalogued variant for that drug at allele fraction >= t, subject to
the depth and quality bars Snippy already used. That call is scored against
MIC > boundary.

The target is vMER -- a resistant isolate called susceptible is the error that
matters. But minimising vMER *alone* is degenerate: calling more variants can only
turn false-susceptible into true-resistant, so the lowest threshold on the grid
always wins and the sweep would recommend 0.05 whatever the data said. The
criterion is therefore the lowest vMER among thresholds that do not cost
specificity against what the pipeline does today, with --mer-tolerance to price an
explicit relaxation. The unconstrained minimum is printed too, labelled as the
floor it is.

One threshold, not one per drug: there is a single consensus sequence per isolate,
so a per-drug optimum would be a number no pipeline can implement.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import random
import sys

DRUGS_WITH_BOUNDARY = {"RIF": -1.0, "INH": -3.0, "EMB": 2.0, "ETH": 2.0,
                       "KAN": 2.0, "LEV": 0.0, "MXF": 0.0, "BDQ": -2.0,
                       "DLM": -3.0, "AMI": 0.0, "CFZ": 2.0, "LZD": 0.0}


def norm(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    i = 0
    while i < len(ref) - 1 and i < len(alt) - 1 and ref[i] == alt[i]:
        i += 1
    return pos + i, ref[i:], alt[i:]


def atomize(pos: int, ref: str, alt: str) -> list[tuple[int, str, str]]:
    """Equal-length haplotype calls split into their substitutions (section 7)."""
    if len(ref) != len(alt) or len(ref) == 1:
        return []
    return [(pos + i, r, a) for i, (r, a) in enumerate(zip(ref, alt)) if r != a]


def load_catalogue(path: pathlib.Path) -> dict[tuple[int, str, str], set[str]]:
    """(pos, ref, alt) -> drugs, from the per-variant table already produced."""
    cat: dict[tuple[int, str, str], set[str]] = {}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                pos = int(r["genome_pos"])
            except (KeyError, TypeError, ValueError):
                continue
            ref, alt = r["ref"].upper(), r["alt"].upper()
            drugs = {d for d in r["drugs"].split(";") if d}
            if not (pos and ref and alt and drugs):
                continue
            for key in [norm(pos, ref, alt), *atomize(pos, ref, alt)]:
                cat.setdefault(key, set()).update(drugs)
    return cat


def isolate_fractions(path: pathlib.Path, wanted: set[int], min_qual: float,
                      min_dp: int) -> dict[tuple[int, str, str], float]:
    """Allele fraction at the catalogue positions this isolate has data for."""
    out: dict[tuple[int, str, str], float] = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 10)
            if len(f) < 10:
                continue
            try:
                pos = int(f[1])
            except ValueError:
                continue
            if pos not in wanted:
                continue
            try:
                if float(f[5]) < min_qual:
                    continue
            except ValueError:
                continue
            fmt, vals = f[8].split(":"), f[9].rstrip("\n").split(":")
            d = dict(zip(fmt, vals))
            try:
                dp = int(d.get("DP", 0))
            except ValueError:
                continue
            if dp < min_dp:
                continue
            alts = f[4].upper().split(",")
            aos = [x for x in d.get("AO", "").split(",") if x not in ("", ".")]
            if len(aos) != len(alts):
                continue
            ref = f[3].upper()
            for alt, ao in zip(alts, aos):
                try:
                    frac = int(ao) / dp
                except (ValueError, ZeroDivisionError):
                    continue
                for key in [norm(pos, ref, alt), *atomize(pos, ref, alt)]:
                    out[key] = max(out.get(key, 0.0), frac)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snippy-runs", required=True)
    ap.add_argument("--catalogue", required=True,
                    help="call_sources_by_variant.tsv from compare_call_sources_vs_mic.py")
    ap.add_argument("--mic-tsv", required=True)
    ap.add_argument("--mic-suffix", default="_int")
    ap.add_argument("--thresholds", default="0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95")
    ap.add_argument("--min-qual", type=float, default=100.0)
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--vmer-weight", type=float, default=2.0,
                    help="how many major errors one very major error is worth; the "
                         "recommendation minimises w*vMER + MER")
    ap.add_argument("--reference-threshold", type=float, default=0.9,
                    help="what the pipeline does today; the recommendation is not "
                         "allowed to be worse than this on MER")
    ap.add_argument("--mer-tolerance", type=float, default=0.0,
                    help="how much mean MER may rise above the reference in "
                         "exchange for lower vMER; 0 means none")
    ap.add_argument("--sample-n", type=int, default=0,
                    help="isolates to read; 0 reads all of them")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default="", help="per-threshold per-drug metrics")
    args = ap.parse_args()

    thresholds = [float(t) for t in args.thresholds.split(",")]
    cat = load_catalogue(pathlib.Path(args.catalogue))
    positions = {k[0] for k in cat}
    drugs_seen = sorted({d for v in cat.values() for d in v} & set(DRUGS_WITH_BOUNDARY))
    print(f"[INFO] catalogue: {len(cat):,} keys at {len(positions):,} positions, "
          f"{len(drugs_seen)} scorable drugs")

    mic: dict[str, dict[str, float]] = {}
    with open(args.mic_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            vals = {}
            for d in DRUGS_WITH_BOUNDARY:
                raw = r.get(d + args.mic_suffix, "")
                if raw not in ("", "NA"):
                    vals[d] = float(raw)
            if vals:
                mic[r["sample_id"]] = vals
    print(f"[INFO] MIC for {len(mic):,} samples")

    runs = pathlib.Path(args.snippy_runs)
    dirs = [d for d in sorted(runs.iterdir())
            if d.name in mic and (d / "snps.raw.vcf").exists()]
    print(f"[INFO] {len(dirs):,} isolates with a raw VCF and a MIC")
    if args.sample_n and args.sample_n < len(dirs):
        dirs = random.Random(args.seed).sample(dirs, args.sample_n)
        print(f"[INFO] reading a random {len(dirs):,}")
    if not dirs:
        raise SystemExit("[ERROR] no isolates to read")

    # per drug and threshold: tp, fp, fn, tn
    counts = {(d, t): [0, 0, 0, 0] for d in drugs_seen for t in thresholds}
    for i, d in enumerate(dirs, 1):
        if i % 500 == 0:
            print(f"[INFO] {i:,}/{len(dirs):,}", flush=True)
        fracs = isolate_fractions(d / "snps.raw.vcf", positions,
                                  args.min_qual, args.min_dp)
        best: dict[str, float] = collections.defaultdict(float)
        for key, frac in fracs.items():
            for drug in cat.get(key, ()):
                if frac > best[drug]:
                    best[drug] = frac
        for drug in drugs_seen:
            if drug not in mic[d.name]:
                continue
            true_r = mic[d.name][drug] > DRUGS_WITH_BOUNDARY[drug]
            top = best.get(drug, 0.0)
            for t in thresholds:
                pred_r = top >= t
                c = counts[(drug, t)]
                c[0 if (true_r and pred_r) else 1 if (pred_r and not true_r)
                  else 2 if true_r else 3] += 1

    rows = []
    print(f"\n{'drug':<5}{'thr':>6}{'sens':>8}{'spec':>8}{'bal':>8}"
          f"{'vMER':>8}{'MER':>8}{'F1':>8}{'n R':>7}")
    for drug in drugs_seen:
        for t in thresholds:
            tp, fp, fn, tn = counts[(drug, t)]
            if tp + fn == 0 or tn + fp == 0:
                continue
            sens, spec = tp / (tp + fn), tn / (tn + fp)
            row = {"drug": drug, "threshold": t, "sens": sens, "spec": spec,
                   "bal": (sens + spec) / 2, "vMER": fn / (tp + fn),
                   "MER": fp / (tn + fp),
                   "F1": 2 * tp / max(2 * tp + fp + fn, 1),
                   "n_resistant": tp + fn}
            rows.append(row)
            print(f"{drug:<5}{t:>6.2f}{sens:>8.3f}{spec:>8.3f}{row['bal']:>8.3f}"
                  f"{row['vMER']:>8.3f}{row['MER']:>8.3f}{row['F1']:>8.3f}"
                  f"{tp+fn:>7}")
        print()

    # One global threshold: there is one consensus per isolate, so the choice
    # cannot be per drug. Drugs are weighted equally rather than by prevalence,
    # or isoniazid and rifampicin would decide it alone.
    #
    # vMER is the error to minimise, but minimising it alone is degenerate: the
    # lowest threshold on the grid always wins, because calling more variants can
    # only reduce false-susceptible. So the criterion is the lowest vMER among
    # thresholds that do not cost specificity against what the pipeline does today.
    per_t = {}
    for t in thresholds:
        sel = [r for r in rows if r["threshold"] == t]
        if sel:
            per_t[t] = (sum(r["bal"] for r in sel) / len(sel),
                        sum(r["vMER"] for r in sel) / len(sel),
                        sum(r["MER"] for r in sel) / len(sel))
    if not per_t:
        raise SystemExit("[ERROR] no threshold produced scorable drugs")

    ref = min(per_t, key=lambda t: abs(t - args.reference_threshold))
    ref_mer = per_t[ref][2]
    budget = ref_mer + args.mer_tolerance
    print(f"Pooled across {len(drugs_seen)} drugs, each weighted equally.")
    print(f"Reference is threshold {ref:g} (what the pipeline does today): "
          f"mean MER {ref_mer:.4f}.")
    print(f"Admissible = mean MER <= {budget:.4f}.\n")
    print(f"cost = {args.vmer_weight:g}*vMER + MER\n")
    print(f"{'thr':>6}{'mean bal':>10}{'mean vMER':>11}{'mean MER':>10}"
          f"{'cost':>9}{'vMER saved':>12}{'':>5}")
    admissible = []
    cost_of = {}
    for t in sorted(per_t):
        bal, vm, me = per_t[t]
        cost_of[t] = args.vmer_weight * vm + me
        ok = me <= budget
        if ok:
            admissible.append((vm, t))
        print(f"{t:>6.2f}{bal:>10.3f}{vm:>11.3f}{me:>10.3f}{cost_of[t]:>9.4f}"
              f"{per_t[ref][1] - vm:>+12.4f}{'  ok' if ok else '  --':>5}")

    best_t = min(cost_of, key=lambda t: (cost_of[t], t))
    bal, vm, me = per_t[best_t]
    print(f"\nRecommended threshold: {best_t:g}  cost {cost_of[best_t]:.4f}")
    print(f"  vMER {vm:.4f} against {per_t[ref][1]:.4f} today, "
          f"MER {me:.4f} against {ref_mer:.4f}, balanced {bal:.3f}")
    if best_t == ref:
        print("  The current cut already minimises the cost — nothing to change.")

    if admissible:
        strict_vm, strict_t = min(admissible)
        if strict_t != best_t:
            print(f"\nCross-check: the lowest threshold that costs no specificity at "
                  f"all is {strict_t:g}\n(vMER {strict_vm:.4f}, MER "
                  f"{per_t[strict_t][2]:.4f}). It differs from the cost-minimiser, so "
                  f"the\nchoice between them is the weight — how much MER the extra "
                  f"vMER is worth.")
    floor_vm, floor_t = min((v[1], t) for t, v in per_t.items())
    print(f"\nFor reference, vMER alone is lowest at {floor_t:g} ({floor_vm:.4f}) "
          f"— the bottom of\nthe grid by construction, since calling more variants "
          f"can only reduce\nfalse-susceptible calls. It is the floor, not a "
          f"recommendation.")
    print("\nAnd read this against the genome-wide picture before acting: it "
          "measures only\nthe benefit at catalogue positions. The same threshold "
          "applies to every site, so\naudit_raw_heterozygosity.py gives what it "
          "costs in artefacts written into the\nsequences the models embed.")

    if args.out:
        p = pathlib.Path(args.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"\n[OK] {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
