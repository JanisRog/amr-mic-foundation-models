#!/usr/bin/env python3
"""Every isolate's catalogue S/R call, at each allele fraction, as a table.

`sweep_consensus_threshold.py` answers where the threshold should sit and reports
means over drugs. This answers the per-isolate question underneath it: which
isolates change call, in which direction, and whether the phenotype agrees -- and
it emits the calls themselves, which is the part that matters beyond diagnostics.

The rule bit in the model's feature set *is* one of these calls. It currently comes
from TB-Profiler at its default 0.1 while the embeddings are built from a ~0.9
consensus and the mutation columns from allele fraction 0.10, so the three blocks
have been describing the same isolate differently: the rule bit can say resistant
on a variant the embedding renders as wild type. Emitting the call at each
threshold is what lets all three agree.

Detection is deliberately identical to the sweep -- same catalogue keys, same
QUAL/DP floors, same AO/DP fraction, same haplotype atomization -- by importing
from it rather than reimplementing, so the two can never drift into disagreeing
about what a call is.

An isolate with no raw VCF is absent from the output rather than called
susceptible; a missing measurement is not a negative one.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from sweep_consensus_threshold import (  # noqa: E402
    load_catalogue, isolate_fractions,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snippy-runs", required=True)
    ap.add_argument("--catalogue", required=True,
                    help="call_sources_by_variant.tsv")
    ap.add_argument("--thresholds", default="0.1,0.2,0.9")
    ap.add_argument("--min-qual", type=float, default=100.0)
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--drugs", default="",
                    help="restrict to these drugs, comma separated")
    ap.add_argument("--mic-tsv", default="",
                    help="optional; adds a phenotype column so flips can be "
                         "scored rather than only counted")
    ap.add_argument("--mic-suffix", default="_int")
    ap.add_argument("--sr-boundary", default="isoniazid=-3.0",
                    help="drug=log2 boundary, repeatable with commas; an isolate "
                         "is phenotypically resistant above it")
    ap.add_argument("--reference-threshold", type=float, default=0.9,
                    help="the threshold flips are counted against")
    ap.add_argument("--score-samples", default="",
                    help="isolates whose phenotype may be used to CHOOSE the "
                         "threshold -- a development sample list, or a fold TSV "
                         "with a sample_id column. The call matrices always "
                         "cover every isolate, because they are features; only "
                         "the scoring is restricted. The threshold is a "
                         "hyperparameter and it feeds the embeddings, so picking "
                         "it against held-out phenotypes contaminates every "
                         "downstream number (section 3). Per drug, since each "
                         "drug has its own split.")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    ts = [float(x) for x in args.thresholds.split(",") if x.strip()]
    if args.reference_threshold not in ts:
        ts.append(args.reference_threshold)
    ts.sort()

    cat = load_catalogue(pathlib.Path(args.catalogue))
    positions = {p for p, _, _ in cat}
    want_drugs = {d.strip() for d in args.drugs.split(",") if d.strip()}
    all_drugs = sorted({d for s in cat.values() for d in s
                        if not want_drugs or d in want_drugs})
    print(f"[INFO] catalogue: {len(cat):,} variant keys, {len(positions):,} "
          f"positions, {len(all_drugs)} drugs")

    bounds = {}
    for part in args.sr_boundary.split(","):
        if "=" in part:
            k, v = part.split("=", 1)
            bounds[k.strip()] = float(v)
    # The catalogue names drugs by code (INH), not by name (isoniazid). A
    # boundary that matches nothing used to leave the phenotype columns silently
    # blank, which reads as "no flips were scoreable" rather than as a typo.
    unmatched = sorted(set(bounds) - set(all_drugs))
    if unmatched:
        raise SystemExit(
            f"[ERROR] --sr-boundary names {unmatched} but the catalogue's drugs "
            f"are {all_drugs}. Use the catalogue's codes.")

    mic = {}
    if args.mic_tsv:
        with open(args.mic_tsv) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sid = str(r.get("sample_id", "")).strip()
                if not sid:
                    continue
                vals = {}
                for d in all_drugs:
                    raw = r.get(d + args.mic_suffix, r.get(d[:3].upper() +
                                                           args.mic_suffix, ""))
                    if raw not in ("", "NA", None):
                        try:
                            vals[d] = float(raw)
                        except ValueError:
                            pass
                if vals:
                    mic[sid] = vals

    root = pathlib.Path(args.snippy_runs)
    dirs = sorted(d for d in root.iterdir()
                  if d.is_dir() and (d / "snps.raw.vcf").exists())
    print(f"[INFO] {len(dirs):,} isolates with a raw VCF")

    calls: dict[float, dict[str, dict[str, int]]] = {t: {} for t in ts}
    for i, d in enumerate(dirs, 1):
        if i % 2000 == 0:
            print(f"[INFO] {i:,}/{len(dirs):,}")
        fr = isolate_fractions(d / "snps.raw.vcf", positions,
                               args.min_qual, args.min_dp)
        best: dict[str, float] = collections.defaultdict(float)
        for key, f in fr.items():
            for drug in cat.get(key, ()):
                if not want_drugs or drug in want_drugs:
                    if f > best[drug]:
                        best[drug] = f
        for t in ts:
            calls[t][d.name] = {dr: int(best.get(dr, 0.0) >= t)
                                for dr in all_drugs}

    pre = pathlib.Path(args.out_prefix)
    pre.parent.mkdir(parents=True, exist_ok=True)
    for t in ts:
        p = pre.with_name(pre.name + f"_af{int(round(t*100)):03d}.tsv")
        with p.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sample_id"] + all_drugs)
            for sid in sorted(calls[t]):
                w.writerow([sid] + [calls[t][sid][dr] for dr in all_drugs])
        n_r = {dr: sum(v[dr] for v in calls[t].values()) for dr in all_drugs}
        print(f"[OK] {p}  resistant: " +
              " ".join(f"{dr[:3]}={n_r[dr]}" for dr in all_drugs))

    scoreable: set[str] | None = None
    if args.score_samples:
        p = pathlib.Path(args.score_samples)
        first = p.read_text().splitlines()[0].split("\t")
        if "sample_id" in first:
            i = first.index("sample_id")
            fold_i = first.index("fold") if "fold" in first else None
            scoreable = set()
            for line in p.read_text().splitlines()[1:]:
                f = line.split("\t")
                if len(f) <= i:
                    continue
                if fold_i is not None and len(f) > fold_i:
                    try:
                        if int(f[fold_i]) < 0:      # holdout rows are marked < 0
                            continue
                    except ValueError:
                        pass
                scoreable.add(f[i].strip())
        else:
            scoreable = {l.split("\t")[0].strip()
                         for l in p.read_text().splitlines() if l.strip()}
        print(f"[INFO] threshold scored on {len(scoreable):,} development "
              f"isolates from {p.name}; calls emitted for all "
              f"{len(dirs):,}")

    ref = args.reference_threshold
    rows = []
    for t in ts:
        if t == ref:
            continue
        for dr in all_drugs:
            pool = [s for s in calls[t]
                    if scoreable is None or s in scoreable]
            gained = [s for s in pool
                      if calls[t][s][dr] and not calls[ref][s][dr]]
            lost = [s for s in pool
                    if not calls[t][s][dr] and calls[ref][s][dr]]
            row = {"threshold": t, "drug": dr, "n_gained_R": len(gained),
                   "n_lost_R": len(lost)}
            if mic and dr in bounds:
                b = bounds[dr]
                gt = [s for s in gained if s in mic and dr in mic[s]]
                row["gained_with_mic"] = len(gt)
                row["gained_truly_R"] = sum(1 for s in gt if mic[s][dr] > b)
                row["gained_truly_S"] = sum(1 for s in gt if mic[s][dr] <= b)
            rows.append(row)
    p = pre.with_name(pre.name + "_flips.tsv")
    cols = sorted({k for r in rows for k in r})
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["threshold", "drug"] +
                                      [c for c in cols
                                       if c not in ("threshold", "drug")])
        w.writeheader()
        w.writerows(rows)
    print(f"\n[OK] {p}")
    print(f"\n{'drug':<16}{'t':>6}{'+R':>7}{'-R':>6}{'+R w/MIC':>10}"
          f"{'truly R':>9}{'truly S':>9}")
    for r in rows:
        print(f"{r['drug']:<16}{r['threshold']:>6.2f}{r['n_gained_R']:>7}"
              f"{r['n_lost_R']:>6}{r.get('gained_with_mic',''):>10}"
              f"{r.get('gained_truly_R',''):>9}{r.get('gained_truly_S',''):>9}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
