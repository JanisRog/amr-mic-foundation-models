#!/usr/bin/env python3
"""Carve the untouched holdout, then the repeated-CV folds inside what is left.

Two separations, and the order matters. First 20% of the cohort is set aside and
not looked at again until the architecture and the config are both fixed --
without it, 27 architectures are compared on every sample and the final number is
a maximum over 27 arms rather than a measurement. Then the remaining 80% is cut
into 5 folds, three times over with different partitions, which is what makes the
comparison between architectures stable enough to act on.

Both cuts are stratified on the integer MIC class, because section 3 says bal_OT
is effectively decided by the rare classes -- a fold that happens to hold no
isolates of a rare class scores that class as a free win for every arm.

The repeats are genuinely different partitions, not three runs of the same one.
Repeating an identical partition measures the RF's own seed noise; repeating a
different one measures what we actually care about, which is how much the answer
depends on who landed in which fold.

Output is one TSV with a row per (sample, repeat) giving its fold, so downstream
steps never re-derive a split and every architecture provably sees the same ones.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import random
import statistics
import sys


def stratified_chunks(by_class: dict[float, list[str]], k: int,
                      rng: random.Random) -> list[list[str]]:
    """Deal each MIC class round-robin into k chunks, from a shuffled start.

    Dealing rather than slicing keeps the class counts within one of each other
    in every chunk, including for classes with fewer members than k.
    """
    chunks: list[list[str]] = [[] for _ in range(k)]
    for cls in sorted(by_class):
        members = list(by_class[cls])
        rng.shuffle(members)
        offset = rng.randrange(k)          # or class 1 always starts in fold 0
        for i, s in enumerate(members):
            chunks[(i + offset) % k].append(s)
    return chunks


def sized_chunks(by_class: dict[float, list[str]], sizes: list[int],
                 rng: random.Random) -> list[list[str]]:
    """Stratified deal to a given size profile rather than to equal folds.

    The control for a leave-one-group-out result: reproduce the fold sizes exactly,
    with random membership. Any difference that remains is the grouping, not the
    training-set size that the grouping happened to impose.

    Each class is split across folds in proportion to the target sizes, largest
    remainder first, so the stratification survives the reshaping.
    """
    total = sum(sizes)
    chunks: list[list[str]] = [[] for _ in sizes]
    for cls in sorted(by_class):
        members = list(by_class[cls])
        rng.shuffle(members)
        n = len(members)
        quota = [n * s / total for s in sizes]
        base = [int(q) for q in quota]
        for i in sorted(range(len(sizes)), key=lambda i: -(quota[i] - base[i]))[
                :n - sum(base)]:
            base[i] += 1
        pos = 0
        for i, b in enumerate(base):
            chunks[i].extend(members[pos:pos + b])
            pos += b
    return chunks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mic-from", required=True,
                    help="sample_splits_*.tsv used as the sample + MIC source")
    ap.add_argument("--restrict-samples", default="",
                    help="cohort file, e.g. rif_susceptible_samples.txt")
    ap.add_argument("--holdout-frac", type=float, default=0.2)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--holdout-force", default="",
                    help="isolates that must be in the holdout — e.g. an earlier "
                         "holdout being carried forward, so nothing already used "
                         "for a design decision can leak into the new one")
    ap.add_argument("--holdout-exclude", default="",
                    help="isolates that must NOT be in the holdout, because they "
                         "were used to choose something")
    ap.add_argument("--allele-tsv", action="append", default=[],
                    help="sample_id<TAB>seq_id for a focus target; whole allele "
                         "groups are then assigned to a fold together, so every "
                         "validation isolate carries a sequence absent from that "
                         "fold's training set (section 4)")
    ap.add_argument("--pin-train-above", type=int, default=0,
                    help="with --allele-tsv, allele groups larger than this are "
                         "never held out; one wild-type allele can be half the "
                         "cohort and no K-way split can balance it")
    ap.add_argument("--fold-sizes", default="",
                    help="comma-separated target sizes, e.g. 2019,1130,251,223 — "
                         "random folds shaped like a leave-one-group-out run, so "
                         "the two differ only in grouping and not in training size")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    with open(args.mic_from) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    col = "log2_mic" if "log2_mic" in rows[0] else list(rows[0])[-1]
    mic = {r["sample_id"]: float(r[col]) for r in rows
           if r.get(col) not in (None, "", "nan")}
    print(f"[INFO] {len(mic):,} samples with MIC")

    if args.restrict_samples:
        keep = {l.strip() for l in open(args.restrict_samples) if l.strip()}
        mic = {s: v for s, v in mic.items() if s in keep}
        print(f"[INFO] {len(mic):,} after the cohort restriction")
    if not mic:
        raise SystemExit("[ERROR] no samples left")

    rng = random.Random(args.seed)
    forced = ({l.strip() for l in open(args.holdout_force) if l.strip()} & set(mic)
              if args.holdout_force else set())
    barred = ({l.strip() for l in open(args.holdout_exclude) if l.strip()}
              if args.holdout_exclude else set())
    if forced:
        print(f"[INFO] {len(forced):,} isolates forced into the holdout")
    if barred:
        print(f"[INFO] {len(barred & set(mic)):,} isolates barred from it")
    overlap = forced & barred
    if overlap:
        raise SystemExit(f"[ERROR] {len(overlap):,} isolates are both forced into "
                         f"and barred from the holdout")

    # holdout: one chunk of a 1/frac-way stratified deal over what is still free
    free = {s: v for s, v in mic.items() if s not in forced and s not in barred}
    by_class: dict[float, list[str]] = collections.defaultdict(list)
    for s, v in free.items():
        by_class[v].append(s)
    want = max(0, round(len(mic) * args.holdout_frac) - len(forced))
    n_chunks = max(2, round(len(free) / want)) if want else 10 ** 9
    chunks = stratified_chunks(by_class, min(n_chunks, max(2, len(free))), rng) \
        if want else [[]]
    holdout = sorted(forced | set(chunks[0]))
    dev = sorted(s for s in mic if s not in set(holdout))
    print(f"[INFO] holdout {len(holdout):,} ({len(holdout)/len(mic):.0%}), "
          f"development {len(dev):,}")

    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "holdout_samples.txt").write_text("\n".join(holdout) + "\n")
    (out / "dev_samples.txt").write_text("\n".join(dev) + "\n")

    # The final evaluation needs one fixed feature set, not a per-fold one: the
    # projection and the mRMR ranking are fitted on all of development and applied
    # unchanged to the holdout. That is exactly this file -- every development
    # isolate sits at fold -1, which no fold index ever matches, so it is training
    # in the single partition; the holdout is fold 0 and is scored once.
    with (out / "cv_folds_final.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "repeat", "fold", "log2_mic"])
        for x in sorted(dev):
            w.writerow([x, 0, -1, f"{mic[x]:g}"])
        for x in sorted(holdout):
            w.writerow([x, 0, 0, f"{mic[x]:g}"])
    print(f"[OK] {out/'cv_folds_final.tsv'} — {len(dev):,} train, "
          f"{len(holdout):,} held out, one partition")
    print("     Use it with --only-fold 0:0 and a single-config JSON. Scoring more "
          "than one\n     config on the holdout is the max-over-configs bias "
          "(section 3) at the last step.")

    # allele grouping: a "unit" is one isolate, or one allele group when
    # --allele-tsv is given, and units are dealt to folds intact
    unit_of = {s: s for s in dev}
    pinned: set[str] = set()
    if args.allele_tsv:
        maps = []
        for spec in args.allele_tsv:
            # PATH or PATH#column — the group file may be a bare two-column map or
            # a wider table with a header, like sample_lineage.tsv
            path, _, want = spec.partition("#")
            with open(path) as fh:
                rd = list(csv.reader(fh, delimiter="\t"))
            if not rd:
                raise SystemExit(f"[ERROR] {path} is empty")
            head, body = rd[0], rd
            col = 1
            if head and head[0].strip().lower() in ("sample_id", "sample"):
                body = rd[1:]
                if want:
                    if want not in head:
                        raise SystemExit(f"[ERROR] no column {want!r} in {path}; "
                                         f"columns are {head}")
                    col = head.index(want)
                print(f"[INFO] {pathlib.Path(path).name}: header found, grouping "
                      f"on column {head[col]!r}")
            elif want:
                raise SystemExit(f"[ERROR] {path} has no header, so #{want} cannot "
                                 f"be resolved")
            d = {r[0].strip(): r[col].strip() for r in body
                 if len(r) > col and r[0].strip()}
            maps.append(d)
            counts = collections.Counter(d.values())
            print(f"[INFO] {pathlib.Path(path).name}: {len(d):,} samples, "
                  f"{len(counts):,} groups; largest {counts.most_common(1)[0][1]:,}")
        covered = [s for s in dev if all(s in d for d in maps)]
        print(f"[INFO] {len(covered):,} of {len(dev):,} development isolates have "
              f"every focus allele")
        groups: dict[tuple, list[str]] = collections.defaultdict(list)
        for s in covered:
            groups[tuple(d[s] for d in maps)].append(s)
        for g, members in groups.items():
            for s in members:
                unit_of[s] = "grp::" + "|".join(g)
        if args.pin_train_above:
            big = {g for g, m in groups.items() if len(m) > args.pin_train_above}
            pinned = {s for g in big for s in groups[g]}
            print(f"[INFO] pinning {len(big)} allele groups ({len(pinned):,} "
                  f"isolates) to training in every fold")

    rotated = [s for s in dev if s not in pinned]
    units: dict[str, list[str]] = collections.defaultdict(list)
    for s in rotated:
        units[unit_of[s]].append(s)
    # a unit is placed by its modal MIC class, so the stratification still holds
    unit_class = {u: statistics.median(mic[s] for s in m) for u, m in units.items()}
    by_unit_class: dict[float, list[str]] = collections.defaultdict(list)
    for u, c in unit_class.items():
        by_unit_class[c].append(u)

    # Leave-one-group-out: when there are no more units than folds, the intent is
    # one group per fold, not a class-balanced deal. Dealing by MIC class would put
    # two lineages in one fold and leave another holding only the 16 unclassified
    # isolates, where most MIC classes are empty and bal_OT is undefined.
    one_per_fold = 0 < len(units) <= args.folds and not args.fold_sizes
    if one_per_fold:
        print(f"[INFO] {len(units)} groups for {args.folds} folds — assigning one "
              f"group per fold (leave-one-group-out)")

    assign: list[tuple[str, int, int, float]] = []
    for rep in range(args.repeats):
        rep_rng = random.Random(args.seed + 1000 + rep)
        if args.fold_sizes:
            want = [int(x) for x in args.fold_sizes.split(",") if x.strip()]
            if len(want) != args.folds:
                raise SystemExit(f"[ERROR] --fold-sizes gives {len(want)} sizes for "
                                 f"{args.folds} folds")
            fold_units = sized_chunks(by_unit_class, want, rep_rng)
        elif one_per_fold:
            big = sorted(units, key=lambda u: -len(units[u]))
            fold_units = [[u] for u in big] + [[] for _ in range(args.folds - len(big))]
        else:
            fold_units = stratified_chunks(by_unit_class, args.folds, rep_rng)
        member_of = {}
        for f, us in enumerate(fold_units):
            for u in us:
                member_of[u] = f
        for s in dev:
            f = -1 if s in pinned else member_of[unit_of[s]]
            assign.append((s, rep, f, mic[s]))
        sizes = [sum(1 for s, r, ff, _ in assign if r == rep and ff == f)
                 for f in range(args.folds)]
        print(f"[INFO] repeat {rep}: fold sizes {sizes}"
              + (f"  (+{len(pinned):,} pinned to training)" if pinned else ""))
        if args.fold_sizes:
            want = [int(x) for x in args.fold_sizes.split(",") if x.strip()]
            if sizes != want:
                print(f"[INFO]   requested {want}; each MIC class is split "
                      f"independently, so small folds land within a couple of "
                      f"isolates rather than exactly")
        if one_per_fold:
            for f, us in enumerate(fold_units):
                if us:
                    print(f"[INFO]   fold {f} = {us[0].replace('grp::','')}"
                          f" ({sizes[f]:,})")
        empty = [f for f, n in enumerate(sizes) if n == 0]
        if empty:
            print(f"[WARN] folds {empty} are empty and will be skipped downstream")

    with (out / "cv_folds.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "repeat", "fold", "log2_mic"])
        for row in sorted(assign):
            w.writerow(row)

    # every sample appears exactly once per repeat, or a fold is not a partition
    per_rep = collections.Counter((s, r) for s, r, _, _ in assign)
    if set(per_rep.values()) != {1}:
        raise SystemExit("[ERROR] a sample lands in two folds of one repeat")
    print(f"\n[OK] {out/'cv_folds.tsv'} — {len(assign):,} rows, "
          f"{args.repeats} x {args.folds} partitions")

    print(f"\n{'log2 MIC':>10}{'dev':>8}{'holdout':>9}{'holdout %':>11}"
          f"{'min fold':>10}")
    rep0 = collections.defaultdict(lambda: collections.Counter())
    for s, r, f, v in assign:
        if r == 0:
            rep0[v][f] += 1
    for cls in sorted(by_class):
        nd = sum(1 for s in dev if mic[s] == cls)
        nh = sum(1 for s in holdout if mic[s] == cls)
        tot = nd + nh
        smallest = min(rep0[cls].get(f, 0) for f in range(args.folds)) if nd else 0
        print(f"{cls:>10.1f}{nd:>8}{nh:>9}{(nh/tot if tot else 0):>10.1%}"
              f"{smallest:>10}")
    print("\n'min fold' is the smallest per-fold count of that class in repeat 0.\n"
          "A zero there means bal_OT is undefined for that class in that fold and\n"
          "the class silently drops out of the average.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
