#!/usr/bin/env python3
"""That drug's catalogued-mutation columns, optionally with the rule bit.

`cv_run_architecture.py` takes a single `--extra-feature-tsv`, so the combined
catalogue arm declared in §17 has to arrive as one file.

The declared arm (§17) is **columns only, no rule bit**: the catalogue's
observations without any part of its verdict. `--rule-tsv` is therefore optional.
Passing it produces the 76-feature variant instead, which §17 keeps available
under the same reporting terms.

The two blocks are not redundant, which is why dropping one is a real choice
rather than a simplification. §11 measured the asymmetry: the rule bit carries the
catalogue's **grading** -- which variant actually confers resistance -- and the
columns carry **which allele**, which is what places an isolate on the MIC ladder.
Presence/absence alone cannot recover the grading, because the variants deciding
vMER have too few carriers each to be learned from labels. Dropping the rule bit
is expected to cost accuracy on the S/R decision and gain it on the ladder.

Columns are selected for the drug from the long table's `drugs` field, never by
string-matching the column name. A column is named for its gene and change, not
for the drug it acts on, so name-matching would silently return nothing -- the
same position-instead-of-name mistake this project has made three times.

Mutations are likewise identified by `feature_class == "resistance_mutation"`
rather than by a `mutation::` name prefix. The first version of this script
assumed the prefix was `mut::` and selected nothing; a declared class is a fact
about the row, while a prefix is a fact about how the row happened to be named.
The long table's other class, `amr_gene_hit`, is gene presence and is deliberately
excluded -- the arm is sequence plus catalogued *mutations*.

No carrier floor by default. A floor would remove exactly the rare alleles the
arm exists to carry; the carrier distribution is printed instead so the choice is
visible rather than implicit.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys


def read_tsv(path: pathlib.Path):
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit(f"[ERROR] {path} is empty")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rule-tsv", default="",
                    help="rule_bit_<code>_tbprofiler.tsv. Optional: the declared "
                         "arm omits it, so that the features carry none of the "
                         "catalogue's verdict")
    ap.add_argument("--mutation-matrix", required=True,
                    help="tbprofiler_resistance_mutation_binary.tsv (panel-wide)")
    ap.add_argument("--long", required=True,
                    help="tbprofiler_amrfinder_style_long.tsv, for the drug map")
    ap.add_argument("--drug", required=True,
                    help="drug name as TB-Profiler writes it in the drugs field")
    ap.add_argument("--samples", default="",
                    help="a folds TSV; every one of its isolates must be covered")
    ap.add_argument("--feature-class", default="resistance_mutation",
                    help="the long table's feature_class to take. The default is "
                         "catalogued mutations; 'amr_gene_hit' would take gene "
                         "presence instead, which is a different arm")
    ap.add_argument("--min-carriers", type=int, default=0,
                    help="0 keeps every column, which is the default on purpose")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    drug = args.drug.strip().lower()
    long_rows = read_tsv(pathlib.Path(args.long))
    wanted, seen_drugs, seen_classes = set(), set(), set()
    for r in long_rows:
        cls = (r.get("feature_class") or "").strip()
        seen_classes.add(cls)
        drugs = {d.strip().lower() for d in (r.get("drugs") or "").split(";")
                 if d.strip() and d.strip() != "-"}
        seen_drugs |= drugs
        if drug in drugs and cls == args.feature_class:
            wanted.add(r["feature"])
    if not wanted:
        raise SystemExit(
            f"[ERROR] no {args.feature_class!r} rows tagged {args.drug!r} in "
            f"{args.long}.\n"
            f"        drugs present: {', '.join(sorted(seen_drugs)[:20])}\n"
            f"        feature classes present: {', '.join(sorted(seen_classes))}")

    mat = read_tsv(pathlib.Path(args.mutation_matrix))
    have = set(mat[0]) - {"sample_id"}
    cols = sorted(wanted & have)
    missing = wanted - have
    if missing:
        print(f"[WARN] {len(missing)} column(s) tagged for this drug are absent "
              f"from the matrix, e.g. {sorted(missing)[:3]}")

    mut = {r["sample_id"]: r for r in mat}
    if args.rule_tsv:
        rule_rows = read_tsv(pathlib.Path(args.rule_tsv))
        rule_cols = [c for c in rule_rows[0] if c != "sample_id"]
        clash = set(rule_cols) & set(cols)
        if clash:
            raise SystemExit(f"[ERROR] column name collision between the rule bit "
                             f"and the mutation matrix: {sorted(clash)}")
        rule = {r["sample_id"]: r for r in rule_rows}
        ids = sorted(set(rule) & set(mut))
    else:
        rule_cols, rule = [], {}
        ids = sorted(mut)
    if args.samples:
        need = {r["sample_id"] for r in read_tsv(pathlib.Path(args.samples))}
        gap = need - set(ids)
        if gap:
            raise SystemExit(
                f"[ERROR] {len(gap)} isolate(s) in {args.samples} are missing from "
                f"the rule bit or the mutation matrix, e.g. {sorted(gap)[:5]}.\n"
                f"        The combined arm must cover the same isolates as the arm "
                f"it is compared against, or the comparison is between cohorts.")
        ids = sorted(need)

    carriers = {c: sum(1 for i in ids if mut[i].get(c) not in ("", "0", None))
                for c in cols}
    kept = [c for c in cols if carriers[c] >= args.min_carriers]
    dropped = len(cols) - len(kept)

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", *rule_cols, *kept])
        for i in ids:
            w.writerow([i, *[rule[i][c] for c in rule_cols],
                        *[mut[i].get(c, "0") or "0" for c in kept]])
    if not rule_cols:
        print("[INFO] no rule bit: these features carry the catalogue's "
              "observations and none of its verdict (§17)")

    singles = sum(1 for c in kept if carriers[c] == 1)
    band = [sum(1 for c in kept if lo <= carriers[c] < hi)
            for lo, hi in ((1, 2), (2, 5), (5, 20), (20, 10 ** 9))]
    print(f"[INFO] {args.drug}: {len(rule_cols)} rule column(s) + {len(kept)} "
          f"mutation column(s) over {len(ids):,} isolates"
          + (f" ({dropped} below --min-carriers {args.min_carriers})" if dropped else ""))
    print(f"[INFO] carriers per mutation column: {band[0]} with 1, {band[1]} with "
          f"2-4, {band[2]} with 5-19, {band[3]} with 20+")
    if singles:
        print(f"[INFO] {singles} column(s) have a single carrier. Kept "
              f"deliberately: §11 found the variants that decide vMER are exactly "
              f"the ones too rare to learn from labels, so dropping them would "
              f"remove what the arm is for.")
    print(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
