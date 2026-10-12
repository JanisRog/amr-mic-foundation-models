#!/usr/bin/env python3
"""Supplementary table: which features the M. tuberculosis model needs.

Six feature sets on identical development folds, crossing three questions --
embeddings or not, TB-Profiler's call or not, catalogued mutation columns or not:

  emb_rule      embeddings + TB-Profiler call                 (the reported set)
  emb_snp       embeddings + mutation columns
  emb_snp_rule  embeddings + mutation columns + TB-Profiler call
  emb           embeddings only
  snp_rule      mutation columns + TB-Profiler call            (no embeddings)
  snp           mutation columns only                          (no embeddings)

plus TB-Profiler's call scored on its own, with no model, from the baseline rows.

Embedding arms are read on the declared primary encoder (Evo 2); the two arms
without embeddings have no encoder. Every arm is its 1-SE-selected configuration,
so no arm is credited with its best draw over 1,000.

The verdict rule is §20's, fixed before any of the four new arms was run:

  "embeddings + TB-Profiler call is the best of the six" is printed only if it has
  the lowest median development cost across the drugs with a verdict AND no other
  arm has lower cost in a paired Wilcoxon test at P < 0.05. Otherwise the table
  names the arm that did better. Either way the reported model does not change:
  it was selected on isoniazid's development folds before this comparison, and its
  held-out set has been scored, so switching now would be a second look.

Bedaquiline and delamanid are shown and excluded from the panel summary, as in
every other panel comparison: the catalogue misses most of their resistant
isolates and their costs are on another scale.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics
import sys

ARMS = [
    ("emb_rule", "Embeddings + TB-Profiler call", "evo2",
     ["both_final1k", "tbrule1k", "final1k"]),
    ("emb_snp", "Embeddings + mutation columns", "evo2", ["both_arm75", "arm75"]),
    ("emb_snp_rule", "Embeddings + mutation columns + TB-Profiler call", "evo2",
     ["both_arm76", "arm76"]),
    ("emb", "Embeddings only", "evo2", ["both_embonly1k", "embonly1k"]),
    ("snp_rule", "Mutation columns + TB-Profiler call", "catalogue", ["snprule1k"]),
    ("snp", "Mutation columns only", "catalogue", ["snponly1k"]),
]
REF = "emb_rule"
METRICS = ("cost", "vMER", "MER", "SR_f1", "bal_OT_common", "MAE")
LOWER = {"cost", "vMER", "MER", "MAE"}
NO_VERDICT = {"bedaquiline", "delamanid"}
SUMMARIES = {
    "isoniazid_final": "cv_final_summaries", "rifampicin": "cv_rif_summaries",
    "ethambutol": "cv_emb_summaries", "ethionamide": "cv_eth_summaries",
    "kanamycin": "cv_kan_summaries", "amikacin": "cv_ami_summaries",
    "levofloxacin": "cv_lev_summaries", "moxifloxacin": "cv_mxf_summaries",
    "bedaquiline": "cv_bdq_summaries", "delamanid": "cv_dlm_summaries",
}


def pretty(d):
    return d[:-len("_final")] if d.endswith("_final") else d


def read_arm(root, drug, model, suffixes):
    for suf in suffixes:
        f = root / f"cv_{drug}" / "runs" / f"{model}_{suf}" / "folds.tsv"
        if not f.exists():
            continue
        rows = [r for r in csv.DictReader(open(f), delimiter="\t")
                if r["aggregation"] == "median"]
        if not rows:
            continue
        out = {"n_features": int(float(rows[0].get("n_features", "nan") or "nan"))
               if rows[0].get("n_features") else None}
        for k in METRICS:
            if k == "cost":
                v = [2 * float(r["vMER"]) + float(r["MER"]) for r in rows]
            else:
                v = [float(r[k]) for r in rows if r.get(k) not in ("", "NA", None)]
            out[k] = statistics.mean(v) if v else None
        return out
    return None


def tbprofiler(summ, drug):
    p = summ / SUMMARIES[drug] / "baseline_rows.tsv"
    if not p.exists():
        return None
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["label"] == "TB-Profiler 6.3.0" and r["split"] == "development_cv":
            return {k: float(r[k]) for k in ("cost", "vMER", "MER", "SR_f1")}
    return None


def wilcoxon(d):
    try:
        from scipy.stats import wilcoxon as w
        if len(d) >= 5 and any(x != 0 for x in d):
            return float(w(d).pvalue)
    except Exception:
        pass
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="slim_arms (selected configs)")
    ap.add_argument("--summaries", required=True, help="dir holding cv_*_summaries/")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    root, summ = pathlib.Path(a.root), pathlib.Path(a.summaries)

    drugs = [d for d in SUMMARIES if (root / f"cv_{d}").is_dir()]
    cell = {(d, k): read_arm(root, d, m, s) for d in drugs for k, _, m, s in ARMS}
    tbp = {d: tbprofiler(summ, d) for d in drugs}
    complete = [d for d in drugs if all(cell[(d, k)] for k, *_ in ARMS)]
    panel = [d for d in complete if d not in NO_VERDICT]
    missing = {d: [k for k, *_ in ARMS if not cell[(d, k)]] for d in drugs
               if d not in complete}

    # ---- long table --------------------------------------------------------
    with open(f"{a.out_prefix}_by_drug.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["drug", "arm", "label", "basis", "n_features", *METRICS])
        for d in drugs:
            for k, label, model, _ in ARMS:
                v = cell[(d, k)]
                w.writerow([pretty(d), k, label, model if model != "catalogue" else "no encoder",
                            (v or {}).get("n_features", "")]
                           + [f"{v[m]:.4f}" if v and v.get(m) is not None else "" for m in METRICS])
            t = tbp[d]
            w.writerow([pretty(d), "tbprofiler", "TB-Profiler call alone (no model)",
                        "catalogue", 1]
                       + [f"{t[m]:.4f}" if t and m in t else "" for m in METRICS])

    if not panel:
        print("[INFO] no drug has all six arms yet; nothing to summarise.")
        for d, v in missing.items():
            print(f"   {pretty(d):<14} missing {', '.join(v)}")
        print(f"[OK] {a.out_prefix}_by_drug.tsv written with what exists")
        return 0

    # ---- panel summary ------------------------------------------------------
    lines = []
    say = lines.append
    say(f"Feature-set ablation, M. tuberculosis, development folds, Evo 2 for "
        f"embedding arms")
    say(f"Panel: {len(panel)} drugs with all six arms and a verdict "
        f"({', '.join(pretty(d) for d in panel)})")
    if missing:
        say("Incomplete, left out of the panel: "
            + "; ".join(f"{pretty(d)} lacks {', '.join(v)}" for d, v in missing.items()))
    say("")
    hdr = f"{'arm (median over panel)':<52}" + "".join(f"{m:>17}" for m in METRICS)
    say(hdr)
    med = {}
    for k, label, *_ in ARMS:
        med[k] = {m: statistics.median([cell[(d, k)][m] for d in panel
                                        if cell[(d, k)][m] is not None])
                  for m in METRICS} if panel else {}
        say(f"{label:<52}" + "".join(f"{med[k].get(m, float('nan')):>17.4f}" for m in METRICS))
    tmed = {m: statistics.median([tbp[d][m] for d in panel if tbp[d] and m in tbp[d]])
            for m in ("cost", "vMER", "MER", "SR_f1")} if panel else {}
    say(f"{'TB-Profiler call alone (no model)':<52}"
        + "".join(f"{tmed[m]:>17.4f}" if m in tmed else f"{'—':>17}" for m in METRICS))

    say("")
    say(f"Paired against '{dict((k, l) for k, l, *_ in ARMS)[REF]}', across the panel:")
    say(f"{'arm':<52}{'Δcost med':>11}{'lower':>7}{'P':>9}{'Δbal_OT_c':>11}{'P':>9}")
    beats = []
    for k, label, *_ in ARMS:
        if k == REF:
            continue
        dc = [cell[(d, k)]["cost"] - cell[(d, REF)]["cost"] for d in panel]
        db = [cell[(d, k)]["bal_OT_common"] - cell[(d, REF)]["bal_OT_common"] for d in panel]
        pc, pb = wilcoxon(dc), wilcoxon(db)
        lower = sum(1 for x in dc if x < 0)
        if pc is not None and pc < 0.05 and statistics.median(dc) < 0:
            beats.append(label)
        say(f"{label:<52}{statistics.median(dc):>+11.4f}{f'{lower}/{len(dc)}':>7}"
            f"{(f'{pc:.4f}' if pc is not None else '—'):>9}"
            f"{statistics.median(db):>+11.4f}{(f'{pb:.4f}' if pb is not None else '—'):>9}")

    say("")
    order = sorted((k for k, *_ in ARMS), key=lambda k: med[k]["cost"])
    best = order[0]
    labels = {k: l for k, l, *_ in ARMS}
    say("Ranked by median development cost: "
        + " < ".join(f"{labels[k]} ({med[k]['cost']:.4f})" for k in order))
    say("")
    say("Verdict (§20 rule, fixed before the new arms were run):")
    if best == REF and not beats:
        say(f"  Embeddings + TB-Profiler call has the lowest median development cost of "
            f"the six, and no arm has significantly lower cost (n = {len(panel)}).")
    else:
        if best != REF:
            say(f"  '{labels[best]}' has the lowest median development cost, not the "
                f"reported set.")
        for b in beats:
            say(f"  '{b}' has significantly lower cost than the reported set.")
        say("  The reported model is unchanged: it was selected on isoniazid's development "
            "folds before this comparison and its held-out set has been scored.")
    say(f"  Second-lowest median cost: '{labels[order[1]]}'.")

    text = "\n".join(lines)
    pathlib.Path(f"{a.out_prefix}_summary.txt").write_text(text + "\n")
    print(text)
    print(f"\n[OK] {a.out_prefix}_by_drug.tsv, {a.out_prefix}_summary.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
