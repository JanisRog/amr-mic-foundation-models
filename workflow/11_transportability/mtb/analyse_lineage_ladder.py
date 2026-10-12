#!/usr/bin/env python3
"""§21: does the embeddings' advantage on the MIC ladder survive leaving a lineage out?

On development folds, embeddings + TB-Profiler call place isolates on the MIC
ladder better than mutation columns + the same call, at equal resistance-call
performance (provisional, six drugs: bal_OT_common +0.009, MAE -0.020, both 6/6).
Whole-gene embeddings also carry each isolate's lineage background, which the
catalogue's mutation columns do not, and §2 found the ladder to be the
lineage-fragile quantity. If part of the advantage is lineage recognition, it
shrinks when a whole lineage is held out.

Both arms are refits of their development-selected configuration on the same
folds: run_panel_lineage.sbatch (embeddings + call, as §19's baseline) and
run_lineage_catalogue.sbatch (mutation columns + call). Each lineage split is
read against random folds of identical sizes, so the comparison is population
shift, not training size.

The rule, fixed in §21 before either catalogue refit was run:

  advantage  = (embeddings + call) - (mutation columns + call), on bal_OT_common
  extra drop = how much more embeddings + call loses than mutation columns + call,
               going from size-matched to lineage folds

  SURVIVES   median advantage on lineage folds > 0, and median extra drop <= 0.005
  WEAKENED   median advantage on lineage folds > 0, but extra drop > 0.005
  DOES NOT   median advantage on lineage folds <= 0

over the drugs with a verdict (bedaquiline and delamanid excluded, as in every
panel comparison). MAE is reported the same way, as a second view.

Read once, with §19, on all ten drugs. Without them it prints the gaps and stops.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics
import sys

DRUGS = ["isoniazid_final", "rifampicin", "ethambutol", "ethionamide", "kanamycin",
         "amikacin", "levofloxacin", "moxifloxacin", "bedaquiline", "delamanid"]
NO_VERDICT = {"bedaquiline", "delamanid"}
RUNS = {  # arm -> (lineage run, size-matched run)
    "emb": ("runs/evo2_lolo_baseline", "sizematched/runs/evo2_sm_baseline"),
    "snp": ("runs/catalogue_lolo_snprule", "sizematched/runs/catalogue_sm_snprule"),
}
TOL = 0.005


def read(path: pathlib.Path):
    f = path / "folds.tsv"
    if not f.exists():
        return None
    rows = [r for r in csv.DictReader(open(f), delimiter="\t")
            if r["aggregation"] == "median"]
    if not rows:
        return None
    out = {}
    for k in ("bal_OT_common", "MAE"):
        v = [float(r[k]) for r in rows if r.get(k) not in ("", "NA", None)]
        out[k] = statistics.mean(v) if v else None
    out["cost"] = statistics.mean(2 * float(r["vMER"]) + float(r["MER"]) for r in rows)
    return out


def pretty(d):
    return d[:-len("_final")] if d.endswith("_final") else d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lineage-root", required=True,
                    help="directory holding cv_<drug>/lineage/...")
    ap.add_argument("--out", default="")
    ap.add_argument("--i-accept-reading-early", action="store_true",
                    help="report on an incomplete panel, spending the single read")
    a = ap.parse_args()
    root = pathlib.Path(a.lineage_root)

    cell, missing = {}, []
    for d in DRUGS:
        L = root / f"cv_{d}" / "lineage"
        for arm, (lin, sm) in RUNS.items():
            cell[(d, arm, "lin")] = read(L / lin)
            cell[(d, arm, "sm")] = read(L / sm)
        gaps = [f"{arm}/{s}" for arm in RUNS for s in ("lin", "sm") if not cell[(d, arm, s)]]
        if gaps:
            missing.append(f"{pretty(d)} ({', '.join(gaps)})")
    if missing and not a.i_accept_reading_early:
        print("[STOP] §21 is read once, with §19, on all ten drugs. Missing:")
        for m in missing:
            print(f"   {m}")
        return 1

    panel = [d for d in DRUGS if d not in NO_VERDICT
             and all(cell[(d, arm, s)] for arm in RUNS for s in ("lin", "sm"))]
    rows, lines = [], []
    say = lines.append
    say(f"§21 — ladder advantage of embeddings over mutation columns (both with "
        f"TB-Profiler's call), {len(panel)} drugs")
    if missing:
        say("READ EARLY on an incomplete panel; missing: " + "; ".join(missing))
    for metric, higher in (("bal_OT_common", True), ("MAE", False)):
        sign = 1 if higher else -1      # orient so that positive = embeddings better
        say("")
        say(f"{metric}   (positive advantage = embeddings better)")
        say(f"{'drug':<14}{'adv size-mat':>13}{'adv lineage':>13}{'emb drop':>10}"
            f"{'snp drop':>10}{'extra drop':>12}")
        adv_sm, adv_lin, extra = [], [], []
        for d in panel:
            es, el = cell[(d, "emb", "sm")][metric], cell[(d, "emb", "lin")][metric]
            ss, sl = cell[(d, "snp", "sm")][metric], cell[(d, "snp", "lin")][metric]
            a_sm, a_lin = sign * (es - ss), sign * (el - sl)
            drop_e, drop_s = sign * (es - el), sign * (ss - sl)   # positive = got worse
            x = drop_e - drop_s
            adv_sm.append(a_sm); adv_lin.append(a_lin); extra.append(x)
            say(f"{pretty(d):<14}{a_sm:>+13.4f}{a_lin:>+13.4f}{drop_e:>+10.4f}"
                f"{drop_s:>+10.4f}{x:>+12.4f}")
            rows.append([d, metric, es, el, ss, sl, a_sm, a_lin, x])
        m_sm, m_lin, m_x = (statistics.median(v) for v in (adv_sm, adv_lin, extra))
        pos = sum(1 for v in adv_lin if v > 0)
        say(f"{'median':<14}{m_sm:>+13.4f}{m_lin:>+13.4f}{'':>20}{m_x:>+12.4f}"
            f"   (advantage positive on lineage folds in {pos}/{len(panel)})")
        if metric == "bal_OT_common":
            if m_lin > 0 and m_x <= TOL:
                v = "SURVIVES — the ladder advantage holds when a lineage is held out."
            elif m_lin > 0:
                v = (f"WEAKENED — the advantage keeps its sign but embeddings lose "
                     f"{m_x:+.4f} more than mutation columns across lineages; part of "
                     f"it is lineage-borne.")
            else:
                v = ("DOES NOT SURVIVE — on lineage folds the embeddings no longer "
                     "place isolates better; the development advantage was lineage "
                     "recognition.")
            verdict = v
    say("")
    say(f"Verdict (§21 rule, on bal_OT_common): {verdict}")
    text = "\n".join(lines)
    print(text)
    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["drug", "metric", "emb_sizematched", "emb_lineage",
                        "snp_sizematched", "snp_lineage", "adv_sizematched",
                        "adv_lineage", "extra_drop_emb"])
            for r in rows:
                w.writerow([r[0], r[1]] + [f"{x:.6f}" for x in r[2:]])
        print(f"\n[OK] {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
