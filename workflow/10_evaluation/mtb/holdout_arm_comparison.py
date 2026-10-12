#!/usr/bin/env python3
"""Held-out comparison of three M. tuberculosis models and the catalogue (§27).

  current   embeddings of the original target set + TB-Profiler call (scored Sept.)
  wider     embeddings of the wider target set + TB-Profiler call    (the adopted model)
  snp       catalogued mutation columns + TB-Profiler call, no embeddings
  TBP       TB-Profiler's call alone, from the baseline rows

Everything is recomputed from per-isolate held-out predictions with one metric
function, so the three models are compared on identical isolates and identical
definitions. Before reporting, the script requires that

  * the three prediction files of a drug cover exactly the same isolates, and
    those are the held-out isolates of the final fold file;
  * the current model's recomputed cost reproduces its reported value;
  * the recomputed metrics of the new runs reproduce what the cluster wrote to
    their folds.tsv.

Paired tests are across the drugs with a verdict (bedaquiline and delamanid shown,
not tested), Wilcoxon signed-rank, as fixed in §27.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics
import sys

import numpy as np

DRUGS = {  # drug: (boundary, fold dir, summaries dir, reported cost of current model)
    "isoniazid": (-3.0, "cv_isoniazid_final", "cv_final_summaries", 0.0996),
    "rifampicin": (-1.0, "cv_rifampicin", "cv_rif_summaries", 0.0760),
    "ethambutol": (2.0, "cv_ethambutol", "cv_emb_summaries", 0.2738),
    "ethionamide": (2.0, "cv_ethionamide", "cv_eth_summaries", 0.4079),
    "kanamycin": (2.0, "cv_kanamycin", "cv_kan_summaries", 0.3122),
    "amikacin": (0.0, "cv_amikacin", "cv_ami_summaries", 0.2038),
    "levofloxacin": (0.0, "cv_levofloxacin", "cv_lev_summaries", 0.1938),
    "moxifloxacin": (0.0, "cv_moxifloxacin", "cv_mxf_summaries", 0.1312),
    "bedaquiline": (-2.0, "cv_bedaquiline", "cv_bdq_summaries", None),
    "delamanid": (-3.0, "cv_delamanid", "cv_dlm_summaries", None),
}
NO_VERDICT = {"bedaquiline", "delamanid"}
KEYS = ("cost", "vMER", "MER", "SR_f1", "exact", "OT", "bal_OT", "bal_OT_common", "MAE")
LOWER = {"cost", "vMER", "MER", "MAE"}
MIN_CLASS_N = 10


def read_preds(path, aggregation):
    rows = [r for r in csv.DictReader(open(path), delimiter="\t") if r["aggregation"] == aggregation]
    return {r["sample_id"]: (float(r["log2_mic"]), float(r["pred"])) for r in rows}


def score(pairs, boundary):
    y = np.array([p[0] for p in pairs]); p = np.array([p[1] for p in pairs])
    err = np.abs(y - p)
    tr, pr = y > boundary, p > boundary
    tp, fp, fn = int((tr & pr).sum()), int((~tr & pr).sum()), int((tr & ~pr).sum())
    n_r, n_s = int(tr.sum()), int((~tr).sum())
    vmer = fn / n_r if n_r else float("nan")
    mer = fp / n_s if n_s else float("nan")
    classes = sorted(set(y.tolist()))
    rec = {c: float(np.mean(err[y == c] <= 1)) for c in classes}
    common = [c for c in classes if int((y == c).sum()) >= MIN_CLASS_N]
    return {"cost": 2 * vmer + mer, "vMER": vmer, "MER": mer,
            "SR_f1": 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else float("nan"),
            "exact": float(np.mean(err == 0)), "OT": float(np.mean(err <= 1)),
            "bal_OT": float(np.mean(list(rec.values()))),
            "bal_OT_common": float(np.mean([rec[c] for c in common])),
            "MAE": float(np.mean(err)), "n": len(y), "n_res": n_r}


def cluster_row(run_dir, aggregation):
    f = run_dir / "folds.tsv"
    if not f.exists():
        return None
    rows = [r for r in csv.DictReader(open(f), delimiter="\t") if r["aggregation"] == aggregation]
    return rows[0] if rows else None


def wilcoxon(d):
    try:
        from scipy.stats import wilcoxon as w
        return float(w(d).pvalue) if any(x != 0 for x in d) else float("nan")
    except Exception:
        return float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True, help="download holding predictions/ and cv_*/runs/*_HOLDOUT")
    ap.add_argument("--current-preds", required=True, help="dir with <drug>_evo2_holdout_predictions*.tsv")
    ap.add_argument("--folds-root", required=True, help="holds cv_<drug>/cv_folds_final.tsv")
    ap.add_argument("--summaries", required=True)
    ap.add_argument("--model", default="evo2")
    ap.add_argument("--aggregation", default="median")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    new, cur, froot, summ = map(pathlib.Path, (a.new, a.current_preds, a.folds_root, a.summaries))

    res, problems = {}, []
    for drug, (bound, fdir, sdir, reported) in DRUGS.items():
        cur_file = cur / (f"{drug}_{a.model}_holdout_predictions_tbrule.tsv" if drug == "isoniazid"
                          else f"{drug}_{a.model}_holdout_predictions.tsv")
        P = {"current": read_preds(cur_file, a.aggregation),
             "wider": read_preds(new / "predictions" / f"{drug}_{a.model}_holdout_predictions_gwfrule1k.tsv", a.aggregation),
             "snp": read_preds(new / "predictions" / f"{drug}_catalogue_holdout_predictions_snprule1k.tsv", a.aggregation)}
        hold = {r["sample_id"] for r in csv.DictReader(open(froot / fdir / "cv_folds_final.tsv"), delimiter="\t")
                if int(r["fold"]) >= 0}
        for k, v in P.items():
            if set(v) != hold:
                problems.append(f"{drug}/{k}: {len(v)} isolates predicted, {len(hold)} held out, "
                                f"{len(set(v) ^ hold)} differ")
        if any(drug in p for p in problems):
            continue
        ids = sorted(hold)
        truth = {k: [v[i][0] for i in ids] for k, v in P.items()}
        if not (truth["current"] == truth["wider"] == truth["snp"]):
            problems.append(f"{drug}: the three files disagree on the measured MICs")
            continue
        M = {k: score([v[i] for i in ids], bound) for k, v in P.items()}
        if reported is not None and abs(M["current"]["cost"] - reported) > 5e-4:
            problems.append(f"{drug}: current cost {M['current']['cost']:.4f} != reported {reported:.4f}")
        for k, run in (("wider", f"{a.model}_gwasfiltered_gwfrule1k_HOLDOUT"), ("snp", "catalogue_snprule1k_HOLDOUT")):
            row = cluster_row(new / fdir / "runs" / run, a.aggregation)
            if row is None:
                problems.append(f"{drug}/{k}: no cluster folds.tsv"); continue
            for key in ("vMER", "MER", "MAE", "bal_OT", "bal_OT_common", "OT"):
                if key in row and abs(float(row[key]) - M[k][key]) > 5e-4:
                    problems.append(f"{drug}/{k}: {key} recomputed {M[k][key]:.4f} != cluster {float(row[key]):.4f}")
        tbp = None
        for r in csv.DictReader(open(summ / sdir / "baseline_rows.tsv"), delimiter="\t"):
            if r["label"] == "TB-Profiler 6.3.0" and r["split"] == "holdout":
                tbp = {k: float(r[k]) for k in ("cost", "vMER", "MER", "SR_f1")}
        M["tbp"] = tbp
        res[drug] = M

    if problems:
        print("[ERROR] inputs do not check out; nothing is reported:")
        for p in problems:
            print("  ", p)
        return 1

    def table(metric):
        print(f"\n{metric}  ({'lower' if metric in LOWER else 'higher'} is better)")
        print(f"{'drug':<14}{'n':>6}{'current':>10}{'wider':>10}{'SNP+call':>10}{'TB-Profiler':>13}"
              f"{'wider-current':>15}{'wider-SNP':>11}")
        for d, M in res.items():
            t = M["tbp"].get(metric) if M["tbp"] else None
            print(f"{d:<14}{M['current']['n']:>6}{M['current'][metric]:>10.4f}{M['wider'][metric]:>10.4f}"
                  f"{M['snp'][metric]:>10.4f}{(f'{t:.4f}' if t is not None else '—'):>13}"
                  f"{M['wider'][metric] - M['current'][metric]:>+15.4f}"
                  f"{M['wider'][metric] - M['snp'][metric]:>+11.4f}" + ("   (no verdict)" if d in NO_VERDICT else ""))

    print(f"Held-out sets, {a.model} for embedding models, {a.aggregation} aggregation. "
          f"All inputs reproduce their reference values.")
    for k in ("cost", "SR_f1", "bal_OT_common", "MAE"):
        table(k)

    panel = [d for d in res if d not in NO_VERDICT]
    print(f"\nPaired across the {len(panel)} drugs with a verdict (Wilcoxon signed-rank):")
    print(f"{'comparison':<34}{'metric':<16}{'median Δ':>10}{'better':>9}{'P':>9}")
    for label, x, y in (("wider - current", "wider", "current"), ("wider (embeddings) - SNP+call", "wider", "snp"),
                        ("current (embeddings) - SNP+call", "current", "snp")):
        for k in ("cost", "bal_OT_common", "MAE"):
            d = [res[g][x][k] - res[g][y][k] for g in panel]
            better = sum(1 for v in d if v != 0 and ((v < 0) == (k in LOWER)))
            print(f"{label:<34}{k:<16}{statistics.median(d):>+10.4f}{f'{better}/{len(d)}':>9}{wilcoxon(d):>9.4f}")
    print("\nAgainst TB-Profiler on cost (held-out, point estimates):")
    for k in ("current", "wider", "snp"):
        w = sum(1 for g in panel if res[g][k]["cost"] < res[g]["tbp"]["cost"])
        print(f"  {k:<8} lower cost than TB-Profiler in {w} of {len(panel)}")

    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["drug", "model", "n", "n_resistant", *KEYS])
            for d, M in res.items():
                for k in ("current", "wider", "snp"):
                    w.writerow([d, k, M[k]["n"], M[k]["n_res"]] + [f"{M[k][x]:.6f}" for x in KEYS])
                if M["tbp"]:
                    w.writerow([d, "tbprofiler", "", ""] + [f"{M['tbp'][x]:.6f}" if x in M["tbp"] else "" for x in KEYS])
        print(f"\n[OK] {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
