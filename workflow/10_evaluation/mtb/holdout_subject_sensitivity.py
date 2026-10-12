#!/usr/bin/env python3
"""Do held-out results depend on isolates whose subject is also in development?

The CRyPTIC collection is not one isolate per subject: 10,228 isolates come from
9,991 subjects, and partitions were drawn by isolate. So a held-out isolate can
share a subject with a development isolate -- 2.5% of held-out measurements do.
If those isolates were easier (a near-identical genome from the same patient in
training), the held-out scores would flatter the model.

This rescoring removes them and compares. Nothing is refit: it reads the saved
per-isolate held-out predictions of the reported model and the final fold file
(fold < 0 development, fold >= 0 held out), and uses the same metric function as
the figures. The full-set cost must reproduce the reported value, which is checked
before anything else is printed.

A subject is read from the sample identifier, site.<s>.subj.<id>.lab.<l>.iso.<n>.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from plot_mic_scatter import metrics  # noqa: E402

SUBJ = re.compile(r"^(site\.[^.]+\.subj\.[^.]+)\.lab\.")
# drug -> (S/R boundary, fold directory, predictions file, reported held-out cost)
DRUGS = {
    "isoniazid": (-3.0, "cv_isoniazid_final", "isoniazid_evo2_holdout_predictions_tbrule.tsv", None),
    "rifampicin": (-1.0, "cv_rifampicin", "rifampicin_evo2_holdout_predictions.tsv", 0.0760),
    "ethambutol": (2.0, "cv_ethambutol", "ethambutol_evo2_holdout_predictions.tsv", 0.2738),
    "ethionamide": (2.0, "cv_ethionamide", "ethionamide_evo2_holdout_predictions.tsv", 0.4079),
    "kanamycin": (2.0, "cv_kanamycin", "kanamycin_evo2_holdout_predictions.tsv", 0.3122),
    "amikacin": (0.0, "cv_amikacin", "amikacin_evo2_holdout_predictions.tsv", 0.2038),
    "levofloxacin": (0.0, "cv_levofloxacin", "levofloxacin_evo2_holdout_predictions.tsv", 0.1938),
    "moxifloxacin": (0.0, "cv_moxifloxacin", "moxifloxacin_evo2_holdout_predictions.tsv", 0.1312),
    "bedaquiline": (-2.0, "cv_bedaquiline", "bedaquiline_evo2_holdout_predictions.tsv", None),
    "delamanid": (-3.0, "cv_delamanid", "delamanid_evo2_holdout_predictions.tsv", None),
}
KEYS = ("cost", "vMER", "MER", "SR_f1", "OT", "bal_OT", "MAE")


def score(rows, boundary):
    y = np.array([float(r["log2_mic"]) for r in rows])
    p = np.array([float(r["pred"]) for r in rows])
    m = metrics(y, p, boundary)
    m["cost"] = 2 * m["vMER"] + m["MER"]
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions-dir", required=True)
    ap.add_argument("--folds-root", required=True, help="holds cv_<drug>/cv_folds_final.tsv")
    ap.add_argument("--aggregation", default="median")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    pdir, froot = pathlib.Path(a.predictions_dir), pathlib.Path(a.folds_root)

    out_rows, problems = [], []
    for drug, (bound, fdir, pfile, reported) in DRUGS.items():
        preds = [r for r in csv.DictReader(open(pdir / pfile), delimiter="\t")
                 if r["aggregation"] == a.aggregation]
        folds = list(csv.DictReader(open(froot / fdir / "cv_folds_final.tsv"), delimiter="\t"))
        dev = {r["sample_id"] for r in folds if int(r["fold"]) < 0}
        hold = {r["sample_id"] for r in folds if int(r["fold"]) >= 0}
        ids = [r["sample_id"] for r in preds]
        if set(ids) != hold or len(ids) != len(set(ids)):
            problems.append(f"{drug}: predictions cover {len(set(ids))} isolates, "
                            f"fold file holds out {len(hold)}")
            continue
        dev_subj = {m.group(1) for s in dev if (m := SUBJ.match(s))}
        shared = [r for r in preds if (m := SUBJ.match(r["sample_id"])) and m.group(1) in dev_subj]
        clean = [r for r in preds if r not in shared]
        full, cl = score(preds, bound), score(clean, bound)
        if reported is not None and abs(full["cost"] - reported) > 5e-4:
            problems.append(f"{drug}: recomputed held-out cost {full['cost']:.4f} "
                            f"does not reproduce the reported {reported:.4f}")
        sh = score(shared, bound) if len(shared) >= 5 else None
        out_rows.append((drug, len(preds), len(shared), full, cl, sh))

    if problems:
        print("[ERROR] not reporting a sensitivity on inputs that do not check out:")
        for p in problems:
            print("  ", p)
        return 1

    print(f"Held-out scores with and without isolates whose subject is also in "
          f"development ({a.aggregation} aggregation, reported model)\n")
    print(f"{'drug':<14}{'n':>6}{'shared':>8}   "
          + "".join(f"{k + ' all':>12}{'excl.':>9}{'Δ':>9}" for k in ("cost", "SR_f1", "bal_OT", "MAE")))
    for drug, n, ns, full, cl, _ in out_rows:
        print(f"{drug:<14}{n:>6}{ns:>8}   "
              + "".join(f"{full[k]:>12.4f}{cl[k]:>9.4f}{cl[k] - full[k]:>+9.4f}"
                        for k in ("cost", "SR_f1", "bal_OT", "MAE")))
    for k in ("cost", "SR_f1", "bal_OT", "MAE"):
        d = [cl[k] - full[k] for _, _, _, full, cl, _ in out_rows]
        print(f"  Δ{k:<7} median {np.median(d):+.4f}   range {min(d):+.4f} to {max(d):+.4f}")

    print("\nThe shared-subject isolates on their own (small n; one-tier accuracy and MAE only):")
    print(f"{'drug':<14}{'n':>6}{'OT shared':>11}{'OT others':>11}{'MAE shared':>12}{'MAE others':>12}")
    for drug, n, ns, full, cl, sh in out_rows:
        if sh:
            print(f"{drug:<14}{ns:>6}{sh['OT']:>11.3f}{cl['OT']:>11.3f}{sh['MAE']:>12.3f}{cl['MAE']:>12.3f}")

    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["drug", "n_holdout", "n_shared_subject"]
                       + [f"{k}_{s}" for k in KEYS for s in ("all", "excluding_shared")])
            for drug, n, ns, full, cl, _ in out_rows:
                w.writerow([drug, n, ns] + [f"{v[k]:.6f}" for k in KEYS for v in (full, cl)])
        print(f"\n[OK] {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
