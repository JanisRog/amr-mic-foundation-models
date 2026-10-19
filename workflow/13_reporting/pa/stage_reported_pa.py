#!/usr/bin/env python3
"""Stage the reported P. aeruginosa model's results under the names the figure scripts read.

The reported model is arm CR (embeddings + ARDaP's call) for the eight antimicrobials
ARDaP calls and arm D (embeddings + AMRFinderPlus presence) for levofloxacin and
aztreonam (decision of 2026-10-19, EVALUATION_PROTOCOL_PSEAER.md). The workbook builder
and the figure drivers read one run directory name, one configuration name and one
predictions name per drug, so this copies each drug's chosen arm to the canonical
names (`<model>_reported1k`, `<model>_reported_HOLDOUT`, `final_config_<drug>_<model>_reported.json`,
`<drug>_<model>_holdout_predictions.tsv`) in a fresh directory laid out like the earlier
`pa_strict` download. Every staged file is listed in STAGING_MANIFEST.tsv with its source
path and arm, so a figure built from the staged tree can always be traced back.

Nothing is renamed in place and nothing is mixed silently: a drug-model pair missing any
required file stops the run (use --models to stage fewer models while the Carbon and
DNABERT-2 runs are pending).
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import shutil
import sys

ARDAP_DRUGS = ("amikacin", "cefepime", "ceftazidime", "ciprofloxacin", "imipenem", "meropenem",
               "piperacillin-tazobactam", "tobramycin")
ALL_DRUGS = ("amikacin", "aztreonam", "cefepime", "ceftazidime", "ciprofloxacin", "imipenem",
             "levofloxacin", "meropenem", "piperacillin-tazobactam", "tobramycin")
RUN_FILES = ("folds.tsv", "per_class.tsv", "summary.tsv")


def arm_of(drug: str) -> str:
    return "CR" if drug in ARDAP_DRUGS else "D"


def find(srcs, *candidates):
    """First existing path among candidates (relative to each source, in order)."""
    for s in srcs:
        for c in candidates:
            p = s / c
            if p.exists():
                return p
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", action="append", required=True, type=pathlib.Path,
                    help="download directory (slim/, predictions_strict/, final_configs_strict/, "
                         "shap_strict/, loso_strict/); repeatable, searched in the order given")
    ap.add_argument("--out", required=True, type=pathlib.Path)
    ap.add_argument("--models", default="evo2,carbon,dnabert")
    a = ap.parse_args()
    models = a.models.split(",")
    if a.out.exists() and any(a.out.iterdir()):
        sys.exit(f"[ERROR] {a.out} is not empty; stage into a fresh directory")
    manifest, problems = [], []

    def put(src: pathlib.Path, dst: pathlib.Path, drug, model, arm, kind):
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        manifest.append([drug, model, arm, kind, str(src), str(dst.relative_to(a.out))])

    for d in ALL_DRUGS:
        arm = arm_of(d)
        f = find(a.src, f"slim/cv_{d}/cv_folds_final.tsv", f"slim_strict/cv_{d}/cv_folds_final.tsv")
        if f is None:
            problems.append(f"{d}: cv_folds_final.tsv"); continue
        put(f, a.out / f"slim_strict/cv_{d}/cv_folds_final.tsv", d, "-", "-", "folds_final")
        b = find(a.src, f"slim/cv_{d}/baseline_rows.tsv", f"slim_strict/cv_{d}/baseline_rows.tsv",
                 f"cv_{d}/baseline_rows.tsv")
        if b is not None:
            put(b, a.out / f"cv_{d}/baseline_rows.tsv", d, "-", "-", "baseline_rows")
        for m in models:
            sfx = "_armCR" if arm == "CR" else ""
            need = {
                "dev": find(a.src, f"slim/cv_{d}/runs/{m}_arm{arm}1k/folds.tsv",
                            f"slim_strict/cv_{d}/runs/{m}_arm{arm}1k/folds.tsv"),
                "hold": find(a.src, f"slim/cv_{d}/runs/{m}_arm{arm}_HOLDOUT/folds.tsv",
                             f"slim_strict/cv_{d}/runs/{m}_arm{arm}_HOLDOUT/folds.tsv"),
                "config": find(a.src, f"final_configs_strict/final_config_{d}_{m}_arm{arm}.json",
                               f"scripts/final_configs_strict/final_config_{d}_{m}_arm{arm}.json"),
                "predictions": find(a.src, f"predictions_strict/{d}_{m}_holdout_predictions{sfx}.tsv"),
            }
            missing = [k for k, v in need.items() if v is None]
            if missing:
                problems.append(f"{d}/{m} (arm {arm}): {', '.join(missing)}"); continue
            for kind, runname in (("dev", f"{m}_reported1k"), ("hold", f"{m}_reported_HOLDOUT")):
                srcdir = need[kind].parent
                for fn in RUN_FILES:
                    if (srcdir / fn).exists():
                        put(srcdir / fn, a.out / f"slim_strict/cv_{d}/runs/{runname}/{fn}", d, m, arm, kind)
            put(need["config"], a.out / f"scripts/final_configs_strict/final_config_{d}_{m}_reported.json",
                d, m, arm, "config")
            put(need["predictions"], a.out / f"predictions_strict/{d}_{m}_holdout_predictions.tsv",
                d, m, arm, "predictions")
            for kind in ("shap_long", "target_summary", "predictions"):
                s = find(a.src, f"shap_strict/{d}_{m}{sfx}_{kind}.tsv", f"predictions_strict/{d}_{m}{sfx}_{kind}.tsv")
                if s is not None and not (kind == "predictions" and s.parent.name == "predictions_strict"
                                          and "holdout" in s.name):
                    put(s, a.out / f"predictions_strict/{d}_{m}_{kind}.tsv", d, m, arm, f"shap:{kind}")
            lo = [find(a.src, f"loso_strict/{d}_{m}{sfx}_{k}_predictions.tsv") for k in ("loso", "cv")]
            if all(lo):
                for k, s in zip(("loso", "cv"), lo):
                    put(s, a.out / f"loso_strict/{d}_{m}_{k}_predictions.tsv", d, m, arm, f"loso:{k}")
            for extra in (f"loso_strict/{d}_studies.tsv",):
                s = find(a.src, extra)
                if s is not None and not (a.out / extra).exists():
                    put(s, a.out / extra, d, "-", "-", "loso:studies")
    if problems:
        shutil.rmtree(a.out, ignore_errors=True)
        sys.exit("[ERROR] nothing staged; missing:\n  " + "\n  ".join(problems))
    with open(a.out / "STAGING_MANIFEST.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["drug", "model", "arm", "kind", "source", "staged_as"])
        w.writerows(manifest)
    arms = {d: arm_of(d) for d in ALL_DRUGS}
    print(f"[OK] staged {len(manifest)} files for models {models} into {a.out}")
    print("     arm per drug: " + ", ".join(f"{d.split('-')[0][:6]}={arms[d]}" for d in ALL_DRUGS))
    n_shap = sum(1 for r in manifest if r[3] == "shap:shap_long")
    n_loso = sum(1 for r in manifest if r[3] == "loso:loso")
    print(f"     SHAP files for {n_shap} drug-model pairs, study-held-out for {n_loso}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
