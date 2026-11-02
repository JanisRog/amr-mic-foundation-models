#!/usr/bin/env python3
"""Supplementary isolate table: every isolate analysed, both organisms.

One row per isolate, with its public accession and, per antimicrobial, the log2
MIC the models were trained on, the partition it was assigned to, and the
held-out prediction of the reported model. Read from the final fold files
(fold < 0 development, fold >= 0 held out) and the saved held-out predictions,
never recomputed.

Before writing, the script requires that each prediction file covers exactly the
held-out isolates of its fold file and carries the same measured MIC. A table
that disagreed with the partition it documents would be worse than none.

  M. tuberculosis   reported model: Evo 2, wider target set + TB-Profiler call
  P. aeruginosa     reported model: Evo 2 embeddings + ARDaP call (arm CR); for levofloxacin and aztreonam,
                    which ARDaP does not call, + AMRFinderPlus presence (arm D). Read from the tree
                    stage_reported_pa.py writes (--pa-folds/--pa-preds) and the reported workbook.
"""

from __future__ import annotations

import argparse
import csv
import os
import pathlib
import re
import sys

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

MTB = {  # antimicrobial: (fold dir, log2 MIC above which resistant)
    "isoniazid": ("cv_isoniazid_final", -3.0), "rifampicin": ("cv_rifampicin", -1.0),
    "ethambutol": ("cv_ethambutol", 2.0), "ethionamide": ("cv_ethionamide", 2.0),
    "kanamycin": ("cv_kanamycin", 2.0), "amikacin": ("cv_amikacin", 0.0),
    "levofloxacin": ("cv_levofloxacin", 0.0), "moxifloxacin": ("cv_moxifloxacin", 0.0),
    "bedaquiline": ("cv_bedaquiline", -2.0), "delamanid": ("cv_delamanid", -3.0),
}
PA = {  # antimicrobial: (log2 MIC above which resistant, MIC column in the metadata)
    "ciprofloxacin": (-1.0, "Ciprofloxacin"), "levofloxacin": (1.0, "Levofloxacin"),
    "tobramycin": (1.0, "Tobramycin"), "amikacin": (4.0, "Amikacin"),
    "imipenem": (2.0, "Imipenem"), "meropenem": (3.0, "Meropenem"),
    "ceftazidime": (3.0, "Ceftazidime"), "cefepime": (3.0, "Cefepime"),
    "aztreonam": (4.0, "Aztreonam"), "piperacillin-tazobactam": (4.0, "Piperacillin-tazobactam"),
}
PA_BIOPROJECT = {
    "Stanton et al. (2022)": "PRJNA288601", "Khaledi et al. (2021)": "PRJNA526797",
    "Alfred Hospital (current report)": "PRJNA1220180", "Cortes-Lara et al. (2021)": "PRJEB40140",
    "Barrio-Tofiño et al. (2019)": "PRJEB31047", "Sastre-Femenia et al. (2023)": "PRJEB61879",
    "Barrio-Tofiño et al. (2017)": "PRJEB21341",
}
SUBJ = re.compile(r"^(site\.[^.]+\.subj\.[^.]+)\.lab\.")


def rd(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def num(x):
    f = float(x)
    return int(f) if f == int(f) else f


def load_drug(folds_path, preds_path, label, problems):
    folds = rd(folds_path)
    if len({r["sample_id"] for r in folds}) != len(folds):
        problems.append(f"{label}: an isolate appears twice in {folds_path}")
    mic = {r["sample_id"]: float(r["log2_mic"]) for r in folds}
    part = {r["sample_id"]: ("held-out" if int(r["fold"]) >= 0 else "development") for r in folds}
    pred = {r["sample_id"]: (float(r["log2_mic"]), float(r["pred"]))
            for r in rd(preds_path) if r["aggregation"] == "median"}
    held = {s for s, p in part.items() if p == "held-out"}
    if set(pred) != held:
        problems.append(f"{label}: {len(pred)} predictions, {len(held)} held-out isolates, "
                        f"{len(set(pred) ^ held)} differ")
    elif any(abs(pred[s][0] - mic[s]) > 1e-9 for s in held):
        problems.append(f"{label}: prediction file and fold file disagree on a measured MIC")
    return mic, part, {s: v[1] for s, v in pred.items()}


def dev_folds(path, part, label, problems):
    """Fold of each development isolate in each repeat, as 'f/f/f' in repeat order."""
    by = {}
    for r in rd(path):
        if int(r["fold"]) >= 0:
            by.setdefault(r["sample_id"], {})[int(r["repeat"])] = int(r["fold"])
    dev = {s for s, p in part.items() if p == "development"}
    if set(by) != dev:
        problems.append(f"{label}: {path} assigns folds to {len(by)} isolates, the final fold file "
                        f"has {len(dev)} development isolates, {len(set(by) ^ dev)} differ")
    reps = {tuple(sorted(v)) for v in by.values()}
    if len(reps) != 1:
        problems.append(f"{label}: isolates differ in the repeats they appear in ({sorted(reps)[:3]})")
    return {s: "/".join(str(v[k]) for k in sorted(v)) for s, v in by.items()}


def write_sheet(wb, title, header, rows):
    ws = wb.create_sheet(title)
    ws.append(header)
    for c in ws[1]:
        c.font = Font(bold=True)
    for r in rows:
        ws.append(r)
    ws.freeze_panes = "B2"
    for i, h in enumerate(header, 1):
        ws.column_dimensions[get_column_letter(i)].width = max(12, min(44, len(str(h)) + 2))
    return ws


def counts_rows(drugs, data, bounds, extra=None):
    out = []
    for d in drugs:
        mic, part = data[d][:2]
        row = [d]
        if extra:
            row.append(extra[d])
        for p in ("development", "held-out"):
            ids = [s for s in mic if part[s] == p]
            res = sum(1 for s in ids if mic[s] > bounds[d])
            row += [len(ids), res, len(ids) - res]
        out.append(row[:1 + bool(extra)] + [len(mic)] + row[1 + bool(extra):] + [bounds[d]])
    return out


def main() -> int:
    home = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtb-folds", default=str(home / "mtb_wider_figure_inputs"))
    ap.add_argument("--mtb-preds", default=str(home / "mtb_holdout_20261012/predictions"))
    ap.add_argument("--mtb-ena", default=str(home / "Documents/New project/mtb_download_ena_assemblies/"
                                                    "mtb_metadata_ena_samples.tsv"))
    ap.add_argument("--mtb-devfolds", default="", help="holds cv_<drug>/cv_folds.tsv; adds the fold columns")
    ap.add_argument("--pa-devfolds", default="", help="holds cv_<drug>/cv_folds.tsv; adds the fold columns")
    ap.add_argument("--mtb-lineage", default="", help="sample_lineage.tsv; adds the main lineage column")
    ap.add_argument("--pa-folds", default=str(home / "pa_strict/slim_strict"))
    ap.add_argument("--pa-preds", default=str(home / "pa_strict/predictions_strict"))
    ap.add_argument("--pa-meta", default=str(home / "Documents/New project/_tmp_pseaer_labels_all/"
                                                    "metadata_flat.tsv"))
    ap.add_argument("--pa-workbook", default=str(home / "AMRPredict_pseaer_strict_results.xlsx"))
    ap.add_argument("--out", default=str(home / "AMRPredict_isolates_used.xlsx"))
    a = ap.parse_args()
    problems = []

    # ---- M. tuberculosis
    mtb = {}
    for d, (fdir, _) in MTB.items():
        mtb[d] = load_drug(pathlib.Path(a.mtb_folds) / fdir / "cv_folds_final.tsv",
                           pathlib.Path(a.mtb_preds) / f"{d}_evo2_holdout_predictions_gwfrule1k.tsv",
                           f"MTB/{d}", problems)
        if a.mtb_devfolds:
            mtb[d] += (dev_folds(pathlib.Path(a.mtb_devfolds) / fdir / "cv_folds.tsv", mtb[d][1],
                                 f"MTB/{d}", problems),)
    ena = {r["Sample_IDs"]: r["ENA_SAMPLE"] for r in rd(a.mtb_ena)}
    mtb_ids = sorted({s for d in mtb for s in mtb[d][0]})
    missing = [s for s in mtb_ids if not ena.get(s)]
    if missing:
        problems.append(f"MTB: {len(missing)} isolates without an ENA sample accession")

    # ---- P. aeruginosa: reported model = lowest development cost
    ws = load_workbook(a.pa_workbook, read_only=True)["summary"]
    rows = list(ws.iter_rows(values_only=True))
    ix = {k: i for i, k in enumerate(rows[0])}
    dev = {}
    for r in rows[1:]:
        if r[ix["model"]] in ("evo2", "carbon", "dnabert") and r[ix["split"]] == "development_cv":
            dev.setdefault(r[ix["drug"]], []).append((r[ix["cost"]], r[ix["model"]]))
    reported = {d: "evo2" for d in dev}
    pa = {}
    for d in PA:
        if d not in reported:
            problems.append(f"PA/{d}: no development rows in {a.pa_workbook}")
            continue
        pa[d] = load_drug(pathlib.Path(a.pa_folds) / f"cv_{d}" / "cv_folds_final.tsv",
                          pathlib.Path(a.pa_preds) / f"{d}_{reported[d]}_holdout_predictions.tsv",
                          f"PA/{d}", problems)
        if a.pa_devfolds:
            pa[d] += (dev_folds(pathlib.Path(a.pa_devfolds) / f"cv_{d}" / "cv_folds.tsv", pa[d][1],
                                f"PA/{d}", problems),)
    meta = {r["Accession number"]: r for r in rd(a.pa_meta)}
    pa_ids = sorted({s for d in pa for s in pa[d][0]})
    for s in pa_ids:
        if s not in meta:
            problems.append(f"PA: {s} absent from the metadata")
        elif meta[s]["Study"] not in PA_BIOPROJECT:
            problems.append(f"PA: {s} belongs to an unmapped study {meta[s]['Study']!r}")

    if problems:
        print("[ERROR] inputs do not check out; no table written:")
        for p in problems[:30]:
            print("  ", p)
        return 1

    wb = Workbook()
    readme = wb.active
    readme.title = "README"
    for line in (
        ["Isolates analysed in this study"],
        [],
        ["MTB_isolates", "One row per M. tuberculosis isolate (CRyPTIC collection) with an MIC for at "
                         "least one of the ten antimicrobials."],
        ["PA_isolates", "One row per P. aeruginosa isolate."],
        ["MTB_counts, PA_counts", "Isolates per antimicrobial and partition, with resistant and susceptible counts."],
        [],
        ["Columns per antimicrobial"],
        ["<antimicrobial> log2 MIC", "log2 of the MIC in mg/L, as used for training and evaluation. For "
                                    "P. aeruginosa, censored values are after the rule described in Methods."],
        ["<antimicrobial> set", "development (model fitting and selection) or held-out (scored once). "
                                "Empty: no MIC for this antimicrobial. Partitions were drawn per antimicrobial."],
        ["<antimicrobial> folds", "Development isolates only; see 'How to read the folds column' below."],
        ["<antimicrobial> predicted log2 MIC", "Prediction of the reported model; held-out isolates only."],
        ["<antimicrobial> reported MIC", "P. aeruginosa only: the MIC as reported by the source study (mg/L)."],
        [],
        ["Other columns"],
        ["subject", "MTB: the site and subject part of the CRyPTIC identifier."],
        ["main lineage", "MTB: main lineage assigned by TB-Profiler; 'mixed' where more than one was detected. "
                         "Isolates with a mixed lineage were excluded from the lineage-held-out analysis only."],
        ["resistant above (log2 MIC)", "Counts sheets: an isolate is resistant if its log2 MIC exceeds this value."],
        ["reported model", "PA_counts: the model whose held-out predictions are listed: Evo 2 embeddings with ARDaP's "
                           "call, or, for levofloxacin and aztreonam (no ARDaP call), with AMRFinderPlus presence."],
        [],
        ["How to read the folds column"],
        ["Scheme", "The development set of each antimicrobial was split into five folds, and this was "
                   "repeated three times with a different random split (five-fold cross-validation, "
                   "three repeats). In each repeat an isolate is the validation isolate in exactly one fold "
                   "and a training isolate in the other four."],
        ["Format", "Three numbers separated by slashes: the fold in which the isolate was the validation "
                   "isolate in repeat 1, repeat 2 and repeat 3. Folds are numbered 0 to 4."],
        ["Example", "3/0/4 means: validation isolate in fold 3 of repeat 1, fold 0 of repeat 2 and fold 4 of "
                    "repeat 3; training isolate in every other fold."],
        ["Empty", "The isolate is in the held-out set for this antimicrobial (it was never used for fitting "
                  "or selection and has a prediction instead), or it has no MIC for this antimicrobial."],
        ["Per antimicrobial", "Folds were drawn separately for each antimicrobial, so the same isolate can "
                              "have different folds, or be held out, for different antimicrobials."],
    ):
        readme.append(line)
    readme["A1"].font = Font(bold=True, size=13)
    for cell in ("A7", "A14", "A20"):
        readme[cell].font = Font(bold=True)
    readme.column_dimensions["A"].width = 36
    readme.column_dimensions["B"].width = 110

    def per_drug(data, drugs, s, raw=None):
        out = []
        for d in drugs:
            mic, part, pred = data[d][:3]
            if raw is not None:
                out.append(raw(d) if s in mic else None)
            out += [num(mic[s]) if s in mic else None, part.get(s)]
            if len(data[d]) > 3:
                out.append(data[d][3].get(s))
            out.append(num(pred[s]) if s in pred else None)
        return out

    lineage = {}
    if a.mtb_lineage:
        lineage = {r["sample_id"]: r["main_lineage"].replace("lineage", "lineage ") for r in rd(a.mtb_lineage)}
        gone = [s for s in mtb_ids if s not in lineage]
        if gone:
            print(f"[ERROR] {len(gone)} M. tuberculosis isolates absent from {a.mtb_lineage}")
            return 1
    head = ["isolate", "ENA sample accession", "subject"] + (["main lineage"] if lineage else []) \
        + ["antimicrobials with MIC"]
    for d in MTB:
        head += [f"{d} log2 MIC", f"{d} set"] + ([f"{d} folds"] if a.mtb_devfolds else []) \
            + [f"{d} predicted log2 MIC"]
    write_sheet(wb, "MTB_isolates", head, (
        [s, ena[s], (m.group(1) if (m := SUBJ.match(s)) else None)] + ([lineage[s]] if lineage else [])
        + [sum(1 for d in mtb if s in mtb[d][0])]
        + per_drug(mtb, MTB, s) for s in mtb_ids))

    chead = ["isolates", "development", "development resistant", "development susceptible",
             "held-out", "held-out resistant", "held-out susceptible", "resistant above (log2 MIC)"]
    write_sheet(wb, "MTB_counts", ["antimicrobial"] + chead,
                counts_rows(MTB, mtb, {d: v[1] for d, v in MTB.items()}))

    head = ["run accession", "study", "BioProject", "antimicrobials with MIC"]
    for d in PA:
        head += [f"{d} reported MIC", f"{d} log2 MIC", f"{d} set"] \
            + ([f"{d} folds"] if a.pa_devfolds else []) + [f"{d} predicted log2 MIC"]
    write_sheet(wb, "PA_isolates", head, (
        [s, meta[s]["Study"], PA_BIOPROJECT[meta[s]["Study"]], sum(1 for d in pa if s in pa[d][0])]
        + per_drug(pa, PA, s, raw=lambda d, s=s: meta[s].get(f"{PA[d][1]}__MIC") or None) for s in pa_ids))
    write_sheet(wb, "PA_counts", ["antimicrobial", "reported model"] + chead,
                counts_rows(PA, pa, {d: v[0] for d, v in PA.items()},
                           extra={d: "Evo 2 + AMRFinderPlus presence" if d in ("levofloxacin", "aztreonam")
                                  else "Evo 2 + ARDaP call" for d in reported}))

    wb.save(a.out)
    subjects = {m.group(1) for s in mtb_ids if (m := SUBJ.match(s))}
    print(f"[OK] {a.out}")
    print(f"  M. tuberculosis: {len(mtb_ids):,} isolates, {len(subjects):,} subjects, "
          f"{len(set(ena[s] for s in mtb_ids)):,} ENA samples")
    print(f"  P. aeruginosa:   {len(pa_ids):,} isolates; reported model "
          + ", ".join(f"{d} {m}" for d, m in reported.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
