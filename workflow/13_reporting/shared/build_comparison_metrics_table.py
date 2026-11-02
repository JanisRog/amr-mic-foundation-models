#!/usr/bin/env python3
"""Supplementary table: metrics of every comparison reported, both organisms.

Nothing is refit. Each sheet is read from a result file the analysis already
produced, or recomputed from saved per-isolate predictions with the metric
functions the figures use:

  MTB_results        development folds and held-out sets, one row per antimicrobial,
                     model and set: holdout_arm_comparison.tsv (Evo 2 models,
                     TB-Profiler), Carbon and DNABERT-2 scored here from their saved
                     predictions, development folds from the results workbook
  MTB_feature_sets   feature_ablation_FINAL_by_drug.tsv (development folds)
  MTB_lineage        s21_lineage_ladder.tsv and s19_lineage_arms.tsv
  PA_results         summary sheet of the P. aeruginosa results workbook
  PA_study_heldout   per-study metrics from saved study-held-out and
                     cross-validation predictions
  Tests              paired two-sided Wilcoxon signed-rank tests across antimicrobials
"""

from __future__ import annotations

import csv
import os
import pathlib
import statistics
import sys

import numpy as np
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font
from scipy.stats import wilcoxon

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "mtb_tbprofiler"))
import holdout_arm_comparison as H  # noqa: E402
from plot_mic_scatter import metrics  # noqa: E402
sys.path.insert(0, str(HERE))
import fda_criteria as FDA  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
HOLD = HOME / "mtb_holdout_20261012"
ABL = HOME / "mtb_slim_final_20261012/feature_ablation_FINAL_by_drug.tsv"
LIN = HOME / "mtb_verdict_20261012"
PAWB = HOME / "AMRPredict_pseaer_reported_results.xlsx"      # reported model: arm CR (8), arm D (2)
PAWB_D = HOME / "AMRPredict_pseaer_strict_results.xlsx"      # arm D, the original specification
REP_PA = HOME / "pa_reported3_20261019"                     # stage_reported_pa.py output
PRESCALL = HOME / "pa_prescall_20261019"                    # presence + ARDaP call, when downloaded
MTBWB = HOME / "AMRPredict_panel_results_wider.xlsx"
LOSO = REP_PA / "loso_strict"
PAMETA = HOME / "pseaer_metadata_flat.tsv"
# Feature-set and same-loci comparisons on the reported target set (protocol §32-§33):
# slimmed development runs (selected configuration only) and held-out predictions.
NEW_MTB = HOME / "mtb_slim_20261016"
NEW_PA = HOME / "pa_strict_20261020"
NO_VERDICT = {"bedaquiline", "delamanid"}
KEYS = ["cost", "vMER", "MER", "SR_f1", "OT", "bal_OT", "bal_OT_common", "exact", "MAE"]
NAMES = {"cost": "Cost (2 x vMER + MER)", "vMER": "vMER", "MER": "MER", "SR_f1": "F1 (resistant)",
         "OT": "Essential agreement", "bal_OT": "Balanced essential agreement",
         "bal_OT_common": "Balanced essential agreement (classes with >= 10 isolates)",
         "exact": "Exact accuracy", "MAE": "Mean absolute error (doubling dilutions)"}
PA_BOUND = {"ciprofloxacin": -1.0, "levofloxacin": 1.0, "tobramycin": 1.0, "imipenem": 2.0, "meropenem": 3.0,
            "ceftazidime": 3.0, "cefepime": 3.0, "aztreonam": 4.0, "amikacin": 4.0, "piperacillin-tazobactam": 4.0}


def rd(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


MTB_SHAP = HOME / "mtb_shap_gwfrule1k/shap_gwfrule1k"
ISOLATES = HOME / "AMRPredict_isolates_used.xlsx"      # Supplementary Table 1


def coresistance_rows():
    """Figure 5: the AUROC matrices behind a and b, and the numbers quoted in the text for a-d.

    AUROC: how well the log2 MIC of the antimicrobial in the column separates the isolates resistant to the
    antimicrobial in the row from the susceptible ones, among isolates with both MICs and at least ten resistant.
    c, d use the isolates tested for all ten antimicrobials of the organism.
    """
    from make_figure5 import CLASSES, FIRST, SECOND, THIRD, auroc, phenotypes
    wbk = load_workbook(ISOLATES, read_only=True)
    long_rows, summ = [], []
    for k, org in (("mtb", "M. tuberculosis"), ("pseaer", "P. aeruginosa")):
        drugs, mic, res = phenotypes(wbk, k)
        M, cnt = {}, {}
        for da in drugs:
            for db in drugs:
                if da == db:
                    continue
                ids = [x for x in res[da] if x in mic[db]]
                pos = [mic[db][x] for x in ids if res[da][x]]
                neg = [mic[db][x] for x in ids if not res[da][x]]
                if len(pos) >= 10 and neg:
                    M[(da, db)] = auroc(pos, neg)
                    cnt[(da, db)] = (len(pos), len(neg))
                    long_rows.append([org, da, db, round(100 * M[(da, db)], 1), len(pos), len(neg)])
        inside = {(x, y) for cls in CLASSES[k] for x in cls for y in cls if x != y}
        ins = [100 * v for q, v in M.items() if q in inside]
        btw = {q: 100 * v for q, v in M.items() if q not in inside}
        summ += [[org, "AUROC within a drug class, median (%)", round(float(np.median(ins)), 1), "", len(ins),
                  "ordered pairs"],
                 [org, "AUROC between classes, median (%)", round(float(np.median(list(btw.values()))), 1), "",
                  len(btw), "ordered pairs"],
                 [org, "Between-class pairs with AUROC >= 80%", sum(v >= 80 for v in btw.values()), "", len(btw),
                  "ordered pairs"]]
        pairs = [(a_, b_) for i, a_ in enumerate(drugs) for b_ in drugs[i + 1:] if (a_, b_) not in inside
                 and (a_, b_) in M and (b_, a_) in M]
        diffs = [abs(100 * (M[(a_, b_)] - M[(b_, a_)])) for a_, b_ in pairs]
        summ += [[org, "Asymmetry of a between-class pair, median |AUROC(A from B) - AUROC(B from A)| (points)",
                  round(float(np.median(diffs)), 1), "", len(diffs), "unordered pairs, both directions defined"],
                 [org, "... pairs with a difference of 10 points or more", sum(d_ >= 10 for d_ in diffs), "",
                  len(diffs), "unordered pairs"]]
        if k == "mtb":
            later = [d for d in SECOND + THIRD if d in drugs]
            f2l = [100 * M[(a_, b_)] for a_ in later for b_ in FIRST if (a_, b_) in M]
            l2f = [100 * M[(a_, b_)] for a_ in FIRST for b_ in later if (a_, b_) in M]
            summ += [[org, "Resistance to a later drug from a first-line MIC, median AUROC (%)",
                      round(float(np.median(f2l)), 1), "", len(f2l), "ordered pairs"],
                     [org, "Resistance to a first-line drug from a later drug's MIC, median AUROC (%)",
                      round(float(np.median(l2f)), 1), "", len(l2f), "ordered pairs"]]
        # panels c, d: isolates tested for all ten antimicrobials
        ids = [x for x in {x for d in drugs for x in res[d]} if all(x in res[d] for d in drugs)]
        n = len(ids)
        count = {x: sum(res[d][x] for d in drugs) for x in ids}
        summ += [[org, "Isolates tested for all ten antimicrobials (panels c, d)", n, "", "", ""],
                 [org, "... susceptible to all", sum(c == 0 for c in count.values()), "", n, ""],
                 [org, "... resistant to more than half (six or more of ten)", sum(c > len(drugs) / 2 for c in
                                                                                  count.values()), "", n, ""],
                 [org, "... resistant to all ten", sum(c == len(drugs) for c in count.values()), "", n, ""]]
        for d in drugs:
            summ.append([org, f"Resistant to {d} (%)", round(100 * np.mean([res[d][x] for x in ids]), 1), "", n, ""])
        if k == "mtb":
            R = lambda x, ds: any(res[d][x] for d in ds if d in res)  # noqa: E731
            first_r = [x for x in ids if R(x, ("isoniazid", "rifampicin"))]
            first_s = [x for x in ids if not R(x, ("isoniazid", "rifampicin"))]
            sec_r = [x for x in ids if R(x, SECOND)]
            new_r = [x for x in ids if R(x, THIRD)]
            a_, b_ = sum(R(x, SECOND) for x in first_r), sum(R(x, SECOND) for x in first_s)
            summ += [[org, "Resistant to isoniazid or rifampicin", len(first_r), "", n, ""],
                     [org, "... of which resistant to a second-line drug (ethionamide, aminoglycosides, fluoroquinolones)",
                      a_, round(100 * a_ / len(first_r), 1), len(first_r), "numerator, %, denominator"],
                     [org, "Susceptible to both isoniazid and rifampicin, resistant to a second-line drug", b_,
                      round(100 * b_ / len(first_s), 1), len(first_s), "numerator, %, denominator"],
                     [org, "Resistant to a second-line drug, also resistant to isoniazid or rifampicin",
                      sum(R(x, ("isoniazid", "rifampicin")) for x in sec_r),
                      round(100 * np.mean([R(x, ("isoniazid", "rifampicin")) for x in sec_r]), 1), len(sec_r),
                      "numerator, %, denominator"],
                     [org, "Resistant to bedaquiline or delamanid", len(new_r), "", n, ""],
                     [org, "... also resistant to a first- or second-line drug",
                      sum(R(x, FIRST + SECOND) for x in new_r),
                      round(100 * np.mean([R(x, FIRST + SECOND) for x in new_r]), 1), len(new_r),
                      "numerator, %, denominator"]]
    # Extended Data Figure 3: the same AUROC within each collection site (M. tuberculosis) or study (P. aeruginosa)
    from make_extended_figure3 import pair_table
    rng = np.random.default_rng(3)
    centre_rows = []
    for k, org in (("mtb", "M. tuberculosis"), ("pseaer", "P. aeruginosa")):
        pairs, widths = pair_table(wbk, k, rng)
        spreads = []
        for da, db, pooled, per in pairs:
            v = [100 * a_ for a_, *_ in per.values()]
            spreads.append(max(v) - min(v))
            for g, (au, npos, nneg) in per.items():
                centre_rows.append([org, da, db, g, round(100 * au, 1), npos, nneg, round(100 * pooled, 1)])
        summ += [[org, "Between-class pairs with pooled AUROC >= 80% and at least three centres with >= 20 resistant and "
                       ">= 20 susceptible isolates (Extended Data Fig. 3)", len(pairs), "", "", ""],
                 [org, "... spread of the AUROC between centres, median (points)", round(float(np.median(spreads)), 1),
                  "", len(pairs), "best minus worst centre"],
                 [org, "... spread, smallest and largest (points)", f"{min(spreads):.0f} to {max(spreads):.0f}", "",
                  len(pairs), ""],
                 [org, "Width of the 95% bootstrap interval of one centre's AUROC, median (points)",
                  round(float(np.median(widths)), 1), "", len(widths), "200 resamples per centre and pair"]]
    return long_rows, summ, centre_rows


def call_role_rows():
    """How the reported model uses the tool's call (Results, Figure 4): per antimicrobial and for both organisms.

    The tool's call per isolate is read from the sign of the call feature's SHAP value (positive = resistant);
    this reproduces the tools' held-out vMER and MER exactly for M. tuberculosis and within 2.6 points for
    P. aeruginosa (intermediate ARDaP calls). Antimicrobials without a single call column (levofloxacin and
    aztreonam, P. aeruginosa) are not listed.
    """
    import pandas as pd
    from scipy.stats import spearmanr
    sys.path.insert(0, os.environ.get("MTB_SCRIPTS") or str(pathlib.Path(__file__).resolve().parent.parent / "mtb_tbprofiler"))
    from make_panel_figures import load as wb_load, pick
    from plot_metric_bars import PANELS
    books = {"mtb": wb_load(MTBWB)[0], "pseaer": wb_load(PAWB)[0]}
    rows, tot = [], {}
    for k, org, shap_dir, bounds in (
            ("mtb", "M. tuberculosis", MTB_SHAP, {d: v[0] for d, v in H.DRUGS.items()}),
            ("pseaer", "P. aeruginosa", REP_PA / "predictions_strict", PA_BOUND)):
        base = PANELS[k]["baseline"]
        for d, b in bounds.items():
            sh = pd.read_csv(shap_dir / f"{d}_evo2_shap_long.tsv", sep="\t")
            call = sh[sh["target_id"].str.startswith("rule")]
            if call.empty or call["target_id"].nunique() > 1:
                continue
            pr = pd.read_csv(shap_dir / f"{d}_evo2_predictions.tsv", sep="\t")
            if "median" in set(pr["aggregation"]):
                pr = pr[pr["aggregation"] == "median"]
            pr = pr.drop_duplicates("sample_id").set_index("sample_id")
            emb = sh[~sh["target_id"].str.startswith("rule")].groupby("sample_id")["shap"].sum().rename("emb")
            df = (call.set_index("sample_id")[["shap"]].join(pr[["log2_mic", "raw"]], how="inner")
                  .join(emb, how="left"))
            c_r, t_r, p_r = df["shap"] > 0, df["log2_mic"] > b, df["raw"] > b
            within = df.groupby(c_r)["raw"].var().mul(df.groupby(c_r).size() - 1).sum() / (len(df) - 1)
            rho = [spearmanr(g["emb"], g["log2_mic"]).correlation for _, g in df.groupby(c_r)
                   if len(g) >= 20 and g["log2_mic"].nunique() > 1]
            ov = p_r != c_r
            right, wrong = int((ov & (p_r == t_r)).sum()), int((ov & (p_r != t_r)).sum())
            ts = pd.read_csv(shap_dir / f"{d}_evo2_target_summary.tsv", sep="\t")
            call_share = ts.loc[ts["target_id"].str.startswith("rule"), "mean_abs_shap"].sum() / ts["mean_abs_shap"].sum()
            rows.append([org, d, len(df), round(100 * call_share, 1), round(100 * (1 - call_share), 1),
                         round(100 * (1 - within / df["raw"].var()), 1), round(float(np.mean(rho)), 3),
                         right + wrong, round(100 * (right + wrong) / len(df), 1), right,
                         round(100 * right / max(right + wrong, 1), 1),
                         f(pick(books[k], d, base, "development_cv", "SR_f1")),
                         f(pick(books[k], d, base, "holdout", "SR_f1"))])
            t = tot.setdefault(org, [0, 0, 0])
            t[0] += len(df); t[1] += right + wrong; t[2] += right
    foot = [[], ["All antimicrobials listed"]]
    for org, (n, o, rt) in tot.items():
        foot.append([org, "", n, "", "", "", "", o, round(100 * o / n, 1), rt, round(100 * rt / o, 1), "", ""])
    sp = lambda col: spearmanr([r[col] for r in rows], [r[8] for r in rows])  # noqa: E731
    for label, col in (("development folds", 11), ("held-out", 12)):
        rho_, p_ = sp(col)
        foot.append(["Spearman correlation across antimicrobials: override rate against the tool's F1 (" + label + ")",
                     "", len(rows), "", "", "", "", "", "", "", "", "", ""])
        foot[-1][3] = round(float(rho_), 3)
        foot[-1][4] = "P = " + format(p_, ".3f")
    return rows, foot


def f(x):
    try:
        return round(float(x), 4)
    except (TypeError, ValueError):
        return None


def sheet(wb, title, header, rows):
    ws = wb.create_sheet(title)
    ws.append(header)
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="top")
    for r in rows:
        ws.append(r)
    ws.column_dimensions["A"].width = 24
    ws.column_dimensions["B"].width = 46
    ws.freeze_panes = "C2"
    return ws


def wtest(d):
    return float(wilcoxon(d).pvalue) if any(x != 0 for x in d) else float("nan")


def dev_from_folds(path):
    """Development values of a slimmed run: mean over the fifteen folds of the one selected configuration."""
    rows = [r for r in rd(path) if r["aggregation"] == "median"]
    if len(rows) != 15 or len({r["config"] for r in rows}) != 1:
        sys.exit(f"[ERROR] {path}: expected the selected configuration on 15 folds, got {len(rows)} rows")
    mean = lambda k: statistics.mean(float(r[k]) for r in rows)  # noqa: E731
    v = {k: mean(k) for k in KEYS if k != "cost"}
    v["cost"] = statistics.mean(2 * float(r["vMER"]) + float(r["MER"]) for r in rows)
    return v, f(rows[0]["n_features"])


# The MIC measures rest on every held-out isolate, so they use all antimicrobials; the resistance
# call (cost, vMER, MER, F1) is thin for bedaquiline and delamanid (12 and 15 resistant isolates) and
# keeps the eight with a verdict.
MIC_KEYS = {"OT", "bal_OT", "bal_OT_common", "exact", "MAE"}


def add_tests(tests, setting, vals, drugs, a, b, la, lb, keys, all_drugs=None):
    """Paired two-sided Wilcoxon across antimicrobials; difference is first minus second.

    `drugs` are the antimicrobials with a verdict; with `all_drugs` the MIC measures use those instead.
    """
    for k in keys:
        ds = all_drugs if (all_drugs and k in MIC_KEYS) else drugs
        dlt = [float(vals[d][a][k]) - float(vals[d][b][k]) for d in ds]
        better = sum(1 for v in dlt if abs(v) > 5e-4 and ((v < 0) == (k in H.LOWER)))   # ties as in the figures: |difference| <= 0.0005
        tests.append([setting, la, lb, NAMES[k], len(dlt), round(statistics.median(dlt), 4),
                      f"{better}/{len(dlt)}", round(wtest(dlt), 4)])


def main() -> int:
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path("AMRPredict_comparison_metrics.xlsx"))
    wb = Workbook()
    readme = wb.active
    readme.title = "README"
    tests = []

    # ---- M. tuberculosis, held-out
    label = {"wider": "Reported model (Evo 2; target loci + TB-Profiler call)",
             "current": "Evo 2; TB-Profiler-annotated loci only + TB-Profiler call",
             "snp": "Mutation columns + TB-Profiler call (no embeddings)",
             "tbprofiler": "TB-Profiler call alone"}
    held = {}
    rows = []
    for r in rd(HOLD / "holdout_arm_comparison.tsv"):
        held.setdefault(r["drug"], {})[r["model"]] = r
    for d, (bound, fdir, _, _) in H.DRUGS.items():
        for m in ("wider", "current", "snp", "tbprofiler"):
            r = held[d][m]
            rows.append([d, label[m], f(r["n"]), f(r["n_resistant"])] + [f(r.get(k)) for k in KEYS])
        ref = H.read_preds(HOLD / "predictions" / f"{d}_evo2_holdout_predictions_gwfrule1k.tsv", "median")
        for enc, name in (("carbon", "Carbon"), ("dnabert", "DNABERT-2")):
            p = H.read_preds(HOLD / "predictions" / f"{d}_{enc}_holdout_predictions_gwfrule1k.tsv", "median")
            if set(p) != set(ref):
                sys.exit(f"[ERROR] {d}/{enc}: held-out isolates differ from the Evo 2 model's")
            s = H.score([p[i] for i in sorted(p)], bound)
            held[d][enc] = {k: s[k] for k in KEYS}
            rows.append([d, f"{name}; target loci + TB-Profiler call", s["n"], s["n_res"]] + [f(s[k]) for k in KEYS])
    panel = [d for d in H.DRUGS if d not in NO_VERDICT]
    # Embeddings only (three models) and variants only, same loci: held-out from saved
    # predictions, development from the slimmed runs. Same isolates as the reported model.
    newm = [("emb_evo2", "evo2_gwasfiltered_embonly1k", "evo2_holdout_predictions_embonly1k"),
            ("emb_carbon", "carbon_gwasfiltered_embonly1k", "carbon_holdout_predictions_embonly1k"),
            ("emb_dnabert", "dnabert_gwasfiltered_embonly1k", "dnabert_holdout_predictions_embonly1k"),
            ("var", "variants_gwasfiltered_varonly1k", "varonly_holdout_predictions_gwasfiltered")]
    newdev = {}
    for d, (bound, fdir, _, _) in H.DRUGS.items():
        ref = set(H.read_preds(HOLD / "predictions" / f"{d}_evo2_holdout_predictions_gwfrule1k.tsv", "median"))
        for k, run, stem in newm:
            p = H.read_preds(NEW_MTB / "predictions" / f"{d}_{stem}.tsv", "median")
            if set(p) != ref:
                sys.exit(f"[ERROR] {d}/{k}: held-out isolates differ from the reported model's")
            s = H.score([p[i] for i in sorted(p)], bound)
            held[d][k] = {x: s[x] for x in KEYS}
            newdev[(d, k)] = dev_from_folds(NEW_MTB / "slim_arms" / fdir / "runs" / run / "folds.tsv")
    for a in ("emb_evo2", "emb_carbon", "emb_dnabert"):
        add_tests(tests, "M. tuberculosis, held-out", held, panel, a, "var",
                  {"emb_evo2": "Evo 2; embeddings only", "emb_carbon": "Carbon; embeddings only",
                   "emb_dnabert": "DNABERT-2; embeddings only"}[a],
                  "Variants only, same loci (no embeddings)", ("OT", "bal_OT_common", "MAE", "SR_f1"),
                  all_drugs=list(H.DRUGS))
    # the three foundation models against each other, embeddings alone (Fig. 1c)
    for b_, lb_ in (("emb_carbon", "Carbon; embeddings only"), ("emb_dnabert", "DNABERT-2; embeddings only")):
        add_tests(tests, "M. tuberculosis, held-out", held, panel, "emb_evo2", b_, "Evo 2; embeddings only", lb_,
                  ("OT", "bal_OT_common", "MAE", "SR_f1"), all_drugs=list(H.DRUGS))
    add_tests(tests, "M. tuberculosis, held-out", held, panel, "emb_carbon", "emb_dnabert", "Carbon; embeddings only",
              "DNABERT-2; embeddings only", ("OT", "bal_OT_common", "MAE", "SR_f1"), all_drugs=list(H.DRUGS))
    # The model restricted to TB-Profiler-annotated loci is a design alternative
    # (Supplementary Table 4); its tests are no longer listed here.
    for a, b in (("wider", "snp"), ("wider", "carbon"), ("wider", "dnabert")):
        for k in ("cost", "bal_OT_common", "MAE"):
            dlt = [float(held[d][a][k]) - float(held[d][b][k]) for d in (list(H.DRUGS) if k in MIC_KEYS else panel)]
            better = sum(1 for v in dlt if abs(v) > 5e-4 and ((v < 0) == (k in H.LOWER)))   # ties as in the figures: |difference| <= 0.0005
            tests.append(["M. tuberculosis, held-out", label["wider"],
                          label.get(b, {"carbon": "Carbon; target loci + TB-Profiler call",
                                        "dnabert": "DNABERT-2; target loci + TB-Profiler call"}.get(b)),
                          NAMES[k], len(dlt), round(statistics.median(dlt), 4), f"{better}/{len(dlt)}",
                          round(wtest(dlt), 4)])

    # ---- M. tuberculosis, development folds and held-out sets in one sheet, laid out
    # like PA_results: one row per antimicrobial, model and set.
    msrc = list(load_workbook(MTBWB, read_only=True)["summary"].iter_rows(values_only=True))
    mx = {k: i for i, k in enumerate(msrc[0])}
    wbkey = {"evo2": "wider", "carbon": "carbon", "dnabert": "dnabert", "TB-Profiler 6.3.0": "tbprofiler"}
    devm, feat = {}, {}
    for r in msrc[1:]:
        d, m, split = r[mx["drug"]], r[mx["model"]], r[mx["split"]]
        if m not in wbkey:
            continue
        k = wbkey[m]
        if split == "development_cv":
            devm[(d, k)] = [f(r[mx[x]]) for x in KEYS]
            feat[(d, k)] = f(r[mx["n_features"]])
        elif k != "tbprofiler" and abs(f(r[mx["cost"]]) - float(held[d][k]["cost"])) > 5e-4:
            # The workbook and the held-out comparison must describe the same models.
            sys.exit(f"[ERROR] {d}/{m}: workbook held-out cost {f(r[mx['cost']])} differs from "
                     f"the held-out comparison's {float(held[d][k]['cost']):.4f}")
    # The two comparison models were selected in the feature-set runs; their
    # development rows come from that table, which carries six of the nine measures.
    for r in rd(ABL):
        k = {"emb_rule": "current", "snp_rule": "snp"}.get(r["arm"])
        if k:
            devm[(r["drug"], k)] = [f(r.get(x)) for x in KEYS]
            feat[(r["drug"], k)] = f(r["n_features"])
    for (d, k), (v, nf) in newdev.items():
        devm[(d, k)] = [f(v[x]) for x in KEYS]
        feat[(d, k)] = nf
    mlabel = [("wider", "Evo 2 (reported)"), ("carbon", "Carbon"), ("dnabert", "DNABERT-2"),
              ("emb_evo2", "Evo 2; embeddings only"), ("emb_carbon", "Carbon; embeddings only"),
              ("emb_dnabert", "DNABERT-2; embeddings only"),
              ("var", "Variants only, same loci (no embeddings)"),
              ("snp", "Mutation columns + TB-Profiler call (no embeddings)"),
              ("tbprofiler", "TB-Profiler")]
    mrow = []
    for d in H.DRUGS:
        for k, name in mlabel:
            if (d, k) in devm:
                mrow.append([d, name, "development folds", feat.get((d, k))] + devm[(d, k)])
            mrow.append([d, name, "held-out", feat.get((d, k))] + [f(held[d][k].get(x)) for x in KEYS])
    sheet(wb, "MTB_results", ["Antimicrobial", "Model", "Set", "Features"] + [NAMES[k] for k in KEYS], mrow)

    # ---- M. tuberculosis, development feature sets on the reported target loci (Evo 2),
    # from the slimmed runs; the reported row must reproduce the results workbook.
    fsets = [("emb_rule", "Embeddings + TB-Profiler call (reported)", "evo2_gwasfiltered_gwfrule1k"),
             ("emb_snp_rule", "Embeddings + mutation columns + TB-Profiler call", "evo2_gwasfiltered_arm76"),
             ("emb_snp", "Embeddings + mutation columns", "evo2_gwasfiltered_gwfsnp1k"),
             ("emb", "Embeddings only", "evo2_gwasfiltered_embonly1k"),
             ("var", "Variants only, same loci", "variants_gwasfiltered_varonly1k"),
             ("snp_rule", "Mutation columns + TB-Profiler call", "catalogue_snprule1k"),
             ("snp", "Mutation columns only", "catalogue_snponly1k")]
    fk = ("cost", "vMER", "MER", "SR_f1", "OT", "bal_OT_common", "MAE")
    by, frows = {}, []
    for d, (_, fdir, _, _) in H.DRUGS.items():
        for arm, name, run in fsets:
            v, nf = dev_from_folds(NEW_MTB / "slim_arms" / fdir / "runs" / run / "folds.tsv")
            by.setdefault(d, {})[arm] = v
            frows.append([d, name, nf] + [f(v[k]) for k in fk])
        if abs(by[d]["emb_rule"]["cost"] - devm[(d, "wider")][KEYS.index("cost")]) > 5e-4:
            sys.exit(f"[ERROR] {d}: reported development cost from the slimmed run differs from the workbook")
        tb = devm.get((d, "tbprofiler"))
        if tb:
            frows.append([d, "TB-Profiler call alone (no model)", None] + [tb[KEYS.index(k)] for k in fk])
    sheet(wb, "MTB_feature_sets", ["Antimicrobial", "Feature set", "Features"] + [NAMES[k] for k in fk], frows)
    lab = {arm: name for arm, name, _ in fsets}
    for arm in ("emb_snp_rule", "emb_snp", "emb", "var", "snp_rule", "snp"):
        add_tests(tests, "M. tuberculosis, development folds", by, panel, "emb_rule", arm,
                  lab["emb_rule"], lab[arm], ("cost", "bal_OT_common", "MAE"), all_drugs=list(H.DRUGS))

    # ---- M. tuberculosis, lineage
    lrows = []
    for r in rd(LIN / "s21_lineage_ladder.tsv"):
        lrows.append([r["drug"].replace("_final", ""), "Reported-type model (embeddings + call) against mutation columns + call",
                      NAMES.get(r["metric"], r["metric"]), f(r["emb_sizematched"]), f(r["emb_lineage"]),
                      f(r["snp_sizematched"]), f(r["snp_lineage"])])
    for r in rd(LIN / "s19_lineage_arms.tsv"):
        lrows.append([r["drug"].replace("_final", ""), "Target loci against TB-Profiler-annotated loci only",
                      NAMES.get(r["metric"], r["metric"]), f(r["rand_gwfrule"]), f(r["lin_gwfrule"]),
                      f(r["rand_baseline"]), f(r["lin_baseline"])])
    sheet(wb, "MTB_lineage", ["Antimicrobial", "Comparison", "Measure", "First model, size-matched random folds",
                              "First model, lineage folds", "Second model, size-matched random folds",
                              "Second model, lineage folds"], lrows)

    # ---- P. aeruginosa, development and held-out
    ws = load_workbook(PAWB, read_only=True)["summary"]
    src = list(ws.iter_rows(values_only=True))
    ix = {k: i for i, k in enumerate(src[0])}
    dev = {}
    for r in src[1:]:
        if r[ix["model"]] in ("evo2", "carbon", "dnabert") and r[ix["split"]] == "development_cv":
            dev.setdefault(r[ix["drug"]], []).append((float(r[ix["cost"]]), r[ix["model"]]))
    reported = {d: "evo2" for d in dev}
    pname = {"evo2": "Evo 2", "carbon": "Carbon", "dnabert": "DNABERT-2",
             "ARDaP (I=R)": "ARDaP (intermediate counted as resistant)",
             "ARDaP (I=S)": "ARDaP (intermediate counted as susceptible)"}
    prow = []
    for r in src[1:]:
        m, d = r[ix["model"]], r[ix["drug"]]
        if m not in pname or f(r[ix["cost"]]) is None:
            continue
        prow.append([d, pname[m] + (" (reported)" if reported.get(d) == m else ""),
                     "held-out" if r[ix["split"]] == "holdout" else "development folds", f(r[ix["n_features"]])]
                    + [f(r[ix[k]]) if k in ix else None for k in KEYS])
    srcd = list(load_workbook(PAWB_D, read_only=True)["summary"].iter_rows(values_only=True))
    ixd = {k: i for i, k in enumerate(srcd[0])}
    for r in srcd[1:]:
        m, d = r[ixd["model"]], r[ixd["drug"]]
        if m not in ("evo2", "carbon", "dnabert") or d in ("aztreonam", "levofloxacin") or f(r[ixd["cost"]]) is None:
            continue
        prow.append([d, pname[m] + "; embeddings + AMRFinderPlus presence (original specification)",
                     "held-out" if r[ixd["split"]] == "holdout" else "development folds", f(r[ixd["n_features"]])]
                    + [f(r[ixd[k]]) if k in ixd else None for k in KEYS])
    # Embeddings only (arm C), variants only and presence only: held-out from saved
    # predictions, development from the slimmed runs. The reported model's held-out
    # values recomputed here must reproduce the workbook's.
    pnew = [("Evo 2; embeddings only", "evo2_armC1k", "evo2_holdout_predictions_armC"),
            ("Carbon; embeddings only", "carbon_armC1k", "carbon_holdout_predictions_armC"),
            ("DNABERT-2; embeddings only", "dnabert_armC1k", "dnabert_holdout_predictions_armC"),
            ("Variants only, same loci + acquired-gene presence (no embeddings)", "varonly",
             "varonly_holdout_predictions"),
            ("AMRFinderPlus presence only (no embeddings)", "base_presence", "presence_holdout_predictions"),
            ("Evo 2; chromosomal embeddings only", "evo2_armA1k", "evo2_holdout_predictions_armA"),
            ("Carbon; chromosomal embeddings only", "carbon_armA1k", "carbon_holdout_predictions_armA"),
            ("DNABERT-2; chromosomal embeddings only", "dnabert_armA1k", "dnabert_holdout_predictions_armA"),
            ("Chromosomal variants only (no embeddings)", "varchromo", "varchromo_holdout_predictions"),
            # exploratory, the 8 antimicrobials ARDaP covers
            ("Evo 2; embeddings + AMRFinderPlus presence + ARDaP call (exploratory)", "evo2_armDR1k",
             "evo2_holdout_predictions_armDR")]
    if (PRESCALL / "predictions_strict").is_dir():
        pnew.append(("AMRFinderPlus presence + ARDaP call (no embeddings)", "base_prescall", "prescall_holdout_predictions"))
    else:
        print("  [INFO] presence + ARDaP call not downloaded yet; its rows and tests are left out")
    pheld = {}
    for d, bound in PA_BOUND.items():
        ref = H.read_preds(REP_PA / "predictions_strict" / f"{d}_evo2_holdout_predictions.tsv", "median")
        rep = H.score([ref[i] for i in sorted(ref)], bound)
        wbrow = [r for r in prow if r[0] == d and r[1] == "Evo 2 (reported)" and r[2] == "held-out"][0]
        if abs(rep["OT"] - wbrow[4 + KEYS.index("OT")]) > 5e-4:
            sys.exit(f"[ERROR] {d}: reported held-out essential agreement differs from the workbook")
        pheld[d] = {"reported": rep}
        for name, run, stem in pnew:
            if "ARDaP" in name and d in ("aztreonam", "levofloxacin"):
                continue                      # no ARDaP call for these two
            base = PRESCALL if run == "base_prescall" else NEW_PA
            p = H.read_preds(base / "predictions_strict" / f"{d}_{stem}.tsv", "median")
            if set(p) != set(ref):
                sys.exit(f"[ERROR] {d}/{name}: held-out isolates differ from the reported model's")
            s = H.score([p[i] for i in sorted(p)], bound)
            pheld[d][name] = s
            v, nf = dev_from_folds(base / "slim" / f"cv_{d}" / "runs" / run / "folds.tsv")
            prow.append([d, name, "development folds", nf] + [f(v[k]) for k in KEYS])
            prow.append([d, name, "held-out", nf] + [f(s[k]) for k in KEYS])
    pd_all = list(PA_BOUND)
    for a in ("Evo 2; embeddings only", "Carbon; embeddings only", "DNABERT-2; embeddings only"):
        add_tests(tests, "P. aeruginosa, held-out", pheld, pd_all, a, pnew[3][0], a, pnew[3][0],
                  ("OT", "bal_OT_common", "MAE", "SR_f1"))
    for b_ in ("Carbon; embeddings only", "DNABERT-2; embeddings only"):      # the models against each other (Fig. 1c)
        add_tests(tests, "P. aeruginosa, held-out", pheld, pd_all, "Evo 2; embeddings only", b_,
                  "Evo 2; embeddings only", b_, ("OT", "bal_OT_common", "MAE", "SR_f1"))
    add_tests(tests, "P. aeruginosa, held-out", pheld, pd_all, "Carbon; embeddings only", "DNABERT-2; embeddings only",
              "Carbon; embeddings only", "DNABERT-2; embeddings only", ("OT", "bal_OT_common", "MAE", "SR_f1"))
    add_tests(tests, "P. aeruginosa, held-out", pheld, pd_all, "reported", pnew[4][0], "Evo 2 (reported)",
              pnew[4][0], ("OT", "bal_OT_common", "MAE", "SR_f1"))
    for a in (pnew[5][0], pnew[6][0], pnew[7][0]):        # chromosome only, both sides
        add_tests(tests, "P. aeruginosa, held-out", pheld, pd_all, a, pnew[8][0], a, pnew[8][0],
                  ("OT", "bal_OT_common", "MAE", "SR_f1"))
    d8 = [d for d in pd_all if d not in ("aztreonam", "levofloxacin")]
    # what the ARDaP call adds to the embeddings, on the eight antimicrobials it calls
    add_tests(tests, "P. aeruginosa, held-out (8 antimicrobials ARDaP calls)", pheld, d8, "reported", pnew[0][0],
              "Evo 2 (reported)", pnew[0][0], ("OT", "bal_OT_common", "MAE", "vMER", "SR_f1"))
    add_tests(tests, "P. aeruginosa, held-out (exploratory)", pheld, d8, pnew[9][0], "reported", pnew[9][0],
              "Evo 2 (reported)", ("OT", "bal_OT_common", "MAE", "vMER", "SR_f1"))
    if len(pnew) > 10:     # embeddings against catalogued determinants, both with the ARDaP call
        add_tests(tests, "P. aeruginosa, held-out (8 antimicrobials ARDaP calls)", pheld, d8, "reported", pnew[10][0],
                  "Evo 2 (reported)", pnew[10][0], ("OT", "bal_OT_common", "MAE", "cost", "SR_f1"))
    # Same row order as MTB_results: antimicrobial, then model, development before held-out.
    porder = ["Evo 2", "Carbon", "DNABERT-2"] + [f"{m}; embeddings + AMRFinderPlus presence (original specification)"
                                                 for m in ("Evo 2", "Carbon", "DNABERT-2")] + [n for n, _, _ in pnew] + [
              "ARDaP (intermediate counted as resistant)", "ARDaP (intermediate counted as susceptible)"]
    dorder = list(PA_BOUND)
    prow.sort(key=lambda r: (dorder.index(r[0]), porder.index(r[1].replace(" (reported)", "")),
                             r[2] != "development folds"))
    sheet(wb, "PA_results", ["Antimicrobial", "Model", "Set", "Features"] + [NAMES[k] for k in KEYS], prow)

    # ---- P. aeruginosa, study-held-out
    study = {r["Accession number"]: r["Study"] for r in rd(PAMETA)}
    srow, pending = [], []
    for d in PA_BOUND:
        m = reported.get(d)
        lf, cf = LOSO / f"{d}_{m}_loso_predictions.tsv", LOSO / f"{d}_{m}_cv_predictions.tsv"
        if not (lf.exists() and cf.exists()):
            pending.append(d)
            continue
        lo = [r for r in rd(lf) if r["aggregation"] == "median"]
        cv = [r for r in rd(cf) if r["aggregation"] == "median"]
        scoreable = {r["study"] for r in rd(LOSO / f"{d}_studies.tsv") if not r["excluded_because"]}
        for st in sorted(scoreable):
            ids = {r["sample_id"] for r in lo if study.get(r["sample_id"]) == st and not r["fold"].endswith("f-1")}
            if not ids:
                continue
            def sc(rows):
                y = np.array([float(r["log2_mic"]) for r in rows]); p = np.array([float(r["pred"]) for r in rows])
                x = metrics(y, p, PA_BOUND[d]); x["cost"] = 2 * x["vMER"] + x["MER"]; return x
            a = sc([r for r in lo if r["sample_id"] in ids])
            reps = {}
            for r in cv:
                if r["sample_id"] in ids:
                    reps.setdefault(r["fold"].split("_")[1], []).append(r)
            b = [sc(v) for v in reps.values()]
            for name, x in (("Study held out", a), ("Standard cross-validation, same isolates",
                                                    {k: float(np.mean([y[k] for y in b])) for k in ("cost", "vMER", "MER", "OT", "bal_OT", "MAE")})):
                srow.append([d, st, name, len(ids)] + [f(x[k]) for k in ("cost", "vMER", "MER", "OT", "bal_OT", "MAE")])
    sheet(wb, "PA_study_heldout", ["Antimicrobial", "Study", "Evaluation", "Isolates"] +
          [NAMES[k] for k in ("cost", "vMER", "MER", "OT", "bal_OT", "MAE")], srow)

    # Extended Data Figure 2: accuracy of the reported model by true MIC class, held-out
    sys.path.insert(0, os.environ.get("MTB_SCRIPTS") or str(pathlib.Path(__file__).resolve().parent.parent / "mtb_tbprofiler"))
    from make_panel_figures import load as wb_load
    cls_rows = []
    for org, wbp in (("M. tuberculosis", MTBWB), ("P. aeruginosa", PAWB)):
        for r_ in wb_load(wbp)[1]:
            if r_["model"] == "evo2" and r_["split"] == "holdout":
                cls_rows.append([org, r_["drug"], int(float(r_["log2_mic"])), int(float(r_["n"])),
                                 round(100 * float(r_["exact"]), 1), round(100 * float(r_["OT"]), 1)])
    sheet(wb, "MIC_class_accuracy", ["Organism", "Antimicrobial", "True log2 MIC class", "Held-out isolates",
                                     "Exact class (%)", "Essential agreement (%)"], cls_rows)
    co_long, co_summ, co_centre = coresistance_rows()
    sheet(wb, "Co_resistance", ["Organism", "Measure", "Value", "Percent", "n (denominator)", "Note"], co_summ)
    sheet(wb, "Co_resistance_by_centre", ["Organism", "Resistance to", "MIC of", "Collection site (M. tuberculosis) or "
                                           "study (P. aeruginosa)", "AUROC in the centre (%)", "Resistant isolates",
                                           "Susceptible isolates", "AUROC, all isolates (%)"], co_centre)
    sheet(wb, "Co_resistance_AUROC", ["Organism", "Resistance to (row of Fig. 5a,b)", "MIC of (column)", "AUROC (%)",
                                      "Resistant isolates", "Susceptible isolates"], co_long)
    cr, foot = call_role_rows()
    sheet(wb, "Call_roles", ["Organism", "Antimicrobial", "Held-out isolates",
                             "Share of attribution on the tool's call (%)", "Share on the embeddings (%)",
                             "Variation in predicted MIC explained by the call alone (%)",
                             "Within-call Spearman correlation, embeddings' contribution with measured MIC",
                             "Isolates where the framework's S/R call differed from the tool's",
                             "... as % of isolates", "... of which matched the measured phenotype",
                             "... % matched", "Tool F1, development folds", "Tool F1, held-out"], cr + foot)
    sheet(wb, "FDA_criteria", FDA.HEADER, FDA.fda_rows())
    sheet(wb, "Tests", ["Setting", "First model", "Second model", "Measure", "Antimicrobials (pairs)",
                        "Median difference (first minus second)", "First model better in", "P (two-sided Wilcoxon signed-rank)"], tests)

    for line in (
        ["Metrics of the comparisons reported"],
        [],
        ["MTB_results", "M. tuberculosis, development folds and held-out sets: the reported model (Evo 2; target "
                        "loci + TB-Profiler call), the same feature set with Carbon and DNABERT-2, the three models on "
                        "embeddings only, a random forest on every variant in the same loci (variants carried by at "
                        "least 10 development isolates, one binary column each), the model on catalogued mutation "
                        "columns + call, and TB-Profiler alone. Development values are means over the fifteen folds "
                        "of the selected configuration; empty cells were not computed."],
        ["MTB_feature_sets", "Development folds, M. tuberculosis: seven feature sets on the reported target loci "
                             "(embedding sets use Evo 2) and TB-Profiler alone."],
        ["MTB_lineage", "Development set, M. tuberculosis: lineage-held-out folds against random folds of identical sizes."],
        ["PA_results", "P. aeruginosa, development folds and held-out sets, same layout as MTB_results: the three "
                       "embedding models (reported feature set: embeddings of the target loci and acquired genes + "
                       "ARDaP's call, S/I/R as one ordered column; for levofloxacin and aztreonam, which ARDaP does not "
                       "call, + AMRFinderPlus presence instead), the original specification (embeddings + AMRFinderPlus "
                       "presence for every antimicrobial), the three models on embeddings only, every variant in the same "
                       "chromosomal loci plus acquired-gene presence (curated point mutations removed), "
                       "AMRFinderPlus presence alone, the same comparison on chromosomal loci only (no acquired "
                       "genes on either side), AMRFinderPlus presence + ARDaP's call without embeddings, one exploratory "
                       "model (original specification + ARDaP's call; 8 antimicrobials, Evo 2), and ARDaP."],
        ["PA_study_heldout", "P. aeruginosa, development set: each study held out, against standard cross-validation on the same isolates."],
        ["Call_roles", "How the reported model uses the tool's call, held-out isolates (Results, Figure 4). The tool's call "
                       "per isolate is read from the sign of the call feature's SHAP value; this reproduces the tools' "
                       "held-out error rates exactly for M. tuberculosis and within 2.6 points for P. aeruginosa. "
                       "Variation explained: share of the variance of the predicted log2 MIC lying between the isolates "
                       "called resistant and those called susceptible by the tool. Within-call correlation: mean over the "
                       "two call groups (at least 20 isolates) of the Spearman correlation between the summed SHAP of "
                       "the embeddings and the measured log2 MIC. Override: the predicted MIC fell on the other side of "
                       "the breakpoint from the tool's call. Levofloxacin and aztreonam (P. aeruginosa) have no single "
                       "call column and are not listed."],
        ["MIC_class_accuracy", "Extended Data Figure 2: the reported model on the held-out sets, by true log2 MIC class: "
                               "isolates in the class, the share predicted exactly and the share within one doubling "
                               "dilution (classes with fewer than 30 isolates are shown outlined in the figure)."],
        ["Co_resistance", "Figure 5: the numbers quoted in the text. AUROC within and between drug classes, the asymmetry of "
                          "between-class pairs, the first-line to later-drug association in M. tuberculosis, and the counts of "
                          "panels c, d (isolates tested for all ten antimicrobials; the two collections were assembled "
                          "differently and their overall resistance levels are not compared)."],
        ["Co_resistance_by_centre", "Extended Data Figure 3: the AUROC of the between-class pairs with a pooled AUROC of at "
                                    "least 80%, within each collection site (M. tuberculosis) or study (P. aeruginosa) "
                                    "with at least 20 resistant and 20 susceptible isolates, one row per dot."],
        ["Co_resistance_AUROC", "Figure 5a,b: AUROC of every pair with at least ten resistant isolates, with the numbers of "
                                "resistant and susceptible isolates."],
        ["FDA_criteria", "Figure 3c and the Results: held-out error counts against the FDA acceptance criteria for AST systems "
                         "(guidance of 28 August 2009, section X.A and Table 8). Major errors: at most 3% of the susceptible "
                         "isolates. Very major errors: the exact binomial 95% interval of the true rate has an upper limit of at "
                         "most 7.5% and a lower limit of at most 1.5%; the 1.5% is a limit of the interval, not a ceiling for the "
                         "rate. The guidance concerns commercial short-incubation MIC systems and excludes anti-mycobacterial "
                         "agents and genotype-based devices; it is applied here by analogy. Error counts are the rates of the "
                         "held-out rows times the held-out isolates of the class (ARDaP's amikacin rate implies 36 resistant "
                         "isolates scored). Bedaquiline and delamanid (12 and 15 resistant) are listed but not scored. "
                         "The distance to the targets of Figure 3c uses the stricter point targets vMER 1.5% and MER 3%."],
        ["Tests", "Paired two-sided Wilcoxon signed-rank tests across antimicrobials. M. tuberculosis: all ten for "
                  "the MIC measures (essential agreement, balanced essential agreement, exact, mean absolute error), which rest on every "
                  "held-out isolate; the eight with a verdict for cost, vMER, MER and F1, since bedaquiline and "
                  "delamanid have only 12 and 15 resistant held-out isolates. P. aeruginosa: all ten (eight for the "
                  "comparisons with ARDaP). The column 'Antimicrobials (pairs)' gives the number used. No adjustment "
                  "for multiple comparisons."],
        [],
        ["Measures", "vMER: resistant isolates called susceptible, as a fraction of resistant isolates. MER: susceptible isolates "
                     "called resistant, as a fraction of susceptible isolates. Essential agreement (also called one-tier accuracy): predictions within one doubling "
                     "dilution. Balanced: averaged over MIC classes."],
    ):
        readme.append(line)
    readme["A1"].font = Font(bold=True, size=13)
    readme.column_dimensions["A"].width = 22
    readme.column_dimensions["B"].width = 120
    wb.save(out)
    print(f"[OK] {out}")
    for n in wb.sheetnames[1:]:
        print(f"  {n}: {wb[n].max_row - 1} rows")
    single = {"cefepime", "levofloxacin"}
    print("  PA study-held-out predictions missing for the reported model:",
          [d for d in pending if d not in single] or "none")
    return 0


if __name__ == "__main__":
    sys.exit(main())
