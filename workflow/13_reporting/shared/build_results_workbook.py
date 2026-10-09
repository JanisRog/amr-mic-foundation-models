#!/usr/bin/env python3
"""One workbook for the panel: development CV, the holdout, and the incumbent.

The numbers that matter for a drug end up spread across a development run
directory, a holdout run directory and a separately computed catalogue baseline,
in three different shapes. Assembling them by hand is how a transcription error
gets into a paper, and it has to be done once per drug -- so `--append` grows one
file across the panel, with all three models under each antibiotic, rather than
leaving twelve files to be merged at the end.

Three sheets:

  **summary**   one row per model per split, plus the baseline rows. Development
                is the *chosen* config only -- not the best of the search, which
                would be the max-over-configs bias section 3 names as the largest
                single bias in this project. Holdout is the same config, scored
                once.
  **per_class** ordinal accuracy by MIC class, for the chosen config on both
                splits. Isoniazid's headline bal_OT of 0.778 averages classes
                running from 0.99 to 0.15, and that is invisible in any summary
                row, so it gets a sheet rather than a footnote.
  **configs**   the selected hyperparameters and how they were selected, so the
                workbook says what was chosen and why without needing the
                protocol open beside it.

The catalogue baseline is read from a TSV rather than recomputed here, because it
needs the MIC table and the fold file, which live on the cluster. Generating it
there and reading it here keeps the numbers off the clipboard.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import pathlib
import sys

METRICS = ["bal_OT", "bal_OT_common", "OT", "exact", "MAE", "RMSE", "spearman",
           "vMER", "MER", "SR_f1", "SR_acc", "SR_AUROC_raw"]


def cost(r: dict) -> float:
    return 2 * float(r["vMER"]) + float(r["MER"])


def read_run(d: pathlib.Path, config: str, aggregation: str = "median") -> dict:
    """The chosen config's row, averaged over whatever partitions the run has."""
    rows = [r for r in csv.DictReader(open(d / "folds.tsv"), delimiter="\t")
            if r["aggregation"] == aggregation and r["config"] == config]
    if not rows:
        raise SystemExit(f"[ERROR] {config} / {aggregation} not in {d}/folds.tsv")
    out = {"n_folds": len(rows)}
    for m in METRICS + ["n_features"]:
        vals = [float(r[m]) for r in rows if r.get(m) not in (None, "")]
        out[m] = sum(vals) / len(vals) if vals else ""
    out["cost"] = sum(cost(r) for r in rows) / len(rows)
    return out


def read_per_class(d: pathlib.Path, config: str,
                   aggregation: str = "median") -> list[dict]:
    out = []
    for r in csv.DictReader(open(d / "per_class.tsv"), delimiter="\t"):
        if r["aggregation"] == aggregation and r["config"] == config:
            out.append({"log2_mic": float(r["log2_mic"]), "n": float(r["n"]),
                        "exact": float(r["exact"]), "OT": float(r["OT"])})
    return sorted(out, key=lambda r: r["log2_mic"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", required=True)
    ap.add_argument("--model", action="append", default=[],
                    metavar="NAME=DEV_DIR:HOLDOUT_DIR:CONFIG_TAG",
                    help="repeatable. Optional: a drug can be added with only a "
                         "--baseline, so the comparator is recorded before any "
                         "model exists and the model rows are appended later. "
                         "Re-running the drug replaces its rows, so nothing is "
                         "duplicated by doing it in two passes")
    ap.add_argument("--baseline", default="",
                    help="TSV of catalogue rows: label, split, and any of the "
                         "metric columns")
    ap.add_argument("--selection-note", default="",
                    help="how the config was chosen, recorded on the configs sheet")
    ap.add_argument("--append", action="store_true",
                    help="add this drug to an existing workbook instead of "
                         "starting one. The panel is finished a drug at a time "
                         "over weeks, and one file per drug is how a comparison "
                         "table stops being assembled and starts being retyped. "
                         "Re-running a drug already present replaces its rows "
                         "rather than duplicating them")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    try:
        from openpyxl import Workbook, load_workbook
        from openpyxl.styles import Font
    except ImportError:
        raise SystemExit("[ERROR] needs openpyxl (pip install openpyxl)")

    head = ["drug", "model", "split", "n_folds", "n_features"] + METRICS + ["cost"]
    out = pathlib.Path(args.out)
    if args.append and out.exists():
        wb = load_workbook(out)
        ws, ws2, ws3 = wb["summary"], wb["per_class"], wb["configs"]
        # a re-run of the same drug should correct it, not append a second copy
        for sheet in (ws, ws2, ws3):
            drop = [c.row for c in sheet["A"]
                    if c.row > 1 and c.value == args.drug]
            for r in reversed(drop):
                sheet.delete_rows(r)
        if drop_note := [c.row for c in ws3["A"] if c.value == "generated"]:
            for r in reversed(drop_note):
                ws3.delete_rows(r)
        print(f"[INFO] appending {args.drug} to {out}")
    else:
        wb = Workbook()
        ws = wb.active
        ws.title = "summary"
        ws.append(head)
        ws2 = wb.create_sheet("per_class")
        ws2.append(["drug", "model", "split", "log2_mic", "n", "exact", "OT"])
        ws3 = wb.create_sheet("configs")
        ws3.append(["drug", "model", "config", "selection_rule"])
        for sheet in (ws, ws2, ws3):
            for c in sheet[1]:
                c.font = Font(bold=True)

    if not args.model and not args.baseline:
        raise SystemExit("[ERROR] nothing to write: pass --model and/or --baseline")
    per_class_rows, config_rows = [], []
    for spec in args.model:
        name, _, rest = spec.partition("=")
        dev_s, hold_s, tag = rest.split(":")
        for split, path in (("development_cv", dev_s), ("holdout", hold_s)):
            if not path:
                continue
            d = pathlib.Path(path)
            r = read_run(d, tag)
            ws.append([args.drug, name, split, r["n_folds"], r.get("n_features", "")]
                      + [r.get(m, "") for m in METRICS] + [r["cost"]])
            for pc in read_per_class(d, tag):
                per_class_rows.append([args.drug, name, split, pc["log2_mic"],
                                       pc["n"], pc["exact"], pc["OT"]])
        config_rows.append([args.drug, name, tag, args.selection_note])

    if args.baseline:
        for r in csv.DictReader(open(args.baseline), delimiter="\t"):
            ws.append([args.drug, r.get("label", "catalogue"),
                       r.get("split", ""), "", ""]
                      + [r.get(m, "") for m in METRICS]
                      + [r.get("cost", "")])

    for row in per_class_rows:
        ws2.append(row)
    for row in config_rows:
        ws3.append(row)
    ws3.append(["generated",
                datetime.datetime.now().isoformat(timespec="seconds")])

    for sheet in (ws, ws2, ws3):
        for col in sheet.columns:
            width = max(len(str(c.value)) for c in col if c.value is not None)
            sheet.column_dimensions[col[0].column_letter].width = min(width + 2, 30)
        sheet.freeze_panes = "A2"

    out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out)
    print(f"[OK] {out}")
    print(f"     summary {ws.max_row - 1} rows, per_class {ws2.max_row - 1}, "
          f"configs {ws3.max_row - 3}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
