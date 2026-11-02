#!/usr/bin/env python3
"""Supplementary table: target loci and resistance determinants, both organisms.

M. tuberculosis targets are read from the per-target summaries of the reported
model (one row per target the model was fitted on), so the table is the model's
own target list. A target is labelled as TB-Profiler-annotated if it is also in
the catalogue-derived model's target list, and otherwise as coming from the
CRyPTIC association study. P. aeruginosa chromosomal loci are read from
pseaer/build_pseaer_target_sets.py; AMRFinderPlus determinants from the per-
antimicrobial presence matrices the models were fitted on.
"""

from __future__ import annotations

import csv
import os
import pathlib
import sys

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "pseaer"))
import build_pa_target_loci_table as PA  # noqa: E402
import build_pseaer_target_sets as T  # noqa: E402

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
WIDER = HOME / "mtb_shap_gwfrule1k/shap_gwfrule1k"
EARLIER = HOME / "Downloads"
AMR = HOME / "pa_strict/labels_imputed_strict/amrfinder/by_drug_tsv"
MTB = ["isoniazid", "rifampicin", "ethambutol", "ethionamide", "kanamycin", "amikacin",
       "levofloxacin", "moxifloxacin", "bedaquiline", "delamanid"]
REGION = {"gene": "gene body", "upstream": "300 bp upstream"}
# Gene symbols and products of H37Rv, from the NCBI feature table of NC_000962.3
# (efetch, db=nuccore, rettype=ft; retrieved 2026-10-13).
# Two loci the feature table describes without a product line.
SPECIAL = {"EBG00000313325": "16S ribosomal RNA", "Rv2427A": "pseudogene (oxyR')"}
ANNOT = {r["locus_tag"]: r for r in csv.DictReader(open(HERE / "h37rv_NC_000962.3_genes.tsv"), delimiter="\t")}


def targets(path):
    return [r for r in csv.DictReader(open(path), delimiter="\t") if not r["target_id"].startswith("rule")]


def sheet(wb, title, header, rows, widths):
    ws = wb.create_sheet(title)
    ws.append(header)
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="top")
    for r in rows:
        ws.append(r)
    for i, w in enumerate(widths):
        ws.column_dimensions[chr(65 + i)].width = w
    ws.freeze_panes = "A2"
    return ws


def main() -> int:
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path("AMRPredict_target_loci.xlsx"))
    wb = Workbook()
    readme = wb.active
    readme.title = "README"

    rows, counts, mismatch = [], [], []
    for d in MTB:
        wide = targets(WIDER / f"{d}_evo2_target_summary.tsv")
        early = {r["target_id"] for r in targets(EARLIER / f"{d}_evo2_target_summary.tsv")}
        ids = {r["target_id"] for r in wide}
        if not early <= ids:
            sys.exit(f"[ERROR] {d}: {len(early - ids)} catalogue-derived targets are absent from the reported set")
        for r in sorted(wide, key=lambda r: (r["locus"], r["region"])):
            a = ANNOT.get(r["locus"], {})
            if r["gene"] and a.get("gene") and r["gene"] != a["gene"]:
                mismatch.append((d, r["locus"], r["gene"], a["gene"]))
            product = a.get("product") or SPECIAL.get(r["locus"]) or "not annotated in NC_000962.3"
            rows.append([d, r["locus"], r["gene"] or a.get("gene") or None, product,
                         REGION.get(r["region"], r["region"]),
                         "TB-Profiler database" if r["target_id"] in early else "CRyPTIC association study"])
        counts.append([d, len(wide), len({r["locus"] for r in wide}), len(early), len(wide) - len(early)])
    sheet(wb, "MTB_targets", ["Antimicrobial", "Locus", "Gene", "Product", "Region", "Source"], rows,
          (16, 12, 10, 58, 18, 28))
    named = sum(1 for r in rows if r[2])
    print(f"  MTB targets with a gene symbol: {named} of {len(rows)}; with a product: "
          f"{sum(1 for r in rows if r[3] != 'not annotated in NC_000962.3')}; symbol differs from NCBI for {len(mismatch)}: {mismatch[:8]}")
    sheet(wb, "MTB_counts", ["Antimicrobial", "Targets", "Loci", "Targets from the TB-Profiler database",
                             "Targets from the CRyPTIC association study"], counts, (16, 10, 8, 22, 24))

    prow = []
    for mech, genes in PA.MECHANISM:
        for g in genes:
            prow.append([g, mech] + [("x" if g in T.DRUGS[d] else None) for d in PA.ORDER])
    ws = sheet(wb, "PA_chromosomal_loci", ["Locus", "Mechanism"] + PA.ORDER, prow, (10, 52) + (14,) * 10)
    ws.append([])
    ws.append(["Targets (gene body and 300 bp upstream)", None] + [2 * len(T.DRUGS[d]) for d in PA.ORDER])

    arow, acount = [], []
    for d in PA.ORDER:
        cols = open(AMR / f"{d}.tsv").readline().rstrip("\n").split("\t")[1:]
        for c in cols:
            name = c[4:] if c.startswith("amr_") else c
            arow.append([d, name])
        acount.append([d, len(cols)])
    sheet(wb, "PA_AMRFinderPlus", ["Antimicrobial", "Determinant (presence column)"], arow, (24, 34))
    sheet(wb, "PA_AMRFinderPlus_counts", ["Antimicrobial", "Presence columns"], acount, (24, 18))

    for line in (
        ["Target loci and resistance determinants"],
        [],
        ["MTB_targets", "Every target of the reported M. tuberculosis model, per antimicrobial. A target is the "
                        "gene body of a locus or the 300 bp upstream of it. Gene symbols and products are those "
                        "of the H37Rv reference annotation (NC_000962.3); loci without a symbol have none in "
                        "that annotation. Source: the TB-Profiler database "
                        "(tbdb, commit 72ef6fa) or the CRyPTIC genome-wide association study."],
        ["MTB_counts", "Targets and loci per antimicrobial, by source."],
        ["PA_chromosomal_loci", "The 31 curated chromosomal loci of P. aeruginosa, their mechanism, and the "
                                "antimicrobials for which each is a target (x)."],
        ["PA_AMRFinderPlus", "Acquired genes and resistance mutations detected by AMRFinderPlus and used as "
                             "presence columns, per antimicrobial (carried by at least ten development isolates)."],
        ["PA_AMRFinderPlus_counts", "Number of presence columns per antimicrobial."],
    ):
        readme.append(line)
    readme["A1"].font = Font(bold=True, size=13)
    readme.column_dimensions["A"].width = 28
    readme.column_dimensions["B"].width = 110
    wb.save(out)
    print(f"[OK] {out}")
    print("  MTB targets per antimicrobial:", ", ".join(f"{c[0]} {c[1]} ({c[2]} loci; {c[3]}+{c[4]})" for c in counts))
    print("  PA presence columns:", ", ".join(f"{c[0]} {c[1]}" for c in acount))
    return 0


if __name__ == "__main__":
    sys.exit(main())
