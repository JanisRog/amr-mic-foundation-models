#!/usr/bin/env python3
"""Supplementary table: curated chromosomal target loci for P. aeruginosa.

Read from pseaer/build_pseaer_target_sets.py, which is what the pipeline ran, so
the table cannot drift from the target sets. The mechanism label is the grouping
the loci were declared under there. The reference column is left empty: the script
records no source for the list, and it must be filled by the author.
"""

from __future__ import annotations

import os
import pathlib
import sys

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "pseaer"))
import build_pseaer_target_sets as T  # noqa: E402

MECHANISM = [
    ("Fluoroquinolone targets (DNA gyrase, topoisomerase IV)", ["gyrA", "gyrB", "parC", "parE"]),
    ("Carbapenem porin", ["oprD"]),
    ("AmpC cephalosporinase and its regulation", T.AMPC),
    ("MexAB-OprM efflux and its regulators", T.MEXAB),
    ("MexXY efflux and its repressor", T.MEXXY),
    ("MexCD-OprJ efflux and its regulator", ["mexC", "mexD", "oprJ", "nfxB"]),
    ("MexEF-OprN efflux and its regulators", ["mexE", "mexF", "oprN", "mexT", "mexS"]),
    ("Elongation factor G", ["fusA1"]),
    ("Two-component system ParRS", ["parR", "parS"]),
]
ORDER = ["ciprofloxacin", "levofloxacin", "tobramycin", "amikacin", "imipenem", "meropenem",
         "ceftazidime", "cefepime", "aztreonam", "piperacillin-tazobactam"]


def main() -> int:
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1
                       else pathlib.Path("AMRPredict_pa_target_loci.xlsx"))
    listed = [g for _, genes in MECHANISM for g in genes]
    assert sorted(listed) == sorted(T.RESOLVE), "mechanism groups do not cover the curated loci exactly"
    wb = Workbook()
    ws = wb.active
    ws.title = "PA_target_loci"
    ws.append(["Locus", "Mechanism"] + ORDER + ["Antimicrobials (n)", "Reference"])
    for mech, genes in MECHANISM:
        for g in genes:
            used = [d for d in ORDER if g in T.DRUGS[d]]
            ws.append([g, mech] + [("x" if d in used else None) for d in ORDER] + [len(used), None])
    ws.append([])
    ws.append(["Loci per antimicrobial", None] + [len(T.DRUGS[d]) for d in ORDER])
    ws.append(["Targets per antimicrobial (gene body and 300 bp upstream)", None]
              + [2 * len(T.DRUGS[d]) for d in ORDER])
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 12
    ws.column_dimensions["B"].width = 52
    ws.freeze_panes = "C2"
    wb.save(out)
    print(f"[OK] {out}: {len(listed)} loci; targets per antimicrobial "
          + ", ".join(f"{d} {2 * len(T.DRUGS[d])}" for d in ORDER))
    return 0


if __name__ == "__main__":
    sys.exit(main())
