#!/usr/bin/env python3
"""Supplementary table: design alternatives evaluated and not adopted.

The rows are a curated record, not a computation: each states what was adopted,
what was tried instead, the data it was tried on, and what happened. Numbers are
copied from EVALUATION_PROTOCOL.md (sections given per row) except the tree
aggregation row, recomputed from the development folds of the reported model.

"TO COMPLETE" in the last column marks a row that must not go to print as it is.
"""

from __future__ import annotations

import os
import pathlib
import sys

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill

HEADER = ["Organism", "Step", "Adopted", "Alternative evaluated", "Evaluated on",
          "Outcome", "Reason not adopted", "Source / status"]
ROWS = [
    ["MTB", "Embedding layer", "Final normalisation layer",
     "Normalised intermediate layers of Evo 2 (blocks 8, 16, 24 and 28)",
     "(a) Correlation of the projected embedding with log2 MIC at individual resistance loci. "
     "(b) Isoniazid, rifampicin-susceptible isolates of an earlier train/test partition (3,651 "
     "training, 898 test), one random forest configuration",
     "(a) Lower at every locus, rising with depth: katG |r| 0.295, 0.399, 0.419 and 0.739, against "
     "0.809 for the final layer. (b) Cost 0.240, 0.256, 0.253 and 0.213, against 0.174 for the final "
     "layer (vMER 0.112, 0.118, 0.118 and 0.099, against 0.079); balanced essential agreement 0.769, "
     "0.749, 0.770 and 0.788, against 0.781; mean absolute error 0.600, 0.628, 0.618 and 0.586, "
     "against 0.585",
     "Higher cost: every intermediate layer missed more resistant isolates",
     "Protocol §7; benchmark 03_mtb_extraction_benchmark (blocks_*_post_norm, norm), "
     "summary.tsv read 2026-10-13"],
    ["MTB", "Collapsing position vectors to one vector per sequence", "Mean over positions",
     "Differences from the most frequent allele of the target, collapsed by the maximum absolute "
     "difference or by the k most different positions",
     "Isoniazid, rifampicin-susceptible isolates of an earlier train/test partition (3,651 training, "
     "898 test), one random forest configuration",
     "Maximum: cost 0.178, balanced essential agreement 0.780, mean absolute error 0.605. Top-k: cost "
     "0.252 (vMER 0.118), balanced essential agreement 0.795, mean absolute error 0.605. Mean: cost "
     "0.174 (vMER 0.079), balanced essential agreement 0.781, mean absolute error 0.585",
     "No gain (maximum), or higher cost through more resistant isolates missed (top-k)",
     "Benchmark 03_mtb_extraction_benchmark (delta_maxabs, delta_topk, norm), summary.tsv read "
     "2026-10-13"],
    ["MTB", "Target selection", "Fixed target loci per antimicrobial",
     "Minimum-redundancy maximum-relevance selection within each fold over 8,262 genome-wide targets",
     "Isoniazid development folds",
     "The rifampicin-resistance gene rpoB had the highest relevance for isoniazid (|r| = 0.616), ahead "
     "of every isoniazid locus; adjusting for lineage left it first (0.594). Selecting 256 targets in "
     "place of 64 changed balanced essential agreement by at most 0.007",
     "Selects markers of co-resistance to other antimicrobials", "Protocol §1"],
    ["MTB", "Target set",
     "Loci that TB-Profiler annotates for the antimicrobial together with the loci reported for it in the "
     "CRyPTIC genome-wide association study, minus loci annotated only for other antimicrobials (4 to 72 "
     "targets)",
     "Loci that TB-Profiler annotates for the antimicrobial only (narrower set)",
     "Development folds, then held-out sets (eight antimicrobials with a verdict; Evo 2)",
     "Development folds: balanced essential agreement +0.027 for the wider set in 7 of 8 antimicrobials, cost "
     "tied. Held-out: cost -0.0004 (4 of 8, P = 0.60), balanced essential agreement +0.015 (6 of 8, P = 0.063), "
     "mean absolute error -0.008 (6 of 8, P = 0.091). With whole lineages held out, the wider set had no "
     "advantage",
     "Gain on random folds only, not significant on the held-out sets. The wider set was adopted on "
     "development evidence after the narrower set had been scored on the held-out set",
     "Protocol §19, §27, §28 (held-out results: holdout_arm_comparison.tsv); added 2026-10-20"],
    ["MTB", "Projection", "One projection shared across targets, four components",
     "Separate principal components per target; components masked by model importance",
     "Isoniazid development folds, fixed targets, three embedding models",
     "All differences in balanced essential agreement within 0.005; cost differed significantly in 2 of "
     "12 comparisons",
     "No gain, and more features (74 to 121 against 49 for isoniazid)", "Protocol §9-§10"],
    ["MTB", "Learner", "Random forest regression",
     "Histogram gradient boosting; extremely randomised trees; ridge regression",
     "Isoniazid development folds",
     "Gradient boosting within 0.002 of random forests on average but more variable across "
     "configurations; extremely randomised trees five to ten times more variable; ridge regression "
     "about 0.3 lower in balanced essential agreement",
     "Random forests were the most consistent across configurations, which matters when a single "
     "configuration is selected", "Protocol §9"],
    ["MTB", "Aggregation of tree predictions", "Median", "Mean",
     "Development folds of the reported model, eight antimicrobials, same forests",
     "Evo 2: cost equal (median difference +0.001), balanced essential agreement higher with the median "
     "in 7 of 8 (+0.002), mean absolute error lower in 6 of 8 (-0.006). Carbon and DNABERT-2: same "
     "pattern (balanced essential agreement higher in 7 and 6 of 8)",
     "Slightly worse placement on the MIC scale at equal cost",
     "Recomputed 2026-10-13 from slim_arms/*/runs/*_gwasfiltered_gwfrule1k/folds.tsv"],
    ["P. aeruginosa", "Target loci", "Curated chromosomal loci and AMRFinderPlus determinants",
     "Loci from a genome-wide association analysis (pyseer v1.3.11, linear mixed models), used alone "
     "or added to AMRFinderPlus-derived targets",
     "Ciprofloxacin; an earlier partition of the collection",
     "Higher cost when used alone (median over configurations 0.38 to 0.50, against 0.33 to 0.36 for "
     "AMRFinderPlus-derived targets) and no clear gain when added (0.30 to 0.34)",
     "No clear gain",
     "evo2_local/02_pseaer/antibiotics/ciprofloxacin/results/holdout_results/sweep_summary_*.tsv, "
     "read 2026-10-13; wording agreed with the author"],
    ["P. aeruginosa", "Catalogue input of the reported model",
     "ARDaP's call as one feature for the eight antimicrobials ARDaP calls (AMRFinderPlus presence columns "
     "for levofloxacin and aztreonam)",
     "AMRFinderPlus presence columns for all ten antimicrobials (original specification, declared in "
     "advance)",
     "Development folds and held-out sets, the eight antimicrobials ARDaP calls (Evo 2)",
     "Development folds: lower cost with the adopted input in 8 of 8 antimicrobials (median -0.107, "
     "P = 0.0078) and higher F1 in 8 of 8 (+0.027, P = 0.0078). Held-out: lower cost in 7 of 8 (-0.088, "
     "P = 0.023), fewer very major errors in 8 of 8 (-0.029, P = 0.0078), F1 +0.027 (6 of 8, P = 0.055), "
     "essential agreement +0.019 (6 of 8, P = 0.11), mean absolute error -0.007 (4 of 8, P = 0.46). "
     "Ciprofloxacin essential agreement was higher under the original specification (85.5% against 82.7%)",
     "Higher cost and lower F1 on the development folds in every antimicrobial, and consistency with the "
     "M. tuberculosis model, which takes the tool's call as a feature. The change followed held-out "
     "scoring of the original specification",
     "Supplementary Table 3 (PA_results, original specification rows), differences recomputed "
     "2026-10-20; P. aeruginosa protocol, decision of 2026-10-19; added 2026-10-20"],
    ["P. aeruginosa", "Censored MICs", "Strict bounds moved one dilution outward, inclusive bounds kept",
     "All bounds moved one dilution outward", "All MICs",
     "Imputed values 12.4% against 4.9%; resistant and susceptible assignments identical",
     "An inclusive bound names a value that is consistent with the measurement",
     "P. aeruginosa protocol"],
]


def main() -> int:
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1
                       else pathlib.Path("AMRPredict_design_alternatives.xlsx"))
    wb = Workbook()
    ws = wb.active
    ws.title = "design_alternatives"
    # The last column (source / status) documents where each number comes from. It
    # stays in this file and is left out of the workbook, which is the printed table.
    ws.append(HEADER[:-1])
    for c in ws[1]:
        c.font = Font(bold=True)
    for r in ROWS:
        ws.append(r[:-1])
    for col, w in zip("ABCDEFG", (14, 26, 34, 44, 40, 62, 48)):
        ws.column_dimensions[col].width = w
    flag = PatternFill("solid", fgColor="FFE699")
    for row, src in zip(ws.iter_rows(min_row=1), [HEADER] + ROWS):
        for c in row:
            c.alignment = Alignment(wrap_text=True, vertical="top")
        if "TO COMPLETE" in src[-1]:
            for c in row:
                c.fill = flag
    ws.freeze_panes = "C2"
    wb.save(out)
    todo = sum("TO COMPLETE" in r[-1] for r in ROWS)
    print(f"[OK] {out}: {len(ROWS)} rows, {todo} marked TO COMPLETE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
