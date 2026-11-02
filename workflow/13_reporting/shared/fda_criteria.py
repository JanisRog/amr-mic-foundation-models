#!/usr/bin/env python3
"""The held-out error rates against the FDA acceptance criteria for antimicrobial susceptibility test systems.

Source: FDA, Class II Special Controls Guidance Document: Antimicrobial Susceptibility Test (AST) Systems
(28 August 2009), section X.A (p. 18) and Table 8 (p. 41):
  major error rate (MER)       <= 3% of the susceptible isolates;
  very major error rate (vMER) exact binomial 95% interval of the true rate with upper limit <= 7.5% and
                               lower limit <= 1.5% (the 1.5% is a limit of the interval, not a ceiling for the estimate).
The guidance is written for commercial short-incubation MIC systems and excludes anti-mycobacterial agents and
genotype-based devices; here it is applied by analogy.

Counts: held-out resistant and susceptible isolates per antimicrobial from the isolate workbook (Supplementary
Table 1); errors = rate x isolates of the class, from the held-out rows of the results workbooks.

Run alone to print the summary; build_comparison_metrics_table.py imports fda_rows() for the FDA_criteria sheet.
"""

from __future__ import annotations

import os
import pathlib
import sys

from openpyxl import load_workbook
from scipy.stats import beta

HOME = pathlib.Path(os.environ.get("AMR_INPUT_DIR", "inputs"))           # staged inputs (README): set AMR_INPUT_DIR
ISOLATES = HOME / "AMRPredict_isolates_used.xlsx"
WORKBOOKS = {"mtb": HOME / "AMRPredict_panel_results_wider.xlsx",
             "pseaer": HOME / "AMRPredict_pseaer_reported_results.xlsx"}
ORG = {"mtb": "M. tuberculosis", "pseaer": "P. aeruginosa"}
BASELINE = {"mtb": ("TB-Profiler 6.3.0", "TB-Profiler"), "pseaer": ("ARDaP (I=R)", "ARDaP")}
NO_VERDICT = {"bedaquiline", "delamanid"}          # 12 and 15 resistant held-out isolates; not scored in the comparisons
MER_LIMIT, VMER_UPPER, VMER_LOWER = 0.03, 0.075, 0.015
FLOOR_N = 48                                        # fewest resistant isolates for which zero errors can pass (Table 8)


def exact_ci(k, n, level=0.95):
    """Clopper-Pearson interval, the 'exact confidence interval based on the binomial distribution' of Table 8."""
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else beta.ppf(a, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a, k + 1, n - k)
    return float(lo), float(hi)


def read(path, name):
    ws = load_workbook(path, read_only=True)[name]
    rows = list(ws.iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def fda_rows():
    """One row per organism, antimicrobial and system (reported model, tool)."""
    out = []
    for k in ("mtb", "pseaer"):
        counts = {r["antimicrobial"]: r for r in read(ISOLATES, "MTB_counts" if k == "mtb" else "PA_counts")}
        summ = read(WORKBOOKS[k], "summary")
        base, base_label = BASELINE[k]
        for r in read(WORKBOOKS[k], "summary"):
            if r["split"] != "holdout" or r["model"] not in ("evo2", base) or r["vMER"] in (None, ""):
                continue
            d = r["drug"]
            c = counts.get(d)
            if c is None:
                continue
            nr, ns = int(c["held-out resistant"]), int(c["held-out susceptible"])
            vm, me = round(float(r["vMER"]) * nr), round(float(r["MER"]) * ns)
            lo, hi = exact_ci(vm, nr)
            v_ok = hi <= VMER_UPPER and lo <= VMER_LOWER
            m_ok = me / ns <= MER_LIMIT
            out.append([ORG[k], d, "Reported model" if r["model"] == "evo2" else base_label, nr, vm,
                        round(100 * vm / nr, 2), round(100 * lo, 2), round(100 * hi, 2),
                        "yes" if v_ok else "no", ns, me, round(100 * me / ns, 2), "yes" if m_ok else "no",
                        "yes" if v_ok and m_ok else "no", "no verdict (too few resistant)" if d in NO_VERDICT else ""])
    return sorted(out, key=lambda x: (x[0] != "M. tuberculosis", x[1], x[2] != "Reported model"))


HEADER = ["Organism", "Antimicrobial", "System", "Resistant held-out isolates", "Very major errors",
          "vMER (%)", "Exact 95% interval, lower (%)", "Exact 95% interval, upper (%)",
          "vMER criterion met (upper <= 7.5%, lower <= 1.5%)", "Susceptible held-out isolates", "Major errors",
          "MER (%)", "MER criterion met (<= 3%)", "Both criteria met", "Note"]


def main() -> int:
    rows = fda_rows()
    scored = [r for r in rows if not r[-1]]
    for org in ORG.values():
        for system in ("Reported model", "tool"):
            sel = [r for r in scored if r[0] == org and (r[2] == "Reported model") == (system == "Reported model")]
            if not sel:
                continue
            n = len(sel)
            print(f"[INFO] {org} {sel[0][2]}: both {sum(r[13] == 'yes' for r in sel)}/{n}; "
                  f"vMER {sum(r[8] == 'yes' for r in sel)}/{n}; MER {sum(r[12] == 'yes' for r in sel)}/{n}")
    for r in rows:
        print("  ", r[0][:2], r[1][:14].ljust(14), r[2][:14].ljust(14), f"R={r[3]:>4} vm={r[4]:>3} ({r[5]:5.2f}%, "
              f"{r[6]:5.2f}-{r[7]:5.2f}) {r[8]:>3}  S={r[9]:>4} maj={r[10]:>3} ({r[11]:5.2f}%) {r[12]:>3}  {r[14]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
