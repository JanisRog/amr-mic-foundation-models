#!/usr/bin/env python3
"""Render the panel figures and their numbers as one reviewable page.

The table is built from the workbook in the same pass that names the figures, so
the numbers under a figure cannot drift from the numbers in it. Every drug x
model x split row is present: three of the series colors sit below 3:1 contrast
on a white ground, and the rule for that is a visible table, not a promise that
the marks are distinguishable.

The figures are light-mode PNGs. They are therefore mounted on an explicit white
plate that stays white in dark mode -- inverting a figure would invert its
palette, and the palette is the part that was validated.
"""

from __future__ import annotations

import argparse
import html
import pathlib
import sys

# Ordered by how much of the ladder the drug actually has: the wide-range drugs
# first, the two newest last, which is also the order the caveats accumulate in.
SCATTER_ORDER = ["isoniazid", "rifampicin", "ethambutol", "ethionamide",
                 "kanamycin", "amikacin", "levofloxacin", "moxifloxacin",
                 "bedaquiline", "delamanid"]

FIELDS = [("cost", "cost"), ("vMER", "vMER"), ("MER", "MER"), ("SR_f1", "F1"),
          ("OT", "OT"), ("exact", "exact"), ("SR_AUROC_raw", "AUROC")]

FIGURES = [
    ("fig1_panel.png", "Figure 1",
     "The panel, held-out set",
     "Cost and F1 for all three models against TB-Profiler, eight drugs with "
     "enough resistant isolates to carry a verdict. Lower cost on 4, higher F1 "
     "on 7. The two metrics disagree because the model shifts toward "
     "specificity, which F1 rewards and a 2:1 very-major-error weighting "
     "punishes; neither is reported alone."),
    ("fig2_tradeoff.png", "Figure 2",
     "Where each drug moved in the error plane",
     "The same result as figure 1, drawn as a direction rather than a margin. "
     "Ethionamide and kanamycin cross iso-cost lines away from the origin while "
     "gaining F1 — that single picture is the whole disagreement between the "
     "two metrics."),
    ("fig3_ladder.png", "Figure 3",
     "The ordinal ladder, and where it fails",
     "The comparator emits one bit per isolate and has no entry in either "
     "panel. The heatmap is the caveat: the headline is carried by the dense "
     "MIC classes, and where a distribution is bimodal the model cannot place "
     "the gap between the modes."),
    ("fig4_newdrugs.png", "Figure 4",
     "Bedaquiline and delamanid",
     "Twelve and fifteen resistant isolates in the holdout. Shown because "
     "omitting them would be selective, separated because ranking them beside "
     "the rest would imply a precision they do not have."),
]

SUPPLEMENTARY = [
    ("figS1_stability.png", "Figure S1",
     "Development against holdout",
     "One point per drug × model, 30 per panel. The only departures from the "
     "equality line are the two underpowered drugs."),
]

CSS = """
:root{
  --ground:#fbfbf9; --surface:#ffffff; --ink:#111110; --ink-2:#4a4945;
  --muted:#84837c; --rule:#e4e3dc; --accent:#1f62b0; --plate:#ffffff;
  --shadow:0 1px 2px rgba(17,17,16,.05),0 8px 24px -12px rgba(17,17,16,.18);
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    --ground:#121211; --surface:#1a1a19; --ink:#f4f3ee; --ink-2:#c3c2b7;
    --muted:#8d8c84; --rule:#2c2c2a; --accent:#6da7ec; --plate:#ffffff;
    --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 28px -14px rgba(0,0,0,.7);
  }
}
:root[data-theme="dark"]{
  --ground:#121211; --surface:#1a1a19; --ink:#f4f3ee; --ink-2:#c3c2b7;
  --muted:#8d8c84; --rule:#2c2c2a; --accent:#6da7ec; --plate:#ffffff;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 28px -14px rgba(0,0,0,.7);
}
*{box-sizing:border-box}
body{
  margin:0; background:var(--ground); color:var(--ink);
  font:400 16px/1.65 Charter,"Iowan Old Style","Source Serif 4",Georgia,serif;
  -webkit-font-smoothing:antialiased;
}
.wrap{max-width:62rem;margin:0 auto;padding:4rem 1.5rem 6rem}
.meta{
  font-family:system-ui,-apple-system,"Segoe UI",sans-serif;
  font-size:.72rem;letter-spacing:.09em;text-transform:uppercase;
  color:var(--muted);
}
h1{font-size:2.15rem;line-height:1.15;margin:.5rem 0 0;font-weight:600;
   letter-spacing:-.015em;text-wrap:balance;max-width:24ch}
.lede{color:var(--ink-2);font-size:1.075rem;margin:1rem 0 0;max-width:62ch}
hr{border:0;border-top:1px solid var(--rule);margin:3rem 0}
figure{margin:0 0 3.5rem}
.plate{
  background:var(--plate);border:1px solid var(--rule);border-radius:3px;
  box-shadow:var(--shadow);padding:.75rem;overflow-x:auto;
}
.plate img{display:block;width:100%;height:auto;min-width:34rem}
figcaption{margin-top:.95rem;max-width:64ch}
.fignum{
  font-family:system-ui,-apple-system,"Segoe UI",sans-serif;
  font-size:.7rem;letter-spacing:.1em;text-transform:uppercase;
  color:var(--accent);font-weight:600;
}
.figttl{font-weight:600;margin:.2rem 0 .35rem;font-size:1.05rem}
figcaption p{margin:0;color:var(--ink-2);font-size:.95rem;line-height:1.6}
h2{font-size:1.4rem;font-weight:600;margin:0 0 .35rem;letter-spacing:-.01em}
h2+p{margin:0 0 1.5rem;color:var(--ink-2);max-width:64ch;font-size:.95rem}
.tablewrap{overflow-x:auto;border:1px solid var(--rule);border-radius:3px;
           background:var(--surface)}
table{border-collapse:collapse;width:100%;min-width:52rem;
      font-family:system-ui,-apple-system,"Segoe UI",sans-serif;font-size:.8rem}
th,td{padding:.42rem .7rem;text-align:right;white-space:nowrap;
      font-variant-numeric:tabular-nums}
th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left;
      font-variant-numeric:normal}
thead th{
  position:sticky;top:0;background:var(--surface);color:var(--muted);
  font-weight:600;font-size:.68rem;letter-spacing:.06em;text-transform:uppercase;
  border-bottom:1px solid var(--rule);
}
tbody tr:not(:last-child) td{border-bottom:1px solid var(--rule)}
tr.grp td{border-top:2px solid var(--rule)}
td.name{color:var(--ink);font-weight:500}
tr.base td{color:var(--ink-2)}
tr.base td.name::before{content:"○ ";color:var(--muted)}
.swatch{display:inline-block;width:.55rem;height:.55rem;border-radius:50%;
        margin-right:.45rem;vertical-align:baseline}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(21rem,1fr));
      gap:2.25rem 1.75rem;margin-top:1.5rem}
.grid.wide{grid-template-columns:1fr;gap:2.5rem}
.grid figure{margin:0}
.grid .plate{padding:.5rem}
.grid .plate img{min-width:0}
.grid figcaption{margin-top:.6rem}
.gcap{font-weight:600;font-size:.95rem}
.gnum{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;font-size:.75rem;
      color:var(--ink-2);margin-top:.15rem;font-variant-numeric:tabular-nums}
.note{color:var(--muted);font-size:.82rem;margin-top:.8rem;
      font-family:system-ui,-apple-system,"Segoe UI",sans-serif;max-width:70ch}
code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.85em;
     background:var(--surface);border:1px solid var(--rule);border-radius:3px;
     padding:.08em .34em}
"""

SERIES = {"evo2": "#2a78d6", "carbon": "#eb6834", "dnabert": "#1baf7a"}
LABEL = {"evo2": "Evo 2 7B", "carbon": "Carbon 8B", "dnabert": "DNABERT-2"}
BASELINE = "TB-Profiler 6.3.0"


def table_html(summary, split: str) -> str:
    drugs = sorted({r["drug"] for r in summary})
    out = ['<div class="tablewrap"><table><thead><tr>',
           "<th>Drug</th><th>Model</th>"]
    out += [f"<th>{h}</th>" for _, h in FIELDS]
    out.append("</tr></thead><tbody>")
    for d in drugs:
        rows = [r for r in summary if r["drug"] == d and r["split"] == split]
        models = [m for m in ("evo2", "carbon", "dnabert") if
                  any(r["model"] == m for r in rows)]
        seq = models + [m for m in sorted({r["model"] for r in rows})
                        if m not in models]
        for i, m in enumerate(seq):
            r = next(rr for rr in rows if rr["model"] == m)
            cls = []
            if i == 0:
                cls.append("grp")
            if m not in SERIES:
                cls.append("base")
            dot = (f'<span class="swatch" style="background:{SERIES[m]}"></span>'
                   if m in SERIES else "")
            name = LABEL.get(m, m)
            out.append(f'<tr class="{" ".join(cls)}">')
            out.append(f'<td class="name">{html.escape(d.capitalize()) if i == 0 else ""}</td>')
            out.append(f'<td class="name">{dot}{html.escape(name)}</td>')
            for f, _ in FIELDS:
                v = r.get(f)
                out.append("<td>" + ("—" if v in ("", None)
                                     else f"{float(v):.3f}") + "</td>")
            out.append("</tr>")
    out.append("</tbody></table></div>")
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workbook", required=True)
    ap.add_argument("--fig-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sys.path.insert(0, str(pathlib.Path(__file__).parent))
    from make_panel_figures import load
    summary, _ = load(pathlib.Path(args.workbook))

    figdir = pathlib.Path(args.fig_dir)
    parts = [
        "<title>Ten-Drug MTB Panel</title>",
        f"<style>{CSS}</style>",
        '<div class="wrap">',
        '<div class="meta">M. tuberculosis · MIC prediction from '
        'foundation-model embeddings</div>',
        "<h1>Ten-drug panel, figure set</h1>",
        '<p class="lede">Embeddings of the catalogue\'s own loci, scored against '
        'TB-Profiler 6.3.0 on a held-out set of roughly 1,300–1,800 isolates per '
        'drug. Every number below is read from the results workbook by the same '
        'script that draws the figures.</p>',
        "<hr>",
    ]
    def render(fname, num, title, cap):
        if not (figdir / fname).exists():
            raise SystemExit(f"[ERROR] missing figure {figdir/fname}")
        return [
            "<figure>",
            f'<div class="plate"><img src="{fname}" alt="{html.escape(title)}"></div>',
            "<figcaption>",
            f'<div class="fignum">{num}</div>',
            f'<div class="figttl">{html.escape(title)}</div>',
            f"<p>{html.escape(cap)}</p>",
            "</figcaption></figure>",
        ]

    for spec in FIGURES:
        parts += render(*spec)

    def g(drug, field):
        for r in summary:
            if (r["drug"] == drug and r["model"] == "evo2"
                    and r["split"] == "holdout"):
                v = r.get(field)
                return None if v in ("", None) else float(v)
        return None

    parts += ["<hr>",
              '<div class="fignum">Figure 5</div>',
              "<h2>MIC placement, every drug</h2>",
              "<p>Evo 2 on the held-out set, one panel per drug. Each figure is "
              "redrawn from its own saved predictions and refuses to render "
              "unless the metrics recomputed from them reproduce that drug's row "
              "in the table below.</p>",
              '<div class="grid">']
    for drug in SCATTER_ORDER:
        f = f"scatter_{drug}.png"
        if not (figdir / f).exists():
            raise SystemExit(f"[ERROR] missing {figdir/f}")
        stats = "  ·  ".join(
            f"{lab} {g(drug, fld):.3f}" for fld, lab in
            (("exact", "exact"), ("OT", "OT"), ("vMER", "vMER"))
            if g(drug, fld) is not None)
        parts += ["<figure>",
                  f'<div class="plate"><img src="{f}" alt="{drug} scatter"></div>',
                  "<figcaption>",
                  f'<div class="gcap">{drug.capitalize()}</div>',
                  f'<div class="gnum">{stats}</div>',
                  "</figcaption></figure>"]
    parts += ["</div>", "<hr>"]
    parts += render("sheet_scatters.png", "Figure 5 — sheet",
                    "All ten confusion grids on one page",
                    "The same panels as above, at a size that lets the panel be "
                    "read as a panel. Counts inside the bubbles are dropped — "
                    "they do not fit — so this version carries the size key the "
                    "individual figures do not need.")
    parts += ["<hr>",
              '<div class="fignum">Figure 6</div>',
              "<h2>Class sizes, and where accuracy follows them</h2>",
              "<p>The same held-out predictions counted by MIC class, against "
              "the development cohort behind each class. The averages every "
              "other figure reports are taken over these very unequal groups — "
              "<code>bal_OT</code> weights them equally, <code>OT</code> does "
              "not, and the difference between the two is visible here as the "
              "line dropping wherever the grey bar is short.</p>",
              '<div class="grid wide">']
    for drug in SCATTER_ORDER:
        f = f"classdist_{drug}.png"
        if not (figdir / f).exists():
            raise SystemExit(f"[ERROR] missing {figdir/f}")
        parts += ["<figure>",
                  f'<div class="plate"><img src="{f}" alt="{drug} class distribution"></div>',
                  "<figcaption>",
                  f'<div class="gcap">{drug.capitalize()}</div>',
                  f'<div class="gnum">bal_OT {g(drug, "bal_OT"):.3f}  ·  '
                  f'OT {g(drug, "OT"):.3f}</div>',
                  "</figcaption></figure>"]
    parts += ["</div>", "<hr>"]
    parts += ['<div class="fignum">Figure 7</div>',
              "<h2>What the model used</h2>",
              "<p>TreeSHAP on each reported model, one point per held-out "
              "isolate, summed over a region's four principal components. Every "
              "drug's top locus is the one its catalogue names — and for "
              "isoniazid and ethionamide the signal sits in the "
              "<em>promoter</em> windows while the gene bodies stay at zero, "
              "which is the mechanism, recovered from sequence without the "
              "model being told where to look.</p>",
              '<div class="grid wide">']
    for drug in SCATTER_ORDER:
        f = f"shap_{drug}.png"
        if not (figdir / f).exists():
            raise SystemExit(f"[ERROR] missing {figdir/f}")
        parts += ["<figure>",
                  f'<div class="plate"><img src="{f}" alt="{drug} SHAP"></div>',
                  "<figcaption>",
                  f'<div class="gcap">{drug.capitalize()}</div>',
                  "</figcaption></figure>"]
    parts += ["</div>", "<hr>"]
    parts += render("sheet_classdist.png", "Figure 6 — sheet",
                    "All ten class distributions on one page",
                    "Read down the column: the accuracy line falls wherever the "
                    "grey bar is short, on every drug. That is the gap between "
                    "OT and bal_OT — the first counts isolates, the second "
                    "counts classes.")
    parts += ["<hr>"]

    for spec in SUPPLEMENTARY:
        parts += render(*spec)

    parts += [
        "<hr>",
        "<h2>Held-out set, every row</h2>",
        "<p>The table is the accessible form of the figures: three of the series "
        "colours sit below a 3:1 contrast ratio on white, so identity is never "
        "left to colour alone.</p>",
        table_html(summary, "holdout"),
        '<p class="note">vMER = very major error rate (1 − sensitivity); '
        'MER = major error rate (1 − specificity); cost = 2·vMER + MER; '
        'OT = within one dilution; AUROC is computed on the continuous '
        'pre-snapping prediction. Catalogue rows have no MIC and therefore no '
        'OT, exact or AUROC.</p>',
        "<hr>",
        "<h2>Development, 3×5-fold CV</h2>",
        "<p>The same rows on the development folds, for the comparison figure S1 "
        "draws.</p>",
        table_html(summary, "development_cv"),
        '<p class="note">Regenerate with '
        '<code>python3 mtb_tbprofiler/make_panel_figures.py --out-dir '
        'figures/mtb</code>, then '
        '<code>build_figure_page.py</code>. PDFs of every figure sit beside the '
        'PNGs in the same directory.</p>',
        "</div>",
    ]
    out = pathlib.Path(args.out)
    out.write_text("\n".join(parts), encoding="utf-8")
    print(f"[OK] {out}  ({len(out.read_text()):,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
