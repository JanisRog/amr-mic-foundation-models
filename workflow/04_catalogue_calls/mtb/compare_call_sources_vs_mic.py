#!/usr/bin/env python3
"""TB-Profiler's calls or Snippy's consensus — which one gets the MIC right?

Two variant call sets exist for every isolate and they disagree. TB-Profiler
reports down to an allele fraction of 0.1; Snippy keeps only what freebayes
genotypes as homozygous alt under a ploidy-2 model, and its consensus is what
every embedding in this project is built from. Section 7
measured the disagreement at 3.77% of catalogued-mutation carriers but never
asked which side is right, and that is the question that decides whether a
targeted re-call is worth building.

The test is the phenotype. For each drug, an isolate is called resistant by
TB-Profiler if it carries any catalogued mutation for that drug, and by Snippy if
that same mutation also survives into Snippy's call set. The isolates where the
two differ are the whole experiment:

  * if TB-Profiler-only isolates are mostly truly resistant, the low-fraction call
    is real, the consensus destroys it, and the embeddings are structurally blind
    to those isolates -- the re-call and a mixedness feature are justified
  * if they look like the wild-type group, Snippy's threshold is filtering noise
    or a minority population without phenotypic effect, and the consensus is the
    better representation -- nothing needs rebuilding

Matching is on (position, ref, alt) after both sides are reduced to the minimal
VCF representation, because the two callers pad indels differently and an
unmatched representation must not be reported as consensus loss. Anything still
unmatched is retried as a positional overlap and counted separately. The three
unmatched-by-any-window columns from section 7 are large indels and are expected
to show up here as genuinely missing from Snippy too.

Runs on the login node in a few minutes for a subset; use the sbatch for all
10,170 samples.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import csv
import json
import math
import pathlib
import statistics
import sys

# TB-Profiler drug names -> the metadata workbook's abbreviations
DRUG_ABBREV = {
    "amikacin": "AMI", "bedaquiline": "BDQ", "clofazimine": "CFZ",
    "delamanid": "DLM", "ethambutol": "EMB", "ethionamide": "ETH",
    "isoniazid": "INH", "kanamycin": "KAN", "levofloxacin": "LEV",
    "linezolid": "LZD", "moxifloxacin": "MXF", "rifabutin": "RFB",
    "rifampicin": "RIF",
}

# Section 3 plus the three the RF500 CV package carried and section 3 did not
# (amikacin >0, clofazimine >2, linezolid >0), on the integer-rounded log2 ladder.
BOUNDARY = {"RIF": -1.0, "INH": -3.0, "EMB": 2.0, "ETH": 2.0, "KAN": 2.0,
            "LEV": 0.0, "MXF": 0.0, "BDQ": -2.0, "DLM": -3.0, "RFB": -3.0,
            "AMI": 0.0, "CFZ": 2.0, "LZD": 0.0}


def first(d: dict, *keys, default=None):
    for k in keys:
        if k in d and d[k] not in (None, ""):
            return d[k]
    return default


def variant_drugs(v: dict) -> set[str]:
    """TB-Profiler has used several shapes for the drug annotation."""
    out: set[str] = set()
    raw = first(v, "drugs", "gene_associated_drugs", default=[])
    if isinstance(raw, str):
        raw = [x for x in raw.replace(";", ",").split(",") if x.strip()]
    for d in raw:
        name = d.get("drug") if isinstance(d, dict) else d
        if name:
            out.add(str(name).strip().lower())
    return out


def atomize(pos: int, ref: str, alt: str) -> list[tuple[int, str, str]]:
    """Split an equal-length MNP into the substitutions it is made of.

    freebayes emits one haplotype call for SNPs close enough together, so an
    isolate with gyrA D94G (7582) and the common lineage polymorphism S95T (7585)
    gets `7582 complex ACAG>GCAC` while TB-Profiler reports `7582 A>G`. Those are
    the same variant written at different granularity, and without splitting them
    every codon-94 carrier reads as a Snippy miss -- which is what made Snippy's
    fluoroquinolone sensitivity look like 0.31.

    Only equal-length ref/alt are split, where each differing offset is an
    independent substitution. Anything involving an indel is left whole, since
    its decomposition is not unique.
    """
    if len(ref) != len(alt) or len(ref) == 1:
        return []
    return [(pos + i, r, a) for i, (r, a) in enumerate(zip(ref, alt)) if r != a]


def norm(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Minimal VCF representation: trim the shared suffix, then the shared prefix.

    Without this, an insertion TB-Profiler writes as 4326100 A>AG and Snippy
    writes as 4326099 CA>CAG never matches, and a representational difference
    gets counted as consensus loss. Both collapse to 4326100 A>AG here.
    """
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    i = 0
    while i < len(ref) - 1 and i < len(alt) - 1 and ref[i] == alt[i]:
        i += 1
    return pos + i, ref[i:], alt[i:]


def load_snippy(path: pathlib.Path) -> tuple[set[tuple[int, str, str]], list[tuple[int, int]], int]:
    """(pos, ref, alt) triples, their spans, and the widest span."""
    exact: set[tuple[int, str, str]] = set()
    spans: list[tuple[int, int]] = []
    if path.name.endswith(".tab"):
        with path.open() as fh:
            rd = csv.reader(fh, delimiter="\t")
            head = next(rd, None)
            if not head:
                return exact, spans, 0
            col = {c.upper(): i for i, c in enumerate(head)}
            ip, ir, ia = col.get("POS"), col.get("REF"), col.get("ALT")
            if None in (ip, ir, ia):
                raise SystemExit(f"[ERROR] unexpected snps.tab header in {path}: {head}")
            for r in rd:
                if len(r) <= max(ip, ir, ia):
                    continue
                try:
                    pos = int(r[ip])
                except ValueError:
                    continue
                ref, alt = r[ir].upper(), r[ia].upper()
                spans.append((pos, pos + max(len(ref), 1) - 1))
                exact.add(norm(pos, ref, alt))
                exact.update(atomize(pos, ref, alt))
    else:
        with path.open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                try:
                    pos = int(f[1])
                except ValueError:
                    continue
                ref = f[3].upper()
                for alt in f[4].upper().split(","):
                    exact.add(norm(pos, ref, alt))
                    exact.update(atomize(pos, ref, alt))
                spans.append((pos, pos + max(len(ref), 1) - 1))
    spans.sort()
    return exact, spans, max((e - s for s, e in spans), default=0)


def overlaps(spans: list[tuple[int, int]], widest: int, lo: int, hi: int) -> bool:
    """Any called variant touching [lo, hi]. Spans are sorted by start.

    A long deletion starts far to the left of the position being asked about, so
    the scan begins one max-span back rather than at `lo` — with a fixed lookback
    a 1.6 kb ethA deletion would be reported as missing when Snippy did call it.
    """
    if not spans:
        return False
    i = bisect.bisect_left(spans, (lo - widest - 1, 0))
    for s, e in spans[i:]:
        if s > hi:
            return False
        if e >= lo:
            return True
    return False


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar on the discordant pairs only.

    The concordant cells carry no information about which caller is better, so
    the test asks whether b of b+c disagreements favouring one side is more
    lopsided than a coin would give.
    """
    n = b + c
    if n == 0:
        return float("nan")
    lo = min(b, c)
    tail = sum(math.comb(n, i) for i in range(lo + 1)) / (2 ** n)
    return min(1.0, 2 * tail)


def metrics(pred_r: list[bool], true_r: list[bool]) -> dict:
    tp = sum(p and t for p, t in zip(pred_r, true_r))
    fp = sum(p and not t for p, t in zip(pred_r, true_r))
    fn = sum((not p) and t for p, t in zip(pred_r, true_r))
    tn = sum((not p) and (not t) for p, t in zip(pred_r, true_r))
    sens = tp / max(tp + fn, 1)
    spec = tn / max(tn + fp, 1)
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "sens": sens, "spec": spec,
            "vMER": fn / max(tp + fn, 1), "MER": fp / max(tn + fp, 1),
            "bal_acc": (sens + spec) / 2,
            "f1": 2 * tp / max(2 * tp + fp + fn, 1)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tbprofiler-json-dir", required=True,
                    help=".../results/tbprofiler/reads/gather/results")
    ap.add_argument("--snippy-runs", required=True,
                    help=".../results/gwas/snippy_h37rv/runs")
    ap.add_argument("--mic-tsv", required=True, help="mtb_all_drug_log2_mic.tsv")
    ap.add_argument("--mic-suffix", default="_int",
                    help="'_int' for the rounded ladder the runs use, '' for raw")
    ap.add_argument("--boundary", action="append", default=[], metavar="DRUG=VAL")
    ap.add_argument("--drug-alias", action="append", default=[], metavar="DRUG=FROM",
                    help="score DRUG using FROM's calls. TB-Profiler annotates no "
                         "rifabutin variants in this catalogue, so RFB=RIF reads "
                         "rifabutin off the rpoB calls -- an approximation, since "
                         "some rpoB alleles keep rifabutin susceptibility")
    ap.add_argument("--min-freq", type=float, default=0.0,
                    help="ignore TB-Profiler calls below this allele fraction")
    ap.add_argument("--limit", type=int, default=0, help="first N samples, for a smoke test")
    ap.add_argument("--out-calls", default="", help="per-sample call table")
    ap.add_argument("--out-variants", default="",
                    help="per-variant match status, for chasing specific mutations")
    ap.add_argument("--out-discordant", default="",
                    help="isolates where exactly one caller gets the S/R call right")
    ap.add_argument("--skip-drugs", default="RFB",
                    help="comma-separated; this TB-Profiler catalogue annotates no "
                         "rifabutin variants, so RFB is zero against zero")
    ap.add_argument("--snippy-id-map", default="",
                    help="TSV mapping a known sample id to the Snippy run directory "
                         "name, when the run dirs are keyed by something the metadata "
                         "does not carry (an ENA run accession, say)")
    args = ap.parse_args()

    skip = {d.strip().upper() for d in args.skip_drugs.split(",") if d.strip()}
    bounds = dict(BOUNDARY)
    for spec in args.boundary:
        k, _, v = spec.partition("=")
        bounds[k.strip().upper()] = float(v)
    aliases = {}
    for spec in args.drug_alias:
        k, _, v = spec.partition("=")
        aliases[k.strip().upper()] = v.strip().upper()
    if aliases:
        print(f"[INFO] drug aliases: {aliases}")

    # MIC table, indexed by both id forms so the JSON naming does not have to be guessed
    mic: dict[str, dict[str, float]] = {}
    alias: dict[str, str] = {}
    names: dict[str, list[str]] = {}      # sample -> every name its Snippy dir might use
    with open(args.mic_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid = r["sample_id"]
            vals = {}
            for d in DRUG_ABBREV.values():
                raw = r.get(d + args.mic_suffix, "")
                if raw not in ("", "NA"):
                    vals[d] = float(raw)
            mic[sid] = vals
            alias[sid] = sid
            names.setdefault(sid, [sid])
            # subject is deliberately not an alias: several isolates share one
            # subject, so it would map different samples onto the same row
            if r.get("ENA_SAMPLE"):
                alias.setdefault(r["ENA_SAMPLE"], sid)
                names[sid].append(r["ENA_SAMPLE"])
    if args.snippy_id_map:
        extra = 0
        with open(args.snippy_id_map) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2:
                    continue
                sid = alias.get(f[0].strip())
                if sid:
                    names[sid].append(f[1].strip())
                    extra += 1
        print(f"[INFO] {extra:,} directory names from {args.snippy_id_map}")
    print(f"[INFO] MIC for {len(mic):,} samples, {len(alias):,} id aliases")

    jdir = pathlib.Path(args.tbprofiler_json_dir)
    sdir = pathlib.Path(args.snippy_runs)
    jsons = sorted(jdir.glob("*.results.json"))
    if args.limit:
        jsons = jsons[:args.limit]
    print(f"[INFO] {len(jsons):,} TB-Profiler results")

    # Pass 1: what TB-Profiler called, and the working catalogue.
    #
    # Snippy's side has to be scored independently or the comparison is rigged:
    # walking only TB-Profiler's list can never find a variant Snippy called and
    # TB-Profiler did not, so "Snippy-only" would be zero by construction rather
    # than by evidence. The union of every drug-annotated variant TB-Profiler
    # reports anywhere in the cohort is the catalogue both callers are scored
    # against. It is a lower bound on the real catalogue -- a variant nobody
    # carries is not in it -- which is the right bound here, since a variant no
    # isolate carries cannot change any isolate's call.
    dr_by_sample: dict[str, list[tuple]] = {}
    catalogue: dict[tuple[int, str, str], set[str]] = {}
    no_mic = 0
    for n, jf in enumerate(jsons, 1):
        if n % 2000 == 0:
            print(f"[INFO] pass 1 {n:,}/{len(jsons):,}", flush=True)
        stem = jf.name[: -len(".results.json")]
        sid = alias.get(stem)
        if sid is None:
            no_mic += 1
            continue
        try:
            obj = json.loads(jf.read_text())
        except (OSError, json.JSONDecodeError) as e:
            print(f"[WARN] unreadable {jf.name}: {e}")
            continue
        recs = []
        for v in obj.get("dr_variants", []):
            drugs = {DRUG_ABBREV[d] for d in variant_drugs(v) if d in DRUG_ABBREV}
            if not drugs:
                continue
            try:
                freq = float(first(v, "freq", "frequency", default=1.0))
            except (TypeError, ValueError):
                freq = 1.0
            if freq < args.min_freq:
                continue
            try:
                pos = int(first(v, "genome_pos", "pos", "position", default=0))
            except (TypeError, ValueError):
                pos = 0
            ref = str(first(v, "ref", "ref_nt", default="")).upper()
            alt = str(first(v, "alt", "alt_nt", default="")).upper()
            gene = str(first(v, "gene", "gene_name", "locus_tag", default=""))
            change = str(first(v, "change", "variant", default=""))
            if pos and ref and alt:
                pos, ref, alt = norm(pos, ref, alt)
            recs.append((pos, ref, alt, freq, gene, change, frozenset(drugs)))
            if pos and ref and alt:
                catalogue.setdefault((pos, ref, alt), set()).update(drugs)
                for key in atomize(pos, ref, alt):
                    catalogue.setdefault(key, set()).update(drugs)
        dr_by_sample[stem] = recs
    print(f"[INFO] catalogue: {len(catalogue):,} distinct drug-annotated variants")

    # Pass 2: score both callers per sample, reading each Snippy file once.
    called: dict[str, dict[str, list[bool]]] = collections.defaultdict(dict)
    match_stat = collections.Counter()
    freq_by_status: dict[str, list[float]] = collections.defaultdict(list)
    snippy_only_var = collections.Counter()
    overlap_only: dict[str, set[str]] = collections.defaultdict(set)
    var_rows = []
    no_snippy = 0

    for n, (stem, recs) in enumerate(dr_by_sample.items(), 1):
        if n % 2000 == 0:
            print(f"[INFO] pass 2 {n:,}/{len(dr_by_sample):,}", flush=True)
        sid = alias[stem]
        snp_file = None
        cands = []
        for nm in dict.fromkeys([stem, *names.get(sid, [])]):
            for leaf in ("snps.tab", "snps.vcf"):
                cands.append(sdir / nm / leaf)
        for cand in cands:
            if cand.exists():
                snp_file = cand
                break
        if snp_file is None:
            no_snippy += 1
            if no_snippy == 1:
                print(f"[WARN] no Snippy run for {sid}. Tried:")
                for c in cands[:6]:
                    print(f"         {c}")
                have = sorted(x.name for x in list(sdir.iterdir())[:2000] if x.is_dir())
                print(f"[WARN] {sdir} holds {len(have):,}+ directories, e.g.:")
                for h in have[:5]:
                    print(f"         {h}")
                print("[WARN] if those names match neither form, pass --snippy-id-map.")
            continue
        exact, spans, widest = load_snippy(snp_file)

        for pos, ref, alt, freq, gene, change, drugs in recs:
            if pos and ref and alt and (pos, ref, alt) in exact:
                status = "exact"
            elif pos and overlaps(spans, widest, pos, pos + max(len(ref), 1) - 1):
                status = "overlap"
            elif not pos:
                status = "no_position"
            else:
                status = "missing"
            match_stat[status] += 1
            freq_by_status[status].append(freq)
            if args.out_variants:
                var_rows.append([sid, gene, change, pos, ref, alt, f"{freq:.3f}",
                                 status, ";".join(sorted(drugs))])
            for ab in drugs:
                called[ab].setdefault(sid, [False, False])[0] = True
                if status == "overlap":
                    # matched only by position, not by representation: not counted
                    # for Snippy, but tracked so the caveat can be sized
                    overlap_only[ab].add(sid)

        # Snippy's own reading of the catalogue, and nothing else. Crediting
        # Snippy for variants TB-Profiler flagged would make its call partly
        # derived from TB-Profiler's, which is the one thing a head-to-head
        # cannot afford.
        own = {(p, r, a) for p, r, a, *_ in recs}
        for key in exact & catalogue.keys():
            if key not in own:
                snippy_only_var[f"{key[0]} {key[1]}>{key[2]}"] += 1
            for ab in catalogue[key]:
                called[ab].setdefault(sid, [False, False])[1] = True

    print(f"\n[INFO] {no_mic:,} results with no MIC row, {no_snippy:,} with no Snippy run")
    tot = sum(match_stat.values())
    print(f"\nCatalogued calls matched against Snippy ({tot:,} calls)")
    print(f"{'status':<12}{'n':>9}{'share':>8}{'median freq':>13}")
    for st in ("exact", "overlap", "missing", "no_position"):
        if not match_stat[st]:
            continue
        print(f"{st:<12}{match_stat[st]:>9,}{match_stat[st]/tot:>8.1%}"
              f"{statistics.median(freq_by_status[st]):>13.3f}")
    print("'overlap' is a representational difference (indel padding), not a loss.\n"
          "'missing' is the consensus loss: TB-Profiler called it, Snippy did not.")

    if snippy_only_var:
        n_so = sum(snippy_only_var.values())
        print(f"\nThe other direction: {n_so:,} calls where Snippy carries a catalogued\n"
              f"variant TB-Profiler did not report for that sample "
              f"({len(snippy_only_var):,} distinct).")
        print(f"{'variant':<28}{'samples':>9}")
        for v, c in snippy_only_var.most_common(15):
            print(f"{v:<28}{c:>9,}")

    # ---- per drug -------------------------------------------------------
    print(f"\n{'drug':<6}{'n':>7}{'%R':>7}{'TBP+':>7}{'SNP+':>7}{'TBPonly':>9}"
          f"{'SNPonly':>9}  | TB-Profiler sens/spec/bal | Snippy sens/spec/bal")
    cells: dict[str, dict[str, list[float]]] = {}
    correctness: dict[str, tuple] = {}
    for ab in sorted(DRUG_ABBREV.values()):
        if ab in skip:
            continue
        b = bounds.get(ab)
        samples = [s for s in mic if ab in mic[s]]
        if not samples:
            continue
        src = aliases.get(ab, ab)
        tbp = [called[src].get(s, [False, False])[0] for s in samples]
        snp = [called[src].get(s, [False, False])[1] for s in samples]
        tbp_only = sum(t and not s_ for t, s_ in zip(tbp, snp))
        snp_only = sum(s_ and not t for t, s_ in zip(tbp, snp))
        if b is None:
            print(f"{ab:<6}{len(samples):>7}{'—':>7}{sum(tbp):>7}{sum(snp):>7}"
                  f"{tbp_only:>9}{snp_only:>9}  | no S/R boundary — cells only")
            cells[ab] = None
            continue
        true_r = [mic[s][ab] > b for s in samples]
        correctness[ab] = ([t == r for t, r in zip(tbp, true_r)],
                           [x == r for x, r in zip(snp, true_r)],
                           samples, true_r, list(zip(tbp, snp)))
        mt = metrics(tbp, true_r)
        ms = metrics(snp, true_r)
        print(f"{ab:<6}{len(samples):>7}{sum(true_r)/len(samples):>7.1%}"
              f"{sum(tbp):>7}{sum(snp):>7}{tbp_only:>9}{snp_only:>9}  |"
              f" {mt['sens']:.3f} {mt['spec']:.3f} {mt['bal_acc']:.3f}          |"
              f" {ms['sens']:.3f} {ms['spec']:.3f} {ms['bal_acc']:.3f}")
        cells[ab] = {"both": [], "tbp_only": [], "snp_only": [], "neither": []}
        for s, t, sn in zip(samples, tbp, snp):
            key = ("both" if t and sn else "tbp_only" if t else
                   "snp_only" if sn else "neither")
            cells[ab][key].append(mic[s][ab])

    # ---- correctness head to head ---------------------------------------
    print("\nWho is right, isolate by isolate. Each caller is scored against the\n"
          "phenotype (MIC > boundary); the concordant cells carry no information\n"
          "about which caller is better, so the test is on the two disagreements.")
    print(f"\n{'drug':<5}{'n':>7}{'TBP acc':>9}{'SNP acc':>9}{'both ok':>9}"
          f"{'TBP only':>10}{'SNP only':>10}{'neither':>9}{'p':>10}")
    disc_rows = []
    tot = collections.Counter()
    for ab in sorted(DRUG_ABBREV.values()):
        if ab in skip or ab not in correctness:
            continue
        t_ok, s_ok, ids, true_r_by, call_by = correctness[ab]
        a = int(sum(t and x for t, x in zip(t_ok, s_ok)))
        bb = int(sum(t and not x for t, x in zip(t_ok, s_ok)))
        cc = int(sum((not t) and x for t, x in zip(t_ok, s_ok)))
        dd = int(sum((not t) and (not x) for t, x in zip(t_ok, s_ok)))
        n = a + bb + cc + dd
        tot["a"] += a; tot["b"] += bb; tot["c"] += cc; tot["d"] += dd
        print(f"{ab:<5}{n:>7}{(a+bb)/n:>9.3f}{(a+cc)/n:>9.3f}{a:>9}"
              f"{bb:>10}{cc:>10}{dd:>9}{mcnemar_exact(bb, cc):>10.2g}")
        for i, sid in enumerate(ids):
            if t_ok[i] == s_ok[i]:
                continue
            disc_rows.append([sid, ab, f"{mic[sid][ab]:g}", f"{bounds[ab]:g}",
                              int(true_r_by[i]), int(call_by[i][0]),
                              int(call_by[i][1]),
                              "tbprofiler" if t_ok[i] else "snippy"])
    n = sum(tot.values())
    print(f"\n{'pooled':<5}{n:>7}{(tot['a']+tot['b'])/n:>9.3f}"
          f"{(tot['a']+tot['c'])/n:>9.3f}{tot['a']:>9}{tot['b']:>10}"
          f"{tot['c']:>10}{tot['d']:>9}{mcnemar_exact(tot['b'], tot['c']):>10.2g}")
    print("\n'TBP only' = TB-Profiler right and Snippy wrong; 'SNP only' the reverse.")
    n_ov = sum(len(v) for k, v in overlap_only.items() if k not in skip)
    if n_ov:
        print(f"Caveat: {n_ov:,} isolate-drug pairs where Snippy has a variant at the\n"
              f"right position but written differently (indel padding). Those are not\n"
              f"credited to Snippy here, so its column is a lower bound by at most that.")
    print("Caveat: the catalogue both callers are scored against is TB-Profiler's own\n"
          "output pooled over the cohort, so a variant TB-Profiler never calls anywhere\n"
          "is invisible to both. The comparison is not perfectly symmetric.")

    # ---- the decisive table --------------------------------------------
    print("\nMIC by agreement cell — the isolates the two callers disagree on are the test.")
    print("If TB-Profiler-only sits with 'both', the low-fraction call is real and the")
    print("consensus loses it. If it sits with 'neither', Snippy is filtering noise.")
    print(f"\n{'drug':<6}{'cell':<10}{'n':>7}{'median':>9}{'%R':>8}")
    for ab in sorted(cells):
        if cells[ab] is None:
            continue
        b = bounds[ab]
        for key in ("both", "tbp_only", "snp_only", "neither"):
            v = cells[ab][key]
            if not v:
                continue
            print(f"{ab:<6}{key:<10}{len(v):>7}{statistics.median(v):>9.1f}"
                  f"{sum(1 for x in v if x > b)/len(v):>8.1%}")

    if args.out_discordant:
        p = pathlib.Path(args.out_discordant)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sample_id", "drug", "log2_mic", "boundary", "true_resistant",
                        "tbprofiler_call", "snippy_call", "correct_caller"])
            w.writerows(disc_rows)
        print(f"\n[OK] {len(disc_rows):,} discordant isolate-drug pairs -> {p}")

    if args.out_calls:
        p = pathlib.Path(args.out_calls)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            cols = sorted(DRUG_ABBREV.values())
            w.writerow(["sample_id"] + [f"{c}_{k}" for c in cols for k in ("tbp", "snippy")])
            for s in sorted(mic):
                row = [s]
                for c in cols:
                    t, sn = called[c].get(s, [False, False])
                    row += [int(t), int(sn)]
                w.writerow(row)
        print(f"\n[OK] per-sample calls -> {p}")

    if args.out_variants:
        p = pathlib.Path(args.out_variants)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sample_id", "gene", "change", "genome_pos", "ref", "alt",
                        "freq", "snippy_status", "drugs"])
            w.writerows(var_rows)
        print(f"[OK] {len(var_rows):,} variant rows -> {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
