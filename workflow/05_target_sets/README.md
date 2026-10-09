# 05 — Target sets

Each locus contributes two targets: the gene body and a 300 bp strand-aware
upstream interval.

## *M. tuberculosis*

1. `build_mtb_all_gene_target_clumps.py` — every gene plus upstream flank in the
   H37Rv annotation (8,243 targets).
2. `build_panel_target_sets.py` — per-drug candidate sets from the published
   association hits (stage 01) and the catalogue's own annotation (`tbdb.bed`):

   | set | contents |
   |---|---|
   | `curatedmech` | association hits the catalogue annotates for this drug |
   | `catgenes` | loci the catalogue annotates for this drug |
   | **`both`** | `curatedmech` ∪ `catgenes` — **the reported set** |
   | `gwasfiltered` | `both` ∪ unannotated hits — the wider set evaluated |
   | `other` | hits annotated for another drug only — never used |

   A locus annotated for a pharmacological relative (inhA for isoniazid and
   ethionamide) is not counted as foreign. An unannotated locus sharing a curated
   intergenic entry with, or within 500 bp of, a mechanism locus is folded into it
   (this keeps *fabG1*). `--max-intergenic-span` (default 20 loci) stops a
   malformed intergenic entry from spanning the chromosome.
3. `subset_target_clumps.py` — cut the all-gene clumps to the union of the panel's
   sets, so extraction runs once.

`read_curated_drug_genes.py` reads one drug's entries without openpyxl.

## *P. aeruginosa*

`build_pseaer_target_sets.py` resolves the curated loci against the PAO1
annotation. Each locus is declared with an explicit rule (gene name, prefix,
product pattern or locus tag), and a rule matching other than exactly one coding
sequence is an error — required for *ampC* (annotated *blaPDC-1*), *mexZ* and
*oprD*. `filter_amrseq_targets.py` drops acquired-gene targets whose gene is
already a curated target.
