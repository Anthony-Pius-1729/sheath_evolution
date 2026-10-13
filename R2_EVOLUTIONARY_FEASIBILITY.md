# R2 pyocin — evolutionary feasibility report

**Date** 2026-10-13 · **Reference** UniProt Q9S574 / PDB 3J9Q, *Pseudomonas
aeruginosa*, 386 aa · **Structural sites** the 143 reference positions carrying a
frozen `R_i` value in `sheath_pilot/qc_output/R_r2pyocin.csv`

## Verdict

| | verdict |
|---|---|
| Orthology boundary | **UNRESOLVED** |
| Evolutionary feasibility for site-rate estimation | **CONDITIONAL** — fails as currently constituted; one concrete test could change that |

**Orthology: unresolved.** Four independent lines of evidence were examined and
three cannot discriminate at all. Domain architecture, protein names and
genomic-neighbourhood marker categories are effectively constant across the
whole retrieval set — every clade, including the most distant, looks like a
contractile-tail sheath in a contractile-tail locus. Only taxonomy and tree
support discriminate, and they place the boundary at a *Pseudomonas*-restricted
clade of 124 sequences. That boundary is not defensible as an ortholog set
either, because 39 of the 57 same-species paralogs that dereplication discarded
fall **inside** it, interleaved rather than monophyletic.

**Feasibility: conditional, and the two candidate routes fail for opposite
reasons.** No boundary is simultaneously clean and informative:

- **Narrow (124–155 sequences).** Alignment is clean and occupancy is perfect,
  but the set is one genus (1.05–2.57 effective genera, 79–99% *Pseudomonas*)
  and nearly invariant at the sites that matter: 43–46 of 143 structural sites
  have a single amino acid, 73–82 have one residue in ≥95% of sequences, median
  2 distinct residues per site, median column entropy 0.19–0.29 bits. This
  route is a **hard fail** and cannot be rescued by adding sequences from the
  same clade — the clade *is* one genus.
- **Wide (293–355 sequences).** Substitution signal is ample (6 and 1 invariant
  sites, median 9–11 residues, entropy 1.50–1.79 bits) and taxonomy is
  independent (21–35 effective genera), but every discarded paralog tested falls
  inside the clade, the alignment inflates to 2.2–2.7 columns per reference
  residue, and mean patristic distance reaches 2.08–2.56 substitutions/site with
  a maximum of 5.0–6.7 — saturation territory for site-rate estimation.

The condition is specific and testable: **positional orthology assigned from
locus gene order**, rather than from the generic marker-category score used
now. If a wide, taxonomically diverse set can be reduced to one positionally
orthologous copy per genome, R2 passes. If gene order cannot separate the
copies, R2 fails and the problem is unlikely to be R2-specific.

## 1. Clade structure and what each evidence line says

Nested sister groups outwards from the reference in the 968-taxon preliminary
tree (`orthology/r2pyocin_groups.csv` has all columns):

| group | n | cumulative | support | identity median | genera | top genera | arch. identical to reference | synteny supported | insertion residues (median) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 7 | 7 | 0.889 | 0.795 | 1 | *Pseudomonas* | 100% | 7/7 | 0 |
| 1 | 112 | 119 | 0.926 | 0.751 | 1 | *Pseudomonas* | 30% | 39/42 | 0 |
| 2 | 5 | 124 | 0.896 | 0.686 | 2 | *Pseudomonas*, *Frateuria* | 100% | 5/5 | 0 |
| 3 | 22 | 146 | **0.401** | 0.714 | 5 | *Burkholderia*, *Ralstonia* | 59% | 19/22 | 0 |
| 4 | 9 | 155 | 0.681 | 0.684 | 6 | *Pseudoduganella*, *Chromobacterium* | 89% | 5/9 | 1 |
| 5 | 3 | 158 | 1.000 | 0.695 | 3 | *Stutzerimonas*, … | 67% | 3/3 | 0 |
| 6 | 3 | 161 | 0.850 | 0.322 | 3 | *Sinobacterium*, *Vibrio* | 0% | 1/3 | 2 |
| 7 | 132 | 293 | 0.887 | 0.422 | 74 | *Pseudoalteromonas*, *Vibrio* | 61% | 76/84 | 7 |
| 8 | 5 | 298 | 0.534 | 0.344 | 3 | *Vibrio*, *Saliniradius* | 0% | 5/5 | 5 |
| 9 | 57 | 355 | 1.000 | 0.334 | 46 | *Bartonella*, *Chromobacterium* | 19% | 50/56 | 14 |
| 10 | 612 | 967 | — | 0.345 | 174 | *Xenorhabdus*, *Pseudomonas* | 91% | 153/178 | 0 |

**Architecture does not discriminate.** No candidate in any group lacks a shared
Pfam signature with the reference. The fraction carrying the reference's exact
architecture (PF22671 + PF04984 + PF17482) moves *non-monotonically* with
phylogenetic distance — 100% in group 0, 30% in group 1, 91% in the most distant
group 10 — which tracks whether the newer PF22671 model has been propagated to
an entry, not biology. Architecture cannot be used as an orthology criterion
here.

**Protein names do not discriminate.** Searching all 967 names for T6SS
(`vipA/vipB/tssB/tssC`), AFP, PVC and pyocin vocabulary returns **zero** matches
in every group. The entries are uniformly "phage tail sheath protein" or
"hypothetical protein". Name-based lineage assignment is not available.

**Locus context does not discriminate.** 406 candidates were put through the
unchanged `synteny.py` scoring (250 from the frozen sample plus 156 added by
`--extend-synteny` to cover the inner groups, which the frozen genus-spread
sample had left nearly empty). 355 resolved to a genomic locus and **all 355
scored ≥2 marker categories**, median 4:

| marker categories | 0 (unresolved) | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| candidates | 51 | 3 | 23 | 146 | 134 | 46 | 3 |

Per-group support spans 33% to 100% of *scored* candidates (group 6: 1/3; group
4: 5/9; groups 0, 2, 5, 8: 100%), but that spread is entirely a locus-resolution
artifact, not a biological difference: **every candidate whose locus resolved
scored ≥2 categories, in every group without exception** (363/363). The low
fractions are groups where IPG returned no location — 4 of group 4's 9 scored
candidates and 2 of group 6's 3 — and an unresolved locus is recorded as 0
categories. Among the groups with ≥20 scored candidates the range is 86–93%,
again limited only by resolution rate.

So the discriminating power is nil either way. Finding a baseplate, tube,
tape-measure and spike gene beside a sheath gene confirms it is a contractile
machine; it says nothing about *which* contractile machine, which is precisely
the distinction orthology needs.

**Taxonomy and tree support do discriminate**, and they agree with each other.
Groups 0–2 are a single genus at 0.65–0.83 identity with supports 0.889, 0.926
and 0.896. The first expansion beyond them (group 3) both crosses the weakest
node in the profile (**0.401**) and jumps taxonomic order into *Burkholderiales*.
On this evidence the natural boundary is level 2 — 124 sequences, 99%
*Pseudomonas*.

## 2. The alignment-width question, resolved

The 1310-column preliminary MSA over a 386-residue reference is **not** evidence
of family mixing in the ortholog set; it is an artifact of aligning the entire
967-sequence retrieval pool:

- 390 columns have ≥50% occupancy — essentially the reference's own 386 columns;
- 898 of 1310 columns are **below 5% occupancy**;
- the reference itself contributes **zero** residues to those sparse columns;
- residues in sparse columns concentrate in the distant, diverse groups
  (median 7 in group 7 and 14 in group 9, against 0 in groups 0–3).

Restricted to the candidate ortholog clade the alignment is 412 columns for a
386-residue reference — **1.07 columns per reference residue**. The width grows
with the boundary: 1.07× at level 2, 1.10× at level 4, 2.21× at level 7, 2.67×
at level 9. So alignment quality is itself one of the costs of widening.

## 3. Duplicated sheath genes — the decisive finding

Paralogy has to be read **before** dereplication. The frozen rule keeps one
representative per species by bitscore, so a genome carrying several prophage
sheaths contributes exactly one row and its duplication is invisible downstream.

- **221 of the 967 species carry more than one gated sheath candidate** (299
  extra sequences), all with distinct accessions and distinct sequences.
  *Escherichia coli* contributes 9, *Salmonella typhimurium* 6, *Chitinimonas
  arctica* 6.
- The bitscore margin between the kept representative and the runner-up has
  median **16 bits**, and **117 of 221 species have a margin below 20 bits** —
  many are 0 to 3. For over half the multi-paralog species, "best bitscore" is
  close to arbitrary.
- No tandem duplication: of 144 extended loci, 142 have exactly one sheath gene
  within ±15 kb and 2 have two. The copies are dispersed — consistent with
  independent prophage insertions rather than local gene duplication.

**The paralog test.** Putting the discarded paralogs back into the alignment and
tree for the level-2 clade (124 members + 57 same-species paralogs from 39
species):

- the clade is recovered intact (124/124 members);
- **39 of 57 discarded paralogs fall inside it**; 18 fall outside;
- they are **not monophyletic** — their MRCA spans 156 of the tree's 182 leaves
  (124 clade members + 57 paralogs + the reference);
- only **3 of 39** are the sister leaf of their own species' clade member;
- median patristic distance from an in-paralog to its clade member is 0.257
  against a clade-wide mean of 0.418.

That pattern is an **ancient duplication whose copies are interleaved through the
clade**, not a set of recent species-specific duplications that could be pruned.
Both copies of each pair are inside the clade, and the tree does not tell us
which one occupies the R-type pyocin position. Widening the boundary does not
help — at levels 7 and 9 **every** discarded paralog falls inside — though that
test loses power there, since the clade spans 30–37% of the whole retrieval pool.

![annotated clade tree]({{artifact:art_35b687a9-c2cf-47f9-a6a3-8371417e490a}})

## 4. Diagnostics at the 143 structural sites

| measurement | L2 (124) | L3 (146) | L4 (155) | L7 (293) | L9 (355) |
|---|---|---|---|---|---|
| distinct orthologs | 124 | 146 | 155 | 293 | 355 |
| identity min / median | 0.65 / 0.75 | 0.65 / 0.75 | 0.64 / 0.75 | 0.30 / 0.69 | 0.30 / 0.45 |
| effective genera | 1.05 | 1.95 | 2.57 | 21.3 | 34.5 |
| largest genus share | 99% | 84% | 79% | 42% | 35% |
| sites with occupancy ≥80% (frozen) | **143** | **143** | **143** | **143** | 142 |
| median site occupancy | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| median ±5-column window occupancy | 1.000 | 1.000 | 1.000 | 0.999 | 0.922 |
| sites with ±5-column window occupancy <80% | 1 | 2 | 2 | 27 | 36 |
| invariant sites | 46 | 43 | 43 | 6 | 1 |
| sites with one residue in ≥95% | 82 | 81 | 73 | 22 | 13 |
| median distinct residues per site | 2 | 2 | 2 | 9 | 11 |
| median column entropy (bits) | 0.185 | 0.246 | 0.292 | 1.503 | 1.794 |
| tree length | 7.8 | 9.6 | 11.2 | 69.5 | 106.7 |
| mean / max patristic distance | 0.42 / 1.01 | 0.47 / 0.96 | 0.50 / 1.05 | 2.08 / 5.02 | 2.56 / 6.71 |
| alignment columns per reference residue | 1.07 | 1.08 | 1.10 | 2.21 | 2.67 |
| outlier tip branches (>Q3 + 3·IQR) | 8 | 9 | 8 | 12 | 16 |

![feasibility trade-off]({{artifact:art_e7d0bee1-72bb-446c-bc31-b45e04cc1ac5}})

**The frozen 80% occupancy threshold is not the binding constraint and never
was.** It passes at 143/143 sites at every boundary except level 9. As
anticipated, occupancy says nothing about whether the aligned residues are
homologous: at level 7 all 143 sites clear 80% occupancy while 27 of them sit in
a ±5-column window that is below 80% occupied, i.e. in locally gappy regions
where column-to-column homology is least trustworthy.

## 5. Disclosure — an outcome-adjacent quantity was computed

While computing per-site substitution signal I also computed the Spearman
correlation between column entropy and the frozen `R_i` at the 143 sites:
**ρ = 0.32 at level 2 and ρ = 0.33 at level 7**. This is **not** a rate estimate,
is **not** covariate-adjusted, and is **not** the frozen primary analysis — but
it is adjacent to it, and it has now been seen. Recording it here so the
pre-registration status of the eventual rate-versus-remodeling test is accurately
documented rather than quietly assumed intact. Nothing in this report's verdict
depends on it, and no rates were estimated.

## 6. Is this specific to R2, or structural to the design?

Probably not specific to R2. The retrieval set is dominated by prophage-derived
sheaths, a gene family that is multi-copy and horizontally mobile — 23% of the
species sampled here carry more than one copy within the frozen identity window.
Any family whose reference sits in a prophage-like locus should be expected to
show the same hazard, so SalCIS, diffocin and the two AFP proteins should be
assumed to have it until shown otherwise. AFP may fare better: its paired-locus
structure gives a positional constraint that a single-gene family does not, and
Afp2/Afp3 must be evaluated together for exactly that reason.

What would be fatal to the design, rather than to R2, is if positional orthology
cannot be established for any of the five families. In that case site-rate
estimation per family would rest on sets of mixed paralogous lineages, and the
`R_i`-versus-rate test could not be interpreted as a within-ortholog comparison.

## 7. Recommended next step (not executed)

1. Attempt **positional orthology from locus gene order** on R2: extract the
   ordered CDS products around each candidate and each discarded paralog, and
   test whether the R-type pyocin gene order distinguishes one copy per genome.
   This is a new discriminator, not a change to any frozen rule.
2. If it resolves the copies, rebuild the ortholog set at a wide boundary
   restricted to positionally orthologous copies, and re-run these diagnostics.
   That set would have both the substitution signal of level 7 and a defensible
   orthology claim.
3. If it does not resolve them, report R2 as failing the evolutionary
   feasibility gate and re-examine whether the study can support per-family
   site-rate estimation at all.

Scaling to Afp2, Afp3, SalCIS and diffocin should wait for the outcome of step 1
or 3. Rate4Site remains unrun, and the inference specification (normalization,
spatial dependence, uncertainty propagation, whether four systems support any
cross-system claim) remains open and un-prejudged.

## Files

| file | content |
|---|---|
| `orthology/r2pyocin_groups.csv` | per-group architecture, names, taxonomy, synteny, insertion load |
| `orthology/r2pyocin_boundary_sweep.csv` | the level-by-level decision table above |
| `orthology/r2pyocin_L*_feasibility.json` | full diagnostics per boundary |
| `orthology/r2pyocin_L*_site_diagnostics.csv` | per-site occupancy, entropy, window occupancy |
| `orthology/r2pyocin_L2_paralog_test.csv` + `_inparalog_pairs.csv` | paralog test result and pair topology |
| `orthology/r2pyocin_L*_ortholog_set.{fasta,afa,nwk}` | candidate sets, alignments, trees |
| `orthology/r2pyocin_duplicates.csv` | the 221 multi-paralog species with bitscore margins |
| `synteny/r2pyocin_synteny_extra.csv` | the 156 added loci, same frozen scoring |
| `orthology_r2pyocin.py`, `ortholog_set_r2pyocin.py` | the two analyses, rerunnable |
