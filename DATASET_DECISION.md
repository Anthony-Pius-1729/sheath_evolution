# Definitive list of evolutionary datasets for the hypothesis test

**Date:** 2026-10-13. Closes the family-by-family feasibility assessment.
**Hypothesis under test:** does residue-level contact remodeling during contraction predict
evolutionary rates independently of solvent accessibility, packing and flexibility?

---

## 1. The list

| system | dataset | n seq | effective genera | structural sites | informative sites | median entropy | status |
|---|---|---|---|---|---|---|---|
| Diffocin (8V3X) | frozen, offset 0, ambiguous genome removed | 33 | 16.71 | 58 | **50 (86%)** | 1.444 | **primary** |
| Diffocin (8V3X) | frozen, offset 0, as screened | 35 | 17.71 | 58 | 50 (86%) | 1.444 | sensitivity |
| AFP (6RBN e3) | Afp3 position | 31 | 5.43 | 82 | **65 (79%)** | 0.823 | **primary** |
| AFP (6RBN e2) | Afp2 position | 27 | 5.49 | 84 | **51 (61%)** | 0.554 | **primary, reduced power** |
| Diffocin (8V3X) | collinear frame | 66 | 31.27 | 58 | 53 (91%) | 1.635 | sensitivity |
| R2 pyocin (3J9Q) | — | — | — | 143 | — | — | **excluded** |
| SalCIS (9R9H) | — | — | — | 68 | — | — | **excluded** |

**Three datasets, from two systems, 223 of 435 primary structural sites fitted**
(58 + 84 + 81; one Afp3 site falls below 80% column occupancy).

Effective-genera values for the diffocin sets are post-deduplication: WP_121586442.1 reached
the screened set through two UniProt accessions for one *Anaerotruncus massiliensis* replicon
and is counted once (`DIFFOCIN_FEASIBILITY.md` section 7b).

![dataset decision]({{artifact:art_f610c25b-d0d1-4b1b-b724-01ccbf2e6405}})

## 2. Why each system landed where it did

**Diffocin — PASS (`DIFFOCIN_FEASIBILITY.md`).** Single sheath gene in the reference locus
and in **each of the 142 sampled loci** across 80 genera — one detected sheath gene per
locus, though some genomes carried more than one candidate locus (7 of 137 replicons); 67
loci anchored on non-sheath flanking
homology at median collinearity rho = 1.00; 34 orthologs at the reference slot, position and
sequence agreeing 36/36. All 58 structural sites at >=80% occupancy, 8 invariant, deepest
pair 2.44 substitutions/site. Declared: one ambiguous species (*Peptacetobacter hiranonis*
contributes two mutually exclusive assignments), 5% of replicons carrying a divergent second
copy, orthology resting on locus position rather than monophyly.

**AFP — PASS on restricted per-subfamily sets (`AFP_FEASIBILITY.md`, `AFP_CHECK23_VERDICT.md`).**
Afp2 and Afp3 are 38.1% identical to each other, so sequence similarity cannot separate the
two tandem genes — best-hit labelling disagrees with position for 15 of 22 genes at the Afp2
slot, and 9 of 29 nodes mixing the labels carry support >=0.9. Positional assignment inside
the reference subfamily works (95% agreement, 35 of 42 loci carrying the full three-gene
run). Declared: a single enterobacterial radiation (5.4 effective genera), 39% invariant
structural sites in Afp2, deep-pair saturation in Afp3, and **Afp2 and Afp3 are not
independent** — same loci, same genomes, one particle.

**R2 pyocin — EXCLUDED (`R2_FINAL_DECISION.md`).** 967 species-dereplicated homologs, but
221 species carry >=2 gated candidates, 39 of 57 discarded paralogs fall inside the candidate
clade, the frozen marker vocabulary scores the true reference locus at zero categories, and
in 11 of 35 testable pairs the discarded paralog reconstructs the reference gene order
*better* than the retained copy. The best single-copy set (136 sequences) passes every
summary diagnostic yet has within-genome copy pairs more divergent than between-genome pairs
(2.674 vs 2.168). Operative claim: the coexistence of divergent sheath homologs within
individual genomes, together with their interleaved phylogenetic placement and conserved
locus architecture, prevents reliable assignment of a single orthologous lineage using the
current sequence- and context-based criteria.

**SalCIS — EXCLUDED (`SALCIS_DECISION.md`).** Of 4,490 blastp hits, exactly one is
*Salmonella* — the reference itself. The gated set is dominated by Bacteroidota
(*Chryseobacterium*, *Flavobacterium*, *Algoriphagus*...) at 30-45% identity to a
Gammaproteobacterial reference. The reference locus is unrecoverable in 84% of sampled loci,
and 28 of 55 loci carry two tandem sheath genes against the reference's one. Maximum
achievable positional ortholog count: 9.

## 3. What this does to the manuscript

**Result 2 justifies the dataset; it is not a second research question.** **Four systems
were evaluated for evolutionary inclusion — R2 pyocin, AFP, SalCIS and diffocin — and two
failed.** (The five structural references comprise four systems: AFP contributes two
proteins, Afp2 and Afp3, not two systems.) The reason for both failures is the same and is
worth reporting: contractile sheath genes occur in multi-copy, mobile, recombining elements,
and the copies within a genome are frequently as divergent as copies between genomes. These
are methodological findings whose role is to justify which datasets enter the test — not an
independent line of inquiry.

**The hypothesis test becomes within-protein, across-site, with three proteins from two
systems.**
`R_i` is a per-residue property, the rate is a per-site estimate, and the association is
tested across sites within each protein — n = 58, 82, 84 sites. The three datasets give
three independent replications of that test (two of which share a locus history), not a
five-system comparative analysis.

**What is no longer claimable.** Any cross-system statement about how contact remodeling
constrains sheath evolution *in general*. Two systems, one of them contributing two
non-independent proteins, cannot support a phylogenetically distributed generalisation. The
manuscript must test the hypothesis per protein, report concordance or discordance across
the three, and treat consistency across them as suggestive rather than as a tested
cross-system claim.

**Result 1 is unaffected.** The frozen `R_i` metric stands for all five families (435
in-domain sites), including R2 pyocin and SalCIS. Reproducible quantification of contact
remodeling does not require an ortholog set.

## 4. What happens next, in order

1. **Freeze the statistical analysis** before any Rate4Site run: primary model, covariate
   definitions (relative solvent accessibility, packing/contact number, B-factor or
   equivalent flexibility proxy), handling of spatial dependence among neighbouring
   residues, propagation of per-site rate uncertainty, and the exact sensitivity analyses
   (diffocin collinear-frame set; diffocin without *P. hiranonis*; Afp2 with and without its
   33 invariant sites).
2. **Run Rate4Site** on the three primary alignments — and only then.
3. **Test, and report the direction obtained.**

## 5. Disclosure carried forward

Spearman correlation between alignment-column entropy and frozen `R_i` at the 143 R2
structural sites was computed as an alignment diagnostic before this decision: **rho = 0.32
(narrow set), 0.33 (wide set)**. It is not covariate-adjusted and not a rate estimate, but
it is outcome-adjacent and it has been seen. The hypothesis cannot be described as
completely unexamined, and the manuscript carries the three-tier disclosure: (1) exploratory
analyses already performed, (2) decisions fixed before any Rate4Site run, (3) genuinely
held-out tests. R2 pyocin, where that correlation was computed, is excluded from the test.

## 6. Files

`orthology/dataset_decision.csv` · `orthology/dataset_decision.png` ·
`AFP_FEASIBILITY.md` · `DIFFOCIN_FEASIBILITY.md` · `SALCIS_DECISION.md` ·
`R2_FINAL_DECISION.md` · `AFP_CHECK23_VERDICT.md` ·
analysis-ready alignments and trees: `orthology/diffocin_frozen_offset0.{afa,nwk}`,
`orthology/afp_sf1_{afp2_like,afp3_like}.{afa,nwk}`, plus
`orthology/diffocin_collinear_frame.{afa,nwk}` for sensitivity.
