# R2 pyocin — copy-number audit and final feasibility decision

**Date** 2026-10-13 (claims 1 and 4 corrected after review) · **Reference**
Q9S574 / NP_249313.1 · **Scope** bounded audit from data already retrieved, plus
one bounded tree test

## Decision

**R2 is excluded from the manuscript's primary evolutionary test.**

> The coexistence of divergent sheath homologs within individual genomes,
> together with their interleaved phylogenetic placement and conserved locus
> architecture, prevents reliable assignment of a single orthologous lineage
> using the current sequence- and context-based criteria.

That is the whole of the claim. It does not require reconstructing an
unobserved duplication history, and this document does not attempt one.

**This is an exclusion for indefinability, not for inconvenience.** The dataset
is not being rejected because it lacks variation — the restricted set has ample
variation (section 3) — and not because it gave an unwelcome answer; no rate was
ever estimated. It is rejected because the evolutionary comparison the analysis
would make cannot be *defined* with confidence. That distinction belongs in the
manuscript's methods and supplement, stated in those terms.

| component | decision |
|---|---|
| structural remodeling measurement | **PASS** |
| sequence retrieval | **PASS** |
| alignment feasibility | **PASS** |
| original orthology assignment (species-dereplicated) | **FAIL** |
| restricted single-copy assignment | **UNRESOLVED** |
| primary evolutionary-rate analysis | **EXCLUDE** |

R2's structural contribution stands unchanged: its `R_i` values are measured
from deposited coordinates, have no evolutionary component, and remain in
Result 1. What it does not contribute is Result 3.

Proceed to **AFP (Afp2 + Afp3 jointly), SalCIS and diffocin**, applying the
audit below as a gate *before* investing in their orthology — it is cheap and
runs off the retrieval output.

## 1. The audit, and exactly what it counts

No new retrieval was needed. The 8,035-record UniProtKB Pfam-family scan behind
the R2 candidate set is upstream of every gate: the identity window, coverage
and length filters are applied after it, so family copies the gates excluded are
still in the record. Counting records per UniProt proteome therefore gives copy
number with no identity filter.

**What it is not.** This is a count of **Pfam-detected copies per UniProt
proteome**, not a definitive genome-wide copy number. It is conditional on the
three Pfam signatures being matched and on the completeness and annotation state
of each proteome. The inference it supports is therefore **asymmetric**:

- **two or more detected copies is strong evidence** that the genome is
  multi-copy for this family;
- **one detected copy does not establish** that no further, more divergent
  sheath homolog exists in that genome — only that none was detected by these
  signatures in that proteome.

Every "single-copy" statement below should be read as "no further Pfam-detected
copy". This asymmetry strengthens rather than weakens the case against treating
the 136 sequences as verified single-copy orthologs: the 136 is an upper bound
on how many of them are genuinely single-copy.

**Control.** *P. aeruginosa* PAO1 (UP000002438) returns exactly **one** record,
the R-type pyocin sheath itself (G3XD39, 386 aa, full PF22671+PF04984+PF17482
architecture). The census is specific to the phage-tail-sheath family and does
not pick up PAO1's three T6SS sheaths, which use different Pfam models — the
right scope for an orthology question about this family.

**Census.** Across the **5,133 proteomes** represented in the family scan:

| Pfam-detected copies | 1 | 2 | 3 | 4 | 5 | 6+ |
|---|---|---|---|---|---|---|
| proteomes | 3,341 | 1,178 | 432 | 99 | 38 | 45 |

65.1% with a single detected copy, median 1, mean 1.56, maximum 20. By the
asymmetry above, the 35% multi-copy figure is a **lower bound**.

## 2. The single-copy restriction does not hold up

Option A from the previous decision was to keep only members from species
contributing exactly one *gated* candidate — 211 of the 293 members at the wide
boundary. Audited against the unfiltered scan, **only 136 of those 211 (64.5%)
have no further detected copy in their proteome**:

| further Pfam-detected copies in the same proteome | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| sequences | **136** | 52 | 18 | 2 | 3 |

Of the 75 that are not single-copy, **16 have an additional copy inside the
30–90% identity window** that was excluded on coverage or length — a plausible
ortholog candidate removed by a tractability gate rather than by biology. None
has a copy above 90% identity to the reference, so there are no recent
near-identical duplicates inflating the count.

![copy-number audit]({{artifact:art_6afe47b8-676b-4a8f-bb96-b8506e637e2f}})

## 3. The restricted set passes every diagnostic — and still is not an ortholog set

Restricting to the 136 sequences with no further detected copy gives, on the
face of it, the best dataset this project has produced for R2:

| measurement | value |
|---|---|
| sequences | 136 |
| genera / effective genera | 57 / 22.3 |
| largest genus share | 37% (*Pseudomonas*) |
| identity to reference, min / median | 0.30 / 0.46 |
| structural sites with occupancy ≥80% | **143 of 143** |
| invariant structural sites | 8 |
| sites with one residue in ≥95% | 17 |
| median distinct residues per site | 7 |
| median column entropy | 1.54 bits |
| sites in a locally gappy window | 14 of 143 |
| alignment columns per reference residue | 1.90 |
| tree length / mean patristic | 34.9 / 1.83 |

Every quantity the feasibility report identified as binding is satisfied. If the
question were "does this set have enough substitutions and taxonomic
independence", the answer would be yes.

**It is not a defensible ortholog set, and the reason is not visible in any of
those numbers.** One-detected-copy-per-proteome addresses *within-genome*
paralogy among detected copies. It does not establish that the retained copies
across genomes correspond to one another — that different genomes retained the
same lineage rather than different ones. Nothing in the diagnostics above tests
correspondence; they test variation, occupancy and taxonomic spread.

## 4. The divergence test, and what it does and does not show

Built one tree from the 136 restricted sequences plus both copies from 25
multi-copy genomes inside the same boundary, and asked how far apart the two
copies of a single genome are relative to the spread of the backbone. Patristic
distance is used rather than clade membership because it does not depend on
rooting.

| | median | 90th percentile | max |
|---|---|---|---|
| between genomes, single-copy pairs (n = 9,180) | 2.168 | 2.808 | 3.943 |
| within one genome, copy pairs (n = 25) | **2.674** | — | 3.562 |

Two copies in one genome are a median **1.23×** as divergent as two single-copy
sequences from different genera. **13 of 25** within-genome pairs exceed the
between-genome median; **10 of 25** exceed its 90th percentile. The most
divergent pair, at 3.562, is close to the maximum distance found anywhere in the
backbone.

**What this shows.** Homologs co-occurring in one genome can be as divergent as,
or more divergent than, sequences drawn from different genera. Co-occurrence
therefore carries no implication of close relatedness in this family, and a
restriction to one copy per genome does not make the retained copies
correspondent across genomes.

**What this does not show.** It does **not** establish that the duplications
predate the diversification of the sampled taxa. Horizontal transfer,
lineage-specific duplication, differing evolutionary rates among copies, and
other gene-history processes can all produce large within-genome divergence, and
this analysis does not distinguish them. An earlier draft of this document drew
that inference; it has been withdrawn. Reconstructing the duplication history is
neither attempted here nor needed for the decision.

Two caveats on the comparison itself. The 9,180 between-genome pairs are **not
independent observations** — they are all pairs drawn from 136 sequences on one
tree, so the between-genome distribution should be read as a descriptive
reference spread, not as a sample supporting a significance test; no p-value is
computed from it and none should be. And the 25 within-genome distances are
bimodal, from 0.018 to 3.562, so "divergent co-occurring homologs" is a
statement about part of the sample, not all of it.

What remains is the conjunction stated in the decision: divergent homologs
co-occur within genomes, they are interleaved in the gene tree rather than
separated, and their locus architecture is conserved to the point that gene
order cannot tell them apart. Those three observations together leave no
criterion among the ones in use that assigns a single orthologous lineage.

## 5. What would be needed, and why it is out of scope

Assigning a single orthologous lineage here would require an explicit model of
the family's gene history — duplication, loss and transfer — fitted and then
carried into the rate estimates' uncertainty. That is a methodological
contribution in its own right, and it is not what the manuscript set out to
establish: its question is whether contact remodeling explains site-specific
evolutionary rate beyond established structural determinants, and R2 is one of
five proteins available to answer it.

Spending the project's remaining effort on that machinery for one family would
invert the priority. R2 is excluded because the comparison cannot be defined
with the criteria in use, not because its data are poor and not because a
rate was estimated and disliked.

## 6. Carry forward to AFP, SalCIS and diffocin

Run the audit in this order, before any orthology work:

1. **Pfam-detected copy census from the retrieval output** — free, no fetches.
   Multi-copy detection is informative; a single detection is not a guarantee.

   **Measured limit, found on AFP 2026-10-13 while applying this step.** The
   census saw **none of the three sheath genes in the characterised AFP
   reference locus**. Afp2 (Q6HAD7) and Afp3 (Q6HAD6) are in the retrieval pool
   but carry no `Proteomes` value at all — they are plasmid-borne pADAP entries
   — so they contribute nothing to a per-proteome count; the third in-locus
   sheath gene (WP_169558217.1) has no entry in the pool whatsoever; and the
   *Serratia entomophila* reference proteome UP001056873 reports a single copy
   that is none of the three. A per-proteome census is therefore a
   **positive-only screen with no negative predictive value**, and it can be
   blind to an entire plasmid-borne locus. Use it to raise suspicion, never to
   clear a family.
2. **Within-genome divergence against the between-genome spread** — one
   alignment and one tree. If within-genome copy pairs reach past the
   between-genome range, co-occurrence carries no relatedness information in
   that family and a one-copy-per-genome restriction will not establish
   correspondence.
3. Only then the clade-boundary and locus work.

AFP is the most promising of the three and should be first: Afp2 and Afp3 are
distinct proteins of the *same* locus, so a candidate pair occupying
corresponding positions in a homologous AFP locus is a positional constraint no
single-gene family offers — and the frozen protocol already treats paired-locus
recovery as an output. It is also the family where the audit matters most,
because Afp2 and Afp3 are themselves related, and a mis-assigned copy could be
the *other* protein of the same machine.

Two defects found while characterising the R2 reference locus apply to all three
and should be fixed before their synteny results are interpreted: the frozen
marker vocabulary scores the true R2 reference locus at **zero** categories, and
the frozen ±15 kb window excludes the pyocin-specific landmarks at slots −10 to
−13. Both are documented in `ORTHOLOGY_DECISION_r2pyocin.md`.

## 7. Status of R2's outputs

| output | status |
|---|---|
| `sheath_pilot/qc_output/R_r2pyocin.csv` — 143 structural sites, frozen `R_i` | **stands**, Result 1 |
| `filtered/r2pyocin_filtered.csv` — 967 retrieved candidates | stands as a retrieval result |
| the 39 interleaved paralogs and the 221 multi-copy species | **retained as evidence**, not pruned |
| R2 as a single-copy ortholog set for rate estimation | **withdrawn** |
| Rate4Site on R2 | not run, and not to be run |

## Files

| file | content |
|---|---|
| `orthology/r2pyocin_genome_copies.csv` | copy count for all 5,133 censused genomes |
| `orthology/r2pyocin_copy_audit_L7single.csv` | per-member audit of the 211 |
| `orthology/r2pyocin_copy_audit_all.csv` | the same for all 967 |
| `orthology/r2pyocin_singlecopy_feasibility.json` | diagnostics for the 136-sequence set |
| `orthology/r2pyocin_singlecopy_site_diagnostics.csv` | per-site occupancy, entropy for the 136 |
| `orthology/r2pyocin_within_genome_distance.csv` | the 25 within-genome copy distances |
| `orthology/r2pyocin_straddle_test.csv` | the rooting-dependent clade-span version |
| `orthology/r2pyocin_singlecopy_genome.{afa,nwk}` | alignment and tree of the 136 |
| `orthology/r2pyocin_straddle.{afa,nwk}` | the 187-taxon test tree |
