# Diffocin orthology and feasibility decision

**Reference:** PDB 8V3X entity 1, *Clostridioides difficile* sheath CD1363, 354 aa,
UniProt A0A9Q7ZU73 / RefSeq WP_003438414.1. **Date:** 2026-10-13.
**Verdict: PASS — the strongest evolutionary dataset in the project, with one declared
ambiguity.**

---

## 1. Reference identity and locus

UniProt sequence **matches reference_sequences.json exactly** (354 aa). Pfam signatures
PF04984 + PF17482 — the same pair as Afp2/Afp3, so diffocin and AFP start from an identical
8,005-record pool and are separated only by sequence and locus evidence.

Reference locus (NZ_AP031492.1, assembly GCF_040436335.1): **51 CDS in ±25 kb, exactly one
sheath-family gene** (slot 0, 703 bits). The neighbourhood is a canonical prophage tail
morphogenesis cassette, fully annotated: phage tail terminator (−1), **tail tube (+1)**,
tail assembly chaperone (+2), phage tail family protein (+4), baseplate protein (+5),
XkdQ/YqbQ (+6), DUF2577 (+7), DUF2634 (+8); upstream a prophage regulatory module (XRE-family
and helix-turn-helix regulators, ImmA/IrrE metallo-endopeptidase) at slots −3 to −5.

Single reference sheath gene means **the AFP labelling problem cannot arise**: there is no
second tandem copy to confuse an assignment with.

## 2. What the retrieval returned

8,005 records → **930 blastp hits** → 192 through the frozen gates → **165
species-dereplicated across 92 genera**; identity median 0.406, range 0.301–**0.646**,
coverage median 0.994.

Diffocin is the most isolated of the five references: 930 hits from the same pool that gave
SalCIS 4,490, and **no homolog closer than 65%**. Of 8 *Clostridioides* hits, 4 are ≥90%
identical and excluded by the ceiling; one *Clostridioides* member survives gating. The
retained set is therefore cross-genus by construction, which is what makes it informative.

## 3. Locus architecture across the family

142 loci characterised from the entire screened pool (151 candidates at ≥200 bits, 141
resolved + cache):

- **Sheath-run length = 1 in every one of the 142 loci.** Each sampled locus contains
  exactly one detected sheath gene, across 80 genera — no family in this project has shown
  this. **This is locus-level, not genome-level:** some genomes contain more than one
  candidate locus (§5, 7 of 137 replicons), each with its single sheath gene. Single-copy
  architecture per locus must not be read as single-copy orthology per genome.
- **67 of 142 anchored** on non-sheath flanking homology (median 9 anchors, median modal
  anchor share 0.727, **median collinearity ρ = 1.00**).
- Modal offsets among anchored loci: **−1 (22), 0 (36), +1 (8), +2 (1)**.
- Of 109 unassigned sheath genes, **78 are unassigned because their locus could not be
  anchored** — the position test could not be run — and only 31 because the gene fell
  outside the mapped slot. Unassigned does not mean positionally rejected.
- Among all assignments, **position and sequence best hit agree 36/36**.

### The offset rule, stated plainly

The frozen rule labels a sheath gene only when it sits at the reference slot *exactly*
(modal offset 0), giving **34 independent orthologs + reference = 35**. Twenty-two further
loci are anchored with median collinearity ρ = 1.00 and differ from the reference frame by
**one flanking gene** (offset −1). Those are the same locus with a single insertion or
deletion upstream, not a different locus.

I am **not** relaxing the rule to gain sequences. Both sets are reported, the frozen set is
primary, and the larger collinear-frame set (**n = 66**) is a declared sensitivity analysis.
The difference between them is one indel in the flank, and that is the whole of it.

## 4. Feasibility at the 58 structural sites

| | **frozen (offset 0)** | collinear frame (sensitivity) |
|---|---|---|
| sequences (incl. reference) | **35** | 66 |
| species / genera | 30 / 21 | 60 / 39 |
| **effective genera** | **17.71** | 31.27 |
| largest genus share | 0.118 | 0.092 |
| alignment columns | 535 | 565 |
| columns per reference residue | 1.51 | 1.60 |
| structural sites scored | 58 | 58 |
| **sites at ≥80% occupancy** | **58 / 58** | **58 / 58** |
| **invariant structural sites** | **8 (14%)** | 5 (9%) |
| sites with ≥95% majority residue | 9 | 10 |
| median per-site entropy (bits) | **1.444** | 1.635 |
| maximum possible entropy at this n | 5.13 | 6.04 |
| tree length | 16.18 | 27.83 |
| mean / max patristic distance | 1.592 / 2.444 | 1.710 / 3.114 |

**Informative sites: 50 of 58 (86%).** Every structural site is scored at full occupancy;
taxonomic breadth is 17.71 effective genera with no genus above 11.4%; the deepest pair is
2.44 substitutions/site, short of saturation. On every axis this is the best set in the
project — compare AFP's 5.4 effective genera and 21–39% invariant sites, and R2's narrow
set at 1.05 effective genera and 32% invariant.

## 5. The ambiguity that must be declared

Seven of 137 replicons (5%) carry a second recovered sheath gene. Where they do, the two
copies are as divergent as copies from different genomes: within-replicon patristic median
**2.268** against between-replicon **1.942** (1.17×, n = 7 pairs; 4 of 7 above the
between-genome median, 3 of 7 above its 90th percentile, range 0.006–4.070). This is the
same signature that defeated R2 pyocin — there at 23% of species, here at 5% of replicons.

Three of the seven pairs touch the retained sets, and one is a genuine internal
contradiction: ***Peptacetobacter hiranonis* (NZ_CP036523.1) contributes two frozen-set
members**, WP_006439294.1 and WP_006440377.1 — two separate collinear loci each with one
sheath gene at the reference position. A single genome cannot hold two positional orthologs
of one reference gene; this genome has a duplicated *locus*, not a duplicated gene.

Handling, fixed now: both sequences are **retained** in the dataset and in
`diffocin_sets.json`, and *P. hiranonis* is declared an **ambiguous species**. The primary
analysis will be reported with it and repeated without it as a named sensitivity test
(4 of 34 members sit on replicons with a second copy; 2 of those 4 are this pair). Nothing
is pruned to tidy the tree.

## 6. Coherence, and what it does and does not show

In a joint tree of all 142 recovered genes plus the reference, the MRCA of the 35 frozen-set
members spans 142 of 145 tips (purity 0.254): the set is **not** a monophyletic clade. That
is expected here and is not the AFP or R2 situation. The interleaved sequences are
overwhelmingly loci that **could not be anchored** (78 of 109), i.e. sequences of unknown
position, not sequences shown to occupy a different position; and there are no tandem
paralogs anywhere in the family to interleave. The positional claim is per-locus and does
not require the set to be monophyletic.

What this does mean: diffocin's orthology rests on **locus position and architecture**, not
on tree topology, and the manuscript must say so.

## 7. Verdict

**PASS.** The frozen 35-sequence set (34 orthologs + reference) is suitable for
site-specific rate estimation at all 58 structural sites: positional orthology established
per locus, one sheath gene per locus throughout the family, 17.71 effective genera,
full occupancy, 86% of structural sites informative, no saturation.

Declared limits: one ambiguous species (*P. hiranonis*, two assignments), 5% of replicons
carrying a divergent second copy, orthology resting on locus evidence rather than
monophyly, and n = 35 (the frozen target was 30–50).

## 7b. Correction applied 2026-10-13 (after the verdict, before any rate was estimated)

Two records in the screened set resolve to **one** protein: WP_121586442.1 reached through
two UniProt accessions (A0ABR7ACX6, A0A498D022) for the same *Anaerotruncus massiliensis*
replicon NZ_JACOIH010000005.1. It is one gene, not two, and is now counted once; the
alignments were already correct (identical names collapse), so only the composition counts
change. Effective genera 17.71 (not 17.85) for the frozen set and 31.27 (not 31.19) for the
collinear-frame set.

**Primary dataset, as used in the hypothesis test.** On instruction, *both* ambiguous
*P. hiranonis* copies are **removed** from the primary set rather than retained, because
retaining two mutually exclusive positional assignments contradicts the orthology principle
used to exclude R2 pyocin. The primary diffocin dataset is therefore **32 orthologs +
reference = 33 sequences**, 29 species, 20 genera, 16.71 effective genera, largest genus
12.5%. The 35-sequence set as screened is retained as a prespecified sensitivity analysis,
as is the 66-sequence collinear-frame set. The exclusion is justified by unresolved
positional orthology and was decided before any evolutionary rate was computed.

## 8. Files

`orthology/diffocin_reference_track.json` ·
`orthology/diffocin_all_{loci,sheath_runs,locus_tracks}` ·
`orthology/diffocin_feasibility.csv` · `orthology/diffocin_sets.json` ·
`orthology/diffocin_{frozen_offset0,collinear_frame,joint}.{afa,nwk}` ·
`orthology/diffocin_within_genome_distance.csv` · `filtered/diffocin_filtered.csv`

Reproduce: `python afp_loci.py --n 160 --dedup species --out-prefix diffocin_all
--ref-track diffocin_reference_track.json --cache diffocin_locus_cache.json
--sheath-slots "0:diffocin_like" --candidates diffocin_candidates.tsv
--prescreen diffocin_blastp_hits.csv`
