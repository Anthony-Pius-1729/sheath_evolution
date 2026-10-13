# R2 pyocin — positional orthology investigation and orthology decision

**Date** 2026-10-13 · **Reference** Q9S574 / NP_249313.1, *Pseudomonas
aeruginosa* PAO1, NC_002516.2 ~681 kb · **Scope** targeted locus comparison, no
new retrieval, no new tree search

## Decision

**Orthology: UNRESOLVED at single-copy resolution, and not forced.** The
duplication is of the *entire pyocin-like cluster*, regulatory cassette
included, so gene order cannot say which copy occupies the position the
reference occupies. In 25 of the 35 testable pairs both copies recover ≥15 of
the reference locus's slots, and in **11 of 35 the discarded paralog matches the
reference locus better than the copy dereplication retained**.

**R2 should not enter the manuscript's primary evolutionary analysis as a
single-copy ortholog set.** It can enter on either of two bases, and that choice
is yours (section 6): a single-copy-species restriction, which is defensible and
leaves 211 sequences across 75 genera at the wide boundary, or an explicitly
paralogy-aware treatment. What it cannot do is enter as "124 orthologs" or "293
orthologs" chosen by best bitscore, because that selection is contradicted by
locus evidence in a third of the cases where the question arises.

## 1. The reference locus, and a finding about the frozen synteny stage

The R2 reference sits in a 57-CDS neighbourhood (±25 kb) with this structure:

| slot | gene | note |
|---|---|---|
| −13 | anthranilate synthase component I (*trpE*) | the documented left boundary of the *P. aeruginosa* pyocin region |
| −12 | transcriptional regulator **PrtN** (NP_249301.1) | pyocin-specific activator |
| −11 | HTH transcriptional regulator **PrtR** (NP_249302.1) | pyocin-specific repressor |
| −10 | repressor **PtrB** (NP_249303.1) | pyocin regulatory cassette |
| −9 … −1 | "hypothetical protein", "bacteriophage protein" | the tail-gene block |
| **0** | **"bacteriophage protein"** (NP_249313.1) | **the sheath itself** |
| +1 … +22 | "bacteriophage protein", "hypothetical protein" | rest of the cluster |

Two consequences, and the first is a methodological problem with the frozen
stage rather than with R2:

**The frozen marker vocabulary scores the true reference locus at zero.** Not one
of the nine frozen marker categories (`tube`, `baseplate`, `tape_measure`,
`spike`, `terminator`, `fibre`, `atpase`, `phage_tail`, `sheath`) matches any
product name in this locus — the sheath gene itself is annotated "bacteriophage
protein". Scored by `synteny.py` the R2 reference locus would read **0
categories, "sparse/uninformative annotation"**, exactly the label the frozen
code reserves for uninformative loci. The earlier result that 363 of 363 resolved
candidate loci score ≥2 categories therefore measures *annotation vocabulary*,
not locus identity: the candidates "pass" a test their own reference fails.

**The frozen ±15 kb window excludes the only pyocin-specific landmarks.** PrtN,
PrtR and *trpE* sit at slots −11 to −13, outside ±15 kb of the sheath. The
discriminating information was never in the window.

## 2. What was compared, defined before any result was seen

Product text cannot be the discriminator — it would score the reference against
itself as a non-match. The comparison is therefore gene-order conservation by
**sequence homology of the neighbours**, the standard positional-orthology
criterion and immune to annotation vocabulary:

1. ordered CDS neighbourhood within ±25 kb, oriented so the sheath anchor reads
   on the + strand, slots signed from the anchor;
2. every neighbour protein blastp'd against the 57 reference-locus proteins;
3. each query gene assigned the reference slot of its best hit (E ≤ 1e-3);
4. per locus: slots recovered, genes at the *same* slot as their counterpart,
   Spearman collinearity of query slot against reference slot, and presence of
   PrtN/PrtR/PtrB homologs — the cassette named *from* the reference, not
   hardcoded.

No threshold was defined and nothing was excluded by this analysis.

## 3. Result: both copies are in pyocin-like clusters

39 in-paralog pairs; 35 testable (4 had no resolvable locus for one copy:
*P. oryzicola*, *P. syringae*, *P. frederiksbergensis*, *P. versuta*).

| | retained copy | discarded paralog |
|---|---|---|
| reference slots recovered (median of 57) | 19 | 17 |
| genes at the identical slot (median) | 12 | 10 |
| regulators present (median of 3) | 2 | 1 |
| all three regulators present | 8 loci | 8 loci |
| collinearity ρ ≥ 0.9 | 22 of 35 | 14 of 35 |
| sheath anchor homologous to the reference sheath | 29/29 | 35/35 |

Pair-by-pair, by slots recovered: the retained copy is better in **18**, the
paralog is better in **11**, and 6 are tied. Separation is clear (|Δ| ≥ 5 slots)
in only **11 of 35** pairs, and those cut both ways.

![locus comparison]({{artifact:art_0ed07897-3aa0-4349-a421-5f8cd74facc3}})

The cases that settle the question are the ones where locus evidence contradicts
the sequence-similarity choice:

| species | retained copy | discarded paralog | reading |
|---|---|---|---|
| *P. soli* | 8 slots, 0 regulators, ρ=1.00 | **17 slots, 1 regulator** | dereplication kept the weaker positional match |
| *P. caspiana* | 9 slots, 0 regulators | **15 slots, 2 regulators** | same |
| *P. chlororaphis* | 20 slots, 3 regulators | 21 slots, 3 regulators | two complete, indistinguishable clusters |
| *P. lurida* | 19 slots, 2 regulators | 20 slots, 3 regulators (×2 copies) | ≥3 pyocin-like clusters in one genome |
| *P. guariconensis* | 18 slots, 1 regulator | 4 slots, 0 regulators | gene order *does* support the retained copy |

*P. chlororaphis* is the general case: both copies carry the complete PrtN–PrtR–PtrB
cassette and recover the same 13 exact slots. There is nothing in the locus to
choose between them. *P. lurida* is worse — the genome carries at least three
such clusters, two of which score higher than the retained copy.

## 4. Why this is a definitional gap, not a data problem

The frozen orthology rule is tree placement **combined with** locus context, as a
four-way classification. Both components are satisfied by **both** copies of each
pair: both fall inside the reference clade, and both sit in a pyocin-like locus.
The four-way classification has no provision for two copies from the same genome
that both satisfy it — the only thing separating them is the dereplication rule,
"one representative per species, best bitscore", which:

- decides by a margin under 20 bits in 117 of 221 multi-copy species;
- is contradicted by locus evidence in 11 of the 35 pairs tested here;
- was never intended as an orthology criterion. It is described in
  `filter_hits.py` as dereplication of "heavy strain-level redundancy in
  RefSeq" — a tractability measure for the synteny stage, not a statement about
  positional orthology.

That is the finding the 39 interleaved sequences were pointing at, and it is why
they must stay in the record. They are not noise to be pruned so the tree looks
monophyletic; they are 39 independent demonstrations that a selection rule
written for redundancy control is being asked to do orthology assignment it
cannot do.

## 5. What this does NOT establish

- **It does not show the retained copies are wrong.** Gene order favours them
  on average (18 pairs against 11, median 19 slots against 17, collinearity
  ρ ≥ 0.9 in 22 against 14). It shows the evidence is not decisive per species.
- **It does not establish how many sheath genes each genome has.** "Single-copy"
  below means one copy inside the frozen 30–90% identity window; a second copy
  more similar than 90% or more distant than 30% would not have been retrieved.
  Settling that needs a per-genome sheath census, which was not run.
- **No rates were estimated and Rate4Site was not run.**

## 6. Options for the manuscript, with what each costs

**Option A — single-copy-species restriction.** Define the ortholog set as the
members from species contributing exactly one gated candidate, and report the
multi-copy species separately as the paralogy evidence they are. Nothing is
deleted; the duplicated species are described rather than silently dropped.

| boundary | members | from single-copy species | genera | identity |
|---|---|---|---|---|
| L2 | 124 | 85 (69%) | 2 | 0.65–0.83 |
| L4 | 155 | 106 (68%) | 11 | 0.64–0.83 |
| L7 | 293 | **211 (72%)** | **75** | 0.30–0.83 |
| L9 | 355 | 269 (76%) | 112 | 0.30–0.83 |

L7-restricted is the only combination in this study that has both a defensible
orthology claim and the substitution signal the site-rate analysis needs (the
wide-boundary diagnostics: 6 invariant sites of 143, median entropy 1.50 bits,
21 effective genera). Its cost is the wide boundary's alignment inflation
(2.21 columns per reference residue) and the unverified assumption that
single-copy-in-window means single-copy-in-genome.

**Option B — paralogy-aware treatment.** Keep both copies and model the
duplication: gene-tree/species-tree reconciliation to label co-orthologs and
in-paralogs explicitly, with rates estimated on the reconciled tree. Scientifically
the stronger answer, and it uses all 39 interleaved sequences as data. It is also
a new methodological component for the manuscript and would need its own
validation.

**Option C — R2 does not enter the primary evolutionary analysis.** Its `R_i`
values are structural and stand regardless; R2 would contribute to Result 1
(reproducible contact remodeling) and not to Result 3. This is the honest
fallback if neither A nor B is acceptable, and it is not a failure of the
structural work.

**Not an option:** retaining "124 orthologs" or "293 orthologs" selected by best
bitscore and describing them as an ortholog set.

## 7. Recommended next step

Option A at L7, restricted to single-copy species, with the per-genome sheath
census run first to test the single-copy assumption — that census is the one
cheap check that would convert Option A from defensible to verified. Option B if
you want the duplication modelled rather than avoided. Either way the same
question must be asked of SalCIS, diffocin and the AFP pair before they are
treated as ortholog sets, since 23% of the species sampled here carry multiple
copies within the identity window.

## Files

| file | content |
|---|---|
| `locus_architecture.py` | ordered neighbourhood extraction, slot indexing, comparison primitives |
| `positional_orthology_r2pyocin.py` | the analysis driver |
| `orthology/r2pyocin_positional_loci.csv` | per-locus slots recovered, exact matches, collinearity, regulators |
| `orthology/r2pyocin_positional_pairs.csv` | the 39 pairs side by side |
| `orthology/r2pyocin_positional_tracks.json` | every gene of every locus with its assigned reference slot |
| `orthology/r2pyocin_reference_track.json` | the reference locus, 57 CDS |
| `orthology/locus_cache.json` | fetched neighbourhoods, so a rerun needs no refetch |
