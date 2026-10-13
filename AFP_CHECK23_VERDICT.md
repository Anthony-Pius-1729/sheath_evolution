# AFP — checks 2 and 3, and the orthology verdict

**Date** 2026-10-13 · **Anchor** reference pADAP locus, Afp2 slot 0, Afp3 slot
+1, third sheath gene slot +2 · **Scope** strictly diagnostic: no ortholog set
selected, no MSA QC, no rates, nothing frozen modified

## Verdict

**AFP orthology: PASS on a restricted subset.** Positions assigned from locus
architecture recover coherent lineages, but only where the locus can be
anchored. By the agreed criteria this is the second outcome — *some confidently
assignable loci, others ambiguous → restrict to the supported subset, document
exclusions* — and the exclusions are substantial and recorded, not quietly
dropped.

| | |
|---|---|
| loci characterised | 103 |
| anchored on non-sheath flanking homology | **57 (55%)** |
| sheath genes positionally assigned | **82** (22 Afp2-position, 53 Afp3-position, 7 third-position) |
| unassigned | 61 (58 in unanchored loci, 3 anchored but outside the mapped slots) |
| nearest-labelled-neighbour concordance in the joint tree | **94.1%** (5 discordant of 85) |
| loci with both an Afp2 and an Afp3 position assigned | **21**, from 21 genera |

The decisive contrast with R2: there, every criterion available failed to
separate the copies. Here position and phylogeny agree at 94%, and the
disagreements are individually interpretable.

## Check 2 — locus characterisation

Candidates were enriched by a screen symmetric over all three reference sheath
genes (bitscore to the best of the three ≥ 200), sampled across genera, and
resolved to loci. Each locus was anchored using **only non-sheath flanking
genes** with homology to the reference locus: a modal candidate-slot minus
reference-slot offset, requiring ≥3 anchors, ≥60% agreeing on the offset and
collinearity ρ ≥ 0.8. Anchored loci have a median of 9 non-sheath anchors and
median collinearity ρ = 1.00.

**Sheath-run length among anchored loci:** 35 have one sheath gene, 16 have two,
6 have three. Among unanchored loci: 41 have one, 2 have two, 2 have three, 1
has five. So the tandem run is real but not universal — roughly 40% of anchored
AFP-like loci carry more than one sheath gene.

**Why best-hit labelling would have failed.** Of 5,470 pool members with a hit
to any of the three references, **4,928 best-hit Afp3** against 440 for Afp2 and
102 for the third gene. Afp3 is an attractor for the whole family. Measured
against the positional assignment:

| | seq best hit → Afp2 | → Afp3 | → third |
|---|---|---|---|
| **position = Afp2** | 7 | **15** | 0 |
| **position = Afp3** | 0 | 52 | 1 |
| **position = third** | 0 | 0 | 7 |

**16 of 82 disagree (20% overall), and the error is one-directional and
concentrated almost entirely on Afp2: 15 of the 22 genes at the Afp2 position —
68% of them — would have been labelled Afp3-like by sequence, against 1 of 53
at the Afp3 position and 0 of 7 at the third.** The 20% figure is the pooled
rate and understates the Afp2-specific risk by more than threefold.
That is the silent interchange the brief warned about, quantified. Length does
not rescue it either — genes at the Afp2 position have median length 477 aa
against the *S. entomophila* Afp2's 354 aa, so the reference is itself atypical.

## Check 3 — one joint phylogeny

A single tree of reference Afp2, reference Afp3, the third reference gene and
all 142 recovered sheath copies, assigned and unassigned alike. Positional
labels are annotation on the tips; the tree knows nothing about them.

![AFP joint tree]({{artifact:art_f38394f1-7603-46a1-9400-c8d016ee6539}})

**Nearest-labelled-neighbour concordance 94.1%** — Afp3-position 96.3%,
Afp2-position 91.3%, third-position 87.5%. Within-label patristic distance is
below between-label in every class: Afp2 3.36 against 4.55, Afp3 3.08 against
5.06, third 3.02 against 6.44.

**Global monophyly is absent, and that is the expected result, not a failure.**
No label is monophyletic across the whole tree — the MRCA of the 23
Afp2-position tips spans 133 tips. But in tree order the labels fall into a
small number of large contiguous blocks (Afp2-position: 8, 7, 7, 1;
Afp3-position: 33, 8, 6, 4, 2, 1; third: 7, 1), and those blocks pair up **by
taxonomy**:

| block | label | n | genera |
|---|---|---|---|
| 3 | third | 7 | *Chania, Erwinia, Pararobbsia, **Photorhabdus**, Serratia, Xenorhabdus, Yersinia* |
| 4 | Afp3 | 8 | identical to block 3 |
| 5 | Afp2 | 8 | same, except ***Dyella*** replaces *Photorhabdus* |
| 9 | Afp2 | 7 | *Algoriphagus, Dyadobacter, Fulvivirga, Mangrovibacterium, Marinigracilibium, **Pedobacter**, Pseudozobellia* |
| 10 | Afp3 | 6 | identical to block 9 less *Pedobacter* (a strict subset) |

The genus composition is **close but not identical**: blocks 3, 4 and 5 share
6 of the 8 genera in their union, with *Photorhabdus* present in the third- and
Afp3-position blocks but not the Afp2-position block, and *Dyella* the reverse.
Blocks 9 and 10 are a strict subset relation. Some of this is sampling — one
sheath gene per locus can fail to anchor while its neighbour succeeds, so a
genus can appear at one position and not another — and the point the pairing
supports is that the blocks are drawn from the *same taxonomic group*, not that
they are genus-identical.

The gene tree groups first by **AFP subfamily** — enterobacterial, and at least
two Bacteroidota groups — and within each subfamily the Afp2, Afp3 and third
positions each form their own block drawn from the same taxonomic group
(near-identical genus sets; see the caveat under the table). That is the
signature of a conserved multi-gene locus evolving vertically within a
subfamily. The Afp2 position in a Bacteroidota locus is not the closest relative
of the Afp2 position in *Serratia*, and there is no reason it should be: the
whole locus diverged. **Correspondence is well defined within a subfamily and
not across the family.**

**The five discordant tips, investigated rather than removed.** Two are a
reciprocal pair from one locus (WP_131557722 / WP_131557724) in which the Afp2-
and Afp3-position genes are each other's nearest labelled neighbour — a recent
within-locus duplication, where the two positions are filled by close paralogs
rather than by the two ancient lineages. Two more are a reciprocal Afp3/third
pair. One is a singleton Afp2-position gene nearest an Afp3-position gene. None
is a case of the positional method producing an incoherent assignment at scale.

**The third sheath gene** is a distinct lineage, not a duplicate of either: its
nearest neighbour is another third-position gene, the third-position block is
the purest of the three (0.89 among labelled tips), and its within-label
distance is far below its between-label distance (3.02 against 6.44). It does
not need its history reconstructed to be set aside — it is a separate position
with its own lineage and should simply not be merged into Afp2 or Afp3.

## What this does and does not license

**Licensed.** Within an AFP subfamily, Afp2-position and Afp3-position sequences
can be identified without arbitrary sequence selection, by locus architecture,
with phylogeny concurring at 94%.

**Not licensed.** A single family-wide "Afp2 ortholog set" pooling all 22
Afp2-position sequences. They are not one lineage, and pooling them across
subfamilies would repeat the R2 error in a new form.

**The open question for the next gate** — not answered here, and not to be
pre-judged: the reference's own subfamily contributes 8 Afp2-position and 8
Afp3-position sequences from 7 genera. If the analysis must be per-subfamily to
be defensible, the per-family sample size for rate estimation will be small. The
largest Bacteroidota block offers 7 + 6. Whether that supports site-rate
estimation at the 84 structural sites of Afp2 and 82 of Afp3 is the feasibility
question, and it should be asked on the restricted subsets, not on the pooled
set.

## Manuscript disclosure, carried forward

The entropy–`R_i` correlation computed earlier (ρ ≈ 0.32 on R2) means the
evolutionary hypothesis cannot be described as completely unexamined. Sequence
entropy is not Rate4Site, but it is an evolutionary variability measure related
to the intended outcome. The manuscript must distinguish three tiers:

1. **exploratory analyses already performed** — the structural pilot, the R2
   retrieval and orthology investigation, and the entropy–`R_i` correlation;
2. **analysis decisions fixed before any Rate4Site run** — the frozen `R_i`
   definition and domain, the 80% occupancy threshold, the orthology criteria
   established here, and the inference specification still to be settled;
3. **genuinely held-out evolutionary tests** — whatever is run on families and
   sites not used to make those decisions.

A Rate4Site analysis remains possible and legitimate; it must be reported inside
that structure rather than as an untouched prediction.

## Files

| file | content |
|---|---|
| `afp_loci.py` | locus characterisation and positional assignment |
| `afp_tree.py` | the joint phylogeny and its three concordance statistics |
| `orthology/afp_loci.csv` | one row per locus: run length, strands, intergenic distances, anchors, offset, collinearity |
| `orthology/afp_sheath_runs.csv` | one row per sheath gene: positional label, reason, and the sequence best hit recorded separately |
| `orthology/afp_locus_tracks.json` | every CDS of every locus with its reference mapping |
| `orthology/afp_check3.json` | concordance, purity and distance statistics |
| `orthology/afp_tree_tips.csv` | per-tip nearest labelled neighbour |
| `orthology/afp_joint.{afa,nwk}` | the joint alignment and tree |
| `orthology/afp_prescreen.csv` | the symmetric three-reference screen used for sampling |
