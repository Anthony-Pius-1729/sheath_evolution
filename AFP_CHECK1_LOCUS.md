# AFP — check 1: reference locus architecture

**Date** 2026-10-13 · **References** Afp2 = Q6HAD7 / WP_010895804.1 (354 aa),
Afp3 = Q6HAD6 / WP_010895805.1 (451 aa), both PDB 6RBN entities 2 and 3,
*Serratia entomophila* · **Locus** NC_002523.5 (pADAP), 36 CDS in ±25 kb

## Correspondence to the structural references: verified

Both RefSeq proteins are **100% identical** over their full length to the
sequences in `reference_sequences.json`, which were taken from the 6RBN
entities — Afp2 727 bits over 354 aa, Afp3 916 bits over 451 aa. The structural
and genomic references are the same molecules.

## The locus

Afp2 and Afp3 are **adjacent, same strand, 20 bp apart** (75,363–76,427 and
76,447–77,802). Slots are signed from Afp2:

| slot | protein | category | product |
|---|---|---|---|
| −6 | WP_010895797.1 | other | phage holin, lambda family |
| −4 … −2 | WP_169558216.1, WP_010895801.1, WP_010895802.1 | other | lysis protein, holin, lysozyme |
| −1 | WP_010895803.1 | phage_tail | phage tail protein (149 aa) |
| **0** | **WP_010895804.1** | phage_tail | **Afp2**, phage tail sheath family |
| **+1** | **WP_010895805.1** | phage_tail | **Afp3**, phage tail sheath family |
| **+2** | **WP_169558217.1** | phage_tail | **a third sheath-family protein, 415 aa** |
| +3 | WP_380740010.1 | phage_tail | phage tail protein (150 aa) |
| +5 | WP_010895809.1 | tube | CIS tube protein |
| +6 | WP_010895810.1 | baseplate | phage baseplate assembly protein V |
| +8 | WP_010895812.1 | spike | contractile injection system spike tip |
| +9, +10 | WP_010895813.1, WP_010895814.1 | baseplate | baseplate J/gp47, baseplate wedge |
| +12 | WP_169558219.1 | tape_measure | CIS tape measure protein |
| +13, +14 | WP_010895817.1, WP_169558220.1 | other | ATP-binding protein, Pvc16 family |

**This locus is well annotated**, in sharp contrast to the R2 reference, whose
sheath gene is called "bacteriophage protein" and whose neighbourhood matches
none of the frozen marker categories. Here `tube`, `baseplate`, `spike` and
`tape_measure` are all present and correctly named, so the frozen synteny score
would behave as intended on AFP. The gene-order-by-homology comparison built for
R2 (`locus_architecture.py`) still applies and remains the method of record —
annotation text helps interpretation, it does not define orthology.

## Two findings that reshape check 2

**1. There are three tandem sheath-family genes, not two.** WP_169558217.1 at
slot +2 is a 415 aa sheath-family protein whose best hits are
**140 bits to Afp3 and 136 bits to Afp2** — essentially tied. Any assignment of
this protein to "the Afp2 lineage" or "the Afp3 lineage" by sequence similarity
alone is a coin flip, in the reference locus itself. The frozen protocol's
treatment of AFP as a two-protein paired locus is therefore incomplete before
any candidate is retrieved.

**2. Afp2 and Afp3 fall inside each other's retrieval window.** They are
**38.1% identical** over a 452-residue alignment (267 bits). The frozen
retrieval window is 30–90% identity. So:

- a search seeded with Afp2 will retrieve Afp3-like proteins as candidates, and
  vice versa;
- the frozen gates cannot separate them, because separating them is exactly the
  30–90% band the gates are designed to keep;
- the third in-locus paralog sits at 29.5–38.8% to both, i.e. on the boundary of
  the same window.

Both reference entries carry the **same Pfam signatures** (PF04984 + PF17482)
and the same InterPro sets, so `uniprot_retrieve.py` generates an identical
candidate query for the two families — one candidate pool, three in-locus
sheath lineages inside it, mutually within the retrieval window.

Note that Afp2 and Afp3 carry PF04984 + PF17482 but **not** PF22671, which R2's
Q9S574 does carry. The AFP candidate pool is therefore defined by a broader
two-signature query than R2's.

## The copy-number gate is blind to this locus

The cheap census proposed in `R2_FINAL_DECISION.md` was run on the AFP pool
(8,005 records from `(PF04984 OR PF17482) AND Bacteria`, 5,126 proteomes). Its
distribution is indistinguishable from R2's — 65.2% single detected copy, median
1, max 19 against R2's 65.1%, median 1, max 20 — because the two queries return
very nearly the same pool.

More importantly, **the census sees none of this locus's three sheath genes**:

- Afp2 (Q6HAD7) and Afp3 (Q6HAD6) are in the pool but have **no `Proteomes`
  value** — they are plasmid-borne pADAP entries, and 27 of 8,005 records are in
  that state;
- the third in-locus sheath gene WP_169558217.1 has **no entry in the pool at
  all**;
- the *S. entomophila* reference proteome UP001056873 reports **one** detected
  copy, a 389 aa protein that is none of the three.

So the census is a positive-only screen with no negative predictive value, and
on a plasmid-borne system it can be blind to the entire locus. That limit has
been written back into the R2 carry-forward protocol. For AFP the census is
simply uninformative, and multi-copy detection would not have meant what it
means for R2 anyway: the AFP locus is *expected* to carry several sheath genes,
and treating that as a paralogy hazard would be the error the brief warns
against.

## Consequence for the protocol

Check 2 cannot be run as two independent single-reference searches with
best-hit assignment — that procedure would mix the three lineages by
construction, and it is the same failure that ended R2, arriving here by a
different route. The assignment has to come from **position within a homologous
locus**: find candidate loci that contain a tandem run of sheath-family genes,
establish the run's homologous correspondence to the reference run
(Afp2 at 0, Afp3 at +1, third at +2) by neighbour homology and gene order, and
read the lineage assignment off the position rather than off the similarity
score.

That is a stronger test than anything available for R2 — R2's single sheath gene
had no positional constraint to exploit — and it is the reason AFP was worth
taking first. It is also a sharper risk: a mis-assignment here does not merely
mix paralogs, it swaps one protein of the machine for another.

The critical question for check 3 is unchanged and now precisely posed: **can
corresponding Afp2 and Afp3 lineages be identified across homologous AFP loci by
locus position, without arbitrary sequence selection?** Check 2 will establish
how often the tandem run is present and intact; check 3 will test whether
position and phylogeny agree.

## Files

| file | content |
|---|---|
| `orthology/afp_reference_track.json` | the 36-CDS reference locus with slots, strands, products |
