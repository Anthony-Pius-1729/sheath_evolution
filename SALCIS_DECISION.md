# SalCIS orthology and feasibility decision

**Reference:** PDB 9R9H entity 1, *Salmonella enterica* subsp. *salamae*, 519 aa,
UniProt A0A701XPN5. **Date:** 2026-10-13. **Verdict: FAIL — no ortholog set is
definable. SalCIS is excluded from the evolutionary test.**

---

## 1. Reference identity and locus

The deposited structure sequence is a **strain variant** of A0A701XPN5: same length,
**4 substitutions, 99.229% identical** (G41C, G175D, H215N, S452G), of which one (position
41) falls inside the 68 in-domain `R_i` sites. The UniProt entry is used only to define the
Pfam family query (PF17482); the blastp query and the `R_i` mapping both use the PDB
sequence, so the structural analysis is unaffected. The substitutions are recorded in
`retrieval_uniprot/salcis_provenance.json`.

The reference protein has **no RefSeq cross-reference** and no UniProt proteome. It sits on
a TPA whole-genome-shotgun contig (DAAMGA010000009.1, *S. enterica* subsp. *salamae* strain
3940-62, 192 kb, 162 CDS), which had to be fetched and parsed directly from INSDC.

The reference locus itself is good: **42 CDS in ±25 kb, exactly one sheath gene**, flanked
downstream by phage tail proteins (+1, +2), a peptidoglycan-binding protein (+3), **VgrG
(+4)** and a T6SS protein (+5), and upstream by conserved housekeeping genes (YchF,
ribose-phosphate pyrophosphokinase, DauA) that mark a clean element boundary. On
architecture alone this is the most favourable reference in the project.

## 2. What the retrieval actually returned

Frozen pipeline: 7,592 PF17482 bacterial records → 4,490 blastp hits → 805 through the
frozen gates → **667 species-dereplicated across 246 genera**, identity median 0.413
(0.300–0.879), coverage median 0.994.

The composition is decisive:

- **Of 4,490 blastp hits, exactly one is a *Salmonella* sequence — the reference itself**
  (99.2%, the UniProt copy, correctly excluded by the 90% ceiling). No second *Salmonella*
  sheath homolog is detected at any identity. The absence of congeners is not an artifact
  of the identity ceiling.
- The gated set's largest genera are **Bacteroidota**: *Chryseobacterium* (37),
  *Flavobacterium* (32), *Algoriphagus* (30), *Chitinophaga* (30), *Pedobacter* (24),
  *Mucilaginibacter* (22). The reference is a Gammaproteobacterium.
- Only **33 of 1,070** screened candidates are Enterobacterales at all, and none is
  *Salmonella*.

The frozen criteria therefore return a cloud of distant homologs from another phylum, not a
set of orthologs of the reference gene.

## 3. Positional orthology

Homology-anchored, ordered neighbourhood comparison (the corrected synteny method; same
code and the same constants used for AFP and diffocin):

| | family-wide (55 loci) | Enterobacterales only (29 loci) |
|---|---|---|
| anchored on non-sheath flanking homology | **9 / 55 (16%)** | **9 / 29 (31%)** |
| sheath-run length 1 | 27 | 19 |
| sheath-run length 2 | **28** | 3 |
| sheath-run length 3 | 0 | **7** |
| positional assignments at the reference slot | **3** | **0** |
| modal offsets among anchored loci | 0 (×3), +1 (×4), +2 (×2) | +1 (×9) |

Two independent failures:

1. **The locus is not recoverable.** In 84% of sampled loci, fewer than three of the
   reference's non-sheath neighbours have a detectable homolog in the candidate
   neighbourhood, so no frame can be established. For AFP the same procedure anchored 74%
   of loci inside the reference subfamily.
2. **The architecture differs where it can be read.** The reference carries one sheath
   gene; **28 of 55** sampled loci carry two tandem sheath genes, and 7 of the 29
   Enterobacterales loci carry three. The copy-number structure of the candidate loci is
   not that of the reference.

The three family-wide assignments and the nine Enterobacterales loci at a uniform +1 offset
are too few to build a dataset under any reading of the rule: relaxing the offset
requirement to the collinear frame (ρ ≥ 0.8, as reported for diffocin) still yields at most
**9 positional orthologs**, against a frozen minimum of 30.

## 4. Verdict

**FAIL.** No ortholog set of the SalCIS reference sheath can be defined with the
sequence- and context-based criteria in use. The failure is not marginal and not a sampling
shortfall: the reference has no detected congeneric homolog, the recovered homolog cloud
belongs to a different phylum, the reference locus is unrecoverable in 84% of candidate
loci, and half the candidate loci have a different sheath copy number.

SalCIS is **excluded from the evolutionary test** for the same reason as R2 pyocin —
indefinability of an orthologous lineage, not inconvenience. Its structural `R_i` values
(68 in-domain sites) remain valid and stand in Result 1.

No rescue was attempted beyond the clade-restricted re-test, which failed harder than the
family-wide run.

## 5. Files

`orthology/salcis_reference_track.json` · `orthology/salcis_{loci,sheath_runs,locus_tracks}`
· `orthology/salcis_ent_{loci,sheath_runs,locus_tracks}` · `filtered/salcis_filtered.csv` ·
`retrieval_uniprot/salcis_{candidates.tsv,blastp_hits.csv,provenance.json}`
