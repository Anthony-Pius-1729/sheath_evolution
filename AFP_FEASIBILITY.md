# AFP evolutionary feasibility decision

**Family:** *Serratia entomophila* antifeeding prophage (AFP), PDB 6RBN entities 2 (Afp2, 354 aa)
and 3 (Afp3, 451 aa). **Date:** 2026-10-13. **Verdict: PASS — two usable per-subfamily
datasets, counting as one independent system.**

This closes the AFP branch. It uses only sequences already retrieved (the 8,005-record
UniProt pool from the shared Pfam query); no new database retrieval was performed.

---

## 1. What the orthology work had established

`AFP_CHECK23_VERDICT.md` licensed **per-subfamily** Afp2 and Afp3 sets and refused
pooled family-wide sets. Two facts forced that:

- Afp2 and Afp3 are only **38.1% identical to each other**, inside the frozen 30–90%
  retrieval window, so a single family-wide query returns both and sequence similarity
  cannot separate them. Best-hit labelling disagreed with positional assignment for
  **15 of 22 (68%)** genes at the Afp2 position.
- In the joint 145-tip tree no label was monophyletic, and of the 29 nodes inside the
  afp2 MRCA that mix afp2 with non-afp2 tips, **9 carry support ≥ 0.9** (median 0.79).
  Pooling is contradicted by well-supported nodes, not merely unresolved.

Positional orthology, by contrast, worked: nearest-labelled-neighbour concordance in the
joint tree was 94.1%, and the three tandem genes occupy reproducible slots in an
anchored, well-annotated locus.

## 2. Correction to the first feasibility pass

A first pass scored the reference-containing subfamily at **n = 8** per gene. That number
was an artifact of my own sampling, not a property of AFP: the family-wide locus screen
dereplicated candidates **one per genus**, and the reference subfamily occupies only eight
genera. The already-retrieved pool contains **279 screened members (bitscore ≥ 200) across
45 species** in those genera.

Re-running the identical locus-characterisation pipeline restricted to those genera and
dereplicated **by species** gives the real set size. The n = 8 rows in
`orthology/afp_feasibility.csv` are superseded by `orthology/afp_sf1_feasibility.csv`.

## 3. The reference subfamily locus set

42 loci characterised, **31 anchored** on non-sheath flanking homology (74%, against 55%
family-wide). Sheath-run lengths: **35 of 42 loci carry a run of three** tandem sheath
genes — the Afp2/Afp3/third architecture of the reference is the subfamily norm, not a
*Serratia* peculiarity. Positional assignment yielded 27 Afp2-position, 31 Afp3-position
and 29 third-position genes, with 45 unassigned.

Within this subfamily, position and sequence best-hit agree for **55 of 58** genes (95%),
against 66 of 82 (80%) family-wide and 7 of 22 (32%) at the Afp2 position family-wide.
The labelling failure documented in check 2 is a cross-subfamily phenomenon; inside the
reference subfamily the two lines of evidence converge. Positional assignment remains the
operative criterion — the agreement is recorded, not used as a rule.

## 4. Feasibility at the structural sites

Each set aligned with its own reference included (FAMSA), tree by FastTree `-lg -gamma`,
scored at the in-domain `R_i` sites of the corresponding 6RBN entity.

| | **Afp2** | **Afp3** |
|---|---|---|
| sequences (= species) | 27 | 31 |
| genera / effective genera | 8 / 5.49 | 8 / 5.43 |
| largest genus share | 0.333 | 0.323 |
| alignment columns | 415 | 708 |
| columns per reference residue | 1.17 | 1.57 |
| structural sites scored | 84 | 82 |
| sites at ≥80% occupancy | 84 / 84 | 81 / 82 |
| **invariant structural sites** | **33 (39%)** | **17 (21%)** |
| sites with ≥95% majority residue | 35 | 24 |
| median per-site entropy (bits) | 0.554 | 0.823 |
| maximum possible entropy at this n | 4.75 | 4.95 |
| tree length | 6.62 | 17.02 |
| mean / max patristic distance | 0.914 / 1.683 | 1.937 / 5.550 |

**Informative sites: 51 of 84 (Afp2), 65 of 82 (Afp3).**

### Reading these numbers

*Both sets clear the gates that R2 failed.* Positional orthology is established, the
alignment is near-unambiguous (1.17 and 1.57 columns per reference residue, against 2.2–2.7
in R2's wide sets), occupancy is complete at essentially every structural site, and the
sampling is one sequence per species.

*Afp3 is the stronger set.* 31 sequences, 79% of structural sites variable, median entropy
0.823 bits. Its max patristic distance of 5.55 substitutions/site signals that the deepest
pairs approach saturation; the 1.57 columns per reference residue reflects indel-rich
regions outside the structural sites, which the ≥80%-occupancy filter already handles.

*Afp2 is usable but thinner.* 27 sequences and **39% of its structural sites carry no
substitution at all**. Those 33 sites cannot contribute to a rate–`R_i` association: they
will all receive the lower bound of the rate scale with wide posteriors. The association
test for Afp2 effectively runs on 51 sites, and with a tree length of 6.62 the rate
estimates at the remaining sites are correspondingly imprecise.

*Taxonomic breadth is the real weakness of both.* 5.4–5.5 effective genera, with a third of
each set in one genus. R2's wide sets reached 21–35 effective genera. AFP's rates therefore
describe substitution within one enterobacterial radiation, and the shared species history
is not averaged out.

## 5. Independence

Afp2 and Afp3 come from **the same loci in the same genomes**. They share the species tree,
the locus history and the selective regime of the same particle. They are two proteins, not
two samples: usable as two per-protein replicates of the within-protein site-level
association, but **not as two independent systems** for any cross-system claim. AFP
contributes **one** system to the comparative design.

## 6. Verdict

**PASS, on two per-subfamily datasets, with declared limits.**

| | |
|---|---|
| Afp3-position set, n = 31 | **suitable** for site-specific rate estimation |
| Afp2-position set, n = 27 | **suitable with reduced power** — 51 informative sites of 84 |
| pooled family-wide sets | **not licensed** (§1); present in `afp_feasibility.csv` labelled as such, for the record only |
| independence | the two sets count as one system |

Limits to carry into the manuscript: a single enterobacterial radiation (5.4 effective
genera), 39% invariant structural sites in Afp2, deep-pair saturation in Afp3, and
positional rather than sequence-based orthology assignment.

## 7. What was deliberately not done

No new retrieval, no additional structural analysis, no exploratory correlation between any
AFP rate proxy and `R_i`, no gene-tree reconciliation, no attempt to rescue the pooled sets
or to widen the subfamily by relaxing the anchoring or bitscore criteria. The 45 unassigned
genes and the 29 third-position genes are retained in `afp_sf1_sheath_runs.csv` and were not
pruned.

## 8. Files

`orthology/afp_sf1_feasibility.csv` (the table above) · `orthology/afp_sf1_loci.csv` ·
`orthology/afp_sf1_sheath_runs.csv` · `orthology/afp_sf1_locus_tracks.json` ·
`orthology/afp_sf1_{afp2_like,afp3_like}.{afa,nwk}` (the two analysis-ready alignments and
trees) · superseded: the `SF1_*` n = 8 rows of `orthology/afp_feasibility.csv`.

Reproduce: `python afp_loci.py --genera
"Serratia,Xenorhabdus,Photorhabdus,Yersinia,Erwinia,Chania,Pararobbsia,Dyella" --n 60
--out-prefix afp_sf1`
