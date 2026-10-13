# Single-family retrieval validation — R2 pyocin

**Date** 2026-10-13 · **Family** `r2pyocin` · **Reference** UniProt Q9S574,
PDB 3J9Q, *Pseudomonas aeruginosa*, 386 aa
**Route** UniProtKB Pfam-family scan → local blastp → frozen gates
**Status** validation of the retrieval route for ONE family. Not an orthology
decision, not a phylogeny, no rates. The other four families have not been run.

## Verdict against the required evidence

| requirement | evidence | where |
|---|---|---|
| Real homologs retrieved | 8,035 bacterial candidates → 4,593 blastp subjects at E ≤ 1e-5 → 1,266 through the frozen gates → **967 species-dereplicated sequences across 280 genera** | `filtered/r2pyocin_retained.fasta`, `filtered/r2pyocin_filtered.csv` |
| Reference recovered | Q9S574 **retrieved** in the candidate set; self-hit 100.0% identity over 100% of the reference; **excluded by the 90% ceiling**, which is the gate working, not a miss | `retrieval_uniprot/r2pyocin_search_provenance.json` → `reference_recovery` |
| Sequence diversity | identity median **0.352**, range **0.300–0.826**; reference coverage median **0.984**, min **0.806** | `evolution_prelim/r2pyocin_identity_coverage.png` |
| Genomic provenance | **938/967** carry a RefSeq accession; **211/250** sampled candidates resolved to a genomic locus through NCBI Identical Protein Groups | `synteny/r2pyocin_synteny.csv` |
| Candidate orthology | preliminary FastTree LG+Γ over 968 taxa; reference sits in a nested series of clades, **211/250** scored loci carry ≥2 contractile-machine marker categories | `evolution_prelim/r2pyocin_{prelim.nwk,clade_profile.csv,summary.json}` |
| Reproducibility | UniProt release, exact query, Pfam signatures, BLAST+ version and parameters, frozen gate constants, file checksums all recorded | `retrieval_uniprot/r2pyocin_{provenance,search_provenance}.json` |

## Why the route changed

The NCBI BLAST URL API is not returning. A single blastp of this reference
against `refseq_protein` sat at `Status=WAITING` for 131 min across 17 polls on
2026-10-12, and a control of the same query against `swissprot` — small enough
to answer in under a minute under normal load — also stayed `WAITING` for 22 min
across 37 polls. `blast_retrieve.py` is retained as the fallback and is correct;
the endpoint is the blocker.

Nothing in the scientific design moved. The identity window, coverage floor,
length window and dereplication rule are **imported** by `local_search.py` from
`filter_hits.py`, so there is exactly one definition of each and the frozen
protocol cannot drift by being re-typed.

## What was run

```
uniprot_retrieve.py --families r2pyocin     # candidate set from UniProtKB
local_search.py     --families r2pyocin     # local blastp + frozen gates
synteny.py                                  # unchanged frozen stage 4
prelim_tree.py      --family  r2pyocin      # preliminary QC tree
```

**Retrieval.** The family definition is taken from the reference's own UniProt
entry, never hardcoded: Q9S574 carries Pfam PF04984, PF17482 and PF22671, giving
the query `(xref:pfam-PF04984 OR xref:pfam-PF17482 OR xref:pfam-PF22671) AND
taxonomy_id:2` against UniProt release **2026_03 (02-September-2026)** → 8,035
records, of which 7,489 carry a RefSeq cross-reference and 8,032 an EMBL one.
The script verifies that the UniProt sequence is byte-identical to
`reference_sequences.json` before retrieving anything, and aborts if not.

**Search.** BLAST+ 2.17.0+ `makeblastdb` + `blastp`, E ≤ 1e-5, `-seg no`,
`-max_target_seqs 20000`, BLOSUM62. Identity and coverage are computed exactly
as `filter_hits.parse` computes them from BLAST XML — identity = nident/align_len,
coverage = (qend−qstart+1)/reference_length, best HSP per subject by bitscore —
so the gates act on the same quantities as they would have on the NCBI route.

**Gates.** 4,593 → 1,266 (30–90% identity, ≥80% reference coverage, 0.75–1.30×
length) → 967 after one-representative-per-species on best bitscore.

## Locus evidence

250 candidates (the frozen `MAX_PER_FAMILY`, spread across genera and identity
deciles by the unchanged `select()`) went to the synteny stage. 211 resolved to
a genomic locus; all 211 scored ≥2 marker categories, and the 39 that scored
zero are exactly the 39 whose locus did not resolve — no candidate resolved to a
locus and then failed to look like a contractile machine.

| marker categories | 0 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| candidates | 39 | 2 | 9 | 78 | 90 | 30 | 2 |

The most common neighbourhoods are `baseplate; phage_tail; tape_measure; tube`
(62) and the same plus `fibre` (51) — the expected contractile-tail architecture
around the sheath gene.

## Preliminary tree

FAMSA alignment of 968 sequences (1,310 columns), FastTree LG+Γ, midpoint
rooted. The nested clades from the reference outwards:

| level | clade size | stem | support | identity min–median–max | synteny-supported / scored |
|---|---|---|---|---|---|
| 0 | 7 | 0.072 | 0.89 | 0.782 / 0.795 / 0.826 | 0 / 0 |
| 1 | 119 | 0.060 | 0.93 | 0.704 / 0.753 / 0.826 | 2 / 2 |
| 3 | 146 | 0.056 | 0.40 | 0.648 / 0.749 / 0.826 | 4 / 6 |
| 6 | 161 | 0.106 | 0.85 | 0.317 / 0.746 / 0.826 | 7 / 12 |
| 7 | 293 | 0.093 | 0.89 | 0.301 / 0.693 / 0.826 | 45 / 56 |
| 9 | 355 | 0.537 | 1.00 | 0.301 / 0.451 / 0.826 | 66 / 80 |
| 10 | 967 | — | — | 0.300 / 0.352 / 0.826 | 219 / 258 |

(full table in `r2pyocin_clade_profile.csv`; the four-way
clade × synteny classification is tabulated at every level in
`r2pyocin_summary.json` under `four_way_by_level_PROVISIONAL`.)

**No orthology cut has been made.** The frozen rule is tree placement combined
with locus context as a four-way classification, and the clade boundary is a
decision to take on this evidence — so the profile is reported and nothing is
fixed. The long stem at level 5 (1.157) and the identity drop between levels 6
and 7 (median 0.746 → 0.693, minimum 0.317 → 0.301) are the two features that
boundary decision will have to address.

## Caveats carried forward

1. **Different recall from BLAST.** A Pfam-family scan retrieves entries
   carrying the reference's signatures. A genuine ortholog annotated without any
   of them is invisible to this route, where BLAST would have found it on
   sequence similarity; conversely this route reaches entries BLAST ranked below
   its hit-list cutoff. The two are different retrieval routes, not substitutes.
   Any dataset must be reported with the route that produced it.
2. **All 967 are unreviewed (TrEMBL) entries.** Expected for bacterial
   prophage-locus proteins, but no manual curation stands behind the
   annotations; the synteny evidence, not the protein name, is what identifies
   the locus.
3. **Synteny covers 250 of 967**, by the frozen cap. The 84% locus-resolution
   rate is a property of that sample, not of the whole set.
4. **`accession` is not unique in the filtered table.** Four RefSeq accessions
   are shared by 2–5 UniProt entries from different species (e.g.
   `WP_036782415.1` ×5), so the 250 scored loci attach to 258 leaves. Identical
   proteins in distinct species survive species-level dereplication
   independently; downstream code must not assume one row per accession.
5. **The tree is a QC screen.** FastTree LG+Γ, midpoint-rooted, one run, no
   bootstrap beyond FastTree's local supports. Model selection and rigorous ML
   inference remain open, as frozen.

## Next

Scale the identical workflow to `afp_Afp2`, `afp_Afp3`, `salcis`, `diffocin`,
then inspect all five datasets together before any rate estimation. Rate4Site
stays unrun.
