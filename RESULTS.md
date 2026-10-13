# Results

**Question.** Does contact-network remodeling during contraction explain site-specific
evolutionary-rate variation beyond solvent accessibility, packing and intrinsic flexibility?

**Answer.** Remodeling carries additional adjusted rate information in AFP, with evidence
strongest in Afp3; Afp2 is unresolved on variable sites alone. No association was detected in
diffocin. The analysis
used a written specification before rate estimation (`FROZEN_ANALYSIS.md`), which remains
unchanged. The manuscript review discloses a covariate-state implementation mismatch and
parallel variable-site checks extending the literal Afp2-only written clause.

---

## Result 1 — Contact remodeling is quantifiable and heterogeneous

The frozen metric `R_i = 1 − |C_E ∩ C_C| / |C_E ∪ C_C|` was computed for all five
structural references (435 in-domain sites: R2 pyocin 143, AFP Afp2 84, AFP Afp3 82,
SalCIS 68, diffocin 58). Details and QC in `sheath_pilot/QC_RESULTS.md`. This result does
not depend on any evolutionary dataset and therefore covers all five references, including
the two excluded below.

## Result 2 — Evolutionary datasets could be built for two of four systems

Four systems were evaluated for evolutionary inclusion. Two yielded defensible positional
ortholog sets and two did not (`DATASET_DECISION.md`):

- **Diffocin** — **each of the 142 sampled loci contained exactly one detected sheath gene**
  across 80 genera, though some genomes contained more than one candidate locus (7 of 137
  replicons carried a second). Locus-level single-copy status is not genome-level
  single-copy orthology, and the two should not be conflated. 67 loci anchored on non-sheath
  flanking homology at median collinearity ρ = 1.00; position and sequence best hit agree
  36/36. **PASS.**
- **AFP** — Afp2 and Afp3 are 38.1% identical to each other, so sequence similarity cannot
  separate the tandem genes (best-hit labelling disagrees with position for 15 of 22 genes
  at the Afp2 slot). Positional assignment inside the reference subfamily works: 35 of 42
  loci carry the full three-gene run, 95% position/sequence agreement. **PASS on
  per-subfamily sets.**
- **R2 pyocin** — 221 of 967 species carry ≥2 gated candidates; 39 of 57 discarded paralogs
  fall inside the candidate clade; in 11 of 35 testable pairs the discarded paralog
  reconstructs the reference gene order better than the retained copy. **FAIL.**
- **SalCIS** — of 4,490 blastp hits exactly one is *Salmonella* (the reference itself); the
  gated set is dominated by Bacteroidota at 30–45% identity to a Gammaproteobacterial
  reference; the reference locus is unrecoverable in 84% of sampled loci. **FAIL.**

These failures reflect divergent coexisting copies and ambiguous locus assignment in R2
pyocin, and inadequate reference-locus recovery in SalCIS. Their role is to justify the
dataset rather than to establish a shared evolutionary history or open a second research
question.

**Final datasets.** Diffocin 33 sequences / 58 sites; AFP Afp2 27 / 84; AFP Afp3 31 / 81
(one site below 80% column occupancy). Trees: IQ-TREE with ModelFinder — diffocin
LG+F+I+G4, Afp2 WAG+I+G4, Afp3 WAG+I+G4 — 1000 ultrafast bootstrap and 1000 SH-aLRT
replicates. Rates: Rate4Site empirical Bayes, WAG+Γ, 16 categories, one run per dataset.

## Result 3 — Remodeling is associated with rate in AFP, not detectably in diffocin

Per protein, across sites, all variables z-scored within protein:

`z(r) = b0 + bR·z(R) + bA·z(RSA) + bP·z(WCN) + bF·z(GNM) + e`

| | **Diffocin** | **AFP Afp2** | **AFP Afp3** |
|---|---|---|---|
| sites | 58 | 84 | 81 |
| sequences | 33 | 27 | 31 |
| **β_R** | **+0.035** | **+0.435** | **+0.425** |
| SE | 0.135 | 0.118 | 0.114 |
| OLS 95% CI (independent sites) | [−0.237, +0.306] | [+0.201, +0.669] | [+0.198, +0.651] |
| **block-bootstrap 95% CI** | **[−0.173, +0.244]** | **[+0.229, +0.636]** | **[+0.281, +0.544]** |
| p (OLS, independent sites) | 0.80 | 0.0004 | 0.00036 |
| partial R² of `R_i` | 0.001 | 0.148 | 0.155 |
| model R² | 0.121 | 0.226 | 0.209 |
| β_RSA | +0.134 | +0.211 | +0.195 |
| β_WCN | −0.562 | +0.265 | +0.329 |
| β_GNM | −0.545 | +0.137 | +0.279 |
| VIF: `R_i` | 1.11 | 1.41 | 1.24 |
| VIF: RSA / WCN / GNM | 3.17 / 8.39 / 5.06 | 2.06 / 4.99 / 4.45 | 2.12 / 4.69 / 3.95 |

In both AFP proteins a one-SD increase in remodeling corresponds to roughly **0.43 SD
faster evolution** after adjustment. The partial R² of about 0.15 is **not** the proportion
of total rate variance explained: it is the fraction of the variation *left unexplained by
the covariate-only model* that adding `R_i` accounts for (RSS_without − RSS_full) / RSS_without.
Total variance explained by the full four-predictor model is 0.226 and 0.209. In diffocin no
association was detected and the partial R² is 0.001.
The two AFP estimates are nearly identical (+0.435, +0.425) but are **not independent
replicates** — same loci, same genomes, one particle.

### `R_i` is not strongly linearly predictable from the conventional descriptors

Its variance inflation factor is 1.11–1.41 across all three proteins, while RSA, WCN and
GNM are mutually collinear (WCN up to 8.39). This establishes that remodeling is not well
predicted by a linear combination of the included covariates **in these datasets** — which
is what licenses reading β_R as incremental explanatory value rather than a restatement of
burial or packing. It does **not** establish that remodeling is mechanistically independent
of them; low collinearity in three proteins is a property of these samples and this linear
model, not a statement about the underlying structural biology. Note also that **none of the three conventional descriptors reaches
significance in any of the three proteins** (all p ≥ 0.065). Small site counts and
collinearity limit the precision of these estimates; their coefficients are interpreted
as adjustment rather than separate findings.

## Complete archived sensitivity results (`results/sensitivity.csv`)

The written plan names Afp2 variable-site restriction, alternative diffocin sets and occupancy
relaxation. Parallel Afp3/diffocin variable-site checks are present in the archived results
but extend the literal written variable-site clause. They are disclosed as archived checks.

| analysis | sites | β_R | block-bootstrap CI |
|---|---|---|---|
| Diffocin — set as screened (n = 35) | 58 | +0.071 | [−0.152, +0.302] |
| Diffocin — collinear-frame set (n = 66) | 58 | +0.121 | [−0.118, +0.361] |
| Diffocin — variable sites only | 50 | −0.059 | [−0.339, +0.234] |
| AFP Afp2 — variable sites only | 51 | +0.275 | [−0.023, +0.545] |
| AFP Afp3 — variable sites only | 64 | +0.355 | [+0.189, +0.490] |
| Diffocin — occupancy ≥50% | 58 | +0.035 | [−0.173, +0.244] |
| AFP Afp2 — occupancy ≥50% | 84 | +0.435 | [+0.229, +0.636] |
| AFP Afp3 — occupancy ≥50% | 82 | +0.419 | [+0.277, +0.541] |

Three things follow, and all belong in the Results rather than the Discussion.

1. **No association was detected in diffocin under any tested specification.** This holds
   across the set as screened, the 66-sequence collinear-frame set with its own tree and
   rate run, the variable-site restriction and the relaxed occupancy threshold. Every
   interval spans zero. The intervals are wide enough to contain effects of potentially
   meaningful magnitude (the primary interval reaches +0.244, and the collinear-frame
   interval +0.361), so this is a failure to detect an association, **not** evidence that
   remodeling has no evolutionary relationship in diffocin.
2. **The occupancy threshold does not matter.** Relaxing it to 50% adds no sites for
   diffocin or Afp2 and one for Afp3, where β_R moves from +0.425 to +0.419. No conclusion
   turns on the 80% rule.
3. **The Afp2 evidence weakens when invariant sites are omitted.** Thirty-three of Afp2's
   84 structural sites carry no substitution in a 27-sequence alignment.
   Restricted to its 51 variable sites, β_R falls from +0.435 to +0.275 with a CI including
   zero ([−0.023, +0.545]). **Afp2 is unresolved on variable sites alone.** The restriction
   changes the population of sites represented by the estimate; it does not establish that
   the primary association is artifactual. Afp3, with 64 variable sites, retains +0.355
   [+0.189, +0.490].

The defensible statement is therefore narrower than the primary table suggests: the
association is supported in **AFP Afp3**, directionally consistent in Afp2 but not
resolved there on substitution-bearing sites alone, and not detected in diffocin.

## Interpretation, bounded

Within the two systems examined, residue-level contact remodeling during contraction carries
information about site-specific evolutionary rate beyond the included frozen monomer
monomer accessibility, packing and flexibility — in AFP, most clearly in Afp3. The same test in diffocin
detects no association.

**Why the contrast might arise.** Diffocin's dataset has the better taxonomic breadth
(16.71 effective genera against AFP's 5.4) and more substitution signal (86% of structural
sites variable against 61–79%). These summaries do not suggest a simple shortage of
sequence diversity, but do not rule out sampling or rate-estimation limitations. Also,
statistical power here does not follow from sequence diversity alone, and two other
quantities differ in the direction of lower power for diffocin:

- **Less predictor variation.** Diffocin is the least-remodeled of the five references
  (median `R_i` 0.236, against 0.50 for both AFP proteins and 0.80 for R2 pyocin). At the
  fitted sites its `R_i` spread is narrower (SD 0.305 vs 0.383 and 0.393) and 19 of 58 sites
  sit at `R_i` = 0.
- **Noisier rate estimates relative to the signal.** The median Rate4Site posterior SD is
  42% of the between-site rate SD in diffocin, against 31% in Afp2 and 27% in Afp3.

A biological difference between the systems remains equally consistent with the data — AFP is
a plasmid-borne antifeeding prophage particle with three tandem sheath genes, diffocin a
chromosomal prophage-derived bacteriocin with one. The analysis cannot distinguish a reduced
dynamic range of remodeling, noisier rate estimation, and a genuine contextual difference,
and with three proteins from two systems it should not try.

## What cannot be claimed

- No general statement about contractile injection systems. Two systems, one contributing
  two non-independent proteins.
- No cross-system generalisation from the AFP/diffocin contrast; the contrast is reported,
  not explained.
- Nothing about R2 pyocin or SalCIS evolution — their `R_i` values stand, their rates were
  never estimated.
- The conventional descriptors' coefficients are not findings at this sample size.

## Limitations

- **Covariate-state reporting.** The retained structural producer averages both monomer
  states, whereas the written specification says extended-state only. Separate state values
  are not archived. The manuscript reports this mismatch and uses the saved predictors
  without refitting or retrospectively choosing states.

- **Dependence between sites is addressed along the sequence, not in three dimensions.**
  The block bootstrap resamples **contiguous 10-residue windows of the reference numbering**
  (diffocin 17 blocks, Afp2 26, Afp3 23; 1–10 analysed sites per block, median 2–3). It
  allows local sequence dependence, but does not necessarily account for dependence between
  residues far apart in sequence but adjacent in the folded sheath. No explicit spatial or
  3D-neighbourhood model was fitted.
- **Bootstrap scaling.** Predictors and the response are re-standardised **within each
  bootstrap replicate**, so the interval is for the standardised slope re-estimated on each
  resample and absorbs resample-to-resample changes in predictor spread; it is not an
  interval for the slope in fixed original units.
- **Rate uncertainty is not propagated, and β_R should not be read as a precisely estimated
  biological effect.** The regression treats estimated rates as observed outcomes, while
  Rate4Site reports a substantial posterior SD per site (median 0.401 for diffocin, 0.218
  and 0.214 for Afp2 and Afp3 — 27–42% of the between-site rate SD). Inverse-variance
  weighting is not an obvious remedy: Rate4Site's posterior SD depends on conservation and
  hence on the rate itself, so naive weighting could bias the association it is meant to
  stabilise. No measurement-error model was prespecified, and none was introduced afterwards
  to strengthen the result.
- **Substitution model.** Rate4Site 2.01 does not implement LG; rates use WAG+Γ in all three
  primary datasets, while ModelFinder's overall choice for diffocin was LG+F+I+G4. WAG was
  the best-fitting *available* matrix in each of the three primary datasets. The two diffocin
  sensitivity sets selected **Q.PFAM+F+G4** — a model family Rate4Site also does not
  implement — and were run under the same frozen WAG+Γ rule, so their topologies come from a
  model the rate estimation does not share. Both return the same null as the primary set, so
  nothing in the conclusion rests on that mismatch; it is recorded because the primary
  diffocin dataset and its sensitivity sets were not fitted under the same selected model.
- **Taxonomic breadth of AFP** is narrow: 5.4 effective genera in a single enterobacterial
  radiation.
- **Afp2's invariant sites** limit what that protein can contribute, as above.
- **Prior exposure.** An entropy–`R_i` Spearman correlation (ρ = 0.32/0.33) was computed at
  the 143 R2 pyocin structural sites as an alignment diagnostic before this analysis was
  specified. It was not covariate-adjusted and was not a rate estimate, and R2 pyocin is in
  none of the datasets tested here — but the hypothesis was not unexamined when the
  specification was written, and this analysis is prespecified only from that point forward.

## Files

Specification `FROZEN_ANALYSIS.md` · code `hypothesis_test.py` · primary results
`results/hypothesis_test_primary.{csv,json}` · sensitivity `results/sensitivity.{csv,json}`
· per-site tables `results/sitetable_*.csv` · figure `results/hypothesis_result.png` ·
trees `phylo/` · rates `rates/` (with the exact command line per dataset in `rates/*.cmd`).
