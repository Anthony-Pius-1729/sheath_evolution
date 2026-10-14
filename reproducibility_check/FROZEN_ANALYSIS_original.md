# Frozen analysis specification

**Written and frozen 2026-10-13, before any evolutionary rate was estimated.**
Nothing below may be changed on the basis of a fitted result. Deviations, if any become
unavoidable, are to be appended to §8 with a reason and a date — never silently applied.

**Question.** Does contact-network remodeling during contraction explain site-specific
evolutionary-rate variation beyond solvent accessibility, packing and intrinsic flexibility?

---

## 1. Datasets (fixed)

| protein | alignment | n sequences | structural sites |
|---|---|---|---|
| Diffocin (8V3X entity 1) | `orthology/diffocin_primary.afa` | 33 | 58 |
| AFP Afp2 (6RBN entity 2) | `orthology/afp_sf1_afp2_like.afa` | 27 | 84 |
| AFP Afp3 (6RBN entity 3) | `orthology/afp_sf1_afp3_like.afa` | 31 | 82 |

Each alignment contains its own structural reference as the row `REF`, whose ungapped
sequence equals the PDB SEQRES exactly — verified in
`orthology/primary_alignment_verification.csv`. Reference residue numbering is therefore
mapped to alignment columns without ambiguity, and `R_i` values join on `seqres_pos`.

The diffocin primary set **excludes both *Peptacetobacter hiranonis* copies**
(WP_006439294.1, WP_006440377.1, replicon NZ_CP036523.1), which received two mutually
exclusive positional assignments. The exclusion is justified by unresolved positional
orthology and was decided before any rate was computed. R2 pyocin and SalCIS are excluded
entirely (`R2_FINAL_DECISION.md`, `SALCIS_DECISION.md`).

## 2. Phylogenies (fixed procedure)

Final trees are inferred per dataset with IQ-TREE: ModelFinder (`-m MFP`) over the standard
amino-acid substitution models with rate heterogeneity, full ML branch-length optimisation,
1000 ultrafast bootstrap replicates and 1000 SH-aLRT replicates. The selected model and the
ML tree are recorded in `phylo/<dataset>.{iqtree,treefile,log}`.

The FastTree topologies used during orthology screening are **not** used for rate
estimation.

Selected models (ModelFinder, BIC): diffocin **LG+F+I+G4**, Afp2 **WAG+I+G4**, Afp3
**WAG+I+G4**.

## 3. Rate estimation (fixed)

Rate4Site (Pupko et al. 2002; Mayrose et al. 2004), **empirical Bayes** rates (`-ib`,
16 rate categories), one run per dataset, on each primary alignment with its IQ-TREE
topology supplied as a fixed tree (`-t`) and the reference row named as the query sequence
(`-a REF`) so rates are reported per reference residue. The `.res` output, the log and the
exact command line are archived under `rates/`.

Three implementation facts, recorded because they constrain the analysis:

- **Binary provenance.** No Rate4Site package exists for this platform in conda or pip. It
  was built from the ConSurf-bundled source tree (`rate4site.3.2.source_fast/sourceMar09`,
  Rostlab/ConSurf), which fails to compile under libc++; three pre-C++11 stream idioms were
  patched (`*_errorOut != cerr` → address comparison in `errorMsg.cpp`; `stream == NULL` →
  `!stream.is_open()` in `someUtil.cpp` and `rate4siteOptions.cpp`). The patches touch only
  error-reporting and file-open checks, no numerical code. **The resulting binary
  self-reports version 2.01**, not 3.2; the source directory label and the binary disagree,
  and 2.01 is what is being run. Source and patches are kept under `tools/src/`.
- **Substitution model.** Rate4Site 2.01 implements JTT, WAG, Dayhoff, REV, cpREV and JC
  only — not LG. The model is therefore fixed to **WAG (`-Mw`)** for all three datasets,
  which is the best-fitting *available* matrix by ModelFinder BIC in every one of them
  (diffocin: WAG+F+I+G4 at 32067.46, the best supported alternative to the overall choice
  LG+F+I+G4; Afp2: WAG+I+G4, also the overall choice; Afp3: WAG+F+I+G4, the best supported,
  overall choice WAG+I+G4). The IQ-TREE run therefore supplies the topology and the model
  ranking, and the rate estimation runs under WAG+Γ throughout.
- **Branch lengths.** The IQ-TREE ML topology is fixed, and branch lengths are
  **re-optimised by Rate4Site under its own gamma model** (the default, and the ConSurf
  convention), so that lengths and rates are on one model scale rather than importing
  lengths fitted under a different matrix. The IQ-TREE ML branch lengths are retained
  unchanged in `phylo/*.treefile` for the record.

**Rate direction and normalisation.** Rate4Site reports relative rates with the alignment
mean at 1, where **larger values mean faster evolution**; the sign convention is verified
against the per-site conservation scores in the same output before use. No transformation
is applied other than standardisation. Rates are then **z-scored within each protein**, so
β coefficients are in units of within-protein standard deviations.

Structural sites whose alignment column falls below 80% occupancy are dropped before
fitting; the count dropped is reported per protein.

Occupancy recomputed on the three **final** primary alignments (not on the sets as
screened): diffocin 58/58 sites, minimum occupancy 0.970 on 33 sequences; Afp2 84/84,
minimum 1.000 on 27; Afp3 81/82, minimum 0.742 on 31 — one Afp3 site dropped. These are
the counts the fitted models use (`results/hypothesis_test_primary.json`,
fields `sites_dropped_occupancy` / `sites_fitted`).

## 4. The model (fixed)

Per protein, ordinary least squares across sites:

```
z(r_i) = b0 + bR·z(R_i) + bA·z(RSA_i) + bP·z(WCN_i) + bF·z(GNM_i) + e_i
```

- `R_i` — frozen contact-remodeling score, `sheath_pilot/qc_output/R_*.csv`.
- `RSA_i`, `WCN_i`, `GNM_i` — the **already frozen monomer-based** descriptors in the same
  tables (`rsa`, `wcn`, `gnm` columns). Monomer, extended state, as frozen. No selection
  among extended-state, contracted-state or averaged variants.
- All predictors and the response are z-scored **within protein**.

**The primary scientific quantity is `bR`** — three protein-specific estimates, reported
separately. No pooled or meta-analytic estimate is the primary result.

Fixed exclusions: no interaction terms, no nonlinear terms, no variable selection, no
additional covariates. Four predictors plus intercept, on 58–84 sites.

## 5. Uncertainty (fixed)

Sites in a protein structure are spatially dependent, so independent-site OLS intervals
understate uncertainty. Both are reported, and the second is the one interpreted:

1. OLS 95% CI, labelled as assuming independent sites.
2. **Spatial block bootstrap**, 10,000 resamples, blocks defined by contiguous runs of
   reference residue numbering of fixed length 10, resampled with replacement to the
   original site count. Block length and replicate count are fixed here and are not tuned.

Reported per protein: `bR`, its OLS CI, its block-bootstrap CI, the partial R² of `R_i`,
and the model R². Variance inflation factors for the four predictors are reported as
descriptive collinearity information.

## 6. Prespecified sensitivity analyses (the complete list)

1. **Diffocin as screened** — the 35-sequence frozen set including both *P. hiranonis*
   copies.
2. **Diffocin collinear-frame set** — the 66-sequence set (`diffocin_collinear_frame.afa`),
   requiring a new IQ-TREE run and a new Rate4Site run under the same specification.
3. **Afp2 invariant sites** — the fit with and without its invariant structural sites, which
   carry no substitution information.
4. **Occupancy threshold** — structural sites at ≥50% instead of ≥80% column occupancy.

Nothing else. No additional sensitivity analysis may be introduced after seeing `bR`.

Sensitivity runs follow this specification unchanged, including the fixed WAG+Γ rate model.
Recorded for the record: ModelFinder selected **Q.PFAM+F+G4** for both diffocin sensitivity
alignments (as screened, and collinear frame), a family Rate4Site 2.01 does not implement;
their topologies therefore come from a model the rate estimation does not share, exactly as
for diffocin's primary LG+F+I+G4 topology.

## 7. Interpretation rules (fixed in advance)

- The result is reported in whichever direction it falls, including a null.
- Any conclusion is scoped to **three proteins from two systems** (AFP and diffocin). AFP
  Afp2 and Afp3 come from the same loci in the same genomes and are **not independent
  system-level replicates**.
- Concordance across the three proteins is described as consistency, not as a tested
  cross-system generalisation.
- No claim is made about contractile injection systems in general.

## 8. Disclosure, and deviations

**Prior exposure.** A Spearman correlation between alignment-column entropy and frozen `R_i`
at the 143 R2 pyocin structural sites was computed earlier as an alignment diagnostic:
ρ = 0.32 (narrow set), 0.33 (wide set). It is not covariate-adjusted and is not a rate
estimate, but it is outcome-adjacent and it has been seen. The analysis specified here is
prespecified **from this point forward**; it is not independent of prior sequence-variability
observations and must not be described as such. R2 pyocin, the protein in which that
correlation was computed, is not in any dataset above.

**Deviations from this specification:** none at the time of freezing.
