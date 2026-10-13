# Contact-network remodeling and bacterial sheath protein evolution

This repository contains the code and retained analysis inputs and outputs for
a study of how intermolecular contacts change during contraction of bacterial injection-system
sheaths. The central question is whether residue-level contact remodeling is associated with
site-specific evolutionary rate after adjustment for monomer solvent accessibility, packing
and flexibility averaged across the extended and contracted states.

The computational analysis is closed. The repository preserves the reported estimates and
separates the primary analysis from historical feasibility checks and manuscript preparation.
Manuscript sources and submission documents are kept locally and excluded from this repository.
The repository supports inspection of the evidence and reproduction from retained downstream inputs;
missing atomic-coordinate provenance limits exact reconstruction of the original structural
calculations.

**Study snapshot:** 13 October 2026

## Main finding and scope

Contact remodeling shows a positive adjusted linear association with evolutionary rate in
AFP Afp2 and Afp3. Evidence is most consistent in Afp3: its sequence-block bootstrap interval
remains positive when invariant alignment sites are omitted. Afp2's interval includes zero
under that restriction. An association was not detected in diffocin.

| Primary dataset | Sequences | Fitted sites | Variable sites | Standardized remodeling coefficient | Sequence-block bootstrap 95% interval |
|---|---:|---:|---:|---:|---|
| AFP Afp2 | 27 | 84 | 51 | +0.435 | [+0.229, +0.636] |
| AFP Afp3 | 31 | 81 | 64 | +0.425 | [+0.281, +0.544] |
| Diffocin | 33 | 58 | 50 | +0.035 | [−0.173, +0.244] |

On variable sites, the coefficients and intervals are +0.275 [−0.023, +0.545] for Afp2,
+0.355 [+0.189, +0.490] for Afp3, and −0.059 [−0.339, +0.234] for diffocin.
Full-precision values are retained in [`results/hypothesis_test_primary.json`](results/hypothesis_test_primary.json)
and [`results/sensitivity.json`](results/sensitivity.json).

The three fitted datasets comprise **223 sites from two systems**. Afp2 and Afp3 belong to
one AFP system and are not independent system-level replicates. Detecting an association
in AFP and not detecting one in diffocin does not establish a difference between their
underlying effects. Diffocin's narrower remodeling range and greater relative rate uncertainty
prevent confident attribution of the contrast to biology. These are conditional, in-sample
linear associations, not causal effects or a test of adaptive evolution.

## Read the study

| File | Purpose |
|---|---|
| [`RESULTS.md`](RESULTS.md) | Numerical findings and their interpretation. |
| [`INDEX.md`](INDEX.md) | Map from each result to its supporting files. |
| [`FROZEN_ANALYSIS.md`](FROZEN_ANALYSIS.md) | Original written specification, preserved unchanged. |
| [`DATASET_DECISION.md`](DATASET_DECISION.md) | Evolutionary inclusion/exclusion decisions. |
| [`REPRODUCIBILITY_AUDIT.md`](REPRODUCIBILITY_AUDIT.md) | Input completeness, deterministic structural selection and unresolved provenance. |

The retained publication figures display structural distributions, dataset selection and
adjusted associations. Figure 3 shows added-variable plots for Afp2, Afp3 and diffocin; its
lines use the saved full-model slopes and are checked against the stored residuals. Full
coefficient and sensitivity tables remain under `supplementary_data/`.

## Structural datasets and remodeling

The structural comparison covers **435 eligible sites in five proteins from four systems**.
R2 pyocin and SalCIS contribute structural results but were excluded from evolutionary-rate
analysis because the available evidence did not support reliable reference-lineage assignment.

| System | Protein | Extended PDB | Contracted PDB | Structurally eligible sites | Rate analysis |
|---|---|---|---|---:|---|
| R2 pyocin | Sheath | 3J9Q | 3J9R | 143 | Excluded: reference lineage ambiguous |
| AFP | Afp2 | 6RBN | 6RC8 | 84 | Retained |
| AFP | Afp3 | 6RBN | 6RC8 | 82 | Retained; one site fails primary occupancy |
| SalCIS | Sheath | 9R9H | 9R9I | 68 | Excluded: reference locus poorly recovered |
| Diffocin | CD1363 | 8V3X | 8V3Y | 58 | Retained |

For residue $i$, remodeling compares its extended and contracted intermolecular contact sets:

```math
R_i = 1 - \frac{\lvert C_i^E \cap C_i^C\rvert}{\lvert C_i^E \cup C_i^C\rvert}.
```

Contacts use minimum heavy-atom distance strictly below 5 Å and a partner identity containing
protein family, relative ring, relative azimuth and canonical residue position. Eligibility
requires at least three identities in the contact union. An index of zero denotes an unchanged
contact set; one denotes no shared identities. The index does not measure force, energy,
displacement or functional importance.

AFP uses deposited biological-assembly operators in both states. Operator expansion is
distinguished from guessed ring completion; none of the five retained references required
the latter. The earlier categorical interface-persistence pilot is historical development
work, not the primary continuous-remodeling analysis. Its files remain under
[`sheath_pilot/`](sheath_pilot/README.md).

## Evolutionary and statistical approach

Candidate homologs were retrieved from bacterial UniProtKB release 2026_03 using reference
Pfam signatures. Sequence similarity and species dereplication were insufficient to establish
orthology for duplicated sheath genes. Non-sheath flanking homology, conserved gene order,
orientation and positional assignment were evaluated before the primary rate models.

Final alignments used FAMSA through pyfamsa 0.7.0. IQ-TREE 3.1.4 selected models with
ModelFinder and supplied the final topologies. Rate4Site self-reported version 2.01 and used
empirical Bayes estimation under WAG+Γ with 16 rate categories. Larger normalized scores
indicate faster relative evolution. Exact retained commands are under [`rates/`](rates/)
and in [`supplementary_data/archived_commands.txt`](supplementary_data/archived_commands.txt).

Each protein has a separate standardized linear model:

```math
\begin{aligned}
z(r_i) = \beta_0 &+ \beta_R z(R_i) + \beta_A z(\mathrm{RSA}_i) \\
&+ \beta_P z(\mathrm{WCN}_i) + \beta_F z(\mathrm{GNM}_i) + \epsilon_i.
\end{aligned}
```

RSA, WCN and GNM are the **archived monomer descriptors averaged across available states**.
The remodeling coefficient describes rate SD per remodeling SD conditional on those inputs.
Partial $R^2$ measures the fraction of variation left unexplained by the covariate-only model
that adding remodeling accounts for.

Uncertainty uses 10,000 resamples of nonempty, ten-reference-position sequence windows,
with variable eligible-site counts, truncation to the original sample size and
restandardization within each replicate. These intervals retain some local sequence dependence.
Their coverage under remaining three-dimensional dependence was not assessed, and rate
posterior uncertainty was not propagated into the regression.

## Repository organization

```text
sheath_pilot/          Structural producer, contact metric, covariates and retained QC
retrieval_uniprot/     Recorded candidate sequence pools and retrieval provenance
filtered/             Similarity/architecture screening outputs
synteny/              Genomic-neighbourhood evidence
evolution_prelim/     Preliminary feasibility outputs
orthology/            Positional assignments, diagnostics and final alignments
phylo/                Final IQ-TREE outputs and Rate4Site-compatible topologies
rates/                Rate estimates, uncertainty, logs and exact commands
results/              Fitted site tables, model JSON/CSV and saved partial residuals
figures/              Current publication figures as PNG/PDF
supplementary_data/   Full-precision presentation tables and selected source copies
```

Root Python scripts implement retrieval, filtering, locus analysis, phylogenetic assessment,
regression and manuscript preparation. Their roles are indexed in [`INDEX.md`](INDEX.md).
Historical alternatives under `orthology/` and preliminary directories are retained as
selection evidence. Their counts and diversity summaries should not be substituted for
those of the final primary datasets.

The final primary alignments are:

- [`orthology/afp_sf1_afp2_like.afa`](orthology/afp_sf1_afp2_like.afa)
- [`orthology/afp_sf1_afp3_like.afa`](orthology/afp_sf1_afp3_like.afa)
- [`orthology/diffocin_primary.afa`](orthology/diffocin_primary.afa)

## Inspect and regenerate presentation

Run commands from the repository root. The following commands use retained outputs and do
not recompute structural contacts, infer phylogenies, estimate rates or fit regressions:

```sh
python3 prepare_supplementary.py
python make_figures.py
```

`prepare_supplementary.py` uses the Python standard library. Figure generation requires
NumPy, pandas and Matplotlib. Preparation regenerates table CSVs, selected source copies and
a local supplementary draft. Independent source manifests remain part of the repository.
The LaTeX conversion and manuscript-validation scripts additionally require local manuscript
files, which are intentionally excluded from Git.

These presentation commands rewrite generated supplementary tables, a local draft and figure files.
They should be run in a working copy when preserving a particular manuscript snapshot.

To prepare a local code/data bundle with a SHA-256 inventory:

```sh
python3 build_reproducibility_archive.py
```

This produces `reproducibility_archive.zip` and its metadata report. The local bundle may
include manuscript files when they are present in the working directory; the ZIP is excluded
from Git. Generated ZIP packages,
extracted duplicate bundles, caches, downloaded article XML and local third-party tool builds
are excluded from Git. The retained source and study data are versioned directly.

### Local manuscript preparation

Manuscript drafts, LaTeX sources, title pages, submission documents and generated manuscript
ZIPs stay in the local working directory and are excluded by `.gitignore`. They are not
available in a fresh clone. The authoring scripts are retained, but their document inputs
must be supplied locally before running:

```sh
python3 write_manuscript_latex.py
python3 validate_manuscript.py
```

The validation script checks manuscript syntax, citations, saved numerical estimates and
package integrity when its local document inputs exist. Manuscript editing and compilation
are separate from inspection of the versioned code and scientific outputs.

## Reproduction boundaries

The retained downstream inputs and outputs are available for inspection. An exact end-to-end
reproduction claim would exceed what has been preserved:

- **Covariate definition:** the written specification names extended-state covariates, while
  the retained producer writes means across states. The reported models condition on the
  actual means. An extended-state-only adjusted association has not been estimated.
- **Structural provenance:** original mmCIF bytes/checksums, selected contact-reference keys,
  representative monomer-chain keys, candidate tie scores and separate descriptor arrays
  were not retained. The selection algorithm is documented, but exact historical selected
  subunits cannot be verified from the saved outputs alone.
- **Sensitivity scope:** the written variable-site clause explicitly names Afp2. Parallel
  Afp3 and diffocin restrictions are archived checks beyond that literal clause; their timing
  is not invented. All eight saved sensitivity results are reported.
- **Dependence and estimated outcomes:** sequence-block intervals do not fully represent
  spatial dependence or uncertainty in rate point estimates, trees, structures or orthology.
- **Environment:** [`requirements.txt`](requirements.txt) records packages from the original
  environment but is not a verified complete lockfile. It omits the statsmodels version used
  for regression. External tools include BLAST+, FastTree, IQ-TREE and the patched Rate4Site
  build; recorded versions/build notes are in the specification and logs. Local compiled
  binaries and bundled third-party source are not distributed through this repository.

`hypothesis_test.py` refits the models and rewrites analysis outputs; structural and retrieval
scripts likewise initiate substantive processing or external requests. They are retained for
scientific inspection and downstream replay, rather than included in the presentation commands
above. No analyses were rerun to prepare this repository or README.

## Citation, availability and declarations

The manuscript title is **Contact-Network Remodeling During Bacterial Sheath Contraction Is
Associated with Site-Specific Evolutionary Rates in an AFP System**.
The manuscript is a local research draft and is not distributed in this repository; no
publication DOI or archival release DOI has been assigned. Cite the manuscript and the specific repository commit or archived release used,
and retain the primary experimental/data citations associated with the analysis.

A source repository and a persistent scholarly data deposit serve different purposes.
Public visibility and an archived release/DOI must be established before replacing the
manuscript's pending public-deposition statement. Local ZIP files are not persistent deposits.

The author reports no funding and no competing interests. The study reanalyzes public protein
structure and sequence data and involves no new human-participant or animal experiments.
Editorial and presentation assistance is disclosed in the manuscript acknowledgments.
A software license has not yet been assigned; third-party resources retain their own terms.
