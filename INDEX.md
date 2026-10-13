# Index — where every result lives

Project root: `/Users/anthony/Documents/sheath_evolution`. Written 2026-10-13.
**Read `RESULTS.md` first** for the findings and their limitations. Manuscript and submission
files listed below belong to the local authoring project and are excluded from GitHub.
The versioned files provide the scientific inputs, outputs, code and decision evidence.

---

## 1. The answer to the question

| | |
|---|---|
| **`MANUSCRIPT.md`** | Scientific draft with focused Introduction, three Results subsections and 23 verified references. |
| **`manuscript.tex`** | Author/contact/declarations, native equations and three main tables/figures; includes `supplementary.tex`. |
| **`manuscript_latex.zip`** | Editable LaTeX package with selected figures, supplementary data and build instructions. |
| **`SUPPLEMENTARY.md` / `supplementary.tex`** | Detailed methods/results, Tables S1–S5 and Figures S1–S5. |
| `supplementary_data/` | Full-precision tables, residue mappings, commands, selected archived outputs and source checksums. |
| `reproducibility_archive.zip` | Local code/data archive prepared for deposition, with complete file/hash inventory; no public URL/DOI assigned. |
| `REPRODUCIBILITY_AUDIT.md` / `structural_provenance_audit.json` | Ordered structural selection procedure, missing input/key audit and downstream replay limits. |
| `REVIEW_RESPONSE.md` | Item-by-item response to the supplied manuscript review. |
| `title_page.tex` / `manuscript_for_review.tex` | Separate identifying page and optional review copy; the full author version remains `manuscript.tex`. |
| `SUBMISSION_CHECKLIST.md` / `COVER_LETTER_DRAFT.md` | Journal preparation and remaining publication/build/deposition facts. |
| `validate_manuscript.py` / `build_reproducibility_archive.py` | Source/numerical reporting checks and local archive preparation; neither fits models. |
| `LATEX_README.md` | Compilation instructions; manuscript PDF build remains unverified. |
| **`MANUSCRIPT_REVIEW.md`** | Scientific/citation review, implementation reporting discrepancies and remaining submission items. |
| **`RESULTS.md`** | Saved scientific results, complete archived sensitivities and bounded interpretation. |
| `figures/Figure1.{png,pdf}` | Remodeling distributions across 435 structural sites. |
| `figures/Figure2.{png,pdf}` | Four-system evolutionary inclusion decisions. |
| `figures/Figure3.{png,pdf}` | Three added-variable panels: Afp2, Afp3 and diffocin, with exact saved slopes. |
| `figures/FigureS1*`–`FigureS5*` | Structural QC, orthology evidence and descriptive rate uncertainty. |
| `prepare_supplementary.py` | Prepare presentation tables and copy selected provenance; no scientific fits. |
| `write_manuscript_latex.py` | Convert reviewed Markdown, renumber citations and create the ZIP/extracted package. |
| `make_figures.py` | Render archived values and explicitly hypothetical schematics; no analysis reruns. |
| `results/hypothesis_test_primary.csv` / `.json` | Saved primary coefficients, uncertainty, variance and predictor statistics. |
| `results/sensitivity.csv` / `.json` | All eight archived sensitivity checks, including the two parallel variable-site restrictions beyond the literal written clause. |

**Headline numbers** (adjusted remodeling effect, SD of rate per SD of `R_i`, sequence
block-bootstrap 95% CI): diffocin **+0.035 [−0.173, +0.244]**, AFP Afp2
**+0.435 [+0.229, +0.636]**, AFP Afp3 **+0.425 [+0.281, +0.544]**. On variable sites only,
Afp2 falls to +0.275 [−0.023, +0.545] and Afp3 holds at +0.355 [+0.189, +0.490].

## 2. How the test was specified and run

| | |
|---|---|
| **`FROZEN_ANALYSIS.md`** | The original written specification, preserved unchanged; state-definition and sensitivity-scope discrepancies are disclosed in the manuscript and review |
| `hypothesis_test.py` | The retained analysis code; includes the rate-direction check and the bootstrap over contiguous 10-residue sequence windows |
| `phylo/*.iqtree`, `phylo/*.treefile` | IQ-TREE model selection and ML trees (5 datasets). `phylo/*.r4s.nwk` = the same trees with support labels stripped, which Rate4Site requires |
| `rates/*.res` | Rate4Site output, per-site rates against the reference numbering |
| `rates/*.cmd` | The exact command line used for each dataset — reproduce from these |
| `rates/*.unnorm`, `rates/*.log` | Un-normalised rates (used for the direction check) and run logs |
| `results/sitetable_primary_*.csv` | The fitted per-site table per protein: rate, rate SD, occupancy, `R_i`, RSA, WCN, GNM |
| `results/partial_residuals.csv` | Covariate-adjusted residuals behind the figure's three panels |
| `tools/rate4site`, `tools/src/` | The Rate4Site binary built here, and the patched source |

## 3. Which datasets were used, and why those

| | |
|---|---|
| **`DATASET_DECISION.md`** | The definitive dataset list, per-system reasoning, what the manuscript can and cannot claim |
| `orthology/dataset_decision.png` | Figure: taxonomic independence vs informative-site fraction, by inclusion status |
| `orthology/dataset_decision.csv` | The same table as data |
| `orthology/diffocin_primary.afa` | **Primary** diffocin alignment (33 seq) — reference row `REF` |
| `orthology/afp_sf1_afp2_like.afa` | **Primary** AFP Afp2 alignment (27 seq) |
| `orthology/afp_sf1_afp3_like.afa` | **Primary** AFP Afp3 alignment (31 seq) |
| `orthology/diffocin_frozen_offset0.afa` | Sensitivity: diffocin as screened (35 seq) |
| `orthology/diffocin_collinear_frame.afa` | Sensitivity: diffocin collinear frame (66 seq) |
| `orthology/primary_alignment_verification.csv` | Proof each alignment's `REF` row equals the PDB SEQRES exactly |

## 4. Per-system feasibility and orthology verdicts

| system | verdict | document |
|---|---|---|
| Diffocin | PASS | `DIFFOCIN_FEASIBILITY.md` |
| AFP | PASS on per-subfamily sets | `AFP_FEASIBILITY.md`, `AFP_CHECK23_VERDICT.md`, `AFP_CHECK1_LOCUS.md` |
| SalCIS | FAIL | `SALCIS_DECISION.md` |
| R2 pyocin | FAIL | `R2_FINAL_DECISION.md`, `ORTHOLOGY_DECISION_r2pyocin.md`, `R2_EVOLUTIONARY_FEASIBILITY.md` |
| retrieval validation | — | `RETRIEVAL_VALIDATION_r2pyocin.md` |

Supporting data: `orthology/` (locus tracks, sheath runs, feasibility tables, trees),
`retrieval_uniprot/` (candidate pools + provenance per family), `filtered/` (gated hits),
`synteny/`.

## 5. Structural analysis (Result 1) and the superseded pilot

| | |
|---|---|
| `sheath_pilot/qc_output/R_*.csv` | The frozen `R_i` values and the `rsa`/`wcn`/`gnm` covariates, per family — the structural input to everything above |
| `sheath_pilot/QC_RESULTS.md` | Structural QC for the continuous metric |
| `sheath_pilot/README.md` | The historical categorical pilot: 0 of 3 families pass. **Supplementary / methodological history only** — read its warnings before reusing any number |
| `sheath_pilot/output/archive_symmetry_inclusive_2026-10-10/` | The superseded symmetry-inclusive run, preserved unmodified |

## 6. Project state

`README.md` — current state at the top, then layout, environment, the five structural
references, and the methodological commitments carried through the project.

## 7. Writing and review status

The Introduction, references and integrated scientific review are complete in
`MANUSCRIPT.md`. `MANUSCRIPT_REVIEW.md` records the author and submission information still
needed, including a public data/code location and final TeX compilation. The expanded main text is approximately 4,950 words excluding tables.

The computational phase is closed. No further proteins, regression variants or evolutionary
measurements are to be added.
