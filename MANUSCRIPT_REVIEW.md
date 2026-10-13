# Integrated manuscript review

Reviewed 2026-10-13. Scope: scientific reporting, writing, primary references and presentation.
The computational analysis remains closed: no contact calculations, rate estimates, tree
inference, model fitting or additional experiments were run during this revision.

## Focused review revision

The supplied review prompted a second methods/provenance pass. The fitted estimand is now
explicitly conditional on monomer descriptors averaged across states in the abstract,
Methods, Results, Figure 3 caption, Discussion and Conclusion. The intended extended-state-only
question remains untested. The justification explains what the averages measure without
claiming equivalence to the written state definition.

Sequence-block interval coverage under spatial dependence is unassessed. Structures,
orthology assignments, trees and rate point estimates are held fixed by this bootstrap.
All interval labels/interpretations preserve that conditional scope. Original point estimates
and bootstrap endpoints remain unchanged.

The provenance audit documents model/chain/operator order, interior-ring preference,
strict <90 Å centroid scoring, self-centroid inclusion and first-in-order ties. It also
separates representative monomer chains from contact-reference subunits. Original mmCIF
bytes, selected keys and candidate scores were not retained; exact historical reference
subunits therefore cannot be verified. `REPRODUCIBILITY_AUDIT.md` and its machine-readable
record retain that limitation rather than substituting example chain IDs.

The three main figures/tables remain. Figure 1's text connects partner exchange to the
connected sheath and does not imply spatial mapping or forces. Figure 3's footer names
the averaged monomer adjustment. `REVIEW_RESPONSE.md` tracks each review item.

A separate local deposition archive is prepared with code, recorded inputs, outputs and
SHA-256 inventory. It has not been pushed or assigned a public URL/DOI. Separate title-page
and review-copy sources, a cover-letter draft and a submission checklist accompany the
full author LaTeX version. The current numerical environment is not retrospectively claimed
to be the original environment; the captured requirements omit the fitting statsmodels version.

## Manuscript structure and declarations

The revised manuscript follows the supplied checklist: five Introduction paragraphs, three
Results subsections, three main figures and three main tables. Its main narrative is about
4,900 words excluding abstract, tables, captions and references; the abstract is about 200
words. Twenty-three references appear in citation order. Native LaTeX equations include the
boxed standardized regression and the partial R² definition. Supplementary Methods/Results
provide Tables S1–S5 and Figures S1–S5, with full-precision CSVs, residue mappings, retained
commands, selected outputs and source checksums. The dedicated archive additionally retains recorded
sequence pools and intermediate decision evidence.

The author is Anthony Yua Ior, Istanbul Technical University, Maslak 34867, Istanbul, Turkiye.
The user confirmed ior23@itu.edu.tr, no funding, no competing interests, and all stated roles:
study design, software/computational analysis, interpretation and manuscript preparation.
Ethics/consent are not applicable to the reported reanalysis of public protein data.

## Scientific interpretation preserved

Afp2 and Afp3 have positive primary adjusted slopes. The abstract, main Table 3, Discussion
and Conclusion retain Afp2's weaker variable-site evidence and Afp3's positive variable-site
interval. Diffocin is reported as an association not detected; narrower remodeling and
greater relative rate uncertainty prevent confident attribution of the AFP–diffocin contrast
to biology. No formal between-system effect test was performed. The two AFP proteins are
not independent system-level replicates, and 223 fitted sites do not provide 223 system
replicates. Partial R² measures variance left unexplained by the reduced model, rather than
fraction of total variance. Results remain in-sample linear associations, not causal or
adaptive-evolution tests.

The bootstrap uses ten-reference-position nonempty sequence windows, variable eligible
row counts, sampling to the original count with final-block truncation, and restandardization
in all 10,000 replicates. Its intervals allow local sequence dependence but do not explicitly
represent sequence-distant three-dimensional dependence. Rate posterior uncertainty remains
unpropagated. Invariant sites receive empirical Bayes estimates, not one common fixed floor.

## Reporting discrepancies found during the implementation audit

| Issue | Evidence | Reporting action |
|---|---|---|
| Covariate state differs from written specification | `FROZEN_ANALYSIS.md` specifies extended-state monomers; `sheath_pilot/structural_qc.py` computes both states and writes `nanmean` across them. Separate state values are not retained. | Methods and limitations disclose the mismatch; no states were selected and no regression was rerun. |
| Engineering exclusion belongs to an earlier producer | The 8 Å shell is in `sheath_pilot/pilot.py`, not the continuous structural producer. | No engineering shell is claimed for the continuous analysis. |
| AFP operator expansion was conflated with ring completion | Stored QC reports 24/48 operator-expanded subunits; generated-contact flags identify guessed ring-completion keys. | Deposited operators and guessed completion are distinguished. No retained reference needs ring completion. |
| Selected contact subunit key is missing | The algorithm selects interior/max-neighbour subunits, but their keys are not written to retained tables. | Main chains are explicitly examples of deposited family members; the unlogged key is disclosed. |
| Variable-site plan is narrower than saved checks | The written clause explicitly names Afp2; saved results include Afp3 and diffocin restrictions. | Parallel restrictions are identified as archived checks beyond that literal clause; their timing is not invented. |
| A saved JSON invariance field is misleading | `invariant_sites_in_fit` counts zero posterior SD; site-table `n_aa` is occupancy, not diversity. | Counts of 51/64/50 variable primary sites were verified directly against cached alignments; original exports remain unchanged. |
| Screening stages differ from final datasets | Early diffocin locus file has 83/40 total/anchored, expanded file 142/67. Historical diversity values describe the screened 35-sequence set. | Figure S4 uses the expanded table; primary 33-sequence and historical 35-sequence statistics are distinguished. |
| Stored cutoff sweeps are not numerically invariant | AFP mean remodeling changes over the stored 4–6 Å sweep. | Figure S1 displays actual saved summaries; blanket stability claims are omitted. |

These are reporting corrections and disclosed provenance limitations. The original frozen
specification, analysis code and raw outputs remain intact. Exact adherence to every written
implementation detail is not claimed. The manuscript is a reviewable scientific draft;
uncertainty in unarchived structural-copy and state-specific provenance remains visible.

## Figure and table selection

| Item | Role |
|---|---|
| Figure 1 / Table 1 | Structural distributions and paired experimental references for all 435 eligible sites. |
| Figure 2 / Table 2 | Compact inclusion decisions and three primary protein datasets. |
| Figure 3 / Table 3 | Three added-variable panels ordered Afp2, Afp3, diffocin, plus primary and variable-site estimates. |
| Figure S1 / Table S1 | Hypothetical metric explanation, saved structural QC, chain sets and complete residue mappings. |
| Figure S2 / Table S2 | Existing R2 copy/locus diagnostics and recorded retrieval/orthology funnels. |
| Figure S3 / Table S3 | Cached AFP locus architecture and diagnostic joint tree; alignment, occupancy, variability and diversity. |
| Figure S4 / Table S4 | Saved SalCIS/diffocin locus evidence; full exported primary coefficients/SEs/VIFs. |
| Figure S5 / Table S5 | Existing descriptive remodeling/rate uncertainty plot; all eight saved sensitivity specifications. |

No validated atomic rendering or residue-coloured structure was archived, so such panels
were not invented to fill Figure 1. Figures S2 and S3b reuse existing images. New drawings
and graphs display saved values or explicitly hypothetical geometry. Figure S5 is descriptive,
not a new residual diagnostic, power test or measurement-error correction.

## Validation and practical limits

All primary estimates and supplementary coefficient tables were checked against full-precision
saved JSONs. The lines in Figure 3 use those exact coefficients; their slopes match the
stored residual identity to numerical tolerance without refitting. Counts of sites, sequences,
variable columns and bootstrap windows agree with cached inputs. The original-results hash
snapshot is retained outside the project, and the supplementary manifest identifies copied
source inputs. Reference numbers and DOI links, TeX labels/braces/environments, figure paths
and ZIP integrity were checked. All figure PNGs were visually reviewed, with geometric text
overlap checks on regenerated vector plots.

No TeX compiler is installed locally. The manuscript PDF and final pagination therefore remain
unverified; compile twice in Overleaf or a TeX environment using `manuscript.tex` as the main
file. Supplementary LaTeX is included automatically. The LaTeX ZIP contains selected processed
outputs and provenance files. The separate deposition ZIP includes retained sequence and
decision inputs, but omits external tool binaries and cannot restore missing atomic inputs.
Older figure files in the project remain historical artifacts and are not included in the
current three-figure package.

## Source verification


The numbers below correspond to the manuscript reference list. Bibliographic information
and the claims supported were checked against primary papers or their authoritative records.

| Ref. | Claim or method supported | Verification source |
|---|---|---|
| 1 | R2 pyocin extended/contracted architecture and sheath mechanics | [Ge et al., PubMed](https://pubmed.ncbi.nlm.nih.gov/25822993/) |
| 2 | AFP architecture, Afp2/Afp3, and intersubunit rearrangement | [Desfosses et al., paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6817355/); [6RBN](https://www.rcsb.org/structure/6RBN); [6RC8](https://www.rcsb.org/structure/6RC8) |
| 3 | Diffocin architecture and contraction | [Cai et al., published paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11329794/); [8V3Y](https://www.rcsb.org/structure/8V3Y) |
| 4 | SalCIS architecture and conformational pair | Published JATS source `g.xml`; [9R9H](https://www.rcsb.org/structure/9R9H); [9R9I](https://www.rcsb.org/structure/9R9I) |
| 5 | Solvent accessibility, packing and interface constraint | [Franzosa and Xia, PubMed](https://pubmed.ncbi.nlm.nih.gov/19597162/) |
| 6 | Packing versus accessibility as rate predictors; WCN | [Yeh et al., paper](https://doi.org/10.1155/2014/572409) |
| 7 | Packing, flexibility and mutational stress; nonlinear rate–flexibility relation | [Huang et al., paper](https://link.springer.com/article/10.1186/1471-2148-14-78) |
| 8 | Rate4Site algorithm | [Pupko et al., PubMed](https://pubmed.ncbi.nlm.nih.gov/12169533/) |
| 9 | Empirical Bayes estimation of site-specific rates | [Mayrose et al., PubMed](https://pubmed.ncbi.nlm.nih.gov/15201400/) |
| 10 | Pfam family annotations | [Paysan-Lafosse et al., paper](https://academic.oup.com/nar/article/53/D1/D523/7900195) |
| 11 | UniProtKB sequences and annotations | [UniProt Consortium, PubMed](https://pubmed.ncbi.nlm.nih.gov/39552041/) |
| 12 | BLAST+ | [Camacho et al., PubMed](https://pubmed.ncbi.nlm.nih.gov/20003500/) |
| 13 | FAMSA alignment algorithm | [Deorowicz et al., PubMed](https://pubmed.ncbi.nlm.nih.gov/27670777/) |
| 14 | IQ-TREE 3 | [Wong et al., paper](https://doi.org/10.1093/molbev/msag117); [official repository citation](https://github.com/iqtree/iqtree3) |
| 15 | ModelFinder | [Kalyaanamoorthy et al., PubMed](https://pubmed.ncbi.nlm.nih.gov/28481363/) |
| 16 | Ultrafast bootstrap approximation | [Hoang et al., PubMed](https://pubmed.ncbi.nlm.nih.gov/29077904/) |
| 17 | SH-like approximate likelihood-ratio support | [Guindon et al., PubMed](https://pubmed.ncbi.nlm.nih.gov/20525638/) |
| 18 | WAG amino-acid substitution matrix | [Whelan and Goldman, PubMed](https://pubmed.ncbi.nlm.nih.gov/11319253/) |
| 19 | Solvent-accessible area algorithm | [Shrake and Rupley, PubMed](https://pubmed.ncbi.nlm.nih.gov/4760134/) |
| 20 | Theoretical maximum solvent-accessible areas | [Tien et al., paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635) |
| 21 | Gaussian network fluctuation model | [Bahar et al., PubMed](https://pubmed.ncbi.nlm.nih.gov/9218955/) |
| 22 | Bootstrap motivation; not an exact description of the implemented window scheme | [Efron, original article scan](https://blogs.helsinki.fi/bk-club/files/2012/05/Efron_Bootstrap_AS1979.pdf); [publisher DOI](https://doi.org/10.1214/aos/1176344552) |
| 23 | Resampling dependent observations; stationary-block assumptions are not claimed | [Künsch, publisher record](https://doi.org/10.1214/aos/1176347265); [JSTOR original record](https://www.jstor.org/stable/2241719) |

Some direct publisher/PMC requests returned access challenges. Indexed primary-paper text,
PubMed bibliographic records, author-hosted primary articles and PDB primary-citation
records were used as appropriate; inaccessible pages were not treated as read full texts.
The gap is framed as information absent from the included monomer descriptor summaries, rather
than as a claim that all prior work on conformational dynamics and protein evolution is
absent. Rate4Site's estimates are relative amino-acid substitution rates, not dN/dS or a
direct test of positive selection.


The two statistical citations ground bootstrap motivation; the implementation's exact mechanics
are stated from its code rather than attributed to those papers. Publisher access challenges
were not treated as full-text access. Public reference sequences and structures are distinct
from the processed study files included in the package.

## Remaining submission work

Compile and review the revised manuscript PDF, confirm portal upload formats, and supply a
persistent public data/code repository plus archived release identifier. The current target
guidance is recorded in `SUBMISSION_CHECKLIST.md`; `title_page.tex` and the review sources
provide alternative identifying-page arrangements.
Neither a public accession nor a completed PDF build has been invented. The implementation
mismatches above require explicit scientific review as disclosed limitations; this revision
does not reopen computational analysis or silently rewrite its specification.
