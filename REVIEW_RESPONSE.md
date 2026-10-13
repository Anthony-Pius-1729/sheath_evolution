# Response to the manuscript review

This revision addresses the supplied methods/reproducibility review. The scientific analysis
remains closed: no structural calculations, rate estimates, tree inference or new regression
specifications were run. All primary and archived sensitivity estimates remain unchanged.

| Review item | Revision and disposition |
|---|---|
| Covariate-state discrepancy | Methods now lead with the implemented available-state arithmetic means, distinguish descriptor averaging from coordinate averaging, and distinguish representative monomer chains from contact-reference subunits. The abstract, estimand, results, Figure 3 caption, Discussion and Conclusion identify the actual adjustment. The extended-state-only question remains untested; the mismatch is disclosed, not resolved by claiming equivalence. |
| Statistical uncertainty | Interval labels identify sequence-block bootstrap estimates. Methods/Discussion state that coverage under remaining spatial dependence was not assessed, and that structures, orthology, trees and rate posteriors are held fixed. Positive AFP intervals are interpreted within those resampling assumptions. |
| Selected structural reference subunits | A source/input audit documents ordered chain/operator construction, interior-ring preference, <90 Å centroid scoring and first-in-order tie handling. The original mmCIF bytes, reference keys and tie scores are missing. Exact historical selected subunits therefore cannot be verified from retained outputs. This limitation remains explicit; no example chain is substituted. |
| Figure 1 biological communication | Results and caption connect the measured distributions to exchange of intermolecular partner identities while the sheath remains connected. They identify the lack of spatial localization and avoid force/energy interpretations. No unvalidated structure panel or pipeline rerun was introduced. |
| Figure 2 priority | Retained in the main text as a concise explanation of attrition from four systems to two systems/three protein datasets. The target journal currently has no fixed figure/table cap. Figure 3 retains the central adjusted-association role. |
| Local code/data archive | `reproducibility_archive.zip` is prepared for deposition with code, recorded inputs, outputs, decision evidence and SHA-256 inventory. It is local only: no public URL or DOI has been assigned. Its README identifies missing structural input provenance and external/environment requirements. |
| Journal preparation | Declaration headings are aligned with the target guidance. Separate title-page and review-copy LaTeX sources are supplied; the full author version remains the editable master. A draft cover letter and submission checklist identify remaining factual/build/deposition items. |

The strongest contribution remains the adjusted linear association in AFP, most convincingly
in Afp3; Afp2 is unresolved on variable sites, and an association was not detected in diffocin.
Detectability differences do not establish a between-system effect difference. Greater relative
rate uncertainty and narrower diffocin remodeling remain unresolved alternatives to biology.

The journal assessment uses the current [Evolutionary Bioinformatics instructions](https://journals.sagepub.com/author-instructions/evb).
The manuscript is revised for review, not declared ready for submission: public deposition,
compiled-PDF/pagination review and unresolved structural provenance remain outstanding.
