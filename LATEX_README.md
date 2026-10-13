# LaTeX manuscript package

Anthony Yua Ior — Istanbul Technical University, Maslak 34867, Istanbul, Turkiye.
Correspondence: ior23@itu.edu.tr. Funding: none. Competing interests: none.
Author contributions: study design, software/analysis, interpretation and manuscript preparation.

The source contains three main figures, three main tables and 23 verified references.
`supplementary.tex` contains supplementary methods/results, Tables S1–S5 and Figures S1–S5.
It is included automatically by `manuscript.tex`; keep both files and `figures/` together.
The complete-precision table CSVs, mappings, exact saved commands, selected original outputs
and source checksums are in `supplementary_data/`.

Upload `manuscript_latex.zip` to Overleaf, set `manuscript.tex` as the main document and
select pdfLaTeX. Bibliography is embedded; BibTeX/Biber is not needed. For local compilation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

The article class and common packages support editing; no publisher class is assumed.
`title_page.tex` is the separate identifying page. `manuscript_for_review.tex` includes
`supplementary_for_review.tex` and omits author identity from document text and PDF-author
metadata. Use these alternative files according to the submission system's requirements;
`manuscript.tex` remains the complete author version. The journal prefers Word and does
not mandate a template; LaTeX upload handling still needs confirmation in the portal.
`SUBMISSION_CHECKLIST.md` and `COVER_LETTER_DRAFT.md` record submission preparation.
`REVIEW_RESPONSE.md` describes the changes made in response to the supplied review.
`REPRODUCIBILITY_AUDIT.md` documents what the archive can and cannot reproduce.
The separate `reproducibility_archive.zip` is local and prepared for deposition, not public.
To produce the main paper without the supplement, comment out `\input{supplementary.tex}`.
The supplement depends on the main preamble and is not a standalone document.

The scientific analysis remains frozen. Presentation scripts read saved outputs:
`prepare_supplementary.py` prepares supplementary tables and provenance,
`make_figures.py` exports presentation figures, and `write_manuscript_latex.py` converts
reviewed Markdown and packages selected assets. The included source can compile without
running those scripts. Regeneration scripts require the complete original project paths.

Reporting distinguishes operator-expanded AFP coordinates from guessed ring completion,
discloses averaged monomer covariates versus the written extended-state definition,
identifies unlogged reference-copy keys, and identifies parallel variable-site checks
beyond the literal Afp2-only clause. None of these issues was corrected by refitting.
Afp2's weaker variable-site evidence, diffocin non-detection, the sequence dependence
limitation and unpropagated rate uncertainty remain visible.

Source consistency, references, saved numbers, archive integrity and figure images have
been checked. No TeX compiler is installed locally, so manuscript PDF compilation and
final page-layout verification remain pending. Public repository deposition is also pending;
no repository accession has been invented. The bundle contains selected processed outputs,
not the full retrieval cache, external tools or complete runnable analysis environment.
