#!/usr/bin/env python3
"""Convert reviewed Markdown to editable LaTeX and package saved presentation assets."""
from pathlib import Path
import re, zipfile, shutil
root=Path(__file__).resolve().parent
source=(root/'MANUSCRIPT.md').read_text()
supp_source=(root/'SUPPLEMENTARY.md').read_text()
# Keep numbered references in order of first appearance, including their DOI links.
body,ref_section=source.split('## References\n',1)
refs=dict(re.findall(r'^(\d+)\. (.+)$',ref_section,re.M))
order=[]
for n in re.findall(r'<sup>(.*?)</sup>',body,re.S):
    for i in re.findall(r'\[(\d+)\]',n):
        if i not in order:order.append(i)
assert set(order)==set(refs) and len(order)==23
renumber={old:str(i) for i,old in enumerate(order,1)}
def renumber_citations(t):
    return re.sub(r'<sup>(.*?)</sup>',lambda m:'<sup>'+re.sub(r'\[(\d+)\]',lambda n:'['+renumber[n[1]]+']',m[1])+'</sup>',t,flags=re.S)
new_refs=[(str(i),refs[n]) for i,n in enumerate(order,1)]
source=renumber_citations(body)+'## References\n\n'+'\n\n'.join(n+'. '+r for n,r in new_refs)+'\n'
(root/'MANUSCRIPT.md').write_text(source)
supp_source=renumber_citations(supp_source);(root/'SUPPLEMENTARY.md').write_text(supp_source)
key_for={n:'ref'+n for n,r in new_refs}
fig_labels={'1':'fig:remodeling','2':'fig:datasets','3':'fig:association',**{f'S{i}':f'fig:supp{i}' for i in range(1,6)}}
table_labels={'1':'tab:structures','2':'tab:datasets','3':'tab:primary',**{f'S{i}':f'tab:supp{i}' for i in range(1,6)}}
def inline(t):
    protected=[]
    def hold(v):
        token=f'ZZLATEXPLACEHOLDER{len(protected)}ZZ';protected.append(v);return token
    t=re.sub(r'<sup>(.*?)</sup>',lambda m:hold(r'\cite{'+','.join(key_for[n] for n in re.findall(r'\[(\d+)\]',m[1]))+'}'),t,flags=re.S)
    # References in prose link to native figure/table labels.
    t=re.sub(r'(?:\*\*)?(?:Supplementary )?Figure (S[1-5]|[1-3])(?:\*\*)?(?!\d)',lambda m:hold('Figure~'+r'\ref{'+fig_labels[m[1]]+'}'),t)
    t=re.sub(r'(?:\*\*)?(?:Supplementary )?Table (S[1-5]|[1-3])(?:\*\*)?(?!\d)',lambda m:hold('Table~'+r'\ref{'+table_labels[m[1]]+'}'),t)
    t=t.replace('|`C_E` ∪ `C_C`| < 3',hold(r'\(\lvert C_E \cup C_C\rvert < 3\)'))
    def code(m):
        value=m[1]
        if value in ['R_i','r_i','C_E','C_C']:return hold(r'\('+value+r'\)')
        return hold(r'\path{'+value+'}')
    t=re.sub(r'`([^`]+)`',code,t)
    t=re.sub(r'\[([^\]]+)\]\((https://(?:[^()]|\([^()]*\))+)\)',lambda m:hold(r'\href{'+m[2]+'}{'+inline(m[1])+'}'),t)
    t=t.replace('ior23@itu.edu.tr',hold(r'\href{mailto:ior23@itu.edu.tr}{ior23@itu.edu.tr}'))
    for old,tex in [('β_R',r'\(\beta_R\)'),('β_0',r'\(\beta_0\)'),('ε_i',r'\(\epsilon_i\)'),('β_WCN',r'\(\beta_{\mathrm{WCN}}\)'),('R²',r'\(R^2\)'),('10⁻⁵',r'\(10^{-5}\)'),('10⁻⁶',r'\(10^{-6}\)')]:t=t.replace(old,hold(tex))
    escape={'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}'}
    symbols={'±':r'\ensuremath{\pm}','²':r'\ensuremath{^2}','·':r'\ensuremath{\cdot}','Å':r'\AA{}','×':r'\ensuremath{\times}','ś':r"\'{s}",'ü':r'\"{u}','Γ':r'\ensuremath{\Gamma}','Δ':r'\ensuremath{\Delta}','α':r'\ensuremath{\alpha}','β':r'\ensuremath{\beta}','ε':r'\ensuremath{\varepsilon}','ρ':r'\ensuremath{\rho}','–':'--','—':'---','−':r'\ensuremath{-}','∩':r'\ensuremath{\cap}','∪':r'\ensuremath{\cup}','≈':r'\ensuremath{\approx}','≤':r'\ensuremath{\leq}','≥':r'\ensuremath{\geq}'}
    t=''.join(escape.get(c,symbols.get(c,c)) for c in t)
    t=re.sub(r'\*\*(.+?)\*\*',lambda m:r'\textbf{'+m[1]+'}',t,flags=re.S)
    t=re.sub(r'\*([^*]+)\*',lambda m:r'\emph{'+m[1]+'}',t,flags=re.S)
    for i,value in enumerate(protected):t=t.replace(f'ZZLATEXPLACEHOLDER{i}ZZ',value)
    assert 'ZZLATEXPLACEHOLDER' not in t
    return t
main_table_captions=[
'Structural datasets. Extended (E) and contracted (C) PDB entities identify the five sheath references. Example author chain identifiers are case-sensitive deposited family members; selected contact-reference subunit keys were not archived.',
'Primary evolutionary datasets. Variable sites contain at least two distinct non-gap amino acids. The 223 fitted sites sum protein-specific samples from two systems, with two AFP proteins in one system.',
'Adjusted remodeling associations and variable-site restrictions. Standardized coefficients describe rate SD per remodeling SD conditional on monomer RSA, WCN and GNM averaged across states. Intervals use ten-reference-position sequence windows; all intervals are sequence-block bootstrap estimates and variable-site documentation scope is explained in the Methods and Table S5.']
supp_table_captions=['Structural references, cross-state modeled sequences and eligibility. E/C denotes extended/contracted; copy counts include both AFP families where deposited operators are applied. Full author chain sets and residue mappings accompany the table.',
'Retrieval funnels and orthology decisions. AFP uses a joint prescreen rather than a directly comparable single-query gated funnel. Counts are screening-stage observations.',
'Alignment metadata, fitted-site occupancy and variability. Values come from saved alignments and site tables; the alternate diffocin sets are sensitivity datasets.',
'All exported non-intercept primary coefficients, independent-site OLS standard errors, intervals and p-values, and variance inflation factors. Intercept inferential statistics were not exported.',
'Complete archived sensitivity results. An asterisk marks Afp3 and diffocin variable-site checks beyond the literal Afp2-only written clause. Other specifications are named in the written sensitivity plan. Bootstrap intervals follow the primary sequence-window scheme.']
counts={'main':0,'supp':0}
def table(lines,kind):
    rows=[[c.strip() for c in line.strip().strip('|').split('|')] for line in lines]
    rows=[r for r in rows if not all(re.fullmatch(r':?-+:?',c) for c in r)]
    assert all(len(r)==len(rows[0]) for r in rows)
    counts[kind]+=1;n=counts[kind];cap=(main_table_captions if kind=='main' else supp_table_captions)[n-1]
    label=table_labels[str(n) if kind=='main' else f'S{n}']
    header=' & '.join(r'\textbf{'+inline(c)+'}' for c in rows[0])+r' \\'
    cells=[' & '.join(inline(c) for c in row)+r' \\' for row in rows[1:]]
    if kind=='main':
        columns={1:r'@{}>{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}p{0.10\textwidth}>{\raggedright\arraybackslash}p{0.15\textwidth}>{\raggedright\arraybackslash}p{0.15\textwidth}>{\raggedright\arraybackslash}p{0.16\textwidth}r@{}',2:r'@{}>{\raggedright\arraybackslash}Xrrr@{}',3:r'@{}>{\raggedright\arraybackslash}X*{3}{>{\centering\arraybackslash}p{0.205\textwidth}}@{}'}[n]
        return '\n'.join([r'\begin{table}[htbp]',r'\centering',r'\caption{'+cap+'}',r'\label{'+label+'}',r'\footnotesize',r'\setlength{\tabcolsep}{3pt}',r'\renewcommand{\arraystretch}{1.25}',r'\begin{tabularx}{\textwidth}{'+columns+'}',r'\toprule',header,r'\midrule',*cells,r'\bottomrule',r'\end{tabularx}',r'\end{table}'])
    weights={1:[.19,.25,.15,.13,.13,.15],2:[.16,.08,.15,.13,.12,.36],3:[.26,.13,.14,.16,.12,.19],4:[.16,.10,.10,.10,.23,.10,.11],5:[.22,.09,.07,.08,.07,.09,.18,.10,.10]}[n]
    columns=''.join(r'>{\raggedright\arraybackslash}p{\dimexpr '+str(w)+r'\linewidth-2\tabcolsep\relax}' for w in weights)
    pre=[r'\clearpage',r'\begin{landscape}'] if n==5 else []
    post=[r'\end{landscape}'] if n==5 else []
    return '\n'.join(pre+[r'\begingroup',r'\footnotesize',r'\setlength{\tabcolsep}{3pt}',r'\renewcommand{\arraystretch}{1.2}',r'\begin{longtable}{'+columns+'}',r'\caption{'+cap+r'}\label{'+label+r'}\\',r'\toprule',header,r'\midrule',r'\endfirsthead',r'\multicolumn{'+str(len(weights))+r'}{l}{\tablename~\thetable{} (continued)}\\',r'\toprule',header,r'\midrule',r'\endhead',r'\bottomrule',r'\endfoot',*cells,r'\end{longtable}',r'\endgroup']+post)
def blocks(t,kind='main'):
    lines=t.strip().splitlines();out=[];i=0;listing=False
    def close():
        nonlocal listing
        if listing:out.append(r'\end{enumerate}');listing=False
    while i<len(lines):
        line=lines[i]
        if not line.strip():i+=1;continue
        if line.startswith('## '):
            close();out.append(r'\section{'+inline(line[3:])+'}');i+=1;continue
        if line.startswith('### '):
            close();heading=re.sub(r'^Table S[1-5]\. ','',line[4:]) if kind=='supp' else line[4:]
            out.append(r'\subsection{'+inline(heading)+'}');i+=1;continue
        if line.startswith('|'):
            close();group=[]
            while i<len(lines) and lines[i].startswith('|'):group.append(lines[i]);i+=1
            out.append(table(group,kind));continue
        if line.strip() == '```math':
            close(); i+=1; math=[]
            while i<len(lines) and lines[i].strip()!='```':
                math.append(lines[i]); i+=1
            assert i<len(lines), 'Unclosed math block'
            i+=1; content='\n'.join(math).strip()
            label=('eq:regression' if 'z(r_i)' in content else
                   'eq:remodeling' if content.startswith('R_i') else
                   'eq:partialr2' if content.startswith('R^2') else None)
            equation=[r'\begin{equation}',content]
            if label:equation.append(r'\label{'+label+'}')
            equation.append(r'\end{equation}');out.append('\n'.join(equation));continue
        if line.startswith('    R_i ='):
            out.append(r'''\begin{equation}
R_i=1-\frac{\lvert C_i^E\cap C_i^C\rvert}{\lvert C_i^E\cup C_i^C\rvert}.
\label{eq:remodeling}
\end{equation}''');i+=1;continue
        if line.startswith('    z(r_i) ='):
            out.append(r'''\begin{equation}
\boxed{\begin{aligned}
z(r_i)=\beta_0 &+\beta_R z(R_i)+\beta_A z(\mathrm{RSA}_i)\\
&+\beta_P z(\mathrm{WCN}_i)+\beta_F z(\mathrm{GNM}_i)+\varepsilon_i .
\end{aligned}}
\label{eq:regression}
\end{equation}''');i+=1;continue
        if line.startswith('    R²_partial'):
            out.append(r'''\begin{equation}
R^2_{\mathrm{partial},R}=\frac{R^2_{\mathrm{full}}-R^2_{\mathrm{reduced}}}{1-R^2_{\mathrm{reduced}}}.
\label{eq:partialr2}
\end{equation}''');i+=1;continue
        item=re.match(r'^\d+\.\s+(.*)',line)
        if item:
            if not listing:out.append(r'\begin{enumerate}');listing=True
            group=[item[1]];i+=1
            while i<len(lines) and lines[i].strip() and not re.match(r'^\d+\.\s',lines[i]):group.append(lines[i].strip());i+=1
            out.append(r'\item '+inline(' '.join(group)));continue
        close();group=[line.strip()];i+=1
        while i<len(lines) and lines[i].strip() and not re.match(r'^(##|\||\d+\.\s|    )',lines[i]):group.append(lines[i].strip());i+=1
        out.append(inline(' '.join(group)))
    close();return '\n\n'.join(out)
def figure(images,caption,label):
    return '\n'.join([r'\begin{figure}[p]',r'\centering',*[(r'\includegraphics[width=\textwidth,height='+height+r'\textheight,keepaspectratio]{'+path+'}'+r'\par\smallskip') for path,height in images],r'\caption{'+inline(caption)+'}',r'\label{'+label+'}',r'\end{figure}'])
title=source.splitlines()[0][2:]
abstract_section=source.split('## Abstract\n',1)[1].split('## Introduction\n',1)[0];abstract,keywords=abstract_section.split('**Keywords:**',1)
main='## Introduction\n'+source.split('## Introduction\n',1)[1].split('## Figures\n',1)[0]
fig_section=source.split('## Figures\n',1)[1].split('## Data and code availability\n',1)[0]
availability=source.split('## Data and code availability\n',1)[1].split('## References\n',1)[0]
preamble=r'''% Generated from MANUSCRIPT.md and SUPPLEMENTARY.md; saved numerical outputs unchanged.
% Compile twice with pdfLaTeX, or use latexmk -pdf manuscript.tex.
% Bibliography is embedded: no BibTeX/Biber step is required.
\documentclass[11pt,a4paper]{article}
\usepackage[margin=25mm]{geometry}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern,microtype}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs,tabularx,array,longtable,pdflscape}
\usepackage{caption}
\usepackage[super,numbers,sort&compress]{natbib}
\usepackage[hyphens]{url}
\usepackage[hidelinks]{hyperref}
\hypersetup{pdftitle={Contact-network remodeling during bacterial sheath contraction},pdfauthor={Anthony Yua Ior}}
\captionsetup{font=small,labelfont=bf}
\setlength{\parindent}{1.2em}
\setlength{\parskip}{0.35em}
\setlength{\emergencystretch}{2em}
\renewcommand{\bibsection}{\section*{References}}
\setlength{\bibsep}{0.5em}
'''
parts=[preamble,r'\title{'+inline(title)+'}',r'''\author{Anthony Yua Ior\thanks{Correspondence: \href{mailto:ior23@itu.edu.tr}{ior23@itu.edu.tr}.}\\[0.5em]
\small Istanbul Technical University\\
\small Maslak 34867, Istanbul, Turkiye}
\date{}
\begin{document}
\maketitle
\begin{abstract}''',blocks(abstract),r'\end{abstract}',r'\noindent\textbf{Keywords:} '+inline(' '.join(keywords.split())),blocks(main),r'\section*{Data and code availability}',blocks(availability),r'\clearpage',r'\section*{Figures}']
figure_files=[]
for num,cap in re.findall(r'\*\*Figure ([1-3])\.\*\*\s*(.*?)(?=\n!\[)',fig_section,re.S):
    cap=' '.join(re.split(r'\s+File:',cap)[0].split());path=f'figures/Figure{num}.pdf';figure_files.append(path)
    parts.append(figure([(path,'0.70')],cap,fig_labels[num]))
parts += [r'\clearpage',r'\begin{thebibliography}{99}']
for n,r in new_refs:parts.append(r'\bibitem{'+key_for[n]+'}\n'+inline(r))
parts += [r'\end{thebibliography}',r'\clearpage',r'\input{supplementary.tex}',r'\end{document}']
supp_body=supp_source.split('## Supplementary Figures\n',1)[0].split('\n',1)[1]
supp_parts=[r'''% Included after the main bibliography; supplementary citations use its reference keys.
\section*{Supplementary Methods and Results}
\setcounter{section}{0}
\setcounter{subsection}{0}
\setcounter{table}{0}
\setcounter{figure}{0}
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\theHsection}{supp.\arabic{section}}
\renewcommand{\theHtable}{supp.\arabic{table}}
\renewcommand{\theHfigure}{supp.\arabic{figure}}
''',blocks(supp_body,'supp'),r'\clearpage',r'\section*{Supplementary Figures}']
supp_images={'S1':[('figures/FigureS1a.pdf','0.34'),('figures/FigureS1b.pdf','0.37')],'S2':[('figures/FigureS2a.png','0.30'),('figures/FigureS2b.png','0.38')],'S3':[('figures/FigureS3a.pdf','0.14'),('figures/FigureS3b.png','0.58')],'S4':[('figures/FigureS4.pdf','0.70')],'S5':[('figures/FigureS5.pdf','0.65')]}
for num,cap in re.findall(r'\*\*Figure (S[1-5])\.\*\*\s*(.*?)(?=\n\n\*\*Figure|\Z)',supp_source.split('## Supplementary Figures\n',1)[1],re.S):
    cap=' '.join(re.split(r'\s+File:',cap)[0].split());supp_parts.append(figure(supp_images[num],cap,fig_labels[num]));figure_files += [p for p,h in supp_images[num]]
assert counts=={'main':3,'supp':5},counts
for name,parts_to_write in [('manuscript.tex',parts),('supplementary.tex',supp_parts)]:
    tex='\n\n'.join(parts_to_write)+'\n';assert tex.isascii(),sorted(set(c for c in tex if ord(c)>127));(root/name).write_text(tex)
# Separate title page and review files follow the journal's identifying-page guidance.
full_tex=(root/'manuscript.tex').read_text()
review_tex=re.sub(r'\\author\{.*?\}\n\\date\{\}',lambda m:r'\author{}'+'\n'+r'\date{}',full_tex,flags=re.S)
review_tex=review_tex.replace('pdfauthor={Anthony Yua Ior}','pdfauthor={}')
review_tex=re.sub(r'\\section\{Acknowledgments\}.*?(?=\\section\{Statements and Declarations\})','',review_tex,flags=re.S)
review_tex=re.sub(r'\\subsection\{Corresponding author\}.*?(?=\\clearpage)','',review_tex,flags=re.S)
review_tex=review_tex.replace(r'\input{supplementary.tex}',r'\input{supplementary_for_review.tex}')
review_supp=re.sub(r'^Anthony Yua Ior,[^\n]+\n+','',(root/'supplementary.tex').read_text(),flags=re.M)
# Repository-owner URLs belong on the identifying title page, not the review variants.
def remove_identifying_repo_links(tex):
    pattern=r'\\href\{https://github\.com/Anthony-Pius-1729/sheath_evolution[^}]*\}\{([^}]+)\}'
    return re.sub(pattern,lambda m:'source repository (details on title page)' if 'repository' in m[1] else m[1],tex)
review_tex=remove_identifying_repo_links(review_tex).replace("so this link identifies the source location", "so these details identify the source location")
review_supp=remove_identifying_repo_links(review_supp)
for name,tex in [('manuscript_for_review.tex',review_tex),('supplementary_for_review.tex',review_supp)]:
    assert not any(v in tex for v in ['Anthony Yua Ior','ior23@itu.edu.tr','Istanbul Technical University','Maslak 34867','Anthony-Pius-1729'])
    (root/name).write_text(tex)
declaration_markdown='## Acknowledgments\n'+source.split('## Acknowledgments\n',1)[1].split('## References\n',1)[0]
title_page=preamble+'\n'+r'\title{'+inline(title)+'}\n'+r'''\author{Anthony Yua Ior\\[0.5em]
\small Istanbul Technical University\\
\small Maslak 34867, Istanbul, Turkiye\\
\small \href{mailto:ior23@itu.edu.tr}{ior23@itu.edu.tr}}
\date{}
\begin{document}
\maketitle
\section*{Data and code availability}
Code and retained study outputs are hosted at
\href{https://github.com/Anthony-Pius-1729/sheath_evolution}{the GitHub repository},
which is currently private. A local code/data archive has been prepared; public access and
an archived release/DOI remain pending. Original coordinate files and selected reference
keys are absent; the provenance audit documents these limits.
'''+blocks(declaration_markdown)+'\n'+r'\end{document}'+'\n'
assert title_page.isascii()
(root/'title_page.tex').write_text(title_page)

readme='''# LaTeX manuscript package

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
The Overleaf ZIP contains one standalone TeX document: `manuscript.tex`. It includes
`supplementary.tex` automatically. Separate title-page and review documents remain in the
local project and are deliberately excluded from this ZIP so the two-page title document
cannot be selected instead of the paper. `manuscript.tex` is the complete author version.

If an existing project shows only acknowledgments and declarations, open `manuscript.tex`,
choose File > Settings > Compiler > Main document > manuscript.tex, and Recompile.
Alternatively upload this ZIP as a NEW Overleaf project; uploading into an old project does
not remove its previously uploaded title-page source. The journal prefers Word and does
not mandate a template; LaTeX upload handling still needs confirmation in the portal.
`SUBMISSION_CHECKLIST.md` and `COVER_LETTER_DRAFT.md` record submission preparation.
`REVIEW_RESPONSE.md` describes the changes made in response to the supplied review.
`REPRODUCIBILITY_AUDIT.md` documents what the archive can and cannot reproduce.
The separate `reproducibility_archive.zip` is local and prepared for deposition, not public.
To produce the main paper without the supplement, comment out `\\input{supplementary.tex}`.
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
final page-layout verification remain pending. Code and retained data have been pushed to the private GitHub repository
https://github.com/Anthony-Pius-1729/sheath_evolution. Manuscript files remain local.
Public access and an archived release/DOI remain pending. The bundle contains selected processed outputs,
not the full retrieval cache, external tools or complete runnable analysis environment.
'''
(root/'LATEX_README.md').write_text(readme)
paths=['manuscript.tex','supplementary.tex','MANUSCRIPT.md','SUPPLEMENTARY.md','MANUSCRIPT_REVIEW.md','REVIEW_RESPONSE.md','REPRODUCIBILITY_AUDIT.md','structural_provenance_audit.json','SUBMISSION_CHECKLIST.md','COVER_LETTER_DRAFT.md','prepare_supplementary.py','validate_manuscript.py','build_reproducibility_archive.py','make_figures.py','write_manuscript_latex.py']+figure_files
if (root/'MANUSCRIPT_VALIDATION.json').exists():paths.append('MANUSCRIPT_VALIDATION.json')
if (root/'REPRODUCIBILITY_ARCHIVE.json').exists():paths.append('REPRODUCIBILITY_ARCHIVE.json')
with zipfile.ZipFile(root/'manuscript_latex.zip','w',zipfile.ZIP_DEFLATED) as z:
    z.write(root/'LATEX_README.md','README.md')
    for p in paths:assert (root/p).exists(),p;z.write(root/p,p)
    for p in sorted((root/'supplementary_data').rglob('*')):
        if p.is_file():z.write(p,str(p.relative_to(root)))
folder=root/'manuscript_latex';folder.mkdir(exist_ok=True)
# Preserve superseded extracted copies as non-TeX backups; original local sources stay intact.
for name in ['title_page.tex','manuscript_for_review.tex','supplementary_for_review.tex']:
    old=folder/name
    if old.exists():
        backup=folder/'optional_source_backups';backup.mkdir(exist_ok=True)
        dest=backup/(name+'.txt');index=1
        while dest.exists():dest=backup/(name+'.'+str(index)+'.txt');index+=1
        old.rename(dest)
with zipfile.ZipFile(root/'manuscript_latex.zip') as z:z.extractall(folder)
with zipfile.ZipFile(root/'manuscript_latex.zip') as z:
    standalone=[n for n in z.namelist() if n.endswith('.tex') and r'\documentclass' in z.read(n).decode()]
    assert standalone==['manuscript.tex'], standalone
print('Created LaTeX source and ZIP: three main figures, three main tables, five supplementary figures/tables, 23 references.')
