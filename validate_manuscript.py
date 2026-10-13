#!/usr/bin/env python3
"""Verify manuscript reporting/source integrity against archived outputs; no model fitting."""
from pathlib import Path
import csv,json,re,hashlib,zipfile
ROOT=Path(__file__).resolve().parent
def rows(p):
    with (ROOT/p).open() as f:return list(csv.DictReader(f))
def syntax(name):
    t=(ROOT/name).read_text();assert t.isascii(),name
    depth=0
    for i,c in enumerate(t):
        if c not in '{}':continue
        j=i-1;n=0
        while j>=0 and t[j]=='\\':n+=1;j-=1
        if n%2:continue
        depth+=1 if c=='{' else -1;assert depth>=0,(name,i)
    assert depth==0,(name,depth)
    stack=[]
    for action,env in re.findall(r'\\(begin|end)\{([^}]+)\}',t):
        if action=='begin':stack.append(env)
        else:assert stack and stack.pop()==env,(name,env)
    assert not stack,(name,stack)
    for p in re.findall(r'\\includegraphics\[[^]]*\]\{([^}]+)\}',t):assert (ROOT/p).exists(),p
    for p in re.findall(r'\\input\{([^}]+)\}',t):assert (ROOT/p).exists(),p
    return t
tex={n:syntax(n) for n in ['manuscript.tex','supplementary.tex','manuscript_for_review.tex','supplementary_for_review.tex','title_page.tex']}
for main,supp in [('manuscript.tex','supplementary.tex'),('manuscript_for_review.tex','supplementary_for_review.tex')]:
    a=tex[main];b=tex[supp];both=a+b
    assert a.count(r'\begin{table}')==3 and b.count(r'\begin{longtable}')==5
    assert a.count(r'\begin{figure}')==3 and b.count(r'\begin{figure}')==5
    labels=re.findall(r'\\label\{([^}]+)\}',both);assert len(labels)==len(set(labels))
    assert set(re.findall(r'\\ref\{([^}]+)\}',both))<=set(labels)
    keys=re.findall(r'\\bibitem\{([^}]+)\}',both);assert len(keys)==len(set(keys))==23
    used={k for group in re.findall(r'\\cite\{([^}]+)\}',both) for k in group.split(',')};assert used==set(keys)
for p in ['manuscript_for_review.tex','supplementary_for_review.tex']:
    assert not any(v in tex[p] for v in ['Anthony Yua Ior','ior23@itu.edu.tr','Istanbul Technical University','Maslak 34867'])
s=(ROOT/'MANUSCRIPT.md').read_text();body,refs=s.split('## References\n');reference=dict(re.findall(r'^(\d+)\. (.+)$',refs,re.M))
first=[]
for group in re.findall(r'<sup>(.*?)</sup>',body,re.S):
    for n,doi in re.findall(r'\[(\d+)\]\((https://doi.org/[^)]*)\)',group):
        assert doi in reference[n]
        if n not in first:first.append(n)
assert first==[str(i) for i in range(1,24)]
assert len(re.findall(r'^### ',body.split('## Results\n')[1].split('## Discussion')[0],re.M))==3
primary={r['dataset']:r for r in json.loads((ROOT/'results/hypothesis_test_primary.json').read_text())}
name={'AFP Afp2':'afp_Afp2','AFP Afp3':'afp_Afp3','Diffocin':'diffocin'}
for row in rows('supplementary_data/TableS4_full_coefficients.csv'):
    m=primary[name[row['dataset']]];p=row['predictor']
    for col,key in [('beta','beta'),('SE','se'),('VIF','vif'),('OLS_p','p')]:assert float(row[col])==m[key][p]
    assert [float(row['OLS_CI_lo']),float(row['OLS_CI_hi'])]==m['ci'][p]
sens=json.loads((ROOT/'results/sensitivity.json').read_text());table=rows('supplementary_data/TableS5_sensitivity_results.csv');assert len(table)==len(sens)==8
for row in table:
    m=sens[row['analysis']]
    assert float(row['beta_R'])==m['beta']['R'] and float(row['SE_R'])==m['se']['R']
    assert [float(row['boot_CI_lo']),float(row['boot_CI_hi'])]==m['block_bootstrap']['ci95']
    assert float(row['partial_R2'])==m['partial_r2_R'] and float(row['full_R2'])==m['r2']
meta=rows('supplementary_data/TableS3_alignment_metadata.csv')
assert [(r['sequences'],r['fitted_sites'],r['variable_sites']) for r in meta[:3]]==[('27','84','51'),('31','81','64'),('33','58','50')]
# Check the main comparison table against independent full-precision saved results.
main_table=body.split('| Quantity | AFP Afp2 | AFP Afp3 | Diffocin |')[1].split('\n\n')[0]
for k in ['afp_Afp2','afp_Afp3','diffocin']:
    m=primary[k]
    assert f"{m['beta']['R']:+.3f}" in main_table
    lo,hi=m['block_bootstrap']['ci95'];assert f'[{lo:+.3f}, {hi:+.3f}]' in main_table.replace('−','-')
manifest=json.loads((ROOT/'supplementary_data/source_sha256.json').read_text())
for p,h in manifest.items():
    assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert hashlib.sha256((ROOT/'supplementary_data/archived'/p).read_bytes()).hexdigest()==h,p
baseline=Path('/private/tmp/sheath_methods_review_before/analysis_sha256.json')
if baseline.exists():
    orig=json.loads(baseline.read_text());assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in orig.items())
    protected=len(orig)
else:protected='session-specific baseline unavailable; source manifest checked'
# Confidence wording remains conditional on the actual scheme and adjustment inputs.
assert 'an extended-state-only adjusted association has not been tested' in body
assert 'coverage' in body and 'not spatially' in body and 'not assessed' in body
assert all(k in body for k in ['+0.435','+0.425','+0.035','−0.023','+0.189','−0.173'])
assert not list(ROOT.glob('sheath_pilot/structures/*.cif'))
audit=json.loads((ROOT/'structural_provenance_audit.json').read_text());assert audit['exact_historical_reference_subunits_verified'] is False
for p,h in audit['source_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
with zipfile.ZipFile(ROOT/'manuscript_latex.zip') as z:
    assert z.testzip() is None
    for p in ['manuscript.tex','supplementary.tex']:assert z.read(p)==(ROOT/p).read_bytes()
    standalone=[n for n in z.namelist() if n.endswith('.tex') and r'\documentclass' in z.read(n).decode()]
    assert standalone==['manuscript.tex'],standalone
    assert not {'title_page.tex','manuscript_for_review.tex','supplementary_for_review.tex'}&set(z.namelist())
    for p in ['MANUSCRIPT.md','SUPPLEMENTARY.md','REVIEW_RESPONSE.md']:assert z.read(p)==(ROOT/p).read_bytes()
def words(t):
    t=re.sub(r'<sup>.*?</sup>','',t,flags=re.S);return len('\n'.join(l for l in t.splitlines() if not l.startswith(('|','#','    '))).split())
abstract=body.split('## Abstract\n')[1].split('**Keywords:**')[0];main=body.split('## Introduction\n')[1].split('## Figures\n')[0]
assert words(abstract)<300 and 4000<=words(main)<=6000
report=dict(original_analysis_hashes_unchanged=protected,main_figures=3,main_tables=3,supplementary_figures=5,supplementary_tables=5,references=23,primary_fitted_sites=223,sensitivity_specifications=8,main_words=words(main),abstract_words=words(abstract),TeX_sources_checked=5,Overleaf_main_document='manuscript.tex',Overleaf_standalone_documents=1,TeX_source_syntax='passed static checks',review_source_identity='removed from text and PDF-author metadata',saved_coefficients='match original JSON',citation_order_DOIs='checked',ZIP_integrity='passed',actual_adjustment='monomer descriptors averaged across available states',statistical_intervals='sequence-block bootstrap; spatial coverage unassessed',structural_reference_reconstruction='not verified: original coordinates and selected keys absent',public_archive='local deposition archive; public URL/DOI pending',manuscript_PDF_compilation='not performed: no local TeX compiler')
(ROOT/'MANUSCRIPT_VALIDATION.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
