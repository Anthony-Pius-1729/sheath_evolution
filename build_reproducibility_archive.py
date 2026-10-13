#!/usr/bin/env python3
"""Package retained code/data for release, excluding local manuscript documents."""
from pathlib import Path
import csv,hashlib,io,json,zipfile,subprocess
ROOT=Path(__file__).resolve().parent
EXT={'.py','.csv','.json','.nwk','.afa','.fasta','.fa','.faa','.tsv','.txt','.log','.treefile','.res','.unnorm','.cmd','.stdout','.iqtree','.contree','.bionj','.png','.pdf','.gff','.gff3','.gb','.gbk','.stderr'}
roots=['orthology','phylo','rates','results','retrieval_uniprot','filtered','synteny','evolution_prelim','figures']
files=set()
for folder in roots:
    base=ROOT/folder
    if base.exists():
        files.update(p for p in base.rglob('*') if p.is_file() and p.suffix in EXT and (folder!='figures' or p.parent==base))
files.update(p for p in ROOT.glob('*.py') if p.is_file())
files.update(p for p in (ROOT/'sheath_pilot').glob('*.py') if p.is_file())
for folder in ['sheath_pilot/qc_output','sheath_pilot/output']:
    files.update(p for p in (ROOT/folder).rglob('*') if p.is_file() and p.suffix in EXT)
files.update(p for p in (ROOT/'supplementary_data').glob('*') if p.is_file())
for name in ['README.md','INDEX.md','DATA_DICTIONARY.md','RETAINED_DATA_VALIDATION.json','requirements.txt','reference_sequences.json','references.fasta','pair_meta.json','FROZEN_ANALYSIS.md','DATASET_DECISION.md','RESULTS.md','AFP_FEASIBILITY.md','DIFFOCIN_FEASIBILITY.md','R2_FINAL_DECISION.md','R2_EVOLUTIONARY_FEASIBILITY.md','REPRODUCIBILITY_AUDIT.md','structural_provenance_audit.json','sheath_pilot/README.md','sheath_pilot/QC_RESULTS.md']:
    p=ROOT/name
    if p.exists():files.add(p)
# The exact five final rate/alignment/topology commands must point to packaged inputs.
required=['hypothesis_test.py','sheath_pilot/structural_qc.py','sheath_pilot/lattice.py','sheath_pilot/covariates.py','sheath_pilot/remodeling.py']
for k in ['afp_Afp2','afp_Afp3','diffocin','diffocin_asscreened','diffocin_collinear']:
    import re
    command=(ROOT/f'rates/{k}.cmd').read_text()
    required += [re.search(r'-s (\S+)',command)[1],re.search(r'-t (\S+)',command)[1],f'rates/{k}.res',f'phylo/{k}.log']
assert all(ROOT/p in files for p in required)
rows=[]
for p in sorted(files):
    rel=str(p.relative_to(ROOT));data=p.read_bytes()
    rows.append(dict(path=rel,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
stream=io.StringIO();w=csv.DictWriter(stream,fieldnames=['path','bytes','sha256']);w.writeheader();w.writerows(rows)
readme='''# Sheath evolution: retained code and data, v1.0.0

Study: Contact-Network Remodeling During Bacterial Sheath Contraction Is Associated with
Site-Specific Evolutionary Rates in an AFP System.
Public source: https://github.com/Anthony-Pius-1729/sheath_evolution
Versioned release: https://github.com/Anthony-Pius-1729/sheath_evolution/releases/tag/v1.0.0
Snapshot: 2026-10-14. DOI: not assigned; separate archival-service deposition pending.
The release URL identifies the publication destination; see GitHub for upload status.

Extract into a new empty directory. Paths are relative to the extracted project root.
`ARCHIVE_SHA256.csv` inventories every included original project file and its byte count/hash.
Read `README.md`, `DATA_DICTIONARY.md` and `REPRODUCIBILITY_AUDIT.md` before reuse.
`FROZEN_ANALYSIS.md` is retained unchanged; its state-definition/sensitivity-scope mismatch
with implemented inputs/archived checks is disclosed in the current text. Historical files
are included as decision evidence and may describe earlier stages rather than final sets.

Primary analysis entry points:

- `orthology/afp_sf1_afp2_like.afa`, `orthology/afp_sf1_afp3_like.afa`,
  `orthology/diffocin_primary.afa`: final primary alignments.
- `phylo/*.treefile`, `phylo/*.r4s.nwk`, `phylo/*.log`: saved final topologies and provenance.
- `rates/*.cmd`, `rates/*.res`, `rates/*.unnorm`: exact commands and retained rate outputs.
- `sheath_pilot/qc_output/R_*.csv`: original structural inputs and averaged monomer descriptors.
- `results/sitetable_primary_*.csv`, `results/hypothesis_test_primary.json`: fitted sites/models.
- `results/sensitivity.json`: all eight archived sensitivity checks, including parallel
  variable-site restrictions extending the literal written Afp2-only clause.
- `supplementary_data/TableS*.csv`: full-precision presentation tables and residue mappings.
- `verify_retained_data.py`, `RETAINED_DATA_VALIDATION.json`: saved-file correspondence
  audit covering all 927 recorded structural positions and 223 primary fitted sites.

No analysis was rerun to create this archive. Source scripts, recorded sequence pools,
retrieval/orthology evidence and intermediate outputs are supplied for inspection and
replay from retained inputs. Rate estimates, structure geometry and models are unchanged.
The numerical analyses require scientific Python dependencies and external tools; the
original requirements record omits the statsmodels version used in regression. Do not
assume a current plotting environment is the original fitting environment.

Important completeness limits:

Original mmCIF coordinate bytes/checksums, selected structural-reference keys, candidate
tie scores and state-specific descriptor intermediates are absent. External tool binaries
and full third-party build source are omitted; recorded versions/build notes remain in
`FROZEN_ANALYSIS.md`. Public structure accessions identify deposited sources, but later
coordinate downloads cannot prove original-byte identity. Exact historical structural-copy
reconstruction has not been verified. This archive does not imply that it has been.

Publication figures are included. Manuscript drafts, LaTeX files, submission documents
and the Overleaf ZIP are excluded, even if present locally. Authoring scripts require
their local manuscript inputs; the retained-data verifier works with the code/data archive.
GitHub release hosting is separate from archival-service deposition and DOI registration.
No DOI, complete original numerical environment or verified structural reconstruction is
claimed. Third-party data retain their source terms; this package adds no license grant.
'''
archive=ROOT/'reproducibility_archive.zip'
excluded={'MANUSCRIPT.md','SUPPLEMENTARY.md','MANUSCRIPT_REVIEW.md','MANUSCRIPT_VALIDATION.json','REVIEW_RESPONSE.md','SUBMISSION_CHECKLIST.md','COVER_LETTER_DRAFT.md','LATEX_README.md','manuscript.tex','supplementary.tex','title_page.tex','manuscript_for_review.tex','supplementary_for_review.tex'}
assert not excluded & {str(p.relative_to(ROOT)) for p in files}
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    z.writestr('ARCHIVE_README.md',readme)
    z.writestr('ARCHIVE_SHA256.csv',stream.getvalue())
    for p in sorted(files):z.write(p,str(p.relative_to(ROOT)))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    assert len(z.namelist())==len(set(z.namelist()))
    assert not excluded & set(z.namelist())
    for row in rows:assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256']
try:
    head=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True)
    build_head=head.stdout.strip() if head.returncode==0 else None
except FileNotFoundError:
    build_head=None
summary=dict(archive=archive.name,status='code/data package prepared for versioned GitHub release; DOI deposition pending',version='v1.0.0',repository='https://github.com/Anthony-Pius-1729/sheath_evolution',release_url='https://github.com/Anthony-Pius-1729/sheath_evolution/releases/tag/v1.0.0',build_git_head=build_head,project_files=len(rows),source_bytes=sum(row['bytes'] for row in rows),archive_bytes=archive.stat().st_size,sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),manuscript_documents_included=False,original_atomic_coordinates_included=False,exact_historical_selected_subunits_verified=False)
(ROOT/'REPRODUCIBILITY_ARCHIVE.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
