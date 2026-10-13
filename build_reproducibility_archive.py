#!/usr/bin/env python3
"""Package existing local code/data for deposition, without running scientific analysis."""
from pathlib import Path
import csv,hashlib,io,json,zipfile
ROOT=Path(__file__).resolve().parent
EXT={'.py','.csv','.json','.nwk','.afa','.fasta','.fa','.faa','.tsv','.txt','.log','.treefile','.res','.unnorm','.cmd','.stdout','.iqtree','.contree','.bionj','.png','.pdf','.gff','.gff3','.gb','.gbk','.stderr'}
roots=['orthology','phylo','rates','results','retrieval_uniprot','filtered','synteny','evolution_prelim']
files=set()
for folder in roots:
    base=ROOT/folder
    if base.exists():
        files.update(p for p in base.rglob('*') if p.is_file() and p.suffix in EXT)
files.update(p for p in ROOT.glob('*.py') if p.is_file())
files.update(p for p in (ROOT/'sheath_pilot').glob('*.py') if p.is_file())
for folder in ['sheath_pilot/qc_output','sheath_pilot/output']:
    files.update(p for p in (ROOT/folder).rglob('*') if p.is_file() and p.suffix in EXT)
files.update(p for p in (ROOT/'supplementary_data').glob('*') if p.is_file())
for name in ['requirements.txt','reference_sequences.json','pair_meta.json','FROZEN_ANALYSIS.md','DATASET_DECISION.md','RESULTS.md','AFP_FEASIBILITY.md','DIFFOCIN_FEASIBILITY.md','R2_FINAL_DECISION.md','R2_EVOLUTIONARY_FEASIBILITY.md','MANUSCRIPT.md','SUPPLEMENTARY.md','MANUSCRIPT_REVIEW.md','MANUSCRIPT_VALIDATION.json','REVIEW_RESPONSE.md','REPRODUCIBILITY_AUDIT.md','structural_provenance_audit.json','SUBMISSION_CHECKLIST.md','COVER_LETTER_DRAFT.md','LATEX_README.md','sheath_pilot/QC_RESULTS.md']:
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
readme='''# Sheath evolution: local code and data archive prepared for deposition

Study: Contact-Network Remodeling During Bacterial Sheath Contraction Is Associated with
Site-Specific Evolutionary Rates in an AFP System.
Author: Anthony Yua Ior, Istanbul Technical University, Maslak 34867, Istanbul, Turkiye.
Contact: ior23@itu.edu.tr. Snapshot: 2026-10-13. Public URL/DOI: not assigned.

Extract into a new empty directory. Paths are relative to the extracted project root.
`ARCHIVE_SHA256.csv` inventories every included original project file and its byte count/hash.
Read `MANUSCRIPT.md`, `SUPPLEMENTARY.md` and `REPRODUCIBILITY_AUDIT.md` before reuse.
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

The manuscript LaTeX and publication figure files are delivered separately in
`manuscript_latex.zip`. This archive is local only and is not itself a persistent repository
or a complete verified computational environment. Select repository metadata/licensing and
obtain a release URL/DOI before public deposition and final manuscript submission.
'''
archive=ROOT/'reproducibility_archive.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    z.writestr('README.md',readme)
    z.writestr('ARCHIVE_SHA256.csv',stream.getvalue())
    for p in sorted(files):z.write(p,str(p.relative_to(ROOT)))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for row in rows:assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256']
summary=dict(archive=archive.name,status='local archive prepared for deposition; no public URL/DOI',project_files=len(rows),source_bytes=sum(row['bytes'] for row in rows),archive_bytes=archive.stat().st_size,sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),original_atomic_coordinates_included=False,exact_historical_selected_subunits_verified=False)
(ROOT/'REPRODUCIBILITY_ARCHIVE.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
