#!/usr/bin/env python3
"""Describe archived outputs; never retrieve sequences, recompute contacts or fit models."""
from pathlib import Path
import csv, json, hashlib, re, shutil
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'supplementary_data'; OUT.mkdir(exist_ok=True)
def readcsv(path):
    with (ROOT/path).open() as f: return list(csv.DictReader(f))
def writecsv(name, rows):
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def mdtable(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def fasta(path):
    seqs={};name=None
    for line in (ROOT/path).read_text().splitlines():
        if line.startswith('>'): name=line[1:].split()[0];seqs[name]=''
        elif name:seqs[name]+=line.strip()
    return seqs
names={'afp_Afp2':'AFP Afp2','afp_Afp3':'AFP Afp3','diffocin':'Diffocin'}
primary={r['dataset']:r for r in json.loads((ROOT/'results/hypothesis_test_primary.json').read_text())}
struct=[('r2pyocin','F1','R2 pyocin','3J9Q/1','3J9R/1','A–R, X, a, d, g, j, m','0–5, A–R, a–l',385,383,'deposited copies',24,36),('afp','F1','AFP Afp2','6RBN/2','6RC8/1','C','A',350,348,'deposited operators',24,48),('afp','F2','AFP Afp3','6RBN/3','6RC8/2','D','B',358,355,'deposited operators',24,48),('salcis','F1','SalCIS','9R9H/1','9R9I/1','T–Z, a–k','O–Z, a–f',514,495,'deposited copies',18,18),('diffocin','F1','Diffocin','8V3X/1','8V3Y/1','E, F, P, Q, S–X, d, e, g, n, q–s, z','A–Z, a–d',352,353,'deposited copies',18,30)]
s1=[];mapping=[]
for sys,fam,name,pe,pc,ce,cc,me,mc,coords,ne,nc in struct:
    rows=[r for r in readcsv(f'sheath_pilot/qc_output/R_{sys}.csv') if r['family']==fam]
    eligible=[r for r in rows if r['in_domain']=='True'];identity=float(rows[0]['identity'])
    for r in rows:mapping.append(dict(protein=name,**r))
    s1.append(dict(protein=name,extended_PDB_entity=pe,contracted_PDB_entity=pc,extended_author_chains=ce,contracted_author_chains=cc,modeled_E=me,modeled_C=mc,identity=identity,eligible_sites=len(eligible),eligible_min_SEQRES=min(int(r['seqres_pos']) for r in eligible),eligible_max_SEQRES=max(int(r['seqres_pos']) for r in eligible),coordinate_source=coords,assembly_subunits_E=ne,assembly_subunits_C=nc,selected_contact_subunit_key='not archived'))
writecsv('TableS1_structural_references.csv',s1);writecsv('TableS1_residue_mapping.csv',mapping)
s2=[dict(system='R2 pyocin',pool=8035,blast_hits=4593,gated=1266,species_dereplicated=967,decision='excluded: reference lineage ambiguous'),dict(system='SalCIS',pool=7592,blast_hits=4490,gated=805,species_dereplicated=667,decision='excluded: reference locus poorly recovered'),dict(system='Diffocin',pool=8005,blast_hits=930,gated=192,species_dereplicated=165,decision='retained: primary 33 sequences'),dict(system='AFP joint route',pool=8005,blast_hits=5470,gated='not archived as a comparable single-query funnel',species_dereplicated='not applicable to joint route',decision='retained: Afp2 27 / Afp3 31')]
writecsv('TableS2_retrieval_decisions.csv',s2)
alns={'afp_Afp2':'orthology/afp_sf1_afp2_like.afa','afp_Afp3':'orthology/afp_sf1_afp3_like.afa','diffocin':'orthology/diffocin_primary.afa','diffocin_asscreened':'orthology/diffocin_ortholog_set.afa','diffocin_collinear':'orthology/diffocin_collinear_frame.afa'}
# Alternate alignment paths are taken from the archived Rate4Site commands.
for k in ('diffocin_asscreened','diffocin_collinear'):
    alns[k]=re.search(r'-s (\S+)',(ROOT/f'rates/{k}.cmd').read_text())[1]
s3=[]
for k,path in alns.items():
    seqs=fasta(path);lens={len(s) for s in seqs.values()};assert len(lens)==1
    if k in primary:tab=readcsv(f'results/sitetable_primary_{k}.csv')
    else:tab=readcsv('results/sitetable_S1_diffocin_as_screened_diffocin.csv' if k=='diffocin_asscreened' else 'results/sitetable_S2_diffocin_collinear_frame_diffocin.csv')
    p2c={};pos=0
    for col,aa in enumerate(seqs['REF']):
        if aa!='-':pos+=1;p2c[pos]=col
    variable=sum(len({seq[p2c[int(r['seqres_pos'])]] for seq in seqs.values()}-{'-'})>=2 for r in tab)
    model=re.search(r'Best-fit model: (\S+) chosen according to BIC',(ROOT/f'phylo/{k}.log').read_text())
    s3.append(dict(dataset=names.get(k,k.replace('_',' ')),sequences=len(seqs),alignment_columns=next(iter(lens)),reference_length=len(seqs['REF'].replace('-','')),fitted_sites=len(tab),variable_sites=variable,min_fitted_occupancy=min(float(r['occupancy']) for r in tab),IQTREE_model=model[1] if model else 'see archived log',alignment=path))
writecsv('TableS3_alignment_metadata.csv',s3)
s4=[]
for k in ['afp_Afp2','afp_Afp3','diffocin']:
    r=primary[k]
    for p in ['R','rsa','wcn','gnm']:
        s4.append(dict(dataset=names[k],predictor=p,beta=r['beta'][p],SE=r['se'][p],OLS_CI_lo=r['ci'][p][0],OLS_CI_hi=r['ci'][p][1],OLS_p=r['p'][p],VIF=r['vif'][p]))
writecsv('TableS4_full_coefficients.csv',s4)
sens=json.loads((ROOT/'results/sensitivity.json').read_text());s5=[]
for tag,r in sens.items():
    s5.append(dict(analysis=tag,dataset=r['dataset'],sequences=r['n_sequences'],sites=r['sites_fitted'],beta_R=r['beta']['R'],SE_R=r['se']['R'],OLS_p=r['p']['R'],boot_CI_lo=r['block_bootstrap']['ci95'][0],boot_CI_hi=r['block_bootstrap']['ci95'][1],partial_R2=r['partial_r2_R'],full_R2=r['r2'],status='archived parallel check beyond literal variable-site clause' if tag in ['S3_afp_Afp3_variable_only','S3_diffocin_variable_only'] else 'named in written sensitivity plan'))
writecsv('TableS5_sensitivity_results.csv',s5)
commands=[]
for k in alns:
    log=(ROOT/f'phylo/{k}.log').read_text();commands.extend([k,re.search(r'^Command: (.+)$',log,re.M)[1],re.search(r'^Seed: (.+)$',log,re.M)[1],(ROOT/f'rates/{k}.cmd').read_text().strip(),''])
(OUT/'archived_commands.txt').write_text('\n'.join(commands)+'\n')
# Copy complete outputs and selected inputs into a reviewable, compact supplementary bundle.
copy_paths=['results/hypothesis_test_primary.json','results/sensitivity.json','results/partial_residuals.csv','results/R_all_references.csv','orthology/primary_alignment_verification.csv','orthology/afp_sf1_feasibility.csv','orthology/dataset_decision.csv','orthology/afp_check3.json','orthology/r2pyocin_genome_copies.csv','orthology/r2pyocin_within_genome_distance.csv','orthology/salcis_loci.csv','orthology/salcis_ent_loci.csv','orthology/diffocin_all_loci.csv','sheath_pilot/qc_output/structural_qc_report.txt','FROZEN_ANALYSIS.md','hypothesis_test.py','requirements.txt','sheath_pilot/structural_qc.py','sheath_pilot/covariates.py']
copy_paths += list(alns.values())
copy_paths += ['REPRODUCIBILITY_AUDIT.md','structural_provenance_audit.json','DATA_DICTIONARY.md','verify_retained_data.py','RETAINED_DATA_VALIDATION.json']
copy_paths += ['reference_sequences.json','DATASET_DECISION.md','AFP_FEASIBILITY.md','DIFFOCIN_FEASIBILITY.md','R2_FINAL_DECISION.md','sheath_pilot/lattice.py','sheath_pilot/remodeling.py','afp_loci.py','afp_tree.py','locus_architecture.py','filter_hits.py','local_search.py','uniprot_retrieve.py','prelim_tree.py']
copy_paths += [f'retrieval_uniprot/{fam}_provenance.json' for fam in ['r2pyocin','afp_Afp2','afp_Afp3','salcis','diffocin']]
copy_paths += [f'retrieval_uniprot/{fam}_search_provenance.json' for fam in ['r2pyocin','salcis','diffocin']]
copy_paths += [f'sheath_pilot/qc_output/R_{fam}.csv' for fam in ['r2pyocin','afp','salcis','diffocin']]
copy_paths += ['orthology/afp_reference_track.json','orthology/salcis_reference_track.json','orthology/salcis_diffocin_reference_tracks.json']
for k in alns:copy_paths += [f'phylo/{k}.treefile',f'phylo/{k}.log',f'phylo/{k}.r4s.nwk',f'rates/{k}.res',f'rates/{k}.unnorm',f'rates/{k}.log',f'rates/{k}.cmd']
for p in copy_paths:
    src=ROOT/p;assert src.exists(),p
    dst=OUT/'archived'/p;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
manifest={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in copy_paths}
(OUT/'source_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
(OUT/'README.md').write_text('Presentation of saved analysis outputs. Table CSVs preserve full numerical precision. archived/ contains selected original outputs and provenance files, not a complete analysis environment or a public repository. Paths in archived_commands.txt are relative to the original project root. No scientific fits were run to create these files.\n')
intro='''# Supplementary Methods and Results

Anthony Yua Ior, Istanbul Technical University, Maslak 34867, Istanbul, Turkiye.
Correspondence: ior23@itu.edu.tr.

This supplement reports retained outputs and reproducibility details for the three-protein
analysis. It does not introduce new rate estimates, model fits or structural calculations.
The complete-precision tables, residue mappings, commands and selected archived outputs
are supplied in `supplementary_data/`. Code and retained study outputs are hosted in the
public [GitHub repository](https://github.com/Anthony-Pius-1729/sheath_evolution) and
[v1.0.0 code/data release](https://github.com/Anthony-Pius-1729/sheath_evolution/releases/tag/v1.0.0).
Separate archival-service deposition and DOI registration remain pending. Current manuscript
drafts and submission files are excluded from the release package.

## Supplementary Methods

### Structural coordinates, mappings and eligibility

All five references have both experimental states and complete adjacent-ring shells in the
stored structural QC. Canonical indices are zero-based, consecutively assigned to aligned cross-state residue
pairs (gapped segments receive no index); `seqres_pos`
is one-based extended entity sequence position. The residue-mapping CSV preserves author
residue numbers in both states, canonical indices, contact counts, eligibility and descriptor
values for every recorded position. Gaps and unresolved residues are not assigned coordinates.
The complete export contains 927 recorded positions, including 435 eligible positions.
`DATA_DICTIONARY.md` defines its fields. `RETAINED_DATA_VALIDATION.json` records an audit of
this export, contact-count arithmetic and joins for all 223 primary fitted sites against
the structural tables, alignments and Rate4Site outputs. This verifies saved-file
correspondence, not the original atomic contacts; individual contact tuples were not retained.
The two author chain sets are case-sensitive; operator-derived AFP copies are distinct from
the deposited chain names. A selected contact-reference subunit key was not archived.

The implementation prefers non-terminal family subunits and then the candidate with the
largest number of subunit centroids at strictly less than 90 Å. Python `max` retains the
first candidate in input order when scores tie. The original mmCIF bytes, selected keys,
candidate scores and tie counts are absent, so exact historical subunits cannot be verified
from the archived outputs. The accompanying `REPRODUCIBILITY_AUDIT.md` documents the
ordered-input procedure and separates downstream input replay from structural reconstruction. Residue correspondence uses global
BLOSUM62 alignment with opening −11 and extension −1. Modeled amino acids require Cα;
family matching requires at least 90% identity, calculated as the number of identical
aligned pairs divided by the shorter modeled sequence length. Cross-state mismatches remain aligned when
the family meets this criterion. The continuous producer does not implement an engineering
mask or the earlier pilot's 8 Å exclusion. Thus neither an engineering exclusion nor a
specific unlogged central chain is claimed. Direct PDB chain examples in the main table
identify family membership rather than the computational reference copy.

R2, SalCIS and diffocin use deposited copies; AFP applies deposited biological-assembly
operators in both states. No retained reference needs guessed screw-based ring completion.
The producer's generated-contact flag detects ring-completion keys, rather than all
operator-derived coordinates. Its zero value must not be used to deny AFP operator expansion.
Contact identity includes partner family, relative ring, relative azimuth modulo inferred
ring multiplicity and canonical partner residue. Intramonomer pairs are excluded; another
copy of the same family is eligible. Only matched partner positions contribute. Primary
eligibility requires a contact union of at least three at strictly less than 5 Å.

The archived cutoff and domain sweeps in Figure S1 describe the structural metric only.
They were not used to select a rate-association model. Changing the contact cutoff changes
both eligible-site composition and mean remodeling; for example, Afp2 means range from
0.49 to 0.62 over the saved 4–6 Å sweep. These records do not support blanket numerical
invariance to the cutoff.

### Retrieval and positional assignment

UniProtKB release 2026_03, released 2 September 2026, supplied bacterial Pfam family pools.
Reference UniProt identifiers are listed below. Single-query BLAST+ 2.17.0+ searches used
E ≤ 10⁻⁵, `-seg no` and `-max_target_seqs 20000`. Candidate gates were 30–90% identity,
at least 80% query coverage, and a 0.75–1.30 sequence-length ratio. The best bitscore per
species was retained. These filters describe a candidate set and do not establish orthology.
AFP's 5,470 joint prescreen records and 279 reference-subfamily screened members are a
separate route; missing single-query funnel stages are not replaced by inferred counts.

Annotated genomic neighbourhoods extend ±25 kb around a candidate. Flanking CDS homology
was searched against the reference neighbourhood using BLAST with E ≤ 10⁻⁵,
`-max_target_seqs 10` and `-seg no`. Sheath genes were excluded from anchoring evidence.
Offsets are candidate CDS index minus reference slot. A locus requires at least three
non-sheath anchors, at least 60% support for its modal offset, and collinearity ρ ≥ 0.8
when the collinearity statistic is defined. The retained tables include loci with undefined
ρ, so this last condition is not a claim that every anchored row carries a finite statistic.
Direction was normalized to the sheath orientation, and positional assignment used the
modal non-sheath offset. Unanchored or ambiguous positions were left unassigned.

AFP positions 0, 1 and 2 represent Afp2, Afp3 and a third sheath-like gene. Among 42
reference-subfamily loci, 31 were anchored and 35 contained a three-gene run. Across 58
assigned genes in the comparison, 55 agreed with the sequence label. The diagnostic joint
FastTree phylogeny has 145 tips, including 85 labeled tips and 60 unassigned tips; the
three reference tips are included in those labeled counts. Family-wide labels do not form
strictly monophyletic classes. Assignment is supported within the reference subfamily,
not by a claim that the joint tree alone establishes universal positional correspondence.
The final rate alignments and IQ-TREE topologies are separate from that diagnostic tree.

R2's 5,133 candidate-positive proteomes contain one detected copy in 3,341, two in 1,178,
three in 432, four in 99, five in 38, and at least six in 45; the maximum is 20.
These are copy counts conditional on the Pfam retrieval and must not be read as an
unbiased genome-wide absence screen. Among 25 tested within-genome pairs, ten exceeded
the between-genome 90th-percentile distance. The 9,180 between-genome pair distances are
not independent observations. Reference-lineage assignment remained ambiguous even in
the most diverse single-copy candidate set; no R2 rate regression was performed.

SalCIS's family-wide locus assessment contains 55 loci, nine anchored; 28 contain two
sheath-like genes and 27 one. The restricted enterobacterial assessment has 29 loci,
nine anchored, with 19 single-, three double- and seven triple-sheath loci. Neither route
yields an adequate corresponding reference-position sample. Diffocin has 142 loci with
one local sheath gene each, 67 anchored, and 36 exact-offset assignments. Seven of 137
replicons nevertheless contain two candidate loci. Removing both ambiguous *P. hiranonis*
assignments gives the primary 33-sequence dataset, distinct from screened 35 and
collinear-frame 66 sequence sets.

### Alignments, trees and rate-estimation settings

FAMSA through pyfamsa 0.7.0 produced the saved alignments. A separate embedded FAMSA build
version was not recorded. No final column trimming is documented. A reference-row gap
is skipped when mapping alignment columns to one-based SEQRES positions. Site occupancy
is the fraction of rows with a non-gap amino acid. The primary threshold is 80%; the saved
50% checks add only one Afp3 site. Invariant columns are retained in primary fits and
excluded only in the recorded variable-site checks.

IQ-TREE 3.1.4 used ModelFinder BIC selection, 1,000 ultrafast bootstrap replicates,
1,000 SH-aLRT replicates and four threads. Exact commands and logged random seeds are
in `archived_commands.txt`. Both AFP topologies selected WAG+I+G4; primary diffocin
selected LG+F+I+G4, and both alternative diffocin sets selected Q.PFAM+F+G4.
Earlier FastTree 2.2.0 LG+Γ phylogenies serve orthology/feasibility QC; their summary
distances are not substituted for the final IQ-TREE fit.

Rate4Site self-reports version 2.01 although its ConSurf source-directory label is 3.2.
The retained build patched three pre-C++11 stream idioms in error and file-open handling,
without changes to numerical routines: address comparison in `errorMsg.cpp` and
`!stream.is_open()` checks in `someUtil.cpp` and `rate4siteOptions.cpp`. Source/build
provenance is recorded in `FROZEN_ANALYSIS.md`; the build source is retained in the project,
not fully reproduced in this compact manuscript bundle. All runs used `-Mw -ib -k 16`,
a fixed supplied topology (`-t`), and reference row `-a REF`. Rate4Site reoptimized branch
lengths under WAG+Γ. It implements neither LG nor Q.PFAM, so matrix differences from the
selected diffocin topology models are retained and disclosed. The saved IQ-TREE branch
lengths remain in the original `.treefile` outputs.

Increasing normalized Rate4Site scores mean faster relative evolution. The output also
retains unnormalized posterior rates and uncertainty. Invariant sites receive empirical
Bayes estimates rather than a single assigned floor. The JSON field `invariant_sites_in_fit`
counts sites with zero posterior SD, rather than alignment invariance, and is zero in the
primary models. The actual invariant counts are 33, 17 and eight for fitted
Afp2, Afp3 and diffocin sites. Table S3 checks distinct non-gap amino acids in the saved
alignments at the fitted reference positions. The `n_aa` field reports occupancy counts,
not amino-acid diversity.

### Structural covariates and statistical implementation

Shrake–Rupley uses the Biopython defaults (1.4 Å probe radius and 100 surface points per
atom) on an isolated chain; residue ASA is divided by Tien theoretical maximum ASA.
WCN sums inverse squared Cα distances to all other modeled residues in that chain,
including sequence neighbours. Self-pairs and squared distances ≤ 10⁻⁶ Å² are excluded.
GNM uses an adjacency of 0 < d < 7.3 Å, the corresponding Kirchhoff Laplacian and the
diagonal of its numerical pseudoinverse. No physical fluctuation scale is inferred.
Covariates use the representative chain selected during family matching, rather than the
interior subunit selected for assembly contacts. The representative monomer-chain keys
were not exported. The producer computes both states and writes the available-state
arithmetic mean of descriptor values; it does not average atomic coordinates.
Only the averaged `rsa`, `wcn` and `gnm` columns are retained; state-specific intermediate
values cannot be checked from the tables. This differs from the written extended-state
specification. The manuscript uses the archived inputs and does not choose states or refit.

For each protein, the outcome and four predictors are standardized with sample SD
(`ddof=1`). OLS includes an intercept and uses the full available eligible sites at the
specified occupancy. The partial R² compares the full model with RSA/WCN/GNM alone.
Saved standard errors, OLS intervals and p-values assume independent sites; scientific
interpretation instead emphasizes the sequence-block intervals. The intercept estimate
and its SE were not exported. It is centered by construction, but no unarchived intercept
inferential statistic is invented in Table S4.

The bootstrap seed is 20261013. Block identifiers are `(seqres_pos - 1) // 10`.
Only nonempty windows enter the draw pool; a sampled block contributes its eligible rows,
so its row count need not equal ten. Blocks are repeatedly drawn with replacement until
the original site count is met or exceeded; the last sampled block is truncated to that
count. All variables are re-standardized in each replicate and the coefficient is
estimated by least squares. The retained runs used all 10,000 replicates. Percentiles
2.5 and 97.5 yield the intervals. This scheme retains some local sequence dependence,
not necessarily dependence between sequence-distant spatial neighbours. Neither full
coverage under spatial dependence nor uncertainty from structures, sequence selection,
phylogenies or rate posteriors is estimated by this resampling procedure.

The written variable-site sensitivity clause explicitly names Afp2. Parallel Afp3 and
diffocin restrictions are present in saved results but not expressly named in that clause.
Their timing relative to the clause is not recoverable from the retained record, and they
are identified as archived parallel checks. The alternative diffocin sets and 50%
occupancy checks are named in the written plan. No additional regression diagnostics,
measurement-error models or between-system effect tests were introduced during writing.

## Supplementary Results and Tables

### Table S1. Structural references and residue eligibility

'''
t1=mdtable(['Protein','PDB/entity E/C','Modeled E/C','Identity','Eligible','Copies E/C'],[[r['protein'],r['extended_PDB_entity']+' / '+r['contracted_PDB_entity'],f"{r['modeled_E']}/{r['modeled_C']}",f"{r['identity']:.3f}",r['eligible_sites'],f"{r['assembly_subunits_E']}/{r['assembly_subunits_C']}"] for r in s1])
chains='\n\n'.join(f"**{r['protein']}:** author chains E: {r['extended_author_chains']}; C: {r['contracted_author_chains']}. Eligible reference-position range {r['eligible_min_SEQRES']}–{r['eligible_max_SEQRES']} (not every intervening position is eligible). Coordinates: {r['coordinate_source']}." for r in s1)
refs=[]
for fam in ['r2pyocin','afp_Afp2','afp_Afp3','salcis','diffocin']:
    p=json.loads((ROOT/f'retrieval_uniprot/{fam}_provenance.json').read_text());r=p['reference'];refs.append(f"{fam}: {r['uniprot']} ({r['length']} aa), Pfam "+', '.join(r['pfam_signatures']))
t2=mdtable(['Route','Pool','BLAST/prescreen','Gated','Species','Decision'],[[r['system'],r['pool'],r['blast_hits'],r['gated'] if isinstance(r['gated'],int) else 'not comparable',r['species_dereplicated'] if isinstance(r['species_dereplicated'],int) else 'not comparable',r['decision']] for r in s2])
t3=mdtable(['Dataset','Sequences','Columns','Fitted sites','Variable','Minimum occupancy'],[[r['dataset'],r['sequences'],r['alignment_columns'],r['fitted_sites'],r['variable_sites'],f"{r['min_fitted_occupancy']:.3f}"] for r in s3])
t4=mdtable(['Protein','Predictor','β','SE','OLS 95% CI','OLS p','VIF'],[[r['dataset'],r['predictor'],f"{r['beta']:+.4f}",f"{r['SE']:.4f}",f"[{r['OLS_CI_lo']:+.4f}, {r['OLS_CI_hi']:+.4f}]",f"{r['OLS_p']:.4g}",f"{r['VIF']:.3f}"] for r in s4])
short={'S1_diffocin_as_screened':'Diffocin screened','S2_diffocin_collinear_frame':'Diffocin collinear','S3_afp_Afp2_variable_only':'Afp2 variable','S3_afp_Afp3_variable_only':'Afp3 variable*','S3_diffocin_variable_only':'Diffocin variable*','S4_diffocin_occ50':'Diffocin occupancy 50%','S4_afp_Afp2_occ50':'Afp2 occupancy 50%','S4_afp_Afp3_occ50':'Afp3 occupancy 50%'}
t5=mdtable(['Specification','Sequences','Sites','β_R','SE','OLS p','Sequence-block 95% CI','Partial R²','Full R²'],[[short[r['analysis']],r['sequences'],r['sites'],f"{r['beta_R']:+.4f}",f"{r['SE_R']:.4f}",f"{r['OLS_p']:.4g}",f"[{r['boot_CI_lo']:+.4f}, {r['boot_CI_hi']:+.4f}]",f"{r['partial_R2']:.4f}",f"{r['full_R2']:.4f}"] for r in s5])
text=intro+t1+'\n\n'+chains+'''

Full mappings and deposited author chain lists are in the Table S1 CSVs. Selected
interior keys are not archived; all five reported shells are complete and no ring completion
was required. Paired identities are rounded in the stored outputs.

### Table S2. Retrieval counts and orthology decisions

'''+t2+'\n\nReference sequence identifiers: '+ '; '.join(refs)+'''.

Counts describe recorded screening stages rather than final interchangeable ortholog sets.
The AFP prescreen differs from the three single-query routes. Family-wide AFP diagnosis
and reference-subfamily selection are separate. Detailed locus records are retained in the
project; Figures S2–S4 summarize selected archived evidence.

### Table S3. Final alignment occupancy and variability

'''+t3+'''

Fitted-site counts use the retained primary or alternative site tables. Afp3 has 82
structural sites, of which one has 0.742 occupancy and is excluded at 80%. Its 81-site
primary model has 64 variable sites; the 82-site occupancy sensitivity has 65.
All Afp2 structural sites have full occupancy. Primary diffocin structural sites have
minimum occupancy 0.970. The supplied CSV records exact alignment paths and reference lengths.

Taxonomic breadth (effective genera) is 5.49 for Afp2, 5.43 for Afp3, 16.71 for primary
diffocin, 17.71 for screened diffocin and 31.27 for collinear diffocin. AFP spans eight
genera; the largest genus contributes 33.3% and 32.3%, respectively. Primary diffocin
contains 29 species in 20 genera and its largest genus contributes 12.5%.
These are sampling summaries, not numbers of independent evolutionary replicates.

Earlier AFP feasibility trees have total lengths 6.62 and 17.02 substitutions per site,
mean patristic distances 0.914 and 1.937, and maximum distances 1.683 and 5.550. The
historical screened 35-sequence diffocin tree has total length 16.18, mean distance 1.592
and maximum 2.444. These QC summaries are not final IQ-TREE summaries and the 35-sequence
values must not be relabeled as properties of the 33-sequence primary tree. The retained
R2 best single-copy set had 136 sequences, 143 sites, eight invariant sites, effective
genera 22.25 and a maximum QC patristic distance of 4.207, but was excluded on orthology
rather than diversity grounds. No corresponding SalCIS rate dataset was constructed.

### Table S4. Full exported primary regression coefficients

'''+t4+'''

All four non-intercept coefficients are standardized within protein. Independent-site OLS
statistics are descriptive; primary remodeling inference uses the block-bootstrap intervals
in main Table 3. Intercept inferential statistics were not archived. VIFs indicate appreciable
collinearity among conventional descriptors despite relatively low VIF for remodeling.
Full-precision values are preserved in the supplied CSV and original JSON.

### Table S5. Complete archived sensitivity results

'''+t5+'''

An asterisk identifies parallel variable-site checks beyond the literal Afp2-only clause.
The other rows correspond to named sensitivity specifications. These results are unchanged
from the saved JSON; rounded CSVs are not used as the authoritative coefficient source.
Each model retains all three structural covariates and the same sequence-window bootstrap.
No association is detected in diffocin in these specifications. Afp2's variable-site interval
includes zero, whereas Afp3's parallel restriction retains a positive interval.

## Supplementary Figures

**Figure S1.** Contact-definition and mapping QC. (a) Explanatory sheath/contact-set schematic;
the geometry and contact example are hypothetical and do not represent an atomic rendering.
(b) Stored cutoff and union-domain sweeps show eligible-site counts and mean remodeling for
the five retained references. These are structural QC summaries from the saved report,
not new calculations or evolutionary sensitivity fits. File: `figures/FigureS1a.pdf`.

**Figure S2.** R2 copy-number and paralogy evidence. (a) The existing copy-number audit and
within- versus between-genome distance comparison. Counts are conditional on retrieval;
pair distances share sequences and are not independent. (b) Existing retained/discarded
candidate locus examples illustrate reference-order ambiguity. These diagnostic candidate
sets were not used for a primary rate regression. File: `figures/FigureS2a.png`.

**Figure S3.** AFP locus architecture and phylogenetic assignment. (a) Reference-neighbourhood
gene order from the saved annotation, oriented relative to the sheath. Arrows show CDS
orientation and relative order, not nucleotide lengths; sheath slots 0–2 are highlighted.
(b) Existing 145-tip joint FastTree diagnostic phylogeny with positional classes and
reference tips. It is distinct from the final protein-specific IQ-TREE rate topologies.
File: `figures/FigureS3a.pdf`.

**Figure S4.** SalCIS and diffocin positional evidence. Reference gene order and counts of
anchored/unanchored loci are drawn from saved annotation and locus tables. SalCIS shows
family-wide and restricted enterobacterial assessments; diffocin shows the family-wide
assessment. These screening-stage counts do not equal final rate-dataset sizes. Arrow
lengths are schematic. Diffocin's local single-sheath pattern does not rule out multiple
loci per replicon. File: `figures/FigureS4.pdf`.

**Figure S5.** Descriptive remodeling range and relative posterior rate uncertainty at the
primary fitted sites. Posterior SD is divided by the sample SD of point-estimated rates
within each protein. Points and boxes summarize saved site tables. This is an existing
descriptive diagnostic, not a new power calculation, uncertainty-propagation model or
regression weighting scheme. File: `figures/FigureS5.pdf`.
'''
(ROOT/'SUPPLEMENTARY.md').write_text(text)
print('Created supplementary text, five full-precision tables, mappings, commands and selected source archive.')
