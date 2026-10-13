# Reproducibility and structural-reference audit

This audit reads existing files and source. It did not download coordinates, rebuild
assemblies, recompute contacts, estimate rates or refit models. Original scientific outputs
and the written specification remain unchanged. Machine-readable findings and source-file
SHA-256 checksums are in `structural_provenance_audit.json`.

## What the retained files can establish

Saved structural CSVs contain canonical/author/SEQRES mappings, contact-count summaries,
eligibility, remodeling and the fitted averaged monomer descriptors. The structural report
records assembly expansion, geometry and complete adjacent-ring shells. Rate inputs and
outputs, final alignments/topologies, model outputs and exact commands are retained.
These permit verification of the supplied site inputs and downstream numerical reporting.
They do not identify the historical coordinate bytes or exact selected subunits.

## Ordered structural procedure in the retained implementation

1. The producer expects each original PDB accession's mmCIF under `sheath_pilot/structures/`.
   A missing file would be fetched from the public download endpoint. No copy of these
   original coordinate files, download revision or byte checksum is now present. Accession
   identity alone does not prove byte identity with a later download.
2. `load` parses the first model. Modeled standard amino acids require Cα, and representative
   chains are taken in parser/model order, with the first encountered chain retained for
   each entity. Cross-state entities are matched using the implemented sequence identity
   threshold. Exact parser order is consequently an input to the procedure.
3. `build_subunits` retains modeled-family chains in model order. Where expansion is needed,
   operators are iterated in mmCIF operator-list order, followed by base-chain order. Expanded
   keys use the author chain plus `#` and the zero-based operator index. These are distinct
   from `~` keys for guessed ring completion, which no primary reference required.
4. The coordinate-derived lattice assigns rings and azimuths. Axis orientation is reconciled
   between states, and the source repeats selection after this orientation step. Centroids
   are calculated from modeled Cα atoms, not from all heavy atoms.
5. For a family, `pick_reference` preserves the order of its candidates. It uses candidates
   strictly inside the minimum/maximum ring where any are available; otherwise it falls
   back to all family candidates. Each candidate's score counts all assembly centroids at
   distance strictly less than 90 Å, including its own centroid. Python `max` returns the
   first highest-scoring candidate when counts tie. There is no separate random tie break.
6. Contact-reference selection is repeated independently for each family/state after the
   assembly/lattice operations. The selected reference must have adjacent rings on both
   sides to pass the recorded shell check. Contact sets use that subunit's coordinates.
7. Monomer descriptors are calculated on the family representative chains from step 2,
   rather than the selected interior contact subunits. Descriptor values, not coordinates,
   are averaged across available states. Neither these representative-chain keys nor
   state-specific intermediate descriptor arrays were exported.

The rule is deterministic conditional on ordered coordinate input and compatible numerical
behavior. The retained code alone does not certify the original selected keys. Candidate
scores, tied candidates, lattice labels and selected keys were not serialized. The CSV
mappings and shell checks do not supply the missing identity. No example chain in main
Table 1 is promoted to a verified computational reference.

## Conclusion of the requested provenance check

Exact historical structural-reference reconstruction is **not established by this archive**.
The original coordinate bytes and selected keys are absent. Current PDB downloads could
support an isolated prospective reconstruction, but that would re-execute the structural
pipeline and would not by itself prove original-byte identity. It was outside this revision's
closed-analysis scope and was not performed. The manuscript retains the limitation.

A future reconstruction would need to preserve input mmCIF checksums/revisions, parser and
numerical versions, representative-chain keys, operator order, subunit keys/ring/azimuths,
candidate scores/tie order and the final reference keys in both states. Comparison should
then be made against archived site outputs in an isolated directory. This describes needed
provenance; it is not a report that such a verification has already occurred.

## Statistical and environment boundaries

The fitted estimand conditions on the archived monomer descriptors averaged across states.
It is not an estimate for extended-state-only adjustment. The implementation departure from
that written definition cannot be settled by editorial wording or by assuming the averages
are interchangeable with extended-state values.

Sequence-block intervals resample fixed ten-reference-position windows and restandardize
within each replicate. They condition on selected structures, alignments/orthology, fitted
sites and rate point estimates. They do not resample spatial contact clusters, posterior
rates, sequences or trees; coverage under remaining spatial dependence is unassessed.

`requirements.txt` is an original captured environment record, not a verified complete
lockfile. In particular, the regression code imports statsmodels but its original version
is not recorded there. A current plotting environment must not be relabeled as the original
analysis environment. The retained IQ-TREE logs and Rate4Site commands identify their actual
versions/settings; the Rate4Site directory label and self-reported version differ as
already documented. External tool binaries/build sources are not part of the deposition ZIP.

## Archive and deposition status

`reproducibility_archive.zip` packages local code, recorded sequence inputs, processed
alignments/trees/rates/site tables, decision evidence and documentation with a SHA-256
inventory. It is prepared for deposition; it has no public URL or DOI yet. Original atomic
coordinates and unlogged intermediate provenance are absent and are not replaced by guessed
values. The archive supports inspection and downstream replay from retained site inputs,
with the limits above. Public deposition and complete structural reconstruction are separate
remaining tasks.
