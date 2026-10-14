# Reproducibility and structural-reference audit

This cumulative audit distinguishes the original retained files from the supplied bounded
reconstruction for AFP Afp2, Afp3 and diffocin. The current packaging review verified saved
outputs and source/input hashes; it did not repeat reconstruction or evolutionary analysis.
Retained scientific measurements and models remain unchanged. The original frozen text is
preserved verbatim, followed by a dated retrospective deviation note. Machine-readable findings and source-file
SHA-256 checksums are in `structural_provenance_audit.json`.

## What the retained files can establish

Saved structural CSVs contain canonical/author/SEQRES mappings, contact-count summaries,
eligibility, remodeling and the fitted averaged monomer descriptors. The structural report
records assembly expansion, geometry and complete adjacent-ring shells. Rate inputs and
outputs, final alignments/topologies, model outputs and exact commands are retained.
These permit verification of the supplied site inputs and downstream numerical reporting.
They do not identify the historical coordinate bytes or exact selected subunits.

The independent saved-file audit in `verify_retained_data.py` checks all 927 recorded
positions from the five retained references, including 435 eligible positions. It verifies
the complete supplementary mapping export, contact-count arithmetic, metric eligibility
and values, and joins for all 223 primary fitted positions against structural inputs,
alignment reference/occupancy/variability and saved Rate4Site scores/posterior SDs.
`RETAINED_DATA_VALIDATION.json` records the checks and source hashes; `DATA_DICTIONARY.md`
defines the fields. Individual partner-contact tuples were not serialized. These checks
establish internal consistency of the retained measurements, not their atomic-coordinate
reconstruction or the validity of the selected historical subunits.

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

## Reconstruction agreement and historical provenance

The supplied reconstruction reproduces the retained structural inputs for the three fitted
proteins. All 224 eligible positions are recovered (84 Afp2, 82 Afp3 and 58 diffocin);
223 pass the primary occupancy requirement. Remodeling, RSA and WCN agree to numerical
precision. Initial GNM discrepancies are recorded in `reconstruction_comparison.csv`.
Explicit zero-mode exclusion reproduces retained GNM values with maximum absolute
difference below 1.7e-13. The currently saved reconstructed CSVs contain those agreeing
values; a separate per-site copy of the initial failing output is not present.

`reproducibility_check/verify_package.py` checks the supplied tables, original producer
checksums, coordinate hashes and original specification. The reconstruction input mmCIF
files are packaged, with source URLs and hashes in `coordinate_inputs.json`. The supported
optional reconstruction runner uses those inputs and the separate zero-mode function,
writes new outputs and records selected keys and numerical versions. Its wrapper was
syntax-checked during packaging; the reconstruction itself was not repeated.

Reconstructed contact-reference subunits (extended / contracted):

- Afp2: `C#1` / `A#2`.
- Afp3: `D#0` / `B#0`.
- Diffocin: `g` / `B`.

The compact identifier record preserves ring/azimuth values from the supplied selection
records. Agreement establishes a reconstruction consistent with retained measurements;
it does not recover the unlogged historical execution trail, candidate scores, monomer
chain identities or historical download bytes. The four packaged coordinate files are
new retrievals, whose byte identity with historical downloads cannot be established.
R2 pyocin and SalCIS were not reconstructed. No historical measurement was replaced.

The earlier report's explanation of a particular LAPACK build is unverified. The supported
claim is that explicit zero-mode exclusion removes the discrepancy. The reconstruction
README corrects the earlier statement about NumPy's default pseudoinverse tolerance and
distinguishes original reported explanations from the evidence actually verified.

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

`reproducibility_archive.zip` supplies code, recorded sequence inputs, processed alignments,
trees, rates, site tables, models, decision evidence, figures and documentation. Version
v1.0.1 adds the bounded reconstruction package and its coordinate files, scripts, identifiers
and checksums. The published v1.0.0 remains available unchanged. The source repository is
public; separate archival-service DOI deposition remains pending. Manuscript drafts and
submission documents are excluded from the code/data release. No complete historical
execution environment or original coordinate-download identity is claimed.
