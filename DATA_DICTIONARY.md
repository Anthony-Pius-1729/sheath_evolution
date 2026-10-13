# Retained residue measurements

The tables preserve measurements from the closed analysis. Inspection of these values
does not require downloading coordinates or rerunning a structural producer. The audit
command below checks their completeness, count arithmetic, and correspondence to the
saved alignments and rate outputs; it does not validate the original atomic contacts.

```sh
python3 verify_retained_data.py
```

The command writes `RETAINED_DATA_VALIDATION.json`. It uses the Python standard library
and reads the saved scientific inputs and outputs without modifying them.

## Tables and identifiers

| File | Contents |
|---|---|
| `sheath_pilot/qc_output/R_r2pyocin.csv` | 257 recorded R2 pyocin positions; 143 eligible. |
| `sheath_pilot/qc_output/R_afp.csv` | Afp2 (`F1`): 184 recorded, 84 eligible; Afp3 (`F2`): 190 recorded, 82 eligible. |
| `sheath_pilot/qc_output/R_salcis.csv` | 168 recorded SalCIS positions; 68 eligible. |
| `sheath_pilot/qc_output/R_diffocin.csv` | 128 recorded diffocin positions; 58 eligible. |
| `supplementary_data/TableS1_residue_mapping.csv` | Complete export of these 927 recorded positions, including 435 eligible positions. Adds a readable `protein` label. |
| `results/sitetable_primary_*.csv` | Regression inputs: Afp2 84, Afp3 81, diffocin 58 positions. One eligible Afp3 position fails the 80% occupancy requirement. |
| `results/hypothesis_test_primary.json` | Full-precision model results and sequence-block intervals. |
| `results/sensitivity.json` | All eight saved sensitivity specifications. |

Join structural records by `(system, family, canon)`, and use `(protein, canon)` in the
combined mapping. `canon` is zero-based and identifies an aligned cross-state residue
pair, not a deposited chain. `seqres_pos` is the one-based extended entity reference
position used to join the Rate4Site output and ungapped `REF` alignment sequence.
`resid_ext` and `resid_con` are author residue numbers. The selected subunit and
representative monomer-chain identities are absent; these numbers cannot recover them.

## Structural columns

| Column | Definition |
|---|---|
| `system`, `family` | Producer system and family identifiers. Family labels are local to a system. |
| `n_ext`, `n_con` | Numbers of distinct intermolecular contact identities in each state. |
| `n_inter`, `n_union` | Counts in the intersection and union; union equals extended + contracted − intersection. |
| `in_domain` | Whether the union contains at least three identities. |
| `R` | One minus intersection/union for eligible positions. Empty values are outside the metric domain, not zero remodeling. |
| `rsa` | Monomer solvent accessibility divided by the residue-specific theoretical maximum. |
| `wcn` | Sum of inverse squared Cα distances to other modeled residues of the same monomer, including sequence neighbours. |
| `gnm` | Diagonal of the pseudoinverse of the monomer Kirchhoff matrix using a 7.3 Å Cα cutoff. |
| `identity` | Cross-state modeled family-sequence identity. |
| `contraction_ratio` | Recorded sheath lattice contraction ratio. |

The three descriptor columns are arithmetic means of available extended- and
contracted-state values calculated on isolated monomers. They are not extended-state-only
values. Individual state arrays were not exported. Contact identities were defined at
heavy-atom distance strictly less than 5 Å and encode partner family, relative ring,
relative azimuth and canonical partner residue. The individual contact tuples are absent;
the saved counts allow arithmetic checks of the reported metric, not independent
reconstruction of atomic contacts.

## Rate and alignment columns

| Column | Definition |
|---|---|
| `aa` | Reference amino acid at the mapped position. |
| `rate` | Normalized Rate4Site posterior score; larger values denote faster relative evolution. |
| `rate_std` | Saved posterior SD, not the regression standard error. |
| `n_aa` | Number of aligned sequences containing a non-gap amino acid at the position. This is an occupancy count, not amino-acid diversity. |
| `n_seq` | Number of sequences in the dataset, including `REF`. |
| `occupancy` | `n_aa / n_seq`. |

Variable sites have at least two distinct non-gap amino acids in the final alignment.
The historical JSON field `invariant_sites_in_fit` counts zero posterior-SD entries;
use the alignment-derived variable counts rather than treating that field as a direct
measure of invariant columns.

## Boundaries

Original mmCIF bytes/checksums, selected subunit keys, tie-candidate scores and monomer
chain keys were not retained. Exact historical structural reconstruction remains
unverified. `REPRODUCIBILITY_AUDIT.md` specifies the ordered selection algorithm and
these limits. The bootstrap resamples ten-reference-position sequence windows;
coverage under three-dimensional structural dependence has not been assessed.
