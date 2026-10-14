# Structural reproduction for the three fitted proteins

Updated 2026-10-14. The supplied reconstruction covers AFP Afp2 and Afp3
(6RBN/6RC8) and diffocin (8V3X/8V3Y). R2 pyocin and SalCIS were not reconstructed.
No retained measurements, rates or regression models are replaced.

## Results and scope

Eligible-site sets agree: 84 Afp2, 82 Afp3 and 58 diffocin positions (224 total).
One eligible Afp3 position fails the primary occupancy criterion, leaving 223 fitted
positions. Remodeling, RSA and WCN agree to numerical precision in the initial
reconstruction. Initial GNM outputs differ substantially at every eligible position.
The supplied family-specific outputs with explicit zero-mode exclusion reproduce
retained GNM values with maximum absolute difference below 1.7e-13.

`verification_report.json` is the current saved-evidence validation report.
`reconstruction_comparison.csv` records the initial comparison, including the GNM
failure. `reconstruction_verdict.json` is the original explanatory report; its
historical numerical-backend interpretation is not established by these outputs.

A Kirchhoff matrix has a constant null mode for each connected component.
`covariates_zeromode_fix.py` excludes eigenvalues at or below 1e-10 times the
maximum before forming the pseudoinverse diagonal. The directly supported finding
is that this treatment reproduces the retained values; the historical LAPACK build
and its handling of small eigenvalues were not recovered. No claim of bitwise
identity across platforms is made. NumPy documents a default relative cutoff of
1e-15 for `pinv` when neither `rcond` nor `rtol` is supplied; the earlier report's
statement of max(M,N)*eps as that default is incorrect. See the
[official NumPy documentation](https://numpy.org/doc/stable/reference/generated/numpy.linalg.pinv.html).

## Verify the saved evidence without reconstructing structures

From the extracted project root:

```sh
python3 reproducibility_check/verify_package.py
```

This uses the Python standard library, verifies the package source/input hashes,
and compares the supplied CSVs against retained measurements. It writes only the
derived verification report. Original producer checksums match the historical
audit. No structure, contact, rate or model is recalculated by this command.

## Optional reconstruction command for independent reuse

```sh
python3 reproducibility_check/run_check.py
```

This requires the scientific dependencies recorded in the project requirements.
It checks the bundled coordinate/source hashes, runs only AFP and diffocin, applies
explicit zero-mode exclusion by default, records compact selected-subunit identifiers
and current Python/library versions, and writes a new directory under `reruns/`.
The preserved producer source is unchanged; the alternative GNM function is injected
only at runtime. No evolutionary pipeline is invoked. The updated wrapper was
inspected and syntax-checked during packaging; it was not executed to repeat the
closed structural reconstruction.

`--gnm-mode historical` instead uses the unmodified pseudoinverse implementation
and may fail the retained-value comparison in the current numerical environment.
`run_check_historical.py` preserves the initial runner verbatim for provenance; it
is not the supported package entry point and depends on the earlier manuscript folder.

## Which saved outputs to inspect

- `qc_output/R_afp.csv` and `qc_output/R_diffocin.csv`: saved reconstructed tables;
  their GNM columns agree with the explicitly corrected family-specific outputs.
- `qc_output_zeromode/R_afp_F1.csv`, `R_afp_F2.csv` and `R_diffocin_F1.csv`: the
  corrected eligible-site outputs used in the saved-evidence verification.
- The combined `R_afp.csv` and `R_diffocin.csv` under `qc_output_zeromode/` also
  agree with those corrected outputs. The initial failing GNM calculation is
  preserved in the comparison summary; a separate per-site copy of that failing
  output is not supplied. The verifier checks the saved tables and distinguishes
  them from the reported initial failure.
- `reconstructed_references.json`: compact E/C subunit keys, ring and azimuth.
- `recovered_identifiers.json`: original verbose selection records.
- `coordinate_inputs.json` and `structures/`: the exact reconstruction inputs,
  source URLs, sizes and SHA-256 hashes. These are newly retrieved coordinates,
  not verified historical download bytes.
- `FROZEN_ANALYSIS_original.md`: byte-identical original specification. The root
  specification appends dated retrospective notes without rewriting that text.
- `package_manifest.json`: checksums for the reconstruction package and its
  required root inputs; excludes itself and the derived verification report.

Reconstructed identifiers demonstrate a selection yielding agreement with retained
measurements. They do not recover the unlogged historical execution trail, candidate
tie scores, representative monomer identities or original contact tuples.
