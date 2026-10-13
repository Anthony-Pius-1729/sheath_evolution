# ARCHIVE — original symmetry-inclusive pilot outputs (superseded)

These are the **original, unmodified** verdict artifacts of the categorical
interface-persistence pilot, as written by `pilot.py` on **2026-10-10 23:20**.
They are retained as the evidence of what the pilot actually computed on that
run. **Do not cite them as the pilot's verdict.**

| file | sha256 |
|---|---|
| `summary.csv` | `9551011e2f01cd35f8fe0bda1f09f2bb901fa0562680574ecb0e8480bdaf3a29` |
| `qc_report.txt` | `5efa7b44a337b605ae3544c6db1527ae830d8fe709f9fdff0fe05e12fdc39cef` |

## Why they are superseded

`pilot.py` counted every mapped residue, including residues whose interface
classification depends on contacts to **symmetry-generated** neighbours. 5MXN
(*Vibrio* extended state) contains a single ring, so every inter-ring contact in
that state has to be generated from the deposited helical operator. The counts
below therefore report the symmetry operator rather than the deposition:

- `qc_report.txt` line "protein families passing: **1 / 3**"
- `summary.csv` marks *Vibrio* VipB (`P1_Type VI secretion protein_473aa`)
  **PASS**, smallest interface class 20.

The governing rule for this project is **native coordinates only**: a residue
whose classification depends on generated contacts does not enter the count. On
that rule the pilot fails in **0 of 3** families — see
`../qc_report.txt` and `../summary.csv`, regenerated 2026-10-12 by
`../../regenerate_native_summary.py`.

## What is *not* in this directory, and why

`residues_*.csv` are **not** archived here because they are rule-agnostic: they
carry `uses_symmetry_generated` and `engineering_excluded` per residue, so both
the symmetry-inclusive and the native-only counts are derivable from them. The
residue tables in `../` are the same files this run produced, unmodified
(verified by checksum in the regenerated report).

`qc_report.txt` here is also the **only** record of the per-state structural
diagnostics, the contact-cutoff sweep and the per-class covariate means. Those
require the coordinates to recompute and were not regenerated; where the native
rule would change them, the regenerated report says so rather than restating
them.

## Scope

This is the **historical categorical pilot** (four-way persistence classes,
pre-registered minimum class size 20). It is **not** the frozen continuous
remodeling metric `R_i`, which lives in `../../qc_output/` and is documented in
`../../QC_RESULTS.md`. The two analyses are separate; neither supersedes the
other.
