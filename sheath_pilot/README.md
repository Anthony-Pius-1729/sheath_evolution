# Contractile sheath interface-persistence — structural pilot

Feasibility gate only. **No evolutionary rates are estimated here.** The pilot asks
whether matched extended/contracted structures yield interface-persistence classes
populated and distinguishable enough to support the planned analysis.

```bash
python3 pilot.py        # downloads structures on first run, writes output/
```

Needs `numpy scipy pandas biopython`. Structures are fetched from RCSB and cached
in `structures/`.

## Pre-registered parameters

Fixed before any counts were computed; not revised afterwards.

| parameter | value |
|---|---|
| contact definition | heavy-atom distance < 5.0 Å |
| cutoff sweep | 4.0, 4.5, 5.0, 5.5, 6.0 Å |
| minimum class size | ≥ 20 residues in the smallest interface class |
| partner identity | `(canonical entity, screw index n, C6 index k)` |
| covariates | monomer RSA, monomer WCN, monomer GNM fluctuation |
| covariate state | computed in both states; primary = mean |
| contact scope | cross-state comparable entities only |
| engineering shell | 8.0 Å around internal indels |

## Systems

| | extended | contracted | twist/rise ext | twist/rise con |
|---|---|---|---|---|
| *Algoriphagus machipongonensis* eCIS (primary) | 7AE0 (2.4 Å) | 7AEK (2.5 Å) | 20.54° / 40.80 Å | 32.20° / 18.04 Å |
| *Vibrio cholerae* T6SS (secondary) | 5MXN (3.7 Å) | 3J9G (3.5 Å) | 23.50° / 37.80 Å | 29.40° / 21.80 Å |

Both are C6 plus a screw. Deposited helical parameters are used for lattice
indexing and verified against the coordinates; `lattice.py` halts rather than
guessing if they do not reproduce the observed chain placements.

## Three faults found and fixed during development

Each produced a plausible-looking but wrong result, so they are recorded here.

1. **Arbitrary axis sign.** The helical axis recovered from superposition has no
   canonical direction. An inconsistent choice between states inverts every screw
   index, converting genuinely persistent partners into apparently changed ones.
   Symptom: Algoriphagus `persistent_same_partner` = 1 against
   `persistent_changed_partner` = 89. The sign is now resolved once per assembly
   against the deposited (positive twist, positive rise) convention and imposed on
   indexing *and* symmetry expansion. Corrected values: 86 and 2.

2. **Deposited entity IDs are not cross-state comparable.** VipB is entity 2 in
   5MXN but VipA is entity 1 in 3J9G, so a raw ID in a contact label makes a
   protomer's own partner look like a different molecule in the other state.
   Descriptions are no better — 5MXN calls both VipA and VipB "Type VI secretion
   protein". Entities are now matched across states by sequence identity (≥ 0.90)
   and relabelled to canonical keys. VipB went from 4/177 to 78/100.

3. **Per-entity sublattices.** The inner tube sits at a different radius and
   azimuth from the sheath; indexing it against the sheath reference flagged all
   18 tube chains as tolerance failures. Indices are now computed per entity
   relative to that entity's own reference subunit.

## Result

Pre-registered four-way gate, native coordinates only, engineering shell applied:

| family | same | changed | ext-only | con-only | smallest | verdict |
|---|---|---|---|---|---|---|
| Algoriphagus sheath (692 aa) | 86 | 2 | 45 | 24 | **2** | FAIL |
| Vibrio VipB (473 aa) | 73 | 17 | 4 | 50 | **4** | FAIL |
| Vibrio VipA (159 aa) | 56 | 7 | 4 | 3 | **3** | FAIL |

**0 of 3 families pass.**

**Authoritative source: `output/summary.csv` and `output/qc_report.txt`**, both
regenerated 2026-10-12 by `regenerate_native_summary.py` under the rule above.
Until that date those two files held the original `pilot.py` run, which counted
symmetry-dependent residues and therefore reported VipB as PASS and 1/3 passing;
that run is preserved unmodified in
`output/archive_symmetry_inclusive_2026-10-10/` and must not be cited as the
verdict. The residue tables were not touched by the regeneration — they carry
`uses_symmetry_generated` and `engineering_excluded` per residue, so either
counting rule is derivable from them, and the script asserts that it reproduces
both the archived counts and the native counts before writing anything.

The VipA row reads 56/7/4/3 rather than the 57/7/4/4 reported before 2026-10-12.
The difference is the engineering shell: 4 of VipA's residues lie within 8.0 Å of
an internal indel, and `pilot.py` had always excluded them from its own counts
while the native table quoted here had not. The shell is now applied in both
places. Verdict unchanged; FAIL either way.

Two findings behind that verdict:

- **`persistent_changed_partner` is near-empty because the lattice keeps its
  connectivity through contraction.** Subunits retain the same `(n, k)` neighbours
  and slide and tilt rather than swapping partners. The class is well defined; it
  is simply not populated. This is a property of the mechanism, not of the data.

- **Vibrio's apparent pass was an artifact.** 5MXN contains a single ring, so every
  inter-ring contact in the extended state must be generated from the deposited
  helical operator. 107 of 429 VipB residues depend on generated neighbours, and
  83 of its 100 `persistent_changed_partner` residues do. With generated contacts
  included VipB scores 20 and passes; on native coordinates it scores 4 and fails.
  Any Vibrio result resting on inter-ring extended-state contacts is reporting the
  symmetry operator, not the structure.

Post-hoc three-way collapse (persistent / state-specific / non-interface), native
only — **reported for design purposes only, chosen after seeing the counts**:

| family | persistent | state-specific | non-interface | smallest |
|---|---|---|---|---|
| Algoriphagus sheath | 88 | 69 | 354 | 69 |
| Vibrio VipB | 90 | 54 | 178 | 54 |
| Vibrio VipA | 63 | 7 | 25 | 7 |

## Flexibility confound — confirmed, and large

GNM fluctuation tracks persistence class in the predicted direction in every
family. Vibrio VipA: `persistent_changed_partner` mean 3.70 against
`persistent_same_partner` 0.81. Monomer RSA and WCN separate every interface class
from `non_interface` (RSA ≈ 0.18–0.22 non-interface against 0.33–0.64 interface).

An analysis that omitted these covariates would recover a strong "state-specific
sites evolve faster" result driven by burial and flexibility alone.

## Output

- `output/summary.csv` — **authoritative**, one row per family, native rule.
  Carries `counting_rule`, `symmetry_dependent_excluded` and the superseded
  `archived_verdict` alongside each verdict.
- `output/qc_report.txt` — **authoritative** verdict report, native rule, with
  input checksums and a per-family archived/native/native+engineering comparison.
- `output/archive_symmetry_inclusive_2026-10-10/` — the original `pilot.py` run,
  unmodified, plus `ARCHIVE_README.md` explaining why it is superseded. It is
  also the only record of the per-state structural diagnostics, the
  contact-cutoff sweep and the per-class covariate means: those need the
  coordinates to recompute, were not regenerated, and were computed **without**
  the native filter. Do not quote them as native-rule numbers.
- `output/residues_<system>_<family>.csv` — residue-level annotation: persistence
  class, partner labels in both states, `uses_symmetry_generated`,
  `engineering_excluded`, `dist_to_engineering`, and RSA/WCN/GNM per state and mean

## Files

- `pilot.py` — driver, classification, reporting
- `lattice.py` — helical axis, screw decomposition, lattice indexing, expansion
- `covariates.py` — monomer RSA, WCN, GNM
- `regenerate_native_summary.py` — rewrites `output/summary.csv` and
  `output/qc_report.txt` from the residue tables under the native rule. Needs
  pandas only (no coordinates, no structure download). `--dry-run` prints the
  report without writing. Aborts rather than writing if the recomputed counts
  disagree with either the archived run or the independently verified native
  counts.

## Manuscript status

Supplementary / methodological-history material, not a primary experimental
section. This pilot's failure is what motivated the continuous remodeling
formulation; the primary manuscript is the frozen `R_i` analysis. If any of it
is reported, its exploratory status and the 0/3 gate failure are disclosed with
it.

## Scope — this directory is not the R_i analysis

Everything above is the **historical categorical pilot**: four-way persistence
classes against a pre-registered minimum class size of 20, a feasibility gate
that failed. The **frozen continuous remodeling metric** `R_i` is a separate,
later analysis — `remodeling.py` and `structural_qc.py`, results in `qc_output/`,
documented in `QC_RESULTS.md`, with its own system set and its own held-out
design. The two share a codebase and nothing else. Keep their numbers, their
directories and their verdicts apart; the categorical FAIL says nothing about
`R_i`, and `QC_RESULTS.md` is not a report on the gate above.
