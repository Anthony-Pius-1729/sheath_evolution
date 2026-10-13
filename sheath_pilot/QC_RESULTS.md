# Structural QC of the frozen remodeling metric — FINAL

`R_i = 1 - |C_E ∩ C_C| / |C_E ∪ C_C|`, contact identity
`(partner_family, Δring, Δazimuth, partner_residue)`, primary cutoff 5.0 Å,
primary domain `|C_E ∪ C_C| ≥ 3`.

**No evolutionary rates computed. No homologs retrieved.**

## The primary structural set

Held-out families that pass QC with **zero dependence on generated coordinates**:

| system | family | n | mean R | median | sd | interior¹ |
|---|---|---|---|---|---|---|
| R2 pyocin | sheath | 143 | 0.693 | 0.800 | 0.334 | 0.490 |
| AFP | Afp2 | 84 | 0.551 | 0.500 | 0.383 | 0.500 |
| AFP | Afp3 | 82 | 0.541 | 0.500 | 0.391 | 0.476 |
| SalCIS | sheath | 68 | 0.569 | 0.583 | 0.414 | 0.397 |
| diffocin | CD1363 | 58 | 0.294 | 0.236 | 0.305 | 0.603 |

**435 residues, 5 families, 4 systems.** No ring completion was applied to any of
them; every contact rests on deposited coordinates.

¹ fraction with 0.05 < R < 0.95.

**Afp2 and Afp3 are two proteins of the same machine in the same organism.** They
are structurally distinct but not independent observations. The
consistency-across-systems claim therefore rests on **four** systems, not five
families.

### Exploratory (hypothesis-generating — excluded from the test)

| system | family | n | mean R | interior |
|---|---|---|---|---|
| AlgoCIS | sheath | 79 | 0.509 | 0.532 |
| PlPVC1 | F1 | 68 | 0.427 | 0.471 |
| PlPVC1 | F2 | 72 | 0.643 | 0.375 |
| PlPVC1 | F3 | 27 | 0.568 | 0.370 |

Vibrio crashes on a non-finite rise (single extended ring) and is excluded anyway.

## Two held-out systems excluded, both for documented reasons

**CIS_Sc (8BL4/8BKY) — sparse extended deposition.** Its extended sheath has 692
sheath–sheath heavy-atom contacts for a central reference chain against R2
pyocin's 4401, a six-fold deficit: 371 of 539 residues modelled at 3.9 Å with a
~158-residue unmodelled stretch. 86% of domain residues share no contacts across
states and R saturates at 0.951 (sd 0.146, interior 0.137). The deposition cannot
support contact-network comparison.

**PVCpnf (6J0B/6J0C) — only two deposited rings in each state.** Every subunit is
terminal, so completion is unavoidable, and the persistent contacts then rest on
it: median `f_generated` = 1.000, 58% of residues above 0.5, R shifting 0.601 →
0.792 when generated contacts are removed symmetrically (Δ = +0.191, ρ = 0.246).
Neither version is trustworthy — with completion R is generated-dependent, without
it the contact sets are truncated by where the deposited segment stops. Retain
only as a sensitivity system, never in the primary analysis.

## Generated-coordinate audit

| family | residues using generated coords | mean f_gen | Δ mean R | ρ | verdict |
|---|---|---|---|---|---|
| AFP Afp2 | 0 | 0.000 | 0.000 | 1.000 | clean (no completion) |
| AFP Afp3 | 0 | 0.000 | 0.000 | 1.000 | clean (no completion) |
| PlPVC1 F1–F3 | 0 | 0.000 | 0.000 | 1.000 | clean |
| R2 / SalCIS / diffocin | — | — | — | — | no completion applied |
| PVCpnf Pvc2 | 38 / 66 | 0.576 | +0.191 | 0.246 | **material — excluded** |

AFP's apparent dependence in the previous run was an artifact of a bug in the
completion rule, not a property of the data. Contacts need only the rings at r−1
and r+1 — any family, since neighbouring rings in an alternating stack hold a
different protein — so three rings suffice regardless of the same-family repeat.
The rule required `2*repeat+1`, i.e. 5 rings for AFP, which has 4. The spurious
completion then shifted the "interior ring" boundary so `pick_reference` could
select a subunit adjacent to generated rings, manufacturing the dependence.
Corrected, AFP needs no completion and `f_generated = 0` throughout.

## Stability and covariates (unchanged by the above)

**Cutoff.** Mean R moves ≤0.05 across 4.0–6.0 Å in every family. The metric is
insensitive to the contact definition.

**Domain threshold.** Matters more: R2 pyocin 0.73→0.60 and AFP 0.59→0.39 moving
from ≥2 to ≥5. The pre-registered ≥3 sits mid-range.

**Covariates (Spearman, in domain).** GNM flexibility is consistently *negatively*
correlated with R (−0.15 to −0.37): more remodeled sites are **less** flexible.
WCN positive but weaker; RSA inconsistent in sign and mostly null. No |ρ| exceeds
0.37, so R is not reducible to burial or flexibility — but flexibility must stay
in the model, and its sign is opposite to the naive expectation.

**Diffocin.** Mildest contraction (ratio 0.65 against 0.42–0.46) and lowest
remodeling (R = 0.294). Less contraction producing less remodeling is the
dose-response a valid metric should show. Retained in the primary set; its
distinct regime is a secondary structural observation, not grounds for exclusion.

## Frozen before homolog retrieval

- primary hypothesis is **non-directional**: R is associated with site-specific
  evolutionary rate independently of monomer packing, RSA and ENM flexibility.
  Direction is secondary and explanatory.
- per-family models, no pooled regression over 435 residues; spatial block
  bootstrap for uncertainty; the result is consistency of β_R across systems.
- MSA occupancy ≥80% non-gap for a site to enter; sensitivity at 70% and 90%.
- no further modification of R, the domain, the cutoff, or the system set.
