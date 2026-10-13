# Held-out two-state contractile sheath systems — verified

Every accession below was checked against the RCSB API (exists, state in title,
resolution) and, where the deposition contains ≥2 rings, the extended/contracted
assignment was **independently validated from coordinates** by measuring ring
spacing. Titles were not trusted: the PlPVC1 search showed extended sheath
subunits living inside entries titled "baseplate" and "cap", and 8BKY carries
helical metadata that contradicts its own coordinates.

Ring spacing is the state check. Extended sheaths are strikingly uniform at
38.3–40.8 Å across every system examined; contracted fall to 16–26 Å.

## Tier A — clean, held-out, native multi-ring in both states

100% cross-state sequence identity, zero substitutions, no engineered trap.

| system | organism | extended | contracted | spacing E→C | ratio | sheath proteins | usable res |
|---|---|---|---|---|---|---|---|
| PVCpnf | *Photorhabdus asymbiotica* | 6J0B (2.9 Å) | 6J0C (3.7 Å) | 39.4 → 17.0 Å | 0.43 | 1 (Pvc2, 355 aa) | 352 |
| R2 pyocin | *Pseudomonas aeruginosa* | 3J9Q (3.5 Å) | 3J9R (3.9 Å) | 38.3 → 16.2 Å | 0.42 | 1 (386 aa) | 383 |
| SalCIS | *Salmonella* subsp. *salamae* | 9R9H (2.28 Å) | 9R9I (3.4 Å) | 39.2 → 17.7 Å | 0.45 | 1 (519 aa) | 495 |
| CIS_Sc | *Streptomyces coelicolor* | 8BL4 (3.9 Å) | 8BKY (3.6 Å) | 39.8 → 17.2 Å | 0.43 | 1 (539 aa) | 367 |
| Diffocin | *Clostridioides difficile* | 8V3X (2.2 Å) | 8V3Y (3.6 Å) | 39.7 → 25.6 Å | 0.65 | 1 (CD1363, 354 aa) | 352 |

**Caveats to carry forward**

- **8BKY metadata is wrong.** Its `_em_helical_entity` reports 23.10°/38.50 Å,
  identical to extended 8BL4. The coordinates show genuine contraction
  (39.8 → 17.2 Å). Derive lattice geometry from coordinates, never from the
  deposited helical record. This alone justifies the ring-clustering approach.
- **Diffocin contracts least** (ratio 0.65 against 0.42–0.45 elsewhere). Either a
  milder contraction mechanism or a partially contracted reconstruction. Treat as
  a potential outlier and check before pooling.

## Tier B — usable, but requires symmetric expansion

| system | organism | extended | contracted | sheath proteins | usable res |
|---|---|---|---|---|---|
| AFP | *Serratia entomophila* | 6RBN (2.8 Å) | 6RC8 (3.8 Å) | **2** (Afp2 354 aa, Afp3 451 aa) | 348 / 355 |

Both depositions contain only the asymmetric unit (4 and 2 chains), so the
assembly must be generated in **both** states. That is acceptable because helical
parameters are deposited for both — the expansion is symmetric, unlike Vibrio.
AFP's deposited parameters are unusual (extended 40.27°/78.64 Å; contracted
3.14°/34.97 Å) and need interpreting before use. AFP is the only clean
**heterotypic two-protein sheath** available, so it is the natural test of whether
the metric behaves on non-homotypic sheaths.

## Exploratory only — do NOT use for the held-out test

| system | why |
|---|---|
| AlgoCIS 7AE0/7AEK | used to generate the hypothesis in the structural pilot |
| *V. cholerae* T6SS 5MXN/3J9G, 5MXN/5MYU | engineered VipA insertion; extended state is a single ring with no native inter-ring contacts |
| PlPVC1 9QGL+9QGN/9QGP | extended sheath resolved only as baseplate- and cap-proximal terminal rings, not a helical trunk; 3 paralogs, all wild-type, so retain as a secondary structural case |

## Excluded

| system | reason |
|---|---|
| Phage E217 8FUV/8FVG | extended entry is an asymmetric unit with **no** deposited helical parameters — cannot expand |
| phiTE 9CB9/9CBA | extended entry has no helical parameters |
| *P. aeruginosa* T6SS 5N8N | contracted only, no experimental extended |
| *Myxococcus* 5URW/5URX | 24–28 Å subtomogram averages, not atomic |
| T4 3FOH/3FOI, 3J2M | rigid-body fittings at 15–16 Å, not atomic models |
| *V. cholerae* 5OJQ | a pseudo-atomic model built into the 5MXN map, not an independent experiment |

## Count

**5 clean held-out systems + AFP** as a sixth with a two-protein sheath. Each
contributes 350–500 usable residues. That is sufficient to pre-register the
continuous remodeling metric and test it on data that played no part in
generating the hypothesis.

## Implementation consequence

`lattice.py` currently indexes against deposited twist/rise. Given 8BKY, that path
must become a fallback at most. The primary path should determine the symmetry
axis from coordinates (the singular vector outside the near-degenerate in-plane
pair — validated: recovers 40.8 Å for 7AE0 and 18.0 Å for 7AEK), cluster chains
into rings by axial coordinate, and index by (ring, azimuth). Deposited parameters
become a cross-check that is *reported*, not trusted.
