"""Coordinate-derived lattice reconstruction for contractile sheath assemblies.

Deposited helical metadata is NOT used to build the lattice. 8BKY (Streptomyces
CIS, contracted) reports twist/rise identical to its extended partner 8BL4
(23.10 deg / 38.50 A) while its coordinates show genuine contraction to 17.2 A
ring spacing. Anything built from that record would be wrong. Metadata is
therefore reconstructed independently and reported as a QC comparison only.

Where a deposition contains only the asymmetric unit (AFP 6RBN/6RC8), the
assembly is generated from the DEPOSITED ASSEMBLY OPERATORS -- authoritative
transforms, not a guessed axis -- and the generated assembly then goes through
the same coordinate-derived reconstruction as everything else. Expansion is
applied per state only where that state needs it, and the QC report records it,
because treating the two conformations asymmetrically is what invalidated the
Vibrio result.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

RING_GAP = 8.0          # A; axial gap separating one ring from the next
CLASH_DIST = 2.0        # A; heavy-atom separation below which an expansion is wrong
AXIS_MIN_CHAINS = 6     # fewer chains than this and the axis is not identifiable


class LatticeError(RuntimeError):
    """Raised when the lattice cannot be established. Never caught internally."""


@dataclass
class Subunit:
    key: str
    family: str                      # canonical cross-state family key
    ring: int = 0
    azim: int = 0
    resids: list = field(default_factory=list)      # author residue numbers
    canon: list = field(default_factory=list)       # canonical indices (None if unmatched)
    coords: list = field(default_factory=list)      # heavy-atom coords per residue
    ca: np.ndarray | None = None
    centroid: np.ndarray | None = None


@dataclass
class Lattice:
    axis: np.ndarray
    origin: np.ndarray
    subunits: list
    rise: float
    twist: float
    cn: int
    n_rings: int
    ring_sizes: list
    expanded: bool
    clashes: int
    repeat: int = 1


# ---------------------------------------------------------------------------
# assembly construction
# ---------------------------------------------------------------------------

def assembly_operators(cif: dict):
    """Deposited assembly transforms, as a list of (R, t)."""
    ids = cif.get("_pdbx_struct_oper_list.id", [])
    if not ids:
        return []
    ops = []
    for i in range(len(ids)):
        R = np.array([[float(cif[f"_pdbx_struct_oper_list.matrix[{r}][{c}]"][i])
                       for c in (1, 2, 3)] for r in (1, 2, 3)])
        t = np.array([float(cif[f"_pdbx_struct_oper_list.vector[{r}]"][i])
                      for r in (1, 2, 3)])
        ops.append((R, t))
    return ops


def build_subunits(model, cif, families: dict, want_expansion: bool):
    """Collect per-chain subunits, applying deposited operators when the
    deposition holds only an asymmetric unit."""
    from Bio.PDB.Polypeptide import is_aa

    base = []
    for ch in model:
        fam = families.get(ch.id)
        if fam is None:
            continue
        resids, coords, ca = [], [], []
        for r in ch:
            if not (is_aa(r, standard=True) and "CA" in r):
                continue
            heavy = np.array([a.get_coord() for a in r if a.element != "H"], dtype=float)
            if not len(heavy):
                continue
            resids.append(r.get_id()[1])
            coords.append(heavy)
            ca.append(r["CA"].get_coord())
        if len(resids) >= 20:
            base.append(Subunit(key=ch.id, family=fam, resids=resids,
                                coords=coords, ca=np.array(ca, dtype=float)))

    if not base:
        raise LatticeError("no chains matched the requested families")

    ops = assembly_operators(cif) if want_expansion else []
    if len(ops) <= 1:
        for s in base:
            s.centroid = s.ca.mean(axis=0)
        return base, False

    out = []
    for oi, (R, t) in enumerate(ops):
        for s in base:
            out.append(Subunit(
                key=f"{s.key}#{oi}", family=s.family, resids=list(s.resids),
                coords=[c @ R.T + t for c in s.coords],
                ca=s.ca @ R.T + t))
    for s in out:
        s.centroid = s.ca.mean(axis=0)
    return out, True


def count_clashes(subunits, sample=6):
    """Heavy-atom clashes between distinct subunits.

    A generated assembly built from wrong operators interpenetrates. This is the
    self-consistency check that an expansion-dependent system has to pass before
    its contacts mean anything.
    """
    from scipy.spatial import cKDTree
    if len(subunits) < 2:
        return 0
    pts, owner = [], []
    for i, s in enumerate(subunits[:sample * 4]):
        for c in s.coords:
            pts.append(c)
            owner.extend([i] * len(c))
    if not pts:
        return 0
    P = np.vstack(pts)
    owner = np.array(owner)
    tree = cKDTree(P)
    pairs = tree.query_pairs(r=CLASH_DIST, output_type="ndarray")
    if len(pairs) == 0:
        return 0
    return int(np.sum(owner[pairs[:, 0]] != owner[pairs[:, 1]]))


# ---------------------------------------------------------------------------
# lattice derivation, from coordinates alone
# ---------------------------------------------------------------------------

def symmetry_axis(centroids: np.ndarray):
    """The symmetry axis is the singular direction OUTSIDE the near-degenerate
    in-plane pair. A ring or tube has two equal in-plane extents; the remaining
    direction is the axis. Picking the smallest singular vector unconditionally
    fails for axially elongated particles (9QGL), and picking the largest fails
    for flat ring stacks (7AEK)."""
    ctr = centroids.mean(axis=0)
    _u, s, vt = np.linalg.svd(centroids - ctr)
    idx = min([(abs(s[0] - s[1]), 2), (abs(s[1] - s[2]), 0), (abs(s[0] - s[2]), 1)])[1]
    return vt[idx], ctr, s


def cylindrical(p, axis, origin, e1, e2):
    v = p - origin
    z = float(v @ axis)
    perp = v - z * axis
    return z, float(np.degrees(np.arctan2(perp @ e2, perp @ e1))), float(np.linalg.norm(perp))


def frame(axis):
    e1 = np.cross(axis, [1.0, 0.0, 0.0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.cross(axis, [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(axis, e1)


def derive(subunits, anchor_family: str, axis_hint: np.ndarray | None = None):
    """Cluster subunits into rings by axial coordinate and index them by
    (ring, azimuth). Returns the lattice with MEASURED rise, twist and Cn."""
    cents = np.array([s.centroid for s in subunits])
    if len(cents) < AXIS_MIN_CHAINS:
        raise LatticeError(
            f"only {len(cents)} subunits; the symmetry axis is not identifiable "
            f"from coordinates (minimum {AXIS_MIN_CHAINS})")

    axis, origin, _sv = symmetry_axis(cents)
    if axis_hint is not None and float(axis @ axis_hint) < 0:
        axis = -axis
    e1, e2 = frame(axis)

    # Rings are clustered over ALL families, not the anchor alone. Several of
    # these sheaths are heteromeric with the families ALTERNATING along the axis
    # (AFP: Afp2/Afp3; PlPVC1: three paralogs). Clustering on one family skips
    # every intermediate ring, reports a doubled rise, and silently removes each
    # subunit's immediate axial neighbours.
    anchor = subunits

    z = np.array([cylindrical(s.centroid, axis, origin, e1, e2)[0] for s in anchor])
    order = np.argsort(z)
    rings, cur = [], [order[0]]
    for prev, nxt in zip(order, order[1:]):
        if z[nxt] - z[prev] > RING_GAP:
            rings.append(cur)
            cur = [nxt]
        else:
            cur.append(nxt)
    rings.append(cur)

    ring_sizes = [len(r) for r in rings]
    cn = int(np.bincount(ring_sizes).argmax())
    ring_z = [float(np.mean(z[r])) for r in rings]
    rise = float(np.mean(np.diff(ring_z))) if len(ring_z) > 1 else float("nan")

    # assign every subunit (including non-anchor families) to the nearest ring
    for s in subunits:
        sz, sphi, _ = cylindrical(s.centroid, axis, origin, e1, e2)
        s.ring = int(np.argmin([abs(sz - rz) for rz in ring_z]))

    step = 360.0 / cn if cn else 360.0
    phi_of = {s.key: cylindrical(s.centroid, axis, origin, e1, e2)[1] for s in subunits}

    # Twist between consecutive rings, as a circular mean over every subunit pair
    # reduced modulo the Cn step. Using one ring's minimum azimuth as the estimator
    # is noisy -- it lands on an arbitrary subunit and was off by ~4 deg against
    # deposited values, enough to misplace a completed ring by several Angstrom.
    by_ring = {}
    for s in subunits:
        by_ring.setdefault(s.ring, []).append(s)

    def circ_mean_mod(vals, period):
        if not vals:
            return 0.0
        z = np.mean(np.exp(2j * np.pi * np.asarray(vals) / period))
        a = np.angle(z) * period / (2 * np.pi)
        return float(a)

    per_ring_twist = []
    rkeys = sorted(by_ring)
    for a, b in zip(rkeys, rkeys[1:]):
        diffs = [phi_of[y.key] - phi_of[x.key]
                 for x in by_ring[a] for y in by_ring[b]]
        per_ring_twist.append(circ_mean_mod(diffs, step))
    twist_per_ring = circ_mean_mod(per_ring_twist, step) if per_ring_twist else 0.0

    # Global azimuth index. Subunits sit at phi = phi0 + ring*twist + k*step, so
    # removing the accumulated twist leaves the strand index k. Indexing each ring
    # against its own minimum instead makes d_azimuth incomparable between the two
    # conformations, because they twist by different amounts per ring.
    # The twist branch must be the SAME in both conformations. A circular mean
    # lands in (-step/2, step/2], which puts extended at +19.9 deg and contracted
    # at -28.5 deg -- two branches of one C6 step differing by 60 deg. That shifts
    # the strand index by one per ring between the states and destroys every
    # cross-ring contact match. Pin the branch to [0, step).
    twist_assign = twist_per_ring % step
    r0 = rkeys[0]
    phi0 = circ_mean_mod([phi_of[s.key] for s in by_ring[r0]], step)
    for s in subunits:
        rel = phi_of[s.key] - (s.ring - r0) * twist_assign - phi0
        s.azim = int(round((rel % 360.0) / step)) % max(cn, 1)

    twists = per_ring_twist

    # Measured twist, between consecutive rings OF THE SAME FAMILY. In an
    # alternating heteromeric stack the next ring up holds a different protein,
    # so the like-to-like screw spans `repeat` rings; that screw is also what
    # ring completion must apply, since stepping by one ring would place a copy
    # of the wrong family at that position.
    anchor_rings = {}
    for s in subunits:
        if s.family != anchor_family:
            continue
        phi = cylindrical(s.centroid, axis, origin, e1, e2)[1]
        anchor_rings.setdefault(s.ring, []).append(phi)
    keys = sorted(anchor_rings)
    repeat = int(np.median(np.diff(keys))) if len(keys) > 1 else 1
    repeat = max(1, repeat)
    twist = twist_per_ring

    return Lattice(axis=axis, origin=origin, subunits=subunits, rise=rise,
                   twist=twist, cn=cn, n_rings=len(rings), ring_sizes=ring_sizes,
                   expanded=False, clashes=0, repeat=repeat)


def complete_rings(lat_obj, needed=1):
    """Extend the stack by `needed` same-family repeats at each end.

    Several depositions hold only two or three rings, so every subunit sits in a
    terminal ring and its contact set is truncated by where the deposited segment
    stops rather than by biology. Completion uses the COORDINATE-DERIVED screw --
    never the deposited record -- and is applied to both conformations alike:
    treating the states asymmetrically is what invalidated Vibrio.
    """
    rings = sorted({s.ring for s in lat_obj.subunits})
    lo, hi = rings[0], rings[-1]
    # Only extend when NO subunit can have a full axial shell as deposited.
    # A contact shell needs the rings at r-1 and r+1 -- any family, since
    # neighbouring rings in an alternating stack hold a different protein -- so
    # three rings suffice regardless of the same-family repeat. Requiring
    # 2*repeat+1 wrongly completed AFP (4 rings, repeat 2); the added rings then
    # moved the 'interior' boundary so the reference could be chosen next to
    # generated coordinates, manufacturing a dependence on them.
    if len(rings) >= 3:
        return 0
    axis, origin = lat_obj.axis, lat_obj.origin
    dz = lat_obj.rise * lat_obj.repeat
    dphi = np.radians(lat_obj.twist * lat_obj.repeat)

    added = []
    for k in list(range(-needed, 0)) + list(range(1, needed + 1)):
        ang = dphi * k
        w, c, s = axis, np.cos(ang), np.sin(ang)
        K = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
        R = np.eye(3) * c + s * K + (1 - c) * np.outer(w, w)
        shift = axis * (dz * k)
        for sub in lat_obj.subunits:
            if not (lo <= sub.ring + k * lat_obj.repeat <= hi):
                new = Subunit(
                    key=f"{sub.key}~{k}", family=sub.family,
                    ring=sub.ring + k * lat_obj.repeat, azim=sub.azim,
                    resids=list(sub.resids), canon=list(sub.canon),
                    coords=[(c0 - origin) @ R.T + origin + shift for c0 in sub.coords],
                    ca=(sub.ca - origin) @ R.T + origin + shift)
                new.centroid = new.ca.mean(axis=0)
                added.append(new)
    lat_obj.subunits = lat_obj.subunits + added
    lat_obj.expanded = bool(added)
    return len(added)


def orient_to(ext_lat: Lattice, ext_ref: Subunit, con_ref: Subunit):
    """Return an axis hint orienting the contracted lattice to match the extended one.

    The sign of a coordinate-derived axis is arbitrary. An inconsistent choice
    between states inverts every delta-ring, which silently destroys the
    cross-state contact correspondence the whole metric depends on. The monomer
    fold is conserved through contraction, so superposing the reference monomers
    transfers the extended frame onto the contracted one.
    """
    common = [(i, j) for i, ce in enumerate(ext_ref.canon)
              for j, cc in enumerate(con_ref.canon) if ce is not None and ce == cc]
    if len(common) < 10:
        raise LatticeError("too few cross-state matched residues to orient the axis")
    P = np.array([ext_ref.ca[i] for i, _ in common])
    Q = np.array([con_ref.ca[j] for _, j in common])
    cp, cq = P.mean(axis=0), Q.mean(axis=0)
    U, _s, Vt = np.linalg.svd((P - cp).T @ (Q - cq))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R @ ext_lat.axis


def metadata_qc(lat: Lattice, cif: dict):
    """Compare the coordinate-derived lattice with the deposited helical record."""
    def one(key):
        v = cif.get(key)
        try:
            return float(v[0])
        except (TypeError, ValueError, IndexError):
            return None

    dep_twist = one("_em_helical_entity.angular_rotation_per_subunit")
    dep_rise = one("_em_helical_entity.axial_rise_per_subunit")
    dep_sym = (cif.get("_em_helical_entity.axial_symmetry") or [None])[0]

    rows = []
    if dep_rise is None:
        rows.append(("rise", f"{lat.rise:.1f} A", "not deposited", "n/a"))
    else:
        # a deposited operator may span the whole `repeat`, not one ring
        ok = min(abs(abs(lat.rise) * k - dep_rise) for k in (1, lat.repeat)) \
            <= max(2.0, 0.1 * dep_rise)
        rows.append(("rise", f"{lat.rise:.1f} A", f"{dep_rise:.1f} A",
                     "PASS" if ok else "FAIL"))
    if dep_twist is None:
        rows.append(("twist", f"{lat.twist:.1f} deg", "not deposited", "n/a"))
    else:
        # Twist is only defined modulo the Cn step: a measured -26.9 deg in a C6
        # lattice is the same physical screw as a deposited 33.1 deg. Compare the
        # minimal circular difference, or every Cn system throws a false error.
        step = 360.0 / lat.cn if lat.cn else 360.0
        diff = min(abs((lat.twist * k - dep_twist + step / 2) % step - step / 2)
                   for k in (1, lat.repeat))
        ok = diff <= max(3.0, 0.15 * abs(dep_twist))
        rows.append(("twist", f"{lat.twist:.1f} deg", f"{dep_twist:.1f} deg",
                     "PASS" if ok else "FAIL"))
    if dep_sym is None:
        rows.append(("Cn", f"C{lat.cn}", "not deposited", "n/a"))
    else:
        rows.append(("Cn", f"C{lat.cn}", str(dep_sym),
                     "PASS" if str(dep_sym).upper() == f"C{lat.cn}" else "FAIL"))
    return rows
