"""Frozen contact-network remodeling metric.

For residue i in state s the contact set is

    C_i^s = { (partner_family, d_ring, d_azimuth, partner_residue) }

and the remodeling index is

    R_i = 1 - |C_i^E & C_i^C| / |C_i^E | C_i^C|

R_i = 0 means the intermolecular contact network is unchanged by contraction;
R_i = 1 means it is completely remodeled.

PRE-REGISTERED, fixed before any held-out system was run:
  primary domain          |C_i^E | C_i^C| >= 3
  sensitivity domains     2, 4, 5
  primary contact cutoff  5.0 A heavy atom
  sensitivity cutoffs     4.0, 4.5, 5.5, 6.0 A
  primary contact identity    lattice-aware, the 4-tuple above
  robustness identity         lattice stripped: (partner_family, partner_residue)

IMPLEMENTATION DECISION not covered by the frozen specification, applied
symmetrically to both states and reported in the QC output: a contact whose
PARTNER residue has no cross-state counterpart is excluded from both C_i^E and
C_i^C. Partner residues are addressed by canonical index -- the column of the
extended-to-contracted alignment -- rather than by author numbering, so the two
states are directly comparable. Including unmatched partners would let coverage
differences between depositions masquerade as remodeling.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

PRIMARY_CUTOFF = 5.0
SWEEP_CUTOFFS = [4.0, 4.5, 5.0, 5.5, 6.0]
PRIMARY_MIN_UNION = 3
SWEEP_MIN_UNION = [2, 3, 4, 5]


def contact_sets(ref, subunits, max_cutoff=max(SWEEP_CUTOFFS)):
    """Minimum heavy-atom distance from each reference residue to each labelled
    partner. Distances are retained so the cutoff can be swept without
    recomputation, which keeps the sensitivity analysis exactly consistent with
    the primary run.

    Returns {ref_canonical_index: {(family, d_ring, d_azim, partner_canon): dist}}
    and the number of contacts dropped because the partner was unmatched.
    """
    ref_pts, ref_owner = [], []
    for canon, coords in zip(ref.canon, ref.coords):
        if canon is None:
            continue
        for p in coords:
            ref_pts.append(p)
            ref_owner.append(canon)
    if not ref_pts:
        return {}, 0
    tree = cKDTree(np.array(ref_pts))
    ref_owner = np.array(ref_owner)

    ref_pts = np.array(ref_pts)
    cn = max(1, _cn_of(subunits))
    # Each (family, d_ring, d_azimuth) names exactly one subunit, so whether a
    # label rests on ring-completed coordinates is well defined per label.
    out, dropped, gen_of_label = {}, 0, {}
    for s in subunits:
        if s is ref:
            continue
        d_ring = s.ring - ref.ring
        d_azim = (s.azim - ref.azim) % cn

        is_gen = "~" in s.key
        pts, labels = [], []
        for canon, coords in zip(s.canon, s.coords):
            if canon is None:
                dropped += 1
                continue
            lbl = (s.family, d_ring, d_azim, canon)
            gen_of_label[lbl] = is_gen
            pts.append(coords)
            labels.extend([lbl] * len(coords))
        if not pts:
            continue
        pts = np.vstack(pts)

        for j, idx_list in enumerate(tree.query_ball_point(pts, r=max_cutoff)):
            if not idx_list:
                continue
            lbl = labels[j]
            dists = np.linalg.norm(ref_pts[idx_list] - pts[j], axis=1)
            for i, d in zip(idx_list, dists):
                slot = out.setdefault(int(ref_owner[i]), {})
                if lbl not in slot or d < slot[lbl]:
                    slot[lbl] = float(d)
    return out, dropped, gen_of_label


def _cn_of(subunits):
    sizes = {}
    for s in subunits:
        sizes.setdefault((s.family, s.ring), 0)
        sizes[(s.family, s.ring)] += 1
    if not sizes:
        return 1
    vals = list(sizes.values())
    return int(np.bincount(vals).argmax())


def _threshold(sets, cutoff, strip_lattice):
    out = {}
    for rid, labels in sets.items():
        keep = set()
        for lbl, d in labels.items():
            if d >= cutoff:
                continue
            keep.add((lbl[0], lbl[3]) if strip_lattice else lbl)
        out[rid] = keep
    return out


def remodeling_index(ext_sets, con_sets, cutoff=PRIMARY_CUTOFF,
                     min_union=PRIMARY_MIN_UNION, strip_lattice=False):
    """R_i for every residue present in both states that meets the domain."""
    E = _threshold(ext_sets, cutoff, strip_lattice)
    C = _threshold(con_sets, cutoff, strip_lattice)

    rows = []
    for rid in sorted(set(E) | set(C)):
        ce, cc = E.get(rid, set()), C.get(rid, set())
        union = ce | cc
        if len(union) < min_union:
            rows.append({"canon": rid, "R": np.nan, "n_union": len(union),
                         "n_inter": len(ce & cc), "n_ext": len(ce),
                         "n_con": len(cc), "in_domain": False})
            continue
        rows.append({"canon": rid, "R": 1.0 - len(ce & cc) / len(union),
                     "n_union": len(union), "n_inter": len(ce & cc),
                     "n_ext": len(ce), "n_con": len(cc), "in_domain": True})
    return rows


def distribution_summary(values):
    """Is R_i continuously distributed, or has it collapsed onto 0 and 1?"""
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if len(v) == 0:
        return {}
    interior = float(np.mean((v > 0.05) & (v < 0.95)))
    return {
        "n": int(len(v)), "mean": float(v.mean()), "median": float(np.median(v)),
        "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
        "q10": float(np.quantile(v, 0.10)), "q90": float(np.quantile(v, 0.90)),
        "frac_zero": float(np.mean(v <= 0.05)),
        "frac_one": float(np.mean(v >= 0.95)),
        "frac_interior": interior,
    }


def generated_dependence(ext_sets, con_sets, gen_e, gen_c,
                         cutoff=PRIMARY_CUTOFF, min_union=PRIMARY_MIN_UNION):
    """How much of each residue's PERSISTENT contact set rests on completed rings.

    f_generated_i = (persistent contacts involving generated coordinates)
                    / |C_E & C_C|

    Also recomputes R with every generated-coordinate contact removed from BOTH
    states, so the comparison is symmetric. If R is materially unchanged, the
    completion is provenance rather than the source of the result.
    """
    E = _threshold(ext_sets, cutoff, False)
    C = _threshold(con_sets, cutoff, False)
    rows = []
    for rid in sorted(set(E) | set(C)):
        ce, cc = E.get(rid, set()), C.get(rid, set())
        inter, union = ce & cc, ce | cc
        gen_inter = {l for l in inter if gen_e.get(l) or gen_c.get(l)}

        ce_n = {l for l in ce if not gen_e.get(l)}
        cc_n = {l for l in cc if not gen_c.get(l)}
        un_n, in_n = ce_n | cc_n, ce_n & cc_n

        rows.append({
            "canon": rid,
            "R": 1.0 - len(inter) / len(union) if len(union) >= min_union else np.nan,
            "R_native_only": 1.0 - len(in_n) / len(un_n) if len(un_n) >= min_union else np.nan,
            "n_inter": len(inter), "n_inter_generated": len(gen_inter),
            "f_generated": len(gen_inter) / len(inter) if inter else np.nan,
            "n_union": len(union), "n_union_native": len(un_n),
        })
    return rows
