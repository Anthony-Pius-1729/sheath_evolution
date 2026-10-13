#!/usr/bin/env python3
"""Structural pilot: interface persistence in contractile injection system sheaths.

Answers one question and no more: do the matched conformational structures yield
sufficiently populated, reliably distinguishable interface-persistence classes to
support the planned evolutionary analysis?

No evolutionary rates are estimated here. The pilot either passes the
pre-registered minimum class size or it does not.

Pre-registered before any counts were computed:
  - minimum class size      : >= 20 mapped residues in the smallest compared class
  - contact definition      : heavy-atom distance < 5.0 A
  - partner identity        : (entity, screw index n, C6 index k), never chain ID
  - structural covariates   : monomer RSA, monomer WCN, monomer GNM fluctuation
  - covariate state         : computed in BOTH states; primary is the mean
  - primary contact scope   : sheath-sheath only; sheath-tube reported separately
"""

from __future__ import annotations

import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from scipy.spatial import cKDTree

import covariates as cov
import lattice as lat

HERE = Path(__file__).parent
STRUCT_DIR = HERE / "structures"
OUT_DIR = HERE / "output"

CONTACT_CUTOFF = 5.0        # Angstrom, heavy atom -- PRIMARY, pre-registered
SWEEP_CUTOFFS = [4.0, 4.5, 5.0, 5.5, 6.0]   # sensitivity analysis
MIN_CLASS_SIZE = 20         # pre-registered, not to be revised after seeing counts
ENGINEERING_SHELL = 8.0     # Angstrom, exclusion shell around internal indels
ENTITY_MATCH_IDENTITY = 0.90

INTERFACE_CLASSES = ["persistent_same_partner", "persistent_changed_partner",
                     "extended_only", "contracted_only"]
ALL_CLASSES = INTERFACE_CLASSES + ["non_interface"]

SYSTEMS = {
    "algoriphagus": {
        "label": "Algoriphagus machipongonensis eCIS (primary)",
        "extended": {"pdb": "7AE0", "twist": 20.54, "rise": 40.80},
        "contracted": {"pdb": "7AEK", "twist": 32.20, "rise": 18.04},
        "sheath_min_length": 400,   # sheath protein; the tube is far shorter
    },
    "vibrio": {
        "label": "Vibrio cholerae T6SS (secondary, engineered extended state)",
        "extended": {"pdb": "5MXN", "twist": 23.50, "rise": 37.80},
        "contracted": {"pdb": "3J9G", "twist": 29.40, "rise": 21.80},
        "sheath_min_length": 100,   # both VipA (159) and VipB (473) are sheath
    },
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def build_aligner():
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    return a


ALIGNER = build_aligner()


# --------------------------------------------------------------------------
# structure loading
# --------------------------------------------------------------------------

def fetch(pdb_id: str) -> Path:
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    path = STRUCT_DIR / f"{pdb_id}.cif"
    if not path.exists():
        print(f"    downloading {pdb_id} ...", flush=True)
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{pdb_id}.cif", path)
    return path


def chain_entities(path: Path) -> dict[str, dict]:
    """Map chain ID -> entity record, keyed by ENTITY ID rather than description.

    Descriptions are not unique: 5MXN labels both VipA and VipB "Type VI secretion
    protein". Keying on the description silently merges two distinct proteins.
    """
    d = MMCIF2Dict(str(path))
    eids = d.get("_entity.id", [])
    descs = d.get("_entity.pdbx_description", [])
    out = {}
    for i, eid in enumerate(d.get("_entity_poly.entity_id", [])):
        desc = descs[eids.index(eid)] if eid in eids else f"entity_{eid}"
        seq = d["_entity_poly.pdbx_seq_one_letter_code_can"][i].replace("\n", "")
        for ch in d["_entity_poly.pdbx_strand_id"][i].split(","):
            out[ch.strip()] = {"entity_id": eid, "desc": desc,
                               "seqres": seq, "length": len(seq)}
    return out


def aa_residues(chain):
    return [r for r in chain if is_aa(r, standard=True) and "CA" in r]


def seq_identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    aln = ALIGNER.align(a, b)[0]
    ident = sum(1
                for (e0, e1), (c0, c1) in zip(aln.aligned[0], aln.aligned[1])
                for off in range(e1 - e0)
                if a[e0 + off] == b[c0 + off])
    return ident / min(len(a), len(b))


# --------------------------------------------------------------------------
# per-state structural analysis
# --------------------------------------------------------------------------

def analyse_state(pdb_id, twist, rise, sheath_min_length, report):
    """Lattice-index one conformational state.

    Returns a dict keyed by entity_id holding the reference chain, its residues,
    contact partners labelled by lattice relationship, and monomer covariates.
    """
    path = fetch(pdb_id)
    model = MMCIFParser(QUIET=True).get_structure(pdb_id, str(path))[0]
    emap = chain_entities(path)

    chains = {}
    for ch in model:
        info = emap.get(ch.id)
        if not info:
            continue
        res = aa_residues(ch)
        if len(res) < 20:
            continue
        chains[ch.id] = {**info, "residues": res, "chain_obj": ch}

    by_entity = defaultdict(list)
    for cid, c in chains.items():
        by_entity[c["entity_id"]].append(cid)

    report.append(f"  {pdb_id}: {len(chains)} modelled polymer chains, "
                  f"{len(by_entity)} polymer entities")
    for eid, cids in sorted(by_entity.items()):
        c = chains[cids[0]]
        role = "sheath" if c["length"] >= sheath_min_length else "tube/other"
        report.append(f"      entity {eid} {c['desc'][:40]:42s} x{len(cids):<3d} "
                      f"construct={c['length']:<4d} modelled={len(c['residues']):<4d} [{role}]")

    sheath_ids = [eid for eid, cids in by_entity.items()
                  if chains[cids[0]]["length"] >= sheath_min_length]
    if not sheath_ids:
        raise lat.LatticeError(f"{pdb_id}: no entity meets the sheath length cutoff")

    # ---- helical axis from the most populous sheath entity ---------------
    anchor = max(sheath_ids, key=lambda e: len(by_entity[e]))
    anchor_coords = {
        cid: np.array([r["CA"].get_coord() for r in chains[cid]["residues"]])
        for cid in by_entity[anchor]
    }
    modal = Counter(len(v) for v in anchor_coords.values()).most_common(1)[0][0]
    anchor_coords = {k: v for k, v in anchor_coords.items() if len(v) == modal}
    if len(anchor_coords) < 3:
        raise lat.LatticeError(
            f"{pdb_id}: only {len(anchor_coords)} equal-length chains of the anchor "
            f"entity; cannot derive the helical axis")

    axis, origin, dev = lat.determine_axis(anchor_coords)
    centroids = {cid: np.mean([r["CA"].get_coord() for r in c["residues"]], axis=0)
                 for cid, c in chains.items()}

    def neighbour_count(cid):
        c0 = centroids[cid]
        return sum(1 for o, c in centroids.items()
                   if o != cid and np.linalg.norm(c - c0) < 90.0)

    anchor_ref = max(anchor_coords, key=neighbour_count)

    # ---- per-entity sublattices ------------------------------------------
    # Each entity occupies its own sublattice at a fixed radius and azimuthal
    # offset. Indexing every chain against the anchor's reference would make the
    # tube's offset look like a tolerance failure, so indices are computed per
    # entity relative to that entity's own reference subunit -- the copy nearest
    # the anchor reference, i.e. the one belonging to the same protomer.
    # The axis orientation is a property of the assembly, not of an entity, so the
    # sign is resolved once against the anchor and then imposed everywhere --
    # including on the symmetry expansion, which would otherwise generate rings on
    # the opposite side of the reference from the indices that label them.
    _a, _d, _b, sign = lat.assign_lattice_indices(
        {c: centroids[c] for c in by_entity[anchor]}, axis, origin,
        anchor_ref, twist, rise)
    axis = axis * sign

    assignments, all_diag, failures = {}, [], 0
    entity_ref = {}
    for eid, cids in by_entity.items():
        sub = {cid: centroids[cid] for cid in cids}
        ref = min(cids, key=lambda c: np.linalg.norm(centroids[c] - centroids[anchor_ref]))
        entity_ref[eid] = ref
        a, diag, bad = lat._assign_one(sub, axis, origin, ref, twist, rise)
        assignments.update(a)
        for d in diag:
            d["entity_id"] = eid
        all_diag.extend(diag)
        failures += len(bad)
        if bad:
            worst = max(bad, key=lambda d: d["resid_phi"])
            report.append(f"      WARNING entity {eid}: {len(bad)}/{len(cids)} chains "
                          f"outside lattice tolerance (worst chain {worst['chain']}: "
                          f"resid_z={worst['resid_z']:.1f} A, "
                          f"resid_phi={worst['resid_phi']:.1f} deg)")

    rings = sorted({n for n, _ in assignments.values()})
    report.append(f"      axis deviation (1-|cos|)={dev:.4f}; anchor reference "
                  f"'{anchor_ref}'; axis sign {sign:+.0f}; "
                  f"{len(all_diag) - failures}/{len(all_diag)} chains indexed in tolerance")
    report.append(f"      screw indices present: {rings}")
    if len(rings) < 3:
        report.append(f"      NOTE only {len(rings)} ring(s) deposited; neighbouring "
                      f"rings generated by the deposited helical operator and flagged")

    # ---- contacts for each sheath entity ---------------------------------
    results = {}
    for eid in sheath_ids:
        ref = entity_ref[eid]
        ref_res = chains[ref]["residues"]
        n0, k0 = assignments[ref]

        neighbours = []
        for cid, c in chains.items():
            if cid == ref:
                continue
            n, k = assignments[cid]
            neighbours.append(((c["entity_id"], n - n0, (k - k0) % 6),
                               c["residues"], False))
        if len(rings) < 3:
            for cid, c in chains.items():
                for shift in (-1, 1):
                    n, k = assignments[cid]
                    moved = _shift_residues(c["residues"], axis, origin,
                                            twist, rise, shift)
                    neighbours.append(((c["entity_id"], n - n0 + shift, (k - k0) % 6),
                                       moved, True))

        results[eid] = {
            "ref_chain": ref,
            "desc": chains[ref]["desc"],
            "seqres": chains[ref]["seqres"],
            "residues": ref_res,
            "partners": contact_partners(ref_res, neighbours, max(SWEEP_CUTOFFS)),
            "rsa": cov.monomer_rsa(chains[ref]["chain_obj"], ref_res),
            "wcn": cov.monomer_wcn(ref_res),
            "gnm": cov.monomer_gnm_fluctuation(ref_res),
            "sheath_ids": set(sheath_ids),
            "rings": rings,
            "lattice_failures": failures,
            "axis_sign": sign,
            "n_chains": len(chains),
        }

    catalog = {eid: {"desc": chains[cids[0]]["desc"],
                     "seqres": chains[cids[0]]["seqres"],
                     "is_sheath": eid in sheath_ids}
               for eid, cids in by_entity.items()}
    return results, catalog


def heavy_atoms(residues):
    coords, owner = [], []
    for r in residues:
        for a in r:
            if a.element != "H":
                coords.append(a.get_coord())
                owner.append(r.get_id())
    return np.array(coords, dtype=float), owner


def contact_partners(ref_residues, neighbours, max_cutoff):
    """Minimum heavy-atom distance from each reference residue to each lattice
    partner, out to `max_cutoff`.

    Storing distances rather than a thresholded set lets the contact definition be
    swept afterwards without recomputing anything, so the sensitivity analysis
    costs nothing and cannot drift from the primary run.
    """
    ref_coords, ref_owner = heavy_atoms(ref_residues)
    tree = cKDTree(ref_coords)
    partners = defaultdict(dict)
    for label, residues, generated in neighbours:
        coords, _ = heavy_atoms(residues)
        if len(coords) == 0:
            continue
        for j, idx_list in enumerate(tree.query_ball_point(coords, r=max_cutoff)):
            for i in idx_list:
                d = float(np.linalg.norm(ref_coords[i] - coords[j]))
                cur = partners[ref_owner[i]].get((label, generated))
                if cur is None or d < cur:
                    partners[ref_owner[i]][(label, generated)] = d
    return partners


class _Shifted:
    __slots__ = ("_id", "_coords")

    def __init__(self, rid, coords):
        self._id, self._coords = rid, coords

    def get_id(self):
        return self._id

    def __iter__(self):
        for c in self._coords:
            yield _Atom(c)


class _Atom:
    __slots__ = ("_c",)
    element = "C"

    def __init__(self, c):
        self._c = c

    def get_coord(self):
        return self._c


def _shift_residues(residues, axis, origin, twist, rise, n):
    out = []
    for r in residues:
        coords = np.array([a.get_coord() for a in r if a.element != "H"], dtype=float)
        if len(coords):
            out.append(_Shifted(r.get_id(),
                                lat.expand_screw(coords, axis, origin, twist, rise, n)))
    return out


# --------------------------------------------------------------------------
# cross-state correspondence and classification
# --------------------------------------------------------------------------

def correspond(res_ext, res_con):
    s_ext = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in res_ext)
    s_con = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in res_con)
    aln = ALIGNER.align(s_ext, s_con)[0]

    mapping, mismatches = {}, []
    for (e0, e1), (c0, c1) in zip(aln.aligned[0], aln.aligned[1]):
        for off in range(e1 - e0):
            re_, rc = res_ext[e0 + off], res_con[c0 + off]
            mapping[re_.get_id()] = rc.get_id()
            if re_.get_resname() != rc.get_resname():
                mismatches.append((re_.get_id()[1], re_.get_resname(), rc.get_resname()))

    matched = set(mapping)
    unmatched = [r for r in res_ext if r.get_id() not in matched]
    nums = [r.get_id()[1] for r in res_ext if r.get_id() in matched]
    lo, hi = (min(nums), max(nums)) if nums else (0, 0)
    terminal = [r for r in unmatched if not (lo <= r.get_id()[1] <= hi)]
    internal = [r for r in unmatched if lo <= r.get_id()[1] <= hi]
    return mapping, mismatches, terminal, internal


def classify(ext, con, mapping, cutoff):
    """Classify each cross-state-mapped residue at a given contact cutoff.

    Only contacts to CROSS-STATE COMPARABLE entities (those carrying a canonical
    "P" key) enter the persistence classes. An entity present in one state only --
    the inner tube, which every contracted deposition lacks because it is expelled
    on firing -- would otherwise be scored extended_only by construction, confusing
    tube expulsion with sheath interface rearrangement. Those contacts are recorded
    in separate columns instead.
    """
    def split(partners, rid):
        got = partners.get(rid, {})
        comparable = {lbl for (lbl, _g), d in got.items()
                      if d < cutoff and str(lbl[0]).startswith("P")}
        other = {lbl for (lbl, _g), d in got.items()
                 if d < cutoff and not str(lbl[0]).startswith("P")}
        gen = any(g for (lbl, g), d in got.items() if d < cutoff)
        return comparable, other, gen

    rows = []
    for rid_e, rid_c in mapping.items():
        pe, tube_e, gen_e = split(ext["partners"], rid_e)
        pc, tube_c, gen_c = split(con["partners"], rid_c)

        if not pe and not pc:
            cls = "non_interface"
        elif pe and not pc:
            cls = "extended_only"
        elif pc and not pe:
            cls = "contracted_only"
        elif pe & pc:
            cls = "persistent_same_partner"
        else:
            cls = "persistent_changed_partner"

        rows.append({
            "resid_ext": rid_e[1], "resid_con": rid_c[1],
            "persistence_class": cls,
            "n_sheath_partners_ext": len(pe), "n_sheath_partners_con": len(pc),
            "shared_partners": len(pe & pc),
            "contacts_single_state_entity_ext": len(tube_e) > 0,
            "contacts_single_state_entity_con": len(tube_c) > 0,
            "partners_ext": ";".join(sorted(f"{e}|{n}|{k}" for e, n, k in pe)),
            "partners_con": ";".join(sorted(f"{e}|{n}|{k}" for e, n, k in pc)),
            "uses_symmetry_generated": gen_e or gen_c,
        })
    return rows


def canonical_entities(cat_ext: dict, cat_con: dict, report):
    """Assign cross-state canonical keys to entities by sequence identity.

    Deposited entity IDs and descriptions are both unusable for this. The IDs are
    per-deposition (VipB is entity 2 in 5MXN but VipA is entity 1 in 3J9G, so a raw
    ID carried into a contact label makes one protomer's own partner look like a
    different molecule in the other state). Descriptions collide (5MXN calls both
    VipA and VipB "Type VI secretion protein"). Sequence identity is the only
    reliable key.

    Returns (canon_ext, canon_con, pairs) where the canon maps send entity IDs to
    a shared key.
    """
    canon_ext, canon_con, pairs = {}, {}, []
    used = set()
    for eid_e, E in sorted(cat_ext.items()):
        best, best_id = None, 0.0
        for eid_c, C in cat_con.items():
            if eid_c in used:
                continue
            ident = seq_identity(E["seqres"], C["seqres"])
            if ident > best_id:
                best, best_id = eid_c, ident
        if best is not None and best_id >= ENTITY_MATCH_IDENTITY:
            key = f"P{len(pairs) + 1}"
            canon_ext[eid_e], canon_con[best] = key, key
            used.add(best)
            pairs.append((eid_e, best, best_id, key))
            report.append(f"      {key}: extended entity {eid_e} "
                          f"({E['desc'][:30]}) <-> contracted entity {best} "
                          f"({cat_con[best]['desc'][:30]}) at {best_id:.3f} identity")
        else:
            canon_ext[eid_e] = f"EXTONLY_{eid_e}"
            report.append(f"      extended entity {eid_e} ({E['desc'][:30]}) has no "
                          f"contracted counterpart above {ENTITY_MATCH_IDENTITY} "
                          f"identity (best {best_id:.3f}) -- not cross-state comparable")
    for eid_c in cat_con:
        canon_con.setdefault(eid_c, f"CONONLY_{eid_c}")
    return canon_ext, canon_con, pairs


def remap_partners(partners, canon):
    """Translate raw entity IDs in contact labels to canonical cross-state keys."""
    out = {}
    for rid, labels in partners.items():
        remapped = {}
        for ((e, n, k), gen), dist in labels.items():
            key = ((canon.get(e, f"UNK_{e}"), n, k), gen)
            if key not in remapped or dist < remapped[key]:
                remapped[key] = dist
        out[rid] = remapped
    return out


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def run_system(sys_key, spec):
    report = ["", "=" * 78, f"SYSTEM: {spec['label']}", "=" * 78]
    e, c = spec["extended"], spec["contracted"]
    report.append(f"  extended  : {e['pdb']}  twist={e['twist']} deg  rise={e['rise']} A")
    report.append(f"  contracted: {c['pdb']}  twist={c['twist']} deg  rise={c['rise']} A")

    ext, cat_e = analyse_state(e["pdb"], e["twist"], e["rise"],
                               spec["sheath_min_length"], report)
    con, cat_c = analyse_state(c["pdb"], c["twist"], c["rise"],
                               spec["sheath_min_length"], report)

    report.append("\n  cross-state entity matching (by sequence identity):")
    canon_e, canon_c, pairs = canonical_entities(cat_e, cat_c, report)

    for eid, R in ext.items():
        R["partners"] = remap_partners(R["partners"], canon_e)
        R["sheath_ids"] = {canon_e[i] for i in R["sheath_ids"]}
    for eid, R in con.items():
        R["partners"] = remap_partners(R["partners"], canon_c)
        R["sheath_ids"] = {canon_c[i] for i in R["sheath_ids"]}

    tables, summary = {}, []
    for eid_e, eid_c, ident, canon_key in pairs:
        if eid_e not in ext or eid_c not in con:
            continue
        E, C = ext[eid_e], con[eid_c]
        mapping, mismatches, terminal, internal = correspond(E["residues"], C["residues"])
        # Descriptions are not unique within a deposition, so the canonical key
        # disambiguates both the report and the output filenames.
        name = f"{canon_key}_{E['desc']}_{len(E['seqres'])}aa"

        report.append(f"\n  --- {name} (extended entity {eid_e}) ---")
        report.append(f"      reference chains      : {E['ref_chain']} (ext) / "
                      f"{C['ref_chain']} (con)")
        report.append(f"      modelled residues     : {len(E['residues'])} ext / "
                      f"{len(C['residues'])} con")
        report.append(f"      cross-state mapped    : {len(mapping)} "
                      f"({100 * len(mapping) / max(len(E['residues']), 1):.1f}% of extended)")
        report.append(f"      substitutions         : {len(mismatches)}")
        report.append(f"      unmatched terminal    : {len(terminal)} (coverage loss)")
        report.append(f"      unmatched internal    : {len(internal)}"
                      + (f" -> extended resids {[r.get_id()[1] for r in internal]}"
                         if internal else ""))

        df = pd.DataFrame(classify(E, C, mapping, CONTACT_CUTOFF))

        if internal:
            eng = np.array([r["CA"].get_coord() for r in internal])
            ca = {r.get_id()[1]: r["CA"].get_coord() for r in E["residues"]}
            df["dist_to_engineering"] = [
                float(np.min(np.linalg.norm(eng - ca[r], axis=1))) if r in ca else np.nan
                for r in df["resid_ext"]]
            df["engineering_excluded"] = df["dist_to_engineering"] < ENGINEERING_SHELL
            report.append(f"      engineering shell     : "
                          f"{int(df['engineering_excluded'].sum())} residues within "
                          f"{ENGINEERING_SHELL} A of an internal indel -> excluded")
        else:
            df["dist_to_engineering"] = np.nan
            df["engineering_excluded"] = False

        ids_e = {r.get_id()[1]: r.get_id() for r in E["residues"]}
        ids_c = {r.get_id()[1]: r.get_id() for r in C["residues"]}
        for nm, se, sc in (("rsa", E["rsa"], C["rsa"]), ("wcn", E["wcn"], C["wcn"]),
                           ("gnm", E["gnm"], C["gnm"])):
            ve = [se.get(ids_e.get(r), np.nan) for r in df["resid_ext"]]
            vc = [sc.get(ids_c.get(r), np.nan) for r in df["resid_con"]]
            df[f"{nm}_ext"], df[f"{nm}_con"] = ve, vc
            df[f"{nm}_mean"] = np.nanmean(np.vstack([ve, vc]), axis=0)

        tables[name] = df
        kept = df[~df["engineering_excluded"]]
        counts = df["persistence_class"].value_counts().to_dict()
        counts_kept = kept["persistence_class"].value_counts().to_dict()

        report.append("\n      persistence class (all mapped / after engineering shell):")
        for cls in ALL_CLASSES:
            report.append(f"        {cls:30s} {counts.get(cls, 0):5d}  /  "
                          f"{counts_kept.get(cls, 0):5d}")
        report.append(f"        {'[contacts single-state entity]':30s} "
                      f"{int(df['contacts_single_state_entity_ext'].sum()):5d}"
                      f"   (excluded from the classes above)")

        smallest = min(counts_kept.get(x, 0) for x in INTERFACE_CLASSES)
        verdict = "PASS" if smallest >= MIN_CLASS_SIZE else "FAIL"
        report.append(f"      smallest interface class = {smallest} "
                      f"(threshold {MIN_CLASS_SIZE}) -> {verdict}")

        # Sensitivity of the pre-registered verdict to the contact definition.
        report.append("\n      contact-cutoff sweep (smallest interface class):")
        sweep_line = []
        for cut in SWEEP_CUTOFFS:
            sdf_ = pd.DataFrame(classify(E, C, mapping, cut))
            sc = sdf_["persistence_class"].value_counts().to_dict()
            sweep_line.append(f"{cut:.1f}A:{min(sc.get(x, 0) for x in INTERFACE_CLASSES)}")
        report.append("        " + "   ".join(sweep_line)
                      + f"   (primary {CONTACT_CUTOFF:.1f}A)")

        # POST HOC -- reported for design purposes only. The pre-registered
        # four-way interface split stands as the verdict above; this collapse was
        # chosen after seeing the counts and must not be substituted for it
        # without re-registering on data not used here.
        persistent = counts_kept.get("persistent_same_partner", 0) + \
            counts_kept.get("persistent_changed_partner", 0)
        state_specific = counts_kept.get("extended_only", 0) + \
            counts_kept.get("contracted_only", 0)
        report.append(f"\n      [POST HOC] collapsed three-way scheme: "
                      f"persistent={persistent}, state_specific={state_specific}, "
                      f"non_interface={counts_kept.get('non_interface', 0)} "
                      f"-> smallest={min(persistent, state_specific)}")

        report.append("\n      monomer covariates by class (state-averaged):")
        report.append(f"        {'class':30s} {'n':>5s} {'RSA':>7s} {'WCN':>8s} {'GNM':>9s}")
        for cls in ALL_CLASSES:
            g = kept[kept["persistence_class"] == cls]
            if len(g) == 0:
                continue
            report.append(f"        {cls:30s} {len(g):5d} {g['rsa_mean'].mean():7.3f} "
                          f"{g['wcn_mean'].mean():8.3f} {g['gnm_mean'].mean():9.4f}")

        summary.append({
            "system": sys_key, "entity": name, "identity": round(ident, 3),
            "modelled_ext": len(E["residues"]), "modelled_con": len(C["residues"]),
            "mapped": len(mapping), "substitutions": len(mismatches),
            "unmatched_internal": len(internal),
            "engineering_excluded": int(df["engineering_excluded"].sum()),
            "rings_ext": len(E["rings"]), "rings_con": len(C["rings"]),
            "symmetry_generated_used": bool(df["uses_symmetry_generated"].any()),
            "lattice_failures_ext": E["lattice_failures"],
            "lattice_failures_con": C["lattice_failures"],
            **{f"n_{x}": counts_kept.get(x, 0) for x in ALL_CLASSES},
            "smallest_interface_class": smallest, "verdict": verdict,
        })
    return report, tables, summary


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = [
        "CONTRACTILE SHEATH INTERFACE-PERSISTENCE STRUCTURAL PILOT",
        "",
        "Pre-registered parameters (fixed before any counts were computed):",
        f"  contact cutoff          {CONTACT_CUTOFF} A heavy atom",
        f"  minimum class size      {MIN_CLASS_SIZE} residues",
        f"  engineering shell       {ENGINEERING_SHELL} A",
        f"  GNM cutoff              {cov.GNM_CUTOFF} A",
        f"  WCN exponent            {cov.WCN_POWER}",
        "  partner identity        (entity, screw index n, C6 index k)",
        "  covariates              monomer only; primary = mean of both states",
        "  primary contact scope   sheath-sheath; sheath-tube reported separately",
        "",
        "No evolutionary rates are estimated in this pilot.",
    ]
    all_summary = []
    for sys_key, spec in SYSTEMS.items():
        try:
            rep, tables, summary = run_system(sys_key, spec)
        except lat.LatticeError as exc:
            out += ["", "=" * 78, f"SYSTEM: {spec['label']}", "=" * 78,
                    f"  HALTED: {exc}", "  No classification emitted."]
            continue
        out += rep
        all_summary += summary
        for name, df in tables.items():
            safe = "".join(ch if ch.isalnum() else "_" for ch in name)[:40]
            df.to_csv(OUT_DIR / f"residues_{sys_key}_{safe}.csv", index=False)

    if all_summary:
        sdf = pd.DataFrame(all_summary)
        sdf.to_csv(OUT_DIR / "summary.csv", index=False)
        out += ["", "=" * 78, "OVERALL", "=" * 78,
                sdf[["system", "entity", "mapped"]
                    + [f"n_{x}" for x in INTERFACE_CLASSES]
                    + ["smallest_interface_class", "verdict"]].to_string(index=False),
                "", f"  protein families passing: "
                    f"{int((sdf['verdict'] == 'PASS').sum())} / {len(sdf)}"]

    text = "\n".join(out)
    (OUT_DIR / "qc_report.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    sys.exit(main())
