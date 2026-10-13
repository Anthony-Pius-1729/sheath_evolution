"""Monomer-based structural covariates.

All three covariates are computed on the ISOLATED reference chain, never on the
assembly. That choice is load-bearing: assembly-based burial and packing are near
restatements of interface persistence, so conditioning on them would absorb the
very effect the study estimates. Monomer values describe burial within the
protein's own fold, leaving the intermolecular information entirely to the
persistence annotation.

The GNM is computed identically from every structure with no dependence on
refinement protocol, which makes flexibility comparable across depositions in a
way normalised B-factors are not.
"""

from __future__ import annotations

import numpy as np
from Bio.PDB.SASA import ShrakeRupley

# Tien et al. 2013, theoretical maximum accessible surface area (A^2).
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

GNM_CUTOFF = 7.3   # Angstrom, standard Gaussian network model contact radius
WCN_POWER = 2      # weighted contact number exponent


def _ca_coords(residues):
    out = []
    for r in residues:
        atom = r["CA"] if "CA" in r else None
        out.append(np.array(atom.get_coord(), dtype=float) if atom is not None else None)
    return out


def monomer_rsa(chain, residues) -> dict:
    """Relative solvent accessibility of each residue in the isolated chain."""
    ShrakeRupley().compute(chain, level="R")
    rsa = {}
    for r in residues:
        maxasa = MAX_ASA.get(r.get_resname())
        sasa = getattr(r, "sasa", None)
        rsa[r.get_id()] = (sasa / maxasa) if (maxasa and sasa is not None) else np.nan
    return rsa


def monomer_wcn(residues) -> dict:
    """Intramolecular weighted contact number: sum of 1/d^2 over other CA atoms
    in the same chain. Higher means more densely packed within the fold."""
    coords = _ca_coords(residues)
    ids = [r.get_id() for r in residues]
    wcn = {}
    for i, ci in enumerate(coords):
        if ci is None:
            wcn[ids[i]] = np.nan
            continue
        total = 0.0
        for j, cj in enumerate(coords):
            if i == j or cj is None:
                continue
            d2 = float(np.sum((ci - cj) ** 2))
            if d2 > 1e-6:
                total += 1.0 / (d2 ** (WCN_POWER / 2.0))
        wcn[ids[i]] = total
    return wcn


def monomer_gnm_fluctuation(residues) -> dict:
    """Gaussian network model mean-square fluctuation per residue.

    Kirchhoff matrix from CA atoms within GNM_CUTOFF; fluctuations are the
    diagonal of its pseudo-inverse. Identical algorithm and inputs for every
    structure, so values are comparable across depositions by construction.
    """
    coords = _ca_coords(residues)
    ids = [r.get_id() for r in residues]
    keep = [i for i, c in enumerate(coords) if c is not None]
    if len(keep) < 10:
        return {rid: np.nan for rid in ids}

    X = np.array([coords[i] for i in keep])
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)
    A = ((D < GNM_CUTOFF) & (D > 0)).astype(float)

    kirchhoff = np.diag(A.sum(axis=1)) - A
    fluct = np.diag(np.linalg.pinv(kirchhoff))

    out = {rid: np.nan for rid in ids}
    for pos, i in enumerate(keep):
        out[ids[i]] = float(fluct[pos])
    return out
