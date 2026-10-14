#!/usr/bin/env python3
"""Explicit zero-mode treatment for reproduction of the retained GNM values.

The archived producer uses diag(numpy.linalg.pinv(kirchhoff)). The supplied
initial reconstruction differs from the retained GNM values. The saved corrected
outputs agree within 1.7e-13 when this function excludes eigenvalues at or below
1e-10 times the maximum. This is evidence of numerical agreement, not proof of
the historical numerical backend or a guarantee of bitwise platform identity.
A graph Kirchhoff matrix has a constant null mode for each connected component.

The original producer remains unchanged. To use the explicit treatment:

    import covariates as cov
    from covariates_zeromode_fix import monomer_gnm_fluctuation
    cov.monomer_gnm_fluctuation = monomer_gnm_fluctuation
"""
from __future__ import annotations

import numpy as np

import covariates as cov

ZERO_MODE_REL_TOL = 1e-10   # eigenvalues below this fraction of the largest are modes to drop


def monomer_gnm_fluctuation(residues) -> dict:
    """GNM mean-square fluctuation per residue, null modes excluded explicitly.

    Same algorithm, same cutoff (covariates.GNM_CUTOFF) and same inputs as the
    archived producer; only the pseudo-inverse's treatment of the zero
    eigenvalue is made deterministic.
    """
    coords = cov._ca_coords(residues)
    ids = [r.get_id() for r in residues]
    keep = [i for i, c in enumerate(coords) if c is not None]
    if len(keep) < 10:
        return {rid: np.nan for rid in ids}

    X = np.array([coords[i] for i in keep])
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)
    A = ((D < cov.GNM_CUTOFF) & (D > 0)).astype(float)
    kirchhoff = np.diag(A.sum(axis=1)) - A

    w, V = np.linalg.eigh(kirchhoff)          # eigendecomposition of a symmetric matrix
    inv = np.zeros_like(w)
    nonzero = w > w.max() * ZERO_MODE_REL_TOL
    inv[nonzero] = 1.0 / w[nonzero]
    fluct = np.diag((V * inv) @ V.T)

    out = {rid: np.nan for rid in ids}
    for pos, i in enumerate(keep):
        out[ids[i]] = float(fluct[pos])
    return out
