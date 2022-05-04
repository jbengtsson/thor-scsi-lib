"""Utilities for phase advance calculations
"""
from .. import lib as tslib
import numpy as np
import copy
import logging

logger = logging.getLogger("thor_scsi")



def compute_nu(M):
    """

    Todo:
        How is it related to :func:`compute_dnu`
    """
    tr = M.trace()
    if tr < 2e0:
        nu = np.arccos(tr / 2e0) / (2e0 * np.pi)
        if M[0][1] < 0e0:
            nu = 1e0 - nu
        return nu, True
    else:
        return np.nan, False


def compute_nus(n_dof: int, M: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    """
    nus = np.zeros(n_dof)
    stable = [False, False]
    for k in range(n_dof):
        k2 = 2 * k
        submat = M[k2: k2 + 2, k2: k2 + 2]
        nu, flag = compute_nu(submat)
        nus[k] = nu
        stable[k] = flag
        if not flag:
            print(f"compute_nu: unstable for plane {k}")

    return nus, stable


def compute_nus_for_symplectic_matrix(
    n_dof: int, M: np.ndarray
) -> (np.ndarray, np.ndarray):
    """Compute nu for a general symplectic periodic transport matrix

    i.e., not assuming mid-plane symmetry.

    Todo:
       consider arrangement of variables (n_dof default to 2?)
       consider to return a masked array
       improve error reporting
    """
    n = 2 * n_dof

    nu, tr = np.zeros(2), np.zeros(2)

    # Should be a local copy.
    M1 = copy.copy(M[:n, :n])
    m_id = np.identity(n)
    detp = np.linalg.det(M1 - m_id)
    detm = np.linalg.det(M1 + m_id)
    for i in range(2):
        s = 2 * i
        e = 2 * i + 1
        tr[i] = M[s:e, s:e].trace()

    X_ = tslib.spatial_index.X
    Y_ = tslib.spatial_index.Y
    if tr[X_] > tr[Y_]:
        sgn = 1e0
    else:
        sgn = -1e0

    b = (detp - detm) / 16e0
    c = (detp + detm) / 8e0 - 1e0

    b2mc = b ** 2 - c
    if b2mc < 0e0:
        nu[:] = np.nan
        logger.warning("compute_nus_symplectic_matrix: unstable ...")
        return nu, False

    stable = True
    for i in range(2):
        if i == 0:
            x = -b + sgn * np.sqrt(b2mc)
        else:
            x = -b - sgn * np.sqrt(b2mc)

        if abs(x) <= 1e0:
            nu[i] = np.arccos(x) / (2e0 * np.pi)
            if M1[2 * i][2 * i + 1] < 0e0:
                nu[i] = 1e0 - nu[i]
            else:
                nu[i] = np.nan
                plane = ["hor", "vert"][i]
                txt = (
                    "compute_nus_symplectic_matrix:"
                    f" unstable {plane:%s} plane: x = {x:%10.3e}"
                )
                logger.warning(txt)
                stable = False
                # If the first plane is unstable why not to compute
                # it also for the second plane?
                return nu, False

    return nu, stable


__all__ = ["compute_nus_for_symplectic_matrix"]
