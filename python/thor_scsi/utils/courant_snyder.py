from .. import lib as tslib
from .output import mat2txt

import numpy as np
import logging


logger = logging.getLogger("thor_scsi")

sign = np.sign


# Configuration space coordinates.
X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]
# Phase-space coordinates.
[x_, px_, y_, py_, ct_, delta_] = [
    tslib.phase_space_index_internal.x,
    tslib.phase_space_index_internal.px,
    tslib.phase_space_index_internal.y,
    tslib.phase_space_index_internal.py,
    tslib.phase_space_index_internal.ct,
    tslib.phase_space_index_internal.delta,
]


def compute_A(
        eta: np.ndarray, alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    r"""Compute twiss diagonalisation

    thus A - i.e., which diagonalises the Poincaré map:

    Args:
        eta:     dispersion
        alpha:  :math:`- 1/2 \beta'`
        beta:    beta function
    
    .. math::
         M = A \cdot R \cdot A^{-1}
    
    from the Twiss parameters for each transversal plane x, y:
    #    twiss = [eta[], alpha[], beta[]].
    
    """
    eta   = np.asarray(eta)
    alpha = np.asarray(alpha)
    beta  = np.asarray(beta)
    
    # eta, alpha, beta = twiss
    A = np.identity(6)
    for k in range(2):
        A[2*k, 2*k] = np.sqrt(beta[k])
        A[2*k+1, 2*k] = -alpha[k]/np.sqrt(beta[k])
        A[2*k+1, 2*k+1] = 1e0/np.sqrt(beta[k])
    A[ct_, ct_] = 1e0
    A[delta_, delta_] =1e0

    B = np.identity(6)
    B[x_, delta_] = eta[x_]
    B[px_, delta_] = eta[px_]
    B[ct_, x_] = eta[px_]
    B[ct_, px_] = -eta[x_]
 
    A = B @ A

    # logger.info("\ncompute_A\nA:\n" + mat2txt(A))
    return A


def compute_dnu(n_dof, A):
    """

    a lot of special treatment for the internal representation of the matrix

    Todo:
        Consider reordering the arguments
        comment why case delta / ct needs to be treated specially?
        And why does it need a minus in front

        Compare function name to :func:`phase_advance.compute_nu`
    """
    eps = 1e-15

    ct_    = tslib.phase_space_index_internal.ct
    delta_ = tslib.phase_space_index_internal.delta

    a = np.zeros(n_dof)

    for k in range(n_dof):
        if k < 2:
            k2 = 2 * k
            dnu = np.arctan2(A[k2][k2 + 1], A[k2][k2]) / (2 * np.pi)
        else:
            dnu = -np.arctan2(A[ct_][delta_], A[ct_][ct_]) / (2 * np.pi)
        if dnu < -eps:
            dnu += 1
        a[k] = dnu
    return a


def compute_A_CS(n_dof, A):
    """compute Courant Snyder form of A

    Rotate A to zero a12 (or equivalent for the other planes)

    Todo:
        check if dnu and R should not match in shape
    """

    n = 2 * n_dof
    dnu = np.zeros(n_dof)
    R = np.identity(6)

    nr, nc = R.shape
    assert nr == 6
    assert nc == 6

    nr, nc = A.shape
    assert nr == 6
    assert nc == 6

    dnu = compute_dnu(n_dof, A)

    # Build required rotation submatrices
    for k in range(n_dof):
        arg = 2 * np.pi * dnu[k]
        c = np.cos(arg)
        s = np.sin(arg)
        # Rotation matrix
        k2 = 2 * k
        if True:
            # fmt: off
            tmp = np.array([
                [c, -s],
                [s,  c]
            ])
            R[k2 : k2 + 2, k2 : k2 + 2] = tmp
            # fmt: on
            # continue
        else:
            R[k2, k2], R[k2][k2 + 1] = c, -s
            R[2 * k + 1][2 * k], R[2 * k + 1][2 * k + 1] = s, c

    # print(f"Rotation matrix\n{mat2txt(R)}")
    # Rotate A matrix back
    Ar = A @ R
    return Ar, dnu


def compute_Twiss_A(A):
    """
    """
    n_dof = 2
    n = 2 * n_dof

    [eta, alpha, beta] = [np.zeros(n), np.zeros(n_dof), np.zeros(n_dof)]

    delta_ = tslib.phase_space_index_internal.delta

    for k in range(n_dof):
        k2 = 2 * k
        eta[k2] = A[k2][delta_]
        eta[k2 + 1] = A[k2 + 1][delta_]
        alpha[k] = \
            -(A[k2][k2] * A[k2 + 1][k2] + A[k2][k2 + 1] * A[k2 + 1][k2 + 1])
        beta[k] = A[k2][k2] ** 2 + A[k2][k2 + 1] ** 2
    dnu = compute_dnu(n_dof, A)

    return np.array(eta), np.array(alpha), np.array(beta), np.array(dnu)


def compute_Twiss_A_A_tp(A):
    """

    Todo:
      difference ti compute_Twiss_A?
    """
    n_dof = 2
    n = 2 * n_dof

    eta, alpha, beta = np.zeros(n), np.zeros(n_dof), np.zeros(n_dof)

    delta_ = tslib.phase_space_index_internal.delta
    A_A_tp = np.dot(A[0:n, 0:n], np.transpose(A[0:n, 0:n]))

    return compute_Twiss_A(A)

    for k in range(n_dof):
        k2 = 2 * k
        eta[k2], eta[k2 + 1] = A[k2][delta_], A[k2 + 1][delta_]
        alpha[k] = -A_A_tp[k2][k2 + 1]
        beta[k] = A_A_tp[k2][k2]
    dnu = compute_dnu(n_dof, A)

    return eta, alpha, beta, dnu


def compute_dispersion(M):
    """

    Todo:
        Define which standard M should follow standard rules or
        thor_scsi / tracy internal rules
    """
    n = 4
    I = np.identity(n)
    D = M[:n, tslib.phase_space_index_internal.delta]

    return np.dot(np.linalg.inv(I - M[:n, :n]), D)


def compute_Twiss_M(M):
    n_dof = 2

    alpha, beta, nu = np.zeros(n_dof), np.zeros(n_dof), np.zeros(n_dof)

    eta = compute_dispersion(M)

    stable = [True, True]
    for k in range(n_dof):
        k2 = 2 * k
        cos = M[k2 : k2 + 2, k2 : k2 + 2].trace() / 2e0
        if abs(cos) >= 1e0:
            print("\ncompute_Twiss_M: {:5.3f}\n".format(cos))
            stable[k] = False
        sin = np.sqrt(1e0 - cos ** 2) * sign(M[k2][k2 + 1])
        alpha[k] = (M[k2][k2] - M[k2 + 1][k2 + 1]) / (2e0 * sin)
        beta[k] = M[k2][k2 + 1] / sin
        nu[k] = np.arctan2(sin, cos) / (2e0 * np.pi)
        if nu[k] < 0e0:
            nu[k] += 1e0

    return [eta, alpha, beta, nu, stable]

__all__ = ["compute_A", "compute_dnu", "compute_A_CS", "compute_Twiss_A",
           "compute_Twiss_A_A_tp", "compute_dispersion", "compute_Twiss_M"]
