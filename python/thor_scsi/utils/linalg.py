"""Eigen for symplectic matrices ?

"""
import numpy as np
from typing import Sequence
import copy as _copy
import thor_scsi.lib as tslib
import logging
import functools
from .phase_space_vector import omega_block_matrix as _obm

logger = logging.getLogger("thor_scsi")


@functools.lru_cache(maxsize=3)
def omega_block_matrix(n_dof):
    """Compute the omega matrix

    Required to test the symplectisism for a matrix in block style

    Todo:
       check if it is correctly implemented for the common 2.5
       phase space ordering
    """
    return _obm(n_dof)

    # n = 2 * n_dof
    # S = np.zeros((n, n))
    # for k in range(n_dof):
    #     s = 2 * k
    #     e = 2 * k + 1
    #     S[s, e] = 1e0
    #     S[e, s] = -1e0
    # return S


def swap(w, i, j):
    w[i], w[j] = w[j], w[i]


def swap_imag_part(a: complex, b: complex) -> (complex, complex):
    return a.real + b.imag * 1j, b.real + a.imag * 1j


def swap_imag(w: Sequence, i: int, j: int):
    w[i], w[j] = swap_imag_part(w[i], w[j])


def swap_mat(A, i, j):
    A[:, [i, j]] = A[:, [j, i]]


def swap_mat_imag(A, i, j):
    n = len(A)
    for k in range(n):
        c = A[k][i].imag
        A[k][i] = A[k][i].real + A[k][j].imag * 1j
        A[k][j] = A[k][j].real + c * 1j


def closest(x: float, x1: float, x2: float, x3: float) -> int:
    """Find the value that is closest to x

    Todo:
        check if corner cases need to be checked
        e.g. dx1 == dx2 ...
    """

    dx1, dx2, dx3 = [np.abs(x - tx) for tx in (x1, x2, x3)]
    if dx1 < dx2 and dx1 < dx3:
        k = 1
    elif dx2 < dx1 and dx2 < dx3:
        k = 2
    else:
        k = 3
    return k


def sort_eigen(n_dof, M, w, v):
    """Sort eigen values to be in order x, y, longitudinal

    Todo:
        can n_dof be derived from the shape of M, w or v?

        Check that shape of M, w, v correspond
    """
    sign = np.sign
    n = 2 * n_dof

    sin_M, cos_M = np.zeros(n_dof), np.zeros(n_dof)
    nu1_M, nu2_M, nu1, nu2 = np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

    for i in range(n_dof):
        j = (i + 1) * 2 - 1
        v = (M[j - 1][j - 1] + M[j][j]) / 2
        cos_M[i] = v
        if np.abs(v) > 1e0:
            logger.warning(
                "sort_eigen: unstable |cos_M[nu_{i +1}]-1e0| = {:10.3e}",
                i + 1,
                np.absolute(v - 1e0),
            )
            stable = False
        sin_M[i] = sign(M[j - 1][j]) * np.sqrt(1e0 - cos_M[i] ** 2)
        nu1_M[i] = np.arctan2(sin_M[i], cos_M[i]) / (2 * np.pi)
        if nu1_M[i] < 0.0:
            nu1_M[i] += 1.0
        if nu1_M[i] <= 0.5:
            nu2_M[i] = nu1_M[i]
        else:
            nu2_M[i] = 1.0 - nu1_M[i]
        nu1[i] = np.arctan2(w[j - 1].imag, w[j - 1].real) / (2 * np.pi)
        if nu1[i] < 0.0:
            nu1[i] += 1.0
        if nu1[i] <= 0.5:
            nu2[i] = nu1[i]
        else:
            nu2[i] = 1.0 - nu1[i]

    for i in range(n_dof):
        c = closest(nu2_M[i], nu2[0], nu2[1], nu2[2])
        if c != i + 1:
            j = c * 2 - 1
            k = i * 2 + 1
            swap_mat(v, j - 1, k - 1)
            swap_mat(v, j, k)
            swap(w, j - 1, k - 1)
            swap(w, j, k)
            swap(nu1, i - 1, c - 1)
            swap(nu2, i - 1, c - 1)
    for i in range(n_dof):
        if (0.5 - nu1_M[i]) * (0.5 - nu1[i]) < 0e0:
            j = i * 2 + 1
            swap_mat_imag(v, j - 1, j)
            swap_imag(w, j - 1, j)


def compute_A_inv(v: np.ndarray, *, n_dof=3, copy=True) -> (np.ndarray, np.ndarray):
    """Compute inverse matrix to matrix v

    What are the prerequistes of v?
    Should that be part of a packages that sounds if it were a
    linear algebra package?

    Warning:
        Current code is untested!

    Todo:
       Check if A should be always returned as a 6x6 matrix

    """
    if copy:
        v = _copy.copy(v)

    A_inv = np.identity(6)
    S = omega_block_matrix(n_dof)

    # Iterate over the degree of freedoms
    for dof in range(n_dof):
        i = dof * 2
        z = np.linalg.multi_dot([v[:, i].real, S, v[:, i].imag])
        sgn = np.sign(z)
        z = np.sqrt(np.absolute(1.0 / z))
        tmp = v[:, i].real + sgn * v[:, i].imag * 1j
        v[:, i] = z * tmp

    # Split complex representation into x, px planes and so on
    # This needs only limited data from v
    # * the elements v[i, i] with i = 0, 2, 4
    #   give a sign flip (accessed with array indices)
    #
    # * for each plane split up: corresponding v value
    #   into real and imaginary part as these will be stuffed to
    #   different places of the array
    #
    # Reworked Johan's code: need to refernce to his tech note
    # why it is done in this manner

    dof_idx = np.arange(n_dof) * 2
    indices = dof_idx[:, np.newaxis] * np.ones(n_dof, dtype=np.int_)[np.newaxis, :]
    indices = indices.ravel()
    logger.debug("compute_A_inv: indices %s", indices)

    signs = np.sign(v[indices, indices].real)

    n = 2 * n_dof
    for i in range(n):
        tmp = v[i, dof_idx]
        # Split up values in real and imaginary part: these will
        # be "stuffed" into different dimensions.
        tmp = [[t.real, t.imag] for t in tmp]
        # I use here that ravel runs fastest on the last index
        # (in this case) ... thus I assume that it will first runs
        # over real, imag  for the first value of v ... and so on
        tmp = np.ravel(tmp)
        # Only fill the requested columns
        A_inv[:n, i] = signs * tmp

    return A_inv, v


_delta = tslib.phase_space_index_internal.delta
_x  = tslib.phase_space_index_internal.x
_px = tslib.phase_space_index_internal.px
_ct = tslib.phase_space_index_internal.ct


def compute_A_inv_with_dispersion(
    v: np.ndarray, eta: np.ndarray, *, n_dof=3, copy=True
) -> (np.ndarray, np.ndarray):
    """

    Todo:
       Check: is the  dispersion implented as "small correction"
    """
    A_inv, v = compute_A_inv(v, n_dof=n_dof, copy=copy)

    # Small correction for dispersion ?
    B = np.identity(6)

    B[_x, _delta] = eta[_x]
    B[_px, _delta] = eta[_px]
    B[_ct, _x] = eta[_px]
    B[_ct, _px] = -eta[_x]

    A_inv = np.dot(A_inv, np.linalg.inv(B))
    return A_inv, v


def compute_A_inv_prev(n_dof, eta, v):
    """

    Original version thanks to Johan

    Todo:
       whant is special about this inverse
    """
    sign = np.sign

    n = 2 * n_dof

    # Should be a local copy.
    # v1 = np.array(v)
    v1 = _copy.copy(v)
    A_inv = np.identity(6)
    S = omega_block_matrix(n_dof)

    for i in range(n):
        if (i + 1) % 2:
            z = np.linalg.multi_dot([v1[:, i].real, S, v1[:, i].imag])
            sgn = sign(z)
            z = np.sqrt(np.abs(1e0 / z))
            for j in range(n):
                v1[j, i] = z * (v1[j, i].real + sgn * v1[j, i].imag * 1j)

    for i in range(n):
        A_inv[0, i] = sign(v1[0][0].real) * v1[i][0].real
        A_inv[1, i] = sign(v1[0][0].real) * v1[i][0].imag
        A_inv[2, i] = sign(v1[2][2].real) * v1[i][2].real
        A_inv[3, i] = sign(v1[2][2].real) * v1[i][2].imag
        if n > 4:
            A_inv[4, i] = sign(v1[4][4].real) * v1[i][4].real
            A_inv[5, i] = sign(v1[4][4].real) * v1[i][4].imag

    # Small correction for dispersion ?
    B = np.identity(6)
    delta_ = tslib.phase_space_index_internal.delta
    x_ = tslib.phase_space_index_internal.x
    px_ = tslib.phase_space_index_internal.px
    ct_ = tslib.phase_space_index_internal.ct

    B[x_, delta_], B[px_, delta_] = eta[x_], eta[px_]
    B[ct_, x_], B[ct_, px_] = eta[px_], -eta[x_]

    A_inv = np.dot(A_inv, np.linalg.inv(B))
    return A_inv, v1
