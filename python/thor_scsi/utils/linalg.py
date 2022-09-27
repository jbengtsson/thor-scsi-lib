"""Eigen for symplectic matrices ?

"""
import numpy as np
from typing import Sequence
import copy as _copy
import thor_scsi.lib as tslib
import logging
import functools
from .phase_space_vector import omega_block_matrix as _obm
from .math import minimum_distance_above_threshold
from dataclasses import dataclass

from .courant_snyder import compute_A_CS
from thor_scsi.utils.output import vec2txt, mat2txt

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

    What about argsort (np.absolute(vec - x))
    """

    dx1, dx2, dx3 = [np.abs(x - tx) for tx in (x1, x2, x3)]
    if dx1 < dx2 and dx1 < dx3:
        k = 1
    elif dx2 < dx1 and dx2 < dx3:
        k = 2
    else:
        k = 3
    return k


class PlaneDistance:
    """Compute distance to different planes and compare to other

    Args:
       value:             tune attributed to this plane
       reference_values:  tunes derived from the eigen vectors
       payload:           an object typically the associated eigen vector

    Intended use case: matching eigenvectors to the different planes

    I guess the reference values should be sorted ....
    Just for corner cases
    """

    def __init__(self, *, value, reference_values, payload):
        self.value = value
        self.reference_values = reference_values
        self.payload = payload

        # What's the metric ... distance of this value to the reference values
        # Well here Mannheim ehh Manhattan distance
        self.distances = np.absolute(value - reference_values)
        self.min_distance = np.argmin(self.distances)

    def __lt__(self, other):
        if self.min_distance != other.min_distance:
            return self.min_distance < other.min_distance

        # Both would match to the same plane ... what to do now ...
        # which is closer
        t_d = self.distances[self.min_distance]
        o_d = other.distances[other.min_distance]

        return t_d < o_d

    def __eq__(self, other):
        # Different planes identified ... fine
        if self.min_distance != other.min_distance:
            return False

        # Both would match to the same plane ... what to do now ...
        # which is closer
        t_d = self.distances[self.min_distance]
        o_d = other.distances[other.min_distance]
        return t_d == o_d

    def __repr__(self):
        cls_name = self.__class__.__name__
        v = self.value
        md = self.min_distance
        ds = self.distances
        start = self.payload.start
        txt = f"{cls_name}(start={start},value={v},min_distance={md},distances={ds})"
        return txt


@dataclass
class PlaneInfo:
    #: phase adavance as deduced from Poincare map
    nu1: float
    #: phase adavance as deduced from Poincare map: different representation
    nu2: float

    #: phase advance as deduced from transformation to Floquet space
    w: float
    #: corresponding eigen vector
    v: np.ndarray

    # original index
    start: int


def match_eigenvalues_to_plane(M: np.ndarray, w: np.ndarray, v: np.ndarray, *, n_dof):
    """Sort eigen values to be in order x, y, longitudinal

    Args:
        M: Poincare Matrix
        w: eigen values
        v: eigen vectors

    Todo:
        can n_dof be derived from the shape of M, w or v?
        Check that shape of M, w, v correspond
        Check it on the Poincare monday

    Warning:
        Untested code ... needs to be checked further
        Will only check if some values are quite similar
    """
    sign = np.sign

    # check that it is a matrix
    # Should we check for quadratic matrix ?
    M = np.atleast_2d(M)
    nr, nc = M.shape
    if nr < n_dof:
        raise AssertionError(f"Poincare map: {nr=} < {n_dof=}")
    if nc < n_dof:
        raise AssertionError(f"Poincare map: {nc=} < {n_dof=}")

    # vector of eigen values
    w = np.atleast_1d(w)
    nr, = w.shape
    if nr < n_dof:
        raise AssertionError(f"Poincare map: {nr=} < {n_dof=}")

    # eigen vectors stored in a matrix
    v = np.atleast_2d(v)
    nr, nc = v.shape
    if nr < n_dof:
        raise AssertionError(f"Poincare map: {nr=} < {n_dof=}")
    if nc < n_dof:
        raise AssertionError(f"Poincare map: {nc=} < {n_dof=}")

    # Why these planes now ...
    indices = np.arange(n_dof) * 2 + 1

    cos_M = (M[indices - 1, indices - 1] + M[indices, indices]) / 2

    # Analyse stabilites of the plane ... questionable if one could sort them ?
    # Would that give us insight which plane is stable ?
    # Should one add an epsilon to 1?
    unstable = np.absolute(cos_M) > 1

    if unstable.any():
        n_unstable = unstable.sum()
        vals = cos_M[unstable]
        vals = np.absolute(vals) - 1
        vals_txt = ["{:10.3e}".format(v) for v in vals]
        txt = (
            f"Matching eigenvalues to plane: found {n_unstable:d} unstable values:"
            f" |cos_M[nu]-1e0| = {vals_txt}"
        )
        logger.warning(txt)
        del n_unstable, vals, vals_txt, txt

    # Use the phase advance to match the eigen vectors to the planes
    # In a first step compute the phase advance for the Poincare mao
    sin_M = sign(M[indices - 1, indices]) * np.sqrt(1e0 - cos_M ** 2)
    nu_M = np.arctan2(sin_M, cos_M) / (2 * np.pi)

    # Let's check that no two nu's are too close
    threshold = 0.01
    flag, idx1, idx2 = minimum_distance_above_threshold(nu_M, threshold)
    if not flag:
        logger.error(f"Investigating phase advance (from Poincare Map: {nu_M}")
        txt = f"Phase advances [{idx1}, {idx2}] too close {nu_M[[idx1,idx2]]}"
        raise AssertionError(txt)

    def arange_nus(nus: np.ndarray) -> np.ndarray:
        nu1 = nus.copy()
        nu1[nu1 < 0] += 1
        # sorted from -.5 ... .5 ?
        # Find out the sign
        nu2 = np.where(nu1 < 0.5, nu1, 1.0 - nu1)
        return nu1, nu2

    nu1_M, nu2_M = arange_nus(nu_M)
    del nu_M

    # Now deduce the eigenvalues from the vectors
    nu = np.arctan2(w[indices - 1].imag, w[indices - 1].real) / (2 * np.pi)
    nu1, nu2 = arange_nus(nu)
    del nu

    # I am not sure if I need to sort nu2 ...
    # Somehow I think it is a more stable approach ?
    # currently not done: it is assumed that the different phase
    # advances are far enough apart that the sorting will not take
    # into account that two values have to be inspected at once to
    # find the more probable one
    sorters = [
        PlaneDistance(
            value=ph,
            reference_values=nu2,
            payload=PlaneInfo(nu1=nu1[i], nu2=nu2[i], v=v[:, i], w=w[i], start=i),
        )
        for i, ph in enumerate(nu2_M)
    ]
    # Just to be sure that there is no mess left ....
    del nu1, nu2# , v, w
    sorters.sort()
    logger.warning(f"sorting result {sorters}")
    logger.warning(f"eigenvector values {nu2_M}")

    # construct w and v in the proper order
    wr = [s.payload.w for s in sorters]
    vr = [s.payload.v for s in sorters]
    nu1 = [s.payload.nu1 for s in sorters]

    w = w.copy()
    v = v.copy()

    # Find out which have to be rearranged ...
    # need to understand that ...
    # assert 0

    w = np.array(w)
    # I think I have to transpose this array to get it back to the original
    # dimension order
    v = np.array(v)
    v = v.T

    # derive the correct sign
    # become pairs of [lambda_ k, 1 / (lambda_ k)], in the common metric ?
    # I think that should be done using the data stored in the payload
    # or in the planeinfo class
    for i in range(n_dof):
        if (0.5 - nu1_M[i]) * (0.5 - nu1[i]) < 0e0:
            j = i * 2 + 1
            swap_mat_imag(v, j - 1, j)
            swap_imag(w, j - 1, j)

    return w, v


def match_eigenvalues_to_plane_orig(M, w, v, *, n_dof):
    """Sort eigen values to be in order x, y, longitudinal

    Args:
        M: Poincare Matrix
        w: eigenvalues
        v: eigen vectors

    Todo:
        can n_dof be derived from the shape of M, w or v?
        Check that shape of M, w, v correspond

        Check reuse of symbol "v":
            * in the 1. for loop v is assinged to a scalar
            * 2. and 3. for loop v is not reassigned but seems
              to be used as vector or matrix ....

    Warning:
        raise an exception if a non stable solution was find
    """

    sign = np.sign
    n = 2 * n_dof

    v = v.copy()
    w = w.copy()

    w_ret = w
    # v is reassigned ... look to todo to the function description
    v_ret = v

    sin_M, cos_M = np.zeros(n_dof), np.zeros(n_dof)
    nu1_M, nu2_M, nu1, nu2 = np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

    for i in range(n_dof):
        j = (i + 1) * 2 - 1
        # Is this renaming of v intentional ? v used as a scalar
        try:
            vtmp = (M[j - 1][j - 1] + M[j][j]) / 2
        except IndexError as ie:
            logger.error(f"{i=}, {j=}")
            raise
        cos_M[i] = vtmp
        if np.abs(vtmp) > 1e0:
            txt = f"sort_eigen: unstable |cos_M[nu_{i+1}]-1e0| = {np.absolute(vtmp - 1e0):10.3e}"
            logger.warning(txt)
            stable = False
            # Check if an exception should be raised here?
            raise ValueError(txt)

        sin_M[i] = sign(M[j - 1][j]) * np.sqrt(1e0 - cos_M[i] ** 2)
        #
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

    # v used as matrix
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
            try:
                swap_mat_imag(v, j - 1, j)
            except Exception as exc:
                logger.error(f"v {v}")
                raise exc
            swap_imag(w, j - 1, j)

    return w_ret, v_ret


def compute_A_inv(v: np.ndarray, *, n_dof=3, copy=True) -> (np.ndarray, np.ndarray):
    """Compute inverse matrix to matrix v

    What are the prerequistes of v?
    Should that be part of a packages that sounds if it were a
    linear algebra package?

    Warning:
        Current code is untested!

    Todo:
       Check if A should be always returned as a 6x6 matrix
       Check if special treatment is required due to the 2.5 dimensional
       ordering of phase space dimensions
    """
    if copy:
        v = v.copy()

    A_inv = np.identity(6)
    S = omega_block_matrix(n_dof)

    # Iterate over the degree of freedoms
    for dof in range(n_dof):
        i = dof * 2
        z = [v[:, i].real @ S @ v[:, i].imag]
        # sgn = np.sign(z)
        # tmp = v[:, i].real + sgn * v[:, i].imag * 1j
        zs = np.sqrt(np.absolute(1.0 / z))
        tmp = v[:, i]
        if z < 0:
            tmp = tmp.conj()
        v[:, i] = zs * tmp

    # Split complex representation into x, px planes and so on
    # This needs only limited data from v
    # * the elements v[i, i] with i = 0, 2, 4
    #   give a sign flip (accessed with array indices)
    #
    # * for each plane split up: corresponding v value
    #   into real and imaginary part as these will be stuffed to
    #   different places of the array
    #
    #
    # Reworked Johan's code: need to reference to his tech note
    # why it is done in this manner

    dof_idx = np.arange(n_dof) * 2
    indices = dof_idx[:, np.newaxis] * np.ones(n_dof, dtype=np.int_)[np.newaxis, :]
    indices = indices.ravel()
    logger.debug("compute_A_inv: indices %s", indices)

    signs = np.sign(v[indices, indices].real)

    n = 2 * n_dof
    for i in range(n):
        # Use complex data for each plane (stored in plane 0, 2 and 4)
        tmp = v[i, dof_idx]
        # Split up values in real and imaginary part: these will
        # be "stuffed" to the different phase space dimensions
        # x, px, y, py, ...
        tmp = [[t.real, t.imag] for t in tmp]
        # I use here that ravel runs fastest on the last index
        # (in this case) ... thus I assume that it will first runs
        # over real, imag  for the first value of v ... and so on
        tmp = np.ravel(tmp)
        # Only fill the requested columns
        A_inv[:n, i] = signs * tmp

    return A_inv, v


_delta = tslib.phase_space_index_internal.delta
_x = tslib.phase_space_index_internal.x
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


def compute_A_prev(n_dof, eta, u):
    """

    Original version thanks to Johan

    Todo:
       what's is special about this inverse
    """
    sign = np.sign

    n = 2 * n_dof

    # Should be a local copy.
    # u1 = np.array(u)
    u1 = _copy.copy(u)
    A = np.identity(6)
    S = omega_block_matrix(n_dof)

    # Normalise eigenvectors: A^T.omega.A = omega.
    for i in range(n_dof):
        z = u1[:, 2*i].real @ S @ u1[:, 2*i].imag
        sgn_im = sign(z)
        scl = np.sqrt(np.abs(z))
        sgn_vec = sign(u1[2*i][2*i].real)
        [u1[:, 2*i], u1[:, 2*i+1]] = \
            [sgn_vec * (u1[:, 2*i].real + sgn_im * u1[:, 2*i].imag * 1j) / scl,
             sgn_vec * (u1[:, 2*i+1].real + sgn_im * u1[:, 2*i+1].imag * 1j)
             / scl]

    for i in range(n_dof):
        [A[:n, 2*i], A[:n, 2*i+1]] = [u1[:, 2*i].real, u1[:, 2*i].imag]

    if n_dof == 2:
        # If coasting beam translate to momentum dependent fix point.
        delta_ = tslib.phase_space_index_internal.delta
        x_     = tslib.phase_space_index_internal.x
        px_    = tslib.phase_space_index_internal.px
        ct_    = tslib.phase_space_index_internal.ct

        B = np.identity(6)
        B[x_, delta_], B[px_, delta_] = eta[x_], eta[px_]
        B[ct_, x_], B[ct_, px_]       = eta[px_], -eta[x_]

        A = B @ A

    return A, u1
