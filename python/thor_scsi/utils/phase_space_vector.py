"""Various utils for phase space vector
"""
import thor_scsi.lib as tslib
import numpy as np


def omega_block_matrix(n_dim: int) -> np.ndarray:
    """create omega block matrix required for testing symplectisism

    Args:
       n_dim : dimension (e.g. 3 for a 6 dimensional phase space)
    """
    n = n_dim * 2
    a = np.zeros([n, n])

    # That would be the standard procedure
    #
    # for i in range(0, n, 2):
    #     j = i + 1
    #     a[i, j] = 1
    #     a[j, i] = -1

    # With advance indexing
    idx = np.arange(0, n, 2)
    a[idx, idx + 1] = 1
    a[idx + 1, idx] = -1
    return a


def jac_intern_to_std(a: np.ndarray, copy: bool = True) -> np.ndarray:
    """Swap the intern representation to the more standard one

    thor_scsi tps places use ct as last row and column.
    Standard praxis is to place it as last but one row / column.
    Then the matrix is block symmetric. This functions swaps the
    rows and columns

    """
    n = 6
    a = np.atleast_2d(a)
    nr, nc = a.shape
    if nr != n:
        raise AssertionError("Number row rows {nr} != {n}")
    if nc != n:
        raise AssertionError("Number of cols {nc} != {n}")

    if copy:
        a = a.copy()

    # Swap rows
    a[[4, 5], :] = a[[5, 4], :]
    # Swap cols
    a[:, [4, 5]] = a[:, [5, 4]]

    return a


def ps_intern_to_std(a: np.ndarray, copy: bool = True) -> np.ndarray:
    """Swap phase space vector form internal to more standard representation

    Similar to :func:`jac_intern_to_std`
    """
    n = 7
    a = np.atleast_1d(a)
    nr, = a.shape
    if nr != n:
        raise AssertionError(f"Number row rows {nr} != {n}")
    if copy:
        a = a.copy()

    # Swap elements
    a[[4, 5]] = a[[5, 4]]
    return a


def get_mat(t_map):
    mat = tslib.ss_vect_tps_to_mat(t_map)
    mat = np.array(mat)
    return mat


def ss_vect_tps2ps_jac(
    ps: tslib.ss_vect_tps, std: bool = True
) -> (np.ndarray, np.ndarray):
    """Extract phase space and jacobian from internal ss_vect_tps representation

    Args:
       std: output as standard array (i.e. ct in 4th and delta in 5th column)
    """
    raw = tslib.ss_vect_tps_to_mat(ps)
    tmp = np.array(raw)
    n_jac = 6
    jac = tmp[:n_jac, :n_jac]
    ps = tmp[n_jac, :]

    if std:
        jac = jac_intern_to_std(jac)
        ps = ps_intern_to_std(ps)
    return ps, jac


def array2ss_vect_tps(a_mat: np.ndarray) -> tslib.ss_vect_tps:
    mat = tslib.Matrix(a_mat)
    t_map = tslib.mat_to_ss_vect_tps(mat)
    return t_map


def ps_jac2ss_vect_tps(ps: np.ndarray, jac: np.ndarray) -> tslib.ss_vect_tps:
    """
    """
    tmp = np.zeros([7, 7], np.float)
    tmp[:6, :6] = jac
    tmp[6, :6] = ps
    return array2ss_vect_tps(tmp)


__all__ = [
    "array2ss_vect_tps",
    "jac_intern_to_std",
    "omega_block_matrix",
    "ps_intern_to_std",
    "ps_jac2ss_vect_tps",
    "ss_vect_tps2ps_jac",
]
