"""Functionallity for computing linear optics
"""
import thor_scsi.lib as tslib
from .phase_space_vector import ss_vect_tps2ps_jac, array2ss_vect_tps
from .courant_snyder import compute_A_CS
from .extract_info import accelerator_info
from .accelerator import instrument_with_standard_observers
from .phase_space_vector import map2numpy
from .linalg import compute_A_inv_with_dispersion, match_eigenvalues_to_plane_orig
from .output import mat2txt

import xarray as xr
import numpy as np
import copy
import logging


logger = logging.getLogger("thor_scsi")


def compute_dispersion(M: np.ndarray) -> np.ndarray:
    """

    Todo:
        Define which standard M should follow standard rules or
        thor_scsi / tracy internal rules
    """
    n = 4
    id_mat = np.identity(n)
    D = M[:n, tslib.phase_space_index_internal.delta]

    return np.linalg.inv(id_mat - M[:n, :n]) @ D


#: scale arctan2 (a12/a11) to Floquet coordinates (correct?)
scale_nu = 1.0 / (2 * np.pi)


def jac2twiss(A: np.ndarray) -> (float, float, float):
    """extract twiss parameters from array for one plane

    returns: alpha, beta, dnu
    See XXX eq (172)
    """
    eps = 1e-15

    A = np.atleast_2d(A)
    nr, nc = A.shape
    if nr != 2:
        raise AssertionError(f"expected 2 rows got {nr}")

    if nc != 2:
        raise AssertionError(f"expected 2 cols got {nc}")

    del nr, nc

    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1, 1]

    alpha = -(a11 * a21 + a12 * a22)
    beta = a11 ** 2 + a12 ** 2

    dnu = np.arctan2(a12, a11) * scale_nu
    # epsilon must be negative here ... ensure that the value is
    # not accidentially zero
    if dnu < -eps:
        dnu += 1.0

    return alpha, beta, dnu


def _extract_tps(elem: tslib.Observer) -> tslib.ss_vect_tps:
    """Extract tps data from the elments' observer
    """
    ob = elem.getObserver()
    assert ob.hasTruncatedPowerSeries()
    return ob.getTruncatedPowerSeries()


def transform_matrix_extract_twiss(A: np.ndarray) -> (np.ndarray, np.ndarray):
    """Extract twiss parameters from matrix rotating towards Floquet space

    Returns:
        eta, alpha, beta, nu

    """
    # njac, dnus = compute_A_CS(2, jac)
    tmp = [jac2twiss(A[p : p + 2, p : p + 2]) for p in range(0, 4, 2)]
    twiss = np.array(tmp)

    # only sensible for the first two planes ...
    delta = tslib.phase_space_index_internal.delta
    etas = A[:4, delta]

    return etas, twiss


def tps2twiss(tps: tslib.ss_vect_tps) -> (np.ndarray, np.ndarray):
    """Extract twiss information from turncated power series

    Warning:
        It is thue users responsibility to ensure that the matrix
        is a matrix rotationg towards to Floquet space
    """
    ps, Aj = ss_vect_tps2ps_jac(tps)
    return transform_matrix_extract_twiss(Aj)


def find_phase_space_origin(M: np.ndarray) -> np.ndarray:
    """Transform to energy dependent fix point

    * Diagonalise M = A R A^(-1)

    This gives:
    * Transfrom to the energy dependent fix point
    * Diagonalising in [x, px, y, py]

       R : Block diagonal
       :math:`A^(-1)` : From phase space ellipse to Floquet space circle

    * Transfrom from energy dependent fix point math::`\eta` math::`\delta` to
      origin of phase space


    Todo:
        Revisit naming: fix point instead of origin
        Add reference to Johan's Tech Note

    """

    n_dof = 2
    n = 2 * n_dof

    M_tp = np.transpose(M[:n, :n])

    # Finds eigen values and eigen vectors
    w, v = np.linalg.eig(M_tp)

    # Sort eigen values to match the trace of the plane
    # print(f"before sort w:{mat2txt(w)}  ")
    # print(f"before sort v:{mat2txt(v)}  ")
    w, v = match_eigenvalues_to_plane_orig(M_tp, w, v, n_dof=n_dof)
    # wabs = np.absolute(w)
    # print(f"sort  wabs   {mat2txt(wabs)}")
    # print(f"after sort w:{mat2txt(w)}  ")
    # print(f"after sort v:{mat2txt(v)}  ")

    eta = compute_dispersion(M)
    A_inv, v1 = compute_A_inv_with_dispersion(v, eta, n_dof=n_dof)
    A = np.linalg.inv(A_inv)
    A, _ = compute_A_CS(2, A)

    return A


def compute_map(
    acc: tslib.Accelerator, calc_config: tslib.ConfigType
) -> tslib.ss_vect_tps:
    """Propagate an identity map through the accelerator
    """
    t_map = tslib.ss_vect_tps()
    t_map.set_identity()
    # acc.propagate(calc_config, t_map, 0, len(acc))
    acc.propagate(calc_config, t_map)
    return t_map


def propagate_and_find_phase_space_orgin(acc, calc_config):
    """propagate once around ring. use this map to find phase space origin

    Todo:
         Rename phase space origin fix point
    """
    t_map = compute_map(acc, calc_config)
    M = map2numpy(t_map)

    Aj = find_phase_space_origin(M)
    Atmp = np.zeros([7, 7], dtype=np.float)
    Atmp[:6, :6] = Aj
    A = array2ss_vect_tps(Atmp)

    return A


def compute_twiss_along_lattice(
    acc: tslib.Accelerator,
    calc_config: tslib.ConfigType = None,
    A: tslib.ss_vect_tps = None,
) -> xr.Dataset:
    """

    Args:
        acc :         an :class:`tlib.Accelerator` instance
        calc_config : an :class:`tlib.Config` instance
        A : see :func:`compute_ring_twiss` for output

    returns xr.Dataset

    Todo:
        Consider if tps should be returned too
        Split up functionallity (not everyone wants output as
        xarrays. But what then to use?)

    """
    if not calc_config:
        calc_config = tslib.ConfigType()

    if A is None:
        A = propagate_and_find_phase_space_orgin(acc, calc_config)

    # Not really required ... but used for convienience
    instrument_with_standard_observers(acc)

    # Propagate through the accelerator
    for k in range(len(acc)):
        acc.propagate(calc_config, A, k, k + 1)
        # Extract the first order term from the map representation of A
        # thus called Aj as it is similar to the Jacobian
        _, Aj = ss_vect_tps2ps_jac(A)
        # A will be rotated ... so rotate it back so that the
        # next step starts at a Courant Snyder form
        rjac, _ = compute_A_CS(2, Aj)
        Atmp = np.zeros([7, 7], dtype=np.float)
        Atmp[:6, :6] = rjac
        A = array2ss_vect_tps(Atmp)

    indices = [elem.index for elem in acc]
    tps_tmp = [_extract_tps(elem) for elem in acc]
    data = [tps2twiss(t) for t in tps_tmp]
    tps_tmp = np.array(tps_tmp)
    twiss_pars = [d[1] for d in data]
    disp = [d[0] for d in data]

    # Stuff information into xarrays ..
    dispersion = xr.DataArray(
        data=disp,
        name="dispersion",
        dims=["index", "phase_coordinate"],
        coords=[indices, ["x", "px", "y", "py"]],
    )
    twiss_parameters = xr.DataArray(
        twiss_pars,
        name="twiss",
        dims=["index", "plane", "par"],
        coords=[indices, ["x", "y"], ["alpha", "beta", "dnu"]],
    )
    phase_space_coords_names = ["x", "px", "y", "py", "ct", "delta"]
    tps = xr.DataArray(
        data=tps_tmp,
        name="tps",
        dims=["index", "phase_coordinate"],
        coords=[indices, phase_space_coords_names],
    )
    info = accelerator_info(acc)
    res = info.merge(dict(twiss=twiss_parameters, dispersion=dispersion, tps=tps))
    return res


__all__ = ["compute_twiss_along_lattice", "jac2twiss"]
