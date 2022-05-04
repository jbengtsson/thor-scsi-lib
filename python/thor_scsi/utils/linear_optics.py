"""Functionallity for computing linear optics

"""
import thor_scsi.lib as tslib
from .phase_space_vector import ss_vect_tps2ps_jac
from .extract_info import accelerator_info
from .accelerator import instrument_with_standard_observers
import xarray as xr
import numpy as np


def compute_dispersion(M: np.ndarray) -> np.ndarray:
    """

    Todo:
        Define which standard M should follow standard rules or
        thor_scsi / tracy internal rules
    """
    n = 4
    id_mat = np.identity(n)
    D = M[:n, tslib.phase_space_index_internal.delta]

    return np.dot(np.linalg.inv(id_mat - M[:n, :n]), D)


def jac2twiss(A: np.ndarray) -> (float, float, float):
    """extract twiss parameters from array for one plane

    returns: alpha, beta, dnu
    See XXX eq (172)
    """
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

    alpha = -a11 * a21 - a12 * a22
    beta = a11 ** 2 + a12 ** 2
    dnu = np.arctan2(a12, a11)

    return alpha, beta, dnu


def _extract_tps(elem: tslib.Observer) -> tslib.ss_vect_tps:
    """Extract tps data from the elments' observer
    """
    ob = elem.getObserver()
    assert ob.hasTruncatedPowerSeries()
    return ob.getTruncatedPowerSeries()


def tps2twiss(
    tps: tslib.ss_vect_tps
) -> (np.ndarray, np.ndarray):
    """

    returns eta, alpha, beta, nu
    """
    ps, jac = ss_vect_tps2ps_jac(tps)

    tmp = [jac2twiss(jac[p: p + 2, p: p + 2]) for p in range(2)]
    twiss = np.array(tmp)
    disp = compute_dispersion(jac)
    return disp, twiss


def compute_twiss_along_lattice(
    acc: tslib.Accelerator,
    calc_config: tslib.ConfigType = None,
    A: tslib.ss_vect_tps = None,
) -> xr.Dataset:
    """

    Args:
        acc : an :class:`tlib.Accelerator` instance
        calc_config : an :class:`tlib.Config` instance
        A :


    returns xr.Dataset

    Todo:
        Consider if tps should be returned too
    """
    if not calc_config:
        calc_config = tslib.ConfigType()

    if A is None:
        A = tslib.ss_vect_tps()
        A.set_identity()

    instrument_with_standard_observers(acc)
    # Propagate through the accelerator
    acc.propagate(calc_config, A)

    # Retrieve information
    tps_tmp = [_extract_tps(elem) for elem in acc]
    data = [tps2twiss(t) for t in tps_tmp]
    twiss_pars = [d[1] for d in data]
    disp = [d[0] for d in data]
    indices = [elem.index for elem in acc]

    tps_tmp = np.array(tps_tmp)

    # Stuff information into xarrays ..
    phase_space_coords_names = ["x", "px", "y", "py", "ct", "delta"]

    tps = xr.DataArray(
        data=tps_tmp,
        name="tps",
        dims=["index", "phase_coordinate"],
        coords=[indices, phase_space_coords_names],
    )
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
        coords=[indices,  ["x", "y"], ["alpha", "beta", "dnu"]],
    )
    info = accelerator_info(acc)
    res = info.merge(dict(twiss=twiss_parameters, dispersion=dispersion, tps=tps))
    print(res)
    return res


__all__ = ["compute_twiss_along_lattice", "jac2twiss"]
