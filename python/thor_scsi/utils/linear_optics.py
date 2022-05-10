"""Functionallity for computing linear optics

"""
import thor_scsi.lib as tslib
from .phase_space_vector import ss_vect_tps2ps_jac, array2ss_vect_tps
from .courant_snyder import compute_A_CS
from .extract_info import accelerator_info
from .accelerator import instrument_with_standard_observers
from .output import mat2txt

import xarray as xr
import numpy as np
import copy


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

    dnu_s = 1.0/(2 * np.pi)
    dnu = np.arctan2(a12, a11) * dnu_s
    if dnu < eps:
        # dnu += 1.0
        pass

    return alpha, beta, dnu


def _extract_tps(elem: tslib.Observer) -> tslib.ss_vect_tps:
    """Extract tps data from the elments' observer
    """
    ob = elem.getObserver()
    assert ob.hasTruncatedPowerSeries()
    return ob.getTruncatedPowerSeries()


def jac_extract_twiss(jac: np.ndarray) -> (np.ndarray, np.ndarray):
    """

    Rename:
       jac to A (A is transformation to Floquet space)
    """
    # njac, dnus = compute_A_CS(2, jac)
    njac = jac
    tmp = [jac2twiss(njac[p: p + 2, p: p + 2]) for p in range(0, 4, 2)]
    # Only sensible for the first two planes ...
    delta = tslib.phase_space_index_internal.delta
    etas = njac[:4, delta]
    twiss = np.array(tmp)
    # disp = compute_dispersion(jac)
    return etas, twiss


def tps2twiss(
    tps: tslib.ss_vect_tps
) -> (np.ndarray, np.ndarray):
    """Extract information from turncated power series

    returns eta, alpha, beta, nu
    """
    ps, jac = ss_vect_tps2ps_jac(tps)
    return jac_extract_twiss(jac)


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
    """
    if not calc_config:
        calc_config = tslib.ConfigType()

    assert(A is not None)
    if A is None:
        # Should not be here
        A = tslib.ss_vect_tps()
        A.set_identity()

    # Not really required
    instrument_with_standard_observers(acc)

    # Propagate through the accelerator
    l = []
    indices = []
    for k in range(len(acc)):
        acc.propagate(calc_config, A, k, k+1)
        _, jac = ss_vect_tps2ps_jac(A)
        # A will be rotated
        # phase in A will be 0
        # print(f"jac at {k}\n{mat2txt(jac)}")
        njac = copy.copy(jac)
        l.append(njac)
        # insert compute_twiss
        # eta, alpha, beta, nu = compute_twiss_A(A)

        rjac, _ = compute_A_CS(2, njac)
        Atmp = np.zeros([7, 7], dtype=np.float)
        Atmp[:6, :6] = rjac
        # Atmp[6, :] = ps
        A = array2ss_vect_tps(Atmp)
        indices.append(acc[k].index)
        if False and k in [10, 11, 12, 13, 14, 15]:
            delta = tslib.phase_space_index_internal.delta
            a_mat_ref = np.load("A_mat_save.npy")
            print(f"jac at {k, acc[k].name}\n{mat2txt(jac)}\n etas {jac[:4, delta]}")
            print(f"rotated jac at {k}\n{mat2txt(rjac)}\n")
            print("A (map)")
            print(A)
            diff = njac - a_mat_ref
            print(f"chck at {k}\n{mat2txt(diff)}\n")
        #if k == 15:
        #    break
        # break

    # Retrieve information
    if False:
        tps_tmp = [_extract_tps(elem) for elem in acc]
        tps_tmp = np.array(tps_tmp)
        phase_space_coords_names = ["x", "px", "y", "py", "ct", "delta"]
        tps = xr.DataArray(
            data=tps_tmp,
            name="tps",
            dims=["index", "phase_coordinate"],
            coords=[indices, phase_space_coords_names],
        )
        data = [tps2twiss(t) for t in tps_tmp]
    else:
        data = [jac_extract_twiss(j) for j in l]
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
        coords=[indices,  ["x", "y"], ["alpha", "beta", "dnu"]],
    )
    info = accelerator_info(acc)
    res = info.merge(dict(twiss=twiss_parameters, dispersion=dispersion,
                          # tps=tps
    ))
    return res


__all__ = ["compute_twiss_along_lattice", "jac2twiss"]
