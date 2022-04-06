"""Functionallity for computing linear optics

"""
from thor_scsi.lib import Accelerator, ss_vect_tps, ConfigType, ss_vect_tps_to_mat
import xarray as xr
import numpy as np
from typing import Sequence
import copy


def compute_closed_orbit():
    """
    """


# def ss_vect_tps_to_ps_jac(ps: ss_vect_tps) -> Sequence[[np.ndarray, np.ndarray]]:


def ss_vect_tps2ps_jac(ps: ss_vect_tps) -> Sequence:
    """Extract phase space and jacobian from internal ss_vect_tps representation

    Todo:
        place the function into the appropriate module
    """
    raw = ss_vect_tps_to_mat(ps)
    tmp = np.array(raw)
    n_jac = 6
    jac = tmp[:n_jac, :n_jac]
    ps = tmp[n_jac, :]
    return ps, jac


# def extract_twiss_parameters(A: ss_vect_tps) -> Sequence[[np.ndarray, np.ndarray, np.ndarray]]:


def extract_twiss_parameters(A: ss_vect_tps) -> Sequence:
    """extract twiss parameters form ss_vect_tps representations

    Warning:
        Currently code just for demonstration purposes
    """

    ps, jac = ss_vect_tps2ps_jac(A)

    # Just to see data
    alpha = np.array(ps[:2])
    beta = np.array(ps[2:4])
    nu = np.array(ps[4:6])

    return alpha, beta, nu


def compute_twiss_parameters(
    acc: Accelerator, calc_config: ConfigType = None
) -> xr.DataArray:
    """Computes Twiss parameters using TPSA

    Args:
        acc: an :class:`thor_scsi.lib.Accelerator` instance
    """
    # Used here for simplicity
    A = ss_vect_tps()
    A.set_identity()

    if not calc_config:
        calc_config = ConfigType()

    along_acc = []
    for elem in acc:
        elem.propagate(calc_config, A)
        along_acc.append(copy.copy(A))

    twiss_tmp = [extract_twiss_parameters(A) for A in along_acc]
    twiss_tmp = np.array(twiss_tmp)

    # Computation is done here .... start to resort data to make it
    # easier to access an individual datum

    # store the twiss paramters in a dataarray
    idx = [elem.index for elem in acc]
    twiss_parameters = xr.DataArray(
        data=twiss_tmp,
        dims=["s", "twiss", "transverse"],
        coords=[idx, ["alpha", "beta", "nu"], ["x", "y"]],
        name="twiss_parameters",
    )
    return twiss_parameters


__all__ = ["compute_closed_orbit", "compute_twiss_parameters", "ss_vect_tps2ps_jac"]
