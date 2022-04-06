"""Functionallity for computing linear optics

"""
from thor_scsi.lib import Accelerator, ss_vect_tps, ConfigType
import xarray as xr
import numpy as np
from typing import Sequence
import copy

def compute_closed_orbit():
    """
    """

#def extract_twiss_parameters(A: ss_vect_tps) -> Sequence[[np.ndarray, np.ndarray, np.ndarray]]:

def extract_twiss_parameters(A: ss_vect_tps) -> Sequence:
    """

    numpy arrays returned as this function is only used internally
    """
    alpha = np.zeros(2, np.float) + 1
    beta = np.zeros(2, np.float) + 2
    nu = np.zeros(2, np.float) + 3

    return alpha, beta, nu


def compute_twiss_parameters(acc: Accelerator) -> xr.DataArray:
    """Computes Twiss parameters using TPSA

    Args:
        acc: an :class:`thor_scsi.lib.Accelerator` instance
    """
    # Used here for simplicity
    A = ss_vect_tps()
    A.set_identity()

    calc_config = ConfigType()

    along_acc = []
    for elem in acc:
        elem.propagate(calc_config, A)
        along_acc.append(A)

    twiss_tmp = [extract_twiss_parameters(A) for A in along_acc]

    idx = [elem.index for elem in acc]
    twiss_parameters = xr.DataArray(
        data=twiss_tmp,
        dims=["s", "twiss", "transverse"],
        coords=[idx, ["alpha", "beta", "nu"], ["x", "y"]],
    )
    # return twiss_parameters

    # alternative return

    # Lets keep reference to the actual elements
    t_lattice = [elem for elem in acc]
    # this is not working yet would be better:
    #  state of the elements when the twiss parameters were calculated
    # t_lattice = [copy.deepcopy(elem) for elem in acc]

    extra_info = xr.DataArray(data=t_lattice, dims=["s"], coords=[idx])
    names = xr.DataArray(data=[elem.name for elem in acc], dims=["s"], coords=[idx])

    res = xr.Dataset(dict(twiss_parameters=twiss_parameters, elements=extra_info, names=names))
    return res

__all__ = ["compute_closed_orbit", "compute_twiss_parameters"]
