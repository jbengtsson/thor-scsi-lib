"""
"""
from thor_scsi.lib import Accelerator
import xarray as xr
import numpy as np


def compute_path_length(acc: Accelerator) -> np.ndarray:
    """Add up the length of the individual elements of the accelerator
    """
    ds = [elem.get_length() for elem in acc]
    return np.add.accumulate(ds)


def accelerator_info(acc: Accelerator) -> xr.Dataset:
    """extra info on accelerator to be combined wit datasets returned by the accelerator

    Args:
        acc:         an acccelerator calculation engine
        calc_config: the c

    Returns:
        xr.Dataset
    """
    # stuff extra info in so that one can circle back to what the Twiss parameters
    # relay too

    # Lets keep reference to the actual elements
    t_lattice = [elem for elem in acc]

    # this is not working yet would be better:
    #  state of the elements when the twiss parameters were calculated
    # t_lattice = [copy.deepcopy(elem) for elem in acc]

    idx = np.array([elem.index for elem in acc])
    s = compute_path_length(acc)

    extra_info = xr.DataArray(data=t_lattice, dims=["index"], coords=[idx])
    names = xr.DataArray(data=[elem.name for elem in acc], dims=["index"], coords=[idx])
    s_da = xr.DataArray(data=s, dims=["index"], coords=[idx])
    res = xr.Dataset(dict(elements=extra_info, names=names, s=s_da))
    return res


__all__ = ["accelerator_info", "compute_path_length"]
