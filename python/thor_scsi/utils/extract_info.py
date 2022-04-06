"""
"""
from thor_scsi.lib import Accelerator
import xarray as xr


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

    idx = [elem.index for elem in acc]

    extra_info = xr.DataArray(data=t_lattice, dims=["s"], coords=[idx])
    names = xr.DataArray(data=[elem.name for elem in acc], dims=["s"], coords=[idx])

    res = xr.Dataset(dict(elements=extra_info, names=names))
    return res
