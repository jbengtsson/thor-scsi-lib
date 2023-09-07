"""Utils for calculating with the accelerator

Todo:
    check if Facade should be created hiding all the functionallity
"""
import gtpsa
import thor_scsi.lib as tslib
from .extract_info import accelerator_info
import numpy as np
import xarray as xr
from typing import Sequence
import logging

logger = logging.getLogger("thor-scsi")


def _instrument_sequence_with_radiation_delegate(
    elems: Sequence[tslib.ElemType], *, delegate
) -> Sequence:
    """Instrument a sequence of elements with object delegates"""
    rad_del = [delegate() for elem in elems]
    for a_del, elem in zip(rad_del, elems):
        elem.set_radiation_delegate(a_del)
    return rad_del


def instrument_with_standard_radiators(
    acc: tslib.Accelerator, *, energy
) -> Sequence[tslib.RadiationDelegate]:
    """Instrument a sequence of elements with radiation delegates

    Todo:
      are these rather observers for curly_H?
    """
    rad_del = _instrument_sequence_with_radiation_delegate(
        acc, delegate=tslib.RadiationDelegate
    )
    return rad_del


def instrument_with_standard_radiator_kickers(
    acc: tslib.Accelerator, *, energy
) -> Sequence[tslib.RadiationDelegate]:
    """Instrument all resonable elements with a radiation delegate

    Todo:
        review if a radation delegate should be registered to any
        element unless it refuses to accept one
    """
    elems = [elem for elem in acc if isinstance(elem, (tslib.Mpole,))]
    logger.info(
        "\nRadiators added to:\n({})".format(", ".join([e.name for e in elems]))
    )
    rad_del_kick = _instrument_sequence_with_radiation_delegate(
        elems, delegate=tslib.RadiationDelegateKick
    )
    for a_del in rad_del_kick:
        a_del.set_energy(energy)
    return rad_del_kick


def instrument_with_standard_observers(
    elems: Sequence[tslib.ElemType], *, mapping=None
) -> Sequence[tslib.StandardObserver]:
    """Instrument accelerator with observers

    Returns accelerator list of created observers
    """
    # Instrument it with observers ... I guess I have to keep them here so that
    # internally these are not weak references ... to be tested
    assert mapping
    observers = [tslib.StandardObserver(mapping) for elem in elems]
    # Register them to the elments
    for ob, elem in zip(observers, elems):
        elem.set_observer(ob)
    return observers


def extract_orbit_from_standard_observers(
    observers: Sequence[tslib.StandardObserver], info: xr.Dataset
) -> xr.Dataset:
    """ """
    observers = np.array(observers)

    # Required later to index into the final xarray
    has_observer = np.array([True if ob else False for ob in observers])

    # The ones that are not null
    ob_sub = observers[has_observer]

    # check that all have seen data
    with_observer = np.array([ob.has_truncated_power_series_a() for ob in ob_sub])
    lf = np.sum(with_observer)
    lo = len(ob_sub)
    if lf != lo:
        raise AssertionError(f"Have {lo} observers but only {lf} has (valid) data")

    indices = np.arange(len(with_observer))
    indices = indices[with_observer == True]
    phase_space_coords_names = ["x", "px", "y", "py", "ct", "delta"]
    tps_tmp = [ob.get_truncated_power_series_a() for ob in ob_sub]
    logger.debug(tps_tmp)
    tps = xr.DataArray(
        data=[t.jacobian() for t in tps_tmp],
        name="tps_jacobian",
        dims=["index", "phase_coordinate_row", "phase_coordinate_col"],
        coords=[indices, phase_space_coords_names, phase_space_coords_names],
    )
    ps = xr.DataArray(
        data=[np.array(t.cst().iloc) for t in tps_tmp],
        name="ps",
        dims=["index", "phase_coordinate"],
        coords=[indices, phase_space_coords_names],
    )
    info = info.isel(index=(with_observer == True))
    res = info.merge(dict(tps=tps, ps=ps))
    return res


def extract_orbit_from_accelerator_with_standard_observers(
    acc: Sequence[tslib.ElemType],
) -> xr.Dataset:
    """ """

    observers = [elem.get_observer() for elem in acc]
    info = accelerator_info(acc)
    return extract_orbit_from_standard_observers(observers, info)


__all__ = [
    "instrument_with_radiators",
    "instrument_with_standard_observers",
    "extract_orbit_from_standard_observers",
    "extract_orbit_from_accelerator_with_standard_observers",
]
