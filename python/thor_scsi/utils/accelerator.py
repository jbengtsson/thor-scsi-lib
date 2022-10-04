"""Utils for calculating with the accelerator

Todo:
    check if Facade should be created hiding all the functionallity
"""
import thor_scsi.lib as tslib
from .extract_info import accelerator_info
import numpy as np
import xarray as xr
from typing import Sequence
import logging

logger = logging.getLogger("thor-scsi")


## def instrument_sequence_with_standard_radiators(
##     elems: Sequence[tslib.ElemType]
## ) -> Sequence[tslib.RadiationDelegate]:
##     """Instrument a sequence of elements with radiation delegates
##     """
##     ps_zero = tslib.ss_vect_double()
##     rad_del = [tslib.RadiationDelegate() for elem in elems]
##     for a_del, elem in zip(rad_del, elems):
##         elem.setRadiationDelegate(a_del)
##         # Just use that that the marker knows who is calling him
##         a_del.view(elem, ps_zero, tslib.ObservedState.start, 0)
##     return rad_del


def instrument_with_radiators(
    acc: tslib.Accelerator, *, energy
) -> Sequence[tslib.RadiationDelegate]:
    """Instrument all resonable elements with a radiation delegate

    Todo:
        review if a radation delegate should be registered to any
        element unless it refuses to accept one
    """

    # I think anything derived from a mpole can tak an radiation delegate
    ps_zero = tslib.ss_vect_double()

    ## Untested code below
    ## for type_name in ["Marker"]:
    ##     elems = [elem for elem in acc.elementsWithNameType(type_name)]
    ##
    ## # Markers and similar devices solely store the data ...
    ## rad_del = [tslib.RadiationDelegate() for elem in elems]
    ## for a_del, elem in zip(rad_del, elems):
    ##     elem.setRadiationDelegate(a_del)
    ##     # Just use that that the marker knows who is calling him
    ##     a_del.view(elem, ps_zero, tslib.ObservedState.start, 0)

    elems = [elem for elem in acc if isinstance(elem, (tslib.Mpole,))]
    logger.info(
        "\nRadiators added to:\n({})".
        format(", ".join([e.name for e in elems]))
    )

    rad_del_kick = [tslib.RadiationDelegateKick() for elem in elems]
    for a_del, elem in zip(rad_del_kick, elems):
        elem.setRadiationDelegate(a_del)
        a_del.setEnergy(energy)
        # Just use that that the marker knows who is calling him
        # a_del.view(elem, ps_zero, tslib.ObservedState.start, 0)

    return rad_del_kick  # , rad_del


def instrument_with_standard_observers(
    acc: tslib.Accelerator
) -> Sequence[tslib.StandardObserver]:
    """Instrument accelerator with observers

    Returns accelerator list of created observers
    """
    # Instrument it with observers ... I guess I have to keep them here so that
    # internally these are not weak references ... to be tested
    observers = [tslib.StandardObserver() for elem in acc]
    # Register them to the elments
    for ob, elem in zip(observers, acc):
        elem.setObserver(ob)
    return observers


def extract_orbit_from_standard_observers(
    observers: Sequence[tslib.StandardObserver], info: xr.Dataset
) -> xr.Dataset:
    """
    """
    observers = np.array(observers)

    # Required later to index into the final xarray
    has_observer = np.array([True if ob else False for ob in observers])

    # The ones that are not null
    ob_sub = observers[has_observer]

    # check that all have seen data
    with_observer = np.array([ob.hasTruncatedPowerSeries() for ob in ob_sub])
    lf = np.sum(with_observer)
    lo = len(ob_sub)
    if lf != lo:
        raise AssertionError(f"Have {lo} observers but only {lf} has (valid) data")

    indices = np.arange(len(with_observer))
    indices = indices[with_observer == True]
    phase_space_coords_names = ["x", "px", "y", "py", "ct", "delta"]
    tps_tmp = [ob.getTruncatedPowerSeries() for ob in ob_sub]
    tps = xr.DataArray(
        data=tps_tmp,
        name="tps",
        dims=["index", "phase_coordinate"],
        coords=[indices, phase_space_coords_names],
    )
    ps = xr.DataArray(
        data=[t.cst() for t in tps_tmp],
        name="ps",
        dims=["index", "phase_coordinate"],
        coords=[indices, phase_space_coords_names],
    )
    info = info.isel(index=(with_observer == True))
    res = info.merge(dict(tps=tps, ps=ps))
    return res


def extract_orbit_from_accelerator_with_standard_observers(
    acc: Sequence[tslib.ElemType]
) -> xr.Dataset:
    """
    """

    observers = [elem.getObserver() for elem in acc]
    info = accelerator_info(acc)
    return extract_orbit_from_standard_observers(observers, info)


__all__ = [
    "instrument_with_radiators",
    "instrument_with_standard_observers",
    "extract_orbit_from_standard_observers",
    "extract_orbit_from_accelerator_with_standard_observers",
]
