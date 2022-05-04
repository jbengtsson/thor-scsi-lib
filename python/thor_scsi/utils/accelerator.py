"""Utils for calculating with the accelerator
"""
import thor_scsi.lib as tslib
from typing import Sequence


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
