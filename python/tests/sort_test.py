import pytest

from thor_scsi.utils.linalg import PlaneDistance
import numpy as np
import logging

logger = logging.getLogger("thor_scsi")


def test00_obvious():
    ref_vals = np.array([1, 2, 3])
    vals = [3, 2, 1]
    targets = [2, 1, 0]
    sorters = [
        PlaneDistance(value=v, reference_values=ref_vals, payload=p)
        for v, p in zip(vals, targets)
    ]
    sorters.sort()

    assert sorters[0].payload == 0
    assert sorters[1].payload == 1
    assert sorters[2].payload == 2


def test01_obvious():
    ref_vals = np.array([3, 1, 2])
    vals = [3, 2, 1]
    targets = [0, 2, 1]
    sorters = [
        PlaneDistance(value=v, reference_values=ref_vals, payload=p)
        for v, p in zip(vals, targets)
    ]
    sorters.sort()

    assert sorters[0].payload == 0
    assert sorters[1].payload == 1
    assert sorters[2].payload == 2


def test10_numeric_close():
    ref_vals = np.array([1, 2, 3])
    vals = [2.8, 1.8, 1.5]
    targets = [2, 1, 0]
    sorters = [
        PlaneDistance(value=v, reference_values=ref_vals, payload=p)
        for v, p in zip(vals, targets)
    ]
    sorters.sort()

    assert sorters[0].payload == 0
    assert sorters[1].payload == 1
    assert sorters[2].payload == 2


def test11_numeric_close():
    ref_vals = np.array([1, 2, 3])
    vals = [2.8, 1.8, 1.5]
    targets = [2, 1, 0]
    sorters = [
        PlaneDistance(value=v, reference_values=ref_vals, payload=p)
        for v, p in zip(vals, targets)
    ]
    sorters.sort()

    assert sorters[0].payload == 0
    assert sorters[1].payload == 1
    assert sorters[2].payload == 2


def test20_bit_realistic():
    ref_vals = np.array([0.25, 0.45])
    vals = [0.4, 0.55]
    targets = [0, 1]
    sorters = [
        PlaneDistance(value=v, reference_values=ref_vals, payload=p)
        for v, p in zip(vals, targets)
    ]
    sorters.sort()

    assert sorters[0].payload == 0
    assert sorters[1].payload == 1


def test21_bit_unrealistic():
    """Or is that unrealistic?
    """
    ref_vals = np.array([0.25, 0.45])
    vals = [0.255, 0.24]
    # While one could find an other optimum, in this case
    # the first is closer to the point than the second
    targets = [0, 1]
    sorters = [
        PlaneDistance(value=v, reference_values=ref_vals, payload=p)
        for v, p in zip(vals, targets)
    ]
    sorters.sort()

    logger.info(f"Sorters {sorters}")
    assert sorters[0].payload == 0
    assert sorters[1].payload == 1
