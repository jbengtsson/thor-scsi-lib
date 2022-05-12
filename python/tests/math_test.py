"""Testing math utilities not found any where
"""

from thor_scsi.utils.math import minimum_distance_above_threshold
import numpy as np
import logging

logger = logging.getLogger("thor_scsi")


def test10_distance():
    """Equidistant range sufficient
    """
    x = np.arange(5)
    flag, e1, e2 = minimum_distance_above_threshold(x, 0.5)
    assert flag


def test20_distance_close():
    """Does it find the correct values?
    """
    x = np.array([1, 2, 3, 4, 5, 3.2])
    flag, e1, e2 = minimum_distance_above_threshold(x, 0.5)
    if not flag:
        logger.info(
            f"Found distance between indices [{e1}, {e2}] corresponding to values {x[[e1,e2]]}"
        )
    assert not flag
