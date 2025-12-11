"""Eigen for symplectic matrices ?

"""
import numpy as np
from typing import Sequence
import logging
import functools
from dataclasses import dataclass, field

from .courant_snyder import compute_A_CS
from .phase_space_vector import omega_block_matrix
from thor_scsi.utils.output import vec2txt, mat2txt

logger = logging.getLogger("thor_scsi")


@dataclass
class PlaneDistance:
    """Helper for sorting values by closest distance to a reference array."""

    distance: float = field(init=False, repr=False)
    value: float = field(compare=True)
    reference_values: np.ndarray = field(compare=False)
    payload: object = field(default=None, compare=False)

    def __post_init__(self):
        self.distance = float(np.min(np.abs(self.reference_values - self.value)))

    def __lt__(self, other):
        if not isinstance(other, PlaneDistance):
            return NotImplemented
        if self.distance != other.distance:
            return self.distance < other.distance
        return self.value < other.value


# Re-export for convenience in tests
__all__ = ["PlaneDistance", "omega_block_matrix"]
