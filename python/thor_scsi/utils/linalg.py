"""Eigen for symplectic matrices ?

"""
import numpy as np
from typing import Sequence
import logging
import functools

from .courant_snyder import compute_A_CS
from thor_scsi.utils.output import vec2txt, mat2txt

logger = logging.getLogger("thor_scsi")
