from thor_scsi.factory import accelerator_from_config
import thor_scsi.lib as tslib
from thor_scsi.utils import engineering

import numpy as np
from numpy.random import Generator, PCG64, default_rng
from typing import Sequence
import os.path

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)
mul_prop = engineering.Property(
    element_info=None, set_method_name="setMultipoles", get_method_name="getMultipoles"
)


quad_variance = np.array([0] + [5e-4 + 5e-4j] + [0] * 3 + [1e-4] + [0] * 3 + [1e-4])
quad_sys_mul = tslib.TwoDimensionalMultipoles(quad_variance)
print(repr(quad_sys_mul))


quadrupoles = [elem for elem in acc if isinstance(elem, tslib.Quadrupole)]
print(quadrupoles)

def create_systematic_multipoles_commands(
    element_indices: Sequence[int], multipole: tslib.TwoDimensionalMultipoles
):
    """
    """
