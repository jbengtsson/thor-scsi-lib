"""Use Case:
     Courant & Snyder matrix diagonalisation.
"""

import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="ERROR")
logger = logging.getLogger("thor_scsi")


import math
import os

import numpy as np

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
    radiate as rad, closed_orbit as co


def read_lattice(file_name):
    # Read in & parse lattice file.
    lat = accelerator_from_config(file_name)

    # Set lattice state (Rf cavity on/off, etc.)
    model_state = ts.ConfigType()

    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    return n_dof, lat, model_state


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 1
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

home_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB",
                     "BESSY-III", "ipac_2023")
file_name = os.path.join(home_dir, "b3_cf425cf_thor_scsi.lat")

print("\nlattice file:  \n", file_name)

n_dof, lat, model_state = read_lattice(file_name)

M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=no)

print("\nM:", M)
