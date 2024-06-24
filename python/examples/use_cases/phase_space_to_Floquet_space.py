"""Use Case:
     Courant & Snyder diagonalisation of an arbitrary matrix M from phase-space
     to Floquet space:

       M = A R A^-1
"""

import os

import gtpsa
import thor_scsi.lib as ts
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils import linear_optics as lo
from thor_scsi.utils.output import mat2txt


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

print("\nlattice file:\n ", file_name)

n_dof, lat, model_state = read_lattice(file_name)

M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=no)

M_mat = M.jacobian()[:6, :6]

stable, A_mat, A_inv_mat, alpha_rad = lo.compute_M_diag(n_dof, M_mat)

print("\nM:\n", mat2txt(M_mat))
print("\nA:\n", mat2txt(A_mat))
print("\nR:\n", mat2txt(A_inv_mat @ M_mat @ A_mat))

