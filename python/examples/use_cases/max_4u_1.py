"""Use Case:
     Parametric scans/evaluations for a unit cell.
"""


import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

from dataclasses import dataclass
import os

import numpy as np

from thor_scsi.utils import lattice_properties as lp


def set_phi_rb(lat_prop, phi_rb):
    # Optimum reverse bend angle is:
    #   phi_rb = -0.37
    phi_b  = 3.0
    b_name = ["b0", "b1", "b2", "b3", "b4", "b5"]
    b_scl = np.array(
        [1.094181/phi_b, 0.151199/phi_b, 0.151101/phi_b, 0.101861/phi_b,
         0.001569/phi_b, 0.000089/phi_b]
    )

    dphi = 3.0 - 2.0*phi_rb

    lat_prop.set_phi_fam("qf", phi_rb, True)
    for k in range(len(b_scl)):
        lat_prop.set_phi_fam(b_name[k], b_scl[k]*dphi, True)


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 2
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 3.0e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_4U")
file_name = os.path.join(home_dir, "max_4u.lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".
      format(lat_prop.compute_phi(lat_prop._lattice)))
lat_prop.prt_M()
if not stable:
    assert False
Twiss = lat_prop.get_Twiss(len(lat_prop._lattice)-1)
lat_prop.prt_Twiss_param(Twiss)

# Compute radiation properties.
stable = lat_prop.compute_radiation()
lat_prop.prt_rad()
lat_prop.prt_M_rad()

types = lat_prop.get_types()
lat_prop.prt_Twiss(types)

if False:
    lat_prop.plt_Twiss(types, False)

if not False:
    phi, eps_x, J_x, J_z, alpha_c = \
        lat_prop.unit_cell_rev_bend(15, -0.97, set_phi_rb)

    lat_prop.plt_scan_phi_rb(
        "plt_scan_phi_rb.png", phi, eps_x, J_x, J_z, alpha_c, True)
