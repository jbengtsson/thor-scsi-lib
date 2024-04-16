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

from thor_scsi.utils import lattice_properties as lp, index_class as ind


def compute_scl_fact(lat_prop, bend_list):
    phi_b = 0e0
    for k in range(len(bend_list)):
        phi_b += lat_prop.get_phi_elem(bend_list[k], 0)
    print("\nphi_b = {:7.5f}".format(phi_b))
    bend_scl = []
    for k in range(len(bend_list)):
        phi = lat_prop.get_phi_elem(bend_list[k], 0)
        bend_scl.append(phi/phi_b)
    return np.array(bend_scl)


def set_phi_rb(lat_prop, phi_uc, bend_list, bend_scl, phi_rb):
    dphi = phi_uc/2e0 - phi_rb
    lat_prop.set_phi_fam("qf", phi_rb)
    for k in range(len(bend_scl)):
        lat_prop.set_phi_fam(bend_list[k], bend_scl[k]*dphi)


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

ind = ind.index_class()

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_4U")
lat_name = "max_4u_match"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.get_types()

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))
lat_prop.prt_M()
if not stable:
    print("\ncomp_per_sol - unstable")
    assert False
lat_prop.prt_lat_param()

lat_prop.prt_Twiss("twiss.txt")

if not False:
    lat_prop.plt_Twiss(lat_name+"_twiss.png", not False)

if not False:
    # Compute radiation properties.
    stable = lat_prop.compute_radiation()
    if not stable:
        print("\ncompute_radiation - unstable")
        assert False
    lat_prop.prt_rad()
    lat_prop.prt_M_rad()

if not False:
    phi_uc    = 3.0
    bend_list = ["b0", "b1", "b2", "b3", "b4", "b5"]

    bend_scl = compute_scl_fact(lat_prop, bend_list)

    phi, eps_x, J_x, J_z, alpha_c = \
        lat_prop.unit_cell_rev_bend(
            15, -0.9, set_phi_rb, phi_uc, bend_list, bend_scl)

    lat_prop.plt_scan_phi_rb(
        lat_name+"_phi_rb.png", phi, eps_x, J_x, J_z, alpha_c, True)
