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


def prt_bend_lat(lat_prop, file_name, bend_list):

    outf = open(file_name, 'w')

    def prt_bend(lat_prop, name):
        L = lat_prop.get_L_elem(name, 0)
        phi = lat_prop.get_phi_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:2s}: Bending, L = {:7.5f}, T = {:7.5f}, K = {:8.5f}, T1 = 0.0"
               ", T2 = 0.0,\n    N = nbend;").format(name, L, phi, b_2),
              file=outf)

    prt_bend(lat_prop, "b0")
    for k in range(len(bend_list)):
        prt_bend(lat_prop, bend_list[k])


def compute_scl_fact(lat_prop, phi_b, bend_list):
    bend_scl = []
    for k in range(len(bend_list)):
        phi = lat_prop.get_phi_elem(bend_list[k], 0)
        bend_scl.append(phi/phi_b)
    return np.array(bend_scl)


def set_phi_rb(lat_prop, phi_b, bend_list, bend_scl, phi_rb):
    dphi = phi_b - 2.0*phi_rb
    lat_prop.set_phi_fam("qf", phi_rb, True)
    for k in range(len(bend_scl)):
        lat_prop.set_phi_fam(bend_list[k], bend_scl[k]*dphi, True)


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
lat_name = "max_4u_1"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi()))
lat_prop.prt_M()
if not stable:
    assert False
lat_prop.prt_lat_param()

# Compute radiation properties.
stable = lat_prop.compute_radiation()
lat_prop.prt_rad()
lat_prop.prt_M_rad()

types = lat_prop.get_types()
lat_prop.prt_Twiss()

if False:
    lat_prop.plt_Twiss(lat_name+"_twiss.png", False)

if not False:
    phi_b     = 3.0
    bend_list = ["b0", "b1", "b2", "b3", "b4", "b5"]

    bend_scl = compute_scl_fact(lat_prop, phi_b, bend_list)

    phi, eps_x, J_x, J_z, alpha_c = \
        lat_prop.unit_cell_rev_bend(
            15, -0.9, set_phi_rb, phi_b, bend_list, bend_scl)

    lat_prop.plt_scan_phi_rb(
        lat_name+"_phi_rb.png", phi, eps_x, J_x, J_z, alpha_c, True)

if False:
    set_phi_rb(lat_prop, -0.37)
    stable = lat_prop.compute_radiation()
    lat_prop.prt_rad()
    lat_prop.prt_M_rad()

    lat_prop.plt_Twiss("twiss.png", not False)
    prt_bend_lat(lat_prop, bend_list)
