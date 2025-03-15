"""Use Case:
     Parametric scans/evaluations for a unit cell.
"""


import os
import sys
from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from thor_scsi.utils import lattice_properties as lp, index_class as ind

import gtpsa


@dataclass
class gtpsa_prop:
    # GTPSA properties.
    # Number of variables - phase-space coordinates & 1 for parameter
    #dependence.
    nv: ClassVar[int] = 6 + 1
    # Max order.
    no: ClassVar[int] = 1
    # Number of parameters.
    nv_prm: ClassVar[int] = 0
    # Parameters max order.
    no_prm: ClassVar[int] = 0
    # Index.
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]


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
    lat_prop.set_phi_fam("r3_h2", phi_rb)
    for k in range(len(bend_scl)):
        lat_prop.set_phi_fam(bend_list[k], bend_scl[k]*dphi)


# TPSA max order.
gtpsa_prop.no = 2

cod_eps = 1e-15
E_0     = 3.0e9

ind = ind.index_class()

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")

file_name = os.path.join(home_dir, sys.argv[1]+".lat")

lat_prop = lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncomp_per_sol: unstable")
        raise ValueError

    # Adjust RF phase for sign of alpha_c.
    cav_loc = lat_prop._lattice.find("cav", 0)
    if lat_prop._alpha_c >= 0e0:
        cav_loc.set_phase(0.0)
    else:
        print("  alpha_c = {:10.3e} phi_rf = 180 deg".
              format(lat_prop._alpha_c))
        cav_loc.set_phase(180.0)

   # Compute radiation properties.
    if not lat_prop.compute_radiation():
        print("\ncompute_radiation: unstable")
        raise ValueError
except ValueError:
    exit
else:
    lat_prop.prt_lat_param()
    lat_prop.prt_rad()
    lat_prop.prt_M()
    lat_prop.prt_M_rad()
    if not False:
        lat_prop.plt_Twiss(file_name+".png", not False)

phi_uc    = 2.85890
bend_list = ["d2_h2_sl_df0", "d2_h2_sl_df1", "d2_h2_sl_df2", "d2_h2_sl_df3",
             "d2_h2_sl_df4", "d2_h2_sl_df5", "d2_h2_sl_df6"]

bend_scl = compute_scl_fact(lat_prop, bend_list)

loc = lat_prop._lattice.find("d2_h2_sl_df0", 0).index
print(f"  loc = {loc:d}")
phi, eps_x, J_x, J_z, alpha_c = \
    lat_prop.unit_cell_rev_bend(
        loc, 30, -0.8, set_phi_rb, phi_uc, bend_list, bend_scl)

if not False:
    lat_prop.plt_scan_phi_rb(
        file_name+"_phi_rb.png", phi, eps_x, J_x, J_z, alpha_c, True)
