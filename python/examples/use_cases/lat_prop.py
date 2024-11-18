"""Use Case:
     Compute, print, and display the global linear optics properties for a
     lattice.
"""


import os
import sys
from dataclasses import dataclass
from typing import ClassVar

import gtpsa

import numpy as np

from thor_scsi.utils import lattice_properties as lp, nonlin_dyn as nld_cl, \
    index_class as ind


ind = ind.index_class()


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


def compute_optics(lat_prop):
    try:
        # Compute Twiss parameters along lattice.
        if not lat_prop.comp_per_sol():
            print("\ncomp_per_sol: unstable")
            raise ValueError

        # Adjust RF phase for sign of alpha_c.
        print("\ncompute_optics:")
        print("  alpha_c = {:10.3e}".format(lat_prop._alpha_c))
        cav_loc = lat_prop._lattice.find("cav", 0)
        if lat_prop._alpha_c >= 0e0:
            cav_loc.set_phase(0.0)
        else:
            print("  phi_rf  = 180 deg") 
            cav_loc.set_phase(180.0)

        # Compute radiation properties.
        stable, stable_rad = lat_prop.compute_radiation()
        print("\n", stable, stable_rad)
        if not stable:
            print("\ncompute_radiation: unstable")
            raise ValueError
    except ValueError:
        assert False


# TPSA max order.
gtpsa_prop.no = 6

cod_eps = 1e-15
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

lat_prop.prt_lat("lat_prop_lat.txt")

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

b_3_list = ["s1", "s2", "s3", "s4"]
nld = nld_cl.nonlin_dyn_class(lat_prop, A_max, beta_inj, delta_max, b_3_list)
if False:
    nld.zero_mult(lat_prop, 3)
    nld.zero_mult(lat_prop, 4)

compute_optics(lat_prop)

lat_prop.prt_lat_param()
lat_prop.prt_rad()
lat_prop.prt_M()
lat_prop.prt_M_rad()

if not False:
    lat_prop.prt_Twiss("lat_prop_Twiss.txt")
    lat_prop.plt_Twiss("lat_prop_Twiss.png", not False)
    lat_prop.plt_chrom("lat_prop_chrom.png", not False)

if False:
    nld.compute_map(lat_prop, no)
    nld.compute_nl(lat_prop)
    nld.prt_nl(lat_prop)
