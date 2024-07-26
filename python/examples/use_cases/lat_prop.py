"""Use Case:
     Compute, print, and display the global linear optics properties for a
     lattice.
"""


import os
import sys

import numpy as np

from thor_scsi.utils import lattice_properties as lp, nonlin_dyn as nld_cl, \
    index_class as ind


ind = ind.index_class()


def compute_optics(lat_prop):
    try:
        # Compute Twiss parameters along lattice.
        if not lat_prop.comp_per_sol():
            print("\ncomp_per_sol: unstable")
            raise ValueError

        # Compute radiation properties.
        stable, stable_rad = lat_prop.compute_radiation()
        print(stable, stable_rad)
        if not stable:
            print("\ncompute_radiation: unstable")
            raise ValueError
    except ValueError:
        assert False


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 6
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.prt_lat("lat_prop_lat.txt")

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

compute_optics(lat_prop)

lat_prop.prt_lat_param()
lat_prop.prt_rad()
lat_prop.prt_M()
lat_prop.prt_M_rad()

assert False

if not False:
    lat_prop.plt_Twiss("lat_prop_Twiss.png", not False)
    lat_prop.plt_chrom("lat_prop_chrom.png", not False)

b_3_list = ["s2", "s3"]

nld = nld_cl.nonlin_dyn_class(lat_prop, A_max, beta_inj, delta_max, b_3_list)

if True:
    nld.set_xi(lat_prop, 0e0, 0e0)
else:
    nld.zero_mult(lat_prop, 3)
    nld.zero_mult(lat_prop, 4)

nld.compute_map(lat_prop, no)
nld.compute_nl(lat_prop)
nld.prt_nl()
