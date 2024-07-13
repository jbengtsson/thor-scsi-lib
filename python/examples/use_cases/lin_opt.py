"""Use Case:
     Compute, print, and display the global linear optics properties for a
     lattice.
"""


import os

from thor_scsi.utils import lattice_properties as lp


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
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_iv")
lat_name = "ake_2"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.prt_lat(lat_name+"_lat.txt")

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncomp_per_sol: unstable")
        raise ValueError

    # Compute radiation properties.
    if not lat_prop.compute_radiation():
        print("\ncompute_radiation: unstable")
        raise ValueError
except ValueError:
    exit

lat_prop.prt_lat_param()
lat_prop.prt_rad()
lat_prop.prt_M()
lat_prop.prt_M_rad()

lat_prop.plt_Twiss(lat_name+"_Twiss.png", not False)
lat_prop.plt_chrom(lat_name+"_chrom.png", not False)
