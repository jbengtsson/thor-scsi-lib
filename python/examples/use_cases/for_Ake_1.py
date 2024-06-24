"""Use Case:
     Compute momentum compaction to arbitrary order.
"""


import os

import numpy as np

import gtpsa

from thor_scsi.utils import lattice_properties as lp, index_class as ind

from thor_scsi.utils.output import mat2txt


ind = ind.index_class()


def compute_alpha_c(map):
    if not True:
        print("\nmap[ct]:\n")
        map.ct.print()

    C = lat_prop.compute_circ()
    print(f"\nC [m] = {C:5.3f}")
    index = np.zeros(nv, int)
    alpha_c = np.zeros(no+1)

    print("\nalpha_c:")
    for k in range(1, no+1):
        index[ind.delta] = k
        alpha_c[k] = map.ct.get(index)/C
        print("  {:1d} {:10.3e}".format(k, alpha_c[k]))

    n_alpha_c = len(alpha_c)
    if n_alpha_c >= 3:
        print("\nFixed points to O(3) [%]: {:5.2f}, {:5.2f}".
              format(0e0, -1e2*alpha_c[1]/alpha_c[2]))

    if n_alpha_c >= 4:
        po2 = alpha_c[2]/(2e0*alpha_c[3])
        q = alpha_c[1]/alpha_c[3]
        pm = np.sqrt(po2**2-q)
        print("Fixed points to O(4) [%]: {:5.2f}, {:5.2f}, {:5.2f}".
              format(0e0, -1e2*(po2+pm), -1e2*(po2-pm)))

    return alpha_c


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 3
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 3.0e9

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_4U")
lat_name = "max_iv_sp_matched_2"
file_name = os.path.join(home_dir, lat_name+".lat")

# Read in lattice & instantiate lattice object.
lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

print("\nTotal bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))
print("Circumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))

# Compute Poincaré Taylor map to arbitrary order.
map = gtpsa.ss_vect_tpsa(desc, no)
# Initialise map to identity.
map.set_identity()
# Track one turn.
lat_prop._lattice.propagate(lat_prop._model_state, map)

# Print linear part.
print("\nmap:\n" + mat2txt(map.jacobian()[:6, :6]))

# Print m_56, m_566, m_5666, etc.:
print("\nm_56, m_566, m_5666,...:")
print("  no    m_56...")
index = np.zeros(nv, int)
for k in range(1, no+1):
    index[ind.delta] = k
    print("  {:1d} {:10.3e}".format(k, map.ct.get(index)))

# Print out momentum compaction & fixed points in the longitudinal (aka
# the "ignored") plane.
alpha_c = compute_alpha_c(map)
