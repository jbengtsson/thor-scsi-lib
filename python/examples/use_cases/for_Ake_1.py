"""Use Case:
     Compute momentum compaction to arbitrary order.
"""


import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, closed_orbit as co, \
    index_class as ind

from thor_scsi.utils.output import mat2txt


ind = ind.index_class()


def plt_alpha_c_vs_delta(file_name, delta, alpha_c, plot):
    fig, gr = plt.subplots()

    fig.suptitle("Momentum Compaction")

    gr.set_title(r"$m_{56}(\delta)$")
    gr.set_xlabel(r"$\delta$ [%]")
    gr.set_ylabel("")
    gr.plot(1e2*delta, alpha_c, "b")
    gr.legend()

    fig.tight_layout()

    plt.savefig(file_name)
    print("\nplt_m_56_vs_delta - plot saved as:", file_name)

    if plot:
        plt.show()


def compute_map(no, desc, delta):
    # Compute Poincaré Taylor M to arbitrary order.
    M = gtpsa.ss_vect_tpsa(desc, no)
    # Compute closed orbit.
    cod = co.compute_closed_orbit(
        lat_prop._lattice, lat_prop._model_state, delta=delta, eps=1e-10,
        desc=desc)
    cod.x0.iloc[ind.ct] = 0e0
    cod.x0.iloc[6] = 0e0
    print("compute_map cod = "+
          "[{:9.3e}, {:9.3e}, {:9.3e}, {:9.3e}, {:9.3e}, {:9.3e}]".format
          (cod.x0.x, cod.x0.px, cod.x0.y, cod.x0.py, cod.x0.delta, cod.x0.ct))
    # Initialise map to identity.
    M.set_identity()
    # Add closed orbit.
    M += cod.x0
    # Track one turn.
    lat_prop._lattice.propagate(lat_prop._model_state, M)
    M.ct.set([0, 0, 0, 0, 0, 0, 0], 0e0, 0e0)
    return M


def compute_alpha_c(M):
    if not True:
        print("\nM[ct]:\n")
        M.ct.print()

    C = lat_prop.compute_circ()
    print(f"\nC [m] = {C:5.3f}")
    index = np.zeros(nv, int)
    alpha_c = np.zeros(no+1)

    print("\nalpha_c:")
    for k in range(1, no+1):
        index[ind.delta] = k
        alpha_c[k] = M.ct.get(index)/C
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


def graph_alpha_c(desc, n_step, delta_max, M):
    C = lat_prop.compute_circ()
    delta_buf = []
    alpha_c_buf = []
    ddelta = delta_max/n_step
    delta = -delta_max
    ps = gtpsa.ss_vect_tpsa(desc, 1)
    ps1 = gtpsa.ss_vect_tpsa(desc, 1)
    ps.set_zero()
    print()
    for k in range(-n_step, n_step+1):
        ps.delta = delta
        ps1.compose(M, ps)
        delta_buf.append(delta)
        alpha_c_buf.append(ps1.cst().ct/C)
        print(" {:10.3e} {:10.3e}".
              format(delta_buf[k+n_step], alpha_c_buf[k+n_step]))
        delta += ddelta
    delta_buf = np.array(delta_buf)
    alpha_c_buf = np.array(alpha_c_buf)
    plt_alpha_c_vs_delta("alpha_c_vs_delta", delta_buf, alpha_c_buf, True)


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 5
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 3.0e9

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

# Read in lattice & instantiate lattice object.
lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

print("\nTotal bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))
print("Circumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))

M = compute_map(no, desc, 0e0)
# Print linear part.
print("\nM:", M)

graph_alpha_c(desc, 10, 5e-2, M)

# assert False

# Print m_56, m_566, m_5666, etc.:
print("\nm_56, m_566, m_5666,...:")
print("  no    m_56...")
index = np.zeros(nv, int)
for k in range(1, no+1):
    index[ind.delta] = k
    print("  {:1d} {:10.3e}".format(k, M.ct.get(index)))

# Print out momentum compaction & fixed points in the longitudinal (aka
# the "ignored") plane.
alpha_c = compute_alpha_c(M)

