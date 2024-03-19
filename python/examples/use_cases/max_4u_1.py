"""Use Case:
     Parametric scans/evaluations for a unit cell.
"""


import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")


import os
import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.utils import periodic_structure as ps, radiate as rp, \
    get_set_mpole as gs

from thor_scsi.utils.output import vec2txt


# Configuration space coordinates.
X_, Y_, Z_ = [
    ts.spatial_index.X,
    ts.spatial_index.Y,
    ts.spatial_index.Z
]
# Phase-space coordinates.
[x_, px_, y_, py_, ct_, delta_] = [
    ts.phase_space_index_internal.x,
    ts.phase_space_index_internal.px,
    ts.phase_space_index_internal.y,
    ts.phase_space_index_internal.py,
    ts.phase_space_index_internal.ct,
    ts.phase_space_index_internal.delta,
]


def prt_default_mapping():
    index_map = gtpsa.default_mapping()
    print("\nIndex Mapping:\n"
          "  x     = {:1d}\n  p_x   = {:1d}\n  y     = {:1d}\n  p_y   = {:1d}\n"
          "  delta = {:1d}\n  ct    = {:1d}".
          format(index_map.x, index_map.px, index_map.y, index_map.py,
                 index_map.delta, index_map.ct))
    

def plt_scan_phi_rb(file_name, phi, eps_x, J_x, J_z, alpha_c, plot):
    fig, (gr_1, gr_2) = plt.subplots(2)

    fig.suptitle("Lattice Trade-Offs vs. Reverse Bend Angle")

    gr_1.set_title(
        r"$[\epsilon_x\left(\phi_{RB}\right), \alpha_c\left(\phi_{RB}\right)]$")
    gr_1.set_xlabel(r"$\phi_{RB}$ [$\degree$]")
    gr_1.set_ylabel(r"$\epsilon_x$ [pm$\cdot$rad]")
    gr_1.yaxis.label.set_color("b")
    gr_1.tick_params(axis="y", colors="b")
    gr_1.plot(phi, 1e12*eps_x, color="b")

    gr_1_r = gr_1.twinx()
    gr_1_r.set_ylabel(r"$\alpha_c$ [$10^{-4}$]")
    gr_1_r.yaxis.label.set_color("g")
    gr_1_r.tick_params(axis="y", colors="g")
    gr_1_r.plot(phi, 1e4*alpha_c, color="g", label=r"$\alpha_c$")

    gr_2.set_title(r"$[J_x\left(\phi_{RB}\right), J_z\left(\phi_{RB}\right)]$")
    gr_2.set_xlabel(r"$\phi_{RB}$ [$\degree$]")
    gr_2.set_ylabel(r"$J_x$")
    gr_2.yaxis.label.set_color("b")
    gr_2.tick_params(axis="y", colors="b")
    gr_2.plot(phi, J_x, color="b")

    gr_2_r = gr_2.twinx()
    gr_2_r.set_ylabel(r"$J_z$")
    gr_2_r.yaxis.label.set_color("g")
    gr_2_r.tick_params(axis="y", colors="g")
    gr_2_r.plot(phi, J_z, color="g")

    fig.tight_layout()

    plt.savefig(file_name)
    print("\n - plot saved as:", file_name)

    if plot:
        plt.show()


def unit_cell_rev_bend(lin_opt, rad_prop, get_set, n_step, phi_min, set_phi):
    phi_rb = 0e0
    phi_step = phi_min/n_step
    phi_rb_buf = []
    eps_x_buf = []
    alpha_c_buf = []
    J_x_buf = []
    J_z_buf = []
    print("\n   phi   phi_tot  eps_x    J_x   J_z    alpha_c     eta_x"
          "     nu_x     nu_y"
          "\n  [deg]   [deg]  [nm.rad]                            [m]")
    for k in range(n_step):
        set_phi(get_set, phi_rb)
        stable = lin_opt.comp_per_sol()
        eta_x = lin_opt._Twiss[0][x_]
        lin_opt.compute_alpha_c()
        if lin_opt._alpha_c > 0e0:
            get_set.set_RF_cav_phase(lin_opt._lattice, "cav", 0.0)
        else:
            get_set.set_RF_cav_phase(lin_opt._lattice, "cav", 180.0)
        stable = rad_prop.compute_radiation(lin_opt)
        if stable:
            phi_rb_buf.append(abs(phi_rb))
            eps_x_buf.append(rad_prop._eps[X_])
            J_x_buf.append(rad_prop._J[X_])
            J_z_buf.append(rad_prop._J[Z_])
            alpha_c_buf.append(lin_opt._alpha_c)
            print("{:7.3f}  {:5.3f}    {:5.1f}    {:4.2f} {:5.2f} {:10.3e}"
                  " {:10.3e}  {:7.5f}  {:7.5f}".
                  format(
                      phi_rb, get_set.compute_phi(lin_opt._lattice),
                      1e12*rad_prop._eps[X_], rad_prop._J[X_], rad_prop._J[Z_],
                      lin_opt._alpha_c, eta_x, rad_prop._nu[X_],
                      rad_prop._nu[Y_]))
        else:
            print("  unstable")
        # lin_opt.prt_Twiss_param()
        phi_rb += phi_step
    return \
        np.array(phi_rb_buf), np.array(eps_x_buf), np.array(J_x_buf), \
        np.array(J_z_buf), np.array(alpha_c_buf)


def set_phi_rb_max_4u(get_set, phi_rb):
    # MAX 4U.
    # Optimum reverse bend angle is:
    #   phi_rb = -0.37
    phi_b  = 3.0
    b0_scl = 1.094181/phi_b
    b1_scl = 0.151199/phi_b
    b2_scl = 0.151101/phi_b
    b3_scl = 0.101861/phi_b
    b4_scl = 0.001569/phi_b
    b5_scl = 0.000089/phi_b

    dphi = 3.0 - 2.0*phi_rb

    get_set.set_phi_fam(lin_opt._lattice, "qf", phi_rb, True)
    get_set.set_phi_fam(lin_opt._lattice, "b0", b0_scl*dphi, True)
    get_set.set_phi_fam(lin_opt._lattice, "b1", b1_scl*dphi, True)
    get_set.set_phi_fam(lin_opt._lattice, "b2", b2_scl*dphi, True)
    get_set.set_phi_fam(lin_opt._lattice, "b3", b3_scl*dphi, True)
    get_set.set_phi_fam(lin_opt._lattice, "b4", b4_scl*dphi, True)
    get_set.set_phi_fam(lin_opt._lattice, "b5", b5_scl*dphi, True)


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
set_phi_rb = set_phi_rb_max_4u

prt_default_mapping()

lin_opt = ps.lin_opt_class(nv, no, nv_prm, no_prm, file_name, E_0)
rad_prop = rp.rad_prop_class(lin_opt, cod_eps)
get_set = gs.get_set_mpole_class()

# Compute Twiss parameters along lattice.
stable = lin_opt.comp_per_sol()
print("\nCircumference [m]      = {:7.5f}".format(lin_opt.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".
      format(get_set.compute_phi(lin_opt._lattice)))
lin_opt.prt_M()
if stable:
    lin_opt.prt_Twiss_param()
else:
    assert False

# Compute radiation properties.
stable = rad_prop.compute_radiation(lin_opt)
rad_prop.prt_rad(lin_opt)
rad_prop.prt_M()

types = lin_opt.get_types()

lin_opt.prt_Twiss(types)

if False:
    lin_opt.plt_Twiss(types, not False)

if not False:
    phi, eps_x, J_x, J_z, alpha_c = \
        unit_cell_rev_bend(lin_opt, rad_prop, get_set, 15, -0.97, set_phi_rb)

    plt_scan_phi_rb(
        "plt_scan_phi_rb.png", phi, eps_x, J_x, J_z, alpha_c, True)
