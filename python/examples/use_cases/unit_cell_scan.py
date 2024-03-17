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

from thor_scsi.utils import periodic_structure as ps
from thor_scsi.utils import radiate as rp

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


def plt_eps_x_alpha_c(file_name, phi, eps_x, alpha_c, plot):
    fig, gr_1 = plt.subplots(1)

    # fig.suptitle("Linear Optics")

    gr_1.set_title(r"$\epsilon_x$ & $\alpha_c$ vs. Reverse Bend Angle")
    gr_1.set_xlabel(r"$\phi_{RB}$ [$\degree$]")
    gr_1.set_ylabel(r"$\epsilon_x$ [pm$\cdot$rad]")
    gr_1.plot(
        phi, 1e12*eps_x, color="b", label=r"$\epsilon_x$")
    gr_1.legend(loc="upper left")

    gr_1_r = gr_1.twinx()
    gr_1_r.set_ylabel(r"$\alpha_c$ [$10^{-4}$]")
    gr_1_r.plot(
        phi, 1e4*alpha_c, color="g", label=r"$\alpha_c$")
    gr_1_r.legend(loc="upper right")

    fig.tight_layout()

    plt.savefig(file_name)
    print("\n - plot saved as:", file_name)

    if plot:
        plt.show()


def get_phi_elem(lat, fam_name, n_kid):
    elem = lat.find(fam_name, n_kid)
    return elem.get_length() * elem.get_curvature() * 180e0 / np.pi


def set_phi_elem(lat, fam_name, kid_num, phi, rect_bend):
    b = lat.find(fam_name, kid_num)
    L = b.get_length()
    h = phi * np.pi / (L * 180e0)
    b.set_curvature(h)
    if rect_bend:
        b.set_entrance_angle(phi/2e0)
        b.set_exit_angle(phi/2e0)


def set_phi_fam(lat, fam_name, phi, rect_bend):
    prt = False
    b = lat.find(fam_name, 0)
    L = b.get_length()
    h = phi * np.pi / (L * 180e0)
    if prt:
        print("set_phi_fam: {:8s} {:5.3f} {:6.3f}".format(b.name, L, phi))
    for b in lat.elements_with_name(fam_name):
        b.set_curvature(h)


def compute_phi(lat):
    """Compute the total bend angle.
    """
    prt = False
    phi = 0e0
    for k in range(len(lat)):
        if type(lat[k]) == ts.Bending:
            dphi = get_phi_elem(lat, lat[k].name, 0)
            phi += dphi
            if prt:
                print("{:8s} {:5.3f} {:6.3f}".
                      format(lat[k].name, lat[k].get_length(), dphi))
    return phi


def set_RF_cav_phase(lat, name, phi):
    # Set RF cavity phase for negative alpha_c.
    cav = lat.find(name, 0)
    cav.set_phase(phi)


def unit_cell_rev_bend(lin_opt, rad_prop, n_step, phi_0, phi_1):
    # Dipoles:
    #   [b1(b1a), rb(rba), mb1(mb1a), mqd(mwba)]
    # Dipole bend angles:
    #   rba  = -0.28;
    #   b1a  = 4.25-2.0*rba;
    #   mwba = 0.2;
    #   mb1a = 2.75-mwba;

    phi_rb = phi_0
    dphi = (phi_1-phi_0)/n_step
    phi_rb_buf = []
    eps_x_buf = []
    alpha_c_buf = []
    print("\n   phi   phi_tot  eps_x     alpha_c   J_x   J_z     eta_x"
          "     nu_x     nu_y"
          "\n  [deg]   [deg]  [nm.rad]                            [m]")
    for k in range(n_step):
        set_phi_fam(lin_opt._lattice, "rb", phi_rb, True)
        set_phi_fam(lin_opt._lattice, "b1", 4.25-2.0*phi_rb, True)
        stable = lin_opt.comp_per_sol()
        eta_x = lin_opt._Twiss[0][x_]
        lin_opt.compute_alpha_c()
        if lin_opt._alpha_c > 0e0:
            set_RF_cav_phase(lin_opt._lattice, "cav", 0.0)
        else:
            set_RF_cav_phase(lin_opt._lattice, "cav", 180.0)
        stable = rad_prop.compute_radiation(lin_opt)
        if stable:
            phi_rb_buf.append(abs(phi_rb))
            eps_x_buf.append(rad_prop._eps[X_])
            alpha_c_buf.append(lin_opt._alpha_c)
            print("{:7.3f}  {:5.3f}    {:5.1f}   {:10.3e}  {:4.2f} {:5.2f}"
                  " {:10.3e}  {:7.5f}  {:7.5f}".
                  format(
                      phi_rb, compute_phi(lin_opt._lattice),
                      1e12*rad_prop._eps[X_], lin_opt._alpha_c, rad_prop._J[X_],
                      rad_prop._J[Z_], eta_x, rad_prop._nu[X_],
                      rad_prop._nu[Y_]))
        else:
            print("  unstable")
        # lin_opt.prt_Twiss_param()
        phi_rb += dphi
    return np.array(phi_rb_buf), np.array(eps_x_buf), np.array(alpha_c_buf)


def prt_default_mapping():
    index_map = gtpsa.default_mapping()
    print("\nIndex Mapping:\n"
          "  x     = {:1d}\n  p_x   = {:1d}\n  y     = {:1d}\n  p_y   = {:1d}\n"
          "  delta = {:1d}\n  ct    = {:1d}".
          format(index_map.x, index_map.px, index_map.y, index_map.py,
                 index_map.delta, index_map.ct))
    

# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 2
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

E_0     = 3.0e9
cod_eps = 1e-15

if True:
    home_dir = os.path.join(
        os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "BESSY-III",
        "ipac_2024")
    file_name = os.path.join(home_dir, "2024Mar.lat")
else:
    home_dir = os.path.join(
        os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_4U")
    file_name = os.path.join(home_dir, "max_4u.lat")

prt_default_mapping()

lin_opt = ps.lin_opt_class(nv, no, nv_prm, no_prm, file_name, E_0)
rad_prop = rp.rad_prop_class(lin_opt, cod_eps)

# Compute Twiss parameters along lattice.
lin_opt.comp_per_sol()
print("\nCircumference [m] = {:7.5f}".format(lin_opt.compute_circ()))
lin_opt.prt_M()
lin_opt.prt_Twiss_param()

# Compute radiation properties.
stable = rad_prop.compute_radiation(lin_opt)
rad_prop.prt_rad(lin_opt)
rad_prop.prt_M()

types = lin_opt.get_types()

lin_opt.prt_Twiss(types)

if False:
    lin_opt.plt_Twiss(types, False)

if not False:
    # BESSY III unit cell.
    phi, eps_x, alpha_c = \
        unit_cell_rev_bend(lin_opt, rad_prop, 15, -0.075, -0.97)
    plt_eps_x_alpha_c("plt_eps_x_alpha_c.png", phi, eps_x, alpha_c, True)
