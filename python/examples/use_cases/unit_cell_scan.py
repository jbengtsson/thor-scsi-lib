"""Use Case:
     Parametric scans/evaluations for a unit cell.
"""

import os
import sys
from dataclasses import dataclass
from typing import ClassVar

import gtpsa

from thor_scsi.utils import lattice_properties as lp


def set_phi_rb_bessy_iii(lat_prop, phi_rb):
    # BESSY III.
    # Dipoles:
    #   [b1(b1a), rb(rba), mb1(mb1a), mqd(mwba)]
    # Dipole bend angles:
    #   rba  = -0.28;
    #   b1a  = 4.25-2.0*rba;
    #   mwba = 0.2;
    #   mb1a = 2.75-mwba;
    lat_prop.set_phi_fam("rb", phi_rb, True)
    lat_prop.set_phi_fam("b1", 4.25-2.0*phi_rb, True)


cod_eps = 1e-15
E_0     = 2.5e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB")

lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = lp.lattice_properties_class(file_name, E_0, cod_eps, 2)

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
lat_prop.prt_M()
if not stable:
    assert False
# Twiss = lat_prop.get_Twiss(len(lat_prop._lattice)-1)
# lat_prop.prt_Twiss_param(Twiss)

# Compute radiation properties.
stable = lat_prop.compute_radiation()
lat_prop.prt_rad()
lat_prop.prt_M()

lat_prop.prt_Twiss("lat_prop_Twiss.txt")

if False:
    lat_prop.plt_Twiss("twiss.png", not False)

if not False:
    phi, eps_x, J_x, J_z, alpha_c = \
        lat_prop.unit_cell_rev_bend(15, -0.97, set_phi_rb_bessy_iii)

    lat_prop.plt_scan_phi_rb(
        "plt_scan_phi_rb.png", phi, eps_x, J_x, J_z, alpha_c, True)
