"""Use Case:
     Parametric scans/evaluations for a unit cell.
"""


import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")


import os
from dataclasses import dataclass
from typing import ClassVar

from thor_scsi.utils import lattice_properties as lp

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


# TPSA max order.
gtpsa_prop.no = 2

cod_eps = 1e-15
E_0     = 2.5e9

if True:
    home_dir = os.path.join(
        os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "BESSY-III",
        "ipac_2024")
    file_name = os.path.join(home_dir, "2024Mar.lat")
else:
    home_dir = os.path.join(
        os.environ["HOME"], "accSoft", "thor_scsi-lib", "python", "examples",
        "lattices")
    file_name = os.path.join(home_dir, "2024Mar.lat")

lat_prop = \
    lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

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
