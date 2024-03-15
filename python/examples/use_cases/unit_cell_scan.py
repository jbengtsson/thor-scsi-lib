"""Use Case:
     Parametric scans of the unit cell.
"""


import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")


import os

import gtpsa

from thor_scsi.utils import periodic_structure as ps

from thor_scsi.utils.output import vec2txt


# from thor_scsi.utils.phase_space_vector import map2numpy

# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 2
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

E_0     = 2.5e9
cod_eps = 1e-15

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "BESSY-III",
    "ipac_2024")
file_name = os.path.join(home_dir, "6BA_0401_UC_with_RB_SVBR_tracy.lat")

lin_opt = ps.lin_opt_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

# Compute Twiss parameters along lattice.
lin_opt.comp_per_sol()
lin_opt.prt_M()
lin_opt.prt_Twiss_param()

# Compute radiation properties.
lin_opt.comp_rad()
lin_opt.prt_rad()

lin_opt.prt_M()
print("\nClosed Orbit:\n"+vec2txt(lin_opt._cod))

types = lin_opt.get_types()

if not False:
   lin_opt.prt_Twiss(types)

if not False:
    lin_opt.plt_Twiss(types, False)

