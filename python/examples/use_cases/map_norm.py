"""Use Case:
     Test map_norm.
"""

import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="ERROR")
logger = logging.getLogger("thor_scsi")

import os

import numpy as np

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo
from thor_scsi.utils.output import mat2txt


# Number of phase-space coordinates.
nv = 7
# Max order for Poincaré map.
no = 4
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no)

cod_eps = 1e-15
E_0     = 3.0e9


home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_4u")
lat_name = "max_4u_uc"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

print("\nCircumference [m] = {:7.5f}".format(lat_prop.compute_circ()))

M = lo.compute_map(
    lat_prop._lattice, lat_prop._model_state, desc=desc, tpsa_order=no)

print("M:\n"+mat2txt(M.jacobian()))

if False:
    h = ts.M_to_h_DF(M)
    h_re = gtpsa.tpsa(desc, no)
    h_im = gtpsa.tpsa(desc, no)
    ts.CtoR(h, h_re, h_im)
    h_re.print("h_re", 1e-30);
    h_im.print("h_im", 1e-30);

ts.Map_Norm(M)
