"""Use Case:
     Minimise the driving terms.
     See SLS Tech Note 09/97

     https://ados.web.psi.ch/slsnotes/sls0997.pdf
"""

import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="ERROR")
logger = logging.getLogger("thor_scsi")


from scipy import optimize

import copy as _copy
from dataclasses import dataclass
import math
import os
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    index_class as ind
from thor_scsi.utils.output import mat2txt

# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt


ind = ind.index_class()


class MultipoleIndex(enum.IntEnum):
    quadrupole = 2
    sextupole  = 3


def prt_map(str, map):
    n_dof = 3
    print(str)
    for k in range(2*n_dof):
        map.iloc[k].print()


def prt_twiss_sxt(lat, data, fam_name):
    print("\n  name        s    beta_x beta_y")
    for k in range(len(lat)):
        if (lat[k].name == fam_name):
            print("  {:8s} {:6.3f} {:6.3f} {:6.3f}".
                  format(lat[k].name, data.s.values[k],
                         data.twiss.values[k, 0, 1],
                         data.twiss.values[k, 1, 1]))

    
# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 3
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

print("Circumference [m] = {:7.5f}".format(lat_prop.compute_circ()))

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncomp_per_sol - unstable")
        raise ValueError

    # Compute radiation properties.
    if not lat_prop.compute_radiation():
        print("\ncompute_radiation - unstable")
        raise ValueError
except ValueError:
    exit
else:
    lat_prop._model_state.radiation = False
    lat_prop._model_state.emittance = False
    lat_prop._model_state.Cavity_on = False

lat_prop.prt_lat_param()
lat_prop.prt_rad()
lat_prop.prt_M()
lat_prop.prt_M_rad()
lat_prop.prt_Twiss(lat_name+"_Twiss.txt")

if False:
    lat_prop.plt_Twiss(lat_name+"_Twiss.png", not False)

M = lo.compute_map(
    lat_prop._lattice, lat_prop._model_state, desc=desc, tpsa_order=no)

if not False:
    print("\nM:\n"+mat2txt(M.jacobian()))

h = ts.M_to_h_DF(M)
h_re = gtpsa.tpsa(desc, no)
h_im = gtpsa.tpsa(desc, no)
ts.CtoR(h, h_re, h_im)
h_re.print("h_re", 1e-30);
h_im.print("h_im", 1e-30);
# (h-ts.RtoC(h_re, h_im)).print()

if False:
    prt_twiss_sxt(lat, data, "om_sf")
    prt_twiss_sxt(lat, data, "om_sd")
