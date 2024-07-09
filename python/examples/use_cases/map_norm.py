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


def prt_map(map, str, eps):
    print(str)
    map.x.print("x", eps)
    map.px.print("p_x", eps)
    map.y.print("y", eps)
    map.py.print("p_y", eps)
    map.delta.print("delta", eps)
    map.ct.print("ct", eps)


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

A_0 = gtpsa.ss_vect_tpsa(desc, no)
A_1 = gtpsa.ss_vect_tpsa(desc, no)
R   = gtpsa.ss_vect_tpsa(desc, no)
g   = gtpsa.tpsa(desc, no)

K = ts.Map_Norm(M, A_0, A_1, R, g)

print("\nA_0:\n"+mat2txt(A_0.jacobian()))
print("\nA_1:\n"+mat2txt(A_1.jacobian()))
print("\nR:\n"+mat2txt(R.jacobian()))
g.print("g", 1e-10)
K.print("K", 1e-10)

Id       = gtpsa.ss_vect_tpsa(desc, no)
A        = gtpsa.ss_vect_tpsa(desc, no)
A_inv    = gtpsa.ss_vect_tpsa(desc, no)
t_map    = gtpsa.ss_vect_tpsa(desc, no)
M_Fl     = gtpsa.ss_vect_tpsa(desc, no)
A_nl     = gtpsa.ss_vect_tpsa(desc, no)
A_nl_inv = gtpsa.ss_vect_tpsa(desc, no)

Id.set_identity();

# Prune Poincaré map to no-1.
ts.get_mns(M, 1, no-1, M)

# M_Fl = (A_0 . A_1)^-1 . M . A_0 . A_1
A.compose(A_0, A_1)
A_inv.inv(A)
t_map.compose(M, A)
M_Fl.compose(A_inv, t_map)
print("\nM_Fl:\n"+mat2txt(M_Fl.jacobian()))

if False:
    prt_map(M_Fl, "M_Fl", 1e-10)

# A_nl = exp(:g:)
ts.h_DF_to_M(-1.0*g, Id, 3, no, A_nl)

# R = exp(-g) . M_Fl . exp(:g:)
A_nl_inv.inv(A_nl)
t_map.compose(M_Fl, A_nl)
M_Fl.compose(A_nl_inv, t_map)
prt_map(M_Fl, "M_Fl", 1e-10)
