"""Use Case:
     Test map_norm.
"""

import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="ERROR")
logger = logging.getLogger("thor_scsi")

import os

import copy as _copy

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


def compute_map(lat_prop, no):
    M = lo.compute_map(
        lat_prop._lattice, lat_prop._model_state, desc=lat_prop._desc,
        tpsa_order=no)
    return M


def compute_h(desc, M):
    h = ts.M_to_h_DF(M)
    h_re = gtpsa.tpsa(desc, no)
    h_im = gtpsa.tpsa(desc, no)
    ts.CtoR(h, h_re, h_im)

    h.print("h", 1e-30);


def compute_map_normal_form(desc, M):
    A_0 = gtpsa.ss_vect_tpsa(desc, no)
    A_1 = gtpsa.ss_vect_tpsa(desc, no)
    R   = gtpsa.ss_vect_tpsa(desc, no)
    g   = gtpsa.tpsa(desc, no)

    K = ts.Map_Norm(M, A_0, A_1, R, g)

    g.print("g", 1e-30)
    K.print("K", 1e-30)

    return K, A_0, A_1, g


def compute_M(desc, K, A_0, A_1, g, no):
    Id       = gtpsa.ss_vect_tpsa(desc, no)
    R        = gtpsa.ss_vect_tpsa(desc, no)
    A_0_inv  = gtpsa.ss_vect_tpsa(desc, no)
    A_1_inv  = gtpsa.ss_vect_tpsa(desc, no)
    A_nl     = gtpsa.ss_vect_tpsa(desc, no)
    A_nl_inv = gtpsa.ss_vect_tpsa(desc, no)
    M        = gtpsa.ss_vect_tpsa(desc, no)

    Id.set_identity()

    ts.h_DF_to_M(K, Id, 2, no, R)
    A_0_inv.inv(A_0)
    A_1_inv.inv(A_1)
    ts.h_DF_to_M(g, Id, 3, no, A_nl)
    A_nl_inv.inv(A_nl)

    M.compose(A_0, A_1)
    M.compose(M, A_nl)
    M.compose(M, R)
    M.compose(M, A_nl_inv)
    M.compose(M, A_1_inv)
    M.compose(M, A_0_inv)
    ts.get_mns(M, 1, no-1, M)
    
    print("\nM:", M)
    prt_map(M, "M:", 1e-30)


def compute_R(desc, map, K, A_0, A_1, g, no):
    Id       = gtpsa.ss_vect_tpsa(desc, no)
    M        = gtpsa.ss_vect_tpsa(desc, no)
    R_inv    = gtpsa.ss_vect_tpsa(desc, no)
    A_0_inv  = gtpsa.ss_vect_tpsa(desc, no)
    A_1_inv  = gtpsa.ss_vect_tpsa(desc, no)
    A_nl     = gtpsa.ss_vect_tpsa(desc, no)
    A_nl_inv = gtpsa.ss_vect_tpsa(desc, no)

    Id.set_identity()

    M = _copy.copy(map)

    ts.h_DF_to_M(-1e0*K, Id, 2, no, R_inv)
    ts.get_mns(R_inv, 1, no-1, R_inv)
    A_0_inv.inv(A_0)
    A_1_inv.inv(A_1)
    ts.h_DF_to_M(g, Id, 3, no, A_nl)
    A_nl_inv.inv(A_nl)

    M.compose(M, A_0)
    M.compose(M, A_1)
    M.compose(M, A_nl)
    M.compose(M, R_inv)
    M.compose(A_0_inv, M)
    M.compose(A_1_inv, M)
    M.compose(A_nl_inv, M)
    ts.get_mns(M, 1, no-1, M)
   
    print("\nM:", M)
    prt_map(M, "M:", 1e-10)


# Number of phase-space coordinates.
nv = 7
# Max order for Poincaré map.
no = 4
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

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

M = compute_map(lat_prop, no)
print("\nM:", M)

K, A_0, A_1, g = compute_map_normal_form(desc, M)

if False:
    compute_h(desc, M)

if False:
    compute_M(desc, K, A_0, A_1, g, no)

if not False:
    compute_R(desc, M, K, A_0, A_1, g, no)
