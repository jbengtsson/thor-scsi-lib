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

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    index_class as ind
from thor_scsi.utils.output import mat2txt


ind = ind.index_class()


def prt_map(map, str, *, eps: float=1e-30):
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


def compute_twoJ(A_max, beta_inj):
    twoJ = \
        np.array(
            [A_max[ind.X]**2/beta_inj[ind.X], A_max[ind.Y]**2/beta_inj[ind.Y]])
    return twoJ


def compute_Id_scl(lat_prop, twoJ):
    Id_scl = \
        gtpsa.ss_vect_tpsa(
            lat_prop._desc, lat_prop._no, index_mapping=lat_prop._named_index)
    Id_scl.set_identity()
    for k in range(4):
        Id_scl.iloc[k].set_variable(0e0, k+1, np.sqrt(twoJ[k//2]))
    Id_scl.delta.set_variable(0e0, 5, delta_max)
    return Id_scl


def compose_bs(h, map):
    Id = \
        gtpsa.ss_vect_tpsa(
        lat_prop._desc, lat_prop._no, index_mapping=lat_prop._named_index)
    t_map = \
        gtpsa.ss_vect_tpsa(
        lat_prop._desc, lat_prop._no, index_mapping=lat_prop._named_index)
    t_map.x = h
    t_map.compose(t_map, map)
    return t_map.x 


def compute_h(desc, M):
    h    = gtpsa.tpsa(desc, no)
    h_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    h_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    M.M_to_h_DF(h)
    h.CtoR(h_re, h_im)
    return h_re, h_im


def compute_map_normal_form(desc, M):
    A_0  = gtpsa.ss_vect_tpsa(desc, no)
    A_1  = gtpsa.ss_vect_tpsa(desc, no)
    R    = gtpsa.ss_vect_tpsa(desc, no)
    g    = gtpsa.tpsa(desc, no)
    g_re = gtpsa.tpsa(desc, no)
    g_im = gtpsa.tpsa(desc, no)
    K    = gtpsa.tpsa(desc, no)
    K_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    K_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)

    M.Map_Norm(A_0, A_1, R, g, K)
    K.CtoR(K_re, K_im)
    return A_0, A_1, R, g_re, g_im, K_re, K_im


def compute_M(desc, K, A_0, A_1, g, no):
    Id       = gtpsa.ss_vect_tpsa(desc, no)
    R        = gtpsa.ss_vect_tpsa(desc, no)
    A_0_inv  = gtpsa.ss_vect_tpsa(desc, no)
    A_1_inv  = gtpsa.ss_vect_tpsa(desc, no)
    A_nl     = gtpsa.ss_vect_tpsa(desc, no)
    A_nl_inv = gtpsa.ss_vect_tpsa(desc, no)
    M        = gtpsa.ss_vect_tpsa(desc, no)

    prt_map(M, "M:", eps=1e-30)

    Id.set_identity()

    K.h_DF_to_M(Id, 2, no, R)
    A_0_inv.inv(A_0)
    A_1_inv.inv(A_1)
    g.h_DF_to_M(Id, 3, no, A_nl)
    A_nl_inv.inv(A_nl)

    M.compose(A_0, A_1)
    M.compose(M, A_nl)
    M.compose(M, R)
    M.compose(M, A_nl_inv)
    M.compose(M, A_1_inv)
    M.compose(M, A_0_inv)
    M.get_mns(1, no-1, M)
    
    print("\nM:", M)
    prt_map(M, "M:", eps=1e-30)


def compute_R(desc, map, A_0, A_1, g, K, no):
    Id       = gtpsa.ss_vect_tpsa(desc, no)
    M        = gtpsa.ss_vect_tpsa(desc, no)
    R        = gtpsa.ss_vect_tpsa(desc, no)
    R_inv    = gtpsa.ss_vect_tpsa(desc, no)
    A_0_inv  = gtpsa.ss_vect_tpsa(desc, no)
    A_1_inv  = gtpsa.ss_vect_tpsa(desc, no)
    A_nl     = gtpsa.ss_vect_tpsa(desc, no)
    A_nl_inv = gtpsa.ss_vect_tpsa(desc, no)

    Id.set_identity()

    M = _copy.copy(map)

    K.h_DF_to_M(Id, 2, no, R)
    R_inv.inv(R)
    R_inv.get_mns(1, no-1, R_inv)
    A_0_inv.inv(A_0)
    A_1_inv.inv(A_1)
    g.h_DF_to_M(Id, 3, no, A_nl)
    A_nl_inv.inv(A_nl)

    M.compose(M, A_0)
    M.compose(M, A_1)
    M.compose(M, A_nl)
    M.compose(M, R_inv)
    M.compose(A_0_inv, M)
    M.compose(A_1_inv, M)
    M.compose(A_nl_inv, M)
    M.get_mns(1, no-1, M)
   
    print("\nM:", M)
    prt_map(M, "M:", eps=1e-10)


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

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_4u")
lat_name = "max_4u_g_1"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

print("\nCircumference [m] = {:7.5f}".format(lat_prop.compute_circ()))

twoJ = compute_twoJ(A_max, beta_inj)
Id_scl = compute_Id_scl(lat_prop, twoJ)

M = compute_map(lat_prop, no)
print("\nM:", M)

h_re, h_im = compute_h(desc, M)
A_0, A_1, R, g_re, g_im, K_re, K_im = compute_map_normal_form(desc, M)

h_im = compose_bs(h_im, Id_scl)
K_re = compose_bs(K_re, Id_scl)

h_im.print("h_im")
K_re.print("K_re")

if False:
    compute_M(desc, K, A_0, A_1, g, no)

if False:
    compute_R(desc, M, A_0, A_1, g, K, no)
