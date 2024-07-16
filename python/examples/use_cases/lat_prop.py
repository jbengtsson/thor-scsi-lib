"""Use Case:
     Compute, print, and display the global linear optics properties for a
     lattice.
"""


import os

import numpy as np

import gtpsa

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    index_class as ind


ind = ind.index_class()


def compute_optics(lat_prop):
    try:
        # Compute Twiss parameters along lattice.
        if not lat_prop.comp_per_sol():
            print("\ncomp_per_sol: unstable")
            raise ValueError

        # Compute radiation properties.
        stable, stable_rad = lat_prop.compute_radiation()
        print(stable, stable_rad)
        if not stable:
            print("\ncompute_radiation: unstable")
            raise ValueError
    except ValueError:
        assert False


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


def compute_h(lat_prop, M):
    h    = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    h_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    h_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)

    M.M_to_h_DF(h)
    h.CtoR(h_re, h_im)
    return h_re, h_im


def compute_map_normal_form(lat_prop, M):
    A_0  = gtpsa.ss_vect_tpsa(lat_prop._desc, lat_prop._no)
    A_1  = gtpsa.ss_vect_tpsa(lat_prop._desc, lat_prop._no)
    R    = gtpsa.ss_vect_tpsa(lat_prop._desc, lat_prop._no)
    g    = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    g_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    g_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    K    = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    K_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    K_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)

    M.Map_Norm(A_0, A_1, R, g, K)
    K.CtoR(K_re, K_im)
    return A_0, A_1, R, g_re, g_im, K_re, K_im


h_dict = {
    "h_10002" : [1, 0, 0, 0, 2, 0, 0],
    "h_20001" : [2, 0, 0, 0, 1, 0, 0],
    "h_00201" : [0, 0, 2, 0, 1, 0, 0],

    "h_30000" : [3, 0, 0, 0, 0, 0, 0],
    "h_21000" : [2, 1, 0, 0, 0, 0, 0],
    "h_11100" : [1, 1, 1, 0, 0, 0, 0],
    "h_10200" : [1, 0, 2, 0, 0, 0, 0],
    "h_10020" : [1, 0, 0, 2, 0, 0, 0]
}

K_dict = {
    "K_22000" : [2, 2, 0, 0, 0, 0, 0],
    "K_11110" : [1, 1, 1, 1, 0, 0, 0],
    "K_00220" : [0, 0, 2, 2, 0, 0, 0],

    "K_11002" : [1, 1, 0, 0, 2, 0, 0],
    "K_00112" : [0, 0, 1, 1, 2, 0, 0]
}


def prt_nl(h_im, K_re):
    print()
    for key in h_dict:
        print("    {:s} = {:10.3e}".format(key, h_im.get(h_dict[key])))
        if key == "h_00201":
            print()
    print()
    for key in K_dict:
        print("    {:s} = {:10.3e}".format(key, K_re.get(K_dict[key])))
        if key == "K_00220":
            print()


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 4
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_4u")
lat_name = input("file_name?> ")
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.prt_lat(lat_name+"_lat.txt")

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

compute_optics(lat_prop)

lat_prop.prt_lat_param()
lat_prop.prt_rad()
lat_prop.prt_M()
lat_prop.prt_M_rad()

if not False:
    lat_prop.plt_Twiss(lat_name+"_Twiss.png", not False)
    lat_prop.plt_chrom(lat_name+"_chrom.png", not False)

twoJ = compute_twoJ(A_max, beta_inj)
Id_scl = compute_Id_scl(lat_prop, twoJ)

M = compute_map(lat_prop, no)
print("\nM:", M)
h_re, h_im = compute_h(lat_prop, M)
A_0, A_1, R, g_re, g_im, K_re, K_im = compute_map_normal_form(lat_prop, M)

h_im = compose_bs(h_im, Id_scl)
K_re = compose_bs(K_re, Id_scl)

prt_nl(h_im, K_re)
