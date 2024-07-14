"""Use Case:
     Sextupole response matrix analysis.

     J. Bengtsson 
     𝑇ℎ𝑒 𝑆𝑒𝑥𝑡𝑢𝑝𝑜𝑙𝑒 𝑆𝑐ℎ𝑒𝑚𝑒 𝑓𝑜𝑟 𝑡ℎ𝑒 𝑆𝑤𝑖𝑠𝑠 𝐿𝑖𝑔ℎ𝑡 𝑆𝑜𝑢𝑟𝑐𝑒 (𝑆𝐿𝑆) – 𝐴𝑛 𝐴𝑛𝑎𝑙𝑦𝑡𝑖𝑐 𝐴𝑝𝑝𝑟𝑜𝑎𝑐ℎ
     See SLS Tech Note 09/97

     https://ados.web.psi.ch/slsnotes/sls0997.pdf
"""

import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

from scipy import optimize

import copy as _copy
from dataclasses import dataclass
import math
import os
from typing import Sequence
from typing import Tuple

import numpy as np
import xarray as xr

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    index_class as ind, knobs

from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv

from thor_scsi.utils.extract_info import accelerator_info


# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt


ind = ind.index_class()


class MultipoleIndex(enum.IntEnum):
    quadrupole = 2
    sextupole  = 3


def compute_optics(lat_prop):
    try:
        # Compute Twiss parameters along lattice.
        if not lat_prop.comp_per_sol():
            print("\ncomp_per_sol: unstable")
            raise ValueError

        # Compute radiation properties.
        if not lat_prop.compute_radiation():
            print("\ncompute_radiation: unstable")
            raise ValueError
    except ValueError:
        assert False

    lat_prop.prt_lat_param()
    lat_prop.prt_rad()
    lat_prop.prt_M()
    lat_prop.prt_M_rad()


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


def compute_twoJ(A_max, beta_inv):
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


def prt_nl(h_im, K_re):
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

    print()
    for key in h_dict:
        print("  {:s} = {:10.3e}".format(key, h_im.get(h_dict[key])))
    print()
    for key in K_dict:
        print("  {:s} = {:10.3e}".format(key, K_re.get(K_dict[key])))


def get_phi_elem(lat, fam_name, n_kid):
    elem = lat.find(fam_name, n_kid)
    return elem.get_length() * elem.get_curvature() * 180e0 / np.pi


def compute_phi(lat):
    """Compute total bend angle.
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


def set_b_n_elem(lat, fam_name, kid_num, n, b_n):
    mp = lat.find(fam_name, kid_num)
    mp.get_multipoles().set_multipole(n, b_n)


def set_b_n_fam(lat, fam_name, n, b_n):
    for mp in lat.elements_with_name(fam_name):
        mp.get_multipoles().set_multipole(n, b_n)


corresponding_types = {
    ts.Sextupole:         ts.SextupoleTpsa,
    ts.Quadrupole:        ts.QuadrupoleTpsa,
    ts.HorizontalSteerer: ts.HorizontalSteererTpsa,
    ts.VerticalSteerer:   ts.VerticalSteererTpsa,
}


def convert_magnet_to_knobbable(a_magnet: ts.Mpole) -> ts.MpoleTpsa:
    config = a_magnet.config()
    corresponding_type = corresponding_types[type(a_magnet)]
    return corresponding_type(config)


def mult_prm(elems, mult_family, n, desc):
    print("\nmult_prm ->")
    for k in range(len(mult_family)):
        index = mult_family[k].index
        elem = convert_magnet_to_knobbable(elems[index])
        elems[index] = \
            knobs.make_magnet_knobbable(
                elem, po=1, desc=desc, named_index=named_index,
                multipole_number=n, offset=False
            )
        print("\nmult_prm:", elems[index])
        # While the RHS pointer can be recasted to:
        #   CellVoid
        #   ElemTypeKnobbed
        #   QuadrupoleType
        # the assignment of the LHS only gives:
        #   CellVoid
        # and the lack of:
        #   ElemtypeKnobbed
        # will lead to an assert on line 158 in:
        #   thor_scsi/std_machine/accelerator.cc
        #
    print("-> mult_prm\n")
    return elems


def lat_mult_prm(mult_prm_name, lat, n, desc):
    elems = [elem for elem in lat]
    mult_family = lat.elements_with_name(mult_prm_name)
    elems = mult_prm(elems, mult_family, n, desc)
    return ts.AcceleratorTpsa(elems)


# Work-around for C++ virtual function -> Python mapping issue.
# See function above:
#   mult_prm
def propagate(lat_prop, named_index):
    M = gtpsa.ss_vect_tpsa \
        (lat_prop._desc, lat_prop._no, lat_prop._nv, lat_prop._named_index)
    M.set_identity()
    for k in range(len(lat_prop._lattice)):
        lat_prop._lattice[k].propagate(model_state, M)
    return M


def param_dep():
    if False:
        # Zero sextopoles.
        print("\nZeroing sextupoles.")
        set_b_n_fam(lat, "sf_h", MultipoleIndex.sextupole, 0e0)
        set_b_n_fam(lat, "sd1", MultipoleIndex.sextupole, 0e0)
        set_b_n_fam(lat, "sd2", MultipoleIndex.sextupole, 0e0)

    if False:
        named_index = gtpsa.IndexMapping(
            dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
        desc = gtpsa.desc(nv, no, nv_prm, no_prm)

        print("\n  C [m]     = {:9.7f}".format(C))
        print("  phi [deg] = {:5.3f}".format(compute_phi(lat)))
        print("  nu        = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))
        print("  xi        = [{:7.5f}, {:7.5f}]".format(xi[X_], xi[Y_]))

    named_index = gtpsa.IndexMapping(
        dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6)
    )

    if False:
        nv_prm = 1
        no_prm = no
        desc = gtpsa.desc(nv, no, nv_prm, no_prm)
        if True:
            lat_ptc = lat_mult_prm("sf_h", lat_prop._lattice, 3, desc)
        else:
            lat_ptc = lat_mult_prm("uq3", lat_prop._lattice, 2, desc)
    else:
        nv_prm = 0
        no_prm = 0
        desc = gtpsa.desc(nv, no, nv_prm, no_prm)
        lat_ptc = lat_mult_prm("", lat_prop._lattice, 0, desc)

    M = propagate(lat_ptc, named_index)

    print("\nM:\n", M)
    if False:
        print_map("\nM:", M)

    if False:
        M2 = M
        # Busted.
        M2.getOrder(M2, 2)
        print("\n:\n", M2)

    if False:
        M_inv = ts.inv(M)

        M_M_inv = gtpsa.ss_vect_tpsa(desc, no, nv, index_mapping=named_index)
        M_M_inv.rcompose(M, M_inv)
        print("\nM:", M)
        print("\nM^-1:", M_inv)
        print("\nM*M^-1:", M_M_inv)
        print("\nM.x:\n", M.x)
        print("\nM_M_inv.x:\n", M_M_inv.x)
        # print("\nM*M^-1:", M_M_inv[0])
        assert False


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 4
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

n_dof     = 2
cod_eps   = 1e-15
E_0       = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_4u",
    "")
lat_name = "max_4u_g_0"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.prt_lat(lat_name+"_lat.txt")

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

if not False:
    compute_optics(lat_prop)

    if False:
        lat_prop.plt_Twiss(lat_name+"_Twiss.png", not False)
        lat_prop.plt_chrom(lat_name+"_chrom.png", not False)

M = compute_map(lat_prop, no)
print("\nM:", M)

twoJ = compute_twoJ(A_max, beta_inj)
Id_scl = compute_Id_scl(lat_prop, twoJ)

h_re, h_im = compute_h(lat_prop, M)
A_0, A_1, R, g_re, g_im, K_re, K_im = compute_map_normal_form(lat_prop, M)

if False:
    # The real part of the driving terms are cancelled by mirror symmetry.
    # The imaginary part of tune shift terms are zero.
    h_im.print("h_im", 1e-5)
    K_re.print("K_re", 1e-5)

prt_nl(h_im, K_re)
h_im = compose_bs(h_im, Id_scl)
K_re = compose_bs(K_re, Id_scl)
prt_nl(h_im, K_re)
