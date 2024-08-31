import os
import sys

import copy as _copy

import enum as en
import numpy as np
from typing import ClassVar
from dataclasses import dataclass

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    courant_snyder as cs, index_class as ind
from thor_scsi.utils.output import mat2txt


class MpoleInd(en.IntEnum):
    quad = 2
    sext = 3

ind = ind.index_class()


def new_tpsa(gtpsa_prop):
    return gtpsa.tpsa(gtpsa_prop.desc, gtpsa_prop.no)


def new_ss_vect_tpsa(gtpsa_prop):
    return gtpsa.ss_vect_tpsa(gtpsa_prop.desc, gtpsa_prop.no)


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


def prt_map(map, str, *, eps: float=1e-30):
    print(str)
    map.x.print("x", eps)
    map.px.print("p_x", eps)
    map.y.print("y", eps)
    map.py.print("p_y", eps)
    map.delta.print("delta", eps)
    map.ct.print("ct", eps)


def compute_map(lat_prop, gtpsa_prop):
    M = lo.compute_map(
        lat_prop._lattice, lat_prop._model_state, desc=lat_prop._desc,
        tpsa_order=gtpsa_prop.no)
    return M


def compute_twoJ(A_max, beta_inj):
    twoJ = \
        np.array(
            [A_max[ind.X]**2/beta_inj[ind.X], A_max[ind.Y]**2/beta_inj[ind.Y]])
    return twoJ


def compute_Id_scl(gtpsa_prop, twoJ):
    Id_scl = new_ss_vect_tpsa(gtpsa_prop)
    Id_scl.set_identity()
    for k in range(4):
        Id_scl.iloc[k].set_variable(0e0, k+1, np.sqrt(twoJ[k//2]))
    Id_scl.delta.set_variable(0e0, 5, delta_max)
    return Id_scl


def compose_bs(gtpsa_prop, h, map):
    t_map = new_ss_vect_tpsa(gtpsa_prop)
    t_map.set_zero()
    t_map.x = h
    t_map.compose(t_map, map)
    return t_map.x 


def compute_h(gtpsa_prop, M):
    h = new_tpsa(gtpsa_prop)
    M.M_to_h_DF(h)
    return h


def compute_map_normal_form(gtpsa_prop, M):
    A_0 = new_ss_vect_tpsa(gtpsa_prop)
    A_1 = new_ss_vect_tpsa(gtpsa_prop)
    R   = new_ss_vect_tpsa(gtpsa_prop)
    g   = new_tpsa(gtpsa_prop)
    K   = new_tpsa(gtpsa_prop)
    M.Map_Norm(A_0, A_1, R, g, K)
    return A_0, A_1, R, g, K


def compute_M(gtpsa_prop, K, A_0, A_1, g):
    Id       = new_ss_vect_tpsa(gtpsa_prop)
    R        = new_ss_vect_tpsa(gtpsa_prop)
    A_0_inv  = new_ss_vect_tpsa(gtpsa_prop)
    A_1_inv  = new_ss_vect_tpsa(gtpsa_prop)
    A_nl     = new_ss_vect_tpsa(gtpsa_prop)
    A_nl_inv = new_ss_vect_tpsa(gtpsa_prop)
    M        = new_ss_vect_tpsa(gtpsa_prop)

    prt_map(M, "M:", eps=1e-30)

    Id.set_identity()

    K.h_DF_to_M(Id, 2, gtpsa_prop.no, R)
    A_0_inv.inv(A_0)
    A_1_inv.inv(A_1)
    g.h_DF_to_M(Id, 3, gtpsa_prop.no, A_nl)
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


def compute_R(desc, map, A_0, A_1, g, K, gtpsa_prop):
    Id       = new_ss_vect_tpsa(gtpsa_prop)
    M        = new_ss_vect_tpsa(gtpsa_prop)
    R        = new_ss_vect_tpsa(gtpsa_prop)
    R_inv    = new_ss_vect_tpsa(gtpsa_prop)
    A_0_inv  = new_ss_vect_tpsa(gtpsa_prop)
    A_1_inv  = new_ss_vect_tpsa(gtpsa_prop)
    A_nl     = new_ss_vect_tpsa(gtpsa_prop)
    A_nl_inv = new_ss_vect_tpsa(gtpsa_prop)

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


def compute_s_ijklm(lat_prop, i, j, k, l, m):
    eta = \
        np.array(lat_prop._Twiss.dispersion.sel(phase_coordinate="x").values)
    alpha = \
        np.array([lat_prop._Twiss.twiss.sel(plane="x", par="alpha").values,
                  lat_prop._Twiss.twiss.sel(plane="y", par="alpha").values]
    )
    beta = \
        np.array([lat_prop._Twiss.twiss.sel(plane="x", par="beta").values,
                  lat_prop._Twiss.twiss.sel(plane="y", par="beta").values]
    )
    dnu = np.array([lat_prop._Twiss.twiss.sel(plane="x", par="nu").values,
                    lat_prop._Twiss.twiss.sel(plane="y", par="nu").values])
    nu = np.array([lat_prop._Twiss.twiss.sel(plane="x", par="nu").values[-1],
                   lat_prop._Twiss.twiss.sel(plane="y", par="nu").values[-1]])

    s_ijklm = 0e0
    for n in range(len(lat_prop._lattice)):
        if type(lat_prop._lattice[n]) == ts.Sextupole:
            m_x = i + j
            m_y = k + l
            n_x = i - j
            n_y = k - l
            b3xL = \
                lat_prop._lattice[n].get_multipoles(). \
                get_multipole(MpoleInd.sext).real
            s_ijklm += \
                b3xL*beta[ind.X, n]**(m_x/2e0)*beta[ind.Y, n]**(m_y/2e0) \
                *np.cos(abs(n_x*dnu[ind.X, n]+n_y*dnu[ind.Y, n] \
                         -np.pi*(n_x*nu[ind.X]+n_y*nu[ind.Y]))) \
                /np.sin(np.pi*(n_x*nu[ind.X]+n_y*nu[ind.Y]))
    return s_ijklm


def compute_ampl_dep_orbit(lat_prop):
    s_21000 =  3e0/8e0*compute_s_ijklm(lat_prop, 2, 1, 0, 0, 0)
    s_30000 =  1e0/8e0*compute_s_ijklm(lat_prop, 3, 0, 0, 0, 0)
    s_10110 =  1e0/8e0*compute_s_ijklm(lat_prop, 1, 0, 1, 1, 0)
    s_10200 =  1e0/8e0*compute_s_ijklm(lat_prop, 2, 1, 0, 0, 0)
    s_10020 = -1e0/8e0*compute_s_ijklm(lat_prop, 2, 1, 0, 0, 0)

    print("\n  s_21000 = {:12.5e}".format(s_21000))
    print("  s_30000 = {:12.5e}".format(s_30000))
    print("  s_10110 = {:12.5e}".format(s_10110))
    print("  s_10200 = {:12.5e}".format(s_10200))
    print("  s_10020 = {:12.5e}".format(s_10020))

    return s_21000, s_30000, s_10110, s_10200, s_10020


def chk_Poisson_bracket(no, desc):
    Id = gtpsa.ss_vect_tpsa(desc, no)
    Id.set_identity()

    (1e0*Id.x).poisbra(1e0*Id.px, 6).print("[x, p_x]")
    (1e0*Id.px).poisbra(1e0*Id.x, 6).print("[p_x, x]")
    (1e0*Id.x).poisbra(1e0*Id.x, 6).print("[x, x]")
    (1e0*Id.px).poisbra(1e0*Id.px, 6).print("[p_x, p_x]")

    (1e0*Id.y).poisbra(1e0*Id.py, 6).print("[y, p_y]")
    (1e0*Id.y).poisbra(1e0*Id.y, 6).print("[y, y]")
    (1e0*Id.py).poisbra(1e0*Id.py, 6).print("[p_y, p_y]")

    (1e0*Id.x).poisbra(1e0*Id.y, 6).print("[x, y]")
    (1e0*Id.x).poisbra(1e0*Id.py, 6).print("[x, p_y]")


# def chk_compose(no, desc):

@dataclass
class gtpsa_prop_class:
    # GTPSA properties.
    # Number of phase-space coordinates.
    nv: ClassVar[int] = 7
    # Max order for Poincaré map.
    no: ClassVar[int] = 4
    # Number of parameters.
    nv_prm: ClassVar[int] = 0
    # Parameters max order.
    no_prm: ClassVar[int] = 0
    # Index.
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]


gtpsa_prop = gtpsa_prop_class()
gtpsa_prop_class.no = 4
print(gtpsa_prop.nv, gtpsa_prop.no)
gtpsa_prop_class.desc = gtpsa.desc(gtpsa_prop.nv, gtpsa_prop.no)

cod_eps = 1e-15
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

twoJ = compute_twoJ(A_max, beta_inj)
Id_scl = compute_Id_scl(gtpsa_prop, twoJ)

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

print("\nCircumference [m] = {:7.5f}".format(lat_prop.compute_circ()))

compute_optics(lat_prop)

lat_prop.prt_lat_param()
lat_prop.prt_M()
# lat_prop.prt_rad()
# lat_prop.prt_M_rad()

lat_prop._model_state.radiation = False
lat_prop._model_state.Cavity_on = False

if False:
    chk_Poisson_bracket(no, desc)

if False:
    chk__compose(no, desc)

M = compute_map(lat_prop, gtpsa_prop)
print("\nM:", M)

if False:
    h_re = new_tpsa(gtpsa_prop)
    h_im = new_tpsa(gtpsa_prop)

    h = compute_h(desc, M)
    h.CtoR(h_re, h_im)
    h_re = compose_bs(no, desc, h_re, Id_scl)
    h_im = compose_bs(no, desc, h_im, Id_scl)

if not False:
    A_0, A_1, R, g, K = compute_map_normal_form(gtpsa_prop, M)

if False:
    compute_M(gtpsa_prop, K, A_0, A_1, g)
if False:
    compute_R(gtpsa_prop, M, A_0, A_1, g, K)

g_re = new_tpsa(gtpsa_prop)
g_im = new_tpsa(gtpsa_prop)
K_re = new_tpsa(gtpsa_prop)
K_im = new_tpsa(gtpsa_prop)

g.CtoR(g_re, g_im)
K.CtoR(K_re, K_im)

if not False:
    g_re = compose_bs(gtpsa_prop, g_re, Id_scl)
    g_im = compose_bs(gtpsa_prop, g_im, Id_scl)
    K_re = compose_bs(gtpsa_prop, K_re, Id_scl)
    K_im = compose_bs(gtpsa_prop, K_im, Id_scl)

g_im.print("g_im")

g_im = g_im.get_mns_1(1, gtpsa_prop.no-1)

# g_re is zero.
# g_re.print("g_re")
g_im.print("g_im")
K_re.print("K_re")
# K_im is zero.
# K_im.print("K_im:")
print("\nA_1:", A_1)

compute_ampl_dep_orbit(lat_prop)

Id = new_ss_vect_tpsa(gtpsa_prop)
Id.set_identity()
x_mean = (1e0*g_im).poisbra(1e0*Id.x, 6)
x_mean.print("\n<x>:")

t_map = new_ss_vect_tpsa(gtpsa_prop)
t_map.set_zero()
t_map.x = g_im
assert False
print("\nt_map:", t_map)
t_map = t_map.compose(A_1, t_map)
x_mean = t_map.iloc[ind.x]
x_mean.print("\n<x>")
