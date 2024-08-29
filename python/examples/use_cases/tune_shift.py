import os
import sys

import copy as _copy

import enum
import numpy as np

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    courant_snyder as cs, index_class as ind
from thor_scsi.utils.output import mat2txt


class MpoleInd(enum.IntEnum):
    quad = 2
    sext = 3

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
    h = gtpsa.tpsa(desc, no)
    M.M_to_h_DF(h)
    return h


def compute_map_normal_form(desc, M):
    A_0 = gtpsa.ss_vect_tpsa(desc, no)
    A_1 = gtpsa.ss_vect_tpsa(desc, no)
    R   = gtpsa.ss_vect_tpsa(desc, no)
    g   = gtpsa.tpsa(desc, no)
    K   = gtpsa.tpsa(desc, no)
    M.Map_Norm(A_0, A_1, R, g, K)
    return A_0, A_1, R, g, K


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


# Number of phase-space coordinates.
nv = 7
# Max order for Poincaré map.
no = 4
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

desc = gtpsa.desc(nv, no)

cod_eps = 1e-15
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

print("\nCircumference [m] = {:7.5f}".format(lat_prop.compute_circ()))

compute_optics(lat_prop)

lat_prop.prt_lat_param()
lat_prop.prt_rad()
lat_prop.prt_M()
lat_prop.prt_M_rad()

lat_prop._model_state.radiation = False
lat_prop._model_state.Cavity_on = False

twoJ = compute_twoJ(A_max, beta_inj)
Id_scl = compute_Id_scl(lat_prop, twoJ)

M = compute_map(lat_prop, no)
print("\nM:", M)

h = compute_h(desc, M)
A_0, A_1, R, g, K = compute_map_normal_form(desc, M)

if False:
    compute_M(desc, K, A_0, A_1, g, no)

if False:
    compute_R(desc, M, A_0, A_1, g, K, no)

# h_im = compose_bs(h_im, Id_scl)
# K_re = compose_bs(K_re, Id_scl)

h_re = gtpsa.tpsa(desc, no)
h_im = gtpsa.tpsa(desc, no)
g_re = gtpsa.tpsa(desc, no)
g_im = gtpsa.tpsa(desc, no)
K_re = gtpsa.tpsa(desc, no)
K_im = gtpsa.tpsa(desc, no)

h.CtoR(h_re, h_im)
g.CtoR(g_re, g_im)
K.CtoR(K_re, K_im)

g_im = g_im.get_mns_1(1, no-1)
g_im.print("g_im")
K_re.print("K_re")

Id = gtpsa.ss_vect_tpsa(desc, no)
Id.set_identity()

r1 = Id.x**2
r1.print("\nr1")
r2 = r1.deriv(1)
r2.print("\nr2")
r2 = r2.integ(2)
r2.print("\nr2")

(Id.y**2).deriv(3).print("\nr2")

r3 = Id.x
r4 = Id.px
r3.poisbra(r4, 6).print("\n[]")

compute_ampl_dep_orbit(lat_prop)
