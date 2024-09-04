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


@dataclass
class gtpsa_prop_class:
    # GTPSA properties.
    # Number of phase-space coordinates.
    nv: ClassVar[int] = 7
    # Max order for Poincaré map.
    no: ClassVar[int] = 1
    # Number of parameters.
    nv_prm: ClassVar[int] = 0
    # Parameters max order.
    no_prm: ClassVar[int] = 0
    # Index.
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]


class new():
    def tpsa():
        return gtpsa.tpsa(gtpsa_prop.desc, gtpsa_prop.no)
    def ss_vect_tpsa():
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


def compute_map(lat_prop):
    M = lo.compute_map(
        lat_prop._lattice, lat_prop._model_state, desc=lat_prop._desc,
        tpsa_order=gtpsa_prop.no)
    return M


def compute_twoJ(A_max, beta_inj):
    twoJ = \
        np.array(
            [A_max[ind.X]**2/beta_inj[ind.X], A_max[ind.Y]**2/beta_inj[ind.Y]])
    return twoJ


def compute_Id_scl(twoJ):
    Id_scl = new.ss_vect_tpsa()
    Id_scl.set_identity()
    for k in range(4):
        Id_scl.iloc[k].set_variable(0e0, k+1, np.sqrt(twoJ[k//2]))
    Id_scl.delta.set_variable(0e0, 5, delta_max)
    return Id_scl


def compose_bs(h, map):
    t_map = new.ss_vect_tpsa()
    t_map.set_zero()
    t_map.x = h
    t_map.compose(t_map, map)
    return t_map.x.to_tpsa()


def compute_h(M):
    h = new.tpsa()
    M.M_to_h_DF(h)
    return h


def compute_map_normal_form(M):
    A_0 = new.ss_vect_tpsa()
    A_1 = new.ss_vect_tpsa()
    R   = new.ss_vect_tpsa()
    g   = new.tpsa()
    K   = new.tpsa()
    M.Map_Norm(A_0, A_1, R, g, K)
    return A_0, A_1, R, g, K


def compute_M(K, A_0, A_1, g):
    Id       = new.ss_vect_tpsa()
    R        = new.ss_vect_tpsa()
    A_0_inv  = new.ss_vect_tpsa()
    A_1_inv  = new.ss_vect_tpsa()
    A_nl     = new.ss_vect_tpsa()
    A_nl_inv = new.ss_vect_tpsa()
    M        = new.ss_vect_tpsa()

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
    Id       = new.ss_vect_tpsa()
    M        = new.ss_vect_tpsa()
    R        = new.ss_vect_tpsa()
    R_inv    = new.ss_vect_tpsa()
    A_0_inv  = new.ss_vect_tpsa()
    A_1_inv  = new.ss_vect_tpsa()
    A_nl     = new.ss_vect_tpsa()
    A_nl_inv = new.ss_vect_tpsa()

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


def get_h_ijklm(g, i, j, k, l, m):
    return g.get([i, j, k, l, m, 0, 0])


def compute_h_ijklm(lat_prop, i, j, k, l, m):
    eta = np.array(lat_prop._Twiss.dispersion.sel(phase_coordinate="x").values)
    alpha = np.array(
        [lat_prop._Twiss.twiss.sel(plane="x", par="alpha").values,
         lat_prop._Twiss.twiss.sel(plane="y", par="alpha").values]
    )
    beta = np.array([lat_prop._Twiss.twiss.sel(plane="x", par="beta").values,
                     lat_prop._Twiss.twiss.sel(plane="y", par="beta").values]
    )
    dnu = np.array([lat_prop._Twiss.twiss.sel(plane="x", par="nu").values,
                    lat_prop._Twiss.twiss.sel(plane="y", par="nu").values])
    dmu = 2e0*np.pi*dnu
    nu = np.array([lat_prop._Twiss.twiss.sel(plane="x", par="nu").values[-1],
                   lat_prop._Twiss.twiss.sel(plane="y", par="nu").values[-1]])

    h_re = 0e0
    h_im = 0e0
    for n in range(len(lat_prop._lattice)):
        if type(lat_prop._lattice[n]) == ts.Sextupole:
            m_x = i + j
            m_y = k + l
            n_x = i - j
            n_y = k - l
            L = lat_prop._lattice[n].get_length()
            b3xL = lat_prop._lattice[n].get_multipoles(). \
                get_multipole(MpoleInd.sext).real*L
            if b3xL != 0e0:
                phi = n_x*dmu[ind.X, n-1] + n_y*dmu[ind.Y, n-1]
                h_abs = b3xL*beta[ind.X, n-1]**(m_x/2e0) \
                    *beta[ind.Y, n-1]**(m_y/2e0)
                dh_re = h_abs*np.cos(phi)
                dh_im = h_abs*np.sin(phi)
                h_re += h_abs*np.cos(phi)
                h_im -= h_abs*np.sin(phi)

    return np.array([h_re, h_im])


def compute_hh(lat_prop, h):
    h_re = new.tpsa()
    h_im = new.tpsa()

    h = h.get_mns_1(1, 3)
    h.CtoR(h_re, h_im)
    h_re.print("h_re", 1e-11)
    h_im.print("h_im", 1e-11)

    h_21000 = -1e0/8e0*compute_h_ijklm(lat_prop, 2, 1, 0, 0, 0)
    h_30000 = -1e0/24e0*compute_h_ijklm(lat_prop, 3, 0, 0, 0, 0)
    h_10110 =  1e0/4e0*compute_h_ijklm(lat_prop, 1, 0, 1, 1, 0)
    h_10200 =  1e0/8e0*compute_h_ijklm(lat_prop, 1, 0, 2, 0, 0)
    h_10020 =  1e0/8e0*compute_h_ijklm(lat_prop, 1, 0, 0, 2, 0)

    h_21000_ref = [get_h_ijklm(h_re, 2, 1, 0, 0, 0),
                   get_h_ijklm(h_im, 2, 1, 0, 0, 0)]
    h_30000_ref = [get_h_ijklm(h_re, 3, 0, 0, 0, 0),
                   get_h_ijklm(h_im, 3, 0, 0, 0, 0)]
    h_10110_ref = [get_h_ijklm(h_re, 1, 0, 1, 1, 0),
                   get_h_ijklm(h_im, 1, 0, 1, 1, 0)]
    h_10200_ref = [get_h_ijklm(h_re, 1, 0, 2, 0, 0),
                   get_h_ijklm(h_im, 1, 0, 2, 0, 0)]
    h_10020_ref = [get_h_ijklm(h_re, 1, 0, 0, 2, 0),
                   get_h_ijklm(h_im, 1, 0, 0, 2, 0)]

    print(f"\n  h_21000 = [{h_21000[0]:12.5e}, {h_21000[1]:12.5e}] ([{h_21000_ref[0]/2e0:12.5e}, {h_21000_ref[1]/2e0:12.5e}]))")
    print(f"  h_30000 = [{h_30000[0]:12.5e}, {h_30000[1]:12.5e}] ([{h_30000_ref[0]/2e0:12.5e}, {h_30000_ref[1]/2e0:12.5e}]))")
    print(f"  h_10110 = [{h_10110[0]:12.5e}, {h_10110[1]:12.5e}] ([{h_10110_ref[0]/2e0:12.5e}, {h_10110_ref[1]/2e0:12.5e}]))")
    print(f"  h_10200 = [{h_10200[0]:12.5e}, {h_10200[1]:12.5e}] ([{h_10200_ref[0]/2e0:12.5e}, {h_10200_ref[1]/2e0:12.5e}]))")
    print(f"  h_10020 = [{h_10020[0]:12.5e}, {h_10020[1]:12.5e}] ([{h_10020_ref[0]/2e0:12.5e}, {h_10020_ref[1]/2e0:12.5e}]))")


def compute_g_ijklm(lat_prop, i, j, k, l, m):
    eta = np.array(lat_prop._Twiss.dispersion.sel(phase_coordinate="x").values)
    alpha = np.array(
        [lat_prop._Twiss.twiss.sel(plane="x", par="alpha").values,
         lat_prop._Twiss.twiss.sel(plane="y", par="alpha").values]
    )
    beta = np.array([lat_prop._Twiss.twiss.sel(plane="x", par="beta").values,
                     lat_prop._Twiss.twiss.sel(plane="y", par="beta").values]
    )
    dnu = np.array([lat_prop._Twiss.twiss.sel(plane="x", par="nu").values,
                    lat_prop._Twiss.twiss.sel(plane="y", par="nu").values])
    dmu = 2e0*np.pi*dnu
    nu = np.array([lat_prop._Twiss.twiss.sel(plane="x", par="nu").values[-1],
                   lat_prop._Twiss.twiss.sel(plane="y", par="nu").values[-1]])

    g_re_ijklm = 0e0
    g_im_ijklm = 0e0
    for n in range(len(lat_prop._lattice)):
        if type(lat_prop._lattice[n]) == ts.Sextupole:
            m_x = i + j
            m_y = k + l
            n_x = i - j
            n_y = k - l
            L = lat_prop._lattice[n].get_length()
            b3xL = \
                lat_prop._lattice[n].get_multipoles(). \
                get_multipole(MpoleInd.sext).real*L
            g_abs = b3xL*beta[ind.X, n]**(m_x/2e0)*beta[ind.Y, n]**(m_y/2e0)
            phi = n_x*dnu[ind.X, n] + n_y*dnu[ind.Y, n]
            g_re_ijklm += g_abs*np.cos(phi) \
                /np.sin(np.pi*(n_x*nu[ind.X]+n_y*nu[ind.Y]))
            g_im_ijklm += g_abs*np.sin(phi) \
                /np.sin(np.pi*(n_x*nu[ind.X]+n_y*nu[ind.Y]))
    return g_im_ijklm, g_re_ijklm

def compute_g(lat_prop, g):
    g_re = new.tpsa()
    g_im = new.tpsa()

    g.CtoR(g_re, g_im)

    g_21000 =  1e0/8e0*compute_g_ijklm(lat_prop, 2, 1, 0, 0, 0)[0]
    g_30000 =  1e0/48e0*compute_g_ijklm(lat_prop, 3, 0, 0, 0, 0)[0]
    g_10110 = -1e0/8e0*compute_g_ijklm(lat_prop, 1, 0, 1, 1, 0)[0]
    g_10200 = -1e0/16e0*compute_g_ijklm(lat_prop, 1, 0, 2, 0, 0)[0]
    g_10020 = -1e0/16e0*compute_g_ijklm(lat_prop, 1, 0, 0, 2, 0)[0]

    g_21000_ref = get_h_ijklm(g_im, 2, 1, 0, 0, 0)
    g_30000_ref = get_h_ijklm(g_im, 3, 0, 0, 0, 0)
    g_10110_ref = get_h_ijklm(g_im, 1, 0, 1, 1, 0)
    g_10200_ref = get_h_ijklm(g_im, 1, 0, 2, 0, 0)
    g_10020_ref = get_h_ijklm(g_im, 1, 0, 0, 2, 0)

    print(f"\n  g_21000 = {g_21000:12.5e} ({g_21000_ref:12.5e})")
    print(f"  g_30000 = {g_30000:12.5e} ({g_30000_ref:12.5e})")
    print(f"  g_10110 = {g_10110:12.5e} ({g_10110_ref:12.5e})")
    print(f"  g_10200 = {g_10200:12.5e} ({g_10200_ref:12.5e}")
    print(f"  g_10020 = {g_10020:12.5e} ({g_10020_ref:12.5e}")

    return g_21000, g_30000, g_10110, g_10200, g_10020


def chk_types():
    Id = new.ss_vect_tpsa()
    Id.set_identity()

    print("\nId.x.to_tpsa():\n  ", type(Id.x), " ",
          type(Id.iloc[ind.x]), sep="")
    print(" ", type(Id.x.to_tpsa()))

    f = new.tpsa()
    g = new.tpsa()

    f = Id.x.to_tpsa()
    f.print("f")
    print("\nf.integ(1):\n  ", type(f), sep="")
    g = f.integ(1)
    g.print("g")
    print(" ", type(f), type(g))

    f = new.tpsa()
    g = new.tpsa()

    print("\nf+g:\n  ", type(f), " ", type(g), sep="")
    h = f + g
    print(" ", type(f+g), type(h))

    f = new.tpsa()

    print("\nf.get_mns_1(1, gtpsa_prop.no-1):\n  ", type(f), sep="")
    f = f.get_mns_1(1, gtpsa_prop.no-1)
    print(" ", type(f))


def chk_Poisson_bracket():
    Id = new.ss_vect_tpsa()
    Id.set_identity()

    ps_dim = 6

    Id.x.to_tpsa().poisbra(Id.px.to_tpsa(), ps_dim).print("[x, p_x]")
    Id.px.to_tpsa().poisbra(Id.x.to_tpsa(), ps_dim).print("[p_x, x]")
    Id.x.to_tpsa().poisbra(Id.x.to_tpsa(), ps_dim).print("[x, x]")
    Id.px.to_tpsa().poisbra(Id.px.to_tpsa(), ps_dim).print("[p_x, p_x]")

    Id.y.to_tpsa().poisbra(Id.py.to_tpsa(), ps_dim).print("[y, p_y]")
    Id.y.to_tpsa().poisbra(Id.y.to_tpsa(), ps_dim).print("[y, y]")
    Id.py.to_tpsa().poisbra(Id.py.to_tpsa(), ps_dim).print("[p_y, p_y]")

    Id.x.to_tpsa().poisbra(Id.y.to_tpsa(), ps_dim).print("[x, y]")
    Id.x.to_tpsa().poisbra(Id.py.to_tpsa(), ps_dim).print("[x, p_y]")

    pow(Id.x.to_tpsa(), 3).poisbra(Id.x.to_tpsa(), ps_dim). \
        print("[x^3, x]", 1e-11)
    pow(Id.x.to_tpsa(), 3).poisbra(Id.px.to_tpsa(), ps_dim). \
        print("[x^3, p_x]", 1e-11)
    pow(Id.x.to_tpsa(), 3).poisbra(pow(Id.x.to_tpsa(), 2)
                                   *Id.px.to_tpsa(), ps_dim). \
        print("[x^3,  x*p_x^2]", 1e-11)


def chk_compose():
    f = new.tpsa()
    g = new.tpsa()
    map_1 = new.ss_vect_tpsa()
    map_2 = new.ss_vect_tpsa()

    map_1.set_identity()
    map_2.set_identity()
    print("\n", type(map_1), type(map_2))
    map_1.compose(map_1, map_2)
    print("\n", type(map_1), type(map_2))

    f.set([1, 0, 0, 0, 0, 0, 0], 0e0, 1e0)
    print("\n", type(f), type(f))
    g = compose_bs(f, map_1)
    print("\n", type(f), type(g))

    map_2.exppb(map_1, map_2)
    print("\nmap_2", map_2)


def pb(f, g):
    a = new.tpsa()
    ps_dim = 6
    for k in range(0, 1):
        print("k =", k)
        fp = f.deriv(2*k+2)
        gp = g.deriv(2*k+1)
        print(fp.get_description())
        print(gp.get_description())
        fp.print("fp", 1e-11)
        gp.print("gp", 1e-11)
        fxg = fp*gp
        fxg.print("fxg", 1e-11)
        a += f.deriv(2*k+1)*g.deriv(2*k+2) - f.deriv(2*k+2)*g.deriv(2*k+1)
    return a


def compute_ampl_dep_orbit_tpsa(g):
    x_avg       = new.tpsa()
    x_avg_re    = new.tpsa()
    x_avg_im    = new.tpsa()
    x3_avg      = new.tpsa()
    x3_avg_re   = new.tpsa()
    x3_avg_im   = new.tpsa()
    x_y2_avg    = new.tpsa()
    x_y2_avg_re = new.tpsa()
    x_y2_avg_im = new.tpsa()

    Id          = new.ss_vect_tpsa()

    ps_dim = 6

    Id.set_identity()

    g = g.get_mns_1(1, 3)
    g.print("g", 1e-11)

    x_avg = g.poisbra(Id.x.to_tpsa(), ps_dim)
    x_avg = compose_bs(x_avg, A_1)
    x_avg = x_avg.get_mns_1(1, 2)
    x_avg.CtoR(x_avg_re, x_avg_im)
    x_avg_re.print("<x_re>", 1e-11)
    s_11000 = x_avg_re.get([1, 1, 0, 0, 0, 0, 0])
    s_00110 = x_avg_re.get([0, 0, 1, 1, 0, 0, 0])

    x3_avg = g.poisbra(pow(Id.x.to_tpsa(), 3), ps_dim)
    x3_avg.print("<x3>", 1e-11)
    x3_avg = compose_bs(x3_avg, A_1)
    x3_avg.print("<x3>", 1e-11)
    x3_avg.CtoR(x3_avg_re, x3_avg_im)
    x3_avg_re.print("<x^3_re>")
    s_22000 = x3_avg_re.get([2, 2, 0, 0, 0, 0, 0])
    s_11110_1 = x3_avg_re.get([1, 1, 1, 1, 0, 0, 0])

    x_y2_avg = g.poisbra(Id.x.to_tpsa()*pow(Id.y.to_tpsa(), 2), ps_dim)
    x_y2_avg = compose_bs(x_y2_avg, A_1)
    x_y2_avg.CtoR(x_y2_avg_re, x_y2_avg_im)
    x_y2_avg_re.print("<x*y^2_re>")
    s_11110_2 = x_y2_avg_re.get([1, 1, 1, 1, 0, 0, 0])
    s_00220 = x_y2_avg_re.get([0, 0, 2, 2, 0, 0, 0])

    print(f"\n  s_11000   = {s_11000:10.3e}")
    print(f"  s_00110   = {s_00110:10.3e}")
    print(f"  s_22000   = {s_22000:10.3e}")
    print(f"  s_11110_1 = {s_11110_1:10.3e}")
    print(f"  s_11110_2 = {s_11110_2:10.3e}")
    print(f"  s_00220   = {s_00220:10.3e}")


def compute_ampl_dep_orbit_tpsa_2(g):
    x_avg       = new.tpsa()
    x_avg_re    = new.tpsa()
    x_avg_im    = new.tpsa()
    x3_avg      = new.tpsa()
    x3_avg_re   = new.tpsa()
    x3_avg_im   = new.tpsa()
    x_y2_avg    = new.tpsa()
    x_y2_avg_re = new.tpsa()
    x_y2_avg_im = new.tpsa()

    Id          = new.ss_vect_tpsa()

    Id.set_identity()

    x_avg = pb(g, Id.x.to_tpsa())
    x_avg = compose_bs(x_avg, A_1)
    x_avg.CtoR(x_avg_re, x_avg_im)
    x_avg_re.print("<x_re>", 1e-11)
    s_11000 = x_avg_re.get([1, 1, 0, 0, 0, 0, 0])
    s_00110 = x_avg_re.get([0, 0, 1, 1, 0, 0, 0])

    f = pb(g, pow(Id.x.to_tpsa(), 3))
    f.print("f", 1e-11)
    assert False
    x3_avg = pb(g, pow(Id.x.to_tpsa(), 3))
    x3_avg.print("<x3>", 1e-11)
    x3_avg = compose_bs(x3_avg, A_1)
    x3_avg.print("<x3>", 1e-11)
    x3_avg.CtoR(x3_avg_re, x3_avg_im)
    x3_avg_re.print("<x^3_re>")
    s_22000 = x3_avg_re.get([2, 2, 0, 0, 0, 0, 0])
    s_11110_1 = x3_avg_re.get([1, 1, 1, 1, 0, 0, 0])

    x_y2_avg = pb(g, Id.x.to_tpsa()*pow(Id.y.to_tpsa(), 2))
    x_y2_avg = compose_bs(x_y2_avg, A_1)
    x_y2_avg.CtoR(x_y2_avg_re, x_y2_avg_im)
    x_y2_avg_re.print("<x*y^2_re>")
    s_11110_2 = x_y2_avg_re.get([1, 1, 1, 1, 0, 0, 0])
    s_00220 = x_y2_avg_re.get([0, 0, 2, 2, 0, 0, 0])

    print(f"\n  s_11000   = {s_11000:10.3e}")
    print(f"  s_00110   = {s_00110:10.3e}")
    print(f"  s_22000   = {s_22000:10.3e}")
    print(f"  s_11110_1 = {s_11110_1:10.3e}")
    print(f"  s_11110_2 = {s_11110_2:10.3e}")
    print(f"  s_00220   = {s_00220:10.3e}")


gtpsa_prop = gtpsa_prop_class()
gtpsa_prop_class.no = 5
gtpsa_prop_class.desc = gtpsa.desc(gtpsa_prop.nv, gtpsa_prop.no)

cod_eps = 1e-15
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

twoJ = compute_twoJ(A_max, beta_inj)
Id_scl = compute_Id_scl(twoJ)

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
    chk_types()
    assert False

if False:
    chk_Poisson_bracket()
    assert False

if False:
    chk_compose()
    assert False

M = compute_map(lat_prop)
print("\nM:", M)

if not False:
    h = compute_h(M)
    h.print("h", 1e-20)

if not False:
    A_0, A_1, R, g, K = compute_map_normal_form(M)
    print("\nA_0:", A_0)
    print("\nA_1:", A_1)
    g.print("g", 1e-10)
    K.print("K", 1e-10)

if False:
    compute_M(K, A_0, A_1, g)
if False:
    compute_R(M, A_0, A_1, g, K)

compute_hh(lat_prop, h)
compute_g(lat_prop, g)

# compute_ampl_dep_orbit_tpsa_2(g)

# compute_ampl_dep_orbit(lat_prop)
