import os
import sys

import copy as _copy

import enum as en
from typing import ClassVar
from dataclasses import dataclass
import numpy as np

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    index_class as ind


class MpoleInd(en.IntEnum):
    quad = 2
    sext = 3

ind = ind.index_class()


@dataclass
class gtpsa_prop:
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
    M.get_mns(1, gtpsa_prop.no-1, M)
    
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

    K.h_DF_to_M(Id, 2, tpsa_prop.no, R)
    R_inv.inv(R)
    R_inv.get_mns(1, tpsa_prop.no-1, R_inv)
    A_0_inv.inv(A_0)
    A_1_inv.inv(A_1)
    g.h_DF_to_M(Id, 3, tpsa_prop.no, A_nl)
    A_nl_inv.inv(A_nl)

    M.compose(M, A_0)
    M.compose(M, A_1)
    M.compose(M, A_nl)
    M.compose(M, R_inv)
    M.compose(A_0_inv, M)
    M.compose(A_1_inv, M)
    M.compose(A_nl_inv, M)
    M.get_mns(1, tpsa_prop.no-1, M)
   
    print("\nM:", M)
    prt_map(M, "M:", eps=1e-10)


def get_h_ijklm(h, i, j, k, l, m):
    return h.get([i, j, k, l, m, 0, 0])


def get_cmplx_h_ijklm(h_re, h_im, i, j, k, l, m):
    return complex(get_h_ijklm(h_re, i, j, k, l, m),
                   get_h_ijklm(h_im, i, j, k, l, m))


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

    h= 0e0
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
                h_abs = b3xL*beta[ind.X, n-1]**(m_x/2e0) \
                    *beta[ind.Y, n-1]**(m_y/2e0)
                phi = n_x*dmu[ind.X, n-1] + n_y*dmu[ind.Y, n-1]
                h += h_abs*np.exp(1j*phi)

    return np.array(h)


def prt_h_cmplx_ijklm(symb, i, j, k, l, m, h):
    ind = np.array([i, j, k, l, m], dtype=int)
    str = ""
    for k in range(len(ind)):
        str += chr(ord('0')+ind[k])
    if h.imag >= 0e0:
        sgn_h_im = '+'
    else:
        sgn_h_im = '-'
    h_abs = abs(h)
    h_arg = np.rad2deg(np.angle(h))
    print(f"  {symb:s}_{str:s} = ({h.real:23.16e} {sgn_h_im:s}"
          f" i{abs(h.imag):22.16e})"
          f"  {h_abs:9.3e} |_ {h_arg:6.1f})")


def prt_h_ijklm(symb, i, j, k, l, m, h_re, h_im):
    h = get_cmplx_h_ijklm(h_re, h_im, i, j, k, l, m)
    prt_h_cmplx_ijklm(symb, i, j, k, l, m, h)


def prt_h(symb, h):
    h_re = new.tpsa()
    h_im = new.tpsa()

    h.CtoR(h_re, h_im)

    prt_h_ijklm(symb, 2, 1, 0, 0, 0, h_re, h_im)
    prt_h_ijklm(symb, 3, 0, 0, 0, 0, h_re, h_im)
    prt_h_ijklm(symb, 1, 0, 1, 1, 0, h_re, h_im)
    prt_h_ijklm(symb, 1, 0, 2, 0, 0, h_re, h_im)
    prt_h_ijklm(symb, 1, 0, 0, 2, 0, h_re, h_im)


def prt_h_num():
    # Can't hash a list (mutable) - but can hash a tuple (immutable).
    h = {}
    h[(2, 1, 0, 0, 0)] = -1e0/8e0 *compute_h_ijklm(lat_prop, 2, 1, 0, 0, 0)
    h[(3, 0, 0, 0, 0)] = -1e0/24e0*compute_h_ijklm(lat_prop, 3, 0, 0, 0, 0)
    h[(1, 0, 1, 1, 0)] =  1e0/4e0 *compute_h_ijklm(lat_prop, 1, 0, 1, 1, 0)
    h[(1, 0, 2, 0, 0)] =  1e0/8e0 *compute_h_ijklm(lat_prop, 1, 0, 2, 0, 0)
    h[(1, 0, 0, 2, 0)] =  1e0/8e0 *compute_h_ijklm(lat_prop, 1, 0, 0, 2, 0)

    for key in h:
        prt_h_cmplx_ijklm("h", key[0], key[1], key[2], key[3], key[4], h[key])

        
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

    g = 0e0
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
            phi = n_x*(dmu[ind.X, n]-np.pi*nu[ind.X]) \
                + n_y*(dmu[ind.Y, n]-np.pi*nu[ind.Y])
            g += 1j*g_abs*np.exp(-1j*(phi)) \
                /np.sin(np.pi*(n_x*nu[ind.X]+n_y*nu[ind.Y]))
    return g


def compute_g(lat_prop, g_tps):
    print("\nLie Generators - TPSA:")
    prt_h('g', g_tps)

    g = {}
    g[(2, 1, 0, 0, 0)] =  1e0/16e0*compute_g_ijklm(lat_prop, 2, 1, 0, 0, 0)
    g[(3, 0, 0, 0, 0)] =  1e0/48e0*compute_g_ijklm(lat_prop, 3, 0, 0, 0, 0)
    g[(1, 0, 1, 1, 0)] = -1e0/8e0 *compute_g_ijklm(lat_prop, 1, 0, 1, 1, 0)
    g[(1, 0, 2, 0, 0)] = -1e0/16e0*compute_g_ijklm(lat_prop, 1, 0, 2, 0, 0)
    g[(1, 0, 0, 2, 0)] = -1e0/16e0*compute_g_ijklm(lat_prop, 1, 0, 0, 2, 0)

    print("\nLie Generators - Numerical:")
    for key in g:
        prt_h_cmplx_ijklm("g", key[0], key[1], key[2], key[3], key[4], g[key])


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


gtpsa_prop.no = 5
gtpsa_prop.desc = gtpsa.desc(gtpsa_prop.nv, gtpsa_prop.no)

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
    A     = new.ss_vect_tpsa()
    A_inv = new.ss_vect_tpsa()
    R_inv = new.ss_vect_tpsa()
    M_Fl  = new.ss_vect_tpsa()

    A_0, A_1, R, g, K = compute_map_normal_form(M)

    A.compose(A_0, A_1)
    A_inv.inv(A)
    R_inv.inv(R)
    M_Fl.compose(M, A)
    M_Fl.compose(A_inv, M_Fl)

    h = compute_h(M_Fl)
    print("\nDriving Terms - TPSA:")
    prt_h('h', h)

if False:
    Id      = new.ss_vect_tpsa()
    g_fp    = new.tpsa()
    g_fp_re = new.tpsa()
    g_fp_im = new.tpsa()

    Id.set_identity()
    op = Id - R_inv
    op.pinv(op, [1, 1, 1, 1, 0, 0, 0])
    g_fp = compose_bs(h, op)
    g_fp.get_mns_1(1, 3)
    print("\nLie Generators - TPSA:")
    prt_h('g', g_fp)

if False:
    compute_M(K, A_0, A_1, g)
if False:
    compute_R(M, A_0, A_1, g, K)

print("\nDriving Terms - Numerical:")
prt_h_num()

if False:
    g_re = new.tpsa()
    g_im = new.tpsa()

    g = g.get_mns_1(1, 3)
    g.CtoR(g_re, g_im)
    g_re.print("g_re")
    g_im.print("g_im")

R_inv = new.ss_vect_tpsa()
R_inv.inv(R)

compute_g(lat_prop, g)

# compute_ampl_dep_orbit_tpsa_2(g)

# compute_ampl_dep_orbit(lat_prop)
