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


def trunc(a, no_cut):
    no = a.get_description().maximum_orders(0, 0)
    a.get_description().truncate(no_cut)
    a = 1e0*a
    a.get_description().truncate(no)
    return a


class nonlinear_beam_dynamics:
    # Private.
    def __init__(self, gtpsa_prop, lat_prop):
        # GTPSA parameters.
        self._gtpsa_prop = gtpsa_prop

        # Lattice properties.
        self._lat_prop   = lat_prop

        # Linear optics.
        self._eta = np.array(
            lat_prop._Twiss.dispersion.sel(phase_coordinate="x").values)
        self._alpha = np.array(
            [lat_prop._Twiss.twiss.sel(plane="x", par="alpha").values,
             lat_prop._Twiss.twiss.sel(plane="y", par="alpha").values]
        )
        self._beta = np.array(
            [lat_prop._Twiss.twiss.sel(plane="x", par="beta").values,
             lat_prop._Twiss.twiss.sel(plane="y", par="beta").values]
        )
        self._dnu = np.array(
            [lat_prop._Twiss.twiss.sel(plane="x", par="nu").values,
             lat_prop._Twiss.twiss.sel(plane="y", par="nu").values])
        self._dmu = 2e0*np.pi*self._dnu
        self._nu = np.array(
            [lat_prop._Twiss.twiss.sel(plane="x", par="nu").values[-1],
             lat_prop._Twiss.twiss.sel(plane="y", par="nu").values[-1]])

        # Identity map.
        self._Id         = new.ss_vect_tpsa()

        # Poincaré map.
        self._M          = np.nan

        # Map normal form.
        self._A_0        = new.ss_vect_tpsa()
        self._A_0_inv    = new.ss_vect_tpsa()
        self._A_1        = new.ss_vect_tpsa()
        self._A_1_inv    = new.ss_vect_tpsa()
        self._A          = new.ss_vect_tpsa()
        self._A_inv      = new.ss_vect_tpsa()
        self._A_nl       = new.ss_vect_tpsa()
        self._A_nl_inv   = new.ss_vect_tpsa()
        self._R          = new.ss_vect_tpsa()
        self._R_inv      = new.ss_vect_tpsa()
        self._g          = new.tpsa()
        self._g_re       = new.tpsa()
        self._g_im       = new.tpsa()
        self._K          = new.tpsa()
        self._K_re       = new.tpsa()
        self._K_im       = new.tpsa()

        # Map in Floquet space.
        self._M_Fl       = np.nan

        self._Id.set_identity()

    # Public.

    def compute_map(self):
        self._M = lo.compute_map(
            self._lat_prop._lattice, self._lat_prop._model_state,
            desc=self._lat_prop._desc,
            tpsa_order=self._gtpsa_prop.no)

    def prt_map(self, str, *, eps: float=1e-30):
        print(str)
        self._map.x.print("x", eps)
        self._map.px.print("p_x", eps)
        self._map.y.print("y", eps)
        self._map.py.print("p_y", eps)
        self._map.delta.print("delta", eps)
        self._map.ct.print("ct", eps)

    def get_f_I(self, f, I):
        return f.get([I[0], I[1], I[2], I[3], I[4]])

    def get_cmplx_f_I(self, f_re, f_im, I):
        return complex(self.get_f_I(f_re, I),
                       self.get_f_I(f_im, I))

    def compute_h_tpsa(self):
        h_tpsa = new.tpsa()
        self.compute_M_Fl()
        self._M_Fl.M_to_h_DF(h_tpsa)
        return h_tpsa

    def compute_map_normal_form(self):
        self._M.Map_Norm(self._A_0, self._A_1, self._R, self._g, self._K)
        self._A.compose(self._A_0, self._A_1)
        self._A_inv.inv(self._A)
        self._R_inv.inv(self._R)
        self._g.CtoR(self._g_re, self._g_im)
        self._K.CtoR(self._K_re, self._K_im)

    def compute_M_Fl(self):
        self._M_Fl = new.ss_vect_tpsa()

        self._M_Fl.compose(self._M, self._A)
        self._M_Fl.compose(self._A_inv, self._M_Fl)
        

    def compute_M_from_MNF(self, K, A_0, A_1, g):
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

    def compute_R_from_MNF(self, map, A_0, A_1, g, K):
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

    def compute_h_I(self, I):
        m_x = I[0] + I[1]
        m_y = I[2] + I[3]
        n_x = I[0] - I[1]
        n_y = I[2] - I[3]
        h_I = 0e0
        for n, ele in enumerate(self._lat_prop._lattice):
            if type(ele) == ts.Sextupole:
                L = ele.get_length()
                b3xL = ele.get_multipoles(). \
                    get_multipole(MpoleInd.sext).real*L
                if b3xL != 0e0:
                    h_abs = b3xL*self._beta[ind.X, n-1]**(m_x/2e0) \
                        *self._beta[ind.Y, n-1]**(m_y/2e0)
                    phi = n_x*self._dmu[ind.X, n-1] + n_y*self._dmu[ind.Y, n-1]
                    h_I += h_abs*np.exp(1j*phi)

        return h_I

    def compute_h_fp(self):
        # Can't hash a list (mutable) - but can hash a tuple (immutable).
        h_fp = {}
        h_fp[(2, 1, 0, 0, 0)] = -1e0/8e0 *self.compute_h_I([2, 1, 0, 0, 0])
        h_fp[(3, 0, 0, 0, 0)] = -1e0/24e0*self.compute_h_I([3, 0, 0, 0, 0])
        h_fp[(1, 0, 1, 1, 0)] =  1e0/4e0 *self.compute_h_I([1, 0, 1, 1, 0])
        h_fp[(1, 0, 2, 0, 0)] =  1e0/8e0 *self.compute_h_I([1, 0, 2, 0, 0])
        h_fp[(1, 0, 0, 2, 0)] =  1e0/8e0 *self.compute_h_I([1, 0, 0, 2, 0])
        return h_fp

    def compute_g_I(self, j, I):
        m_x = I[0] + I[1]
        m_y = I[2] + I[3]
        n_x = I[0] - I[1]
        n_y = I[2] - I[3]
        g_fp_I = 0e0
        for k, ele in enumerate(self._lat_prop._lattice):
            if type(ele) == ts.Sextupole:
                L = ele.get_length()
                b3xL = \
                    ele.get_multipoles(). \
                    get_multipole(MpoleInd.sext).real*L
                if (b3xL != 0e0):
                    g_abs = b3xL*self._beta[ind.X, k-1]**(m_x/2e0) \
                        *self._beta[ind.Y, k-1]**(m_y/2e0)
                    phi = n_x*(self._dmu[ind.X, k-1]-self._dmu[ind.X, j-1]) \
                        + n_y*(self._dmu[ind.Y, k-1]-self._dmu[ind.Y, j-1])
                    dnu = np.pi*(n_x*self._nu[ind.X]+n_y*self._nu[ind.Y])
                    g_fp_I -= 1j*g_abs*np.exp(1j*(phi-dnu))/np.sin(dnu)
        return g_fp_I

    def compute_g_fp(self, j):
        g_fp = {}
        g_fp[(2, 1, 0, 0, 0)] =  1e0/16e0*self.compute_g_I(j, [2, 1, 0, 0, 0])
        g_fp[(3, 0, 0, 0, 0)] =  1e0/48e0*self.compute_g_I(j, [3, 0, 0, 0, 0])
        g_fp[(1, 0, 1, 1, 0)] = -1e0/8e0 *self.compute_g_I(j, [1, 0, 1, 1, 0])
        g_fp[(1, 0, 2, 0, 0)] = -1e0/16e0*self.compute_g_I(j, [1, 0, 2, 0, 0])
        g_fp[(1, 0, 0, 2, 0)] = -1e0/16e0*self.compute_g_I(j, [1, 0, 0, 2, 0])
        return g_fp

    def prt_f_cmplx_I(self, symb, I, f):
        str = ""
        for k in I:
            str += chr(ord('0')+k)
        if f.imag >= 0e0:
            sgn_f_im = '+'
        else:
            sgn_f_im = '-'
        f_abs = abs(f)
        f_arg = np.rad2deg(np.angle(f))
        print(f"  {symb:s}_{str:s} = ({f.real:23.16e} {sgn_f_im:s}"
              f" i{abs(f.imag):22.16e})"
              f"  {f_abs:9.3e} |_ {f_arg:6.1f})")

    def prt_f_I(self, symb, I, f_re, f_im):
        f = self.get_cmplx_f_I(f_re, f_im, I)
        self.prt_f_cmplx_I(symb, I, f)

    def prt_f_tpsa(self, title, symb, f_tpsa):
        f_re = new.tpsa()
        f_im = new.tpsa()

        f_tpsa.CtoR(f_re, f_im)

        print(title)
        self.prt_f_I(symb, [2, 1, 0, 0, 0], f_re, f_im)
        self.prt_f_I(symb, [3, 0, 0, 0, 0], f_re, f_im)
        self.prt_f_I(symb, [1, 0, 1, 1, 0], f_re, f_im)
        self.prt_f_I(symb, [1, 0, 2, 0, 0], f_re, f_im)
        self.prt_f_I(symb, [1, 0, 0, 2, 0], f_re, f_im)

    def prt_f_fp(self, title, symb, f_fp):
        print(title)
        for key, value in f_fp.items():
            self.prt_f_cmplx_I(symb, key, value)

    def compute_dq_dJ_tpsa(self):
        x = new.tpsa()
        x_re = new.tpsa()
        x_im = new.tpsa()

        x = self._g.poisbra(self._Id.x.to_tpsa(), 4)
        x.CtoR(x_re, x_im)
        x = complex(x_re.get([1, 1, 0, 0, 0]), x_im.get([1, 1, 0, 0, 0]))
        print(f"\n  x_11000 = {x:21.3e}")
        x = complex(x_re.get([0, 0, 1, 1, 0]), x_im.get([0, 0, 1, 1, 0]))
        print(f"  x_00110 = {x:21.3e}")

        x = self._g.poisbra(self._Id.x.to_tpsa()**3, 4)
        x.CtoR(x_re, x_im)
        x = complex(x_re.get([2, 2, 0, 0, 0]), x_im.get([2, 2, 0, 0, 0]))
        print(f"  x_22000 = {x:21.3e}")

        x = self._g.poisbra(self._Id.x.to_tpsa()*self._Id.y.to_tpsa()**2, 4)
        x.CtoR(x_re, x_im)
        x = complex(x_re.get([1, 1, 1, 1, 0]), x_im.get([1, 1, 1, 1, 0]))
        print(f"  x_11110 = {x:21.3e}")
        x = complex(x_re.get([0, 0, 2, 2, 0]), x_im.get([0, 0, 2, 2, 0]))
        print(f"  x_00220 = {x:21.3e}")

    def compute_dq_dJ_fp_I(self, j, I):
        m_x = I[0] + I[1]
        m_y = I[2] + I[3]
        n_x = I[0] - I[1]
        n_y = I[2] - I[3]
        h_I = 0e0
        dq = 0e0
        for k, ele in enumerate(self._lat_prop._lattice):
            if type(ele) == ts.Sextupole:
                L = ele.get_length()
                b3xL = \
                    ele.get_multipoles(). \
                    get_multipole(MpoleInd.sext).real*L
                if (b3xL != 0e0):
                    dq_abs = b3xL*self._beta[ind.X, k-1]**(m_x/2e0) \
                        *self._beta[ind.Y, k-1]**(m_y/2e0)
                    phi = abs(n_x*(self._dmu[ind.X, k-1]-self._dmu[ind.X, j-1])
                              + n_y*(self._dmu[ind.Y, k-1]
                                     -self._dmu[ind.Y, j-1]))
                    dnu = np.pi*(n_x*self._nu[ind.X]+n_y*self._nu[ind.Y])
                    dq += dq_abs*np.cos((phi-dnu))/np.sin(dnu)
        return dq
        
    def compute_dq_dJ_fp(self, j):
        Id = new.ss_vect_tpsa()
        
        Id.set_identity()

        dq = -np.sqrt(self._beta[ind.X, j-1])*(
            nlbd.compute_dq_dJ_fp_I(1, [2, 1, 0, 0, 0])
            *Id.x.to_tpsa()*Id.px.to_tpsa()
            -nlbd.compute_dq_dJ_fp_I(1, [1, 0, 1, 1, 0])
            *Id.y.to_tpsa()*Id.py.to_tpsa()) \
            /8e0

        dq -= 3e0*np.sqrt(self._beta[ind.X, j-1])*(
            nlbd.compute_dq_dJ_fp_I(1, [3, 0, 0, 0, 0])
            +3e0*nlbd.compute_dq_dJ_fp_I(1, [2, 1, 0, 0, 0])) \
            *Id.x.to_tpsa()**2*Id.px.to_tpsa()**2 \
            /64e0

        dq += np.sqrt(self._beta[ind.X, j-1])*(
            4e0*(nlbd.compute_dq_dJ_fp_I(1, [1, 0, 2, 0, 0])
                 -nlbd.compute_dq_dJ_fp_I(1, [2, 1, 0, 0, 0])
                 -nlbd.compute_dq_dJ_fp_I(1, [1, 0, 0, 2, 0]))
            *Id.x.to_tpsa()*Id.px.to_tpsa()*Id.y.to_tpsa()*Id.py.to_tpsa()

            +(1e0*nlbd.compute_dq_dJ_fp_I(1, [1, 0, 2, 0, 0])
              +4e0*nlbd.compute_dq_dJ_fp_I(1, [1, 0, 1, 1, 0])
              +1e0*nlbd.compute_dq_dJ_fp_I(1, [1, 0, 0, 2, 0]))
            *Id.y.to_tpsa()**2*Id.py.to_tpsa()**2
        )/64e0

        return dq

    def compute_dnu_dJ(self):
        a = np.zeros(3, dtype=complex)
        a_k = np.zeros(3, dtype=complex)
        for j, ele in enumerate(self._lat_prop._lattice):
            if type(ele) == ts.Sextupole:
                L = ele.get_length()
                b3xL = ele.get_multipoles(). \
                    get_multipole(MpoleInd.sext).real*L
                if b3xL != 0e0:
                    g_k = self.compute_g_fp(j)
                    a[0] -= b3xL*self._beta[ind.X, j-1] \
                        *(g_k[2, 1, 0, 0, 0]+g_k[3, 0, 0, 0, 0])
                    a[1] -= b3xL*(
                        self._beta[ind.X, j-1]*g_k[1, 0, 1, 1, 0]
                        -self._beta[ind.Y, j-1]
                        *(g_k[1, 0, 2, 0, 0]+g_k[1, 0, 0, 2, 0]))
                    a[2] += b3xL*(
                        self._beta[ind.Y, j-1]*(
                            g_k[1, 0, 1, 1, 0]+g_k[1, 0, 2, 0, 0]
                            +g_k[1, 0, 0, 2, 0]))

                    a_k[0] += b3xL*(self._beta[ind.Y, j-1]*g_k[1, 0, 1, 1, 0])
                    a_k[1] += b3xL*(self._beta[ind.Y, j-1]*g_k[1, 0, 2, 0, 0])
                    a_k[2] += b3xL*(self._beta[ind.Y, j-1]*g_k[1, 0, 0, 2, 0])

        return a/(2e0*np.pi), a_k/(2e0*np.pi)
        

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


def compose_bs(h, map):
    t_map = new.ss_vect_tpsa()
    t_map.set_zero()
    t_map.x = h
    t_map.compose(t_map, map)
    return t_map.x.to_tpsa()


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

    ps_dim = 4

    f = get_mn([3, 0, 0, 0, 0])
    f.poisbra(get_mn([0, 3, 0, 0, 0]), ps_dim).print("[h^3+_x, h^3-_x]")

    f = get_mn([2, 1, 0, 0, 0])
    f.poisbra(get_mn([1, 2, 0, 0, 0]),
              ps_dim).print("[h^2+_x*h^-_, h^+_x*h^2-_x]")

    assert False

    f = Id.x.to_tpsa()**2*Id.px.to_tpsa()
    g = Id.px.to_tpsa()*Id.y.to_tpsa()*Id.py.to_tpsa()
    f.poisbra(g, ps_dim).print("[f, g]")

    f = Id.x.to_tpsa()*Id.y.to_tpsa()*Id.py.to_tpsa()
    g = Id.px.to_tpsa()*Id.y.to_tpsa()*Id.py.to_tpsa()
    f.poisbra(g, ps_dim).print("[f, g]")

    f = Id.x.to_tpsa()*Id.py.to_tpsa()**2
    g = Id.px.to_tpsa()*Id.y.to_tpsa()**2
    f.poisbra(g, ps_dim).print("[f, g]")

    assert False

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
    pow(Id.x.to_tpsa(), 3).poisbra(pow(Id.x.to_tpsa(), ps_dun)
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

    f.set([1, 0, 0, 0, 0], 0e0, 1e0)
    print("\n", type(f), type(f))
    g = compose_bs(f, map_1)
    print("\n", type(f), type(g))

    map_2.exppb(map_1, map_2)
    print("\nmap_2", map_2)


def pb(f, g):
    n_dof = 2
    a = new.tpsa()
    for k in range(n_dof):
        a += f.deriv(2*k+1)*g.deriv(2*k+2) - f.deriv(2*k+2)*g.deriv(2*k+1)
    return a


def get_mn(I):
    Id = new.ss_vect_tpsa()
    Id.set_identity()
    mn = 1e0
    for k, i in enumerate(I):
        mn *= Id.iloc[k].to_tpsa()**i
    return mn


def chk_terms(nlbd):
    ps_dim = 4;

    x_re = new.tpsa()
    x_im = new.tpsa()

    print("\n<x>:")
    x = nlbd._g.poisbra(get_mn([1, 0, 0, 0, 0]), ps_dim)
    x.CtoR(x_re, x_im)
    x_I  = nlbd.get_f_I(x_re, [1, 1, 0, 0, 0])
    scl = np.sqrt(nlbd._beta[ind.X, 0])/8e0
    dq = -scl*nlbd.compute_dq_dJ_fp_I(1, [2, 1, 0, 0, 0])
    print(f"  s_11000 = {x_I:12.5e}")
    print(f"            {dq:12.5e}")

    x_I  = nlbd.get_f_I(x_re, [0, 0, 1, 1, 0])
    dq = scl*nlbd.compute_dq_dJ_fp_I(1, [1, 0, 1, 1, 0])
    print(f"\n  s_00110 = {x_I:12.5e}")
    print(f"            {dq:12.5e}")

    print("\n<x^3>:")
    x = nlbd._g.poisbra(get_mn([3, 0, 0, 0, 0]), ps_dim)
    x.CtoR(x_re, x_im)
    x_I  = nlbd.get_f_I(x_re, [2, 2, 0, 0, 0])

    # A x4 to correct for cosine & sine terms - for both f & g -> [f, g].
    scl = 4e0/8e0
    g = nlbd._g_im.get([3, 0, 0, 0, 0])*get_mn([3, 0, 0, 0, 0])
    x_1_tps = scl*g.poisbra(get_mn([0, 3, 0, 0, 0]), ps_dim)
    x_1_tps = x_1_tps.get([2, 2, 0, 0, 0])
    g = nlbd._g_im.get([2, 1, 0, 0, 0])*get_mn([2, 1, 0, 0, 0])
    x_2_tps = 3e0*scl*g.poisbra(get_mn([1, 2, 0, 0, 0]), ps_dim)
    x_2_tps = x_2_tps.get([2, 2, 0, 0, 0])
    x_tps = x_1_tps + x_2_tps

    scl = 3e0*np.sqrt(nlbd._beta[ind.X, 0])/64e0
    dq_1_fp = -scl*nlbd.compute_dq_dJ_fp_I(1, [3, 0, 0, 0, 0])
    dq_2_fp = -scl*3e0*nlbd.compute_dq_dJ_fp_I(1, [2, 1, 0, 0, 0])
    dq_fp = dq_1_fp + dq_2_fp

    print(f"  s_22000 = {x_I:12.5e}")
    print(f"            {x_tps:12.5e} = {x_1_tps:+12.5e} {x_2_tps:+12.5e}")
    print(f"            {dq_fp:12.5e} = {dq_1_fp:+12.5e} {dq_2_fp:+12.5e}")

    print("\n<x*y^2>:")
    x = nlbd._g.poisbra(get_mn([1, 0, 2, 0, 0]), ps_dim)
    x.CtoR(x_re, x_im)
    x_I  = nlbd.get_f_I(x_re, [1, 1, 1, 1, 0])

    # A x4 to correct for cosine & sine terms - for both f & g -> [f, g].
    scl = 4e0/8e0
    g = nlbd._g_im.get([1, 0, 2, 0, 0])*get_mn([1, 0, 2, 0, 0])
    x_1_tps = scl*g.poisbra(get_mn([0, 1, 0, 2, 0]), ps_dim)
    x_1_tps = x_1_tps.get([1, 1, 1, 1, 0])
    g = nlbd._g_im.get([2, 1, 0, 0, 0])*get_mn([2, 1, 0, 0, 0])
    x_2_tps = 2e0*scl*g.poisbra(get_mn([0, 1, 1, 1, 0]), ps_dim)
    x_2_tps = x_2_tps.get([1, 1, 1, 1, 0])
    g = nlbd._g_im.get([1, 0, 0, 2, 0])*get_mn([1, 0, 0, 2, 0])
    x_3_tps = scl*g.poisbra(get_mn([0, 1, 2, 0, 0]), ps_dim)
    x_3_tps = x_3_tps.get([1, 1, 1, 1, 0])
    x_tps = x_1_tps + x_2_tps + x_3_tps

    scl = np.sqrt(nlbd._beta[ind.X, 0])/16e0
    dq_1_fp = scl*nlbd.compute_dq_dJ_fp_I(1, [1, 0, 2, 0, 0])
    dq_2_fp = -scl*nlbd.compute_dq_dJ_fp_I(1, [2, 1, 0, 0, 0])
    dq_3_fp = -scl*nlbd.compute_dq_dJ_fp_I(1, [1, 0, 0, 2, 0])
    dq_fp = dq_1_fp + dq_2_fp + dq_3_fp

    print(f"  s_11110 = {x_I:12.5e}")
    print(f"            {x_tps:12.5e} = {x_1_tps:+12.5e} {x_2_tps:+12.5e}"
          f" {x_3_tps:+12.5e}")
    print(f"            {dq_fp:12.5e} = {dq_1_fp:+12.5e} {dq_2_fp:+12.5e}"
          f" {dq_3_fp:+12.5e}")

    x_I  = nlbd.get_f_I(x_re, [0, 0, 2, 2, 0])

    # A x4 to correct for cosine & sine terms - for both f & g -> [f, g].
    scl = 4e0/8e0
    g = nlbd._g_im.get([1, 0, 2, 0, 0])*get_mn([1, 0, 2, 0, 0])
    x_1_tps = scl*g.poisbra(get_mn([0, 1, 0, 2, 0]), ps_dim)
    x_1_tps = x_1_tps.get([0, 0, 2, 2, 0])
    g = nlbd._g_im.get([1, 0, 1, 1, 0])*get_mn([1, 0, 1, 1, 0])
    x_2_tps = 2e0*scl*g.poisbra(get_mn([0, 1, 1, 1, 0]), ps_dim)
    x_2_tps = x_2_tps.get([0, 0, 2, 2, 0])
    g = nlbd._g_im.get([1, 0, 0, 2, 0])*get_mn([1, 0, 0, 2, 0])
    x_3_tps = scl*g.poisbra(get_mn([0, 1, 2, 0, 0]), ps_dim)
    x_3_tps = x_3_tps.get([0, 0, 2, 2, 0])
    x_tps = x_1_tps + x_2_tps + x_3_tps

    scl = np.sqrt(nlbd._beta[ind.X, 0])/64e0
    dq_1_fp = scl*nlbd.compute_dq_dJ_fp_I(1, [1, 0, 2, 0, 0])
    dq_2_fp = scl*4e0*nlbd.compute_dq_dJ_fp_I(1, [1, 0, 1, 1, 0])
    dq_3_fp = scl*nlbd.compute_dq_dJ_fp_I(1, [1, 0, 0, 2, 0])
    dq_fp = dq_1_fp + dq_2_fp + dq_3_fp
    print(f"\n  s_00220 = {x_I:12.5e}")
    print(f"            {x_tps:12.5e} = {x_1_tps:+12.5e} {x_2_tps:+12.5e}"
          f" {x_3_tps:+12.5e}")
    print(f"            {dq_fp:12.5e} = {dq_1_fp:+12.5e} {dq_2_fp:+12.5e}"
          f" {dq_3_fp:+12.5e}")


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
    s_11000 = x_avg_re.get([1, 1, 0, 0, 0])
    s_00110 = x_avg_re.get([0, 0, 1, 1, 0])

    x3_avg = g.poisbra(pow(Id.x.to_tpsa(), 3), ps_dim)
    x3_avg.print("<x3>", 1e-11)
    x3_avg = compose_bs(x3_avg, A_1)
    x3_avg.print("<x3>", 1e-11)
    x3_avg.CtoR(x3_avg_re, x3_avg_im)
    x3_avg_re.print("<x^3_re>")
    s_22000 = x3_avg_re.get([2, 2, 0, 0, 0])
    s_11110_1 = x3_avg_re.get([1, 1, 1, 1, 0])

    x_y2_avg = g.poisbra(Id.x.to_tpsa()*pow(Id.y.to_tpsa(), 2), ps_dim)
    x_y2_avg = compose_bs(x_y2_avg, A_1)
    x_y2_avg.CtoR(x_y2_avg_re, x_y2_avg_im)
    x_y2_avg_re.print("<x*y^2_re>")
    s_11110_2 = x_y2_avg_re.get([1, 1, 1, 1, 0])
    s_00220 = x_y2_avg_re.get([0, 0, 2, 2, 0])

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
    s_11000 = x_avg_re.get([1, 1, 0, 0, 0])
    s_00110 = x_avg_re.get([0, 0, 1, 1, 0])

    f = pb(g, pow(Id.x.to_tpsa(), 3))
    f.print("f", 1e-11)
    assert False
    x3_avg = pb(g, pow(Id.x.to_tpsa(), 3))
    x3_avg.print("<x3>", 1e-11)
    x3_avg = compose_bs(x3_avg, A_1)
    x3_avg.print("<x3>", 1e-11)
    x3_avg.CtoR(x3_avg_re, x3_avg_im)
    x3_avg_re.print("<x^3_re>")
    s_22000 = x3_avg_re.get([2, 2, 0, 0, 0])
    s_11110_1 = x3_avg_re.get([1, 1, 1, 1, 0])

    x_y2_avg = pb(g, Id.x.to_tpsa()*pow(Id.y.to_tpsa(), 2))
    x_y2_avg = compose_bs(x_y2_avg, A_1)
    x_y2_avg.CtoR(x_y2_avg_re, x_y2_avg_im)
    x_y2_avg_re.print("<x*y^2_re>")
    s_11110_2 = x_y2_avg_re.get([1, 1, 1, 1, 0])
    s_00220 = x_y2_avg_re.get([0, 0, 2, 2, 0])

    print(f"\n  s_11000   = {s_11000:10.3e}")
    print(f"  s_00110   = {s_00110:10.3e}")
    print(f"  s_22000   = {s_22000:10.3e}")
    print(f"  s_11110_1 = {s_11110_1:10.3e}")
    print(f"  s_11110_2 = {s_11110_2:10.3e}")
    print(f"  s_00220   = {s_00220:10.3e}")


gtpsa_prop.no = 4
gtpsa_prop.desc = gtpsa.desc(gtpsa_prop.nv, gtpsa_prop.no)

cod_eps = 1e-15
E_0     = 3.0e9

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

nlbd = nonlinear_beam_dynamics(gtpsa_prop, lat_prop)

nlbd.compute_map()
print("\nM:", nlbd._M)

nlbd.compute_map_normal_form()

nlbd._K_re.print("K_re")

if False:
    chk_terms(nlbd)
    assert False

if False:
    compute_M_from_MNF(K, A_0, A_1, g)
if False:
    compute_R_from_MNF(M, A_0, A_1, g, K)

nlbd.compute_dq_dJ_tpsa()

dq = nlbd.compute_dq_dJ_fp(1)
dq.print("dq")
assert False

h_tpsa = nlbd.compute_h_tpsa()
h_tpsa = h_tpsa.get_mns_1(1, 3)
h_tpsa_re = new.tpsa()
h_tpsa_im = new.tpsa()
h_tpsa.CtoR(h_tpsa_re, h_tpsa_im)
h_tpsa_re.print("h_tpsa_re")
nlbd.prt_f_tpsa("\nDriving Terms - TPSA:", 'h', h_tpsa)

h_fp = nlbd.compute_h_fp()
nlbd.prt_f_fp("\nDriving Terms - FP:", 'h', h_fp)

nlbd.prt_f_tpsa("\nLie Generators - TPSA:", 'g', nlbd._g)
g_fp = nlbd.compute_g_fp(1)
nlbd.prt_f_fp("\nLie Generators - FP:", 'g', g_fp)

K = [K_re.get([2, 2, 0, 0, 0]), K_re.get([1, 1, 1, 1, 0]),
     K_re.get([0, 0, 2, 2, 0])]
print(f"\nK:\n [{K[0]:10.3e}, {K[1]:10.3e}, {K[2]:10.3e}]")
assert False
K, K_k = nlbd.compute_dnu_dJ()
print(f"  [{abs(K[0]):10.3e} |_ {np.rad2deg(np.angle(K[0])):6.1f}"
      f", {abs(K[1]):10.3e} |_ {np.rad2deg(np.angle(K[1])):6.1f}"
      f", {abs(K[2]):10.3e}] |_ {np.rad2deg(np.angle(K[2])):6.1f}")
print(f"  [{abs(K_k[0]):10.3e} |_ {np.rad2deg(np.angle(K_k[0])):6.1f}"
      f", {abs(K_k[1]):10.3e} |_ {np.rad2deg(np.angle(K_k[1])):6.1f}"
      f", {abs(K_k[2]):10.3e}] |_ {np.rad2deg(np.angle(K_k[2])):6.1f}")

# compute_ampl_dep_orbit_tpsa_2(g)

# compute_ampl_dep_orbit(lat_prop)
