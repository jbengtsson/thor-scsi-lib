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

from thor_scsi.utils.output import mat2txt


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
    def ctpsa():
        return gtpsa.ctpsa(gtpsa_prop.desc, gtpsa_prop.no)
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

        # Constants
        self._I_21000 = [2, 1, 0, 0, 0]
        self._I_30000 = [3, 0, 0, 0, 0]
        self._I_10110 = [1, 0, 1, 1, 0]
        self._I_10200 = [1, 0, 2, 0, 0]
        self._I_10020 = [1, 0, 0, 2, 0] 

        self._I_03000 = [0, 3, 0, 0, 0]

        self._I_11000 = [1, 1, 0, 0, 0]
        self._I_00110 = [0, 0, 1, 1, 0]
        self._I_22000 = [2, 2, 0, 0, 0]
        self._I_11110 = [1, 1, 1, 1, 0]
        self._I_00220 = [0, 0, 2, 2, 0]

        # Can't hash a list (mutable) - but can hash a tuple (immutable).
        self._I_21000_tuple = (2, 1, 0, 0, 0)
        self._I_30000_tuple = (3, 0, 0, 0, 0)
        self._I_10110_tuple = (1, 0, 1, 1, 0)
        self._I_10200_tuple = (1, 0, 2, 0, 0)
        self._I_10020_tuple = (1, 0, 0, 2, 0) 

        # TPSA parameters.
        self._gtpsa_prop = gtpsa_prop

        # Lattice properties.
        self._lat_prop   = lat_prop

        # Linear optics.
        self._s = np.array(lat_prop._Twiss.s)
        self._eta = np.array(
            [lat_prop._Twiss.dispersion.sel(phase_coordinate="x").values,
            lat_prop._Twiss.dispersion.sel(phase_coordinate="px").values]
        )
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

    # def compute_map(self):
    #     self._M = lo.compute_map(
    #         self._lat_prop._lattice, self._lat_prop._model_state,
    #         desc=self._lat_prop._desc,
    #         tpsa_order=self._gtpsa_prop.no)

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
                b3xL = ele.get_multipoles().get_multipole(MpoleInd.sext).real*L
                if b3xL != 0e0:
                    h_abs = b3xL*self._beta[ind.X, n-1]**(m_x/2e0) \
                        *self._beta[ind.Y, n-1]**(m_y/2e0)
                    phi = n_x*self._dmu[ind.X, n-1] + n_y*self._dmu[ind.Y, n-1]
                    h_I += h_abs*np.exp(1j*phi)

        return h_I

    def compute_h_fp(self):
        # Can't hash a list (mutable) - but can hash a tuple (immutable).
        h_fp = {}
        h_fp[self._I_21000_tuple] = -1e0/8e0 *self.compute_h_I(self._I_21000)
        h_fp[self._I_30000_tuple] = -1e0/24e0*self.compute_h_I(self._I_30000)
        h_fp[self._I_10110_tuple] =  1e0/4e0 *self.compute_h_I(self._I_10110)
        h_fp[self._I_10200_tuple] =  1e0/8e0 *self.compute_h_I(self._I_10200)
        h_fp[self._I_10020_tuple] =  1e0/8e0 *self.compute_h_I(self._I_10020)

        return h_fp

    def compute_g_I(self, j, I):
        m_x = I[0] + I[1]
        m_y = I[2] + I[3]
        n_x = I[0] - I[1]
        n_y = I[2] - I[3]
        g_I = 0e0
        dmu_j = np.zeros(2)
        for k, ele in enumerate(self._lat_prop._lattice):
            if type(ele) == ts.Sextupole:
                L = ele.get_length()
                b3xL = ele.get_multipoles().get_multipole(MpoleInd.sext).real*L
                if (b3xL != 0e0):
                    for l in range(2):
                        # -1 is last element.
                        dmu_j[l] = self._dmu[l, j-1]
                        if j-1 < 0:
                            dmu_j[l] -= 2e0*np.pi*self._nu[l]
                    g_abs = b3xL*self._beta[ind.X, k-1]**(m_x/2e0) \
                        *self._beta[ind.Y, k-1]**(m_y/2e0)
                    phi = n_x*abs(self._dmu[ind.X, k-1]-dmu_j[ind.X]) \
                        + n_y*abs(self._dmu[ind.Y, k-1]-dmu_j[ind.Y])
                    dnu = np.pi*(n_x*self._nu[ind.X]+n_y*self._nu[ind.Y])
                    g_I -= 1j*g_abs*np.exp(1j*(phi-dnu))/np.sin(dnu)
        return g_I

    def compute_g_fp(self, j):
        # Compute the Lie generators.
        g_re = new.tpsa()
        g_im = new.tpsa()
        g = 1e0/16e0*self.compute_g_I(j, self._I_21000)
        g_re.set(self._I_21000, 0e0, g.real)
        g_im.set(self._I_21000, 0e0, g.imag)
        g = 1e0/48e0*self.compute_g_I(j, self._I_30000)
        g_re.set(self._I_30000, 0e0, g.real)
        g_im.set(self._I_30000, 0e0, g.imag)
        g = -1e0/8e0 *self.compute_g_I(j, self._I_10110)
        g_re.set(self._I_10110, 0e0, g.real)
        g_im.set(self._I_10110, 0e0, g.imag)
        g = -1e0/16e0*self.compute_g_I(j, self._I_10200)
        g_re.set(self._I_10200, 0e0, g.real)
        g_im.set(self._I_10200, 0e0, g.imag)
        g = -1e0/16e0*self.compute_g_I(j, self._I_10020)
        g_re.set(self._I_10020, 0e0, g.real)
        g_im.set(self._I_10020, 0e0, g.imag)
        return g_re, g_im

    def prt_ampl_phase(self, str, c, nl="\n"):
        print(f"{str:s}{abs(c):10.3e} \u2220 {np.rad2deg(np.angle(c)):6.1f}",
              end=nl)

    def I2str(self, I):
        str = ""
        for k in I:
            str += chr(ord('0')+k)
        return str

    def prt_f_cmplx_I(self, symb, I, f):
        str = self.I2str(I)
        print(f"  {symb:s}_{str:s} =", end="")
        self.prt_ampl_phase("", f, nl="")

    def prt_f_I(self, symb, I, f_re, f_im):
        f = self.get_cmplx_f_I(f_re, f_im, I)
        self.prt_f_cmplx_I(symb, I, f)

    def prt_f_tpsa(self, title, symb, f_tpsa):
        f_re = new.tpsa()
        f_im = new.tpsa()

        f_tpsa.CtoR(f_re, f_im)

        print(title)
        self.prt_f_I(symb, self._I_21000, f_re, f_im)
        self.prt_f_I(symb, self._I_30000, f_re, f_im)
        self.prt_f_I(symb, self._I_10110, f_re, f_im)
        self.prt_f_I(symb, self._I_10200, f_re, f_im)
        self.prt_f_I(symb, self._I_10020, f_re, f_im)

    def prt_f_fp(self, title, symb, f_fp):
        print(title)
        for key, value in f_fp.items():
            self.prt_f_cmplx_I(symb, key, value)
            print()

    def propagate_g(self, j):
        # Propagate g.
        Id = new.ss_vect_tpsa()
        R = new.ss_vect_tpsa()
        Id.set_identity()
        R.set_identity()
        c = np.zeros(2, dtype=float)
        s = np.zeros(2, dtype=float)
        dmu_j = np.zeros(2, dtype=float)
        for k in range(2):
            # -1 is last element.
            dmu_j[k] = self._dmu[k, j-1]
            if j-1 < 0:
                dmu_j[k] -= 2e0*np.pi*self._nu[k]
        for k, ele in enumerate(c):
            c[k] = np.cos(dmu_j[k])
            s[k] = np.sin(dmu_j[k])
        for k in range(2):
            R.iloc[2*k] = c[k]*Id.iloc[2*k].to_tpsa() \
                - s[k]*Id.iloc[2*k+1].to_tpsa()
            R.iloc[2*k+1] = s[k]*Id.iloc[2*k].to_tpsa() \
                + c[k]*Id.iloc[2*k+1].to_tpsa()
        self._g = compose_bs(self._g, R)

    def compute_dx_dJ_PB(self, j):
        ps_dim = 4

        dx = new.tpsa()
        dx_re = new.tpsa()
        dx_im = new.tpsa()
        x1_re = new.tpsa()
        x3_re = new.tpsa()
        xy2_re = new.tpsa()

        # g = self.propagate_g(j)

        # Imaginary part is zero.

        self._A_1 = compute_A_CS(2, self._A_1)[1]

        # <x>.
        dx = self._g.poisbra(
            compose_bs(self._Id.x.to_tpsa(), self._A_1), ps_dim)
        dx.CtoR(dx_re, dx_im)
        x1_re.set(self._I_11000, 0e0, self._beta[ind.X, j]
                  *dx_re.get(self._I_11000))
        x1_re.set(self._I_00110, 0e0, self._beta[ind.Y, j]
                  *dx_re.get(self._I_00110))

        # <x^3>.
        dx = (self._g.get(self._I_03000)*get_mn(self._I_03000)).poisbra(
            compose_bs(self._Id.x.to_tpsa()**3, self._A_1), ps_dim)
        dx.CtoR(dx_re, dx_im)

        # Divide by J_x.
        x3_re.set(self._I_11000, 0e0, dx_re.get(self._I_22000))
        x3_re.set(self._I_00110, 0e0, dx_re.get(self._I_11110))

        # <x*y^2>.
        dx = self._g.poisbra(
            compose_bs(self._Id.x.to_tpsa(), self._A_1)
            *compose_bs(self._Id.y.to_tpsa()**2, self._A_1), ps_dim)
        dx.CtoR(dx_re, dx_im)

        # Divide by J_y.
        xy2_re.set(self._I_11000, 0e0, -dx_re.get(self._I_11110))
        xy2_re.set(self._I_00110, 0e0, -dx_re.get(self._I_00220))

        return x1_re, x3_re, xy2_re

    def compute_dx_dJ_fp_I(self, j, I):
        dmu_j = np.zeros(2)
        dmu_k = np.zeros(2)

        for l in range(2):
            # -1 is last element.
            dmu_j[l] = self._dmu[l, j-1]
            if j-1 < 0:
                dmu_j[l] -= 2e0*np.pi*self._nu[l]

        m_x = I[0] + I[1]
        m_y = I[2] + I[3]
        n_x = I[0] - I[1]
        n_y = I[2] - I[3]
        h_I = 0e0
        dx_dJ = 0e0
        for k, ele in enumerate(self._lat_prop._lattice):
            if type(ele) == ts.Sextupole:
                L = ele.get_length()
                b3xL = ele.get_multipoles().get_multipole(MpoleInd.sext).real*L
                if (b3xL != 0e0):
                    for l in range(2):
                        # -1 is last element.
                        dmu_k[l] = self._dmu[l, k-1]
                        if k-1 < 0:
                            dmu_k[l] -= 2e0*np.pi*self._nu[l]
                    dx_dJ_abs = b3xL*self._beta[ind.X, k-1]**(m_x/2e0) \
                        *self._beta[ind.Y, k-1]**(m_y/2e0)
                    phi = n_x*abs(dmu_k[ind.X]-dmu_j[ind.X]) \
                        + n_y*abs((dmu_k[ind.Y]-dmu_j[ind.Y]))
                    dnu = np.pi*(n_x*self._nu[ind.X]+n_y*self._nu[ind.Y])
                    dx_dJ += dx_dJ_abs*np.cos(phi-dnu)/np.sin(dnu)

        return dx_dJ

    def compute_dx_dJ_k_fp(self, j):
        dx_dJ_k = np.zeros(9, dtype=float)

        # <x^3>.
        scl = self._beta[ind.X, j-1]**(3e0/2e0)/64e0

        dx_dJ_k[0] = -3e0*scl*nlbd.compute_dx_dJ_fp_I(j, self._I_21000)
        dx_dJ_k[1] = -scl*nlbd.compute_dx_dJ_fp_I(j, self._I_30000)
        dx_dJ_k[2] = 4e0*scl*nlbd.compute_dx_dJ_fp_I(j, self._I_10110)

        # <x*y^2>.
        scl = np.sqrt(self._beta[ind.X, j-1])*self._beta[ind.Y, j-1]/64e0
        dx_dJ_k[3] = 8e0*scl*nlbd.compute_dx_dJ_fp_I(j, self._I_21000)
        dx_dJ_k[4] = -4e0*scl*nlbd.compute_dx_dJ_fp_I(j, self._I_10200)
        dx_dJ_k[5] = 4e0*scl*nlbd.compute_dx_dJ_fp_I(j, self._I_10020)

        dx_dJ_k[6] = -4e0*scl*nlbd.compute_dx_dJ_fp_I(j, self._I_10110)
        dx_dJ_k[7] = -scl*nlbd.compute_dx_dJ_fp_I(j, self._I_10200)
        dx_dJ_k[8] = -scl*nlbd.compute_dx_dJ_fp_I(j, self._I_10020)

        return dx_dJ_k

    def compute_dx_dJ(self, dx_dJ_k):
        x3 = new.tpsa()
        xy2 = new.tpsa()
        x3.set(self._I_11000, 0e0, dx_dJ_k[0]+dx_dJ_k[1])
        x3.set(self._I_00110, 0e0, dx_dJ_k[2])
        xy2.set(self._I_11000, 0e0, dx_dJ_k[3]+dx_dJ_k[4]+dx_dJ_k[5])
        xy2.set(self._I_00110, 0e0, dx_dJ_k[6]+dx_dJ_k[7]+dx_dJ_k[8])
        return x3, xy2

    def compute_dnu_dJ(self):
        K = [new.tpsa(), new.tpsa()]
        for j, ele in enumerate(self._lat_prop._lattice):
            if type(ele) == ts.Sextupole:
                L = ele.get_length()
                b3xL = ele.get_multipoles().get_multipole(MpoleInd.sext).real*L
                if b3xL != 0e0:
                    nlbd._M = compute_map(lat_prop, j)
                    nlbd.compute_map_normal_form()
                    if not True:
                        x3_re, xy2_re = self.compute_dx_dJ_PB(j)
                    else:
                        dx_dJ_k = self.compute_dx_dJ_k_fp(j)
                        x3_re, xy2_re = self.compute_dx_dJ(dx_dJ_k)
                    K[0] -= b3xL*(x3_re+0*xy2_re.get(self._I_11110)
                                  *get_mn(self._I_11110))
                    K[1] -= b3xL*xy2_re
        K[0].print("K")
        K[1].print("K")
       

def compute_map_k(lat_prop, k):
    # Compute Poincaré map at location k.
    map = new.ss_vect_tpsa()
    map.set_identity()
    n = len(lat_prop._lattice)
    lat_prop._lattice.propagate(lat_prop._model_state, map, k, n-k)
    if k != 0:
        lat_prop._lattice.propagate(lat_prop._model_state, map, 0, k)
    return map


def compute_M_k(lat_prop, k):
    # Compute transer map for element k.
    M = new.ss_vect_tpsa()
    M.set_identity()
    n = len(lat_prop._lattice)
    lat_prop._lattice.propagate(lat_prop._model_state, M, k, 1)
    return M


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


def get_cmn(c, I):
    cmn = new.ctpsa()
    cmn.set(I, 0e0, c)
    return cmn


def RtoCmplx(a_re, a_im):
    a = new.ctpsa()
    for k, coeff in enumerate(a_re.get_coefficients()):
        I = coeff[0]
        c = coeff[1]
        a.set(I, 0e0, c)
    for k, coeff in enumerate(a_im.get_coefficients()):
        I = coeff[0]
        c = a_re.get(I) + 1j*coeff[1]
        a.set(I, 0e0, c)
    return a


def tpsa_mul(scl, a):
    b = new.tpsa()
    for k, coeff in enumerate(a.get_coefficients()):
        I = coeff[0]
        c = coeff[1]
        b.set(I, 0e0, scl*c)
    return b


def compute_alpha_beta(A_tpsa):
    A = A_tpsa.jacobian()[:4, :4]
    A_A_tp = A @ A.transpose()
    alpha = []
    beta = []
    for k in range(2):
        alpha.append(-A_A_tp[2*k, 2*k+1])
        beta.append(A_A_tp[2*k, 2*k])
    return np.array(alpha), np.array(beta)


def get_ij(i, j, A):
    I = [0, 0, 0, 0, 0]
    I[j] = 1
    return A.iloc[i].get(I)


def set_ij(i, j, a, A):
    I = [0, 0, 0, 0, 0]
    I[j] = 1
    A.iloc[i].set(I, 0e0, a)
    return A


def compute_dnu(n_dof, A):
    eps = 1e-15;
    dnu = []
    for k in range(n_dof):
        if k < 2:
            dnu.append(np.arctan2(get_ij(2*k, 2*k+1, A), get_ij(2*k, 2*k, A))
                       /(2e0*np.pi))
        elif n == 3:
            dnu.append(-np.arctan2(get_ij(ct, ind.delta, A),
                                   get_ij(ct, ind.ct, A))/(2e0*np.pi))
            if dnu[-1] < -eps:
                dnu[-1] += 1e0;
    return np.array(dnu)


def compute_A_CS(n_dof, A):
    A_CS = new.ss_vect_tpsa()
    R = new.ss_vect_tpsa()
    R.set_identity();
    dnu = compute_dnu(n_dof, A)
    for k in range(n_dof):
        c = np.cos(2e0*np.pi*dnu[k])
        s = np.sin(2e0*np.pi*dnu[k]);
        set_ij(2*k, 2*k, c, R)
        set_ij(2*k, 2*k+1, -s, R)
        set_ij(2*k+1, 2*k, s, R)
        set_ij(2*k+1, 2*k+1, c, R)
        A_CS.compose(A, R)
    return dnu, A_CS;


def compute_A_0_k(j):
    A = new.ss_vect_tpsa()
    A.set_identity()
    set_ij(ind.x,     ind.delta,  nlbd._eta[ind.x,  j], A)
    set_ij(ind.px,    ind.delta,  nlbd._eta[ind.px, j], A)
    set_ij(ind.delta, ind.x,     -nlbd._eta[ind.px, j], A)
    set_ij(ind.delta, ind.px,     nlbd._eta[ind.x,  j], A)
    return A


def compute_A_1_k(j):
    A = new.ss_vect_tpsa()
    A.set_identity()
    for k in range(2):
        set_ij(2*k,   2*k,   np.sqrt(nlbd._beta[k, j]),                    A)
        set_ij(2*k+1, 2*k,   -nlbd._alpha[k, j]/np.sqrt(nlbd._beta[k, j]), A)
        set_ij(2*k+1, 2*k+1, 1e0/np.sqrt(nlbd._beta[k, j]),                A)
    return A


def prt_alpha_beta(lat_prop, k, A):
    name = lat_prop._lattice[k].name
    alpha, beta = compute_alpha_beta(nlbd._A)
    print(f"   {k:2d} {name:10s}{nlbd._s[k]:5.3f} {alpha[ind.X]:8.5f}"
          f" {beta[ind.X]:7.5f}  {alpha[ind.Y]:8.5f} {beta[ind.Y]:7.5f}")


def map_trunc(no_cut, A):
    for k in range(7):
        A.iloc[k].get_description().truncate(no_cut)
    return A


def chk_g(lat_prop, nlbd):
    g_re = new.tpsa()
    g_im = new.tpsa()
    M = new.ss_vect_tpsa()
    M_j_inv = new.ss_vect_tpsa()
    A = new.ss_vect_tpsa()
    A_inv = new.ss_vect_tpsa()
    R = new.ss_vect_tpsa()
    R_inv = new.ss_vect_tpsa()

    # Can't hash a list (mutable) - but can hash a tuple (immutable).
    ind_list = {}
    ind_list[nlbd._I_30000_tuple] = nlbd._I_30000
    ind_list[nlbd._I_21000_tuple] = nlbd._I_21000
    ind_list[nlbd._I_10110_tuple] = nlbd._I_10110
    ind_list[nlbd._I_10200_tuple] = nlbd._I_10200
    ind_list[nlbd._I_10020_tuple] = nlbd._I_10020

    g_0 = nlbd._g
    print("\n")
    # for j, ele in enumerate(nlbd._lat_prop._lattice):
    for j in range(1, len(nlbd._lat_prop._lattice)-1):
        # Compute map at element entrance.
        if not True:
            nlbd._M = compute_map_k(lat_prop, j)
        else:
            M_j = compute_M_k(lat_prop, j)
            M_j_inv.inv(M_j)
            nlbd._M.compose(nlbd._M, M_j_inv)
            nlbd._M.compose(M_j, nlbd._M)

        nlbd.compute_map_normal_form()
        g_1 = nlbd._g

        A_0 = compute_A_0_k(j-1)
        A_1 = compute_A_1_k(j-1)
        A.compose(A_0, A_1)
        A_inv.inv(A)
        R.compose(M_j, A)
        A_0 = compute_A_0_k(j)
        A_1 = compute_A_1_k(j)
        A.compose(A_0, A_1)
        A_inv.inv(A)
        R.compose(A_inv, R)
        R_inv.inv(R)

        nlbd._g = g_0
        nlbd.propagate_g(j)
        g_2 = nlbd._g
        # g = compose_bs(g, R)

        trunc(nlbd._g, 3).print("nlbd._g")
        trunc(g, 3).print("g")

        nlbd._g.CtoR(g_re, g_im)

        if False:
            prt_alpha_beta(lat_prop, j, nlbd._A_1)

        g_re_fp, g_im_fp = nlbd.compute_g_fp(j)

        print(f"\n  {j:2d}"
              f"  nu = [{nlbd._dnu[ind.X, j]:7.5f},"
              f" {nlbd._dnu[ind.Y, j]:7.5f}]"
              f"  mu = [{nlbd._dmu[ind.X, j]:7.5f},"
              f" {nlbd._dmu[ind.Y, j]:7.5f}]"
              f"  alpha = [{nlbd._alpha[ind.X, j]:7.5f},"
              f" {nlbd._alpha[ind.Y, j]:8.5f}]"
              f"  beta = [{nlbd._beta[ind.X, j]:7.5f},"
              f" {nlbd._beta[ind.Y, j]:8.5f}]")

        for key, value in ind_list.items():
            str = nlbd.I2str(ind_list[key])
            c = complex(g_re_fp.get(ind_list[key]), g_im_fp.get(ind_list[key]))
            print(f"  g_{str:s} = ", end="")
            nlbd.prt_ampl_phase("", c, nl="")
            c = complex(g_re.get(ind_list[key]), g_im.get(ind_list[key]))
            nlbd.prt_ampl_phase("\t\t", c, nl="\n")


def chk_it(nlbd):
    x3_re = new.tpsa()
    x3_im = new.tpsa()
    xy2_re = new.tpsa()
    xy2_im = new.tpsa()
    K = [new.tpsa(), new.tpsa()]
    
    I_11000 = [1, 1, 0, 0, 0]
    I_00110 = [0, 0, 1, 1, 0]

    print("\n")
    for j, ele in enumerate(nlbd._lat_prop._lattice):
        if type(ele) == ts.Sextupole:
            L = ele.get_length()
            b3xL = ele.get_multipoles().get_multipole(MpoleInd.sext).real*L
            if b3xL != 0e0:
                nlbd._M = compute_map(lat_prop, j)
                nlbd.compute_map_normal_form()
                # g = nlbd.propagate_g(j)

                x1_re, x3_re, xy2_re = nlbd.compute_dx_dJ_PB(j)

                dx_dJ_k = nlbd.compute_dx_dJ_k_fp(j)
                x3_fp, xy2_fp = nlbd.compute_dx_dJ(dx_dJ_k)
                x_fp_re = x3_fp + xy2_fp

                print(f"  {dx_dJ_k[1]:10.3e}"
                      f"\t {x3_re.get(I_11000):10.3e}")

                # print(f"  {x_fp_re.get(I_11000):10.3e}"
                #       f" ({x_re.get(I_11000):10.3e}"
                #       f" = {x1_re.get(I_11000):10.3e}"
                #       f" {x3_re.get(I_11000):+10.3e}"
                #       f" {xy2_re.get(I_11000):+10.3e})"
                #       f" {x_fp_re.get(I_00110):10.3e}"
                #       f" ({x_re.get(I_00110):10.3e}"
                #       f" = {x1_re.get(I_00110):10.3e}"
                #       f" {x3_re.get(I_00110):+10.3e}"
                #       f" {xy2_re.get(I_00110):+10.3e})")


def compute_dx_dJ_g_k(nlbd, j):
    ps_dim = 4;

    g_re = new.tpsa()
    g_im = new.tpsa()

    dx_dJ_g_k = np.zeros(9, dtype=complex)

    g = nlbd.propagate_g(j)
    g.CtoR(g_re, g_im)

    # A x4 to correct for cosine & sine terms - for both f & g -> [f, g].
    scl = -1j*4e0/8e0

    # <x^3>.

    g_k_re = g_re.get([2, 1, 0, 0, 0])*get_mn([2, 1, 0, 0, 0])
    g_k_im = g_im.get([2, 1, 0, 0, 0])*get_mn([2, 1, 0, 0, 0])
    dx_dJ_g_k[0] = scl*complex(
        g_re.poisbra(get_mn([1, 2, 0, 0, 0]), ps_dim).get([2, 2, 0, 0, 0]),
        g_im.poisbra(get_mn([1, 2, 0, 0, 0]), ps_dim).get([2, 2, 0, 0, 0]))
    g_re_k = g_re.get([3, 0, 0, 0, 0])*get_mn([3, 0, 0, 0, 0])
    g_im_k = g_im.get([3, 0, 0, 0, 0])*get_mn([3, 0, 0, 0, 0])
    dx_dJ_g_k[1] = scl*complex(
        g_re.poisbra(get_mn([0, 3, 0, 0, 0]), ps_dim).get([2, 2, 0, 0, 0])/3e0,
        g_im.poisbra(get_mn([0, 3, 0, 0, 0]), ps_dim).get([2, 2, 0, 0, 0])/3e0)

    g_re_k = g_re.get([1, 0, 1, 1, 0])*get_mn([1, 0, 1, 1, 0])
    g_im_k = g_im.get([1, 0, 1, 1, 0])*get_mn([1, 0, 1, 1, 0])
    dx_dJ_g_k[2] = scl*complex(
        g_re.poisbra(get_mn([1, 2, 0, 0, 0]), ps_dim).get([1, 1, 1, 1, 0]),
        g_im.poisbra(get_mn([1, 2, 0, 0, 0]), ps_dim).get([1, 1, 1, 1, 0]))

    # <x*y^2>.

    g_re_k = g_re.get([2, 1, 0, 0, 0])*get_mn([2, 1, 0, 0, 0])
    g_im_k = g_im.get([2, 1, 0, 0, 0])*get_mn([2, 1, 0, 0, 0])
    dx_dJ_g_k[3] = scl*complex(
        g_re.poisbra(2e0*get_mn([0, 1, 1, 1, 0]), ps_dim).get([1, 1, 1, 1, 0]),
        g_im.poisbra(2e0*get_mn([0, 1, 1, 1, 0]), ps_dim).get([1, 1, 1, 1, 0]))
    g_re_k = g_re.get([1, 0, 2, 0, 0])*get_mn([1, 0, 2, 0, 0])
    g_im_k = g_im.get([1, 0, 2, 0, 0])*get_mn([1, 0, 2, 0, 0])
    dx_dJ_g_k[4] = scl*complex(
        g_re.poisbra(get_mn([0, 1, 0, 2, 0]), ps_dim).get([1, 1, 1, 1, 0]),
        g_im.poisbra(get_mn([0, 1, 0, 2, 0]), ps_dim).get([1, 1, 1, 1, 0]))
    g_re_k = g_re.get([1, 0, 0, 2, 0])*get_mn([1, 0, 0, 2, 0])
    g_im_k = g_im.get([1, 0, 0, 2, 0])*get_mn([1, 0, 0, 2, 0])
    dx_dJ_g_k[5] = scl*complex(
        g_re.poisbra(get_mn([0, 1, 2, 0, 0]), ps_dim).get([1, 1, 1, 1, 0]),
        g_im.poisbra(get_mn([0, 1, 2, 0, 0]), ps_dim).get([1, 1, 1, 1, 0]))

    g_re_k = g_re.get([1, 0, 1, 1, 0])*get_mn([1, 0, 1, 1, 0])
    g_im_k = g_im.get([1, 0, 1, 1, 0])*get_mn([1, 0, 1, 1, 0])
    dx_dJ_g_k[6] = scl*complex(
        g_re.poisbra(2e0*get_mn([0, 1, 1, 1, 0]), ps_dim).get([0, 0, 2, 2, 0]),
        g_im.poisbra(2e0*get_mn([0, 1, 1, 1, 0]), ps_dim).get([0, 0, 2, 2, 0]))
    g_re_k = g_re.get([1, 0, 2, 0, 0])*get_mn([1, 0, 2, 0, 0])
    g_im_k = g_im.get([1, 0, 2, 0, 0])*get_mn([1, 0, 2, 0, 0])
    dx_dJ_g_k[7] = scl*complex(
        g_re.poisbra(get_mn([0, 1, 0, 2, 0]), ps_dim).get([0, 0, 2, 2, 0]),
        g_im.poisbra(get_mn([0, 1, 0, 2, 0]), ps_dim).get([0, 0, 2, 2, 0]))
    g_re_k = g_re.get([1, 0, 0, 2, 0])*get_mn([1, 0, 0, 2, 0])
    g_im_k = g_im.get([1, 0, 0, 2, 0])*get_mn([1, 0, 0, 2, 0])
    dx_dJ_g_k[8] = scl*complex(
        g_re.poisbra(get_mn([0, 1, 2, 0, 0]), ps_dim).get([0, 0, 2, 2, 0]),
        g_im.poisbra(get_mn([0, 1, 2, 0, 0]), ps_dim).get([0, 0, 2, 2, 0]))

    return dx_dJ_g_k


def prt_ampl_phase(str, a, nl="\n"):
    print(f"{str:s}{abs(a):10.3e} \u2220 {np.rad2deg(np.angle(a)):6.1f}",
          end=nl)


def prt_terms(
        nlbd, x3_g_PB, xy2_g_PB, dx_dJ_g_k, x3_g, xy2_g, dx_dJ_k_fp, x3_fp,
        xy2_fp):
    print("\n<x^3>:")
    prt_ampl_phase("  s_22000 = ", x3_g_PB.get([2, 2, 0, 0, 0]))
    prt_ampl_phase("            ", x3_g.get([2, 2, 0, 0, 0]), nl="")
    prt_ampl_phase(" = ", dx_dJ_g_k[0], nl="")
    prt_ampl_phase("", dx_dJ_g_k[1])
    prt_ampl_phase("            ", x3_fp.get([2, 2, 0, 0, 0]), nl="")
    prt_ampl_phase(" = ", dx_dJ_k_fp[0], nl="")
    prt_ampl_phase("", dx_dJ_k_fp[1])

    prt_ampl_phase("\n  s_11110 = ", x3_g_PB.get([1, 1, 1, 1, 0]))
    prt_ampl_phase("            ", x3_g.get([1, 1, 1, 1, 0]))
    prt_ampl_phase("            ", dx_dJ_k_fp[2])

    print("\n<x*y^2>:")
    prt_ampl_phase("  s_11110 = ", xy2_g_PB.get([1, 1, 1, 1, 0]))
    prt_ampl_phase("            ", xy2_g.get([1, 1, 1, 1, 0]), nl="")
    prt_ampl_phase(" = ", dx_dJ_g_k[3], nl="")
    prt_ampl_phase("", dx_dJ_g_k[4], nl="")
    prt_ampl_phase("", dx_dJ_g_k[5])
    prt_ampl_phase("            ", xy2_fp.get([1, 1, 1, 1, 0]), nl="")
    prt_ampl_phase(" = ", dx_dJ_k_fp[3], nl="")
    prt_ampl_phase("", dx_dJ_k_fp[4], nl="")
    prt_ampl_phase("", dx_dJ_k_fp[5])

    prt_ampl_phase("\n  s_00220 = ", xy2_g_PB.get([0, 0, 2, 2, 0]))
    prt_ampl_phase("            ", xy2_g.get([0, 0, 2, 2, 0]), nl="")
    prt_ampl_phase(" = ", dx_dJ_g_k[6], nl="")
    prt_ampl_phase("", dx_dJ_g_k[7], nl="")
    prt_ampl_phase("", dx_dJ_g_k[8])
    prt_ampl_phase("            ", xy2_fp.get([0, 0, 2, 2, 0]), nl="")
    prt_ampl_phase(" = ", dx_dJ_k_fp[6], nl="")
    prt_ampl_phase("", dx_dJ_k_fp[7], nl="")
    prt_ampl_phase("", dx_dJ_k_fp[8])


def chk_terms(nlbd, j):
    x3_g_PB, xy2_g_PB = nlbd.compute_dx_dJ_PB(j)

    dx_dJ_g_k = compute_dx_dJ_g_k(nlbd, j)
    x3_g, xy2_g = nlbd.compute_dx_dJ(dx_dJ_g_k)

    dx_dJ_k_fp = nlbd.compute_dx_dJ_k_fp(j)
    x3_fp, xy2_fp = nlbd.compute_dx_dJ(dx_dJ_k_fp)
    x3_fp.print("x3_fp")
    xy2_fp.print("xy2_fp")

    prt_terms(nlbd, x3_g_PB, xy2_g_PB, dx_dJ_g_k, x3_g, xy2_g, dx_dJ_k_fp,
              x3_fp, xy2_fp)


gtpsa_prop.no = 4
gtpsa_prop.desc = gtpsa.desc(gtpsa_prop.nv, gtpsa_prop.no)

cod_eps = 1e-15
E_0     = 3.0e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

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

nlbd._M = compute_map_k(lat_prop, 0)
print("\nM:", nlbd._M)

nlbd.compute_map_normal_form()

nlbd._K_re.print("K_re")

if False:
    compute_M_from_MNF(K, A_0, A_1, g)
if False:
    compute_R_from_MNF(M, A_0, A_1, g, K)

if False:
    if True:
        loc = 1
    else:
        # loc = 237
        loc = 357
    chk_terms(nlbd, loc)
    assert False

if not False:
    chk_g(lat_prop, nlbd)
    assert False

if False:
    chk_it(nlbd)
    assert False

nlbd.compute_dnu_dJ()
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

K, K_k = nlbd.compute_dnu_dJ()
print(f"  [{abs(K[0]):10.3e}  \u2220 {np.rad2deg(np.angle(K[0])):6.1f}"
      f", {abs(K[1]):10.3e}  \u2220 {np.rad2deg(np.angle(K[1])):6.1f}"
      f", {abs(K[2]):10.3e}]  \u2220 {np.rad2deg(np.angle(K[2])):6.1f}")
print(f"  [{abs(K_k[0]):10.3e}  \u2220 {np.rad2deg(np.angle(K_k[0])):6.1f}"
      f", {abs(K_k[1]):10.3e}  \u2220 {np.rad2deg(np.angle(K_k[1])):6.1f}"
      f", {abs(K_k[2]):10.3e}]  \u2220 {np.rad2deg(np.angle(K_k[2])):6.1f}")
