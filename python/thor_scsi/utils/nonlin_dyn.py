"""Class for nonlinear dynamics analysis.
"""

import enum

import numpy as np
from scipy import linalg as la

from thor_scsi.utils import index_class as ind, linear_optics as lo


ind = ind.index_class()

import gtpsa


class MpoleInd(enum.IntEnum):
    quad = 2
    sext = 3


class nonlin_dyn_class:
    # Private.

    def __init__(self, lat_prop, no, A_max, beta_inj, delta_max, b_3_list):
        self._no        = no
        self._A_max     = A_max
        self._beta_inj  = beta_inj
        self._delta_max = delta_max
        self._twoJ      = np.nan
        self._Id_scl    = np.nan

        self._b_3_list  = b_3_list

        self._M         = np.nan

        self._h_im_rms  = np.nan
        self._K_re_rms  = np.nan
        self._h_re      = np.nan
        self._h_im      = np.nan
        self._K_re      = np.nan
        self._K_im      = np.nan

        self.compute_twoJ()
        self.compute_Id_scl(lat_prop)

    # Public.

    h_re_dict = {
        "h_22000" : [2, 2, 0, 0, 0, 0, 0],
        "h_11110" : [1, 1, 1, 1, 0, 0, 0],
        "h_00220" : [0, 0, 2, 2, 0, 0, 0]
    }

    h_im_dict = {
        "h_10002" : [1, 0, 0, 0, 2, 0, 0],
        "h_20001" : [2, 0, 0, 0, 1, 0, 0],
        "h_00201" : [0, 0, 2, 0, 1, 0, 0],

        "h_30000" : [3, 0, 0, 0, 0, 0, 0],
        "h_21000" : [2, 1, 0, 0, 0, 0, 0],
        "h_10110" : [1, 0, 1, 1, 0, 0, 0],
        "h_10200" : [1, 0, 2, 0, 0, 0, 0],
        "h_10020" : [1, 0, 0, 2, 0, 0, 0]
    }

    h_re_dict = {
        "h_22000" : [2, 2, 0, 0, 0, 0, 0],
        "h_11110" : [1, 1, 1, 1, 0, 0, 0],
        "h_00220" : [0, 0, 2, 2, 0, 0, 0]
    }

    h_im_dict = {
        "h_10002" : [1, 0, 0, 0, 2, 0, 0],
        "h_20001" : [2, 0, 0, 0, 1, 0, 0],
        "h_00201" : [0, 0, 2, 0, 1, 0, 0],

        "h_30000" : [3, 0, 0, 0, 0, 0, 0],
        "h_21000" : [2, 1, 0, 0, 0, 0, 0],
        "h_10110" : [1, 0, 1, 1, 0, 0, 0],
        "h_10200" : [1, 0, 2, 0, 0, 0, 0],
        "h_10020" : [1, 0, 0, 2, 0, 0, 0]
    }

    K_dict = {"K_22000" : [2, 2, 0, 0, 0, 0, 0],
              "K_11110" : [1, 1, 1, 1, 0, 0, 0],
              "K_00220" : [0, 0, 2, 2, 0, 0, 0],

              "K_11002" : [1, 1, 0, 0, 2, 0, 0],
              "K_00112" : [0, 0, 1, 1, 2, 0, 0]}

    def compute_twoJ(self):
        self._twoJ = \
            np.array(
                [self._A_max[ind.X]**2/self._beta_inj[ind.X],
                 self._A_max[ind.Y]**2/self._beta_inj[ind.Y]])

    def compute_Id_scl(self, lat_prop):
        self._Id_scl = \
            gtpsa.ss_vect_tpsa(
                lat_prop._desc, lat_prop._no,
                index_mapping=lat_prop._named_index)
        self._Id_scl.set_identity()
        for k in range(4):
            self._Id_scl.iloc[k].set_variable(
                0e0, k+1, np.sqrt(self._twoJ[k//2]))
        self._Id_scl.delta.set_variable(0e0, 5, self._delta_max)

    def compose_bs(self, lat_prop, h, map):
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

    def compute_map(self, lat_prop, no):
        self._M = lo.compute_map(
            lat_prop._lattice, lat_prop._model_state, desc=lat_prop._desc,
            tpsa_order=no)

    def zero_sext(lat_prop, b_3_list):
        # Zero sextupoles.
        print("\nZeroing sextupoles.")
        for b3_name in b_3_list:
            lat_prop.set_b_n_fam(b3_name, MpoleInd.sext, 0e0)

    def compute_sext_resp_mat_num(self, lat_prop, b_3_list):
        # Uses integrated sextupole strength.
        db_3xL = 1e-3
        n = len(b_3_list)
        A = np.zeros((2, n))
        self.compute_map(lat_prop, 2)
        stable, _, xi = lo.compute_nu_xi(lat_prop._desc, lat_prop._no, self._M)
        if stable:
            for k in range(n):
                b_3xL = \
                    lat_prop.get_b_nxL_elem(b_3_list[k], 0, MpoleInd.sext)
                lat_prop.set_b_nxL_fam(b_3_list[k], MpoleInd.sext, b_3xL-db_3xL)
                self.compute_map(lat_prop, 2)
                stable, _, xi_1 =  \
                    lo.compute_nu_xi(lat_prop._desc, lat_prop._no, self._M)
                if stable:
                    lat_prop.set_b_nxL_fam(
                        b_3_list[k], MpoleInd.sext, b_3xL+db_3xL)
                    self.compute_map(lat_prop, 2)
                    stable, _, xi_2 = \
                        lo.compute_nu_xi(lat_prop._desc, lat_prop._no, self._M)
                    a_ij = (xi_2-xi_1)/(2e0*db_3xL)
                    A[ind.X, k] = a_ij[ind.X]
                    A[ind.Y, k] = a_ij[ind.Y]
                else:
                    break
            return stable, xi, A
        else:
             return False, np.nan, np.nan

    def set_xi(self, lat_prop, xi_x, xi_y):
        n = len(self._b_3_list)

        stable, xi, A = self.compute_sext_resp_mat_num(lat_prop, self._b_3_list)
        if stable:
            A_inv = la.pinv(A)
            db_3xL = A_inv @ (-xi)

            for k in range(n):
                b_3xL = \
                    lat_prop.get_b_nxL_elem(self._b_3_list[k], 0, MpoleInd.sext)
                lat_prop.set_b_nxL_fam(
                    self._b_3_list[k], MpoleInd.sext, b_3xL+db_3xL[k])
                b_3 = lat_prop.get_b_n_elem(self._b_3_list[k], 0, MpoleInd.sext)

                M = gtpsa.ss_vect_tpsa(
                    lat_prop._desc, lat_prop._no, lat_prop._nv,
                index_mapping=lat_prop._named_index)
                self.compute_map(lat_prop, 2)
                stable, _, xi = \
                    lo.compute_nu_xi(lat_prop._desc, lat_prop._no, self._M)

        return stable

    def compute_h(self, lat_prop):
        h    = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        self._h_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        self._h_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)

        self._M.M_to_h_DF(h)
        h.CtoR(self._h_re, self._h_im)

    def compute_map_normal_form(self, lat_prop):
        A_0        = gtpsa.ss_vect_tpsa(lat_prop._desc, lat_prop._no)
        A_1        = gtpsa.ss_vect_tpsa(lat_prop._desc, lat_prop._no)
        R          = gtpsa.ss_vect_tpsa(lat_prop._desc, lat_prop._no)
        g          = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        g_re       = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        g_im       = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        K          = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        self._K_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        self._K_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)

        self._M.Map_Norm(A_0, A_1, R, g, K)
        K.CtoR(self._K_re, self._K_im)

    def compute_rms(self, h, dict):
        var = 0e0
        for key in dict:
            var += h.get(dict[key])**2
        return np.sqrt(var)

    def compute_nl(self, lat_prop):
        # Compute map to order no.
        self.compute_map(lat_prop, self._no)

        self.compute_h(lat_prop)
        self.compute_map_normal_form(lat_prop)

        self._h_re = self.compose_bs(lat_prop, self._h_re, self._Id_scl)
        self._h_im = self.compose_bs(lat_prop, self._h_im, self._Id_scl)
        self._K_re = self.compose_bs(lat_prop, self._K_re, self._Id_scl)
        self._K_im = self.compose_bs(lat_prop, self._K_im, self._Id_scl)

        self._h_im_rms = self.compute_rms(self._h_im, self.h_im_dict)
        self._K_re_rms = self.compute_rms(self._K_re, self.K_dict)

    def prt_nl(self):
        print("\n  h_im rms = {:9.3e}".format(self._h_im_rms))
        print("  K_re rms = {:9.3e}".format(self._K_re_rms))

        print()
        for key in self.h_im_dict:
            print("  {:s}  = {:10.3e}".
                  format(key, self._h_im.get(self.h_im_dict[key])))
            if key == "h_00201":
                print()
        print()
        for key in self.h_re_dict:
            print("  {:s}  = {:10.3e}".
                  format(key, self._h_re.get(self.h_re_dict[key])))
            if key == "h_00201":
                print()
        print()
        for key in self.K_dict:
            print("  {:s}  = {:10.3e}".
                  format(key, self._K_re.get(self.K_dict[key])))
            if key == "K_00220":
                print()


__all__ = [nonlin_dyn_class]
