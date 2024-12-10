"""Class for nonlinear dynamics analysis.
"""

import enum

import numpy as np
from scipy import linalg as la

import thor_scsi.lib as ts
from thor_scsi.utils import index_class as ind, linear_optics as lo


ind = ind.index_class()

import gtpsa


class MpoleInd(enum.IntEnum):
    quad = 2
    sext = 3


class nonlin_dyn_class:
    # Private.

    def __init__(self, lat_prop, A_max, beta_inj, delta_max, b_3_list):
        self._A_max        = A_max
        self._beta_inj     = beta_inj
        self._delta_max    = delta_max
        self._twoJ         = np.nan
        self._Id_scl       = np.nan

        self._b_3_list     = b_3_list

        self._M            = np.nan

        self._h_im_scl_rms = np.nan
        self._K_re_scl_rms = np.nan

        self._h_re         = np.nan
        self._h_im         = np.nan
        self._K_re         = np.nan
        self._K_im         = np.nan

        self._h_re_scl     = np.nan
        self._h_im_scl     = np.nan
        self._K_re_scl     = np.nan
        self._K_im_scl     = np.nan

        if b_3_list != []:
            self.compute_twoJ()
            self.compute_Id_scl(lat_prop)
        else:
            print("\nnonlin_dyn_class: b_3_list is empty")

    # Public.

    h_dict = {"h_10002" : [1, 0, 0, 0, 2, 0, 0],
              "h_20001" : [2, 0, 0, 0, 1, 0, 0],
              "h_00201" : [0, 0, 2, 0, 1, 0, 0],

              "h_30000" : [3, 0, 0, 0, 0, 0, 0],
              "h_21000" : [2, 1, 0, 0, 0, 0, 0],
              "h_10110" : [1, 0, 1, 1, 0, 0, 0],
              "h_10200" : [1, 0, 2, 0, 0, 0, 0],
              "h_10020" : [1, 0, 0, 2, 0, 0, 0]}

    h_scl_dict = {"h_10002" : 1e-2,
                  "h_20001" : 1e-2,
                  "h_00201" : 1e-2,

                  "h_30000" : 1e-2,
                  "h_21000" : 1e-2,
                  "h_10110" : 1e-2,
                  "h_10200" : 1e-2,
                  "h_10020" : 1e-2}

    K_xi_dict = {"K_11001" : [1, 1, 0, 0, 1, 0, 0],
                 "K_00111" : [0, 0, 1, 1, 1, 0, 0]}

    K_geom_dict = {"K_22000" : [2, 2, 0, 0, 0, 0, 0],
                   "K_11110" : [1, 1, 1, 1, 0, 0, 0],
                   "K_00220" : [0, 0, 2, 2, 0, 0, 0],

                   "K_33000" : [3, 3, 0, 0, 0, 0, 0],
                   "K_22110" : [2, 2, 1, 1, 0, 0, 0],
                   "K_11220" : [1, 1, 2, 2, 0, 0, 0],
                   "K_00330" : [0, 0, 3, 3, 0, 0, 0]}

    K_geom_scl_dict = {"K_22000" : 1.0,
                       "K_11110" : 1.0,
                       "K_00220" : 1.0,

                       "K_33000" : 5.0,
                       "K_22110" : 5.0,
                       "K_11220" : 5.0,
                       "K_00330" : 5.0}

    K_chrom_dict = {"K_11002" : [1, 1, 0, 0, 2, 0, 0],
                    "K_00112" : [0, 0, 1, 1, 2, 0, 0],

                    "K_11003" : [1, 1, 0, 0, 3, 0, 0],
                    "K_00113" : [0, 0, 1, 1, 3, 0, 0]}

    K_chrom_scl_dict = {"K_11002" : 1.0,
                        "K_00112" : 1.0,

                        "K_11003" : 5.0,
                        "K_00113" : 5.0}

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

    def zero_mult(self, lat_prop, n):
        # Zero multipoles.
        print("\nZeroing multipoles: n = {:d}".format(n))
        for el in lat_prop._lattice:
            if (type(el) == ts.Bending) or (type(el) == ts.Quadrupole) \
               or (type(el) == ts.Sextupole) or (type(el) == ts.Octupole) \
               or (type(el) == ts.Multipole):
                el.get_multipoles().set_multipole(n, 0e0)

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

            b_3_val = []
            for k in range(n):
                b_3xL = \
                    lat_prop.get_b_nxL_elem(self._b_3_list[k], 0, MpoleInd.sext)
                lat_prop.set_b_nxL_fam(
                    self._b_3_list[k], MpoleInd.sext, b_3xL+db_3xL[k])
                b_3 = \
                    lat_prop.get_b_n_elem(self._b_3_list[k], 0, MpoleInd.sext)
                b_3_val.append(b_3)

            M = gtpsa.ss_vect_tpsa(
                lat_prop._desc, lat_prop._no, lat_prop._nv,
            index_mapping=lat_prop._named_index)
            self.compute_map(lat_prop, 2)
            stable, _, xi = \
                lo.compute_nu_xi(lat_prop._desc, lat_prop._no, self._M)

        return stable, xi, np.array(b_3_val)

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
        self._g_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        self._g_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        K          = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        self._K_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        self._K_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)

        self._M.Map_Norm(A_0, A_1, R, g, K)
        K.CtoR(self._K_re, self._K_im)

    def scl_h(self, h, h_scl_dict, h_dict):
        for key in h_dict:
            h.set(h_dict[key], 0e0, h_scl_dict[key]*h.get(h_dict[key]))

    def compute_rms(self, h, h_dict):
        rms = 0e0
        for key in h_dict:
            rms += h.get(h_dict[key])**2
        return np.sqrt(rms)

    def compute_nl(self, lat_prop):
        self.compute_h(lat_prop)
        self.compute_map_normal_form(lat_prop)

        self._h_re_scl = self.compose_bs(lat_prop, self._h_re, self._Id_scl)
        self._h_im_scl = self.compose_bs(lat_prop, self._h_im, self._Id_scl)
        self._K_re_scl = self.compose_bs(lat_prop, self._K_re, self._Id_scl)
        self._K_im_scl = self.compose_bs(lat_prop, self._K_im, self._Id_scl)

        self.scl_h(self._h_re_scl, self.h_scl_dict, self.h_dict)
        self.scl_h(self._h_im_scl, self.h_scl_dict, self.h_dict)

        self.scl_h(self._K_re_scl, self.K_geom_scl_dict, self.K_geom_dict)
        self.scl_h(self._K_re_scl, self.K_chrom_scl_dict, self.K_chrom_dict)
        self.scl_h(self._K_im_scl, self.K_geom_scl_dict, self.K_geom_dict)
        self.scl_h(self._K_im_scl, self.K_chrom_scl_dict, self.K_chrom_dict)

        self._h_re_scl_rms = self.compute_rms(self._h_re_scl, self.h_dict)
        self._h_im_scl_rms = self.compute_rms(self._h_im_scl, self.h_dict)

        K_geom = self.compute_rms(self._K_re_scl, self.K_geom_dict)
        K_chrom = self.compute_rms(self._K_re_scl, self.K_chrom_dict)

        self._K_re_scl_rms = np.sqrt(K_geom**2+K_chrom**2)

    def prt_nl(self, lat_prop):
        print("\n  nu = [{:7.5f}, {:7.5f}]".format(
            -self._K_re.get([1, 1, 0, 0, 0, 0, 0])/np.pi,
            -self._K_re.get([0, 0, 1, 1, 0, 0, 0])/np.pi))
        print("  xi = [{:5.3f}, {:5.3f}]\n".format(
            -self._K_re.get([1, 1, 0, 0, 1, 0, 0])/np.pi,
            -self._K_re.get([0, 0, 1, 1, 1, 0, 0])/np.pi))
        for key in self.K_xi_dict:
            print("  {:s}  = {:10.3e}".
                  format(key, self._K_re.get(self.K_xi_dict[key])))

        print("\n  h_re rms = {:9.3e}".format(self._h_re_scl_rms))
        print("  h_im rms = {:9.3e}".format(self._h_im_scl_rms))
        print("  K_re rms = {:9.3e}".format(self._K_re_scl_rms))

        print("\nRe{h}:")
        for key in self.h_dict:
            print("  {:s}  = {:10.3e}".
                  format(key, self._h_re_scl.get(self.h_dict[key])))
            if key == "h_00201":
                print()

        print("\nIm{h}:")
        for key in self.h_dict:
            print("  {:s}  = {:10.3e}".
                  format(key, self._h_im_scl.get(self.h_dict[key])))
            if key == "h_00201":
                print()

        print("\nRe{K}:")
        for key in self.K_geom_dict:
            print("  {:s}  = {:10.3e}".
                  format(key, self._K_re_scl.get(self.K_geom_dict[key])))
            if key == "K_00220":
                if lat_prop._no < 6:
                    break
                else:
                    print()

        print()
        for key in self.K_chrom_dict:
            print("  {:s}  = {:10.3e}".
                  format(key, self._K_re_scl.get(self.K_chrom_dict[key])))
            if key == "K_00112":
                if lat_prop._no < 6:
                    break
                else:
                    print()


__all__ = [nonlin_dyn_class]
