"""
"""

import numpy as np


quad = 2


class opt_phi_class:
    # Private

    def __init__(self, lat_prop, phi_max):
        self._phi_max = phi_max

    # Public.



class bend_prm_class:
    # Private

    def __init__(
            self, lat_prop, dip_0, dip_list, b_2_prm, phi_max, b_2_bend_max):
        self._dip_0        = dip_0
        self._phi_0        = np.nan
        self._dip_list     = dip_list
        self._b_2_prm      = b_2_prm

        self._phi_max      = phi_max
        self._b_2_bend_max = b_2_bend_max

    # Public.

    def compute_phi_bend(self, lat_prop):
        # Compute total bend angle.
        phi = 0e0
        for k in range(len(self._dip_list)):
            phi += lat_prop.get_phi_elem(self._dip_list[k], 0)
        phi += lat_prop.get_phi_elem(self._dip_0, 0)
        return phi

    def get_bend_prm(self, lat_prop):
        prm = []
        bounds = []
        for k in range(len(self._dip_list)):
            dphi = lat_prop.get_phi_elem(self._dip_list[k], 0)
            prm.append(dphi)
            bounds.append((-self._phi_max, self._phi_max))
            if self._b_2_prm:
                b_2 = lat_prop.get_b_n_elem(self._dip_list[k], 0, quad)
                prm.append(b_2)
                bounds.append((-self._b_2_bend_max, self._b_2_bend_max))
                self._phi_0 = lat_prop.get_phi_elem(self._dip_0, 0)
        if self._b_2_prm:
            b_2 = lat_prop.get_b_n_elem(self._dip_0, 0, quad)
            prm.append(b_2)
            bounds.append((-self._b_2_bend_max, self._b_2_bend_max))
        return prm, bounds

    def set_bend_prm(self, lat_prop, prm, prm_ind):
        phi_bend = self.compute_phi_bend(lat_prop)
        phi = 0e0
        for k in range(len(self._dip_list)):
            lat_prop.set_phi_fam(self._dip_list[k], prm[prm_ind])
            phi += prm[prm_ind]
            prm_ind += 1
            if self._b_2_prm:
                lat_prop.set_b_n_fam(self._dip_list[k], quad, prm[prm_ind])
                prm_ind += 1
        self._phi_0 = phi_bend - phi
        lat_prop.set_phi_fam(self._dip_0, self._phi_0)
        if self._b_2_prm:
            lat_prop.set_b_n_fam(self._dip_0, quad, prm[prm_ind])
            prm_ind += 1
        return prm_ind


class rev_bend_prm_class:
    # Private

    def __init__(self, lat_prop, rb_list, dip_list, phi_max):
        self._rb_list  = rb_list
        self._dip_list = dip_list
        self._n_rb     = np.nan
        self._n_dip    = np.nan
        self._phi_list = np.nan

        self._phi_max  = phi_max

    # Public.

    def get_rev_bend_prm(self, lat_prop):
        self._n_rb = 0
        for k in range(len(self._rb_list)):
            self._n_rb += lat_prop.get_n_kid(self._rb_list[k])
        self._n_dip = 0
        for k in range(len(self._dip_list)):
            self._n_dip += lat_prop.get_n_kid(self._dip_list[k])

        phi_rb = lat_prop.get_phi_elem(self._rb_list[0], 0)
        bounds = (-self._phi_max, self._phi_max)
        return phi_rb, bounds

    def get_bend_angles(self, lat_prop):
        # Get bend angles without reverse bends.
        phi_rb = \
            lat_prop.get_phi_elem(self._rb_list[0], 0)*self._n_rb/self._n_dip
        phi_list = []
        for k in range(len(self._dip_list)):
            phi_list.append(
                lat_prop.get_phi_elem(self._dip_list[k], 0)+phi_rb)
        return np.array(phi_list)

    def set_rev_bend_prm(self, lat_prop, prm, prm_ind):
        phi_list = self.get_bend_angles(lat_prop)
        for k in range(len(self._rb_list)):
           lat_prop.set_phi_fam(self._rb_list[k], prm[prm_ind])
        dphi = self._n_rb*prm[prm_ind]/self._n_dip
        for k in range(len(self._dip_list)):
           lat_prop.set_phi_fam(self._dip_list[k], phi_list[k]-dphi)
        self._prm_0 = prm[prm_ind]
        prm_ind += 1
        return prm_ind


class prm_class(bend_prm_class, rev_bend_prm_class):
    # Private

    def __init__(self, lat_prop, prm_list, b_2_max):
        # Dictionary of parameter types and corresponding get functions.
        self._get_prm_func_dict = {
            "L":   lat_prop.get_L_elem,
            "L_b": lat_prop.get_L_elem,
            "phi": lat_prop.get_phi_elem,
            "b_2": lat_prop.get_b_2_elem
        }

        # Dictionary of parameter types and corresponding set functions.
        self._set_prm_func_dict = {
            "L":   lat_prop.set_L_fam,
            "L_b": lat_prop.set_L_bend_fam,
            "phi": lat_prop.set_phi_fam,
            "b_2": lat_prop.set_b_2_fam
        }

        self._b_2_max  = b_2_max
        self._prm_list = prm_list

    # Public.

    def get_prm(self, lat_prop):
        prm = []
        bounds = []
        for k in range(len(self._prm_list)):
            if self._prm_list[k][0] == "bend":
                p, b = self._prm_list[k][1].get_bend_prm(lat_prop)
                prm.extend(p)
                bounds.extend(b)
            elif self._prm_list[k][0] == "rev_bend":
                p, b = self._prm_list[k][1].get_rev_bend_prm(lat_prop)
                prm.append(p)
                bounds.append(b)
            else:
                prm.append(
                    self._get_prm_func_dict[self._prm_list[k][1]](
                        self._prm_list[k][0], 0))
                bounds.append((-self._b_2_max, self._b_2_max))
        return np.array(prm), bounds

    def set_prm(self, lat_prop, prm):
        prm_ind = 0
        for k in range(len(self._prm_list)):
            if self._prm_list[k][0] == "bend":
                prm_ind = \
                    self._prm_list[k][1].set_bend_prm(lat_prop, prm, prm_ind)
            elif self._prm_list[k][0] == "rev_bend":
                prm_ind = \
                    self._prm_list[k][1].set_rev_bend_prm(
                        lat_prop, prm, prm_ind)
            else:
                self._set_prm_func_dict[self._prm_list[k][1]](
                    self._prm_list[k][0], prm[prm_ind])
                prm_ind += 1
        
    def prt_prm(self, prm):
        n_prt = 5
        print("\n  prm =")
        for k in range(len(prm)):
            print("  " if k == 0 else "",
                  "{:12.5e}".format(prm[k]),
                  "\n  " if (k != 0) and (k % n_prt == n_prt-1) else "",
                  end="")
        if k % n_prt != n_prt-1:
            print()


def prt_lat(lat_prop, file_name, prm_list):
    outf = open(file_name, 'w')

    def prt_drift(name):
        L = lat_prop.get_L_elem(name, 0)
        print(("{:5s}: Drift, L = {:7.5f};").format(name, L), file=outf)

    def prt_dip(name):
        L = lat_prop.get_L_elem(name, 0)
        phi = lat_prop.get_phi_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:5s}: Bending, L = {:7.5f}, T = {:8.5f}, K = {:8.5f}"
               ", T1 = 0.0, T2 = 0.0,\n       N = n_bend;")
              .format(name, L, phi, b_2), file=outf)

    def prt_bend(bend):
        prt_dip(bend._dip_0)
        for k in range(len(bend._dip_list)):
            prt_dip(bend._dip_list[k])

    def prt_rev_bend(rev_bend):
        for k in range(len(rev_bend._rb_list)):
            prt_dip(rev_bend._rb_list[k])
        for k in range(len(rev_bend._dip_list)):
            prt_dip(rev_bend._dip_list[k])

    def prt_quad(name):
        L = lat_prop.get_L_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:5s}: Quadrupole, L = {:7.5f}, K = {:8.5f}, N = n_quad;")
              .format(name, L, b_2), file=outf)

    # Dictionary of parameter types and corresponding print functions.
    get_prm_func_dict = {
        "L":   prt_drift,
        "L_b": prt_dip,
        "phi": prt_dip,
        "b_2": prt_dip
    }

    for k in range(len(prm_list._prm_list)):
        if prm_list._prm_list[k][0] == "bend":
             prt_bend(prm_list._prm_list[k][1])
        elif prm_list._prm_list[k][0] == "rev_bend":
             prt_rev_bend(prm_list._prm_list[k][1])
        else:
            get_prm_func_dict[prm_list._prm_list[k][1]](
                prm_list._prm_list[k][0])


__all__ = [prm_class, prt_lat]
