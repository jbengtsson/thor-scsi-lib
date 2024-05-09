"""
"""

import numpy as np


quad = 2


class opt_phi_class:
    # Private

    def __init__(self, lat_prop, bend_tweak, bend_list, phi_max):
        self._bend_tweak = bend_tweak
        self._bend_list  = bend_list
        self._phi_tot    = self.compute_phi_bend(lat_prop)
        self._phi_tweak  = lat_prop.get_phi_elem(self._bend_tweak, 0)

        self._phi_max    = phi_max

    # Public.

    def compute_phi_bend(self, lat_prop):
        # Compute total bend angle.
        phi = 0e0
        for k in range(len(self._bend_list)):
            n = lat_prop.get_n_kid(self._bend_list[k])
            phi += n*lat_prop.get_phi_elem(self._bend_list[k], 0)
        n = lat_prop.get_n_kid(self._bend_tweak)
        phi += n*lat_prop.get_phi_elem(self._bend_tweak, 0)
        return phi

    def get_bend_prm(self, lat_prop):
        prt = not False
        prm = []
        bounds = []
        if prt:
            print("\nget_bend_prm:")
        for k in range(len(self._bend_list)):
            if prt:
                print("  {:s}".format(self._bend_list[k]))
            phi = lat_prop.get_phi_elem(self._bend_list[k], 0)
            prm.append(phi)
            bounds.append((-self._phi_max, self._phi_max))
        return prm, bounds

    def set_bend_prm(self, lat_prop, prm, prm_ind):
        for k in range(len(self._bend_list)):
            lat_prop.set_phi_fam(self._bend_list[k], prm[prm_ind])
            prm_ind += 1
        phi = self.compute_phi_bend(lat_prop)
        n = lat_prop.get_n_kid(self._bend_tweak)
        self._phi_tweak += (self._phi_tot - phi) / n
        lat_prop.set_phi_fam(self._bend_tweak, self._phi_tweak)
        return prm_ind


class phi_bend_class:
    # Private

    def __init__(self, lat_prop, bend_tweak, bend_list, phi_max):
        self._bend_tweak = bend_tweak
        self._bend_list  = bend_list
        self._phi_tot    = self.compute_phi_bend(lat_prop)
        self._phi_tweak  = lat_prop.get_phi_elem(self._bend_tweak, 0)
        self._phi_max    = phi_max

    # Public.

    def compute_phi_bend(self, lat_prop):
        # Compute total bend angle.
        phi = 0e0
        for k in range(len(self._bend_list)):
            n = lat_prop.get_n_kid(self._bend_list[k])
            phi += n*lat_prop.get_phi_elem(self._bend_list[k], 0)
        n = lat_prop.get_n_kid(self._bend_tweak)
        phi += n*lat_prop.get_phi_elem(self._bend_tweak, 0)
        return phi

    def get_bend_prm(self, lat_prop):
        prt = not False
        prm = []
        bounds = []
        if prt:
            print("\nget_bend_prm:")
        for k in range(len(self._bend_list)):
            if prt:
                print("  {:s}".format(self._bend_list[k]))
            phi = lat_prop.get_phi_elem(self._bend_list[k], 0)
            prm.append(phi)
            bounds.append((-self._phi_max, self._phi_max))
        return prm, bounds

    def set_bend_prm(self, lat_prop, prm, prm_ind):
        for k in range(len(self._bend_list)):
            lat_prop.set_phi_fam(self._bend_list[k], prm[prm_ind])
            prm_ind += 1
        phi = self.compute_phi_bend(lat_prop)
        n = lat_prop.get_n_kid(self._bend_tweak)
        self._phi_tweak += (self._phi_tot - phi) / n
        lat_prop.set_phi_fam(self._bend_tweak, self._phi_tweak)
        return prm_ind


class phi_bend_class:
    # Private

    def __init__(self, lat_prop, bend_list, phi_max):
        self._bend_list      = bend_list
        self._bend_L         = np.nan
        self._bend_phi       = np.nan
        self._bend_phi_ratio = np.nan
        self._phi_max        = phi_max

        self.compute_bend_L(lat_prop)
        self.compute_bend_phi(lat_prop)
        self.compute_phi_ratios(lat_prop)

    # Public.

    def compute_bend_L(self, lat_prop):
        # Compute total length.
        self._bend_L = 0e0
        for k in range(len(self._bend_list)):
            self._bend_L += lat_prop.get_L_elem(self._bend_list[k], 0)

    def compute_bend_phi(self, lat_prop):
        # Compute total bend angle.
        self._bend_phi = 0e0
        for k in range(len(self._bend_list)):
            print("\ncompute_bend_phi:", self._bend_list[k])
            self._bend_phi += lat_prop.get_phi_elem(self._bend_list[k], 0)

    def compute_phi_ratios(self, lat_prop):
        # Compute phi ratios.
        self._bend_phi_ratio = []
        for k in range(len(self._bend_list)):
            phi = lat_prop.get_phi_elem(self._bend_list[k], 0)
            self._bend_phi_ratio.append(phi/self._bend_phi)
        self._bend_phi_ratio = np.array(self._bend_phi_ratio)

    def get_bend_phi_prm(self):
        prm = self._bend_phi
        bound = (-self._phi_max, self._phi_max)
        return prm, bound

    def set_bend_phi(self, lat_prop, prm):
        prt = False
        self._bend_phixL = prm
        if prt:
            print("\nset_bend_phi\n  prm = {:7.5f}\n".format(prm))
        for k in range(len(self._bend_list)):
            if prt:
                print("  {:5s} {:8.5f}".
                      format(self._bend_list[k], self._bend_phi_ratio[k]*prm))
            lat_prop.set_phi_fam(
                self._bend_list[k], self._bend_phi_ratio[k]*prm)

    def print(self):
        print("  L   = {:8.5f}".format(self._bend_L))
        print("  phi = {:8.5f}\n".format(self._bend_phi))
        for k in range(len(self._bend_list)):
            print("  {:5s} {:8.5f}".
                  format(self._bend_list[k], self._bend_phi_ratio[k],
                         self._bend_phi*self._bend_phi_ratio[k]))


class b_2_bend_class:
    # Private

    def __init__(self, lat_prop, bend_list, b_2_max):
        self._bend_list        = bend_list
        self._bend_L           = np.nan
        self._bend_b_2xL       = np.nan
        self._bend_b_2xL_ratio = np.nan

        self._b_2_max          = b_2_max

        self.compute_bend_L(lat_prop)
        self.compute_bend_b_2xL(lat_prop)
        self.compute_b_2xL_ratios(lat_prop)

    # Public.

    def compute_bend_L(self, lat_prop):
        # Compute total length.
        self._bend_L = 0e0
        for k in range(len(self._bend_list)):
            self._bend_L += lat_prop.get_L_elem(self._bend_list[k], 0)

    def compute_bend_b_2xL(self, lat_prop):
        # Compute total b_2xL.
        self._bend_b_2xL = 0e0
        for k in range(len(self._bend_list)):
            print("\ncompute_bend_b_2xL:", self._bend_list[k])
            self._bend_b_2xL += \
                lat_prop.get_b_nxL_elem(self._bend_list[k], 0, quad)

    def compute_b_2xL_ratios(self, lat_prop):
        # Compute b_2xL ratios.
        self._bend_b_2xL_ratio = []
        for k in range(len(self._bend_list)):
            b_2xL = lat_prop.get_b_nxL_elem(self._bend_list[k], 0, quad)
            self._bend_b_2xL_ratio.append(b_2xL/self._bend_b_2xL)
        self._bend_b_2xL_ratio = np.array(self._bend_b_2xL_ratio)

    def get_bend_b_2_prm(self):
        prm = self._bend_b_2xL
        bound = (-self._b_2_max, self._b_2_max)
        return prm, bound

    def set_bend_b_2(self, lat_prop, prm):
        prt = False
        self._bend_b_2xL = prm
        if prt:
            print("\nset_bend_b_2\n  prm = {:7.5f}\n".format(prm))
        for k in range(len(self._bend_list)):
            if prt:
                print("  {:5s} {:8.5f}".
                      format(self._bend_list[k], self._bend_b_2xL_ratio[k]*prm))
            lat_prop.set_b_n_fam(
                self._bend_list[k], quad, self._bend_b_2xL_ratio[k]*prm)

    def print(self):
        print("  L   = {:8.5f}".format(self._bend_L))
        print("  b_2 = {:8.5f}\n".format(self._bend_b_2xL/self._bend_L))
        for k in range(len(self._bend_list)):
            print("  {:5s} {:8.5f}".
                  format(self._bend_list[k], self._bend_b_2xL_ratio[k],
                         self._bend_b_2xL*self._bend_b_2xL_ratio[k]))


class prm_class(phi_bend_class, b_2_bend_class):
    # Private

    def __init__(self, lat_prop, prm_list, b_2_max):
        self._b_2_max  = b_2_max
        self._prm_list = prm_list

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

    # Public.

    def get_prm(self):
        prt = not False
        prm = []
        bounds = []
        if prt:
            print("\nget_prm:")
        for k in range(len(self._prm_list)):
            if self._prm_list[k][0] == "phi_bend":
                if prt:
                    print("  {:5s}".format(self._prm_list[k][0]))
                p, b = self._prm_list[k][1].get_bend_phi_prm()
                self._prm_list[k][1].print()
                prm.append(p)
                bounds.append(b)
            elif self._prm_list[k][0] == "b_2_bend":
                if prt:
                    print("  {:5s}".format(self._prm_list[k][0]))
                p, b = self._prm_list[k][1].get_bend_b_2_prm()
                self._prm_list[k][1].print()
                prm.append(p)
                bounds.append(b)
            else:
                if prt:
                    print("  {:5s}".format(self._prm_list[k][0]))
                prm.append(
                    self._get_prm_func_dict[self._prm_list[k][1]](
                        self._prm_list[k][0], 0))
                bounds.append((-self._b_2_max, self._b_2_max))
        return np.array(prm), bounds

    def set_prm(self, lat_prop, prm):
        for k in range(len(self._prm_list)):
            if self._prm_list[k][0] == "phi_bend":
                self._prm_list[k][1].set_bend_phi(lat_prop, prm[k])
            elif self._prm_list[k][0] == "b_2_bend":
                self._prm_list[k][1].set_bend_b_2(lat_prop, prm[k])
            else:
                self._set_prm_func_dict[self._prm_list[k][1]](
                    self._prm_list[k][0], prm[k])
        
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
        for k in range(len(bend._bend_list)):
            prt_dip(bend._bend_list[k])

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
        if prm_list._prm_list[k][0] == "phi_bend":
             prt_bend(prm_list._prm_list[k][1])
        elif prm_list._prm_list[k][0] == "b_2_bend":
             prt_bend(prm_list._prm_list[k][1])
        else:
            get_prm_func_dict[prm_list._prm_list[k][1]](
                prm_list._prm_list[k][0])


__all__ = [prm_class, prt_lat]
