"""
"""

import numpy as np


quad = 2


class bend_class:
    # Compound bending magnet.

    # Private

    def __init__(self, lat_prop, bend_list, phi_max, b_2_max):
        self._lat_prop         = lat_prop
        self._bend_list        = bend_list

        self._bend_L           = np.nan
        self._bend_phi_ratio   = np.nan
        self._bend_b_2xL_ratio = np.nan

        self._bend_phi         = np.nan
        self._bend_b_2xL       = np.nan

        self._phi_max          = phi_max
        self._b_2_max          = b_2_max

        self.compute_bend_L()
        self.compute_bend_phi()
        self.compute_bend_b_2xL()
        self.compute_phi_ratios()
        self.compute_b_2xL_ratios()

    # Public.

    def compute_bend_L(self):
        # Compute total length.
        self._bend_L = 0e0
        for k in range(len(self._bend_list)):
            self._bend_L += self._lat_prop.get_L_elem(self._bend_list[k], 0)

    def compute_bend_phi(self):
        # Compute total bend angle.
        self._bend_phi = 0e0
        for k in range(len(self._bend_list)):
            self._bend_phi += self._lat_prop.get_phi_elem(self._bend_list[k], 0)

    def compute_bend_b_2xL(self):
        # Compute total b_2xL.
        self._bend_b_2xL = 0e0
        for k in range(len(self._bend_list)):
            self._bend_b_2xL += \
                self._lat_prop.get_b_nxL_elem(self._bend_list[k], 0, quad)

    def compute_phi_ratios(self):
        # Compute phi ratios.
        self._bend_phi_ratio = []
        for k in range(len(self._bend_list)):
            phi = self._lat_prop.get_phi_elem(self._bend_list[k], 0)
            self._bend_phi_ratio.append(phi/self._bend_phi)
        self._bend_phi_ratio = np.array(self._bend_phi_ratio)

    def compute_b_2xL_ratios(self):
        # Compute b_2xL ratios.
        self._bend_b_2xL_ratio = []
        for k in range(len(self._bend_list)):
            b_2xL = self._lat_prop.get_b_nxL_elem(self._bend_list[k], 0, quad)
            self._bend_b_2xL_ratio.append(b_2xL/self._bend_b_2xL)
        self._bend_b_2xL_ratio = np.array(self._bend_b_2xL_ratio)

    def get_bend_phi(self):
        return self._bend_phi

    def get_bend_phi_prm(self):
        prm = self._bend_phi
        bound = (-self._phi_max, self._phi_max)
        return prm, bound

    def get_bend_b_2_prm(self):
        prm = self._bend_b_2xL
        bound = (-self._b_2_max, self._b_2_max)
        return prm, bound

    def set_bend_phi(self, prm):
        prt = False
        self._bend_phixL = prm
        if prt:
            print("\nset_bend_phi\n  prm = {:7.5f}\n".format(prm))
        for k in range(len(self._bend_list)):
            if prt:
                print("  {:5s} {:8.5f}".
                      format(self._bend_list[k], self._bend_phi_ratio[k]*prm))
            self._lat_prop.set_phi_fam(
                self._bend_list[k], self._bend_phi_ratio[k]*prm)

    def set_bend_b_2(self, prm):
        prt = False
        self._bend_b_2xL = prm
        if prt:
            print("\nset_bend_b_2\n  prm = {:7.5f}\n".format(prm))
        for k in range(len(self._bend_list)):
            if prt:
                print("  {:5s} {:8.5f}".
                      format(self._bend_list[k], self._bend_b_2xL_ratio[k]*prm))
            self._lat_prop.set_b_n_fam(
                self._bend_list[k], quad, self._bend_b_2xL_ratio[k]*prm)

    def print(self):
        print("\n  bend_class:")
        print("    L   = {:8.5f}".format(self._bend_L))
        print("    phi = {:8.5f}".format(self._bend_phi))
        print("    b_2 = {:8.5f}".format(self._bend_b_2xL/self._bend_L))
        print("\n    part   fraction   phi    fraction   b_2")
        for k in range(len(self._bend_list)):
            print("    {:5s} {:8.5f} {:8.5f} {:8.5f} {:8.5f}".
                  format(self._bend_list[k], self._bend_phi_ratio[k],
                         self._bend_phi*self._bend_phi_ratio[k],
                         self._bend_b_2xL_ratio[k],
                         self._bend_b_2xL*self._bend_b_2xL_ratio[k]))


class phi_tot_class():
    # Private

    def __init__(self, lat_prop, phi_list, phi_n_list, phi_max):
        self._lat_prop   = lat_prop
        self._phi_list   = phi_list
        self._phi_n_list = phi_n_list
        self._phi        = np.nan
        self._phi_max    = phi_max

        self.compute_phi()

    # Public.

    def compute_phi(self):
        # Compute total bend angle.
        self._phi = 0e0
        for k in range(len(self._phi_list)):
            if type(self._phi_list[k]) == str:
                self._phi += \
                    self._phi_n_list[k]*self.get_phi_elem(self._phi_list[k], 0)
            elif isinstance(self._phi_list[k], bend_class):
                self._phi += \
                    self._phi_n_list[k]*self._phi_list[k].get_bend_phi()
            else:
                print("\nget_prm - undefined element type:", self._phi_list[k])

    def get_phi_prm(self):
        prt = not False
        prm = []
        bounds = []
        if prt:
            self.print()
        p, b = self._phi_list[0].get_bend_phi_prm()
        return p, b

    def print(self):
        print("  phi_tot_class:")
        print("    phi = {:8.5f}".format(self._phi))
        print("    phi_list:")
        for k in range(len(self._phi_list)):
            if self._phi_list[k] == str:
                print("      {:s}".format(self._phi_list[k]))
            elif isinstance(self._phi_list[k], bend_class):
                print("      bend_class")
            else:
                print("\nphi_tot_class::print() - undefined parameter type.")
                assert False


class prm_class(bend_class):
    # Private

    def __init__(self, lat_prop, prm_list, b_2_max):
        self._lat_prop = lat_prop
        self._b_2_max  = b_2_max
        self._prm_list = prm_list

        # Dictionary of parameter types and corresponding get functions.
        self._get_prm_func_dict = {
            "L":   self._lat_prop.get_L_elem,
            "L_b": self._lat_prop.get_L_elem,
            "phi": self._lat_prop.get_phi_elem,
            "b_2": self._lat_prop.get_b_2_elem
        }

        # Dictionary of parameter types and corresponding set functions.
        self._set_prm_func_dict = {
            "L":   self._lat_prop.set_L_fam,
            "L_b": self._lat_prop.set_L_bend_fam,
            "phi": self._lat_prop.set_phi_fam,
            "b_2": self._lat_prop.set_b_2_fam
        }

    # Public.

    def get_prm(self):
        prt = not False
        prm = []
        bounds = []
        if prt:
            print("\nprm_class::get_prm:")
        for k in range(len(self._prm_list)):
            if prt:
                print("\n  {:s}".format(self._prm_list[k][0]))
            if type(self._prm_list[k][1]) == str:
                p = \
                    self._get_prm_func_dict[self._prm_list[k][1]](
                        self._prm_list[k][0], 0)
                b = (-self._b_2_max, self._b_2_max)
            elif isinstance(self._prm_list[k][1], bend_class):
                if self._prm_list[k][0] == "phi_bend":
                    p, b = self._prm_list[k][1].get_bend_phi_prm()
                elif self._prm_list[k][0] == "b_2_bend":
                    p, b = self._prm_list[k][1].get_bend_b_2_prm()
                else:
                    print("\nprm_class::get_prm - undefined parameter type:",
                          self._prm_list[k][1])
                    assert False
            elif isinstance(self._prm_list[k][1], phi_tot_class):
               p, b = self._prm_list[k][1].get_phi_prm()
            else:
                print("\nprm_class::get_prm - undefined parameter type:",
                      self._prm_list[k][1])
                assert False
            prm.append(p)
            bounds.append(b)
            print("  prm    = ", prm)
            print("  bounds = ", bounds)
        return np.array(prm), bounds

    def set_prm(self, prm):
        for k in range(len(self._prm_list)):
            prm_name, prm_type = [self._prm_list[k][0], self._prm_list[k][1]]
            if type(prm_type) == str:
                self._set_prm_func_dict[prm_type](prm_name, prm[k])
            elif isinstance(prm_type, bend_class):
                if prm_name == "phi_bend":
                    prm_type.set_bend_phi(prm[k])
                elif prm_name == "b_2_bend":
                    prm_type.set_bend_b_2(prm[k])
                else:
                    print("\nprm_class::set_prm - undefined parameter type:",
                          prm_type)
            elif self._prm_list[k][0] == "phi_tot":
                # The 1st element in the list is the bend angle paramameter
                # whereas the 2nd one is used for maintaining the total bend
                # angle.
                prm_type._phi_list[0].set_bend_phi(prm[k])
                phi = \
                    (prm_type._phi-prm_type._phi_n_list[0]*prm[k]) \
                    /prm_type._phi_n_list[1]
                prm_type._phi_list[1].set_bend_phi(phi)
            else:
                print("\nprm_class::set_prm - undefined parameter:", prm_type)
        
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
        if type(prm_list._prm_list[k][1]) == str:
            get_prm_func_dict[prm_list._prm_list[k][1]](
                prm_list._prm_list[k][0])
        elif isinstance(prm_list._prm_list[k][1], bend_class):
            prt_bend(prm_list._prm_list[k][1])
        elif isinstance(prm_list._prm_list[k][1], phi_tot_class):
            # The 1st element in the list is the bend angle paramameter whereas
            # the 2nd one is used for maintaining the total bend angle.
            prt_bend(prm_list._prm_list[k][1]._phi_list[0])
            prt_bend(prm_list._prm_list[k][1]._phi_list[1])
        else:
            print("\nprt_lat - undefined parameter:", prm_list._prm_list[k][1])


__all__ = [prm_class, prt_lat]
