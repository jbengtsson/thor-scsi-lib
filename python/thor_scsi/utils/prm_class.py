"""
"""

import numpy as np


import thor_scsi.lib as ts


quad = 2


class bend_class:
    # Compound bending magnets.

    # Private

    def __init__(self, lat_prop, bend_list):
        self._lat_prop         = lat_prop
        self._bend_list        = bend_list

        self._bend_L_tot       = self.compute_bend_L_tot()

        self._bend_phi         = self.compute_bend_phi()
        self._bend_b_2xL       = self.compute_bend_b_2xL()

        self._bend_phi_ratio   = self.compute_phi_ratios()
        self._bend_b_2xL_ratio = self.compute_b_2_ratios()

        self._phi_min          = np.nan
        self._phi_max          = np.nan
        self._b_2_min          = np.nan
        self._b_2_max          = np.nan

    def __str__(self):
        return "bend_class"

    # Public.

    def compute_bend_L_tot(self):
        # Compute total length.
        L = 0e0
        for k in range(len(self._bend_list)):
            try:
                elem_name = self._bend_list[k]
                dL = self._lat_prop.get_L_elem(elem_name, 0)
            except:
                print(f"\n*** compute_bend_phi - undefined element:"
                      f" {elem_name:s}")
                exit(1)
            else:
                L += dL
        return L

    def compute_bend_phi(self):
        # Compute total bend angle.
        phi = 0e0
        for k in range(len(self._bend_list)):
            phi += self._lat_prop.get_phi_elem(self._bend_list[k], 0)
        return phi

    def compute_bend_b_2xL(self):
        # Compute total b_2xL.
        b_2xL = 0e0
        for k in range(len(self._bend_list)):
            b_2xL += self._lat_prop.get_b_nxL_elem(self._bend_list[k], 0, quad)
        return b_2xL

    def compute_phi_ratios(self):
        # Compute phi ratios.
        self._bend_phi_ratio = []
        for k in range(len(self._bend_list)):
            phi = self._lat_prop.get_phi_elem(self._bend_list[k], 0)
            self._bend_phi_ratio.append(phi/self._bend_phi)
        return np.array(self._bend_phi_ratio)

    def compute_b_2_ratios(self):
        # Compute b_2 ratios.
        self._bend_b_2xL_ratio = []
        for k in range(len(self._bend_list)):
           b_2xL = self._lat_prop.get_b_nxL_elem(self._bend_list[k], 0, quad)
           self._bend_b_2xL_ratio.append(b_2xL/self._bend_b_2xL)
        return np.array(self._bend_b_2xL_ratio)

    def get_bend_phi(self):
        return self._bend_phi

    def get_bend_phi_prm(self):
        prm = self._bend_phi
        bound = (self._phi_min, self._phi_max)
        return prm, bound

    def get_bend_b_2_prm(self):
        prm = self._bend_b_2xL/self._bend_L_tot
        bound = (self._b_2_min, self._b_2_max)
        return prm, bound

    def set_bend_dphi(self, prm):
        self._bend_phi = prm
        n = len(self._bend_list)
        dphi = prm/n
        for k in range(n):
            self._lat_prop.set_dphi_fam(self._bend_list[k], dphi)

    def set_bend_phi(self, prm):
        prt = False
        self._bend_phi = prm
        if prt:
            print("\nset_bend_phi\n  prm = {:7.5f}\n".format(prm))
        for k in range(len(self._bend_list)):
            self._lat_prop.set_phi_fam(
                self._bend_list[k], self._bend_phi_ratio[k]*prm)

    def correct_bend_phi(self):
        # Maintain total bend angle.
        prt = False
        if prt:
            print("\ncorrect_bend_phi:")
        dphi = self._bend_phi - self.compute_bend_phi()
        n = len(self._bend_list)
        dphi /= n
        for k in range(n):
            if prt:
                print("  {:5s}".format(self._bend_list[k]))
            self._lat_prop.set_dphi_fam(self._bend_list[k], dphi)
        if prt:
            print("  {:9.3e} ({:9.3e})".format(
                self.compute_bend_phi(), self._bend_phi))

    def set_bend_b_2(self, prm):
        prt = False
        self._bend_b_2xL = prm*self._bend_L_tot
        if prt:
            print("\nset_bend_b_2\n  prm = {:7.5f} bend_b_2xL = {:7.5f}\n".
                  format(prm, self._bend_b_2xL))
        for k in range(len(self._bend_list)):
            L_k = self._lat_prop.get_L_elem(self._bend_list[k], 0)
            b_2 = self._bend_b_2xL_ratio[k]*self._bend_b_2xL/L_k
            if prt:
                print("  {:15s} {:8.5f}".format(self._bend_list[k], b_2))
            self._lat_prop.set_b_n_fam(self._bend_list[k], quad, b_2)

    def print(self):
        prt_debug = False
        print(" L = {:8.5f} phi = {:8.5f} b_2 = {:8.5f}".
              format(self._bend_L_tot, self._bend_phi,
                     self._bend_b_2xL/self._bend_L_tot))
        print("\n        part            L    phi_k/phi   phi"
              "   (b_2xL)_k/b_2xL  b_2 comp.    b_2")
        if prt_debug:
            phi_ratio = 0e0
            b_2_ratio = 0e0
        for k in range(len(self._bend_list)):
            if prt_debug:
                phi_ratio += self._bend_phi_ratio[k]
                b_2_ratio += self._bend_b_2xL_ratio[k]
            L_k = self._lat_prop.get_L_elem(self._bend_list[k], 0)
            print("    {:15s} {:8.5f} {:8.5f} {:8.5f}    {:8.5f}      {:8.5f}"
                  "  {:8.5f}".format(
                      self._bend_list[k],
                      self._lat_prop.get_L_elem(self._bend_list[k], 0),
                      self._bend_phi_ratio[k],
                      self._bend_phi*self._bend_phi_ratio[k],
                      self._bend_b_2xL_ratio[k],
                      self._bend_b_2xL*self._bend_b_2xL_ratio[k]/L_k,
                      self._lat_prop.get_b_n_elem(self._bend_list[k], 0, quad)))
        if prt_debug:
            print("\n  Sum{{phi_k/phi}} = {:5.3f}"
                  " Sum{{(b_2xL)_k/b_2xL}} = {:5.3f}".format(
                      phi_ratio, b_2_ratio))


class phi_lat_class():
    # Class to maintain the total bend angle for a super period.
    # Private

    def __init__(self, lat_prop, n_bend, dip_prm):
        self._lat_prop = lat_prop
        self._dip_prm  = dip_prm
        self._n_dip    = n_bend
        self._phi_lat  = lat_prop.compute_phi_lat()

        print("\nphi_lat_class: {:5.3f}".format(self._phi_lat), self._dip_prm,
              "{:d}".format(self._n_dip))

    # Public.

    def set_phi_lat(self):
        prt = False
        phi_lat_0 = self._phi_lat
        phi_lat = self._lat_prop.compute_phi_lat()
        dphi = (phi_lat-phi_lat_0)/self._n_dip
        if prt:
            print("\nset_phi_lat:\n  phi_lat     = {:5.3f}"
                  "\n  phi_lat_ref = {:5.3f}\n  dphi_lat    = {:5.3f}".
                  format(phi_lat, phi_lat_0, dphi))
        if type(self._dip_prm) == str:
            self._lat_prop.set_dphi_fam(self._dip_prm, -dphi)
        elif isinstance(self._dip_prm, bend_class):
            if prt:
                print("  set_bend_dphi(-dphi)")
            self._dip_prm.set_bend_dphi(-dphi)
            phi_lat = self._lat_prop.compute_phi_lat()
            dphi = (phi_lat-phi_lat_0)/self._n_dip
            if prt:
                print("\n  phi_lat     = {:5.3f}\n  dphi_lat    = {:5.3f}".
                      format(phi_lat, dphi))
        else:
            print("\nset_phi_lat:\n undefined type: ", type(self._dip_prm))
            assert False

class phi_list_class():
    # Private

    def __init__(self, dip_list):
        self._dip_list = dip_list
        self._phi_bend = compute_phi_list()

        print("\nphi_lat_class: {:5.3f}".format(self._phi))

    # Public.

    def set_phi_lat(self):
        prt = False
        phi_lat_0 = self._phi_lat
        phi_lat = self._lat_prop.compute_phi_lat()
        dphi = (phi_lat-phi_lat_0)/self._n_dip
        if prt:
            print("\nset_phi_lat:\n  phi_lat     = {:5.3f}"
                  "\n  phi_lat_ref = {:5.3f}\n  dphi_lat    = {:5.3f}".
                  format(phi_lat, phi_lat_0, dphi))
        if type(self._dip_prm) == str:
            self._lat_prop.set_dphi_fam(self._dip_prm, -dphi)
        elif isinstance(self._dip_prm, bend_class):
            if prt:
                print("  set_bend_dphi(-dphi)")
            self._dip_prm.set_bend_dphi(-dphi)
            phi_lat = self._lat_prop.compute_phi_lat()
            dphi = (phi_lat-phi_lat_0)/self._n_dip
            if prt:
                print("\n  phi_lat     = {:5.3f}\n  dphi_lat    = {:5.3f}".
                      format(phi_lat, dphi))
        else:
            print("\nset_phi_lat:\n undefined type: ", type(self._dip_prm))
            assert False

class prm_class(bend_class):
    # Private

    def __init__(self, lat_prop, prm_list):
        self._lat_prop = lat_prop
        self._prm_list = prm_list

        # Dictionary of parameter types and corresponding get functions.
        self._get_prm_func_dict = {
            "L":   self._lat_prop.get_L_elem,
            "L_b": self._lat_prop.get_L_elem,
            "phi": self._lat_prop.get_phi_elem,
            "h":   self._lat_prop.get_h_elem,
            "b_2": self._lat_prop.get_b_2_elem,
            "b_3": self._lat_prop.get_b_3_elem,
            "b_4": self._lat_prop.get_b_4_elem}

        # Dictionary of parameter types and corresponding set functions.
        self._set_prm_func_dict = {
            "L":   self._lat_prop.set_L_fam,
            "L_b": self._lat_prop.set_L_bend_fam,
            "phi": self._lat_prop.set_phi_fam,
            "h":   self._lat_prop.set_h_fam,
            "b_2": self._lat_prop.set_b_2_fam,
            "b_3": self._lat_prop.set_b_3_fam,
            "b_4": self._lat_prop.set_b_4_fam}

    # Public.

    def get_prm(self):
        prt = not False
        prm = []
        bounds = []
        if prt:
            print("\nprm_class::get_prm:")
        print("----------------------------------------")
        for k in range(len(self._prm_list)):
            if prt:
                print("  {:10s}".format(self._prm_list[k][1]), end="")
            if type(self._prm_list[k][0]) == str:
                if prt:
                    print(" {:10s}".format(self._prm_list[k][0]), end="")
                p = \
                    self._get_prm_func_dict[self._prm_list[k][1]](
                        self._prm_list[k][0], 0)
                b = (self._prm_list[k][2], self._prm_list[k][3])
                prm.append(p)
                bounds.append(b)
                if prt:
                    print("    {:12.5e} ({:5.2f}, {:5.2f})".
                      format(prm[k], bounds[k][0], bounds[k][1]))
            elif isinstance(self._prm_list[k][0], bend_class):
                if self._prm_list[k][1] == "phi_bend":
                    self._prm_list[k][0]._phi_min = self._prm_list[k][2]
                    self._prm_list[k][0]._phi_max = self._prm_list[k][3]
                    p, b = self._prm_list[k][0].get_bend_phi_prm()
                elif self._prm_list[k][1] == "b_2_bend":
                    self._prm_list[k][0]._b_2_min = self._prm_list[k][2]
                    self._prm_list[k][0]._b_2_max = self._prm_list[k][3]
                    p, b = self._prm_list[k][0].get_bend_b_2_prm()
                else:
                    print("\nprm_class::get_prm - undefined parameter type:",
                          self._prm_list[k][1])
                    assert False
                prm.append(p)
                bounds.append(b)
                if prt:
                    print("", self._prm_list[k][0],
                          "   {:12.5e} ({:5.2f}, {:5.2f})".format(
                              prm[k], bounds[k][0], bounds[k][1]), end="")
                    self._prm_list[k][0].print()
            else:
                print("\nprm_class::get_prm - undefined parameter type:",
                      self._prm_list[k][1])
                assert False
            print("----------------------------------------")
        return np.array(prm), bounds

    def set_prm(self, prm):
        for k in range(len(self._prm_list)):
            prm_name, prm_type = [self._prm_list[k][0], self._prm_list[k][1]]
            if type(prm_name) == str:
                self._set_prm_func_dict[prm_type](prm_name, prm[k])
            elif isinstance(prm_name, bend_class):
                if prm_type == "phi_bend":
                    prm_name.set_bend_phi(prm[k])
                elif prm_type == "b_2_bend":
                    prm_name.set_bend_b_2(prm[k])
                else:
                    print("\nprm_class::set_prm - undefined parameter type:",
                          prm_name)
            else:
                print("\nprm_class::set_prm - undefined parameter:", prm_name)
        
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


def prt_lat(
        lat_prop,
        file_name,
        prm_list,
        *,
        bend_list = None):
    outf = open(file_name, 'w')

    def prt_drift(fam_name, prt_list):
        if prt_list.count(fam_name) == 0:
            prt_list.append(fam_name)
            L = lat_prop.get_L_elem(fam_name, 0)
            print("{:5s}: Drift, L = {:7.5f};".format(fam_name, L), file=outf)
        return prt_list

    def prt_dip(fam_name, prt_list):
        if prt_list.count(fam_name) == 0:
            prt_list.append(fam_name)
            L = lat_prop.get_L_elem(fam_name, 0)
            phi = lat_prop.get_phi_elem(fam_name, 0)
            phi_1 = lat_prop.get_phi_1_elem(fam_name, 0)
            phi_2 = lat_prop.get_phi_1_elem(fam_name, 0)
            b_2 = lat_prop.get_b_n_elem(fam_name, 0, quad)
            print("{:5s}: Bending, L = {:7.5f}".format(fam_name, L), end="",
                  file=outf)
            print(", Phi = {:8.5f}".format(phi), end="", file=outf)
            print(", Phi_1 = {:8.5f}".format(phi_1), end="", file=outf)
            print(",\n    Phi_2 = {:8.5f}".format(phi_2), end="", file=outf)
            if b_2 != 0e0:
                print(", B_2 = {:8.5f}".format(b_2), end="", file=outf)
            print(", N = n_bend;", file=outf)
        return prt_list

    def prt_bend(bend, prt_list):
        for k in range(len(bend._bend_list)):
            prt_list = prt_dip(bend._bend_list[k], prt_list)
        return prt_list

    def prt_quad(fam_name, prt_list):
        if prt_list.count(fam_name) == 0:
            k = lat_prop._lattice.find(fam_name, 0).index
            if type(lat_prop._lattice[k]) == ts.Bending:
                prt_list = prt_dip(fam_name, prt_list)
            else:
                prt_list.append(fam_name)
                L = lat_prop.get_L_elem(fam_name, 0)
                b_2 = lat_prop.get_b_n_elem(fam_name, 0, quad)
                phi = lat_prop.get_phi_elem(fam_name, 0)
                if phi == 0e0:
                    print("{:5s}: Quadrupole, L = {:7.5f}, B_2 = {:8.5f}"
                          ", N = n_quad;"
                          .format(fam_name, L, b_2), file=outf)
                else:
                    print("{:5s}: Bending, L = {:7.5f}, Phi = {:8.5f}"
                          ", B_2 = {:8.5f}, N = n_bend;"
                          .format(fam_name, L, phi, b_2), file=outf)
        return prt_list

    def prt_sext(fam_name, prt_list):
        if prt_list.count(fam_name) == 0:
            prt_list.append(fam_name)
            k = lat_prop._lattice.find(fam_name, 0).index
            L = lat_prop._lattice[k].get_length()
            b_3 = lat_prop._lattice[k].get_multipoles().get_multipole(3).real
            if type(lat_prop._lattice[k]) == ts.Sextupole:
                print("{:5s}: Sextupole, L = {:7.5f}, B_3 = {:10.5f}"
                      ", N = n_sext;"
                      .format(fam_name, L, b_3), file=outf)
            elif type(lat_prop._lattice[k]) == ts.Multipole:
                b_2 = \
                    lat_prop._lattice[k].get_multipoles().get_multipole(2).real
                print("{:5s}: Multipole, L = {:7.5f}, B_3 = {:10.5f}".
                      format(fam_name, L, b_3), end="", file=outf)
                if b_2 != 0e0:
                    print(", B_2 = {:8.5f}".format(b_2), end="", file=outf)
                print(", N = n_sext;", file=outf)
            else:
                print("\nprt_sext - undefined element type: {:s}".
                      format(fam_name))
                assert False
        return prt_list

    def prt_oct(fam_name, prt_list):
        if prt_list.count(fam_name) == 0:
            k = lat_prop._lattice.find(fam_name, 0).index
            L = lat_prop._lattice[k].get_length()
            b_4 = lat_prop._lattice[k].get_multipoles().get_multipole(4).real
            print("{:5s}: Octupole, L = {:7.5f}, B_4 = {:12.5e}, N = n_sext;"
                  .format(fam_name, L, b_4), file=outf)
        return prt_list

    # Dictionary of parameter types and corresponding print functions.
    prt_prm_func_dict = {
        "L":   prt_drift,
        "L_b": prt_dip,
        "phi": prt_dip,
        "h": prt_dip,
        "b_2": prt_quad,
        "b_3": prt_sext,
        "b_4": prt_oct
    }

    # Print family once.
    prt_list = []
    for k in range(len(prm_list._prm_list)):
        if type(prm_list._prm_list[k][0]) == str:
            prt_list = \
                prt_prm_func_dict[prm_list._prm_list[k][1]](
                    prm_list._prm_list[k][0], prt_list)
        elif isinstance(prm_list._prm_list[k][0], bend_class):
            prt_list = \
                prt_bend(prm_list._prm_list[k][0], prt_list)
        else:
            print("\nprt_lat - undefined parameter:", prm_list._prm_list[k][1])

    if (bend_list != None) and (bend_list != []):
        for bend in bend_list:
            if type(bend) == str:
                prt_list = \
                    prt_dip(bend, prt_list)
            elif isinstance(bend, bend_class):
                prt_list = \
                    prt_bend(bend, prt_list)
            else:
                print("\nprt_lat - undefined parameter:", bend._dip_prm)

    outf.close()

__all__ = [prm_class, phi_lat_class, bend_class, prt_lat]
