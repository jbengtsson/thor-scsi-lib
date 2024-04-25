"""Class to get/set multipole magnet parameters.
"""

import enum

import math
import numpy as np

import thor_scsi.lib as ts


class MultipoleIndex(enum.IntEnum):
    quadrupole = 2
    sextupole = 3


class get_set_mpole_class:
    # Private

    def __init__(self):
        pass

    # Public.

    def get_n_kid(self, fam_name):
        return len(self._lattice.elements_with_name(fam_name))

    def set_RF_cav_phase(self, fam_name, phi):
        # Set RF cavity phase for negative alpha_c.
        cav = self._lattice.find(fam_name, 0)
        cav.set_phase(phi)

    def get_L_elem(self, fam_name, kid_num):
        elem = self._lattice.find(fam_name, kid_num)
        return elem.get_length()

    def set_L_elem(self, fam_name, kid_num, L):
        elem = self._lattice.find(fam_name, kid_num)
        elem.set_length(L)

    def set_L_fam(self, fam_name, L):
        n_kid = len(self._lattice.elements_with_name(fam_name))
        for k in range(n_kid):
            self.set_L_elem(fam_name, k, L)

    def set_L_bend_elem(self, fam_name, kid_num, L):
        elem = self._lattice.find(fam_name, kid_num)
        phi = self.get_phi_elem(fam_name, kid_num)
        self.set_L_elem(fam_name, kid_num, L)
        self.set_phi_elem(fam_name, kid_num, phi)

    def set_L_bend_fam(self, fam_name, L):
        n_kid = len(self._lattice.elements_with_name(fam_name))
        for k in range(n_kid):
            self.set_L_bend_elem(fam_name, k, L)

    def get_phi_elem(self, fam_name, kid_num):
        b = self._lattice.find(fam_name, kid_num)
        return math.degrees(b.get_length() * b.get_curvature())

    def set_phi_elem(self, fam_name, kid_num, phi):
        b = self._lattice.find(fam_name, kid_num)
        L = b.get_length()
        h = math.radians(phi) / L
        b.set_curvature(h)

    def set_phi_fam(self, fam_name, phi):
        n_kid = len(self._lattice.elements_with_name(fam_name))
        for k in range(n_kid):
            self.set_phi_elem(fam_name, k, phi)

    def set_dphi_elem(self, fam_name, kid_num, dphi):
        b = self._lattice.find(fam_name, kid_num)
        L = b.get_length()
        # phi [rad].
        phi = b.get_length() * b.get_curvature()
        h = (phi + dphi) / L
        b.set_curvature(h)

    def set_dphi_fam(self, fam_name, dphi):
        n_kid = len(self._lattice.elements_with_name(fam_name))
        for k in range(n_kid):
            self.set_dphi_elem(fam_name, k, dphi)

    def set_phi_rect_elem(self, fam_name, kid_num, phi):
        self.set_phi_elem(fam_name, k, phi)
        b = self._lattice.find(fam_name, kid_num)
        b.set_entrance_angle(phi/2e0)
        b.set_exit_angle(phi/2e0)

    def set_phi_rect_fam(self, fam_name, phi):
        n_kid = len(self._lattice.elements_with_name(fam_name))
        for k in range(n_kid):
            self.set_phi_elem(fam_name, k, phi, rect_bend)

    def compute_phi_lat(self):
        """Compute the total bend angle.
        """
        prt = False
        phi = 0e0
        for k in range(len(self._lattice)):
            if (type(self._lattice[k]) == ts.Bending) \
               or (type(self._lattice[k]) == ts.Quadrupole):
                dphi = self.get_phi_elem(self._lattice[k].name, 0)
                phi += dphi
                if prt:
                    print("{:8s} {:5.3f} {:6.3f}".
                          format(self._lattice[k].name,
                                 self._lattice[k].get_length(), dphi))
        return phi

    def compute_scl_fact(self, bend_list):
        phi_b = 0e0
        for k in range(len(bend_list)):
            phi_b += self.get_phi_elem(bend_list[k], 0)
        print("\ncompute_scl_fact: phi_b = {:7.5f}".format(phi_b))
        bend_scl = []
        for k in range(len(bend_list)):
            phi = self.get_phi_elem(bend_list[k], 0)
            bend_scl.append(phi/phi_b)
        return np.array(bend_scl)

    def get_b_n_elem(self, fam_name, kid_num, n):
        mp = self._lattice.find(fam_name, kid_num)
        return mp.get_multipoles().get_multipole(n).real

    def set_b_n_elem(self, fam_name, kid_num, n, b_n):
        mp = self._lattice.find(fam_name, kid_num)
        mp.get_multipoles().set_multipole(n, b_n)

    def set_b_n_fam(self, fam_name, n, b_n):
        n_kid = len(self._lattice.elements_with_name(fam_name))
        for k in range(n_kid):
            self.set_b_n_elem(fam_name, k, n, b_n)

    def get_b_2_elem(self, fam_name, kid_num):
        return self.get_b_n_elem(fam_name, kid_num, 2)

    def set_b_2_fam(self, fam_name, b_n):
        self.set_b_n_fam(fam_name, 2, b_n)


__all__ = [get_set_mpole_class]
