"""Class to get/set multipole magnet parameters.
"""

import enum

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

    def set_RF_cav_phase(self, name, phi):
        # Set RF cavity phase for negative alpha_c.
        cav = self._lattice.find(name, 0)
        cav.set_phase(phi)

    def get_phi_elem(self, fam_name, n_kid):
        elem = self._lattice.find(fam_name, n_kid)
        return elem.get_length() * elem.get_curvature() * 180e0 / np.pi

    def set_phi_elem(self, fam_name, kid_num, phi, rect_bend):
        b = self._lattice.find(fam_name, kid_num)
        L = b.get_length()
        h = phi * np.pi / (L * 180e0)
        b.set_curvature(h)
        if rect_bend:
            b.set_entrance_angle(phi/2e0)
            b.set_exit_angle(phi/2e0)

    def set_phi_fam(self, fam_name, phi, rect_bend):
        n = len(self._lattice.elements_with_name(fam_name))
        for k in range(n):
            self.set_phi_elem(fam_name, k, phi, rect_bend)

    def set_dphi_elem(self, fam_name, kid_num, dphi, rect_bend):
        phi = self.get_phi_elem(fam_name, kid_num) + dphi
        self.set_phi_elem(fam_name, kid_num, phi, rect_bend)

    def set_dphi_fam(self, fam_name, phi, rect_bend):
        n = len(self._lattice.elements_with_name(fam_name))
        for k in range(n):
            self.set_dphi_elem(fam_name, k, phi, rect_bend)

    def compute_phi(self):
        """Compute the total bend angle.
        """
        phi = 0e0
        for k in range(len(self._lattice)):
            if (type(self._lattice[k]) == ts.Bending) \
               or (type(self._lattice[k]) == ts.Quadrupole):
                dphi = self.get_phi_elem(self._lattice[k].name, 0)
                phi += dphi
        return phi

    def get_b_n_elem(self, fam_name, kid_num, n):
        mp = self._lattice.find(fam_name, kid_num)
        return mp.get_multipoles().get_multipole(n).real

    def set_b_n_elem(self, fam_name, kid_num, n, b_n):
        mp = self._lattice.find(fam_name, kid_num)
        mp.get_multipoles().set_multipole(n, b_n)

    def set_b_n_fam(self, fam_name, n, b_n):
        n = len(self._lattice.elements_with_name(fam_name))
        for k in range(n):
            self.set_b_n_elem(fam_name, k, b_n)

    def get_L_elem(self, fam_name, n_kid):
        elem = self._lattice.find(fam_name, n_kid)
        return elem.get_length()

    def set_L_elem(self, fam_name, kid_num, L):
        elem = self._lattice.find(fam_name, kid_num)
        elem.set_length(L)

    def set_L_fam(self, fam_name, L):
        n = len(self._lattice.elements_with_name(fam_name))
        for k in range(n):
            self.set_L_elem(fam_name, n, L)

    def get_b_2_elem(self, fam_name, kid_num):
        return self.get_b_n_elem(
            fam_name, kid_num, MultipoleIndex.quadrupole)

    def set_b_2_fam(self, fam_name, b_2):
        n = len(self._lattice.elements_with_name(fam_name))
        self.set_b_n_(fam_name, 2, b_2)


__all__ = [get_set_mpole_class]
