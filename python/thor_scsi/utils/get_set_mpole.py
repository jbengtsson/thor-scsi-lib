"""Class to get/set multipole magnet parameters.
"""


import numpy as np

import thor_scsi.lib as ts


class get_set_mpole_class:
        # Private

    def __init__(self):
        pass

    # Public.

    def set_RF_cav_phase(self, lat, name, phi):
        # Set RF cavity phase for negative alpha_c.
        cav = lat.find(name, 0)
        cav.set_phase(phi)

    def get_phi_elem(self, lat, fam_name, n_kid):
        elem = lat.find(fam_name, n_kid)
        return elem.get_length() * elem.get_curvature() * 180e0 / np.pi

    def set_phi_elem(self, lat, fam_name, kid_num, phi, rect_bend):
        b = lat.find(fam_name, kid_num)
        L = b.get_length()
        h = phi * np.pi / (L * 180e0)
        b.set_curvature(h)
        if rect_bend:
            b.set_entrance_angle(phi/2e0)
            b.set_exit_angle(phi/2e0)

    def set_phi_fam(self, lat, fam_name, phi, rect_bend):
        prt = False
        b = lat.find(fam_name, 0)
        L = b.get_length()
        h = phi * np.pi / (L * 180e0)
        if prt:
            print("set_phi_fam: {:8s} {:5.3f} {:6.3f}".format(b.name, L, phi))
        for b in lat.elements_with_name(fam_name):
            b.set_curvature(h)

    def compute_phi(self, lat):
        """Compute the total bend angle.
        """
        prt = False
        phi = 0e0
        for k in range(len(lat)):
            if (type(lat[k]) == ts.Bending) or (type(lat[k]) == ts.Quadrupole):
                dphi = self.get_phi_elem(lat, lat[k].name, 0)
                phi += dphi
                if prt:
                    print("{:8s} {:5.3f} {:6.3f}".
                          format(lat[k].name, lat[k].get_length(), dphi))
        return phi

    def get_b_n_elem(self, lat, fam_name, kid_num, n):
        mp = lat.find(fam_name, kid_num)
        return mp.get_multipoles().get_multipole(n).real

    def set_b_n_elem(self, lat, fam_name, kid_num, n, b_n):
        mp = lat.find(fam_name, kid_num)
        mp.get_multipoles().set_multipole(n, b_n)

    def set_b_n_fam(self, lat, fam_name, n, b_n):
        for mp in lat.elements_with_name(fam_name):
            mp.get_multipoles().set_multipole(n, b_n)

    def get_L_elem(self, lat, fam_name, n_kid):
        elem = lat.find(fam_name, n_kid)
        return elem.get_length()

    def set_L_fam(self, lat, fam_name, L):
        for elem in lat.elements_with_name(fam_name):
            elem.set_length(L)

    def get_b_2_elem(self, lat, fam_name, kid_num):
        return get_b_n_elem(lat, fam_name, kid_num, MultipoleIndex.quadrupole)

    def set_b_2_fam(self, lat, fam_name, b_2):
        set_b_n_fam(lat, fam_name, MultipoleIndex.quadrupole, b_2)


__all__ = [get_set_mpole_class]
