"""Class for linear optics analysis of a periodic structure.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.factory import accelerator_from_config

from . import linear_optics as lo, courant_snyder as cs, index_class as ind
from .output import mat2txt, vec2txt

# from thor_scsi.utils.output import vec2txt


ind = ind.index_class()


class periodic_structure_class:
    # Private.

    def __init__(self, nv, no, nv_prm, no_prm, file_name, E_0):
        self._prt_Twiss_file_name = "twiss.txt"

        self._named_index = []
        self._desc        = []

        self._lattice     = []
        self._model_state = []
        self._no          = no
        self._n_dof       = np.nan

        self._M           = np.nan
        self._alpha_c     = np.nan
        self._A           = np.nan

        self._data        = []
              
        # Initialise.

        self._named_index = \
            gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

        # Descriptor for Truncated Power Series Algebra variables.
        self._desc = gtpsa.desc(nv, no, nv_prm, no_prm)

        # Read in & parse lattice file.
        self._lattice = accelerator_from_config(file_name)
        # Set lattice state (Rf cavity on/off, etc.).
        self._model_state = ts.ConfigType()

        # D.O.F. (Degrees-Of-Freedom) - coasting beam.
        self._n_dof = 2
        self._model_state.Energy = E_0
        self._model_state.radiation = False
        self._model_state.Cavity_on = False

    # Public.

    def compute_circ(self):
        return np.sum([elem.get_length() for elem in self._lattice])

    def compute_alpha_c(self):
        self._alpha_c = \
            self._M.iloc[ind.ct].getv(1)[ind.delta]/self.compute_circ()

    def get_type(self, elem):
        match type(elem):
            case ts.Marker:
                type_code = 0.0
            case ts.Drift:
                type_code = 0.0
            case ts.Bending:
                if elem.get_multipoles().get_multipole(2).real == 0e0:
                    type_code = np.sign(elem.get_curvature())*0.5
                else:
                    type_code = \
                        np.sign(elem.get_multipoles().get_multipole(2).real) \
                        *0.75
            case ts.Quadrupole:
                type_code = \
                    np.sign(elem.get_multipoles().get_multipole(2).real)*1.0
            case ts.Sextupole:
                type_code = \
                    np.sign(elem.get_multipoles().get_multipole(3).real)*1.25
            case ts.Octupole:
                type_code = \
                    np.sign(elem.get_multipoles().get_multipole(4).real)*1.5
            case ts.Cavity:
                type_code = 0.0 
            case _:
                print("\nget_type - undefined type:", elem.name)
                type_code = np.nan
        return type_code

    def get_types(self):
        n = len(self._lattice)
        self._types = np.zeros(n)
        for k in range(n):
            self._types[k] = self.get_type(self._lattice[k])

    def plt_Twiss(self, file_name, plot):
        fig, (gr_1, gr_2) = plt.subplots(2)

        fig.suptitle("Linear Optics")

        gr_1.set_title("Beta Functions")
        gr_1.set_xlabel("s [m]")
        gr_1.set_ylabel(r"$\beta_{x,y}$ [m]")
        gr_1.plot(
            self._data.s, self._data.twiss.sel(plane="x", par="beta"), "b",
            label=r"$\beta_x$")
        gr_1.plot(
            self._data.s, self._data.twiss.sel(plane="y", par="beta"), "r",
            label=r"$\beta_y$")
        gr_1.legend()

        gr_1_r = gr_1.twinx()
        gr_1_r.set_ylim([-2.0, 20.0])
        gr_1_r.set_yticks([])
        gr_1_r.step(self._data.s, self._types, "k")

        gr_2.set_title("Dispersion")
        gr_2.set_xlabel("s [m]")
        gr_2.set_ylabel(r"$\eta_x$ [m]")
        gr_2.plot(
            self._data.s, self._data.dispersion.sel(phase_coordinate="x"), "b")

        gr_2_r = gr_2.twinx()
        gr_2_r.set_ylim([-2.0, 20.0])
        gr_2_r.set_yticks([])
        gr_2_r.step(self._data.s, self._types, "k")

        fig.tight_layout()

        plt.savefig(file_name)
        print("\nplt_Twiss - plot saved as:", file_name)

        if plot:
            plt.show()

    def plt_scan_phi_rb(self, file_name, phi, eps_x, J_x, J_z, alpha_c, plot):
        fig, (gr_1, gr_2) = plt.subplots(2)

        fig.suptitle("Lattice Trade-Offs vs. Reverse Bend Angle")

        gr_1.set_title(
            r"$[\epsilon_x\left(\phi_{RB}\right)$"
            r"$, \alpha_c\left(\phi_{RB}\right)]$")
        gr_1.set_xlabel(r"$\phi_{RB}$ [$\degree$]")
        gr_1.set_ylabel(r"$\epsilon_x$ [pm$\cdot$rad]")
        gr_1.yaxis.label.set_color("b")
        gr_1.tick_params(axis="y", colors="b")
        gr_1.plot(phi, 1e12*eps_x, color="b")

        gr_1_r = gr_1.twinx()
        gr_1_r.set_ylabel(r"$\alpha_c$ [$10^{-4}$]")
        gr_1_r.yaxis.label.set_color("g")
        gr_1_r.tick_params(axis="y", colors="g")
        gr_1_r.plot(phi, 1e4*alpha_c, color="g", label=r"$\alpha_c$")

        gr_2.set_title(
            r"$[J_x\left(\phi_{RB}\right), J_z\left(\phi_{RB}\right)]$")
        gr_2.set_xlabel(r"$\phi_{RB}$ [$\degree$]")
        gr_2.set_ylabel(r"$J_x$")
        gr_2.yaxis.label.set_color("b")
        gr_2.tick_params(axis="y", colors="b")
        gr_2.plot(phi, J_x, color="b")

        gr_2_r = gr_2.twinx()
        gr_2_r.set_ylabel(r"$J_z$")
        gr_2_r.yaxis.label.set_color("g")
        gr_2_r.tick_params(axis="y", colors="g")
        gr_2_r.plot(phi, J_z, color="g")

        fig.tight_layout()

        plt.savefig(file_name)
        print("\n - plot saved as:", file_name)

        if plot:
            plt.show()

    def prt_default_mapping(self):
        index_map = gtpsa.default_mapping()
        print("\nIndex Mapping:\n"
              "  x     = {:1d}\n  p_x   = {:1d}\n  y     = {:1d}"
              "\n  p_y   = {:1d}\n  delta = {:1d}\n  ct    = {:1d}".
              format(index_map.x, index_map.px, index_map.y, index_map.py,
                     index_map.delta, index_map.ct))
    
    def prt_M(self):
        n_dof = 3
        print("\nM:\ntpsa cst:")
        for k in range(2*n_dof):
            print(" {:13.6e}".format(self._M.cst().iloc[k]), end="")
        print("\ntpsa linear:\n"+mat2txt(self._M.jacobian()[:6, :6]))

    def prt_Twiss_param(self, Twiss):
        """
        """
        eta, alpha, beta, nu = Twiss
        print("\nTwiss:")
        print(f"  eta   = [{eta[ind.X]:9.3e}, {eta[ind.Y]:9.3e}]")
        print(f"  alpha = [{alpha[ind.X]:9.3e}, {alpha[ind.Y]:9.3e}]")
        print(f"  beta  = [{beta[ind.X]:5.3f}, {beta[ind.Y]:5.3f}]")


    def prt_lat_param(self):
        eta = np.zeros(2)
        alpha = np.zeros(2)
        beta = np.zeros(2)
        nu = np.zeros(2)
        loc = len(self._lattice)-1
        eta[ind.X] = self._data.dispersion.sel(phase_coordinate="x").values[loc]
        eta[ind.Y] = self._data.dispersion.sel(phase_coordinate="y").values[loc]
        alpha[ind.X] = self._data.twiss.sel(plane="x", par="alpha").values[loc]
        alpha[ind.Y] = self._data.twiss.sel(plane="y", par="alpha").values[loc]
        beta[ind.X] = self._data.twiss.sel(plane="x", par="beta").values[loc]
        beta[ind.Y] = self._data.twiss.sel(plane="y", par="beta").values[loc]
        nu[ind.X] = self._data.twiss.sel(plane="x", par="nu").values[loc]
        nu[ind.Y] = self._data.twiss.sel(plane="y", par="nu").values[loc]
        print("\nTwiss:")
        print(f"  eta     = [{eta[ind.X]:9.3e}, {eta[ind.Y]:9.3e}]")
        print(f"  alpha   = [{alpha[ind.X]:9.3e}, {alpha[ind.Y]:9.3e}]")
        print(f"  beta    = [{beta[ind.X]:5.3f}, {beta[ind.Y]:5.3f}]")
        print(f"  nu      = [{nu[ind.X]:7.5f}, {nu[ind.Y]:7.5f}]")
        print(f"  alpha_c = {self._alpha_c:10.3e}")


    def prt_Twiss(self):
        """
        Print Twiss parameters along the lattice.
        """
        file = open(self._prt_Twiss_file_name, "w")
        nu = np.zeros(2, dtype=float)
        print("\n     Name          s    type  alpha_x   beta_x  nu_x    eta_x"
              "   eta'_x    alpha_y   beta_y  nu_y    eta_y   eta'_y",
              file=file)
        print("                  [m]                    [m]             [m]"
              "     [m]                 [m]             [m]", file=file)
        for k in range(len(self._data.index)):
            print("{:3d} {:10s} {:7.3f} {:4.1f} {:9.5f} {:8.5f} {:7.5f} {:7.5f}"
                  " {:8.5f} {:9.5f} {:8.5f} {:7.5f} {:7.5f} {:8.5f}".
                  format(k, self._lattice[k].name,  self._data.s[k],
                         self._types[k],
                         self._data.twiss.sel(plane="x", par="alpha").values[k],
                         self._data.twiss.sel(plane="x", par="beta").values[k],
                         self._data.twiss.sel(plane="x", par="nu").values[k],
                         self._data.dispersion.sel(
                             phase_coordinate="x").values[k],
                         self._data.dispersion.sel(
                             phase_coordinate="px").values[k],
                         self._data.twiss.sel(plane="y", par="alpha").values[k],
                         self._data.twiss.sel(plane="y", par="beta").values[k],
                         self._data.twiss.sel(plane="y", par="nu").values[k],
                         self._data.dispersion.sel(
                             phase_coordinate="y").values[k],
                         self._data.dispersion.sel(
                             phase_coordinate="py").values[k]),
                  file=file)

        print("\nprt_Twiss - Linear optics saved as:",
              self._prt_Twiss_file_name)

    def get_Twiss(self, loc):
        eta = np.array([
            self._data.dispersion.sel(phase_coordinate="x").values[loc],
            self._data.dispersion.sel(phase_coordinate="px").values[loc],
            self._data.dispersion.sel(phase_coordinate="y").values[loc],
            self._data.dispersion.sel(phase_coordinate="py").values[loc],
        ])
        alpha = np.array([
                self._data.twiss.sel(plane="x", par="alpha").values[loc],
                self._data.twiss.sel(plane="y", par="alpha").values[loc],
            ])
        beta = np.array([
                self._data.twiss.sel(plane="x", par="beta").values[loc],
                self._data.twiss.sel(plane="y", par="beta").values[loc],
            ])
        nu = np.array([
                self._data.twiss.sel(plane="x", par="nu").values[loc],
                self._data.twiss.sel(plane="y", par="nu").values[loc],
            ])
        return eta, alpha, beta, nu,

    def comp_per_sol(self):
        """Compute the periodic solution for a super period.
        Degrees of freedom - RF cavity is off; i.e., coasting beam.
        """
        self._n_dof = 2
        self._model_state.radiation = False
        self._model_state.Cavity_on = False

        stable, self._M, self._nu, self._A = \
            lo.compute_map_and_diag(
                self._n_dof, self._lattice, self._model_state, desc=self._desc)
        if stable:
            self.compute_alpha_c()
            A_map = gtpsa.ss_vect_tpsa(self._desc, self._no)
            A_map.set_jacobian(self._A)
            self._data = \
                lo.compute_Twiss_along_lattice(
                    self._n_dof, self._lattice, self._model_state, A=A_map,
                    desc=self._desc, mapping=self._named_index)
        else:
            self._data = np.nan
            print("\ncomp_per_sol: unstable")

        return stable

    def unit_cell_rev_bend(
            self, n_step, phi_min, set_phi_rb, phi_b, bend_list, bend_scl):
        phi_rb = 0e0
        phi_step = phi_min/n_step
        phi_rb_buf = []
        eps_x_buf = []
        alpha_c_buf = []
        J_x_buf = []
        J_z_buf = []
        print("\n   phi   phi_tot  eps_x    J_x   J_z    alpha_c     eta_x"
              "     nu_x     nu_y"
              "\n  [deg]   [deg]  [nm.rad]                            [m]")
        for k in range(n_step):
            set_phi_rb(self, phi_b, bend_list, bend_scl, phi_rb)
            stable = self.comp_per_sol()
            Twiss = self.get_Twiss(len(self._lattice)-1)
            eta_x = Twiss[0][ind.x]
            self.compute_alpha_c()
            if self._alpha_c > 0e0:
                self.set_RF_cav_phase("cav", 0.0)
            else:
                self.set_RF_cav_phase("cav", 180.0)
            stable = self.compute_radiation()
            if stable:
                phi_rb_buf.append(abs(phi_rb))
                eps_x_buf.append(self._eps[ind.X])
                J_x_buf.append(self._J[ind.X])
                J_z_buf.append(self._J[ind.Z])
                alpha_c_buf.append(self._alpha_c)
                print("{:7.3f}  {:5.3f}    {:5.1f}    {:4.2f} {:5.2f} {:10.3e}"
                      " {:10.3e}  {:7.5f}  {:7.5f}".
                      format(
                          phi_rb, self.compute_phi(), 1e12*self._eps[ind.X],
                          self._J[ind.X], self._J[ind.Z], self._alpha_c, eta_x,
                          self._nu[ind.X], self._nu[ind.Y]))
            else:
                self.compute_alpha_c()
                print("unit_cell_rev_bend: unstable alpha_c = {:10.3e}".
                      format(self._alpha_c))
            phi_rb += phi_step
        return \
            np.array(phi_rb_buf), np.array(eps_x_buf), np.array(J_x_buf), \
            np.array(J_z_buf), np.array(alpha_c_buf)


__all__ = [periodic_structure_class]
