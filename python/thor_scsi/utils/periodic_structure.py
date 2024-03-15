"""Class for linear optics analysis of a periodic structure.
"""

import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
     radiate as rad

from thor_scsi.utils.output import mat2txt, vec2txt


# Configuration space coordinates.
X_, Y_, Z_ = [
    ts.spatial_index.X,
    ts.spatial_index.Y,
    ts.spatial_index.Z
]
# Phase-space coordinates.
[x_, px_, y_, py_, ct_, delta_] = [
    ts.phase_space_index_internal.x,
    ts.phase_space_index_internal.px,
    ts.phase_space_index_internal.y,
    ts.phase_space_index_internal.py,
    ts.phase_space_index_internal.ct,
    ts.phase_space_index_internal.delta,
]

class lin_opt_class:
    # Private

    def __init__(self, nv, no, nv_prm, no_prm, file_name, E_0, cod_eps):
        self._prt_Twiss_file_name = "twiss.txt"
        self._plt_Twiss_file_name = "twiss.png"

        self._named_index = []
        self._desc        = []

        self._lattice     = []
        self._model_state = []
        self._no          = no
        self._E_0         = E_0
        self._n_dof       = np.nan
        self._cod_eps     = cod_eps

        self._M           = np.nan
        self._A           = np.nan

        self._Twiss       = np.nan
        self._data        = []
        
        self._nu          = np.nan

        self._cod         = np.nan
        self._U_0         = np.nan
        self._J           = np.nan
        self._tau         = np.nan
        self._eps         = np.nan
        self._D_rad       = np.nan
        
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
        self._model_state.radiation = False
        self._model_state.Cavity_on = False

    # Public.

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
            case _:
                type_code = np.nan
        return type_code

    def get_types(self):
        n = len(self._lattice)
        types = np.zeros(n)
        for k in range(n):
            types[k] = self.get_type(self._lattice[k])
        return types

    def plt_Twiss(self, types, plot):
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
        gr_1_r.step(self._data.s, types, "k")

        gr_2.set_title("Dispersion")
        gr_2.set_xlabel("s [m]")
        gr_2.set_ylabel(r"$\eta_x$ [m]")
        gr_2.plot(
            self._data.s, self._data.dispersion.sel(phase_coordinate="x"), "b")

        gr_2_r = gr_2.twinx()
        gr_2_r.set_ylim([-2.0, 20.0])
        gr_2_r.set_yticks([])
        gr_2_r.step(self._data.s, types, "k")

        fig.tight_layout()

        plt.savefig(self._plt_Twiss_file_name)
        print("\nplt_Twiss - plot saved as:", self._plt_Twiss_file_name)

        if plot:
            plt.show()

    def prt_Twiss_param(self):
        """
        """
        eta, alpha, beta = self._Twiss
        self._nu = lo.compute_nu_symp(2, self._M.jacobian())
        print("\nTwiss:")
        print(f"  eta   = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
        print(f"  alpha = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
        print(f"  beta  = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")
        print(f"  nu    = [{self._nu[X_]:7.5f}, {self._nu[Y_]:7.5f}]")

    def prt_Twiss(self, types):
        """
        Print Twiss parameters along the lattice.
        """
        file = open(self._prt_Twiss_file_name, 'w')
        s = 0e0
        nu = np.zeros(2, dtype=float)
        print("\n    Name         s    type  alpha_x   beta_x  nu_x    eta_x"
              "   eta'_x    alpha_y   beta_y  nu_y    eta_y   eta'_y",
              file=file)
        print("                [m]                    [m]             [m]"
              "     [m]                 [m]             [m]", file=file)
        for k in range(len(self._data.index)):
            s += self._lattice[k].get_length()
            nu[X_] += self._data.twiss.sel(plane="x", par="dnu").values[k]
            nu[Y_] += self._data.twiss.sel(plane="y", par="dnu").values[k]
            print("{:3d} {:8s} {:7.3f} {:4.1f} {:9.5f} {:8.5f} {:7.5f} {:7.5f}"
                  " {:8.5f} {:9.5f} {:8.5f} {:7.5f} {:7.5f} {:8.5f}".
                  format(k, self._lattice[k].name, s,
                         self.get_type(self._lattice[k]),
                         self._data.twiss.sel(plane="x", par="alpha").values[k],
                         self._data.twiss.sel(plane="x", par="beta").values[k],
                         nu[X_],
                         self._data.dispersion.sel(
                             phase_coordinate="x").values[k],
                         self._data.dispersion.sel(
                             phase_coordinate="px").values[k],
                         self._data.twiss.sel(plane="y", par="alpha").values[k],
                         self._data.twiss.sel(plane="y", par="beta").values[k],
                         nu[Y_],
                         self._data.dispersion.sel(
                             phase_coordinate="y").values[k],
                         self._data.dispersion.sel(
                             phase_coordinate="py").values[k]),
                  file=file)

        print("\nprt_Twiss - Linear optics saved as:",
              self._prt_Twiss_file_name)


    def comp_per_sol(self):
        """
        Compute the periodic solution for a super period.
        Degrees of freedom - RF cavity is off; i.e., coasting beam.
        """
        self._n_dof = 2
        self._model_state.radiation = False
        self._model_state.Cavity_on = False

        stable, self._M, self._A = \
            lo.compute_map_and_diag(
                self._n_dof, self._lattice, self._model_state, desc=self._desc)
        if stable:
            res = cs.compute_Twiss_A(self._A)
            self._Twiss = res[:3]
            A_map = gtpsa.ss_vect_tpsa(self._desc, self._no)
            A_map.set_jacobian(self._A)
            self._data = \
                lo.compute_Twiss_along_lattice(
                    self._n_dof, self._lattice, self._model_state, A=A_map,
                    desc=self._desc, mapping=self._named_index)
        else:
            self._data = np.nan
            print("\ncomp_per_sol: unstable")

    def prt_M(self):
        n_dof = 3
        print("\nM:\ntpsa cst:")
        for k in range(2*n_dof):
            print(" {:13.6e}".format(self._M.cst().iloc[k]), end="")
        print("\ntpsa linear:\n"+mat2txt(self._M.jacobian()[:6, :6]))

    def comp_rad(self):
        stable, self._M, self._cod, self._A, self._U_0, self._J, self._tau, \
            self._eps, self._D_rad = \
                rad.compute_radiation(
                    self._lattice, self._model_state, self._E_0, self._cod_eps,
                    desc=self._desc)

    def prt_rad(self):
        print("\nRadiation Properties:")
        print("  eps_x [m.rad] = [{:9.3e}, {:9.3e}, {:9.3e}]".format(
            self._eps[X_], self._eps[Y_], self._eps[Z_]))
        print("  J             = [{:5.3f}, {:5.3f}, {:5.3f}]".format(
            self._J[X_], self._J[Y_], self._J[Z_]))
        print("  tau [msec]    = [{:5.3f}, {:5.3f}, {:5.3f}]".format(
            1e3*self._tau[X_], 1e3*self._tau[Y_], 1e3*self._tau[Z_]))
        print("  D             = [{:11.5e}, {:11.5e}, {:11.5e}]".format(
            self._D_rad[X_], self._D_rad[Y_], self._D_rad[Z_]))


__all__ = [lin_opt_class]
