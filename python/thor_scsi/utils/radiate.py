"""
"""
import thor_scsi.lib as tslib

import gtpsa
from .accelerator import instrument_with_radiators
import numpy as np
from scipy.constants import c as c0
from dataclasses import dataclass

from . import periodic_structure as ps
from .closed_orbit import compute_closed_orbit
from .linear_optics import compute_M_diag
from .output import mat2txt, vec2txt

import logging
logger = logging.getLogger("thor-scsi")


X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]

x_, px_, y_, py_, ct_, delta_ = [
    tslib.phase_space_index_internal.x,
    tslib.phase_space_index_internal.px,
    tslib.phase_space_index_internal.y,
    tslib.phase_space_index_internal.py,
    tslib.phase_space_index_internal.ct,
    tslib.phase_space_index_internal.delta
]


@dataclass
class RadiationResult:
    #: the different taus
    relaxation_constants: np.ndarray
    emittance: np.ndarray
    #: for consistency checks
    fractional_tunes: np.ndarray


class lat_prop_class(ps.lin_opt_class):
    # Private

    def __init__(self, nv, no, nv_prm, no_prm, file_name, E_0, cod_eps):
        print("\nlat_prop_class - __init__")
        super().__init__(nv, no, nv_prm, no_prm, file_name, E_0)
        self._cod_eps       = cod_eps
        self._rad_del_kicks = []
        self._M             = []       # ss_vect_tpsa.
        self._nu            = np.nan
        self._A             = np.nan
        self._dE            = np.nan
        self._U_0           = np.nan
        self._J             = np.nan
        self._alpha_rad     = np.nan
        self._tau           = np.nan
        self._eps           = np.nan
        self._D_rad         = np.nan

    # Public.

    def compute_circ(self):
        return np.sum([elem.get_length() for elem in self._lattice])

    def compute_diffusion_coefficients(self):
        dD_rad = \
            np.array([rk.get_diffusion_coefficients_increments()
                      for rk in self._rad_del_kicks])
        self._D_rad = np.sum(dD_rad, axis=0)

    def compute_rad_prop(self):
        dof = 3
        self._J = np.zeros(dof)
        self._tau = np.zeros(dof)
        self._eps = np.zeros(dof)
        C = self.compute_circ()
        logger.info("\nC = %5.3f", C)
        self._U_0 = self._model_state.Energy*self._dE
        for k in range(dof):
            self._J[k] = \
                2e0*(1e0+self._M.cst().delta)*self._alpha_rad[k]/self._dE
            self._tau[k] = -C/(c0*self._alpha_rad[k])
            self._eps[k] = -self._D_rad[k]/(2e0*self._alpha_rad[k])

        logger.info(
            "\nE [GeV]     = {:3.1f}\nU0 [keV]    = {:3.1f}\neps         ="
            " {:12.6e} {:12.6e} {:12.6e}\ntau [msec]  = {:8.6f} {:8.6f} {:8.6f}"
            "\nJ           = {:8.6f} {:8.6f} {:8.6f}\nalpha_rad   = {:13.6e}"
            " {:13.6e} {:13.6e}\nD_rad       = {:12.6e} {:12.6e} {:12.6e}".
            format(1e-9*self._model_state.Energy, 1e-3*self._U_0,
                   self._eps[X_], self._eps[Y_], self._eps[Z_],
                   1e3*self._tau[X_], 1e3*self._tau[Y_], 1e3*self._tau[Z_],
                   self._J[X_], self._J[Y_], self._J[Z_], self._alpha_rad[X_],
                   self._alpha_rad[Y_], self._alpha_rad[Z_], self._D_rad[X_],
                   self._D_rad[Y_], self._D_rad[Z_]))

    def compute_radiation(self):
        dof = 3

        self._model_state.radiation = True
        self._model_state.emittance = False
        self._model_state.Cavity_on = True

        # Install radiators that radiation is calculated
        self._rad_del_kicks = \
            instrument_with_radiators(
                self._lattice, energy=self._model_state.Energy)

        r = \
            compute_closed_orbit(
                self._lattice, self._model_state, delta=0e0,
                eps=self._cod_eps)
        # self._M = r.one_turn_map[:6, :6]
        self._M = r.one_turn_map

        logger.info(
            "\nM:\n" + mat2txt(self._M.jacobian())
            + "\n\ncod ="
            + vec2txt(np.array(
                [self._M.cst().x, self._M.cst().px, self._M.cst().y,
                 self._M.cst().py, self._M.cst().delta, self._M.cst().ct]))
        )

        self._model_state.dE = 0e0
        ps = self._M.cst()
        # dE is computed by the RF cavity propagator.
        self._lattice.propagate(self._model_state, ps)
        self._dE = self._model_state.dE

        stable, self._nu, self._A, A_inv, self._alpha_rad = \
            compute_M_diag(dof, self._M.jacobian())

        A_7x7 = np.zeros((7, 7))
        A_7x7[:6, :6] = self._A
        A_7x7[6, 6] = 1e0
        if stable:
            self._model_state.emittance = True

            A_cpy  = gtpsa.ss_vect_tpsa(self._desc, 1)
            A_cpy += self._M.cst()
            A_cpy.set_jacobian(A_7x7)
            self._lattice.propagate(self._model_state, A_cpy)

            self.compute_diffusion_coefficients()

            self.compute_rad_prop()
        else:
            self._U_0 = np.nan
            self._J = np.zeros(3, float)
            self._tau = np.zeros(3, float)
            self._eps = np.zeros(3, float)

        return stable

    def prt_rad(self):
        print("\nRadiation Properties:")
        print("  E [GeV]       = {:5.3f}".
              format(1e-9*self._model_state.Energy))
        print("  U_0 [keV]     = {:5.1f}".
              format(1e-3*self._U_0))
        print("  eps_x [m.rad] = [{:9.3e}, {:9.3e}, {:9.3e}]".format(
            self._eps[X_], self._eps[Y_], self._eps[Z_]))
        print("  J             = [{:5.3f}, {:5.3f}, {:5.3f}]".format(
            self._J[X_], self._J[Y_], self._J[Z_]))
        print("  tau [msec]    = [{:5.3f}, {:5.3f}, {:5.3f}]".format(
            1e3*self._tau[X_], 1e3*self._tau[Y_], 1e3*self._tau[Z_]))
        print("  D             = [{:11.5e}, {:11.5e}, {:11.5e}]".format(
            self._D_rad[X_], self._D_rad[Y_], self._D_rad[Z_]))

    def prt_M(self):
        n_dof = 3
        print("\nM:\ntpsa cst:")
        for k in range(2*n_dof):
            print(" {:13.6e}".format(self._M.cst().iloc[k]), end="")
        print("\ntpsa linear:\n"+mat2txt(self._M.jacobian()[:6, :6]))


# def calculate_radiation(
#     lat: tslib.Accelerator,
#     *,
#     energy,
#     model_state: tslib.ConfigType = None,
#     install_radiators: bool = True,
#     dof=3
# ):
#     """

#     Todo:
#         Rename function
#         Inspect if radiators are installed?
#         check in calc config if radation is requested
#         in this case install the radiators

#     1. calculate fix point and Poincarè Map M with damped system
#        (i.e. radiation on and cavity on (without dispersion in a second case)
#     2. diagonalise M = A $\Gamma$ A$^{-1}$
#     3. eigenvalues:

#         - complex part: tunes,
#         - real part: damping times  (refer equation)

#        use eigen values of symplectic matrix to identify the planes

#     4. propagate A, thin kick will create diffusion coeffs (don't forget to
#        zero them before calculation starts (sum it up afterwards

#     """

#     if install_radiators:
#         logger.debug("Installing radiators")
#         # keep variable as long as you need to do calculations
#         rad_del = instrument_with_radiators(lat, energy=energy)

#     if model_state is None:
#         raise AssertionError
#         model_state = tslib.ConfigType()
#         model_state.radiation = True
#         # is this used anywhere?
#         model_state.emittance = False
#         model_state.Cavity_on = True

#     model_state.Energy = energy
#     logger.debug(
#         f"model_state radiation { model_state.radiation} emmittance {model_state.emittance} Cavity on {model_state.Cavity_on}"
#     )


#     # Diffusion coefficients
#     lat.propagate(model_state, Ap)
#     print("Ap after calc")
#     print(Ap)

#     r = RadiationResult(relaxation_constants=w.real, fractional_tunes=w.imag)


__all__ = [lat_prop_class]
