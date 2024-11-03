"""
"""
import thor_scsi.lib as tslib

import gtpsa
from .accelerator import instrument_with_radiators
import numpy as np
from scipy.constants import c as c0
from dataclasses import dataclass

from . import periodic_structure as ps, get_set_mpole as gs, index_class as ind
from .closed_orbit import compute_closed_orbit
from .linear_optics import compute_M_diag
from .output import mat2txt, vec2txt

import logging
logger = logging.getLogger("thor-scsi")


ind = ind.index_class()


@dataclass
class RadiationResult:
    #: the different taus
    relaxation_constants: np.ndarray
    emittance: np.ndarray
    #: for consistency checks
    fractional_tunes: np.ndarray


class lattice_properties_class(
        ps.periodic_structure_class, gs.get_set_mpole_class):
    # Private.

    def __init__(self, gtpsa_prop, file_name, E_0, cod_eps):
        # Initialise periodic_structure_class.
        super().__init__(gtpsa_prop, file_name, E_0)
        self._cod_eps       = cod_eps
        self._rad_del_kicks = []
        self._M_rad         = []       # ss_vect_tpsa.
        self._nu_rad        = np.nan
        self._A_rad         = np.nan
        self._dE            = np.nan
        self._U_0           = np.nan
        self._J             = np.nan
        self._alpha_rad     = np.nan
        self._tau           = np.nan
        self._eps           = np.nan
        self._D_rad         = np.nan
        self._sigma_s       = np.nan
        self._sigma_delta   = np.nan

    # Public.

    def __str__(self):
        return \
            f"\nlattice_properties_class:\n" \
            +f"  no    = {self._no:d}\n  nv    = {self._nv:d}\n" \
            +f"  n_dof = {self._n_dof:d}"

    def compute_circ(self):
        return np.sum([elem.get_length() for elem in self._lattice])

    def compute_diffusion_coefficients(self):
        dD_rad = \
            np.array([rk.get_diffusion_coefficients_increments()
                      for rk in self._rad_del_kicks])
        self._D_rad = np.sum(dD_rad, axis=0)

    def compute_rad_prop(self):
        dof = 3
        planes = ["x", "y", "z"]
        self._J = np.zeros(dof)
        self._tau = np.zeros(dof)
        self._eps = np.zeros(dof)
        logger.info("\nC = %5.3f", self._circ)
        self._U_0 = self._model_state.Energy*self._dE
        for k in range(dof):
            if self._alpha_rad[k] < 0e0:
                self._J[k] = \
                    2e0*(1e0+self._M_rad.cst().delta)*self._alpha_rad[k] \
                    /self._dE
                self._tau[k] = -self._circ/(c0*self._alpha_rad[k])
                self._eps[k] = -self._D_rad[k]/(2e0*self._alpha_rad[k])
            else:
                # Unstable.
                print("\ncompute_rad_prop:"+
                      " unstable alpha_{:s} >= 0 ({:9.3e})".format(
                          planes[k], self._alpha_rad[k]))
                return False

        # Longitudinal Twiss parameters.

        alpha_z = \
            -self._A_rad[ind.ct][ind.ct]*self._A_rad[ind.delta][ind.ct] \
            - self._A_rad[ind.ct][ind.delta] \
            *self._A_rad[ind.delta][ind.delta];
        beta_z = \
            self._A_rad[ind.ct][ind.ct]**2 \
            + self._A_rad[ind.ct][ind.delta]**2;
        gamma_z = (1e0+alpha_z**2)/beta_z;

        # Bunch size.
        self._sigma_s = np.sqrt(beta_z*self._eps[ind.Z]);
        self._sigma_delta = np.sqrt(gamma_z*self._eps[ind.Z]);

        logger.info(
            "\nE [GeV]     = {:3.1f}\nU0 [keV]    = {:3.1f}\neps         ="
            " {:12.6e} {:12.6e} {:12.6e}\ntau [msec]  = {:8.6f} {:8.6f} {:8.6f}"
            "\nJ           = {:8.6f} {:8.6f} {:8.6f}\nalpha_rad   = {:13.6e}"
            " {:13.6e} {:13.6e}\nD_rad       = {:12.6e} {:12.6e} {:12.6e}".
            format(1e-9*self._model_state.Energy, 1e-3*self._U_0,
                   self._eps[ind.X], self._eps[ind.Y], self._eps[ind.Z],
                   1e3*self._tau[ind.X], 1e3*self._tau[ind.Y],
                   1e3*self._tau[ind.Z], self._J[ind.X], self._J[ind.Y],
                   self._J[ind.Z], self._alpha_rad[ind.X],
                   self._alpha_rad[ind.Y], self._alpha_rad[ind.Z],
                   self._D_rad[ind.X], self._D_rad[ind.Y], self._D_rad[ind.Z]))

        return True

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
        self._M_rad = r.one_turn_map

        logger.info(
            "\nM_rad:\n" + mat2txt(self._M_rad.jacobian())
            + "\n\ncod ="
            + vec2txt(np.array(
                [self._M_rad.cst().x, self._M_rad.cst().px,
                 self._M_rad.cst().y, self._M_rad.cst().py,
                 self._M_rad.cst().delta, self._M_rad.cst().ct]))
        )

        self._model_state.dE = 0e0
        ps = self._M_rad.cst()
        # dE is computed by the RF cavity propagator.
        self._lattice.propagate(self._model_state, ps)
        self._dE = self._model_state.dE

        stable, self._nu_rad, self._A_rad, A_inv, self._alpha_rad = \
            compute_M_diag(dof, self._M_rad.jacobian())

        if stable:
            self._model_state.emittance = True

            A_7x7 = np.zeros((7, 7))
            A_7x7[:6, :6] = self._A_rad
            A_7x7[6, 6] = 1e0
            A_cpy  = gtpsa.ss_vect_tpsa(self._desc, 1)
            A_cpy += self._M_rad.cst()
            A_cpy.set_jacobian(A_7x7)
            self._lattice.propagate(self._model_state, A_cpy)

            self.compute_diffusion_coefficients()

            stable_rad = self.compute_rad_prop()
        else:
            print("\ncompute_radiation: unstable")
            self._U_0 = np.nan
            self._J = np.zeros(3, float)
            self._tau = np.zeros(3, float)
            self._eps = np.zeros(3, float)
            return stable, False

        return stable, stable_rad

    def prt_rad(self):
        print("\nRadiation Properties:")
        print("  E [GeV]       = {:5.3f}".
              format(1e-9*self._model_state.Energy))
        print("  U_0 [keV]     = {:5.1f}".
              format(1e-3*self._U_0))
        print("  eps_x [m.rad] = [{:9.3e}, {:9.3e}, {:9.3e}]".format(
            self._eps[ind.X], self._eps[ind.Y], self._eps[ind.Z]))
        print("  sigma_s [mm]  = {:5.3f}".format(1e3*self._sigma_s))
        print("  sigma_delta]  = {:9.3e}".format(self._sigma_delta))
        print("  J             = [{:5.3f}, {:5.3f}, {:5.3f}]".format(
            self._J[ind.X], self._J[ind.Y], self._J[ind.Z]))
        print("  tau [msec]    = [{:5.3f}, {:5.3f}, {:5.3f}]".format(
            1e3*self._tau[ind.X], 1e3*self._tau[ind.Y], 1e3*self._tau[ind.Z]))
        print("  alpha_rad     = [{:9.3e}, {:9.3e}, {:9.3e}]".format(
            self._alpha_rad[ind.X], self._alpha_rad[ind.Y],
            self._alpha_rad[ind.Z]))
        print("  D             = [{:11.5e}, {:11.5e}, {:11.5e}]".format(
            self._D_rad[ind.X], self._D_rad[ind.Y], self._D_rad[ind.Z]))

    def prt_M_rad(self):
        print("\nM:_rad", self._M_rad, end="")


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


__all__ = [lattice_properties_class]
