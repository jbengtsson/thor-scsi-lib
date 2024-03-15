"""
"""
import thor_scsi.lib as tslib

import gtpsa
from .accelerator import instrument_with_radiators
import numpy as np
from scipy.constants import c as c0
from dataclasses import dataclass

from thor_scsi.utils.closed_orbit import compute_closed_orbit
from thor_scsi.utils.linear_optics import compute_M_diag
from thor_scsi.utils.output import mat2txt, vec2txt

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


def compute_circ(lat):
    return np.sum([elem.get_length() for elem in lat])


def compute_diffusion_coefficients(rad_del_kicks):
    dD_rad = \
        np.array([rk.get_diffusion_coefficients_increments()
                  for rk in rad_del_kicks])
    D_rad = np.sum(dD_rad, axis=0)
    return D_rad


def compute_rad_prop(lat, model_state, x0, dE, alpha_rad, D_rad):
    dof = 3
    J = np.zeros(dof)
    tau = np.zeros(dof)
    eps = np.zeros(dof)
    C = compute_circ(lat)
    logger.info("\nC = %5.3f", C)
    U_0 = model_state.Energy*dE
    for k in range(dof):
        J[k] = 2e0*(1e0+x0.delta)*alpha_rad[k]/dE
        tau[k] = -C/(c0*alpha_rad[k])
        eps[k] = -D_rad[k]/(2e0*alpha_rad[k])

    logger.info("\nE [GeV]     = {:3.1f}\nU0 [keV]    = {:3.1f}\neps         = {:12.6e} {:12.6e} {:12.6e}\ntau [msec]  = {:8.6f} {:8.6f} {:8.6f}\nJ           = {:8.6f} {:8.6f} {:8.6f}\nalpha_rad   = {:13.6e} {:13.6e} {:13.6e}\nD_rad       = {:12.6e} {:12.6e} {:12.6e}".
                format(1e-9*model_state.Energy,
                       1e-3*U_0,
                       eps[X_], eps[Y_], eps[Z_],
                       1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_],
                       J[X_], J[Y_], J[Z_],
                       alpha_rad[X_], alpha_rad[Y_], alpha_rad[Z_],
                       D_rad[X_], D_rad[Y_], D_rad[Z_]))

    return U_0, J, tau, eps



def compute_radiation(
    lat: tslib.Accelerator,
    model_state: tslib.ConfigType,
    E,
    eps,
    *, desc
):

    dof = 3

    model_state.Energy    = E
    model_state.radiation = True
    model_state.emittance = False
    model_state.Cavity_on = True

    # Install radiators that radiation is calculated
    rad_del_kicks = instrument_with_radiators(lat, energy=E)

    r = compute_closed_orbit(lat, model_state, delta=0e0, eps=eps)
    # M = r.one_turn_map[:6, :6]
    M = r.one_turn_map.jacobian()[:6, :6]

    logger.info(
        "\nM:\n" + mat2txt(M)
        + "\n\nx0 ="
        + vec2txt(np.array(
            [r.x0.x, r.x0.px, r.x0.y, r.x0.py, r.x0.delta, r.x0.ct]))
    )

    model_state.dE = 0e0
    ps = r.x0
    lat.propagate(model_state, ps)
    dE = model_state.dE

    stable, A, A_inv, alpha_rad = compute_M_diag(dof, M)

    A_7x7 = np.zeros((7, 7))
    A_7x7[:6, :6] = A
    A_7x7[6, 6] = 1e0
    if stable:
        model_state.emittance = True

        #A_cpy = vec_mat2ss_vect_tps(r.x0, A)
        A_cpy  = gtpsa.ss_vect_tpsa(desc, 1)
        A_cpy += r.x0
        A_cpy.set_jacobian(A_7x7)
        lat.propagate(model_state, A_cpy)

        D_rad = compute_diffusion_coefficients(rad_del_kicks)

        U_0, J, tau, eps = \
            compute_rad_prop(lat, model_state, r.x0, dE, alpha_rad, D_rad)
    else:
        U_0 = np.nan
        J = np.zeros(3, float)
        tau = np.zeros(3, float)
        eps = np.zeros(3, float)

    cod = np.array([r.x0.x, r.x0.px, r.x0.y, r.x0.py, r.x0.delta, r.x0.ct])

    # return stable, M, cod, A, U_0, J, tau, eps, D_rad
    return stable, r.one_turn_map, cod, A, U_0, J, tau, eps, D_rad


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
