"""
"""
import thor_scsi.lib as tslib
from .accelerator import instrument_with_radiators
from .closed_orbit import compute_closed_orbit
from . import linalg
from . import linear_optics as lo
from .output import mat2txt
from .phase_advance import compute_nus_for_symplectic_matrix
from .phase_space_vector import ss_vect_tps2ps_jac, array2ss_vect_tps
import numpy as np
from dataclasses import dataclass

import logging
logger = logging.getLogger("thor-scsi")

@dataclass
class RadiationResult:
    #: the different taus
    relaxation_constants: np.ndarray
    emittance: np.ndarray
    #: for consistency checks
    fractional_tunes: np.ndarray


def calculate_radiation(
    acc: tslib.Accelerator,
    *,
    energy,
    calc_config: tslib.ConfigType = None,
    install_radiators: bool = True,
    dof=3
):
    """

    Todo:
        Rename function
        Inspect if radiators are installed?
        check in calc config if radation is requested
        in this case install the radiators

    1. calculate fix point and Poincarè Map M with damped system
       (i.e. radiation on and cavity on (without dispersion in a second case)
    2. diagonalise M = A $\Gamma$ A$^{-1}$
    3. eigenvalues:

        - complex part: tunes,
        - real part: damping times  (refer equation)

       use eigen values of symplectic matrix to identify the planes

    4. propagate A, thin kick will create diffusion coeffs (don't forget to
       zero them before calculation starts (sum it up afterwards

    """

    if install_radiators:
        logger.debug("Installing radiators")
        # keep variable as long as you need to do calculations
        rad_del = instrument_with_radiators(acc, energy=energy)

    if calc_config is None:
        raise AssertionError
        calc_config = tslib.ConfigType()
        calc_config.radiation = True
        # is this used anywhere?
        calc_config.emittance = False
        calc_config.Cavity_on = True

    calc_config.Energy = energy
    logger.debug(
        f"calc_config radiation { calc_config.radiation} emmittance {calc_config.emittance} Cavity on {calc_config.Cavity_on}"
    )

    # Compute the fixed point ... radation is on
    r = compute_closed_orbit(acc, calc_config, delta=0e0)
    print("r.one_turn_map")
    print(mat2txt(r.one_turn_map))
    # diagonalise M
    n = 2 * dof
    M_tp = np.transpose(r.one_turn_map[:n, :n])

    # Finds eigen values and eigen vectors ... same functionallity as in
    # find_phase_space_origin / find_phase_space_fix_point
    w, v = np.linalg.eig(M_tp)
    print("v")
    print(mat2txt(v))
    w, v = linalg.match_eigenvalues_to_plane(M_tp, w, v, n_dof=dof)
    print("v matched to planes ?")
    print(mat2txt(v))
    # print(mat2txt(M_t))
    eta = np.zeros(6, np.float)
    A_inv, v1 = linalg.compute_A_inv_prev(dof, eta, v)

    A = np.linalg.inv(A_inv)

    Atmp = np.zeros([7, 7], dtype=np.float)
    Atmp[:6, :6] = A
    print("Atmp ")
    print(mat2txt(Atmp))
    # Atmp[[2, 3, 4, 5], :] = Atmp[[4, 5, 2, 3], :]
    print("Atmp resuffled")
    print(mat2txt(Atmp))
    Ap = array2ss_vect_tps(Atmp)

    print("Ap before calc")
    print(Ap)

    # Diffusion coefficients
    acc.propagate(calc_config, Ap)
    print("Ap after calc")
    print(Ap)

    r = RadiationResult(relaxation_constants=w.real, fractional_tunes=w.imag)
