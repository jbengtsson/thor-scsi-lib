"""
"""
import thor_scsi.lib as tslib
from .accelerator import instrument_with_radiators
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


    # Diffusion coefficients
    acc.propagate(calc_config, Ap)
    print("Ap after calc")
    print(Ap)

    r = RadiationResult(relaxation_constants=w.real, fractional_tunes=w.imag)
