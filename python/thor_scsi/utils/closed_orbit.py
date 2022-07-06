import thor_scsi.lib as tslib
from .phase_space_vector import map2numpy
import numpy as np
import logging
from dataclasses import dataclass

logger = logging.getLogger("thor-scsi-lib")


@dataclass
class ClosedOrbitResult:

    #: was a closed orbit found?
    found_closed_orbit: bool

    #: start vector for found closed orbit
    x0: tslib.ss_vect_double

    # one turn map for the closed orbit
    one_turn_map: tslib.ss_vect_tps

    #: last element traced to
    last_element: int


def compute_closed_orbit(
    acc: tslib.Accelerator,
    conf: tslib.ConfigType,
    delta: float,
    max_iter: int = 10,
    eps: float = 1e-10,
) -> ClosedOrbitResult:
    """searches for the closed orbit

    Args:
        max_iter: maximum iteration steps
        delta:

    returns a :class:`ClosedOrbitResult` object. Its member
    :any:`ClosedOrbitResult.found_closed_orbit` is set to true
    if a result is found

    Todo:
        * trace down why delta is distributed to different places
        * x0 as part of the options?
        *
    """

    logger.debug("computing closed orbit")
    conf.dPparticle = delta

    if conf.Cavity_on:
        n = 6
    else:
        n = 4

    x0 = tslib.ss_vect_double()
    if n == 4:
        x0.set_zero()
        x0[tslib.phase_space_index_internal.delta] = delta
    elif n == 6:
        # To be revisited ...
        # if delta != 0 add eta * delta
        x0.set_zero()
        x0[tslib.phase_space_index_internal.delta] = delta
    else:
        raise AssertionError("Only implemented for 4D or 6D phase space")

    logger.debug("x0 %s", x0)

    # create weighting matrix for inverse calculation
    jj = np.zeros(tslib.ss_dim, np.int)
    for k in range(tslib.ss_dim):
        jj[k] = 1 if k < n else 0

    I = tslib.ss_vect_tps()
    I.set_identity()
    M = tslib.ss_vect_tps()

    dx_abs = 1e30

    closed_orbit = False
    n_elements = len(acc)
    # Newton's method for root finding
    for n_iter in range(max_iter):
        if dx_abs < eps:
            closed_orbit = True
            break

        # prepare return map
        M.set_identity()
        M += x0

        next_element = acc.propagate(conf, M)

        if next_element == n_elements:
            # Managed to get around the ring ... good
            # Compute difference of last point to first point
            x1 = M.cst()
            dx = x0 - x1

            # gradient search
            tmp = M - I - x1
            # jj is used to flag if longitudinal calculations are r
            # required i.e. the cavity is on
            gradient = tslib.partialInverse(tmp, jj)
            dx0 = gradient * dx

            # Next start point following line search ?
            x0 += dx0.cst()
            dx_abs = tslib.xabs(n, dx)
        else:
            dx_abs = np.nan
            break

        logger.debug(
            "n_iter %3d, dx_abs %7.1e, eps  %7.1e, x0 %s", n_iter, dx_abs, eps, x0
        )

    else:
        raise AssertionError("Exceeded maximum iterations f{max_iter}")

    # has a closed orbit been reached
    if closed_orbit:
        assert dx_abs < eps

    if closed_orbit:
        t_map = map2numpy(M)
        logger.debug(f"Poincaré Map {M}")
        result = ClosedOrbitResult(
            found_closed_orbit=closed_orbit,
            x0=x0,
            one_turn_map=t_map,
            last_element=next_element,
        )

        # Propagate with the fixed point around the lattice
        # so that observers can memorize this orbit
        acc.propagate(conf, M)

    else:
        result = ClosedOrbitResult(
            found_closed_orbit=False,
            x0=x0,
            one_turn_map=t_map,
            last_element=next_element,
        )
        elem = acc[next_element]
        ob = elem.observer()
        xkm1 = None
        if ob:
            try:
                xkm1 = ob.getPhaseSpace()
            except:
                logger.error(f"Could not retrieve phase space from observer {ob}")
        logger.error(
            f"compute_closed_orbit: failed to converge after {max_iter:d}"
            f"  delta = {delta:12.5e}, particle lost at element {next_element:3d}"
            f"  x_0   = {x0}  x_k-1 = {xkm1}"
            f"  x_k   = {M.cst():13.5e}"
        )

    return result


__all__ = ["ClosedOrbitResult", "compute_closed_orbit"]
