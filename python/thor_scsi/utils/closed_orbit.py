import thor_scsi.lib as tslib
from .phase_space_vector import map2numpy
import numpy as np
import logging
from dataclasses import dataclass
import copy

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
    *,
    delta: float = None,
    x0: tslib.ss_vect_double = None,
    max_iter: int = 10,
    eps: float = 1e-6,
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

    if eps <= 0e0:
        raise AssertionError(f"tolerance {eps} <= 0e0")

    if x0 is not None and not np.isfinite(x0).all():
        raise AssertionError(f"given start {x0} is not finite!")

    if conf.Cavity_on:
        n = 6
    else:
        n = 4

    logger.debug(f" Cavity on ? {conf.Cavity_on} {n=}")

    if x0 is None:
        assert delta is not None
        conf.dPparticle = delta
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
    else:
        if delta is not None:
            raise AssertionError("if x0 is given delta must be None")
        conf.dPparticle = x0[tslib.phase_space_index_internal.delta]
    logger.debug("x0 %s", x0)

    # create weighting matrix for inverse calculation
    jj = np.zeros(tslib.ss_dim, np.int)

    # jj[:n] = 1  # select active phase space coordinates
    for k in range(tslib.ss_dim):
        jj[k] = 1 if k < n else 0

    I = tslib.ss_vect_tps()
    I.set_identity()
    M = tslib.ss_vect_tps()

    dx_abs = 1e30

    closed_orbit = False
    n_elements = len(acc)
    # Newton's method for root finding
    logger.warning(
        "start     , dx_abs %7.1e, eps  %7.1e, x0 %s", dx_abs, eps, x0
    )


    for n_iter in range(max_iter):
        if dx_abs < eps:
            closed_orbit = True
            break

        # prepare return map
        M.set_identity()
        M += x0
        logger.warning(f"{n_iter=},  Start propagation at \n {M.cst()}\n")
        next_element = acc.propagate(conf, M)
        logger.warning(f"{n_iter=},  End propagation at \n {M.cst()}\n")

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
            # dx_aps = np.sqrt(np.sum(dx[:n] ** 2))
            dx_abs = tslib.xabs(n, dx)

        else:
            dx_abs = np.nan
            break

        logger.warning(
            "n_iter %3d, dx_abs %7.1e, eps  %7.1e, x0 %s", n_iter + 1, dx_abs, eps, x0
        )

    else:
        logger.error(
            f"compute closed orbit did not finish within {max_iter} "
            f"iterations. x0 = {x0}, target: eps = {eps} reached {dx_abs}"
        )
        raise AssertionError(f"Exceeded maximum iterations {max_iter}")

    # has a closed orbit been reached
    if closed_orbit:
        assert dx_abs < eps
        M.set_identity()
        M += x0

        logger.debug(f" Fixed point check propagate again start : x0\n{M.cst()}\nMat {M}")
        Mref = copy.copy(M)
        acc.propagate(conf, M)
        logger.debug(f" Fixed point check propagate returned    : x0\n{M.cst()}\nMat {M}")

    if closed_orbit:
        t_map = map2numpy(M)
        logger.debug(f"Poincaré Map\n {M}")
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
