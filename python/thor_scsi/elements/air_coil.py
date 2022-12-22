"""
"""
import thor_scsi.lib as tslib
import numpy as np
import logging

logger = logging.getLogger("thor-scsi")

mu0 = 4 * np.pi * 1e-7


class AirCoilMagneticField(tslib.Field2DInterpolation):
    """Magnetic field created by an 2D air coil magnet


    Warning:
        If a function defined as a virtual function does not
        exist it will crash the interpreter ...


    Precomputes the factor

    .. math::
       \\frac{I \\mu_0}{2 \\pi}
    """

    def __init__(self, *, positions, currents):
        """

        Args:
            positions: transversal position of the individual wires
            currents:  current flowing through the wires
        """
        tslib.Field2DInterpolation.__init__(self)

        self.positions = positions
        self.currents = currents

        #  \mu_0
        # ---------
        #  2*pi

        factor = mu0 / (2 * np.pi)

        self.precomp = factor * self.currents  # * lp_inv

    def field_py(self, pos, field):
        """Calculate field created by the wires

        Internally calculating in the complex plane
        """
        pos = np.asarray(pos)
        pos0 = np.array([self.positions.real, self.positions.imag])
        dpos = pos[:, np.newaxis] - pos0

        # Compute distance to different points
        dpos * dpos
        dz = dpos
        r = np.absolute(dz)
        phi = np.angle(dz)
        Bphi = self.precomp * 1 / r

        # Bphi 90 degree to radius vector
        B = Bphi * np.exp((phi + np.pi / 2) * 1j)
        B = B.sum()

        field[0] = B.imag
        field[1] = B.real
        logger.debug(
            "%s: pos %s -> field %s <- %s", self.__class__.__name__, pos, field, B
        )

    def gradient_py(self, pos, gradient):
        raise NotImplementedError("gradient not (yet) implemented")

    def field(self, x, y, Bx, By):
        raise NotImplementedError("Don't know how to handle this function here ")

    def gradient(self, x, y, Bx, By):
        raise NotImplementedError("Gradient not implemented")

    def __repr__(self):
        cls_name = self.__class__.__name__
        pos = self.positions
        cur = self.currents
        return f"{cls_name}(positions={pos}, currents={cur})"

    def __str__(self):
        return self.__repr__()


class NonlinearKickerField(AirCoilMagneticField):
    """Nonlinear kicker field

    Field created by a classical telephone transmission cable

    Here the interest is within the conductors while the classical's
    application focus was on minimising the field on the outside
    (thus reducing cross talk to other conversations)

    See

    "DEVELOPMENT OF A NON-LINEAR KICKER SYSTEM TO FACILITATE
    A NEW INJECTION SCHEME FOR THE BESSY II STORAGE RING"

    T. Atkinson, M. Dirsat, O. Dressler, P. Kuske, H. Rast,
    Proceedings of IPAC2011, San Sebastián, Spain, THPO024
    2011
    """

    def __init__(self, *, position, current):
        pos = position
        #
        positions = np.array(
            [
                # right upper
                pos,
                # right lower
                pos.conjugate(),
                # left upper
                -pos.conjugate(),
                # left lower
                -pos,
            ]
        )
        currents = np.array([current] * len(positions))
        currents *= [1, -1, -1, 1]
        AirCoilMagneticField.__init__(self, positions=positions, currents=currents)


__all__ = ["AirCoilMagneticField", "NonlinearKickerField"]
