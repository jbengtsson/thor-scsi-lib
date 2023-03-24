import thor_scsi.lib as tslib
from thor_scsi.elements.air_coil import NonlinearKickerField, AirCoilMagneticField
from thor_scsi.pyflame import Config
import numpy as np
import matplotlib.pyplot as plt
from numpy import fft


def plot_field():
    ref_pos = 20e-3 + 10e-3j
    t_current = 7e2

    nlkf = NonlinearKickerField(position=ref_pos, current=t_current)
    left = AirCoilMagneticField(
        positions=np.array([ref_pos, ref_pos.conjugate()]),
        currents=np.array([-t_current, t_current]),
    )
    right = AirCoilMagneticField(
        positions=np.array([-ref_pos.conjugate(), -ref_pos]),
        currents=np.array([t_current, -t_current]),
    )

    x = np.linspace(-30e-3, 30e-3, num=300)
    pos = np.zeros([len(x), 2], dtype=np.float)
    pos[:, 0] = x

    fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
    ax_x, ax_y = axes

    B = np.zeros([len(x), 2], dtype=np.float)
    xs = x * 1e3

    [nlkf.field_py(tp, tB) for tp, tB in zip(pos, B)]
    B = B * 1e3
    ax_x.plot(xs, B[:, 0], "-", label="nlk")
    ax_y.plot(xs, B[:, 1], "-", label="nlk")
    ax_x.set_xlabel("x [mm]")
    ax_y.set_xlabel("y [mm]")
    ax_x.set_ylabel("B$_x$ [mT]")
    ax_y.set_ylabel("B$_y$ [mT]")

    B = np.zeros([len(x), 2], dtype=np.float)
    [left.field_py(tp, tB) for tp, tB in zip(pos, B)]
    B = B * 1e3
    ax_x.plot(xs, B[:, 0], "--", label="left")
    ax_y.plot(xs, B[:, 1], "--", label="left")

    B = np.zeros([len(x), 2], dtype=np.float)
    [right.field_py(tp, tB) for tp, tB in zip(pos, B)]
    B = B * 1e3
    ax_x.plot(xs, B[:, 0], "-.", label="right")
    ax_y.plot(xs, B[:, 1], "-.", label="right")

    fig.savefig("nlk_field_x.pdf")


def plot_field_circle(fignum=None):
    ref_pos = 20e-3 + 10e-3j
    t_current = 7e2

    a = 30e-3
    b = 6e-3

    nlkf = NonlinearKickerField(position=ref_pos, current=t_current)
    # for fft
    N = 256 * 16

    tphi = np.arange(N)
    phi = tphi * 2 * np.pi / N
    pphi = phi / np.pi * 180.0

    z = a * np.cos(phi) + 1j * b * np.sin(phi)
    pos = np.zeros([len(z), 2], dtype=float)
    pos[:, 0] = z.real
    pos[:, 1] = z.imag

    scale = 1000
    fig, ax = plt.subplots(1, 1, figsize=[4, 2], num=fignum, sharex=True, sharey=True)
    ax.plot(nlkf.positions.real * scale, nlkf.positions.imag * scale, "o")
    ax.plot(pos[:, 0] * scale, pos[:, 1] * scale, "-")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.axis("equal")
    fig.tight_layout()
    fig.savefig("nlk_geometry.pdf")

    fignum2 = fignum if fignum is None else fignum + 1
    fig, axes = plt.subplots(1, 2, num=fignum2, figsize=[8/2, 6/2], sharex=True, sharey=True)
    ax_x, ax_y = axes

    B = np.zeros([len(z), 2], dtype=float)

    [nlkf.field_py(tp, tB) for tp, tB in zip(pos, B)]
    B = B * 1e3
    ax_x.plot(pphi, B[:, 0], "-", label="nlk")
    ax_y.plot(pphi, B[:, 1], "-", label="nlk")
    ax_x.set_xticks([0, 180, 360])
    ax_y.set_xticks([0, 180, 360])
    ax_x.set_xlabel("$\phi$ [deg]")
    ax_y.set_xlabel("$\phi$ [deg]")
    ax_x.set_ylabel("B$_x$ [mT]")
    ax_y.set_ylabel("B$_y$ [mT]")
    fig.tight_layout()
    fig.savefig("nlk_field_ellipse.pdf")

    # Elliptical harmonic multipoles
    Bz = B[:, 1] + B[:, 0] * 1j
    spec = fft.fft(Bz)

    fignum3 = fignum if fignum is None else fignum + 2
    fig, axes = plt.subplots(1, 2, num=fignum3, sharex=True, sharey=True)
    ax_x, ax_y = axes

    ax_x.plot(spec.real)
    ax_y.plot(spec.imag)
    ax_x.set_xlabel("N")
    ax_y.set_xlabel("N")
    ax_x.set_ylabel("Re[spec]")
    ax_x.set_ylabel("Im[spec]")
    return Bz, spec


def main():
    ref_pos = 20e-3 + 10e-3j
    t_current = 1e3

    nlkf = NonlinearKickerField(position=ref_pos, current=t_current)
    print(nlkf)

    C = Config()
    C.setAny("L", 0e0)
    C.setAny("name", "nlk")
    C.setAny("N", 1)

    nlk = tslib.FieldKick(C)
    nlk.set_field_interpolator(nlkf)

    config = tslib.ConfigType()

    use_double = False
    print(repr(nlk))

    for x_off in (0, 10e-3):
        if use_double:
            ps = tslib.ss_vect_double()
            ps.set_zero()
        else:
            ps = tslib.ss_vect_tps()
            ps.set_identity()


        ps[0] = x_off
        print(f"propagate offset x={x_off}\n", ps)
        nlk.propagate(config, ps)
        print("pos\n", ps)

    x_off = 10.0e-3
    config = tslib.ConfigType()
    nlk.propagate(config, ps)
    print("pos\n", ps)


if __name__ == "__main__":
    # plot_field()
    # plot_field_circle()
    # plt.show()
    # print(tslib)
    main()
