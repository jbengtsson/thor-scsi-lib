import thor_scsi.lib as tslib
from thor_scsi.elements.air_coil import NonlinearKickerField, AirCoilMagneticField
from thor_scsi.pyflame import Config
import numpy as np
import matplotlib.pyplot as plt
from numpy import fft


# fmt: off
ref_pos1 =  8e-3 +  7e-3j
ref_pos2 = 17e-3 + 15e-3j
# fmt: on
ref_pos = 8 + 25e-3j
t_current = -7e2

# fmt: off
t_current *= 1 - 1 * 0.14 / 2
ref_pos1  *= 1 - 0.14
ref_pos2  *= 1 - 0.14
# fmt: on


def compute_mirror_position_plate(ref_pos, mirror_pos, *, y_plane=True):
    """ """
    assert y_plane
    dy = ref_pos.imag - mirror_pos.imag
    return ref_pos - 2 * dy * 1j


plate_position1 = 5e-3j
mirror_pos1 = compute_mirror_position_plate(ref_pos1, plate_position1)
print(f"{ref_pos1*1e3=} {mirror_pos1*1e3}")


def plot_field(fignum=10):
    binary = True
    if binary:
        inner = (ref_pos1.real, ref_pos1.imag,  t_current)
        outer = (ref_pos2.real, ref_pos2.imag, -t_current)
        mirror = (mirror_pos1.real, mirror_pos1.imag, -t_current * 0.14)
        nlkf = tslib.NonlinearKicker([inner, outer, mirror])
        nlkf.set_scale(1.0)
        print("nlkf", repr(nlkf))
    else:
        nlkf = NonlinearKickerField(positions=[ref_pos1, ref_pos2], current=t_current * (1-0.18))

    nlkfwm = NonlinearKickerField(
        positions=[ref_pos1, ref_pos2],
        current=t_current * (1-0.0),
        mirror_positions=[mirror_pos1, None],
    )
    left = AirCoilMagneticField(
        positions=np.array([ref_pos, ref_pos.conjugate()]),
        currents=np.array([-t_current, t_current]),
    )
    right = AirCoilMagneticField(
        positions=np.array([-ref_pos.conjugate(), -ref_pos]),
        currents=np.array([t_current, -t_current]),
    )

    x = np.linspace(-40e-3, 40e-3, num=300)
    pos = np.zeros([len(x), 2], dtype=np.float)
    pos[:, 0] = x

    if True:
        fig, axes = plt.subplots(1, 2, sharex=True, sharey=True, num=fignum)
        ax_x, ax_y = axes
    else:
        fig = None

    B = np.zeros([len(x), 2], dtype=float)
    Bwm = np.zeros([len(x), 2], dtype=float)
    xs = x * 1e3

    def evaluate_nlkf(pos, nlk=nlkf):
        return  nlk.field(pos)
        # print(r, field)

    B = np.array([nlkf.field(tp) for tp in pos])
    [nlkfwm.field_py(tp, tB) for tp, tB in zip(pos, Bwm)]
    B = B * 1e3
    Bwm = Bwm * 1e3
    if fig:
        ax_x.plot(xs, B[:, 0], "b-", label="nlk")
        ax_y.plot(xs, B[:, 1], "b-", label="nlk")
        ax_x.plot(xs, Bwm[:, 0], "c-", label="nlk with mirror currents")
        ax_y.plot(xs, Bwm[:, 1], "c-", label="nlk with mirror currents")
        ax_x.set_xlabel("x [mm]")
        ax_y.set_xlabel("x [mm]")
        ax_x.set_ylabel("B$_x$ [mT]")
        ax_y.set_ylabel("B$_y$ [mT]")

    Bl = np.zeros([len(x), 2], dtype=np.float)
    [left.field_py(tp, tB) for tp, tB in zip(pos, Bl)]
    Bl = Bl * 1e3
    # ax_x.plot(xs, B[:, 0], "--", label="left")
    # ax_y.plot(xs, B[:, 1], "--", label="left")

    Br = np.zeros([len(x), 2], dtype=np.float)
    [right.field_py(tp, tB) for tp, tB in zip(pos, Br)]
    Br = Br * 1e3
    # ax_x.plot(xs, B[:, 0], "-.", label="right")
    # ax_y.plot(xs, B[:, 1], "-.", label="right")

    # Comparison to Marc Dirsat's field

    t, x, By = np.loadtxt("By_750ns_2d.csv", skiprows=1, delimiter=",").T
    # x in mm need to calculate with it
    x = x / 1000.0
    x = x - x.mean()

    # ploting quantities
    px = x * 1000
    pBy = By * 1000
    if fig:
        ax_y.plot(px, pBy, "g.-")

        fig.savefig("nlk_field_x.pdf")

    fig2, ax = plt.subplots(
        1, 1, figsize=[8, 6], sharex=True, sharey=True, num=fignum + 1
    )
    ax.plot(xs, B[:, 1], "b-", label="nlk")
    ax.plot(xs, Bwm[:, 1], "c:")
    ax.plot(px, pBy, "g.-")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("B$_y$ [mT]")
    ax.set_xlim(-12,12)
    fig2.savefig("nlk_field_y.pdf")

    # Evaluate field at position for difference
    pos = np.zeros([len(x), 2], dtype=np.float)
    B_ref = np.zeros([len(x), 2], dtype=np.float)
    B_ref_wm = np.zeros([len(x), 2], dtype=np.float)
    pos[:, 0] = x

    B_ref = np.array([nlkf.field(tp) for tp in pos])
    # [nlkf.field_py(tp, tB) for tp, tB in zip(pos, B_ref)]
    [nlkfwm.field_py(tp, tB) for tp, tB in zip(pos, B_ref_wm)]
    px_ref = pos[:, 0] * 1000
    pBy_ref = B_ref[:, 1] * 1000
    pBy_ref_wm = B_ref_wm[:, 1] * 1000
    if fig:
        ax_y.plot(px_ref, pBy_ref, "r--")
        ax_y.plot(px_ref, pBy_ref_wm, "m.--")

    fig3, ax = plt.subplots(1, 1, sharex=True, sharey=True, num=fignum + 2)
    ax.plot(px_ref, pBy - pBy_ref, 'g-')
    ax.plot(px_ref, pBy - pBy_ref_wm, 'm--')
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("B$_y$ [mT]")
    ax.set_xlim(-12,12)


def plot_field_circle(fignum=None):
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
    fig, axes = plt.subplots(
        1, 2, num=fignum2, figsize=[8 / 2, 6 / 2], sharex=True, sharey=True
    )
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
    plot_field()
    # plot_field_circle()
    # plt.show()
    # print(tslib)
    # main()
    plt.show()
