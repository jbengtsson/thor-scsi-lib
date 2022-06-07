import thor_scsi.lib as tslib
from thor_scsi.elements.air_coil import NonlinearKickerField, AirCoilMagneticField
from thor_scsi.pyflame import Config
import numpy as np
import matplotlib.pyplot as plt


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
    print(nlkf, left, right)
    x = np.linspace(-30e-3, 30e-3, num=300)
    pos = np.zeros([len(x), 2], dtype=np.float)
    pos[:, 0] = x

    fig, axes = plt.subplots(1, 2, sharex=True)
    ax_x, ax_y = axes
    B = np.zeros([len(x), 2], dtype=np.float)
    [nlkf.field_py(tp, tB) for tp, tB in zip(pos, B)]
    ax_x.plot(x, B[:, 0], "-", label="nlk")
    ax_y.plot(x, B[:, 1], "-", label="nlk")
    ax_x.set_xlabel("x")
    ax_y.set_xlabel("y")

    B = np.zeros([len(x), 2], dtype=np.float)
    [left.field_py(tp, tB) for tp, tB in zip(pos, B)]
    ax_x.plot(x, B[:, 0], "--", label="left")
    ax_y.plot(x, B[:, 1], "--", label="left")

    B = np.zeros([len(x), 2], dtype=np.float)
    [right.field_py(tp, tB) for tp, tB in zip(pos, B)]
    ax_x.plot(x, B[:, 0], "-.", label="right")
    ax_y.plot(x, B[:, 1], "-.", label="right")


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
    nlk.setFieldInterpolator(nlkf)

    print(repr(nlk))
    ps = tslib.ss_vect_double()
    ps.set_zero()
    config = tslib.ConfigType()
    print("propagate center")
    nlk.propagate(config, ps)
    print("pos", ps)

    ps = tslib.ss_vect_double()
    ps.set_zero()
    x_off = 10.0e-3
    ps[0] = x_off
    config = tslib.ConfigType()
    print("propagate offset x={x_off}")
    nlk.propagate(config, ps)
    print("pos", ps)


if __name__ == "__main__":
    # plot_field()
    # plt.show()
    main()
