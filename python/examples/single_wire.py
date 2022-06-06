import thor_scsi.lib as tslib
from thor_scsi.elements.air_coil import AirCoilMagneticField
import numpy as np
import matplotlib.pyplot as plt


class SingleWire(AirCoilMagneticField):
    """A single wire in space

    Used to check the properties of the potential
    """

    def __init__(self, *, position, current):
        pos = position
        positions = np.array([pos])
        currents = np.array([current])
        AirCoilMagneticField.__init__(self, positions=positions, currents=currents)


# position of the wire
w_pos = (20 + 10j) * 1e-3
sgw = SingleWire(position=w_pos, current=700)
x = np.linspace(0, w_pos.real * 2, num=20 * 2 * 10 + 1)
y = np.linspace(0, w_pos.imag * 2, num=10 * 2 * 10 + 1)

X, Y = np.meshgrid(x, y)

pos = np.zeros(X.shape + (2,))
pos[..., 0] = X
pos[..., 1] = Y
B = np.zeros(X.shape + (2,))
rB = B.reshape(-1, 2)
rpos = pos.reshape(-1, 2)

[sgw.field_py(p, rb) for p, rb in zip(rpos, rB)]

print(rB[..., 0].max(), rB[..., 0].max())
print(rB[..., 1].max(), rB[..., 1].max())

B = rB.reshape(X.shape + (2,))

r = 3e-3
w_idx = ((pos[..., 0] - w_pos.real) ** 2 + (pos[..., 1] - w_pos.imag) ** 2) < r ** 2

B[w_idx, 0] = np.nan
B[w_idx, 1] = np.nan
fig, axes = plt.subplots(1, 3)
ax_x, ax_y, ax_a = axes
extent = np.array([x.min(), x.max(), y.min(), y.max()]) * 1e3
print(B[..., 0].max(), B[..., 0].max())
print(B[..., 1].max(), B[..., 1].max())
ax_x.imshow(B[..., 0], extent=extent, origin="lower")
ax_y.imshow(B[..., 1], extent=extent, origin="lower")

Ba = np.absolute(B[..., 0] + B[..., 1] * 1j)
ax_a.imshow(np.log(Ba), extent=extent, origin="lower")


phi = np.linspace(0, 2 * np.pi, num=1000)
R = r
probe = r * np.exp(phi * 1j) + w_pos
probe_pos = np.zeros((len(probe), 2))
probe_pos[:, 0] = probe.real
probe_pos[:, 1] = probe.imag
Bp = np.zeros((len(probe), 2))

pphi = phi * 180 / np.pi


fig, axes = plt.subplots(1, 3)
ax_x, ax_y, ax_a = axes

ax_x.plot(pphi, probe_pos[:, 0])
ax_y.plot(pphi, probe_pos[:, 1])
ax_a.plot(pphi, np.absolute(probe - w_pos))
[sgw.field_py(p, rb) for p, rb in zip(probe_pos, Bp)]

fig, axes = plt.subplots(1, 3)
ax_x, ax_y, ax_a = axes

Ba = Bp[:, 0] * np.sin(phi) + Bp[:, 1] * np.cos(phi)
ax_x.plot(pphi, Bp[:, 0])
ax_y.plot(pphi, Bp[:, 1])
ax_a.plot(pphi, Ba)

plt.show()
