"""Compute closed orbit for a given lattice

Computes the orbit
* for the lattice as given
* for an element with an additional dipole component

"""
from thor_scsi.lib import ConfigType
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.closed_orbit import compute_closed_orbit
from thor_scsi.utils.accelerator import (
    instrument_with_standard_observers,
    extract_orbit_from_standard_observers,
)
from thor_scsi.utils.output import mat2txt, vec2txt
import matplotlib.pyplot as plt
import os

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)
# Standard setup does not come with observers ...
# see why in a second
observers = instrument_with_standard_observers(acc)

conf = ConfigType()

# Observers will then flag if new data has arrived
[ob.reset() for ob in observers]
result = compute_closed_orbit(acc, conf, 0e0, max_iter=10, eps=1e-10)
print("Unperturbated orbit:")
print("closed orbit solution:", result.found_closed_orbit)
print("start vector", vec2txt(result.x0))
print(f"one turn matrix\n{mat2txt(result.one_turn_map)}")
print("-----------------------")

res = extract_orbit_from_standard_observers(acc)
fig, axes = plt.subplots(2, 1, sharex=True)
ax_x, ax_y = axes
ax_x.set_xlabel("s [m]")
ax_y.set_xlabel("s [m]")
ax_x.set_ylabel("x [mm]")
ax_y.set_ylabel("y [mm]")
transverse_scale = 1000
ax_x.plot(
    res.s,
    res.ps.sel(phase_coordinate="x") * transverse_scale,
    label="undisturbed orbit",
)
ax_y.plot(
    res.s,
    res.ps.sel(phase_coordinate="y") * transverse_scale,
    label="undisturbed orbit",
)


elem = acc.find("uq1", 0)
muls = elem.getMultipoles()
muls.setMultipole(1, 1e-3 - 1e-3j)
print("Perturbated orbit:")
print("Changed element", repr(elem))
# Observers will then flag if new data has arrived
[ob.reset() for ob in observers]
result = compute_closed_orbit(acc, conf, 0e0, max_iter=10, eps=1e-10)
print("Perturbated orbit:")
print("closed orbit solution:", result.found_closed_orbit)
print("start vector", vec2txt(result.x0))
print(f"one turn matrix\n{mat2txt(result.one_turn_map)}")
print("-----------------------")
res = extract_orbit_from_standard_observers(acc)
ax_x.plot(
    res.s,
    res.ps.sel(phase_coordinate="x") * transverse_scale,
    label="distributed orbit",
)
ax_y.plot(
    res.s,
    res.ps.sel(phase_coordinate="y") * transverse_scale,
    label="distributed orbit",
)

# mark the elements' position
s_pos = res.s.sel(index=elem.index)
for ax in [ax_x, ax_y]:
    axis = ax.axis()
    ax.plot([s_pos, s_pos], axis[2:], 'k-', linewidth=2)
    ax.axis(axis)
for ax in [ax_x, ax_y]:
    ax.legend()
