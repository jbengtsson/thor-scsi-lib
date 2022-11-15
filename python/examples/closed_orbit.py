"""Compute closed orbit for a given lattice

Computes the orbit
* for the lattice as given
* for an element with an additional dipole component

It illustrates how observers are used to obtain the
orbit data. Furthermore it is shown that observers
can be replaced by new ones. This allows inspecting
the data the observers have seen at a later stage.
"""
from thor_scsi.lib import ConfigType
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.closed_orbit import compute_closed_orbit
from thor_scsi.utils.accelerator import (
    instrument_with_standard_observers,
    extract_orbit_from_standard_observers,
    extract_orbit_from_accelerator_with_standard_observers,
)
from thor_scsi.utils.extract_info import accelerator_info
from thor_scsi.utils.output import mat2txt, vec2txt
import matplotlib.pyplot as plt
import os

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

t_file = os.path.join("lattices", "tme.lat")
# t_file = os.path.join("lattices", "tme_rb.lat")
acc = accelerator_from_config(t_file)
conf = ConfigType()

# Standard setup does not come with observers ...
# see why in a second
observers_non_perturbated = instrument_with_standard_observers(acc)
# Observers will then flag if new data has arrived
[ob.reset() for ob in observers_non_perturbated]
result = compute_closed_orbit(acc, conf, delta=0e0, max_iter=10, eps=1e-10)

fmt = (
    """Orbit: {}
closed orbit solution: {}
start vector\n{}\n
one turn matrix\n{}
"""
    + "\n"
    + "-" * 78
    + "\n"
)
print(
    fmt.format(
        "as defined by lattice",
        result.found_closed_orbit,
        vec2txt(result.x0),
        mat2txt(result.one_turn_map),
    )
)

res = extract_orbit_from_accelerator_with_standard_observers(acc)
fig, axes = plt.subplots(2, 1, sharex=True)
ax_x, ax_y = axes
ax_x.set_xlabel("s [m]")
ax_y.set_xlabel("s [m]")
ax_x.set_ylabel("x [mm]")
ax_y.set_ylabel("y [mm]")
transverse_scale = 1000
line_x, = ax_x.plot(
    res.s,
    res.ps.sel(phase_coordinate="x") * transverse_scale,
    "-",
    label="undisturbed orbit",
)
line_y, = ax_y.plot(
    res.s,
    res.ps.sel(phase_coordinate="y") * transverse_scale,
    "-",
    label="undisturbed orbit",
)


# elem = acc.find("uq1", 0)
elem = acc.find("QF", 0)
muls = elem.get_multipoles()
muls.set_multipole(1, 1e-3 - 1e-3j)
observers_perturbated = instrument_with_standard_observers(acc)
# Observers will then flag if new data has arrived.
# Not strictly necessary here as a new set or observers was
# created and registered.
[ob.reset() for ob in observers_perturbated]
result = compute_closed_orbit(acc, conf, delta=0e0, max_iter=10, eps=1e-10)

print(
    fmt.format(
        f"added dipole component to element {elem.name}, ",
        result.found_closed_orbit,
        vec2txt(result.x0),
        mat2txt(result.one_turn_map),
    )
)

# Get data again
res_dis = extract_orbit_from_accelerator_with_standard_observers(acc)
ax_x.plot(
    res_dis.s,
    res_dis.ps.sel(phase_coordinate="x") * transverse_scale,
    "-",
    label="distributed orbit",
)
ax_y.plot(
    res_dis.s,
    res_dis.ps.sel(phase_coordinate="y") * transverse_scale,
    "-",
    label="distributed orbit",
)

# Add a little artificial offset if you can not distinguish the redrawn lines
offset = 0e0
# Plot again using data from original observers
info = accelerator_info(acc)
res_rep = extract_orbit_from_standard_observers(observers_non_perturbated, info)
ax_x.plot(
    res_rep.s,
    res_rep.ps.sel(phase_coordinate="x") * transverse_scale + offset,
    "--",
    color=line_x.get_color(),
    label="undisturbed orbit (reexported)",
)
line_y, = ax_y.plot(
    res_rep.s,
    res_rep.ps.sel(phase_coordinate="y") * transverse_scale + offset,
    "--",
    color=line_y.get_color(),
    label="undisturbed orbit  (reexported)",
)

# mark the elements' position
s_pos = res_dis.s.sel(index=elem.index)
for ax in [ax_x, ax_y]:
    axis = ax.axis()
    ax.plot([s_pos, s_pos], axis[2:], "k-", linewidth=2)
    ax.axis(axis)
for ax in [ax_x, ax_y]:
    ax.legend()
plt.show()

# If you run it in e.g. ipython have fun inspecting the results
# or the observers. These should still contain the data from them
# different configurations.
