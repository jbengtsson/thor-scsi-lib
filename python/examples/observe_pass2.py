"""
"""
import gtpsa
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.extract_info import accelerator_info
import thor_scsi.lib as tslib
from thor_scsi.utils.output import mat2txt, vec2txt
import matplotlib.pyplot as plt

import xarray as xr
import os.path
from thor_scsi.utils.accelerator import (
    instrument_with_standard_observers,
    extract_orbit_from_standard_observers,
    extract_orbit_from_accelerator_with_standard_observers,
)


# TME lattice
t_file = os.path.join(os.path.dirname(__file__), "lattices", "tme.lat")
acc = accelerator_from_config(t_file)
for elem in acc:
    print(elem)

calc_config = tslib.ConfigType()
ass = instrument_with_standard_observers(acc)

# Add a little distortion to the lattice..
elem = acc.find("QF", 0)
muls = elem.get_multipoles()
muls.set_multipole(1, 1e-3 - 1e-3j)


# Start the calculation ... with the start
ps = gtpsa.ss_vect_double(0e0)
ps.set_zero()
if True:
    # Closed orbit for the TME lattice with
    # calculated by closed orbit
    ps[0] = -7.007412e-04
    ps[1] = 6.045201e-03
    ps[2] = -4.268650e-03
    ps[3] = 2.542382e-04
elif True:
    # Closed orbit for the TME lattice with reverse bend as
    # calculated by closed orbit
    ps[0] = -4.978941e-04
    ps[1] = 2.217995e-03
    ps[2] = -6.858292e-03
    ps[3] = 1.835468e-04
else:
    pass


print("Start vector {}".format(vec2txt(ps)))
acc.propagate(calc_config, ps)
print("End vector   {}".format(vec2txt(ps)))
print("End vector should quite match start vector")


# Following code should be matched to
# :func:`extract_orbit_from_accelerator_with_standard_observers`
ps_tmp = [elem.get_observer().get_phase_space() for elem in acc]
phase_space_coords_names = ["x", "px", "y", "py", "delta", "ct"]
indices = [elem.index for elem in acc]
ps_xr = xr.DataArray(
    data=ps_tmp,
    name="ps",
    dims=["index", "phase_coordinate"],
    coords=[indices, phase_space_coords_names],
)

info = accelerator_info(acc)
orb = info.merge(dict(ps=ps_xr))
# print(orb)

transverse_scale = 1000

fig, ax = plt.subplots(1, 1, sharex=True)
ax.set_xlabel("s [m]")
ax.set_ylabel("x,y [mm]")
for plane in "x", "y":
    ax.plot(
        orb.s, orb.ps.sel(phase_coordinate=plane) * transverse_scale, "-", label=plane
    )
ax.legend()
plt.show()
