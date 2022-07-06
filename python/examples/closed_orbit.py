"""Compute closed orbit for a given lattice
"""
import thor_scsi.lib as tslib
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.closed_orbit import compute_closed_orbit
from thor_scsi.utils.output import mat2txt, vec2txt

import os

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)
conf = tslib.ConfigType()

result = compute_closed_orbit(acc, conf, 0e0, max_iter=10, eps=1e-10)
print("Unperturbated orbit:")
print("closed orbit solution:", result.found_closed_orbit)
print("start vector",  vec2txt(result.x0))
print(f"one turn matrix\n{mat2txt(result.one_turn_map)}")
print("-----------------------")

elem = acc.find("uq1", 1)
muls = elem.getMultipoles()
muls.setMultipole(1, 1e-3 - 1e-3j)
print("Perturbated orbit:")
print("Changed element", repr(elem))
result = compute_closed_orbit(acc, conf, 0e0, max_iter=10, eps=1e-10)
print("Perturbated orbit:")
print("closed orbit solution:", result.found_closed_orbit)
print("start vector",  vec2txt(result.x0))
print(f"one turn matrix\n{mat2txt(result.one_turn_map)}")
print("-----------------------")
