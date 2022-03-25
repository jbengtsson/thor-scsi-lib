"""Reading a lattice file
"""

import os.path
import sys
t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

from thor_scsi.flame import GLPSParser
import thor_scsi
from thor_scsi.lib import Accelerator, ConfigType, ss_vect_tps_type
p = GLPSParser()
with open(t_file) as fp:
    text = fp.read()
C = p.parse_byte(text, t_dir)
m = Accelerator(C)
# print(m)
print("accessing element")
sys.stdout.flush()
elem =  m.find("b_t", 0)


calc_config = ConfigType()
ps = ss_vect_tps_type()
ps.set_identity()
m.propagate(calc_config, ps, 0, 5)
print(ps)
sys.stdout.flush()
