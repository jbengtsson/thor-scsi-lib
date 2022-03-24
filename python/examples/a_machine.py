"""Reading a lattice file
"""

import os.path
t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

from thor_scsi.flame import GLPSParser
import thor_scsi
from thor_scsi.lib import Accelerator
p = GLPSParser()
with open(t_file) as fp:
    text = fp.read()
C = p.parse_byte(text, t_dir)
m = Accelerator(C)
print(m)
elem =  m.find("uq1", 0)
