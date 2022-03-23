"""Reading a lattice file
"""

import os.path
t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

from thor_scsi.flame import GLPSParser
import thor_scsi
from thor_scsi.lib import Accelerator
p = GLPSParser()
C = p.parse_file(t_file, False)
print(Accelerator)

m = Accelerator(C, elements=None)
