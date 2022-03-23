from thor_scsi.flame import Config
from thor_scsi.lib import Drift

C = Config()
C.setAny("L", 2)
C.setAny("name", "a_drift")
drift = Drift(C)
print(drift)
print(repr(drift))
print(drift.propagate.__doc__)
print(dir(drift))
