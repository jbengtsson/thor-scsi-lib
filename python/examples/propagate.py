import thor_scsi.lib as tslib
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.accelerator import instrument_with_standard_observers
import time
import os

# t_file = os.path.join("lattices", "tme.lat")
# t_file = os.path.join("lattices", "tme_rb.lat")

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)

energy = 2.5e9

calc_config = tslib.ConfigType()
calc_config.radiation = True
# is this used anywhere?
calc_config.emittance = False
calc_config.Cavity_on = True
calc_config.Energy = energy

rad_del = instrument_with_radiators(acc, energy=energy)
if False:
    ps = tslib.ss_vect_double()
    ps.set_zero()
else:
    ps = tslib.ss_vect_tps()
    ps.set_identity()

print(ps.cst())
print(ps)
acc.propagate(conf, ps)
print(ps.cst())
print(ps)

# Need to keep reference here: otherwise pybind11 will
# remove them ... I assume
observers = instrument_with_standard_observers(acc)

# # Warm up
# for i in range(1):
#    acc.propagate(conf, ps)

## # Measurement
## n_turns = 200 #* 1000
## tic = time.clock_gettime_ns(0)
## for i in range(n_turns):
##     acc.propagate(conf, ps)
## tac = time.clock_gettime_ns(0)
## dt = (tac - tic) / 1e6
## print(f"Wall time for {n_turns}: {dt} ms ")
## tic = time.clock_gettime_ns(0)
## acc.propagate(conf, ps, n_turns=n_turns * 10000)
## tac = time.clock_gettime_ns(0)
## dt = (tac - tic) / 1e6
## print(f"Wall time for {n_turns}: {dt} ms ")
