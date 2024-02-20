
from thor_scsi.utils.bessy_ii_mml import middle_layer

on_line       = False
rd_only       = not True
physics_units = True

if not on_line:
    home_dir = "/Users/johan/Teresia"
else:
    home_dir = \
        "/net/nfs/srv/MachinePhysics/MachineDevelopment/mcateer/Jupyter/" \
        "olsson-test"

file_name_conv_coeff = {}
file_name_conv_coeff["quad"] = \
    home_dir + "/conversion-factors-quadrupoles.csv"
file_name_conv_coeff["sext"] = \
    home_dir + "/conversion-factors-sextupoles.txt"

# Allocate & initialise the Middle Layer.

mml = middle_layer()

mml.middle_layer_init(file_name_conv_coeff, on_line, rd_only, physics_units)

if not False:
    mml.prt_pwr_supp("sext")
    mml.prt_conv_fact("sext")


pwr_supp = "S3P1T6R"
init_sp = mml.get_pv_sp("sext", pwr_supp)

scaling = 1.01
mml.put_pv_sp("sext", pwr_supp, init_sp*scaling)

mml.get_pv_rb("sext", pwr_supp)

mml.put_pv_sp("sext", pwr_supp, init_sp)
