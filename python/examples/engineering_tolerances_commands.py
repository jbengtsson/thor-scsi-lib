"""just create a few engineering tolerances commands

"""
import numpy as np
from numpy.random import Generator, PCG64, default_rng

from thor_scsi.utils import engineering

indices = np.arange(5)

mul_prop = engineering.Property(
    element_info=None, set_method_name="setMultipoles", get_method_name="getMultipoles"
)


eng_cmd = engineering.ScaleEngineeringDistributionCommand(
    t_property=mul_prop,
    loc=0,
    size=1e-3,
    distribution_name="normal",
    vector_length=2
    # vector_length=20,
)
print(eng_cmd)

# This seems to start with a different seed every time
# rng = default_rng()

rng = Generator(PCG64(seed=355113))
dist_cmds = engineering.create_distribution_commands(indices, eng_cmd)
print("Distribution commands")
for dc in dist_cmds:
    print("\t", dc)
del dc

cmds = engineering.create_commands_for_randomised_state(dist_cmds, rng=rng)
print("EngineeringCommands for a specific set ")
for c in cmds:
    print("\t", c)
del c
