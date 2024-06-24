"""Use Case:
     Translate Tracy-2 lattice to AT by utilising PyAt.
"""


import os
import at


def tracy2at(file_name, harmonic_number):
    ring = at.load_tracy(file_name+".lat", harmonic_number=harmonic_number)
    at.save_m(ring, filename=file_name+".m")


home_dir = os.path.join(
    os.environ["HOME"], os.environ["TRACY_LIB"], "projects", "in", "lattice",
    "MAX_4U")
lat_name = "max_iv_sp_matched_2"
file_name = os.path.join(home_dir, lat_name)

tracy2at(file_name, 176)
