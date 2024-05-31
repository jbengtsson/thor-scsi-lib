"""Use Case:
     Translate Tracy-2 lattice to AT by utilising PyAt.
"""


import os
import at


def tracy2at(file_name, harmonic_number):
    try:
        ring = at.load_tracy(file_name+".lat", harmonic_number=harmonic_number)
    except:
        print("\nload_tracy failed")

    try:
        at.save_m(ring, filename=file_name+".m")
    except:
        print("\nsave_m failed")


home_dir = os.path.join(
    os.environ["HOME"], os.environ["TRACY_LIB"], "projects", "in", "lattice",
    "MAX_4U")
lat_name = "max_iv_ref"
file_name = os.path.join(home_dir, lat_name)

tracy2at(file_name, 176)
