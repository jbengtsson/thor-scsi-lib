"""Reading a lattice file
"""

from thor_scsi.lib import (
    ConfigType,
)
from thor_scsi.factory import accelerator_from_config
import gtpsa
import numpy as np
import os.path

#  Find the file
t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

print(f"Reading lattice file {t_file}")
m = accelerator_from_config(t_file)

print(f"An accelerator with {len(m):d} elements")
i = 42
elem = m[i]
print(f"Element at index {i:03d} = {elem}")
# In this manner one would access an element per name.
# if of the same name, use the number to find which one it is
elem = m.find("b2", 0)
# Element will print out itself
print(elem)
# or represent itself
print(repr(elem))
# it provides the following methods
print(dir(elem))

# Access to the multipoles
muls = elem.get_multipoles()
gradient = muls.get_multipole(2)
coeffs = muls.get_coefficients()
txt = f"""------------------------------------------------------------------
Investigating an element: should be implemented by a visitor
Element
     name             {elem.name}
     index in lattice {elem.index}
     type             {type(elem)}

Geometry:
    length            {elem.get_length(): 8.4f}
    curvature         {elem.get_curvature(): 8.4f}

    curved trajctory  {elem.assuming_curved_trajectory():}

Magnetic field:
    main multipole
         number       {elem.get_main_multipole_number():2d}
         strength     {elem.get_main_multipole_strength(): 8.4f}

    multipoles
         number       {len(coeffs)}
         gradient     {gradient}
         higher order {coeffs[2:]}


Angles:                                        even internally managed
                                               in degrees ?
    bending:          {elem.get_bending_angle(): 8.4f}
    entrance:         {elem.get_entrance_angle(): 8.4f}
    exit:             {elem.get_exit_angle(): 8.4f}


Propagation
   Thick Lens:        {elem.is_thick()}
   Integration method {elem.get_integration_method()}
   Integration steps  {elem.get_number_of_integration_steps()}
   for thick lenses

------------------------------------------------------------------
"""
print(txt)


# Prepare for calculation ... need the calculation configuration
# from thor_scsi not from flame
calc_config = ConfigType()
# create a phase space. It depends on the chosen type which
# 'kind' of calculation will be executed.
desc = gtpsa.desc(6,2)
ps = gtpsa.ss_vect_tpsa(desc, 1)
ps.set_identity()

print(ps)
m.propagate(calc_config, ps, 0, 2000)
print(ps)

res = np.array(ps.jacobian())
np.set_printoptions(precision=4)
print(res)

separator = "-" * 66
print(separator)
# Access to some type
type_name = "Quadrupole"
print(f"All elements of type name '{type_name}' ")
for elem in m.elements_with_name_type(type_name):
    print(elem)
print(separator)

elem_name = "chv"
print(separator)
print(f"All elements with name {elem_name} ")
for elem in m.elements_with_name(elem_name):
    print(repr(elem))
print(separator)
