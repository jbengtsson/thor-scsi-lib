"""Reading a lattice file
"""

from thor_scsi.flame import GLPSParser
import thor_scsi
from thor_scsi.lib import (Accelerator, ConfigType, ss_vect_tps, ss_vect_double,
                           ss_vect_tps_to_lists)
from thor_scsi.utils import ss_vect_to_masked_matrix
import numpy as np
import os.path

#  Find the file
t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

print(f'Reading lattice file {t_file}')
# I had issues reading the lattice file directly thus I preseent
# the string stream to the parser
with open(t_file) as fp:
    text = fp.read()

# Create parser, and create configuration
C = GLPSParser().parse_byte(text, t_dir)

# The machine creates it then based on the configuration
m = Accelerator(C)

print(f"An accelerator with {len(m):d} elements")
i = 42
elem = m[i]
print(f"Element at index {i:03d} = {elem}")
# In this manner one would access an element per name.
# if of the same name, use the number to find which one it is
elem =  m.find("b2", 0)
# Element will print out itself
print(elem)
# or represent itself
print(repr(elem))
# it provides the following methods
print(dir(elem))

# Access to the multipoles
muls = elem.getMultipoles()
gradient = muls.getMultipole(2)
coeffs = muls.getCoefficients()
txt = f"""------------------------------------------------------------------
Investigating an element: should be implemented by a visitor
Element
     name             {elem.name}
     index in lattice {elem.index}
     type             {type(elem)}

Geometry:
    length            {elem.getLength(): 8.4f}
    curvature         {elem.getCurvature(): 8.4f}

    curved trajctory  {elem.assumingCurvedTrajectory():}

Magnetic field:
    main multipole
         number       {elem.getMainMultipoleNumber():2d}
         strength     {elem.getMainMultipoleStrength(): 8.4f}

    multipoles
         number       {len(coeffs)}
         gradient     {gradient}
         higher order {coeffs[2:]}


Angles:                                        even internally managed
                                               in degrees ?
    bending:          {elem.getBendingAngle(): 8.4f}
    entrance:         {elem.getEntranceAngle(): 8.4f}
    exit:             {elem.getExitAngle(): 8.4f}


Propagation
   Thick Lens:        {elem.isThick()}
   Integration method {elem.getIntegrationMethod()}
   Integration steps  {elem.getNumberOfIntegrationSteps()}
   for thick lenses

------------------------------------------------------------------
"""
print(txt)



# Prepare for calculation ... need the calculation configuration
# from thor_scsi not from flame
calc_config = ConfigType()
# create a phase space. It depends on the chosen type which
# 'kind' of calculation will be executed.
ps = ss_vect_tps()
ps.set_identity()

print(ps)
m.propagate(calc_config, ps, 0, 2000)
print(ps)


tmp = ss_vect_tps_to_lists(ps)
res = ss_vect_to_masked_matrix(tmp)
np.set_printoptions(precision=4)
print(res.data)
print(res)
