#!/usr/bin/env python
"""
Compute the constants of the forth order integrator for cross check

c_1 = 1/(2*(2-2^(1/3)))
c_2 = (1-2^(1/3))/(2*(2-2^(1/3)))
d_1 = 1/(2-2^(1/3))
d_2 = -2^(1/3)/(2-2^(1/3))
"""

c_1 = 1/(2*(2-2**(1/3)))
c_2 = (1-2**(1/3))/(2*(2-2**(1/3)))
d_1 = 1/(2-2**(1/3))
d_2 = -2**(1/3)/(2-2**(1/3))

txt = f"""
Constants as computed by python

c_1 = {c_1:20.16f}
c_2 = {c_2:20.16f}

d_1 = {d_1:20.16f}
d_2 = {d_2:20.16f}
"""
print(txt)
