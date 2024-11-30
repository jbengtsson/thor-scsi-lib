"""Use Case - Gaellic Football:
   Apollo 13

   Jim Lovell – 𝐺𝑒𝑛𝑡𝑙𝑒𝑚𝑒𝑛, 𝑖𝑡'𝑠 𝑏𝑒𝑒𝑛 𝑎 𝑝𝑟𝑖𝑣𝑖𝑙𝑒𝑔𝑒 𝑓𝑙𝑦𝑖𝑛𝑔 𝑤𝑖𝑡ℎ 𝑦𝑜𝑢!

   https://youtu.be/aozWYg8xbeU?t=118

   Gene Krantz – 𝐹𝑎𝑖𝑙𝑢𝑟𝑒 𝑖𝑠 𝑛𝑜𝑡 𝑎𝑛 𝑜𝑝𝑡𝑖𝑜𝑛!

   https://youtu.be/Tid44iy6Rjs?t=108

   Gaellic Football – 𝐼𝑓 𝑡ℎ𝑒𝑟𝑒'𝑠 𝑎 𝑑𝑟𝑎𝑤, 𝑦𝑜𝑢'𝑙𝑙 ℎ𝑎𝑣𝑒 𝑡𝑜 𝑟𝑒𝑝𝑙𝑎𝑦 𝑖𝑡?!

   https://youtu.be/6MPukDeUcCQ?t=135
"""

import sympy as smp

x = smp.symbols("x", real=True)
f = smp.sin(x)**3 * smp.exp(-5*x)

F = smp.integrate(f, x)
print("\nF = ", F)
