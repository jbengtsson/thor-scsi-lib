from dataclasses import dataclass
from typing import ClassVar

import gtpsa

@dataclass
class tpsa_class:
    # GTPSA properties.
    # Number of variables - phase-space coordinates & 1 for parameter
    #dependence.
    _nv: ClassVar[int] = 6 + 1
    # Max order.
    _no: ClassVar[int] = 1
    # Truncation order.
    _to: ClassVar[int] = 1
    # Number of parameters.
    _nv_prm: ClassVar[int] = 0
    # Parameters max order.
    _no_prm: ClassVar[int] = 0
    # Index.
    _named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]

    def tpsa(self):
        return gtpsa.tpsa(self._desc, self._no)
    def ctpsa(self):
        return gtpsa.ctpsa(self._desc, self._no)
    def ss_vect_tpsa(self):
        return gtpsa.ss_vect_tpsa(self._desc, self._no)

__all__ = [tpsa_class]
