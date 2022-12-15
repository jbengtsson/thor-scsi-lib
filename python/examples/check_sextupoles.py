import os.path
from thor_scsi.factory import accelerator_from_config
import thor_scsi.lib as tslib
from typing import Sequence
import re

t_dir = os.path.join(os.environ["HOME"], "Devel", "gitlab", "dt4cc", "lattices")
t_file = os.path.join(t_dir, "b2_stduser_beamports_blm_tracy_corr.lat")

lat = accelerator_from_config(t_file)


class PositionName:
    def __init__(self, name):

        self.straight_type_ = None
        self.elem_type = None
        self.magnet_type = None
        self.is_magnet = False
        self.child_number = None

        try:
            self._processName(name)
        except:
           print(f"Failed to process name {name}")
           raise

    def _processName(self, name):
        self._processRing(name[-1])
        self._processSector(name[-2])

        self._processFamily(name[0])
        self._processChildNumber(name[1])

    def _processRing(self, character):
        assert(character == "r")

    def _processChildNumber(self, character):
        chld = int(character)
        assert (child in [12])

    def _processFamily(self, character):
        family = character
        if family == "S":
            self.magnet_type = "sextupole"
            self.is_magnet = True
        elif family == "Q":
            self.magnet_type = "Quadrupole"
            self.is_magnet = True
        elif family == "B":
            self.magnet_type = "Bending"
            self.is_magnet = True
        else:
            pass

    @property
    def elementType(self):
        return self.elem_type

    @property
    def straight_type(self):
        return self.straight_type



def s1_sextupoles(seq: Sequence):
    for elem in seq:
        if elem.name[:2] == "s1":
            yield elem


def s2pdr_sextupoles(seq: Sequence):
    for elem in seq:
        pos = PositionName(elem.name)

def check_main_multipole(seq):
    for elem in s1_sextupoles(seq):
        print(elem.name, elem.get_main_multipole_strength().real * 2)
        # return


def match_name(elem):
    r = component_match.match(elem.name)
    if r:
        # print(elem.name, r.groupdict())
        pass
    elif bpm_match.match(elem.name):
        pass
    elif bpm_match2.match(elem.name):
        pass
    elif drift_match.match(elem.name):
        pass
    else:
        print(f"No match for {elem.name}")


def main():
    # check_main_multipole(s1_sextupoles(lat))
    for cnt, elem in enumerate(lat):
        match_name(elem)
        if cnt  > 100:
            # return
            pass

    elements = list(lat)
    check_main_multipole(s2pdr_sextupoles(elements[1:]))


if __name__ == "__main__":
    main()
