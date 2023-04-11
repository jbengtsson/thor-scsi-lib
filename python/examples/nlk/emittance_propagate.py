"""Simple propagation of emittance

A look to non linear kicker
"""
from pathlib import Path
from thor_scsi.factory import accelerator_from_config
from thor_scsi.pyflame import Config
import thor_scsi.lib as tslib

import numpy as np
import gtpsa
import os

x_, px_ = 0, 2

# t_file = Path(os.environ["HOME"]) / "Devel"/ "gitlab" / "dt4cc"/"lattices" / "b2_stduser_beamports_blm_tracy_corr.lat"
t_dir =  Path(os.environ["HOME"]) / "Devel" / "gitlab" / "dt4acc" / "lattices"
t_file = t_dir / "b2_stduser_beamports_blm_tracy_corr_with_nlk.lat"

def create_nlk_interpolation(nlk_name):
    def compute_mirror_position_plate(ref_pos, mirror_pos, *, y_plane=True):
        """ """
        assert y_plane
        dy = ref_pos.imag - mirror_pos.imag
        return ref_pos - 2 * dy * 1j

    # fmt: off
    ref_pos1 =  8e-3 +  7e-3j
    ref_pos2 = 17e-3 + 15e-3j
    # fmt: on
    t_current = -7e2

    # fmt: off
    t_current *= 1 - 1 * 0.14 / 2
    ref_pos1  *= 1 - 0.14
    ref_pos2  *= 1 - 0.14

    plate_position1 = 5e-3j
    mirror_pos1 = compute_mirror_position_plate(ref_pos1, plate_position1)
    print(f"{ref_pos1*1e3=} {mirror_pos1*1e3}")

    inner = tslib.aircoil_filament(ref_pos1.real, ref_pos1.imag,  t_current)
    outer = tslib.aircoil_filament(ref_pos2.real, ref_pos2.imag, -t_current)
    mirror = tslib.aircoil_filament(mirror_pos1.real, mirror_pos1.imag, -t_current * 0.14)
    nlkf_intp = tslib.NonLinearKickerInterpolation([inner, outer, mirror])

    c = Config()
    c.setAny("L", 0e0)
    c.setAny("name", nlk_name)
    c.setAny("N", 1)
    nlk = tslib.FieldKick(c)
    nlk.set_field_interpolator(nlkf_intp)
    return nlk, nlkf_intp

# Setup function start
print(f"Reading lattice file {t_file}")
acc_orig = accelerator_from_config(t_file)
calc_config = tslib.ConfigType()

elem_names = [elem.name for elem in acc_orig]


# thor scsi non linear kicker names
nlk_name = "pkdnl1kr"
nlk_index = elem_names.index(nlk_name)
nlk = acc_orig[nlk_index]
print(nlk)
nlk.get_field_interpolator().set_scale(1.0)

nv = 6
mo = 1
desc = gtpsa.desc(nv, 1)
ps = gtpsa.ss_vect_tpsa(desc, mo, nv)
ps.set_identity()
print("bessy ii standard machine: start")
print(ps)
acc_orig.propagate(calc_config, ps)
print("propagated\n", ps)
print()

nlkfk, nlkf_intp = create_nlk_interpolation(nlk.name)
assert(nlkfk.name == nlk_name)
# Use set_scale to set it to zero for any following turns ...
nlkf_intp.set_scale(1.0)

# place the thin kicker in the middle of the lattice
C = Config()
C.setAny("name", nlk_name + "_drift_upstream")
C.setAny("L", nlk.get_length()/2.0)
C_upstream = C
C = Config()
C.setAny("name", nlk_name + "_drift_downstream")
C.setAny("L", nlk.get_length()/2.0)
C_downstream = C
del C

# Split up the drift and insert the kicker in the middle of it
nlkf_elements = [tslib.Drift(C_upstream), nlkfk, tslib.Drift(C_downstream)]
print([repr(elem) for elem in nlkf_elements])
# construct the non linear kicker
print(nlk_name, nlk, nlk.get_length())

# make elements a list, these are easier to manipulate with
elements_orig = [elem for elem in acc_orig]
print("Stopping before nlk?", elements_orig[:nlk_index][-1])
print("Starting after  nlk?", elements_orig[nlk_index+1:][0])
print("Elements around nlk", [elem.name for elem in elements_orig[nlk_index - 2: nlk_index + 2]])

elements = elements_orig[:nlk_index] + nlkf_elements + elements_orig[nlk_index + 1:]
elem_names = [elem.name for elem in elements]
nlk_index = elem_names.index(nlk_name)

# round the ring ...
# do it twice once with the kicker on and once without it
acc = tslib.Accelerator(elements)

length_orig = np.sum([elem.get_length() for elem in acc_orig])
length = np.sum([elem.get_length() for elem in acc])
print(f"{length_orig=}, {length=}")

# print(" ".join(elem_names))
mu_x = 1e-3
mu_px = 1e-4

# 70 nm rad
emittance = 70e-9
sigma_x = .25e-9
sigma_px = np.sqrt(emittance**2 - sigma_x**2 )
print(f"{sigma_x=} {sigma_px=}")

nv = 6
mo = 1
desc = gtpsa.desc(nv, 1)
ps = gtpsa.ss_vect_tpsa(desc, mo, nv)
ps.set_identity()

t_x  = ps[x_]
t_x.set(0, mu_x)
t_px = ps[px_]
t_px.set(0, mu_px)
ps[x_] = t_x
ps[px_] = t_px

print(ps)
ps_orig = ps.copy()
print("setup\n", ps_orig)

nlkfk.propagate(calc_config, ps)
print("nlkfk\n", ps)

ps = ps_orig.copy()
ps = gtpsa.ss_vect_tpsa(desc, mo, nv)
ps.set_identity()
acc.propagate(calc_config, ps)
around_ring = ps.copy()

print("maps: nlk on")
print(around_ring)

ps = ps_orig.copy()
# Set kicker off
nlkf_intp.set_scale(0e0)
ps = gtpsa.ss_vect_tpsa(desc, mo, nv)
ps.set_identity()
acc.propagate(calc_config, ps)
print("maps: nlk off")
print(ps)


ps = ps_orig.copy()
print(ps)
ps = gtpsa.ss_vect_tpsa(desc, mo, nv)
ps.set_identity()
acc_orig.propagate(calc_config, ps)
print("bessy ii standard machine")
print(ps)


jac = around_ring.jacobian()
jac_x = jac[:2, :2]
sigma_matrix_x = np.array([
    [sigma_x**2, sigma_x * sigma_px],
    [sigma_x * sigma_px, sigma_px**2],
])

print(sigma_matrix_x)
print(jac_x @ sigma_matrix_x @ jac_x.T)

# Schneide bei 15 mm

# print(f"emittance {emittance*1e9:3f} urad, sigma_x {sigma_x*1e6:3f}, sigma_px {sigma_px*1e6:3f}")
