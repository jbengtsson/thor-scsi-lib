"""Use Case:
     Parametric unit dipole cell properties vs. Reverse Bend angle.
"""


import os
import sys

from thor_scsi.utils import lattice_properties as lp, index_class as ind


def set_phi_rb(lat_prop, rb_name, phi_uc, bend_list, bend_scl, phi_rb):
    dphi = phi_uc/2e0 - phi_rb
    lat_prop.set_phi_fam(rb_name, phi_rb)
    for k in range(len(bend_scl)):
        lat_prop.set_phi_fam(bend_list[k], bend_scl[k]*dphi)


cod_eps = 1e-15
E_0     = 3.0e9

ind = ind.index_class()

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB")

file_name = os.path.join(home_dir, sys.argv[1]+".lat")

lat_prop = lp.lattice_properties_class(file_name, E_0, cod_eps, 2)

bend_list = ["d2_h2_sl_d0a", "d2_h2_sl_d0b", "d2_h2_sl_d0c", "d2_h2_sl_df1",
             "d2_h2_sl_df2", "d2_h2_sl_df3", "d2_h2_sl_df4", "d2_h2_sl_df5"]
# Reverse Bend name.
rb_name = "r2_h2"
# Location for eta_x(phi_RB) - centre of dipole unit cell.
loc = lat_prop._lattice.find("d2_h2_sl_d0a", 0).index
# Total bend angle.
phi_uc = lat_prop.compute_phi_lat()
bend_scl = lat_prop.compute_scl_fact(bend_list)

print(f"\n  phi_uc = {phi_uc:5.3f}")
print(f"  Reverse Bend name: {rb_name:s}")
print(f"  {lat_prop._lattice[loc].name:s} (loc = {loc:d})")

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncomp_per_sol: unstable")
        raise ValueError

    # Adjust RF phase for sign of alpha_c.
    cav_loc = lat_prop._lattice.find("cav", 0)
    if lat_prop._alpha_c >= 0e0:
        cav_loc.set_phase(0.0)
    else:
        print("  alpha_c = {:10.3e} phi_rf = 180 deg".
              format(lat_prop._alpha_c))
        cav_loc.set_phase(180.0)

   # Compute radiation properties.
    if not lat_prop.compute_radiation():
        print("\ncompute_radiation: unstable")
        raise ValueError
except ValueError:
    exit
else:
    lat_prop.prt_lat_param()
    lat_prop.prt_rad()
    lat_prop.prt_M()
    lat_prop.prt_M_rad()
    if not False:
        lat_prop.plt_Twiss(file_name+".png", not False)

phi, eps_x, J_x, J_z, alpha_c = \
    lat_prop.unit_cell_rev_bend(
        rb_name, loc, 25, -0.8, set_phi_rb, phi_uc, bend_list, bend_scl)

if not False:
    lat_prop.plt_scan_phi_rb(
        file_name+"_phi_rb.png", phi, eps_x, J_x, J_z, alpha_c, True)
