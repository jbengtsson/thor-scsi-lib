import thor.lib as scsi

file_name = "tests/lattices/b2_stduser_beamports_blm_tracy_corr"
lat = scsi.LatticeType()
lat.Lat_Read(file_name)
lat.Lat_Init()

for cnt, elem in enumerate(lat.elems):
    print(f"{cnt:04d} {elem}")
