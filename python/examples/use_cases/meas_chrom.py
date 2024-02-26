from thor_scsi.utils.bessy_ii_mml import middle_layer


self._on_line       = False
self._rd_only       = True

if on_line:
    from epics import caget, caput, cainfo


class epics_class:
    # Private.

    def __init__(self):
        self._on_line = False
        self._rd_only = True

    # Public.

    def mml_init(on_line, rd_only):
        self._on_line = on_line
        self._rd_only = rd_only
        

    def get_pv_sp(self, pv):
        if self._on_line: 
            val = caget(pv)
        else:
            val = 0.1234
        print("  {:10s} {:8.5f}".format(pv, val))
        return val

    def put_pv_sp(self, pv, val):
        if self._on_line:
            if not rd_only:
                caput(pv, val)
            else:
                print("\nput_pv_sp – rd_only mode")
        else:
            print("  {:10s} {:8.5f}".format(pv, val))
            
    def get_pv_rb(self, pv):
        if self._on_line:
            val = caget(pv)
        else:
            val = 0.1234
        print("  {:10s} {:8.5f}".format(pv, val))
        return val            
            

    def meas_chrom(n_point, df_RF, pv_sp_RF):
        f_0_RF = get_pv_sp(pv_sp_RF)
        print("\nmeas_chrom:")
        for k in range(-n_point, n_point+1):
            f_RF = f_0_RF + k*df_RF
            put_pv_sp(pv, f_RF)
            print("  {:2d} {:11.5e}".format(k, f_RF))


on_line = False
rd_only = not True

if not on_line:
    home_dir = "/Users/johan/Teresia"
else:
    home_dir = \
        "/net/nfs/srv/MachinePhysics/MachineDevelopment/mcateer/Jupyter/" \
        "olsson-test"

mml = epics_class()

mml.init(on_line, rd_only)

scaling = 1.01
mml.put_pv_sp("sext", init_sp*scaling)

mml.get_pv_rb("sext", )

mml.put_pv_sp("sext", pwr_supp, init_sp)
