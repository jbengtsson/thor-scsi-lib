# Author:
#
#   Johan Bengtsson
#   13/02/24
#
# Middle Layer Class.
#
# Interface module for EPICS for:
#
#   Converting from Engineering to Physics units & vice versa.
#
#   Mapping the magnetic multipole families to the power supply circuits.
#
#   Epics I/O - get/put - to the EPICS control system via the Python object
#   interface.

#-------------------------------------------------------------------------------

class middle_layer:
    # Private.

    def __init__(self):
        self._on_line       = False
        self._rd_only       = True
        self._physics_units = True
        
        self._pwr_supp = {"quad": {
            "Q1PDR": [
                "Q1M1D1R", "Q1M2D1R", "Q1M1D2R", "Q1M2D2R", "Q1M1D3R",
                "Q1M2D3R", "Q1M1D4R", "Q1M2D4R", "Q1M1D5R", "Q1M2D5R",
                "Q1M1D6R", "Q1M2D6R", "Q1M1D7R", "Q1M2D7R", "Q1M1D8R",
                "Q1M2D8R"
            ],
            "Q1PTR": [
                "Q1M1T1R", "Q1M2T1R", "Q1M1T2R", "Q1M2T2R", "Q1M1T3R",
                "Q1M2T3R", "Q1M1T4R", "Q1M2T4R", "Q1M1T5R", "Q1M2T5R",
                "Q1M1T6R", "Q1M2T6R", "Q1M1T7R", "Q1M2T7R", "Q1M1T8R",
                "Q1M2T8R"
            ],
            "Q2PDR": [
                "Q2M1D1R", "Q2M2D1R", "Q2M1D2R", "Q2M2D2R", "Q2M1D3R",
                "Q2M2D3R", "Q2M1D4R", "Q2M2D4R", "Q2M1D5R", "Q2M2D5R",
                "Q2M1D6R", "Q2M2D6R", "Q2M1D7R", "Q2M2D7R", "Q2M1D8R",
                "Q2M2D8R"
            ],
            "Q2PTR": [
                "Q2M1T1R", "Q2M2T1R", "Q2M1T2R", "Q2M2T2R", "Q2M1T3R",
                "Q2M2T3R", "Q2M1T4R", "Q2M2T4R", "Q2M1T5R", "Q2M2T5R",
                "Q2M1T6R", "Q2M2T6R", "Q2M1T7R", "Q2M2T7R", "Q2M1T8R",
                "Q2M2T8R"],
            "Q3PD1R": [
                "Q3M1D1R", "Q3M2D1R"],
            "Q3PD2R": [
                "Q3M1D2R", "Q3M2D2R"],
            "Q3PD3R": [
                "Q3M1D3R", "Q3M2D3R"],
            "Q3PD4R": [
                "Q3M1D4R", "Q3M2D4R"],
            "Q3PD5R": [
                "Q3M1D5R", "Q3M2D5R"],
            "Q3PD6R": [
                "Q3M1D6R", "Q3M2D6R"],
            "Q3PD7R": [
                "Q3M1D7R", "Q3M2D7R"],
            "Q3PD8R": [
                "Q3M1D8R", "Q3M2D8R"],
            "Q3P1T1R": [
                "Q3M1T1R"],
            "Q3P2T1R": [
                "Q3M2T1R"],
            "Q3PT2R": [
                "Q3M1T2R", "Q3M2T2R"],
            "Q3PT3R": [
                "Q3M1T3R", "Q3M2T3R"],
            "Q3PT4R": [
                "Q3M1T4R", "Q3M2T4R"],
            "Q3PT5R": [
                "Q3M1T5R", "Q3M2T5R"],
            "Q3P1T6R": [
                "Q3M1T6R"],
            "Q3P2T6R": [
                "Q3M2T6R"],
            "Q3PT7R": [
                "Q3M1T7R", "Q3M2T7R"],
            "Q3P1T8R": [
                "Q3M1T8R"],
            "Q3P2T8R": [
                "Q3M2T8R"],
            "Q4PD1R": [
                "Q4M1D1R", "Q4M2D1R"],
            "Q4PD2R": [
                "Q4M1D2R", "Q4M2D2R"],
            "Q4PD3R": [
                "Q4M1D3R", "Q4M2D3R"],
            "Q4PD4R": [
                "Q4M1D4R", "Q4M2D4R"],
            "Q4PD5R": [
                "Q4M1D5R", "Q4M2D5R"],
            "Q4PD6R": [
                "Q4M1D6R", "Q4M2D6R"],
            "Q4PD7R": [
                "Q4M1D7R", "Q4M2D7R"],
            "Q4PD8R": [
                "Q4M1D8R", "Q4M2D8R"],
            "Q4P1T1R": [
                "Q4M1T1R"],
            "Q4P2T1R": [
                "Q4M2T1R"],
            "Q4PT2R": [
                "Q4M1T2R", "Q4M2T2R"],
            "Q4PT3R": [
                "Q4M1T3R", "Q4M2T3R"],
            "Q4PT4R": [
                "Q4M1T4R", "Q4M2T4R"],
            "Q4PT5R": [
                "Q4M1T5R", "Q4M2T5R"],
            "Q4P1T6R": [
                "Q4M1T6R"],
            "Q4P2T6R": [
                "Q4M2T6R"],
            "Q4PT7R": [
                "Q4M1T7R", "Q4M2T7R"],
            "Q4P1T8R": [
                "Q4M1T8R"],
            "Q4P2T8R": [
                "Q4M2T8R"],
            "Q5P1T1R": [
                "Q5M1T1R"],
            "Q5P2T1R": [
                "Q5M2T1R"],
            "Q5PT2R": [
                "Q5M1T2R", "Q5M2T2R"],
            "Q5PT3R": [
                "Q5M1T3R", "Q5M2T3R"],
            "Q5PT4R": [
                "Q5M1T4R", "Q5M2T4R"],
            "Q5PT5R": [
                "Q5M1T5R", "Q5M2T5R"],
            "Q5P1T6R": [
                "Q5M1T6R"],
            "Q5P2T6R": [
                "Q5M2T6R"] ,
            "Q5PT7R": [
                "Q5M1T7R", "Q5M2T7R"],
            "Q5P1T8R": [
                "Q5M1T8R"],
            "Q5P2T8R": [
                "Q5M2T8R"],
            "PQIPT6": [
                "QIT6"]
        }}

        self._pwr_supp["sext"] = {
            "S1PR": [
                "S1MD1R", "S1MT1R", "S1MD2R", "S1MT2R", "S1MD3R", "S1MT3R",
                "S1MD4R", "S1MT4R", "S1MD5R", "S1MT5R", "S1MD6R", "S1MT6R",
                "S1MD7R", "S1MT7R", "S1MD8R", "S1MT8R"
            ],
            "S2PDR": [
                "S2M1D1R", "S2M2D1R", "S2M1D2R", "S2M2D2R", "S2M1D3R",
                "S2M2D3R", "S2M1D4R", "S2M2D4R", "S2M1D5R", "S2M2D5R",
                "S2M1D6R", "S2M2D6R", "S2M1D7R", "S2M2D7R", "S2M1D8R",
                "S2M2D8R"
            ],
            "S2PTR": [
                "S2M1T1R", "S2M2T1R", "S2M1T2R", "S2M2T2R", "S2M1T3R",
                "S2M2T3R", "S2M1T4R", "S2M2T4R", "S2M1T5R", "S2M2T5R",
                "S2M1T6R", "S2M2T6R", "S2M1T7R", "S2M2T7R", "S2M1T8R",
                "S2M2T8R"
            ],
            "S3PDR": [
                "S3M1D2R", "S3M2D2R", "S3M1D3R", "S3M2D3R", "S3M1D4R",
                "S3M2D4R", "S3M1D5R", "S3M2D5R", "S3M1D6R", "S3M2D6R",
                "S3M1D7R", "S3M2D7R", "S3M1D8R", "S3M2D8R"
            ],
            "S3PTR": [
                "S3M1T1R", "S3M2T1R", "S3M1T2R", "S3M2T2R", "S3M1T3R",
                "S3M2T3R", "S3M1T4R", "S3M2T4R", "S3M1T5R", "S3M2T5R",
                "S3M1T7R", "S3M2T7R", "S3M1T8R", "S3M2T8R"
            ],
            "S4PDR": [
                "S4M1D2R", "S4M2D2R", "S4M1D3R", "S4M2D3R", "S4M1D4R",
                "S4M2D4R", "S4M1D5R", "S4M2D5R", "S4M1D6R", "S4M2D6R",
                "S4M1D7R", "S4M2D7R", "S4M1D8R", "S4M2D8R"],
            "S4PTR": [
                "S4M1T1R", "S4M2T1R", "S4M1T2R", "S4M2T2R", "S4M1T3R",
                "S4M2T3R", "S4M1T4R", "S4M2T4R", "S4M1T5R", "S4M2T5R",
                "S4M1T7R", "S4M2T7R", "S4M1T8R", "S4M2T8R"],
            "S3PD1R": [
                "S3M1D1R", "S3M2D1R"],
            "S4PD1R": [
                "S4M1D1R", "S4M2D1R"],
            "S3P1T6R": [
                "S3M1T6R"],
            "S3P2T6R": [
                "S3M2T6R"],
            "S4P1T6R": [
                "S4M1T6R"],
            "S4P2T6R": [
                "S4M2T6R"]
        }

        # Unit is [1/m A] for thick quadrupoles.
        self._conv_fact = {"quad":{}}
        # Unit is [1/m^2 A] for thick sextupoles.
        self._conv_fact = {"sext":{}}

        self._pv_sp_get = {"quad":{}, "sext": {}}
        self._pv_sp_put = {"quad":{}, "sext": {}}
        self._pv_rb_get = {"quad":{}, "sext": {}}
      
    # Public.

    def middle_layer_init(self, file_name, on_line, rd_only, physics_units):
        self._on_line = on_line
        self._rd_only = rd_only
        self._physics_units = physics_units

        if on_line:
            import epics

        self.rd_conv_fact(file_name["sext"], "sext")

        if on_line:
            self.epics_init("sext")

    def rd_conv_fact(self, file_name, type):
        with open(file_name) as f:
            for line in f:
                if line[0] != "#":
                    token = line.split()
                    self._conv_fact[type].update({token[0] : float(token[1])})

    def epics_init(self, type):
        prt = not False
        n_prt = 5
        if prt:
            print("\nepics_init")
        n = 0
        for pwr_supp in self._pwr_supp[type]:
            if prt:
                print("  {:10s}".format(pwr_supp), end="")
                n += 1
            if self._on_line:
                self._pv_sp_get[type][pwr_supp] = epics.PV(pwr_supp+":set")
                self._pv_sp_put[type][pwr_supp] = epics.PV(pwr_supp+":set")
                self._pv_rb_get[type][pwr_supp] = epics.PV(pwr_supp+":rdbk")
            else:
                self._pv_sp_get[type][pwr_supp] = pwr_supp+":set"
                self._pv_sp_put[type][pwr_supp] = pwr_supp+":set"
                self._pv_rb_get[type][pwr_supp] = pwr_supp+":rdbk"
            if prt and (n % n_prt == 0):
                print()
                n = 0
        if prt:
            print()

    def get_pv_sp(self, type, pv):
        if self._on_line:
            val = self._pv_sp_get[type][pv].get()
        else:
            val = 0.1234
        if self._physics_units:
            val *= self._conv_fact[type][pv]
        print("  {:10s} {:8.5f}".format(pv, val))
        return val

    def put_pv_sp(self, type, pv, val):
        if self._physics_units:
            val /= self._conv_fact[type][pv]
        if self._on_line:
            if not rd_only:
                self._pv_sp_put[type][pv].put(val)
            else:
                print("\nput_pv_sp – rd_only mode")
        else:
            print("  {:10s} {:8.5f}".format(pv, val))
            
    def get_pv_rb(self, type, pv):
        if self._on_line:
            val = self._pv_rb_get[type][pv].get()
        else:
            val = 0.1234
        if self._physics_units:
            val *= self._conv_fact[type][pv]
        print("  {:10s} {:8.5f}".format(pv, val))
        return val            
            
    def prt_pwr_supp(self, type):
        n_prt = 5
        print("\nPower Supplies:")
        n = 0
        for mag in self._pwr_supp[type]:
            print("  {:10s}".format(mag), end="")
            n += 1
            if n % n_prt == 0:
                print()
                n = 0
        print()

    def prt_conv_fact(self, type):
        print("\nConversion Factors:")
        for pwr_supp in self._conv_fact[type]:
            print("  {:10s} {:8.5f}".
                  format(pwr_supp, self._conv_fact[type][pwr_supp]))


__all__ = ["middle_layer_init", "rd_conv_fact", "epics_init", "get_pv_sp",
           "put_pv_sp", "get_pv_rb", "prt_pwr_supp", "prt_conv_fact"]
