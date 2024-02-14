'''
Author:

  Johan Bengtsson
  13/02/24

Interface module for the EPICS Python interface for/to:

  Convert from Engineering to Physics units & vice versa.

  Map the magnetic/multipole families to the power supply circuits.

  Epics I/O - get/put - to the EPICS control system via the Python object
  interface.
'''


on_line = False

import os
import numpy as np

if on_line:
    from epics import caget, caput, cainfo

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
    radiate as rad, closed_orbit as co

from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt


# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

# Configuration space coordinates.
X_, Y_, Z_ = [
    ts.spatial_index.X,
    ts.spatial_index.Y,
    ts.spatial_index.Z
]
# Phase-space coordinates.
[x_, px_, y_, py_, ct_, delta_] = [
    ts.phase_space_index_internal.x,
    ts.phase_space_index_internal.px,
    ts.phase_space_index_internal.y,
    ts.phase_space_index_internal.py,
    ts.phase_space_index_internal.ct,
    ts.phase_space_index_internal.delta,
]

#-------------------------------------------------------------------------------

class middle_layer:
    def __init__(self):
        self.fam = {"quad": [
            "Q1PDR",   "Q1PTR",   "Q2PDR",   "Q2PTR",   "Q3PD1R",  "Q3PD2R",
            "Q3PD3R",  "Q3PD4R",  "Q3PD5R",  "Q3PD6R",  "Q3PD7R",  "Q3PD8R",
            "Q3P1T1R", "Q3P2T1R", "Q3PT2R",  "Q3PT3R",  "Q3PT4R",  "Q3PT5R",
            "Q3P1T6R", "Q3P2T6R", "Q3PT7R",  "Q3P1T8R", "Q3P2T8R", "Q4PD1R",
            "Q4PD2R",  "Q4PD3R",  "Q4PD4R",  "Q4PD5R",  "Q4PD6R",  "Q4PD7R",
            "Q4PD8R",  "Q4P1T1R", "Q4P2T1R", "Q4PT2R",  "Q4PT3R",  "Q4PT4R",
            "Q4PT5R",  "Q4P1T6R", "Q4P2T6R", "Q4PT7R",  "Q4P1T8R", "Q4P2T8R",
            "Q5P1T1R", "Q5P2T1R", "Q5PT2R",  "Q5PT3R",  "Q5PT4R",  "Q5PT5R",
            "Q5P1T6R", "Q5P2T6R", "Q5PT7R",  "Q5P1T8R", "Q5P2T8R", "PQIPT6"
        ]}

        self.pwr_supp = {"quad": {
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

        self.sp = {"quad_nom": {
            "Q1PDR":    2.44045585,
            "Q1PTR":    2.44045585,
            "Q2PDR":   -1.8536747,
            "Q2PTR":   -1.8536747,
            "Q3PD1R":  -2.02322285,
            "Q3PD2R":  -2.12441276,
            "Q3PD3R":  -2.12608143,
            "Q3PD4R":  -2.1282793,
            "Q3PD5R":  -2.1211438,
            "Q3PD6R":  -2.11223413,
            "Q3PD7R":  -2.11883984,
            "Q3PD8R":  -2.13404738,
            "Q3P1T1R": -2.50764398,
            "Q3P2T1R": -2.46682595,
            "Q3PT2R":  -2.45526041,
            "Q3PT3R":  -2.43119165,
            "Q3PT4R":  -2.44037407,
            "Q3PT5R":  -2.44818682,
            "Q3P1T6R": -2.69386876,
            "Q3P2T6R": -2.32789462,
            "Q3PT7R":  -2.43591598,
            "Q3P1T8R": -2.47132446,
            "Q3P2T8R": -2.51228342,
            "Q4PD1R":   1.40046286,
            "Q4PD2R":   1.4802205,
            "Q4PD3R":   1.48692991,
            "Q4PD4R":   1.4883633,
            "Q4PD5R":   1.48010676,
            "Q4PD6R":   1.48545637,
            "Q4PD7R":   1.47643621,
            "Q4PD8R":   1.49055699,
            "Q4P1T1R":  2.63205252,
            "Q4P2T1R":  2.56505973,
            "Q4PT2R":   2.57722952,
            "Q4PT3R":   2.57917393,
            "Q4PT4R":   2.58038995,
            "Q4PT5R":   2.57768425,
            "Q4P1T6R":  2.25837798,
            "Q4P2T6R":  2.55873747,
            "Q4PT7R":   2.58020271,
            "Q4P1T8R":  2.56384946,
            "Q4P2T8R":  2.64079373,
            "Q5P1T1R": -2.52154146,
            "Q5P2T1R": -2.51058167,
            "Q5PT2R":  -2.5831049,
            "Q5PT3R":  -2.62044465,
            "Q5PT4R":  -2.59546801,
            "Q5PT5R":  -2.58439541,
            "Q5P1T6R": -1.09078314,
            "Q5P2T6R": -2.42521942,
            "Q5PT7R":  -2.60426005,
            "Q5P1T8R": -2.50807154,
            "Q5P2T8R": -2.50807154,
            "PQIPT6":  -1.08082489
        }}

        self.sp["quad_low_alpha"] = {
            "Q1PDR":    2.527701,
            "Q1PTR":    2.527701,

            "Q2PDR":   -1.995355,
            "Q2PTR":   -1.995355,

            "Q3PD1R":  -1.979674,
            "Q3PD2R":  -1.979674,
            "Q3PD3R":  -1.979674,
            "Q3PD4R":  -1.979674,
            "Q3PD5R":  -1.979674,
            "Q3PD6R":  -1.979674,
            "Q3PD7R":  -1.979674,
            "Q3PD8R":  -1.979674,
            "Q3P1T1R": -1.978542,
            "Q3P2T1R": -1.978542,
            "Q3PT2R":  -1.978542,
            "Q3PT3R":  -1.978542,
            "Q3PT4R":  -1.978542,
            "Q3PT5R":  -1.978542,
            "Q3P1T6R": -1.978542,
            "Q3P2T6R": -1.978542,
            "Q3PT7R":  -1.978542,
            "Q3P1T8R": -1.978542,
            "Q3P2T8R": -1.978542,

            "Q4PD1R":   1.343589,
            "Q4PD2R":   1.343589,
            "Q4PD3R":   1.343589,
            "Q4PD4R":   1.343589,
            "Q4PD5R":   1.343589,
            "Q4PD6R":   1.343589,
            "Q4PD7R":   1.343589,
            "Q4PD8R":   1.343589,
            "Q4P1T1R":  1.343675,
            "Q4P2T1R":  1.343675,
            "Q4PT2R":   1.343675,
            "Q4PT3R":   1.343675,
            "Q4PT4R":   1.343675,
            "Q4PT5R":   1.343675,
            "Q4P1T6R":  1.343675,
            "Q4P2T6R":  1.343675,
            "Q4PT7R":   1.343675,
            "Q4P1T8R":  1.343675,
            "Q4P2T8R":  1.343675,

            "Q5P1T1R":  0.0,
            "Q5P2T1R":  0.0,
            "Q5PT2R":   0.0,
            "Q5PT3R":   0.0,
            "Q5PT4R":   0.0,
            "Q5PT5R":   0.0,
            "Q5P1T6R":  0.0,
            "Q5P2T6R":  0.0,
            "Q5PT7R":   0.0,
            "Q5P1T8R":  0.0,
            "Q5P2T8R":  0.0,
            "PQIPT6":   0.0
        }

        self.fam["sext"] = [
            "S1PR",  "S2PDR",  "S2PTR",  "S3PDR",   "S3PTR",   "S4PDR",
            "S4PTR", "S3PD1R", "S4PD1R", "S3P1T6R", "S3P2T6R", "S4P1T6R",
            "S4P2T6R"
        ]

        self.pwr_supp["sext"] = {
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

        self.sp["sext_nom"] = {
            "S1PR":     53.71159807/2,

            "S2PDR":   -44.968873/2,
            "S2PTR":   -44.968873/2,

            "S3PDR":   -47.03/2,

            "S3PTR":   -52.2/2,

            "S4PDR":    42.31/2,

            "S4PTR":    64.87/2,

            "S3PD1R":  -29.55/2,

            "S4PD1R":   28.02/2,

            "S3P1T6R": -52.2/2,
            "S3P2T6R": -52.2/2,

            "S4P1T6R":  64.87/2,
            "S4P2T6R":  64.87/2
        }

        self.conv_fact = {"quad":{}, "sext": {}}
        self.pv_get = {"quad":{}, "sext": {}}
        self.pv_put = {"quad":{}, "sext": {}}


    def middle_layer_init(self, file_name):
        # self.rd_conv_fact("sext", file_name["quad"])
        self.rd_conv_fact("sext", file_name["sext"])

    def rd_conv_fact(self, type, file_name):
        # Unit is [1/A] for thick sextupoles.
        with open(file_name) as f:
            for line in f:
                if line[0] != "#":
                    token = line.split()
                    self.conv_fact[type].update({token[0] : float(token[1])})

    def epics_init(self, type):
        prt = False
        n_prt = 5
        if prt:
            print("\nepics_init")
        for fam in self.fam[type]:
            if prt:
                print("  Family -", fam)
                n = 0
                print("  ", end="")
            for pwr_supp in self.pwr_supp[type][fam]:
                if prt:
                    print(f"  {pwr_supp:10s}", end="")
                    n += 1
                if on_line:
                    self.pv_get[type][pwr_supp] = epics.PV(pwr_supp+":get")
                    self.pv_put[type][pwr_supp] = epics.PV(pwr_supp+":set")
                else:
                    self.pv_get[type][pwr_supp] = pwr_supp+":get"
                    self.pv_put[type][pwr_supp] = pwr_supp+":set"
                if prt and (n % n_prt == 0):
                    print()
                    print("  ", end="")
                    n = 0
            if prt:
                print()

    def get_pv(self, type, pv):
        if on_line:
            val = self.pv_get[type][pv].get()
            print(f"  {pv:10s} {val:8.5f}")
        else:
            val = self.pv_get[type][pv]
            print(f"  {pv:10s} {val:14s}")

    def put_pv(self, type, pv):
        if on_line:
            val = self.pv_put[type][pv].get()
            print(f"  {pv:10s} {val:8.5f}")
        else:
            val = self.pv_put[type][pv]
            print(f"  {pv:10s} {val:14s}")

    def get_pv_fam(self, type, fam):
        print("\nget_pv -", fam, ":")
        for pwr_supp in self.pwr_supp[type][fam]:
            if on_line:
                val = self.pv_get[type][pwr_supp].get()
                print(f"  {pwr_supp:10s} {val:8.5f}")
            else:
                val = self.pv_get[type][pwr_supp]
                print(f"  {pwr_supp:10s} {val:14s}")

    def put_pv_fam(self, type, fam):
        print("\nput_pv -", fam, ":")
        for pwr_supp in self.pwr_supp[type][fam]:
            if on_line:
                val = self.pv_put[type][pwr_supp].get()
                print(f"  {pwr_supp:10s} {val:8.5f}")
            else:
                val = self.pv_put[type][pwr_supp]
                print(f"  {pwr_supp:10s} {val:14s}")

    def print_fam(self, type):
        n_prt = 5
        print(f"\nFamilies - {type:s}:")
        n = 0
        for fam in self.fam[type]:
            print(f"  {fam:10s}", end="")
            n += 1
            if n % n_prt == 0:
                print()
                n = 0
        print()

    def print_pwr_supp(self, type, fam):
        n_prt = 5
        print(f"\nPower Supplies - {fam:s}:")
        n = 0
        for pwr_supp in self.pwr_supp[type][fam]:
            print(f"  {pwr_supp:10s}", end="")
            n += 1
            if n % n_prt == 0:
                print()
                n = 0
        print()

    def print_conv_fact(self, type):
        n_prt = 5
        print(f"\nConversion Factors - {type:s}:")
        n = 0
        for pwr_supp in self.conv_fact[type]:
            print(f"  {pwr_supp:10s} {self.conv_fact[type][pwr_supp]:7.5f}",
                  end="")
            n += 1
            if n % n_prt == 0:
                print()
                n = 0
        print()

#-------------------------------------------------------------------------------

def read_lattice(file_name):
    # Read in & parse lattice file.
    lat = accelerator_from_config(file_name)

    # Set lattice state (Rf cavity on/off, etc.)
    model_state = ts.ConfigType()

    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    return n_dof, lat, model_state


def print_Twiss_param(str, Twiss):
    # eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    # that way I also check that Twiss has exactly three parameters
    eta, alpha, beta = Twiss
    print(str, end="")
    print(f"  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
    print(f"  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
    print(f"  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")


def compute_periodic_solution(lat, model_state, named_index, desc):
    """
    Todo:
        model_state: rename to calculation_configuration or calc_config
    """
    # Compute the periodic solution for a super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    stable, M, A = lo.compute_map_and_diag(n_dof, lat, model_state, desc=desc)
    print("\nM:\n" + mat2txt(M.jacobian()[:6, :6]))
    res = cs.compute_Twiss_A(A)
    Twiss = res[:3]
    print_Twiss_param("\nTwiss:\n", Twiss)
    A_map = gtpsa.ss_vect_tpsa(desc, no)
    A_map.set_jacobian(A)
    ds = \
        lo.compute_Twiss_along_lattice(
            n_dof, lat, model_state, A=A_map, desc=desc, mapping=named_index)

    return M, A, ds

#-------------------------------------------------------------------------------

# Main program.

# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 1
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

# Beam Energy [eV].
E_0 = 1.7e9

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

home_dir = os.path.join(os.environ["HOME"])

#-------------------------------------------------------------------------------
#
# Lattice.

file_name_lat = \
    os.path.join \
    (home_dir, "Nextcloud", "thor_scsi", "JB", "BESSY-II",
     "b2_stduser_beamports_blm_tracy_corr.lat")

# Read in lattice and compute the linear optics.

n_dof, lat, model_state = read_lattice(file_name_lat)

if False:
    M, A, data = \
        compute_periodic_solution(lat, model_state, named_index, desc)

#-------------------------------------------------------------------------------
#
# Middle Layer.

file_name_conv_coeff = {}
file_name_conv_coeff["quad"] = \
    os.path.join(home_dir, "Teresia/conversion-factors-quadrupoles.csv")
file_name_conv_coeff["sext"] = \
    os.path.join(home_dir, "Teresia/conversion-factors-sextupoles.txt")

# Initialise the Middle Layer.

ml = middle_layer()

ml.middle_layer_init(file_name_conv_coeff)

if not False:
    ml.print_fam("quad")
    ml.print_fam("sext")
    ml.print_pwr_supp("sext", ml.fam["sext"][0])
    ml.print_conv_fact("sext")
    print()

# Test of EPICS I/O.

ml.epics_init("sext")

if False:
    for fam in ml.fam["sext"]:
        ml.get_pv_fam("sext", fam)

if not False:
    ml.get_pv("sext", "S1MD1R")
    ml.put_pv("sext", "S1MD1R")
    ml.get_pv("sext", "S1MD1R")

if False:
    ml.get_pv_fam("sext", "S1PR")
    ml.put_pv_fam("sext", "S1PR")
    ml.get_pv_fam("sext", "S1PR")

if False:
    for fam in ml.fam["sext"]:
        ml.get_pv_fam("sext", fam)
