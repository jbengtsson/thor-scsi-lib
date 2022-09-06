"""Configuring and using a mpole
"""
from thor_scsi.pyflame import Config
from thor_scsi.lib import Mpole, ConfigType, Drift
from thor_scsi.lib import ss_vect_double
from flask import Flask
from flask_restful import Api, Resource , reqparse
import json

app = Flask(__name__)
api = Api(app)

#add the coming arguments + their validations
mpole_args = reqparse.RequestParser()
mpole_args.add_argument("name", type = str, help = "name is mandatory")
mpole_args.add_argument("l", type = float)
mpole_args.add_argument("n", type = int)
mpole_args.add_argument("type", type = str)
mpole_args.add_argument("x", type = float)
mpole_args.add_argument("px", type = float)
mpole_args.add_argument("y", type = float)
mpole_args.add_argument("py", type = float)
mpole_args.add_argument("dt", type = float)
mpole_args.add_argument("ct", type = float)
#get the beam positions into the args


class Propagate(Resource):
    def post(self):
        args = mpole_args.parse_args()
        print("x has a value of ")
        print(args.x)
        C = Config()
        C.setAny("L", args.l)
        C.setAny("N", args.n)
        C.setAny("name", args.name)
        # C.setAny("K", args.k)
        # what are other fields for example t , t1....
        # quadrupole = 1
        if args.type == "drift":
            mpole = Drift(C) #element instead of mpole
        else:
            mpole = Mpole(C)
            mpole.getFieldInterpolator().setMultipole(5, 4.2e-4) 
            mul = mpole.getFieldInterpolator()
            print(mul)
            print(repr(mul))
        ps = ss_vect_double()
        # ps.set_zero()
        #get the beam positions as well
        ps[0] = args.x
        ps[1] = args.px
        # ps[2] = args.y
        # ps[3] = args.py
        # ps[4] = args.dt
        # ps[5] = args.ct
        # ps[0] = 1e-3
        # ps[2] = -1e-3
        calc_config = ConfigType()
        mpole.propagate(calc_config, ps)

        # Convert to JSON: elementWithPS
        mpoleWithPhaseSpace =  '{ "x":'+str(ps[0]) +',  "px":'+str(ps[1]) +', "y":'+str(ps[2]) +', "py":'+str(ps[3]) +', "delta":'+str(ps[4]) +',"ct":'+str(ps[5]) +', "name":"'+args.name+'", "type":"' +args.type+'", "l":'+ str(args.l)+',"n":'+str(args.n)+'}'
        #print(mpoleWithPhaseSpace)
        # parse JSON:
        return json.loads(mpoleWithPhaseSpace)
api.add_resource(Propagate, "/propagate")
        

if __name__ == "__main__":
    app.run(debug=True)

