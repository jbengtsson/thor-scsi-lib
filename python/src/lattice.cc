#include <pybind11/stl.h>
#include "thor_py.h"
#include <pybind11/pybind11.h>
#include <thor_scsi/core/lattice.h>

namespace tsc = thor_scsi::core;
namespace py = pybind11;

class LatticeTypePy : public tsc::LatticeType {

public:
  std::tuple<double, double, double, std::vector<double>,
	     std::vector<double>,std::vector<double>>
  get_eps_x_py(const bool prt) {
    double eps_x, sigma_delta, U_0;
    // TODO: allocate appropriate space
    std::vector<double> J(3, 0.0), tau(3, 0.0), I(6, 0.0);

    get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, prt);

    std::tuple<double, double, double, std::vector<double>,
	       std::vector<double>,std::vector<double>>
      tup(make_tuple(eps_x, sigma_delta, U_0, J, tau, I));
    return tup;
  }
};

void py_thor_scsi_init_lattice(py::module_ &m)
{

  	py::class_<LatticeTypePy>(m, "LatticeType")
		// std::vector<ElemFamType>
		.def_readwrite("elemf", &LatticeTypePy::elemf)
		// std::vector<ElemType*>
		.def_readwrite("elems", &LatticeTypePy::elems)
		// ConfigType
		.def_readwrite("conf",  &LatticeTypePy::conf)
		.def(py::init<>())
		.def("Lat_Init",      &LatticeTypePy::Lat_Init)
		.def("prt_fams",      &LatticeTypePy::prt_fams)
		.def("prt_elems",     &LatticeTypePy::prt_elems)
		.def("GetnKid",       &LatticeTypePy::GetnKid)
		.def("ElemIndex",     &LatticeTypePy::ElemIndex)
		.def("Elem_GetPos",   &LatticeTypePy::Elem_GetPos)
		.def("Lat_Read",      &LatticeTypePy::Lat_Read)
		.def("prtmfile",      &LatticeTypePy::prtmfile)
		.def("rdmfile",       &LatticeTypePy::rdmfile)
		.def("prt_lat1",      &LatticeTypePy::prt_lat1)
		.def("prt_lat2",      &LatticeTypePy::prt_lat2)
		.def("prt_lat3",      &LatticeTypePy::prt_lat3)
		.def("prt_lat4",      &LatticeTypePy::prt_lat4)
		.def("getcod",        &LatticeTypePy::getcod)
		.def("Ring_GetTwiss", &LatticeTypePy::Ring_GetTwiss)
		.def("ChamberOff",    &LatticeTypePy::ChamberOff)
		.def("print",         &LatticeTypePy::print)
		.def("get_eps_x",     &LatticeTypePy::get_eps_x_py)
		.def("GetEmittance",  &LatticeTypePy::GetEmittance)
		.def("Cell_Pass1",    &LatticeTypePy::Cell_Pass1)
		.def("Cell_Pass2",    &LatticeTypePy::Cell_Pass2);


}
