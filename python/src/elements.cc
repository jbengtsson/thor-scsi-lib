#include <pybind11/stl.h>
#include "thor_py.h"
#include <pybind11/pybind11.h>
#include <thor_scsi/core/elements.h>

using namespace thor_scsi::elements;
namespace py = pybind11;

void py_thor_scsi_init_elements(py::module_ &m)
{

	py::class_<ElemType>(m, "ElemType")
		// std::string
		.def_readwrite("Name",    &ElemType::Name)
		// bool
		.def_readwrite("Reverse", &ElemType::Reverse)
		// double
		.def_readwrite("PL",      &ElemType::PL)
		// PartsKind
		.def_readwrite("Pkind",   &ElemType::Pkind)
		// Virtual class: no constructor.
		// .def(py::init<>())
		.def("prt_elem", &ElemType::prt_elem)
		.def("repr_elem", &ElemType::repr_elem)
		.def("repr", &ElemType::repr)
		.def("__repr__", &ElemType::repr)
		.def("print",    &ElemType::print);
	//

	py::class_<ElemFamType>(m, "ElemFamType")
		// ElemType
		// int
		.def_readwrite("nKid",    &ElemFamType::nKid)
		.def_readwrite("NoDBN",   &ElemFamType::NoDBN)
		// std::vector<int>
		.def_readwrite("KidList", &ElemFamType::KidList)
		// std::vector<string>
		.def_readwrite("DBNlist", &ElemFamType::DBNlist)
		.def(py::init<>());

	py::class_<DriftType, ElemType>(m, "DriftType")
		.def(py::init<>());


	py::class_<MpoleType, ElemType>(m, "MpoleType")
		// int
		.def_readwrite("Pmethod",  &MpoleType::Pmethod)
		.def_readwrite("PN",       &MpoleType::PN)
		.def_readwrite("Porder",   &MpoleType::Porder)
		.def_readwrite("n_design", &MpoleType::n_design)
		// double
		.def_readwrite("PdTpar",   &MpoleType::PdTpar)
		.def_readwrite("PdTsys",   &MpoleType::PdTsys)
		.def_readwrite("PdTrms",   &MpoleType::PdTrms)
		.def_readwrite("PdTrnd",   &MpoleType::PdTrnd)
		.def_readwrite("PTx1",     &MpoleType::PTx1)
		.def_readwrite("PTx2",     &MpoleType::PTx2)
		.def_readwrite("Pgap",     &MpoleType::Pgap)
		.def_readwrite("Pirho",    &MpoleType::Pirho)
		.def_readwrite("Pc0",      &MpoleType::Pc0)
		.def_readwrite("Pc1",      &MpoleType::Pc1)
		.def_readwrite("Ps1",      &MpoleType::Ps1)
		// std::vector<double>
		.def_readwrite("PdSsys",   &MpoleType::PdSsys)
		.def_readwrite("PdSrms",   &MpoleType::PdSrms)
		.def_readwrite("PdSrnd",   &MpoleType::PdSrnd)
		// MpoleArray
		.def_readwrite("PBpar",    &MpoleType::PBpar)
		.def_readwrite("PBsys",    &MpoleType::PBsys)
		.def_readwrite("PBrms",    &MpoleType::PBrms)
		.def_readwrite("PBrnd",    &MpoleType::PBrnd)
		.def_readwrite("PB",       &MpoleType::PB)
		// pthicktype
		.def_readwrite("Pthick",   &MpoleType::Pthick)
		// std::vector< std::vector<double> >
		.def_readwrite("M_elem",   &MpoleType::M_elem)
		.def(py::init<>());


	py::class_<CavityType,     ElemType>(m, "CavityType")
		.def(py::init<>());

	py::class_<MarkerType,     ElemType>(m, "MarkerType")
		.def(py::init<>());

	py::class_<WigglerType,    ElemType>(m, "WigglerType")
		.def(py::init<>());

	py::class_<InsertionType,  ElemType>(m, "InsertionType")
		.def(py::init<>());

	py::class_<FieldMapType,   ElemType>(m, "FieldMapType")
		.def(py::init<>());

	py::class_<SpreaderType,   ElemType>(m, "SpreaderType")
		.def(py::init<>());

	py::class_<RecombinerType, ElemType>(m, "RecombinerType")
		.def(py::init<>());

	py::class_<SolenoidType,   ElemType>(m, "SolenoidType")
		.def(py::init<>());

    py::class_<MapType,        ElemType>(m, "MapType")
	    .def(py::init<>());

}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */