#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include "thor_scsi.h"
#include <thor_scsi/elements/field_kick.h>
#include <thor_scsi/elements/drift.h>
#include <thor_scsi/elements/marker.h>
#include <thor_scsi/elements/bending.h>
#include <thor_scsi/elements/quadrupole.h>
#include <thor_scsi/elements/sextupole.h>
#include <thor_scsi/elements/octupole.h>
#include <thor_scsi/elements/cavity.h>


namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace py = pybind11;


class PyElemType: public tsc::ElemType {
public:
	using tsc::ElemType::ElemType;

	void pass(thor_scsi::core::ConfigType & conf, ss_vect<double>& ps) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::ElemType, pass, conf, ps);
	}
	void pass(thor_scsi::core::ConfigType & conf, ss_vect<tps>& ps) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::ElemType, pass, conf, ps);
	}
	const char * type_name(void) const override {
		PYBIND11_OVERRIDE_PURE(
			const char *, /* Return type */
			tsc::ElemType,      /* Parent class */
			type_name,          /* Name of function in C++ (must match Python name) */
			/* Argument(s) */
			);
	}
};

class PyClassicalMagnet: public tse::ClassicalMagnet {
	using tse::ClassicalMagnet::ClassicalMagnet;
	int getMainMultipoleNumber(void) const override {
		PYBIND11_OVERRIDE_PURE(
			int, /* Return type */
			tse::ClassicalMagnet,      /* Parent class */
			getMainMultipoleNumber,          /* Name of function in C++ (must match Python name) */
			/* Argument(s) */
			);
	}
};


static const char pass_d_doc[] = "pass the element (state as doubles)";
static const char pass_tps_doc[] = "pass the element (state as tps)";
void py_thor_scsi_init_elements(py::module_ &m)
{

	py::class_<tsc::ElemType, PyElemType>(m, "ElemType")
		.def("__str__",       &tsc::CellVoid::pstr)
		.def("__repr__",      &tsc::CellVoid::repr)
		.def_readonly("name", &tsc::CellVoid::name)
		.def("getLength", &tse::DriftType::getLength)
		.def("setLength", &tse::DriftType::setLength)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<double>&>(&tse::ElemType::pass), pass_d_doc)
		 .def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<tps>&>(&tse::ElemType::pass),    pass_tps_doc)
		.def(py::init<const Config &>());

	py::class_<tse::DriftType, tsc::ElemType>(m, "Drift")
		.def(py::init<const Config &>());

	py::class_<tse::MarkerType, tsc::ElemType>(m, "Marker")
		.def(py::init<const Config &>());

	py::class_<tse::CavityType, tsc::ElemType>(m, "Cavity")
		.def("setFrequency",      &tse::CavityType::setFrequency)
		.def("getFrequency",      &tse::CavityType::getFrequency)
		.def("setVoltage",        &tse::CavityType::setVoltage)
		.def("getVoltage",        &tse::CavityType::getVoltage)
		.def("setPhase",          &tse::CavityType::setPhase)
		.def("getPhase",          &tse::CavityType::getPhase)
		.def("setHarmonicNumber", &tse::CavityType::setHarmonicNumber)
		.def("getHarmonicNumber", &tse::CavityType::getHarmonicNumber)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<double>&>(&tse::CavityType::pass), pass_d_doc)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<tps>&>   (&tse::CavityType::pass),    pass_tps_doc)
		.def(py::init<const Config &>());

	/*
	 * Needs to be defined as shared ptr class as it is returned as shared pointer
	 * by classical magnet's method getMultipoles
	*/
	py::class_<tsc::PlanarMultipoles, std::shared_ptr<tsc::PlanarMultipoles>>(m, "TwoDimensionalMultipoles")
		.def("getMultipole", &tsc::PlanarMultipoles::getMultipole)
		.def("setMultipole", &tsc::PlanarMultipoles::setMultipole)
		.def("__str__",   &tsc::PlanarMultipoles::pstr)
		.def("__repr__",  &tsc::PlanarMultipoles::repr)
		.def(py::init<>());

	py::class_<tse::ClassicalMagnet, PyClassicalMagnet, tsc::ElemType>(m, "ClassicalMagnet")
		//.def("__str__",   &tse::ClassicalMagnet::pstr)
		//.def("__repr__",  &tse::ClassicalMagnet::repr)
		.def("getMultipoles",            &tse::ClassicalMagnet::getMultipoles)
		.def("getMainMultipoleNumber",   &tse::ClassicalMagnet::getMainMultipoleNumber)
		.def("getMainMultipoleStrength", &tse::ClassicalMagnet::getMainMultipoleStrength)
		.def("setMainMultipoleStrength", py::overload_cast<const double>(&tse::ClassicalMagnet::setMainMultipoleStrength))
		.def("setMainMultipoleStrength", py::overload_cast<const tsc::cdbl>(&tse::ClassicalMagnet::setMainMultipoleStrength))
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<double>&>(&tse::ClassicalMagnet::pass), pass_d_doc)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<tps>&>   (&tse::ClassicalMagnet::pass),    pass_tps_doc)
		.def(py::init<const Config &>());

	py::class_<tse::QuadrupoleType, tse::ClassicalMagnet>(m, "Quadrupole")
		.def(py::init<const Config &>());
	py::class_<tse::SextupoleType, tse::ClassicalMagnet>(m, "Sextupole")
		.def(py::init<const Config &>());
	py::class_<tse::OctupoleType, tse::ClassicalMagnet>(m, "Octupole")
		.def(py::init<const Config &>());


#if 0
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
#endif

}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
