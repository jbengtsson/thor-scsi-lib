#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include "thor_scsi.h"
#include <thor_scsi/elements/field_kick.h>
#include <thor_scsi/elements/drift.h>
#include <thor_scsi/elements/marker.h>
#include <thor_scsi/elements/bpm.h>
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

class PyObserver: public tsc::Observer {
public:
	using tsc::Observer::Observer;

	void view(std::shared_ptr<const tsc::CellVoid> elem, const ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, view, elem, ps, os, cnt);
	}
	void view(std::shared_ptr<const tsc::CellVoid> elem, const ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, view, elem, ps, os, cnt);
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

class PyRadDelInt: public tse::RadiationDelegateInterface {
	using tse::RadiationDelegateInterface::RadiationDelegateInterface;
	void view(const tse::ElemType& elem, const ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateInterface, view, elem, ps, os, cnt);
	}
	void view(const tse::ElemType& elem, const ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateInterface, view, elem, ps, os, cnt);
	}
};

class PyRadDel: public tse::RadiationDelegate {
	using tse::RadiationDelegate::RadiationDelegate;
	void view(const tse::ElemType& elem, const ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegate, view, elem, ps, os, cnt);
	}
	void view(const tse::ElemType& elem, const ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegate, view, elem, ps, os, cnt);
	}
};


class PyRadDelKickInt: public tse::RadiationDelegateKickInterface {
	using tse::RadiationDelegateKickInterface::RadiationDelegateKickInterface;
	void view(const tse::FieldKickAPI& elem, const ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
	void view(const tse::FieldKickAPI& elem, const ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
};

class PyRadDelKick: public tse::RadiationDelegateKick {
	using tse::RadiationDelegateKick::RadiationDelegateKick;
	void view(const tse::FieldKickAPI& elem, const ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKick, view, elem, ps, os, cnt);
	}
	void view(const tse::FieldKickAPI& elem, const ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKick, view, elem, ps, os, cnt);
	}
};


static const char pass_d_doc[] = "pass the element (state as doubles)";
static const char pass_tps_doc[] = "pass the element (state as tps)";
/*
 * I tested if the derived classes have to mention that their memory needs to be
 * managed by shared pointers. Worked for me only if it is done by shared pointers
 */
void py_thor_scsi_init_elements(py::module &m)
{

	py::enum_<tsc::ObservedState>(m, "ObservedState")
		.value("event",     tsc::ObservedState::event)
		.value("start",     tsc::ObservedState::start)
		.value("end",       tsc::ObservedState::end)
		.value("undefined", tsc::ObservedState::undefined)
		.value("failure",   tsc::ObservedState::failure);


	py::class_<tsc::Observer,  PyObserver, std::shared_ptr<tsc::Observer>> observer(m, "Observer");
	observer.def(py::init<>());

	py::class_<tsc::ElemType,  PyElemType, std::shared_ptr<tsc::ElemType>> elem_type(m, "ElemType");
	elem_type.def("__str__",       &tsc::ElemType::pstr)
		.def("__repr__",       &tsc::ElemType::repr)
		.def_readonly("name",  &tsc::ElemType::name)
		.def_readonly("index", &tsc::ElemType::index)
		.def("getLength",      &tse::ElemType::getLength)
		.def("setLength",      &tse::ElemType::setLength)
		.def("getObserver",    &tse::ElemType::observer)
		.def("setObserver",    &tse::ElemType::set_observer)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<double>&>(&tse::ElemType::pass), pass_d_doc)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<tps>&>(&tse::ElemType::pass),    pass_tps_doc)
		.def(py::init<const Config &>());


	py::class_<tse::DriftType, std::shared_ptr<tse::DriftType>>(m, "Drift", elem_type)
		.def(py::init<const Config &>());

	py::class_<tse::RadiationDelegateInterface, PyRadDelInt, std::shared_ptr<tse::RadiationDelegateInterface>> rad_del_int(m, "RadiationDelegateInterface");
	rad_del_int
		.def("__repr__", &tse::RadiationDelegateInterface::repr)
		.def(py::init<>());

	py::class_<tse::RadiationDelegate, PyRadDel, std::shared_ptr<tse::RadiationDelegate>>(m, "RadiationDelegate", rad_del_int)
		.def("reset",       &tse::RadiationDelegate::reset)
		.def("getCurlydHx", &tse::RadiationDelegate::getCurlydHx)
		.def(py::init<>());

	py::class_<tse::RadiationDelegateKickInterface, PyRadDelKickInt, std::shared_ptr<tse::RadiationDelegateKickInterface>> rad_del_kick_int(m, "RadiationDelegateKickInterface");
	rad_del_int
		//.def("__repr__", &tse::RadiationDelegateKickInterface::repr)
		.def(py::init<>());

	py::class_<tse::RadiationDelegateKick, PyRadDelKick, std::shared_ptr<tse::RadiationDelegateKick>>(m, "RadiationDelegateKick", rad_del_kick_int)
		.def("reset",                              &tse::RadiationDelegateKick::reset)
		.def("getCurlydHx",                        &tse::RadiationDelegateKick::getCurlydHx)
		.def("getSynchrotronIntegralsIncrements",  &tse::RadiationDelegateKick::getSynchrotronIntegralsIncrement)
		.def("getDiffusionCoefficientsIncrements", &tse::RadiationDelegateKick::getDiffusionCoefficientsIncrement)
		.def("computeDiffusion",                   &tse::RadiationDelegateKick::computeDiffusion)
		.def("isComputingDiffusion",               &tse::RadiationDelegateKick::isComputingDiffusion)
		.def("setEnergy",                          &tse::RadiationDelegateKick::setEnergy)
		.def("getEnergy",                          &tse::RadiationDelegateKick::getEnergy)
		.def("view",
		     py::overload_cast<const tse::FieldKickAPI&, const ss_vect<double> &, const enum tsc::ObservedState, const int>(&tse::RadiationDelegateKick::view))
		.def("view",
		     py::overload_cast<const tse::FieldKickAPI&, const ss_vect<double> &, const enum tsc::ObservedState, const int>(&tse::RadiationDelegateKick::view))
		.def("__repr__",                           &tse::RadiationDelegateKick::repr)
		.def(py::init<>());

	py::class_<tse::MarkerType, std::shared_ptr<tse::MarkerType>>(m, "Marker", elem_type)
		.def("getRadiationDelegate", &tse::MarkerType::getRadiationDelegate)
		.def("setRadiationDelegate", &tse::MarkerType::setRadiationDelegate)
		.def(py::init<const Config &>());

	py::class_<tse::BPMType, std::shared_ptr<tse::BPMType>>(m, "BPM", elem_type)
		.def(py::init<const Config &>());

	//, std::shared_ptr<tse::>
	py::class_<tse::CavityType, std::shared_ptr<tse::CavityType>>(m, "Cavity", elem_type)
		.def("setFrequency",      &tse::CavityType::setFrequency)
		.def("getFrequency",      &tse::CavityType::getFrequency)
		.def("setVoltage",        &tse::CavityType::setVoltage)
		.def("getVoltage",        &tse::CavityType::getVoltage)
		.def("setPhase",          &tse::CavityType::setPhase)
		.def("getPhase",          &tse::CavityType::getPhase)
		.def("setHarmonicNumber", &tse::CavityType::setHarmonicNumber)
		.def("getHarmonicNumber", &tse::CavityType::getHarmonicNumber)
		.def(py::init<const Config &>());

	/*
	 * Needs to be defined as shared ptr class as it is returned as shared pointer
	 * by classical magnet's method getMultipoles
	*/
	py::class_<tsc::PlanarMultipoles, std::shared_ptr<tsc::PlanarMultipoles>>(m, "TwoDimensionalMultipoles")
		.def("getMultipole",     &tsc::PlanarMultipoles::getMultipole)
		.def("setMultipole",     &tsc::PlanarMultipoles::setMultipole)
		.def("applyRollAngle",   &tsc::PlanarMultipoles::applyRollAngle)
		.def("__str__",          &tsc::PlanarMultipoles::pstr)
		.def("__repr__",         &tsc::PlanarMultipoles::repr)
		.def("getCoefficients",  &tsc::PlanarMultipoles::getCoeffsConst)
		.def("applyTranslation", py::overload_cast<const tsc::cdbl>(&tsc::PlanarMultipoles::applyTranslation))
		.def("applyTranslation", py::overload_cast<const double, const double>(&tsc::PlanarMultipoles::applyTranslation))
		.def(py::init<>());

	py::class_<tse::FieldKick, std::shared_ptr<tse::FieldKick>> field_kick(m, "FieldKick", elem_type);
	field_kick
		.def("isThick",                     &tse::FieldKick::isThick)
		.def("asThick",                     &tse::FieldKick::asThick)
		.def("getNumberOfIntegrationSteps", &tse::FieldKick::getNumberOfIntegrationSteps)
		.def("setNumberOfIntegrationSteps", &tse::FieldKick::setNumberOfIntegrationSteps)
		.def("getIntegrationMethod",        &tse::FieldKick::getIntegrationMethod)
		.def("getCurvature",                &tse::FieldKick::getCurvature)
		.def("setCurvature",                &tse::FieldKick::setCurvature)
		.def("assumingCurvedTrajectory",    &tse::FieldKick::assumingCurvedTrajectory)
		.def("getBendingAngle",             &tse::FieldKick::getBendingAngle)
		.def("setBendingAngle",             &tse::FieldKick::setBendingAngle)
		.def("setEntranceAngle",            &tse::FieldKick::setEntranceAngle)
		.def("getEntranceAngle",            &tse::FieldKick::getEntranceAngle)
		.def("setExitAngle",                &tse::FieldKick::setExitAngle)
		.def("getExitAngle",                &tse::FieldKick::getExitAngle)
		.def("getRadiationDelegate",        &tse::FieldKick::getRadiationDelegate)
		.def("setRadiationDelegate",        &tse::FieldKick::setRadiationDelegate)
		.def(py::init<const Config &>());

	py::class_<tse::MpoleType, std::shared_ptr<tse::MpoleType>> mpole_type(m, "Mpole", field_kick);
	mpole_type
		.def(py::init<const Config &>());
	py::class_<tse::ClassicalMagnet, PyClassicalMagnet, std::shared_ptr<tse::ClassicalMagnet>> cm(m, "ClassicalMagnet", mpole_type);
	cm
		.def("getMultipoles",            &tse::ClassicalMagnet::getMultipoles)
		.def("getMainMultipoleNumber",   &tse::ClassicalMagnet::getMainMultipoleNumber)
		.def("getMainMultipoleStrength", &tse::ClassicalMagnet::getMainMultipoleStrength)
		.def("setMainMultipoleStrength", py::overload_cast<const double>(&tse::ClassicalMagnet::setMainMultipoleStrength))
		.def("setMainMultipoleStrength", py::overload_cast<const tsc::cdbl>(&tse::ClassicalMagnet::setMainMultipoleStrength))
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<double>&>(&tse::ClassicalMagnet::pass), pass_d_doc)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<tps>&>   (&tse::ClassicalMagnet::pass),    pass_tps_doc)
		.def(py::init<const Config &>());

	py::class_<tse::QuadrupoleType, std::shared_ptr<tse::QuadrupoleType>>(m, "Quadrupole", cm)
		.def(py::init<const Config &>());
	py::class_<tse::SextupoleType, std::shared_ptr<tse::SextupoleType>> (m, "Sextupole", cm)
		.def(py::init<const Config &>());
	py::class_<tse::OctupoleType, std::shared_ptr<tse::OctupoleType>>(m, "Octupole", cm)
		.def(py::init<const Config &>());
	py::class_<tse::BendingType, std::shared_ptr<tse::BendingType>>(m, "Bending", cm)
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
