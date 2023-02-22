#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <thor_scsi/elements/radiation_delegate.h>

namespace py = pybind11;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

class PyRadDelInt: public tse::RadiationDelegateInterface {
	using tse::RadiationDelegateInterface::RadiationDelegateInterface;
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateInterface, view, elem, ps, os, cnt);
	}
    /*
        void view(const tse::ElemType& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateInterface, view, elem, ps, os, cnt);
	}
    */
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateInterface, view, elem, ps, os, cnt);
	}
};

class PyRadDel: public tse::RadiationDelegate {
	using tse::RadiationDelegate::RadiationDelegate;
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegate, view, elem, ps, os, cnt);
	}
	/*
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegate, view, elem, ps, os, cnt);
	}
	*/
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegate, view, elem, ps, os, cnt);
	}
};


class PyRadDelKickInt: public tse::RadiationDelegateKickInterface {
	using tse::RadiationDelegateKickInterface::RadiationDelegateKickInterface;
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
	/*
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
	*/
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
};

class PyRadDelKick: public tse::RadiationDelegateKick {
	using tse::RadiationDelegateKick::RadiationDelegateKick;
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKick, view, elem, ps, os, cnt);
	}
	/*
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKick, view, elem, ps, os, cnt);
	}
	*/
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKick, view, elem, ps, os, cnt);
	}
};


void py_thor_scsi_init_radiation(py::module &m)
{

	py::class_<tse::RadiationDelegateInterface, PyRadDelInt, std::shared_ptr<tse::RadiationDelegateInterface>> rad_del_int(m, "RadiationDelegateInterface");
	rad_del_int
		.def("__repr__",        &tse::RadiationDelegateInterface::repr)
		.def(py::init<>());

	py::class_<tse::RadiationDelegate, PyRadDel, std::shared_ptr<tse::RadiationDelegate>>(m, "RadiationDelegate", rad_del_int)
		.def("reset",               &tse::RadiationDelegate::reset)
		.def("get_curly_dHx",       &tse::RadiationDelegate::getCurlydHx)
		.def("get_delegator_name",  &tse::RadiationDelegate::getDelegatorName)
		.def("get_delegator_index", &tse::RadiationDelegate::getDelegatorIndex)
		.def("view",
		     py::overload_cast<const tse::ElemType&, const gtpsa::ss_vect<double> &, const enum tsc::ObservedState, const int>(&tse::RadiationDelegate::view))
		.def("view",
		     py::overload_cast<const tse::ElemType&, const gtpsa::ss_vect<double> &, const enum tsc::ObservedState, const int>(&tse::RadiationDelegate::view))
		.def(py::init<>());

	py::class_<tse::RadiationDelegateKickInterface, PyRadDelKickInt, std::shared_ptr<tse::RadiationDelegateKickInterface>> rad_del_kick_int(m, "RadiationDelegateKickInterface");
	rad_del_int
		//.def("__repr__", &tse::RadiationDelegateKickInterface::repr)
		.def(py::init<>());

	// Why did it not end up at Johan?
	py::class_<tse::RadiationDelegateKick, PyRadDelKick, std::shared_ptr<tse::RadiationDelegateKick>>(m, "RadiationDelegateKick", rad_del_kick_int)
		.def("reset",                                 &tse::RadiationDelegateKick::reset)
		.def("get_curly_dHx",                         &tse::RadiationDelegateKick::getCurlydHx)
		.def("get_delegator_name",                    &tse::RadiationDelegateKick::getDelegatorName)
		.def("get_delegator_index",                   &tse::RadiationDelegateKick::getDelegatorIndex)
		.def("get_synchrotron_integrals_increments",  &tse::RadiationDelegateKick::getSynchrotronIntegralsIncrement)
		.def("get_diffusion_coefficients_increments", &tse::RadiationDelegateKick::getDiffusionCoefficientsIncrement)
		.def("compute_diffusion",                     &tse::RadiationDelegateKick::computeDiffusion)
		.def("is_computing_diffusion",                &tse::RadiationDelegateKick::isComputingDiffusion)
		.def("set_energy",                            &tse::RadiationDelegateKick::setEnergy)
		.def("get_energy",                            &tse::RadiationDelegateKick::getEnergy)
		.def("view",
		     py::overload_cast<const tse::FieldKickAPI&, const gtpsa::ss_vect<double> &, const enum tsc::ObservedState, const int>(&tse::RadiationDelegateKick::view))
		.def("view",
		     py::overload_cast<const tse::FieldKickAPI&, const gtpsa::ss_vect<double> &, const enum tsc::ObservedState, const int>(&tse::RadiationDelegateKick::view))
		.def("__repr__",                           &tse::RadiationDelegateKick::repr)
		.def(py::init<>());

}
