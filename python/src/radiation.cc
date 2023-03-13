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


template <class FC>
class PyRadDelKickInt: public tse::RadiationDelegateKickInterfaceKnobbed<FC> {
	using tse::RadiationDelegateKickInterfaceKnobbed<FC>::RadiationDelegateKickInterfaceKnobbed;
	void view(const FC& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
	/*
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
	*/
	void view(const FC& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
};


template <class FC>
class PyRadDelKick: public tse::RadiationDelegateKickKnobbed<FC> {
	using tse::RadiationDelegateKickKnobbed<FC>::RadiationDelegateKickKnobbed;
	void view(const FC& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKickKnobbed<FC>, view, elem, ps, os, cnt);
	}
	/*
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKick, view, elem, ps, os, cnt);
	}
	*/
	void view(const FC& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKickKnobbed<FC>, view, elem, ps, os, cnt);
	}
};


template<typename FC, typename Class>
void add_methods_radiation_field_kick(py::class_<Class> t_kick)
{


    //typedef tse::RadiationDelegateKick<FC> rk_t;


    t_kick
		.def("reset",                                 &Class::reset)
		.def("get_curly_dHx",                         &Class::getCurlydHx)
		.def("get_delegator_name",                    &Class::getDelegatorName)
		.def("get_delegator_index",                   &Class::getDelegatorIndex)
		.def("get_synchrotron_integrals_increments",  &Class::getSynchrotronIntegralsIncrement)
		.def("get_diffusion_coefficients_increments", &Class::getDiffusionCoefficientsIncrement)
		.def("compute_diffusion",                     &Class::computeDiffusion)
		.def("is_computing_diffusion",                &Class::isComputingDiffusion)
		.def("set_energy",                            &Class::setEnergy)
		.def("get_energy",                            &Class::getEnergy)
		.def("view",
		     py::overload_cast<const FC&, const gtpsa::ss_vect<double> &, const enum tsc::ObservedState, const int>(&Class::view))
		.def("view",
		     py::overload_cast<const FC&, const gtpsa::ss_vect<double> &, const enum tsc::ObservedState, const int>(&Class::view))
		.def("__repr__",                           &Class::repr)
	;
}

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

	typedef tse::FieldKickAPIKnobbed<tsc::StandardDoubleType> fka_dt;
	typedef tse::FieldKickAPIKnobbed<tsc::TpsaVariantType> fka_dvt;

	typedef tse::RadiationDelegateKickInterfaceKnobbed<fka_dt>  rdlki_d_t;
	typedef tse::RadiationDelegateKickInterfaceKnobbed<fka_dvt> rdlki_t_t;
	py::class_<rdlki_d_t, PyRadDelKickInt<fka_dt>, std::shared_ptr<rdlki_d_t>> rad_del_kick_int(m, "RadiationDelegateKickInterface");
	rad_del_int
		//.def("__repr__", &tse::RadiationDelegateKickInterface::repr)
		.def(py::init<>());
	py::class_<rdlki_t_t, PyRadDelKickInt<fka_dvt>, std::shared_ptr<rdlki_t_t>> rad_del_kick_int_tpsa(m, "RadiationDelegateKickInterfaceTpsa");
	rad_del_int
		//.def("__repr__", &tse::RadiationDelegateKickInterface::repr)
		.def(py::init<>());

	// Why did it not end up at Johan?
	typedef tse::RadiationDelegateKickKnobbed<fka_dt>  rdlk_d_t;
	typedef tse::RadiationDelegateKickKnobbed<fka_dvt> rdlk_t_t;

	py::class_<rdlk_d_t, PyRadDelKick<fka_dt>, std::shared_ptr<rdlk_d_t>> rad_del_kick(m, "RadiationDelegateKick", rad_del_kick_int);
	rad_del_kick
		.def(py::init<>());
	add_methods_radiation_field_kick<fka_dt, rdlk_d_t>(rad_del_kick);

	py::class_<rdlk_t_t, PyRadDelKick<fka_dvt>, std::shared_ptr<rdlk_t_t>> rad_del_kick_tpsa(m, "RadiationDelegateKickTpsa", rad_del_kick_int);
	rad_del_kick_tpsa
		.def(py::init<>());
	add_methods_radiation_field_kick<fka_dvt, rdlk_t_t>(rad_del_kick_tpsa);

}
