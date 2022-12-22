#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "thor_scsi.h"
#include <thor_scsi/elements/field_kick.h>
#include <thor_scsi/elements/drift.h>
#include <thor_scsi/elements/marker.h>
#include <thor_scsi/elements/bpm.h>
#include <thor_scsi/elements/bending.h>
#include <thor_scsi/elements/quadrupole.h>
#include <thor_scsi/elements/sextupole.h>
#include <thor_scsi/elements/octupole.h>
#include <thor_scsi/elements/corrector.h>
#include <thor_scsi/elements/cavity.h>
#include <thor_scsi/elements/standard_observer.h>
#include <thor_scsi/elements/standard_aperture.h>


namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace py = pybind11;


class PyElemType: public tsc::ElemType {
public:
	using tsc::ElemType::ElemType;

	void propagate(thor_scsi::core::ConfigType & conf, gtpsa::ss_vect<double>& ps) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::ElemType, propagate, conf, ps);
	}
	void propagate(thor_scsi::core::ConfigType & conf, gtpsa::ss_vect<tps>& ps) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::ElemType, propagate, conf, ps);
	}
	void propagate(thor_scsi::core::ConfigType & conf, gtpsa::ss_vect<gtpsa::tpsa>& ps) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::ElemType, propagate, conf, ps);
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

	void view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, view, elem, ps, os, cnt);
	}
	void view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, view, elem, ps, os, cnt);
	}
	void view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, view, elem, ps, os, cnt);
	}

	void show(std::ostream& strm, int level) const override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, strm, level);
	}
};

class Py2DAperture: public tsc::TwoDimensionalAperture {
public:
	using tsc::TwoDimensionalAperture::TwoDimensionalAperture;

	double isWithin(const double x, const double y) const override {
		PYBIND11_OVERRIDE_PURE(double, tsc::TwoDimensionalAperture, x, y);
	};
	const char * type_name(void) const override {
		PYBIND11_OVERRIDE_PURE(const char *, tsc::TwoDimensionalAperture, void);
	};

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
	bool isSkew(void) const override {
		PYBIND11_OVERRIDE_PURE(
			bool, /* Return type */
			tse::ClassicalMagnet,      /* Parent class */
			isSkew,          /* Name of function in C++ (must match Python name) */
			/* Argument(s) */
			);
	}
};

class PyRadDelInt: public tse::RadiationDelegateInterface {
	using tse::RadiationDelegateInterface::RadiationDelegateInterface;
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateInterface, view, elem, ps, os, cnt);
	}
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateInterface, view, elem, ps, os, cnt);
	}
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateInterface, view, elem, ps, os, cnt);
	}
};

class PyRadDel: public tse::RadiationDelegate {
	using tse::RadiationDelegate::RadiationDelegate;
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegate, view, elem, ps, os, cnt);
	}
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegate, view, elem, ps, os, cnt);
	}
	void view(const tse::ElemType& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegate, view, elem, ps, os, cnt);
	}
};


class PyRadDelKickInt: public tse::RadiationDelegateKickInterface {
	using tse::RadiationDelegateKickInterface::RadiationDelegateKickInterface;
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tse::RadiationDelegateKickInterface, view, elem, ps, os, cnt);
	}
};

class PyRadDelKick: public tse::RadiationDelegateKick {
	using tse::RadiationDelegateKick::RadiationDelegateKick;
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKick, view, elem, ps, os, cnt);
	}
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKick, view, elem, ps, os, cnt);
	}
	void view(const tse::FieldKickAPI& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE(void, tse::RadiationDelegateKick, view, elem, ps, os, cnt);
	}
};

class PyField2DInterpolation: public tsc::Field2DInterpolation{
public:
	using tsc::Field2DInterpolation::Field2DInterpolation;

#if 0
	void field_py(const std::array<double, 2> pos, std::array<double, 2> field) const {
		PYBIND11_OVERRIDE_PURE(void, PyField2DInterpolationIntermediate, field_py, pos, field);
	}
#endif
	// clang now accepting this code
	// needs to be tested if it still works as expected (data can be returned)
	// originally made for non linear kicker
	// this way it does not work ....
	virtual void field_py    (const py::array_t<double>  &t_pos, py::array_t<double> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, field_py, t_pos, t_field);
	}
	virtual void field_py    (const std::array<tps    , 2>  &t_pos, std::array<tps   , 2> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, field_py   , t_pos, t_field);
	}
#if 0
	virtual void gradient_py (const std::array<double&, 2> &t_pos, std::array<double& , 2> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, gradient_py, t_pos, t_field);
	}
	virtual void gradient_py (const std::array<tps&   , 2>  &t_pos, std::array<tps&   , 2> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, gradient_py, t_pos, t_field);
	}
	virtual void gradient_py (const std::array<tps&   , 2>  &t_pos, std::array<double&, 2> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, gradient_py, t_pos, t_field);
	}
#endif
private:

	template<typename T>
	inline void _gradient(const T x, const T y, T *Bx, T *By) const {
		throw std::runtime_error("Not implemented");
	}

public:
	void field(const double x, const double y, double *Bx, double *By) const override {
		auto pos = py::array_t<double>{2};
		auto t_field = py::array_t<double>{2};
		double *pos_p = static_cast<double *>(pos.request().ptr),
			*field_p = static_cast<double *>(t_field.request().ptr);
		pos_p[0] = x;
		pos_p[1] = y;
		this->field_py(pos, t_field);
		*Bx = field_p[0];
		*By = field_p[1];

		/*
		std::cerr << "x " << x << ", y "
			  << ", By " << *By <<  ",  Bx " <<  *Bx
			  << std::endl;
		*/
	}

	void field(const tps x, const tps y, tps *Bx, tps *By) const override {
		//this->_field(x, y, Bx, By);
		std::array<tps, 2> pos = {x, y};
		std::array<tps, 2> t_field;
		this->field_py(pos, t_field);
		*Bx = t_field[0];
		*By = t_field[1];

	}
	void field(const gtpsa::tpsa x, const gtpsa::tpsa y, gtpsa::tpsa *Bx, gtpsa::tpsa *By) const override {
		//this->_field(x, y, Bx, By);
		std::array<gtpsa::tpsa, 2> pos = {x, y};
		std::array<gtpsa::tpsa*, 2> t_field = {Bx, By};
		// this->field_py(pos, t_field);
		//*Bx = t_field[0];
		//*By = t_field[1];

	}
	void gradient(const double x, const double y, double *Gx, double *Gy) const  override{
		PYBIND11_OVERRIDE_PURE(void, tsc::Field2DInterpolation, gradient, x, y, Gx, Gy);
	}
	void gradient(const tps x, const tps y, tps *Gx, tps *Gy) const   override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Field2DInterpolation, gradient, x, y, Gx, Gy);
	}
	void gradient(const tps x, const tps y, double *Gx, double *Gy) const  override{
		PYBIND11_OVERRIDE_PURE(void, tsc::Field2DInterpolation, gradient, x, y, Gx, Gy);
	}
	void gradient(const gtpsa::tpsa x, const gtpsa::tpsa y, gtpsa::tpsa *Gx, gtpsa::tpsa *Gy) const   override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Field2DInterpolation, gradient, x, y, Gx, Gy);
	}
	void gradient(const gtpsa::tpsa x, const gtpsa::tpsa y, double *Gx, double *Gy) const  override{
		PYBIND11_OVERRIDE_PURE(void, tsc::Field2DInterpolation, gradient, x, y, Gx, Gy);
	}

	void show(std::ostream& strm, int level) const override {
		const std::string txt = (level <= 1) ? this->pstr() : this->repr();
		strm << txt;
	}

	std::string repr(void) const {
		PYBIND11_OVERRIDE_NAME(std::string, tsc::Field2DInterpolation, "__repr__", repr);
	}

	std::string pstr(void) const {
		PYBIND11_OVERRIDE_NAME(std::string, tsc::Field2DInterpolation, "__str__", pstr);
	}

};


class PyField2DInterpolationProxy: public tsc::Field2DInterpolation
{
};
#if 0
class PyFieldKick: public tse::FieldKick {
	using tse::FieldKick::FieldKick;

	void show(std::ostream& strm, int level) const override {
		const std::string txt = (level <= 1) ? this->pstr() : this->repr();
		strm << txt;
	}

	std::string repr(void) const {
		PYBIND11_OVERRIDE_NAME(std::string, tsc::FieldKick, "__repr__", repr);
	}

	std::string pstr(void) const {
		PYBIND11_OVERRIDE_NAME(std::string, tsc::FieldKick, "__str__", pstr);
	}
};
#endif

static const char pass_d_doc[] = "pass the element (state as doubles)";
static const char pass_tps_doc[] = "pass the element (state as tps)";
static const char pass_tpsa_doc[] = "pass the element (state as gtpsa)";
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

	py::class_<tse::StandardObserver, tsc::Observer, std::shared_ptr<tse::StandardObserver>> std_observer(m, "StandardObserver");
	std_observer
		.def("__str__",                      &tse::StandardObserver::pstr)
		.def("__repr__",                     &tse::StandardObserver::repr)
		.def("get_observed_name",            &tse::StandardObserver::getObservedName)
		.def("get_observed_index",           &tse::StandardObserver::getObservedIndex)
		.def("has_phase_space",              &tse::StandardObserver::hasPhaseSpace)
		.def("get_phase_space",              &tse::StandardObserver::getPhaseSpace)
		.def("has_truncated_power_series",   &tse::StandardObserver::hasTruncatedPowerSeries)
		.def("get_truncated_power_series",   &tse::StandardObserver::getTruncatedPowerSeries)
		.def("has_truncated_power_series_a", &tse::StandardObserver::hasTruncatedPowerSeriesA)
		.def("get_truncated_power_series_a", &tse::StandardObserver::getTruncatedPowerSeriesA)
		.def("reset",                        &tse::StandardObserver::reset)
		.def(py::init<>());

	py::class_<tsc::TwoDimensionalAperture,  Py2DAperture, std::shared_ptr<tsc::TwoDimensionalAperture>> aperture(m, "Aperture");
	aperture.def("__repr__", &tsc::TwoDimensionalAperture::repr)
		.def(py::init<>());

	const char rect_ap_doc[] = "initialise rectangular aperture";
	py::class_<tse::RectangularAperture, tsc::TwoDimensionalAperture, std::shared_ptr<tse::RectangularAperture>> rect_ap(m, "RectangularAperture");
	rect_ap.def(py::init<const double, const double, const double, const double>(), rect_ap_doc,
		    py::arg("width"), py::arg("height"), py::arg("x") = 0, py::arg("y") = 0);

	const char circ_ap_doc[] = "initialise round aperture";
	py::class_<tse::CircularAperture, tsc::TwoDimensionalAperture, std::shared_ptr<tse::CircularAperture>> circ_ap(m, "CircularAperture");
	circ_ap.def(py::init<const double, const double, const double>(), circ_ap_doc,
		    py::arg("radius"), py::arg("x") = 0, py::arg("y") = 0);


	py::class_<tsc::CellVoid, std::shared_ptr<tsc::CellVoid>> cell_void(m, "CellVoid");
	py::class_<tsc::ElemType,  PyElemType, std::shared_ptr<tsc::ElemType>> elem_type(m, "ElemType", cell_void);
	elem_type.def("__str__",       &tsc::ElemType::pstr)
		.def("__repr__",       &tsc::ElemType::repr)
		.def_readonly("name",  &tsc::ElemType::name)
		.def_readonly("index", &tsc::ElemType::index)
		.def("get_length",     &tse::ElemType::getLength)
		.def("set_length",     &tse::ElemType::setLength)
		.def("get_observer",   &tse::ElemType::observer)
		.def("set_observer",   &tse::ElemType::set_observer)
		.def("get_aperture",   &tse::ElemType::getAperture)
		.def("set_aperture",   &tse::ElemType::setAperture)
		.def("propagate", py::overload_cast<tsc::ConfigType&, gtpsa::ss_vect<double>&>(&tse::ElemType::propagate), pass_d_doc)
		.def("propagate", py::overload_cast<tsc::ConfigType&, gtpsa::ss_vect<tps>&>(&tse::ElemType::propagate),    pass_tpsa_doc)
		.def(py::init<const Config &>());


	py::class_<tse::DriftType, std::shared_ptr<tse::DriftType>>(m, "Drift", elem_type)
		.def(py::init<const Config &>());

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

	py::class_<tse::MarkerType, std::shared_ptr<tse::MarkerType>>(m, "Marker", elem_type)
		.def("get_radiation_delegate", &tse::MarkerType::getRadiationDelegate)
		.def("set_radiation_delegate", &tse::MarkerType::setRadiationDelegate)
		.def(py::init<const Config &>());

	py::class_<tse::BPMType, std::shared_ptr<tse::BPMType>>(m, "BPM", elem_type)
		.def(py::init<const Config &>());

	//, std::shared_ptr<tse::>
	py::class_<tse::CavityType, std::shared_ptr<tse::CavityType>>(m, "Cavity", elem_type)
		.def("set_frequency",       &tse::CavityType::setFrequency)
		.def("get_frequency",       &tse::CavityType::getFrequency)
		.def("set_voltage",         &tse::CavityType::setVoltage)
		.def("get_voltage",         &tse::CavityType::getVoltage)
		.def("set_phase",           &tse::CavityType::setPhase)
		.def("get_phase",           &tse::CavityType::getPhase)
		.def("set_harmonic_number", &tse::CavityType::setHarmonicNumber)
		.def("get_harmonic_number", &tse::CavityType::getHarmonicNumber)
		.def(py::init<const Config &>());

	/*
	 * Needs to be defined as shared ptr class as it is returned as shared pointer
	 * by classical magnet's method getMultipoles
	*/


	// py::class_<tsc::Observer,  PyObserver, std::shared_ptr<tsc::Observer>> observer(m, "Observer");
	// observer.def(py::init<>());
	py::class_<tsc::Field2DInterpolation, PyField2DInterpolation,  std::shared_ptr<tsc::Field2DInterpolation>> field2dintp(m, "Field2DInterpolation");
	field2dintp
		// can I write the following lines in a c++-14 compatible way
		//.def("field",        py::overload_cast<const double, const double, double *, double *>(&tsc::Field2DInterpolation::field))
		.def("field",  [](tsc::Field2DInterpolation &intp, const py::array_t<double> t_pos) -> py::array_t<double> {
				py::buffer_info t_pos_buf = t_pos.request();
				if (t_pos_buf.ndim != 1){
					throw std::runtime_error("position must be a vector");
				}
				if (t_pos_buf.size != 2){
					throw std::runtime_error("position must be a vector of size 2");
				}
				auto t_field = py::array_t<double>(t_pos_buf.size);
				py::buffer_info t_field_buf = t_field.request();
				double *pos_p =  static_cast<double *>(t_pos_buf.ptr);
				double *field_p =  static_cast<double *>(t_field_buf.ptr);

				intp.field(pos_p[0], pos_p[1], &field_p[0], &field_p[1]);
				// std::cerr << "Returning field " << field_p[0] << ", " << field_p[1] << "\n";
				return t_field;
			})
		.def("field",  [](tsc::Field2DInterpolation &intp, const std::array<tps, 2> t_pos, std::array<tps, 2> t_field) -> void {
				       tps Bx, By;
				       intp.field(t_pos[0], t_pos[1], &Bx, &By);
				       t_field[0] = Bx;
				       t_field[1] = By;
			       }
		)
		.def(py::init<>());

	/*
	py::class_<PyField2DInterpolation, std::shared_ptr<PyField2DInterpolation>> field2dintp(m, "Field2DInterpolation", _field2dintp);
	field2dintp
		.def("field", static_cast<void (PyField2DInterpolation::*)(const std::array<double, 2>&,  std::array<double, 2>&) const>(&PyField2DInterpolation::field_py))
		.def("field", static_cast<void (PyField2DInterpolation::*)(const std::array<tps, 2>&, std::array<tps, 2>&) const>(&PyField2DInterpolation::field_py))
		.def(py::init<>());
	*/
	py::class_<tsc::TwoDimensionalMultipoles, std::shared_ptr<tsc::TwoDimensionalMultipoles>> multipole(m, "TwoDimensionalMultipoles", field2dintp
													    , py::buffer_protocol()
		);
	multipole
		.def("get_multipole",     &tsc::TwoDimensionalMultipoles::getMultipole)
		.def("set_multipole",     &tsc::TwoDimensionalMultipoles::setMultipole)
		.def("apply_roll_angle",  &tsc::TwoDimensionalMultipoles::applyRollAngle)
		.def("__str__",           &tsc::TwoDimensionalMultipoles::pstr)
		.def("__repr__",          &tsc::TwoDimensionalMultipoles::repr)
		.def("__len__",           &tsc::TwoDimensionalMultipoles::size)
		.def("get_coefficients",  &tsc::TwoDimensionalMultipoles::getCoeffsConst)
		//.def("field", static_cast<void (tsc::Field2DInterpolation::*)(const double, const double, double *, double *) const>(&tsc::TwoDimensionalMultipoles::field))
		//.def("field", static_cast<void (tsc::Field2DInterpolation::*)(const tps, const tps, tps *, tps *) const>(&tsc::TwoDimensionalMultipoles::field))
		.def("apply_translation", py::overload_cast<const tsc::cdbl>(&tsc::TwoDimensionalMultipoles::applyTranslation))
		.def("apply_translation", py::overload_cast<const double, const double>(&tsc::TwoDimensionalMultipoles::applyTranslation))
		.def(py::self += py::self)
		.def(py::self + py::self)
		.def(py::self += double())
		.def(py::self + double())
		.def(py::self *= double())
		.def(py::self * double())
		// .def(double() * py::self)
		.def(py::self *= std::vector<double>())
		.def(py::self * std::vector<double>())
		.def_buffer([](tsc::TwoDimensionalMultipoles &muls) -> py::buffer_info {
				    std::vector<tsc::cdbl_intern> coeffs = muls.getCoeffs();
				    size_t n = coeffs.size();
				    py::buffer_info r;
				    r.ptr = coeffs.data();         /* Pointer to buffer */
				    r.itemsize = sizeof(tsc::cdbl_intern); /* Size of one scalar */
				    r.format = py::format_descriptor<tsc::cdbl_intern>::format(); /* Python struct-style format descriptor */
				    r.ndim = 1;
				    r.shape = { static_cast<py::ssize_t>(n)};/* Number of dimensions */
				    r.strides = {
					    /* Strides (in bytes) for each index */
					    static_cast<py::ssize_t>(sizeof(tsc::cdbl_intern))
				    };
				    return r;
			    })
		.def(py::init<std::vector<tsc::cdbl_intern>>())
#if 0
		.def(py::init([](py::array_t<tsc::cdbl_intern, py::array::c_style|py::array::forcecast> b) {
				      /* Request a buffer descriptor from Python */
				      py::buffer_info info = b.request();

				      /* Some sanity checks ... */
				      if (info.format != py::format_descriptor<tsc::cdbl_intern>::format())
					      throw std::runtime_error("Incompatible format: expected a double array!");

				      if (info.ndim != 1){
					      std::stringstream strm;
					      strm << "Incompatible buffer: expected 1 but got "
						   << info.ndim << "dimensions!";
					      throw std::runtime_error(strm.str());
				      }
				      /*
				      std::vector<tsc mat = arma::mat(static_cast<const double *>(info.ptr),
							   info.shape[0], info.shape[1]);
				      */
				      return mat;
			      }))
#endif
		.def(py::init<const unsigned int>(), "initalise multipoles", py::arg("h_max") = tsc::max_multipole);

#if 0
	py::class_<tsc::BegninTwoDimensionalMultipoles, std::shared_ptr<tsc::BegninTwoDimensionalMultipoles>>(m, "BegninTwoDimensionalMultipoles", multipole)
		.def(py::self += py::self)
		.def(py::self + py::self)
#if 0
		.def(py::self += double())
		.def(py::self + double())
		.def(py::self *= double())
		.def(py::self * double())
#endif
		.def(py::self *= std::vector<double>())
		.def(py::self * std::vector<double>());
#endif


	/*
	 * causes seg fault .... memory management ...
	 *
	py::class_<tsc::PhaseSpaceGalileanPRot2DTransform, std::shared_ptr<tsc::PhaseSpaceGalileanPRot2DTransform>> prtf(m, "PhaseSpaceGalileanPRot2DTransform");
	prtf
		.def("setDx", &tsc::PhaseSpaceGalileanPRot2DTransform::setDx);
	*/
	py::class_<tse::FieldKick, std::shared_ptr<tse::FieldKick>> field_kick(m, "FieldKick", elem_type);
	field_kick
		// .def("getTransform",                &tse::FieldKick::getTransform) // causes segfault
		.def("set_dx",                          [](tse::FieldKick &kick, const double dx){kick.getTransform()->setDx(dx);})
		.def("set_dy",                          [](tse::FieldKick &kick, const double dy){kick.getTransform()->setDy(dy);})
		.def("set_roll",                        [](tse::FieldKick &kick, const double roll){kick.getTransform()->setRoll(roll);})
		.def("get_dx",                          [](tse::FieldKick &kick) -> double {return kick.getTransform()->getDx();})
		.def("get_dy",                          [](tse::FieldKick &kick) -> double {return kick.getTransform()->getDy();})
		.def("get_roll",                        [](tse::FieldKick &kick) -> double {return kick.getTransform()->getRoll();})
		.def("is_thick",                        &tse::FieldKick::isThick)
		.def("as_thick",                        &tse::FieldKick::asThick)
		.def("get_number_of_integration_steps", &tse::FieldKick::getNumberOfIntegrationSteps)
		.def("set_number_of_integration_steps", &tse::FieldKick::setNumberOfIntegrationSteps)
		.def("get_integration_method",          &tse::FieldKick::getIntegrationMethod)
		.def("get_curvature",                   &tse::FieldKick::getCurvature)
		.def("set_curvature",                   &tse::FieldKick::setCurvature)
		.def("assuming_curved_trajectory",      &tse::FieldKick::assumingCurvedTrajectory)
		.def("get_bending_angle",               &tse::FieldKick::getBendingAngle)
		.def("set_bending_angle",               &tse::FieldKick::setBendingAngle)
		.def("set_entrance_angle",              &tse::FieldKick::setEntranceAngle)
		.def("get_entrance_angle",              &tse::FieldKick::getEntranceAngle)
		.def("set_exit_angle",                  &tse::FieldKick::setExitAngle)
		.def("get_exit_angle",                  &tse::FieldKick::getExitAngle)
		.def("get_radiation_delegate",          &tse::FieldKick::getRadiationDelegate)
		.def("set_radiation_delegate",          &tse::FieldKick::setRadiationDelegate)
		.def("get_field_interpolator",          &tse::FieldKick::getFieldInterpolator)
		.def("set_field_interpolator",          &tse::FieldKick::setFieldInterpolator)
		.def(py::init<const Config &>());

	py::class_<tse::MpoleType, std::shared_ptr<tse::MpoleType>> mpole_type(m, "Mpole", field_kick);
	mpole_type
		.def("getFieldInterpolator",            &tse::ClassicalMagnet::getFieldInterpolator)
		.def(py::init<const Config &>());

	py::class_<tse::ClassicalMagnet, PyClassicalMagnet, std::shared_ptr<tse::ClassicalMagnet>> cm(m, "ClassicalMagnet", mpole_type);
	cm
		.def("get_multipoles",            &tse::ClassicalMagnet::getMultipoles)
		//.def("getBegninMultipoles",      &tse::ClassicalMagnet::getBegninMultipoles)
		// .def("setMultipoles",            &tse::ClassicalMagnet::setMultipoles)
		.def("get_main_multipole_number",   &tse::ClassicalMagnet::getMainMultipoleNumber)
		.def("get_main_multipole_strength", &tse::ClassicalMagnet::getMainMultipoleStrength)
		.def("get_main_multipole_strength_component", &tse::ClassicalMagnet::getMainMultipoleStrengthComponent)
		.def("set_main_multipole_strength", py::overload_cast<const double>(&tse::ClassicalMagnet::setMainMultipoleStrength))
		.def("set_main_multipole_strength", py::overload_cast<const tsc::cdbl>(&tse::ClassicalMagnet::setMainMultipoleStrength))
		.def("propagate", py::overload_cast<tsc::ConfigType&, gtpsa::ss_vect<double>&>      (&tse::ClassicalMagnet::propagate), pass_d_doc)
		.def("propagate", py::overload_cast<tsc::ConfigType&, gtpsa::ss_vect<tps>&>         (&tse::ClassicalMagnet::propagate), pass_tps_doc)
		.def("propagate", py::overload_cast<tsc::ConfigType&, gtpsa::ss_vect<gtpsa::tpsa>&> (&tse::ClassicalMagnet::propagate), pass_tps_doc)
		.def(py::init<const Config &>());

	py::class_<tse::QuadrupoleType, std::shared_ptr<tse::QuadrupoleType>>(m, "Quadrupole", cm)
		.def(py::init<const Config &>());
	py::class_<tse::SextupoleType, std::shared_ptr<tse::SextupoleType>> (m, "Sextupole", cm)
		.def(py::init<const Config &>());
	py::class_<tse::OctupoleType, std::shared_ptr<tse::OctupoleType>>(m, "Octupole", cm)
		.def(py::init<const Config &>());
	py::class_<tse::BendingType, std::shared_ptr<tse::BendingType>>(m, "Bending", cm)
		.def(py::init<const Config &>());
	py::class_<tse::HorizontalSteererType, std::shared_ptr<tse::HorizontalSteererType>>(m, "HorizontalSteerer", cm)
		.def(py::init<const Config &>());
	py::class_<tse::VerticalSteererType, std::shared_ptr<tse::VerticalSteererType>>(m, "VerticalSteerer", cm)
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
