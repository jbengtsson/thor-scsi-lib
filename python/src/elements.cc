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
#include <thor_scsi/elements/cavity.h>
#include <thor_scsi/elements/standard_observer.h>
#include <thor_scsi/elements/standard_aperture.h>


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

class PyField2DInterpolation: public tsc::Field2DInterpolation{
public:
	using tsc::Field2DInterpolation::Field2DInterpolation;
#if 0
	void field_py(const py::array<double, 2> pos, std::array<double, 2> field) const {
		PYBIND11_OVERRIDE_PURE(void, PyField2DInterpolationIntermediate, field_py, pos, field);
	}
#endif
	virtual void field_py(const py::array_t<double> &t_pos, py::array_t<double> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, field_py, t_pos, t_field);
	}
	virtual void field_py(const py::array_t<tps> &t_pos, py::array_t<tps> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, field_py, t_pos, t_field);
	}
	virtual void gradient_py(const std::array<double, 2> &t_pos, std::array<double, 2> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, gradient_py, t_pos, t_field);
	}
	virtual void gradient_py(const std::array<tps, 2> &t_pos, std::array<tps, 2> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, gradient_py, t_pos, t_field);
	}
	virtual void gradient_py(const std::array<tps, 2> &t_pos, std::array<double, 2> &t_field) const {
		PYBIND11_OVERRIDE(void, PyField2DInterpolation, gradient_py, t_pos, t_field);
	}
private:
	template<typename T>
	inline void _field(const T x, const T y, T *Bx, T *By) const {
		auto pos = py::array_t<T>({2}), t_field = py::array_t<T>({2});
		T *pos_p = static_cast<T *>(pos.request().ptr),
			*field_p = static_cast<T *>(t_field.request().ptr);
		pos_p[0] = x;
		pos_p[1] = y;
		this->field_py(pos, t_field);
		*Bx = field_p[0];
		*By = field_p[1];

		std::cerr << "By " << *By <<  ",  Bx" <<  *Bx << std::endl;
	}
	template<typename T>
	inline void _gradient(const T x, const T y, T *Bx, T *By) const {
	}

public:
	void field(const double x, const double y, double *Bx, double *By) const override {
		this->_field(x, y, Bx, By);
	}
	void field(const tps x, const tps y, tps *Bx, tps *By) const override {
		//this->_field(x, y, Bx, By);
		throw std::runtime_error("field with tps to python not implemented");
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
		.def("__str__",                  &tse::StandardObserver::pstr)
		.def("__repr__",                 &tse::StandardObserver::repr)
		.def("getObservedName",          &tse::StandardObserver::getObservedName)
		.def("getObservedIndex",         &tse::StandardObserver::getObservedIndex)
		.def("hasPhaseSpace",            &tse::StandardObserver::hasPhaseSpace)
		.def("getPhaseSpace",            &tse::StandardObserver::getPhaseSpace)
		.def("hasTruncatedPowerSeries",  &tse::StandardObserver::hasTruncatedPowerSeries)
		.def("getTruncatedPowerSeries",  &tse::StandardObserver::getTruncatedPowerSeries)
		.def("reset",                    &tse::StandardObserver::reset)
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


	py::class_<tsc::ElemType,  PyElemType, std::shared_ptr<tsc::ElemType>> elem_type(m, "ElemType");
	elem_type.def("__str__",       &tsc::ElemType::pstr)
		.def("__repr__",       &tsc::ElemType::repr)
		.def_readonly("name",  &tsc::ElemType::name)
		.def_readonly("index", &tsc::ElemType::index)
		.def("getLength",      &tse::ElemType::getLength)
		.def("setLength",      &tse::ElemType::setLength)
		.def("getObserver",    &tse::ElemType::observer)
		.def("setObserver",    &tse::ElemType::set_observer)
		.def("getAperture",    &tse::ElemType::getAperture)
		.def("setAperture",    &tse::ElemType::setAperture)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<double>&>(&tse::ElemType::pass), pass_d_doc)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ss_vect<tps>&>(&tse::ElemType::pass),    pass_tps_doc)
		.def(py::init<const Config &>());


	py::class_<tse::DriftType, std::shared_ptr<tse::DriftType>>(m, "Drift", elem_type)
		.def(py::init<const Config &>());

	py::class_<tse::RadiationDelegateInterface, PyRadDelInt, std::shared_ptr<tse::RadiationDelegateInterface>> rad_del_int(m, "RadiationDelegateInterface");
	rad_del_int
		.def("__repr__",        &tse::RadiationDelegateInterface::repr)
		.def(py::init<>());

	py::class_<tse::RadiationDelegate, PyRadDel, std::shared_ptr<tse::RadiationDelegate>>(m, "RadiationDelegate", rad_del_int)
		.def("reset",             &tse::RadiationDelegate::reset)
		.def("getCurlydHx",       &tse::RadiationDelegate::getCurlydHx)
		.def("getDelegatorName",  &tse::RadiationDelegate::getDelegatorName)
		.def("getDelegatorIndex", &tse::RadiationDelegate::getDelegatorIndex)
		.def("view",
		     py::overload_cast<const tse::ElemType&, const ss_vect<double> &, const enum tsc::ObservedState, const int>(&tse::RadiationDelegate::view))
		.def("view",
		     py::overload_cast<const tse::ElemType&, const ss_vect<double> &, const enum tsc::ObservedState, const int>(&tse::RadiationDelegate::view))
		.def(py::init<>());

	py::class_<tse::RadiationDelegateKickInterface, PyRadDelKickInt, std::shared_ptr<tse::RadiationDelegateKickInterface>> rad_del_kick_int(m, "RadiationDelegateKickInterface");
	rad_del_int
		//.def("__repr__", &tse::RadiationDelegateKickInterface::repr)
		.def(py::init<>());

	// Why did it not end up at Johan?
	py::class_<tse::RadiationDelegateKick, PyRadDelKick, std::shared_ptr<tse::RadiationDelegateKick>>(m, "RadiationDelegateKick", rad_del_kick_int)
		.def("reset",                              &tse::RadiationDelegateKick::reset)
		.def("getCurlydHx",                        &tse::RadiationDelegateKick::getCurlydHx)
		.def("getDelegatorName",                   &tse::RadiationDelegateKick::getDelegatorName)
		.def("getDelegatorIndex",                  &tse::RadiationDelegateKick::getDelegatorIndex)
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
		.def("getMultipole",     &tsc::TwoDimensionalMultipoles::getMultipole)
		.def("setMultipole",     &tsc::TwoDimensionalMultipoles::setMultipole)
		.def("applyRollAngle",   &tsc::TwoDimensionalMultipoles::applyRollAngle)
		.def("__str__",          &tsc::TwoDimensionalMultipoles::pstr)
		.def("__repr__",         &tsc::TwoDimensionalMultipoles::repr)
		.def("__len__",          &tsc::TwoDimensionalMultipoles::size)
		.def("getCoefficients",  &tsc::TwoDimensionalMultipoles::getCoeffsConst)
		//.def("field", static_cast<void (tsc::Field2DInterpolation::*)(const double, const double, double *, double *) const>(&tsc::TwoDimensionalMultipoles::field))
		//.def("field", static_cast<void (tsc::Field2DInterpolation::*)(const tps, const tps, tps *, tps *) const>(&tsc::TwoDimensionalMultipoles::field))
		.def("applyTranslation", py::overload_cast<const tsc::cdbl>(&tsc::TwoDimensionalMultipoles::applyTranslation))
		.def("applyTranslation", py::overload_cast<const double, const double>(&tsc::TwoDimensionalMultipoles::applyTranslation))
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
		.def("getFieldInterpolator",        &tse::FieldKick::getFieldInterpolator)
		.def("setFieldInterpolator",        &tse::FieldKick::setFieldInterpolator)

		.def(py::init<const Config &>());

	py::class_<tse::MpoleType, std::shared_ptr<tse::MpoleType>> mpole_type(m, "Mpole", field_kick);
	mpole_type
		.def(py::init<const Config &>());

	py::class_<tse::ClassicalMagnet, PyClassicalMagnet, std::shared_ptr<tse::ClassicalMagnet>> cm(m, "ClassicalMagnet", mpole_type);
	cm
		.def("getMultipoles",            &tse::ClassicalMagnet::getMultipoles)
		//.def("getBegninMultipoles",      &tse::ClassicalMagnet::getBegninMultipoles)
		// .def("setMultipoles",            &tse::ClassicalMagnet::setMultipoles)
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
