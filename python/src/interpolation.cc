#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <thor_scsi/core/field_interpolation.h>
#include <thor_scsi/core/multipoles.h>
#include <complex>

namespace py = pybind11;
namespace tsc = thor_scsi::core;

typedef std::complex<double> cdbl;
typedef std::complex<double> cdbl_intern;

template<class C>
class PyField2DInterpolation: public tsc::Field2DInterpolationKnobbed<C> {
public:
	using base = tsc::Field2DInterpolationKnobbed<C>;

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
	void field(const double& x, const double& y, double *Bx, double *By) const override {
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

	void field(const tps& x, const tps& y, tps *Bx, tps *By) const override {
		//this->_field(x, y, Bx, By);
		std::array<tps, 2> pos = {x, y};
		std::array<tps, 2> t_field;
		this->field_py(pos, t_field);
		*Bx = t_field[0];
		*By = t_field[1];

	}
	void field(const gtpsa::tpsa& x, const gtpsa::tpsa& y, gtpsa::tpsa *Bx, gtpsa::tpsa *By) const override {
		//this->_field(x, y, Bx, By);
		std::array<gtpsa::tpsa, 2> pos = {x, y};
		std::array<gtpsa::tpsa*, 2> t_field = {Bx, By};
		// this->field_py(pos, t_field);
		//*Bx = t_field[0];
		//*By = t_field[1];

	}
	void gradient(const double& x, const double& y, double *Gx, double *Gy) const  override{
		PYBIND11_OVERRIDE_PURE(void, base, gradient, x, y, Gx, Gy);
	}
	void gradient(const tps& x, const tps& y, tps *Gx, tps *Gy) const   override {
		PYBIND11_OVERRIDE_PURE(void, base, gradient, x, y, Gx, Gy);
	}
	void gradient(const tps& x, const tps& y, double *Gx, double *Gy) const  override{
		PYBIND11_OVERRIDE_PURE(void, base, gradient, x, y, Gx, Gy);
	}
	void gradient(const gtpsa::tpsa& x, const gtpsa::tpsa& y, gtpsa::tpsa *Gx, gtpsa::tpsa *Gy) const   override {
		PYBIND11_OVERRIDE_PURE(void, base, gradient, x, y, Gx, Gy);
	}
	void gradient(const gtpsa::tpsa& x, const gtpsa::tpsa& y, double *Gx, double *Gy) const  override{
		PYBIND11_OVERRIDE_PURE(void, base, gradient, x, y, Gx, Gy);
	}

	void show(std::ostream& strm, int level) const override {
		const std::string txt = (level <= 1) ? this->pstr() : this->repr();
		strm << txt;
	}

	std::string repr(void) const {
		PYBIND11_OVERRIDE_NAME(std::string, base, "__repr__", repr);
	}

	std::string pstr(void) const {
		PYBIND11_OVERRIDE_NAME(std::string, base, "__str__", pstr);
	}
};


template<typename Types, typename Class>
void add_methods_interpolation(py::class_<Class> t_mapper)
{
	t_mapper
		.def("__repr__", &Class::repr)
		.def("__str__",  &Class::pstr)
		;
}

template<typename Types, typename Class>
void add_methods_multipoles(py::class_<Class> t_mapper)
{
	t_mapper
		.def("get_multipole",    &Class::getMultipole)
		.def("set_multipole",    &Class::setMultipole)
		.def("apply_roll_angle", &Class::applyRollAngle)
		.def("__str__",          &Class::pstr)
		.def("__repr__",         &Class::repr)
		.def("__len__",          &Class::size)
		.def("get_coefficients", &Class::getCoeffsConst)
		.def(py::self += py::self)
		.def(py::self + py::self)
	    #if 0
		.def(py::self += std::complex<double>())
		.def(py::self + std::complex<double>())
		.def(py::self *= std::complex<double>())
		.def(py::self * std::complex<double>())
	    #endif
		.def(std::complex<double>() * py::self)
		.def(std::complex<double>() + py::self)
		.def(py::self *= std::vector<typename Types::complex_type>())
		.def(py::self * std::vector<typename Types::complex_type>())
		.def("field", [](const Class& inst, typename Types::complex_type& z) {
			return inst._cfield(z);
		})
		;
}

void py_thor_scsi_init_field_interpolation(py::module &m) {

	py::class_<
		tsc::Field2DInterpolation,
		PyField2DInterpolation<tsc::StandardDoubleType>,
		std::shared_ptr<tsc::Field2DInterpolation>
		>
		field2dintp (m, "Field2DInterpolation");
	add_methods_interpolation<tsc::StandardDoubleType, tsc::Field2DInterpolation>(field2dintp);

	py::class_<
		tsc::Field2DInterpolationDependence,
		PyField2DInterpolation<tsc::TpsaVariantType>,
		std::shared_ptr<tsc::Field2DInterpolationDependence>
		> field2dintpvar(m, "Field2DInterpolationDependence");
	add_methods_interpolation<tsc::StandardDoubleType, tsc::Field2DInterpolation>(field2dintpvar);
	field2dintpvar
            .def(py::init<>());


	field2dintp
            // can I write the following lines in a c++-14 compatible way
            //.def("field",        py::overload_cast<const double, const double, double *, double *>(&tsc::Field2DInterpolation::field))
            .def("field", [](tsc::Field2DInterpolation &intp, const py::array_t<double> t_pos) -> py::array_t<double> {
		    py::buffer_info t_pos_buf = t_pos.request();
		    if (t_pos_buf.ndim != 1) {
			throw std::runtime_error("position must be a vector");
		    }
		    if (t_pos_buf.size != 2) {
			throw std::runtime_error("position must be a vector of size 2");
		    }
		    auto t_field = py::array_t<double>(t_pos_buf.size);
		    py::buffer_info t_field_buf = t_field.request();
		    double *pos_p = static_cast<double *>(t_pos_buf.ptr);
		    double *field_p = static_cast<double *>(t_field_buf.ptr);

		    intp.field(pos_p[0], pos_p[1], &field_p[0], &field_p[1]);
		    // std::cerr << "Returning field " << field_p[0] << ", " << field_p[1] << "\n";
		    return t_field;
		})
            .def("field", [](tsc::Field2DInterpolation &intp, const std::array<tps, 2> t_pos,
                             std::array<tps, 2> t_field) -> void {
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
	py::class_<
	    tsc::TwoDimensionalMultipolesKnobbed<tsc::StandardDoubleType>,
	    std::shared_ptr<tsc::TwoDimensionalMultipolesKnobbed<tsc::StandardDoubleType>>
	    >
	    multipoles_knobbed(m, "_TwoDimensionalMultipolesKnobbedDouble", field2dintp //, py::buffer_protocol()
		);
	add_methods_multipoles<tsc::StandardDoubleType, tsc::TwoDimensionalMultipolesKnobbed<tsc::StandardDoubleType>>(multipoles_knobbed);
	py::class_<tsc::TwoDimensionalMultipoles, std::shared_ptr<tsc::TwoDimensionalMultipoles>>
	    multipole(m, "TwoDimensionalMultipoles", multipoles_knobbed //, py::buffer_protocol()
		);

	add_methods_multipoles<tsc::StandardDoubleType, tsc::TwoDimensionalMultipoles>(multipole);
	multipole
	    //.def("field", static_cast<void (tsc::Field2DInterpolation::*)(const double, const double, double *, double *) const>(&tsc::TwoDimensionalMultipolesKnobbed::field))
	    //.def("field", static_cast<void (tsc::Field2DInterpolation::*)(const tps, const tps, tps *, tps *) const>(&tsc::TwoDimensionalMultipolesKnobbed::field))
	    .def("apply_translation", py::overload_cast<const cdbl>(&tsc::TwoDimensionalMultipoles::applyTranslation))
	    .def("apply_translation",
		 py::overload_cast<const double, const double>(&tsc::TwoDimensionalMultipoles::applyTranslation))
	    .def(py::init<std::vector<cdbl_intern>>())
	    .def(py::init<const std::complex<double>, const unsigned int>(), "initalise multipoles",
		     py::arg("default_value"), py::arg("h_max") = tsc::max_multipole);
#if 0
	.def_buffer([](tsc::TwoDimensionalMultipoles &muls) -> py::buffer_info {
	    std::vector<cdbl_intern> coeffs = muls.getCoeffs();
	    size_t n = coeffs.size();
	    py::buffer_info r;
                r.ptr = coeffs.data();         /* Pointer to buffer */
                r.itemsize = sizeof(cdbl_intern); /* Size of one scalar */
                r.format = py::format_descriptor<cdbl_intern>::format(); /* Python struct-style format descriptor */
                r.ndim = 1;
                r.shape = {static_cast<py::ssize_t>(n)};/* Number of dimensions */
                r.strides = {
                        /* Strides (in bytes) for each index */
                        static_cast<py::ssize_t>(sizeof(cdbl_intern))
                };
                return r;
            })
	    .def(py::init([](py::array_t<cdbl_intern, py::array::c_style|py::array::forcecast> b) {
                                  /* Request a buffer descriptor from Python */
                                  py::buffer_info info = b.request();

                                  /* Some sanity checks ... */
                                  if (info.format != py::format_descriptor<cdbl_intern>::format())
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

		py::class_<
		    tsc::TwoDimensionalMultipolesKnobbed<tsc::TpsaVariantType>,
		std::shared_ptr<tsc::TwoDimensionalMultipolesKnobbed<tsc::TpsaVariantType>>
	    > multipoles_tpsa_base(m, "_MultipolesBaseTpsa", field2dintpvar);
	        add_methods_multipoles<tsc::TpsaVariantType, tsc::TwoDimensionalMultipolesKnobbed<tsc::TpsaVariantType>>(multipoles_tpsa_base);

		py::class_<tsc::TwoDimensionalMultipolesTpsa, std::shared_ptr<tsc::TwoDimensionalMultipolesTpsa>>
		    multipoles_tpsa(m, "TwoDimensionalMultipolesTpsa", multipoles_tpsa_base //, py::buffer_protocol()
	    );
	    add_methods_multipoles<tsc::TpsaVariantType, tsc::TwoDimensionalMultipolesTpsa>(multipoles_tpsa);
	    multipoles_tpsa
		    .def(py::init<const std::complex<double>, const unsigned int>(), "initalise multipoles",
			 py::arg("default_value"), py::arg("h_max") = tsc::max_multipole)
		    ;


}
