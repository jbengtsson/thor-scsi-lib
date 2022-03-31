// #include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <tps/tpsa_lin.h>
#include <tps/enums.h>
#include <string>
#include <sstream>
#include <cassert>
#include "thor_scsi.h"

namespace py = pybind11;


static void ss_vect_check_index(const int index)
{
	std::stringstream strm;
	const int max_index = ps_dim - 1;
	if(index < 0){
		strm << "ss_vect[" << index << "] invalid: < 0!";
		throw py::index_error(strm.str());
	}
	else if (index > max_index){
		strm << "ss_vect[" << index << "] invalid: > " << max_index << "!";
		throw py::index_error(strm.str());
	}
}

template<typename T>
static auto declare_ss_vect(py::module &m, const std::string pyclass_name) {
	using Class = ss_vect<T>;

	py::class_<Class> a_class(m, pyclass_name.c_str()
			  //,py::buffer_protocol(), py::dynamic_attr()
		);
	a_class
		.def("set_identity", &ss_vect<T>::set_identity)
		.def("set_zero",     &ss_vect<T>::set_zero)
		//.def("cst",          &ss_vect<T>::cst, "return constant term")
		/*
		.def("__eq__",       [](const ss_vect<T> &self, (const ss_vect<T> &a)){
					     return self.cst() == a.cst();
				     })
		*/
		.def("__repr__",     &ss_vect<T>::repr)
		.def("__str__",      &ss_vect<T>::pstr)
		.def("__len__",      [](const ss_vect<T> &a) -> int {return ps_dim;})
		.def("__getitem__",  [](const ss_vect<T> &a, const int k){
					     ss_vect_check_index(k);
					     return a.at(k);
				     })
		.def("__setitem__",  [](ss_vect<T> &a, const int k, const T &x){
					     ss_vect_check_index(k);
					     a.at(k) = x;
				     })
		/*
		.def(py::self += py::self)
		.def(py::self -= py::self)
		.def(py::self *= double())
		.def(py::self *= tps())
		.def(py::self + py::self)
		.def(py::self - py::self)
		*/
		//
		.def(py::init<>());
	return a_class;
}

static void tps_check_index(const int index)
{
	std::stringstream strm;
	const int max_index = (nv_tps + 1) - 1;
	if(index < 0){
		strm << "ss_vect[" << index << "] invalid: < 0!";
		throw py::index_error(strm.str());
	}
	else if (index > max_index){
		strm << "ss_vect[" << index << "] invalid: > " << max_index << "!";
		throw py::index_error(strm.str());
	}
}

void py_thor_scsi_init_tps(py::module &m)
{

	py::enum_<spatial_ind>(m, "spatial_ind")
		.value("X_", spatial_ind::X_)
		.value("Y_", spatial_ind::Y_)
		.value("Z_", spatial_ind::Z_);

	py::enum_<phase_space_ind>(m, "phase_space_ind")
		.value("x_",     phase_space_ind::x_)
		.value("px_",    phase_space_ind::px_)
		.value("y_",     phase_space_ind::y_)
		.value("py_",    phase_space_ind::py_)
		.value("delta_", phase_space_ind::delta_)
		.value("ct_",    phase_space_ind::ct_);


	// Polymorphic number class.
	// Classes.

	py::class_<tps>(m, "tps")
		.def(py::init<>())
		.def(py::init<const double>())
		.def(py::init<const double, const int>())
		//.def(py::init<const double, const std::array<const int, >>())
		//.def("cst",      &tps::cst, "constant term")
		.def("__str__",  &tps::pstr)
		.def("__repr__", &tps::repr)

		.def("__getitem__", [](const tps &self, const int k){
					    tps_check_index(k);
					    return self[k];
				    })
		.def(py::self -= py::self)
		.def(py::self *= py::self)
		.def(py::self /= py::self)

		.def(py::self += double())
		.def(py::self -= double())
		.def(py::self *= double())
		.def(py::self /= double())

		.def(py::self * py::self)

		.def(py::self + py::self)
		.def(py::self + double())
		.def(double() + py::self)

		.def(py::self - py::self)
		.def(py::self - double())
		.def(double() - py::self)

		.def(py::self * py::self)
		.def(py::self * double())
		.def(double() * py::self)

		.def(py::self / py::self)
		.def(py::self / double())
		.def(double() / py::self)


		;

	auto ss_vect_tps = declare_ss_vect<tps>(m, "ss_vect_tps");
	auto ss_vect_double = declare_ss_vect<double>(m, "ss_vect_double");

	//m.def("lists_to_ss_vect_tps", mattomap_save);
	m.def("ss_vect_tps_to_mat", &maptomat);
	m.def("mat_to_ss_vect_tps", &mattomap_check);

}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
