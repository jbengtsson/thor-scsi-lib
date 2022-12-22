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
		// J.B. 23-05-22: enabled cst.
		.def("cst",          &ss_vect<T>::cst, "return constant term")
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
		.def("__setitem__",  [](ss_vect<T> &a, const int k, const double &x){
					     ss_vect_check_index(k);
					     a.at(k) = x;
				     })
		.def("__copy__",  [](ss_vect<T> &self) {
					  return ss_vect<T>(self);
				  })
		/*
		.def(py::self += py::self)
		.def(py::self -= py::self)
		.def(py::self *= double())
		.def(py::self *= tps())
		.def(py::self + py::self)
		.def(py::self - py::self)
		*/
		// J.B. 23-05-22: added missing ss_vect<> operators.
		.def(py::self += py::self)
		.def(py::self += ss_vect<double>())

		.def(py::self -= py::self)
		.def(py::self -= ss_vect<double>())

		.def(py::self + py::self)
		.def(py::self + ss_vect<double>())
		.def(ss_vect<double>() - py::self)

		.def(py::self - py::self)
		.def(py::self - ss_vect<double>())
		.def(ss_vect<double>() - py::self)
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

/**
 * @brief check that the coordinate is within range
 *
 * Consider moving to ss_vect header ...
 */
static void ss_vect_tps_check_coor(const int coor)
{
	if(coor < 1){
		std::ostringstream strm;
		strm << "Coordinate (first index) must be at least 1 but was " << coor;
		throw py::index_error(strm.str());
	}
	if(coor > ps_dim) {
		std::ostringstream strm;
		strm << "Coordinate (first index) must be at less than "<< ps_dim << " but was " << coor;
		throw py::index_error(strm.str());
	}
}

void py_thor_scsi_init_tps(py::module &m)
{

	py::enum_<spatial_ind>(m, "spatial_index")
		.value("X", spatial_ind::X_)
		.value("Y", spatial_ind::Y_)
		.value("Z", spatial_ind::Z_);

	// does not follow standard convention
	py::enum_<phase_space_ind>(m, "phase_space_index_internal")
		.value("x",     phase_space_ind::x_)
		.value("px",    phase_space_ind::px_)
		.value("y",     phase_space_ind::y_)
		.value("py",    phase_space_ind::py_)
		.value("delta", phase_space_ind::delta_)
		.value("ct",    phase_space_ind::ct_);


	// Polymorphic number class.
	// Classes.
	py::class_<tps>(m, "tps")
		.def(py::init<>())
		.def(py::init<const double>())
		.def(py::init<const double, const int>())
		//.def(py::init<const double, const std::array<const int, >>())
		.def("cst",      &tps::cst, "constant term")
		//.def("ps",       &tps::ps, "phase space")
		.def("__str__",  &tps::pstr)
		.def("__repr__", &tps::repr)
		.def("__getitem__", [](const tps &self, const long int idx){
					    tps_check_index(idx);
					    return self[idx];
				    })
		.def("peek", [](const tps &self, tpsa_index idx){
				     // I think that check is not sufficient
				     // code crashes
				     for(auto k: idx){
					     tps_check_index(k);
				     }
				     return self[idx];
			     }, "get value at this set of 7 indices")
		.def("pook", [](tps &self, tpsa_index idx, double val){
				     // I think that check is not sufficient
				     // code crashes
				     for(auto k: idx){
					     tps_check_index(k);
				     }
				     self.pook(idx, val);
			     }, "set value at this set of 7 indices")
		.def(py::self += py::self)
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

	auto ss_vect_double = declare_ss_vect<double>      (m, "ss_vect_double");
	auto ss_vect_tps    = declare_ss_vect<tps>         (m, "ss_vect_tps"   );
	ss_vect_tps.def(py::self * ss_vect<double>());

	// m.def("partialInverse", &PInv, "partial inverse depending on the numbers of freedoms");
	m.def("partialInverse",  py::overload_cast<const ss_vect<tps>&, const tpsa_index&>(&PInv), "partial inverse depending on the numbers of freedoms");
	m.def("select_subpart", &select_subpart, "support experiment to understand the code");
	m.def("inverse", &Inv, "full inverse");
	m.def("xabs", [](int n, ss_vect<double>&x) -> double {
			      if (n > ss_dim){
				      std::stringstream strm;
				      strm << "Index n " << n
					   << " larger than ss_dim = " << ss_dim;
				      throw std::out_of_range(strm.str());
			      }
			      return xabs(n, x);
		 },
		"computes vector norm", py::arg("n"), py::arg("x"));

#if 0
	ss_vect_tps
		.def("__getitem__", [](ss_vect<tps> &ps, py::sequence &seq) -> double {
					    py::int_ coord= seq[0], term =  seq[1];
					    const int c = coord, t = term;
					    // std::cerr << "accessing element coord " << c <<  " term " << t << std::endl;
					    ss_vect_check_index(c);
					    tps_check_index(t);
					    // I don't know why this reduces c by 1 internally
					    return get_m_ij_save(ps, c + 1, t);
				    })
		.def("__setitem__", [](ss_vect<tps> &ps, py::sequence &seq, const double val) {
					    py::int_ coord= seq[0], term =  seq[1];
					    const int c = coord, t = term;
					    // std::cerr << "accessing element coord " << c <<  " term " << t << std::endl;
					    ss_vect_check_index(c);
					    ss_vect_tps_check_coor(c);
					    tps_check_index(t);
					    put_m_ij_save(ps, c + 1, t, val);
				    });
#endif

	//m.def("lists_to_ss_vect_tps", mattomap_save);
	m.def("ss_vect_tps_to_mat", &maptomat);
	m.def("mat_to_ss_vect_tps", &mattomap_check);
	m.def("vec_mat_to_ss_vect", &vecmattomap);
	m.attr("ps_dim") = ps_dim;
	m.attr("ss_dim") = ss_dim;

}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
