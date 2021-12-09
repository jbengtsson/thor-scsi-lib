#ifndef _THOR_SCSI_PY_H_
#define _THOR_SCSI_PY_H_ 1

#include <pybind11/pybind11.h>
#include <string>
#include <sstream>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>

namespace py = pybind11;

template<typename T>
void declare_field(py::module &scsi, const std::string &typestr) {
  using Class = ss_vect<T>;
  std::string pyclass_name = std::string("ss_vect") + typestr;
  py::class_<Class>(scsi, pyclass_name.c_str(), py::buffer_protocol(),
		    py::dynamic_attr())
    .def(py::init<>())
    // Python print redirect.
    .def("__repr__", [](const ss_vect<T> &a) {
		       std::stringstream stream;
		       stream << a;
		       return stream.str();
		     })
    .def("print",      &Class::print)
    .def("__getitem__", [](const ss_vect<T> &a, const int k) {
			  return a[k];
			})
    .def("__setitem__", [](ss_vect<T> &a, const int k, const T &x) {
			  a[k] = x;
			})
    .def("zero",       &Class::zero)
    .def("identity",   &Class::identity);
}

void py_thor_scsi_init_enums(py::module_ &m);
void py_thor_scsi_init_elements(py::module_ &m);
void py_thor_scsi_init_config_type(py::module_ &m);
void py_thor_scsi_init_lattice(py::module_ &m);

#endif /* _THOR_SCSI_PY_H_ */
