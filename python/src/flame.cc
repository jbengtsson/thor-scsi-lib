#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/iostream.h>
#include <flame/core/config.h>
#include <flame/core/base.h>
#include <ostream>
namespace py = pybind11;

void py_flame_init(py::module_ &m)
{
	py::class_<Config>(m, "Config")
		.def("getAny", &Config::getAny)
		.def("setAny", &Config::setAny)
	        .def("__repr__", &Config::repr)
		.def("__iter__", [](const Config &conf) {
			return py::make_iterator(conf.begin(), conf.end());},
			py::keep_alive<0, 1>())
		.def(py::init<>());

	//GLPSParser parse;

	py::class_<GLPSParser>(m, "GLPSParser")
		//.def("parse_file",  py::overload_cast<const char *, const bool>(&GLPSParser::parse_file))
		//.def("parse_file",  py::overload_cast<const char *, FILE *, const char *>(&GLPSParser::parse_file))
		.def("parse_byte",  py::overload_cast<const std::string&, const char *>(&GLPSParser::parse_byte))
		.def(py::init<>());

}

PYBIND11_MODULE(pyflame, m) {
	m.doc() = "Machinery of FLAME to construct accelerator (but without FLAMES)";
	py_flame_init(m);
}

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
