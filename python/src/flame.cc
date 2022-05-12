#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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
		.def(py::init<>());
		// .def("get", &Config::get)

	//GLPSParser parse;

	py::class_<GLPSParser>(m, "GLPSParser")
		//.def("parse_file",  py::overload_cast<const char *, const bool>(&GLPSParser::parse_file))
		//.def("parse_file",  py::overload_cast<const char *, FILE *, const char *>(&GLPSParser::parse_file))
		.def("parse_byte",  py::overload_cast<const std::string&, const char *>(&GLPSParser::parse_byte))
		.def(py::init<>());

}

PYBIND11_MODULE(flame, m) {
	m.doc() = "Machinery of FLAME to construct accelerator (but without FLAMES)";
	py_flame_init(m);
}

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
