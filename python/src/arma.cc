#include <pybind11/pybind11.h>
#include <armadillo>
#include "thor_scsi.h"

namespace py = pybind11;


void py_thor_scsi_init_arma(py::module &m)
{
	//arma::mat;

	py::class_<arma::mat>(m, "Matrix", py::buffer_protocol())
		.def("__repr__", [](arma::mat &mat) -> std::string {
					 std::stringstream strm;
					 mat.print(strm, "<tpsa array>");
					 return strm.str();
				 })
		//.def("__str__", &PyArma::pstr)
		.def_buffer([](arma::mat &mat) -> py::buffer_info {
				    size_t n_cols = static_cast<size_t>(mat.n_cols);
				    py::buffer_info r;
				    r.ptr = mat.memptr();         /* Pointer to buffer */
				    r.itemsize = sizeof(double); /* Size of one scalar */
				    r.format = py::format_descriptor<double>::format(); /* Python struct-style format descriptor */
				    r.ndim = 2;
				    r.shape = { static_cast<py::ssize_t>(mat.n_rows),
						static_cast<py::ssize_t>(mat.n_cols)
				    };/* Number of dimensions */
				    r.strides = {
					    /* Strides (in bytes) for each index */
					    static_cast<py::ssize_t>(sizeof(double)),
					    static_cast<py::ssize_t>(sizeof(double) * n_cols)};
				    return r;
			    });
}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
