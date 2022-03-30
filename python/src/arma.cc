#include <pybind11/pybind11.h>
#include <armadillo>
#include "thor_scsi.h"

namespace py = pybind11;

void py_thor_scsi_init_arma(py::module &m)
{
	//arma::mat;

  py::class_<arma::mat>(m, "Matrix", py::buffer_protocol())
	  .def_buffer([](arma::mat &mat) -> py::buffer_info {
			py::buffer_info r;
			r.ptr = mat.memptr();         /* Pointer to buffer */
			r.itemsize = sizeof(double); /* Size of one scalar */
			r.format = py::format_descriptor<double>::format(); /* Python struct-style format descriptor */
			r.ndim = 2;
			r.shape = { static_cast<long int>(mat.n_rows),
				    static_cast<long int>(mat.n_cols)
			};/* Number of dimensions */
			r.strides = { sizeof(double) * static_cast<long int>(mat.n_cols),           /* Strides (in bytes) for each index */
				      sizeof(double) };

			return r;
		      });
}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
