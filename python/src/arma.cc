#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
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
			    })
		//.def(py::init([](py::buffer b) {
		.def(py::init([](py::array_t<double, py::array::c_style|py::array::forcecast> b) {
				      /* Request a buffer descriptor from Python */
				      py::buffer_info info = b.request();

				      /* Some sanity checks ... */
				      if (info.format != py::format_descriptor<double>::format())
					      throw std::runtime_error("Incompatible format: expected a double array!");

				      if (info.ndim != 2){
					      std::stringstream strm;
					      strm << "Incompatible buffer: expected 2 but got "
						   << info.ndim << "dimensions!";
					      throw std::runtime_error(strm.str());
				      }

				      bool need_transpose = false;
				      if(info.strides[0] != sizeof(double)){
					      need_transpose = true;
				      }
				      /*
					std::cerr << "Strides [" << info.strides[0] << ", " << info.strides[1] << "]"
					<< ": need_transpose " << std::boolalpha << need_transpose << std::endl;
				      */
				      auto mat = arma::mat(static_cast<const double *>(info.ptr),
							   info.shape[0], info.shape[1]);
				      if(need_transpose){
					      arma::inplace_trans(mat);
				      }
				      return mat;
			      }));
}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
