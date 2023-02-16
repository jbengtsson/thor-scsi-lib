#include <pybind11/pybind11.h>
#include <thor_scsi/elements/standard_aperture.h>

namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace py = pybind11;

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

void py_thor_scsi_init_aperture(py::module &m)
{

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
}
