#include <pybind11/pybind11.h>
#include <thor_scsi/version.h>
#include <thor_scsi/elements/constants.h>
#include "thor_scsi.h"
#include <string>
#include <sstream>

namespace py = pybind11;
namespace tse = thor_scsi::elements;


PYBIND11_MODULE(lib, m) {
    m.doc() = "Self-Consistent Symplectic Integrator (SCSI)";

    py_thor_scsi_init_tps(m);
    py_thor_scsi_init_field_interpolation(m);
    py_thor_scsi_init_custom(m);
    py_thor_scsi_init_elements(m);
    py_thor_scsi_init_aperture(m);
    py_thor_scsi_init_radiation(m);
    py_thor_scsi_init_observers(m);
    py_thor_scsi_init_accelerator(m);
    py_thor_scsi_init_config_type(m);
    // py_thor_scsi_init_lattice(scsi);


    // Beamline class.

    // Constants.
    //scsi.attr("HOMmax") = HOMmax;
    m.attr("c0")     = tse::speed_of_light;
    // scsi.attr("q_e")    = tsc::q_e;
    // scsi.attr("m_e")    = tsc::m_e;
    // scsi.attr("mu_0")   = tsc::mu_0;
    // scsi.attr("eps_0")  = tsc::eps_0;
    // scsi.attr("r_e")    = tsc::r_e;
    // scsi.attr("h_bar")  = tsc::h_bar;
    /*
    */
    // Enums.


}
