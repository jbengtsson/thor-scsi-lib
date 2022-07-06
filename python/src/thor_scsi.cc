#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
// #include <pybind11/operators.h>
#include <vector>
#include <tuple>

#include <thor_scsi/version.h>
#include <tps/enums.h>
// #include <thor_scsi/core/defines.h>
#include <thor_scsi/elements/constants.h>
#include <thor_scsi/core/cells.h>
// #include <thor_scsi/core/elements_enums.h>
#include "thor_scsi.h"
#include <string>
#include <sstream>

namespace py = pybind11;
namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;


// Polymorphic number template class wrapper.



/*
 * get_eps_x: provide the memory and return it to python
 *            TODO: check for memory leaks!
 */


PYBIND11_MODULE(lib, scsi) {
    scsi.doc() = "Self-Consistent Symplectic Integrator (SCSI)";

    py_thor_scsi_init_arma(scsi);
    py_thor_scsi_init_tps(scsi);
    py_thor_scsi_init_elements(scsi);
    py_thor_scsi_init_accelerator(scsi);
    py_thor_scsi_init_config_type(scsi);
    // py_thor_scsi_init_lattice(scsi);


    // Beamline class.

    // Constants.
    //scsi.attr("HOMmax") = HOMmax;
    scsi.attr("c0")     = tse::speed_of_light;
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
