#ifndef _THOR_SCSI_CORE_CONSTANTS_H_
#define _THOR_SCSI_CORE_CONSTANTS_H_ 1

/*
 * Todo:
 *    Review if math.h should be included in global namespace just for MPI
 *
 */
#include <math.h>
namespace thor_scsi{
	namespace core {
/**
 *
 * Physical constants
 *
 * Todo:
 *     Consider if replacing with GSL constants
 */

/**
 *  speed of light in vacuum
 */
		const double  c0    = 2.99792458e8;
/**
 *  electron charge
 */
		const double  q_e   = 1.602e-19;
		const double  m_e   = 0.51099906e6;             ///< electron rest mass [eV/c^2]
		const double  mu_0  = 4.0*M_PI*1e-7;            ///< permittivity of free space
		const double  eps_0 = 1.0/(c0*c0*mu_0);         ///< permeability of free space
		const double  r_e   = q_e/(4.0*M_PI*eps_0*m_e); ///< classical electron radius
		const double  h_bar = 6.58211899e-16;          ////< reduced Planck constant  [eV s]
  }
}
#endif /* _THOR_SCSI_CORE_CONSTANTS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
