#ifndef _THOR_SCSI_ELEMENTS_CONSTANTS_H_
#define _THOR_SCSI_ELEMENTS_CONSTANTS_H_ 1
#include <cmath>
#include <thor_scsi/elements/utils.h>
// for definition of sqr ....
#include <tps/tps_type.h>

namespace thor_scsi::elements{
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
	const double speed_of_light = 2.99792458e8;
	const double m_e   = 0.51099906e6;             ///< electron rest mass [eV/c^2].
	const double h_bar = 6.58211899e-16;           //< reduced Planck constant [eVs].

        /**
	 * See Sands ....
	 *
	 * @todo: remove 1e-9 in calculation of the result
	 *        if so remeve the energy scale in radiation_delegate
	 */
	inline double compute_C_gamma(void)  {
		const double q_e   = 1.602e-19;                ///< electron charge.
		const double mu_0  = 4.0*M_PI*1e-7;            // permittivity of free space.
		const double eps_0 = 1.0/(sqr(speed_of_light)*mu_0);   // permeability of free space.
		const double r_e   = q_e/(4.0*M_PI*eps_0*m_e); ///< classical electron radius.

		double r =  4e0 * M_PI * r_e / (3e0 * cube(m_e));
		return r;
	}

	inline double compute_C_q(void) {
		double c0 = speed_of_light;
		const double C_u = 55e0/(24e0*sqrt(3e0));
		return 3e0*C_u*h_bar*c0/(4e0*m_e);
	}

	const double C_gamma = compute_C_gamma();
	const double C_q = compute_C_q();


}
#endif /* _THOR_SCSI_ELEMENTS_CONSTANTS_H_ */
