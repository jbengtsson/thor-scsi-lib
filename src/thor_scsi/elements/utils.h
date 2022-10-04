#ifndef _THOR_SCSI_ELEMENTS_UTILS_H_
#define _THOR_SCSI_ELEMENTS_UTILS_H_ 1
#include <cmath>

namespace thor_scsi::elements{
	template <typename T>
	inline T cube(T v){ return v * v * v;}

	inline double degtorad(const double val){
		const double scale = M_PI/180e0;
		return val * scale;
	}

	/**
	 * @brief third cubic root solved by pow function
	 *
	 * Currently will only accept values larger than zero.
	 * Could be extended to negative values for a if necessary
	 */
	double thirdroot(const double a);

}
#endif /* _THOR_SCSI_ELEMENTS_UTILS_H_ */
