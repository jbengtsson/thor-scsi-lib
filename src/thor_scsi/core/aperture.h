#ifndef _THOR_SCSI_CORE_APERTURE_H_
#define _THOR_SCSI_CORE_APERTURE_H_ 1
#include <cmath>
#include <ostream>

namespace thor_scsi::core {
	/**
	 * @brief: a geometric boundary condition
	 */
	class TwoDimensionalAperture{
	public:

		/**
		 * @brief similar to a bool: positive: within the aperture, negative outside
		 *
		 * It basically returns the distance from the aperture. It is
		 * doable for a circular aperture, but solvable for other
		 * apertures?
		 * Perhaps a target: warrent the sign, focus in the calculation
		 * to be precise if close to the boundary
		 *
		 * Rationale: support optimisations that want to get the beam
		 * close to a boundary. For centre, optimisation use the center
		 */
		virtual double isWithin(const double x, const double y) const = 0;
		virtual void show(std::ostream&, int level) const;

		/*
		 * Same naming convention as for CellVoid
		 */
		virtual const char * type_name(void) const = 0;
	};

	inline
	std::ostream& operator<<(std::ostream& strm, const TwoDimensionalAperture& aperture)
	{
		aperture.show(strm, 0);
		return strm;
	}


}
#endif /* _THOR_SCSI_APERTURE_H_ */
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
