#ifndef _THOR_SCSI_CORE_FIELD_INTERPOLATION_H_
#define _THOR_SCSI_CORE_FIELD_INTERPOLATION_H_ 1
#include <ostream>
namespace thor_scsi::core {
  	/**
	 * Pure virtual class of multipoles
	 *
	 * External code just requires to calculate multipoles
	 *
	 * Todo: move to API part?
	 */
	struct Field2DInterpolation{
		virtual ~Field2DInterpolation(){};
		/**
		 * @brief interpolate field at position x and y
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * Args:
		 *     x: first coordinate (typically horizontal position)
		 *     y: second coordinate  (typically vertical position)
		 *     Bx: pointer to variable where field component `x` will be stored
		 *     By: pointer to variable where field component `y` will be stored
		 *
		 * Warning:
		 *     Note: * thor_scsi* is using a left hand coordinate system.
		 *     This function here is expecting all coordinates in a left
		 *     handed coordiante system
		 *
		 * \endverbatim
		 */
		virtual inline void field(const double x, const double y, double *Bx, double *By) const = 0;

		/**
		 * @brief interpolate the gradient at the current position
		 *
		 * \verbatim embed:rst:leading-asterisk
		 * \endverbatim
		 */
		virtual void gradient(const double x, const double y, double *Gx, double *Gy) const = 0;
		virtual void show(std::ostream&, int level) const = 0;

	};

	inline
	std::ostream& operator<<(std::ostream& strm, const Field2DInterpolation& intp)
	{
		intp.show(strm, 0);
		return strm;
	}

}
#endif /* _THOR_SCSI_FIELD_INTERPOLATION_H_ */
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
